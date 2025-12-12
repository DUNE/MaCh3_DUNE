#include "SampleHandlerAtm.h"

#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wfloat-conversion"
//Standard Record includes
#include "duneanaobj/StandardRecord/StandardRecord.h"
#pragma GCC diagnostic pop

#include <Eigen/Dense>
#include <iomanip>
#include <openssl/evp.h>

SampleHandlerAtm::SampleHandlerAtm(std::string mc_version_, ParameterHandlerGeneric* xsec_cov_, const std::shared_ptr<OscillationHandler>&  Oscillator_) : SampleHandlerFD(mc_version_, xsec_cov_, Oscillator_) {
  KinematicParameters = &KinematicParametersDUNE;
  ReversedKinematicParameters = &ReversedKinematicParametersDUNE;
  
  Initialise();
}

SampleHandlerAtm::~SampleHandlerAtm() {
}

void SampleHandlerAtm::Init() {
  IsELike = Get<bool>(SampleManager->raw()["SampleOptions"]["IsELike"],__FILE__,__LINE__);
  ExposureScaling = Get<double>(SampleManager->raw()["SampleOptions"]["ExposureScaling"],__FILE__,__LINE__);

  if (SampleManager->raw()["SampleOptions"]["EigenFile"]) {
    EigenInputFile = Get<std::string>(SampleManager->raw()["SampleOptions"]["EigenFile"],__FILE__,__LINE__);
    EigenInputFileMD5Sum = Get<std::string>(SampleManager->raw()["SampleOptions"]["EigenFileMD5Sum"],__FILE__,__LINE__);
  } else {
    EigenInputFile = "";
  }
}

void SampleHandlerAtm::SetupSplines() {
  SplineHandler = nullptr;
}

void SampleHandlerAtm::SetupWeightPointers() {
  for (size_t i = 0; i < dunemcSamples.size(); ++i) {
    MCSamples[i].total_weight_pointers.push_back(&(dunemcSamples[i].flux_w));
    MCSamples[i].total_weight_pointers.push_back(MCSamples[i].osc_w_pointer);
    MCSamples[i].total_weight_pointers.push_back(&(MCSamples[i].xsec_w));
    MCSamples[i].total_weight_pointers.push_back(&(ExposureScaling));
  }  
}

template<class Matrix>
inline void WriteBinaryToFile(const std::string& filename, const Matrix& matrix){
    std::ofstream out(filename, std::ios::out | std::ios::binary | std::ios::trunc);
    if(out.is_open()) {
        typename Matrix::Index rows=matrix.rows(), cols=matrix.cols();
        out.write(reinterpret_cast<char*>(&rows), sizeof(typename Matrix::Index));
        out.write(reinterpret_cast<char*>(&cols), sizeof(typename Matrix::Index));
        out.write(reinterpret_cast<const char*>(matrix.data()), rows*cols*static_cast<typename Matrix::Index>(sizeof(typename Matrix::Scalar)) );
        out.close();
    }
    else {
      MACH3LOG_ERROR("Can not write to file:{}",filename);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
}

template<class Matrix>
inline void ReadBinaryFromFile(const std::string& filename, Matrix& matrix){
    std::ifstream in(filename, std::ios::in | std::ios::binary);
    if (in.is_open()) {
        typename Matrix::Index rows=0, cols=0;
        in.read(reinterpret_cast<char*>(&rows),sizeof(typename Matrix::Index));
        in.read(reinterpret_cast<char*>(&cols),sizeof(typename Matrix::Index));
        matrix.resize(rows, cols);
        in.read(reinterpret_cast<char*>(matrix.data()), rows*cols*static_cast<typename Matrix::Index>(sizeof(typename Matrix::Scalar)) );
        in.close();
    }
    else {
      MACH3LOG_ERROR("Can not open binary matrix file:{}",filename);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
}

// Compute MD5 digest of buffer 'data'. Returns true on success and fills 'digest' (16 bytes).
bool compute_md5(const unsigned char* data, size_t len, unsigned char digest[EVP_MAX_MD_SIZE], unsigned int &digest_len) {
    EVP_MD_CTX *ctx = nullptr;
    bool ok = false;

    ctx = EVP_MD_CTX_new();
    if (!ctx) {
        std::cerr << "EVP_MD_CTX_new failed\n";
        return false;
    }

    // Use EVP_md5() — still available under OpenSSL 3.0's default provider.
    if (1 != EVP_DigestInit_ex(ctx, EVP_md5(), nullptr)) {
        std::cerr << "EVP_DigestInit_ex failed\n";
	throw;
        goto done;
    }

    if (len > 0) {
        if (1 != EVP_DigestUpdate(ctx, data, len)) {
            std::cerr << "EVP_DigestUpdate failed\n";
	    throw;
            goto done;
        }
    }

    if (1 != EVP_DigestFinal_ex(ctx, digest, &digest_len)) {
        std::cerr << "EVP_DigestFinal_ex failed\n";
        throw;
        goto done;
    }

    ok = true;
done:
    EVP_MD_CTX_free(ctx);
    return ok;
}

// Read entire file into vector<char>, returns false on error
bool read_file(const std::string &path, std::vector<unsigned char> &out) {
    std::ifstream ifs(path, std::ios::binary | std::ios::ate);
    if (!ifs) return false;
    std::streamsize size = ifs.tellg();
    ifs.seekg(0, std::ios::beg);
    out.resize(static_cast<size_t>(size));
    if (!ifs.read(reinterpret_cast<char*>(out.data()), size)) return false;
    return true;
}

// Convert binary digest to lowercase hex string
std::string to_hex(const unsigned char* data, size_t len) {
    std::ostringstream oss;
    oss << std::hex << std::setfill('0');
    for (size_t i = 0; i < len; ++i) {
        oss << std::setw(2) << static_cast<int>(data[i]);
    }
    return oss.str();
}

bool CompareMD5Sum(std::string KnownSum, std::string FilePathToCheck) {
  std::vector<unsigned char> data;
  read_file(FilePathToCheck,data);
  
  unsigned char digest[EVP_MAX_MD_SIZE];
  unsigned int digest_len = 0;
  compute_md5(data.data(), data.size(), digest, digest_len);
  std::string CalculatedMD5Sum = to_hex(digest, digest_len);
  
  if (KnownSum != CalculatedMD5Sum) {
    std::cout << "MD5 Sums are different!" << std::endl;
    throw;
  } else {
    std::cout << "MD5 Sums are identical!" << std::endl;
  }

  return true; 
}

int SampleHandlerAtm::ReadFromEigen() {
  Eigen::MatrixXd Matrix;

  MACH3LOG_INFO("Reading Eigen matrix from:{}",EigenInputFile);
  CompareMD5Sum(EigenInputFileMD5Sum,EigenInputFile);
  ReadBinaryFromFile(EigenInputFile.c_str(),Matrix);

  int nEntries = static_cast<int>(Matrix.rows());
  dunemcSamples.resize(nEntries);

  for (size_t iEvent=0;iEvent<dunemcSamples.size();++iEvent) {
    dunemcSamples[iEvent].Target = static_cast<int>(Matrix(iEvent,Variables::Target));
    dunemcSamples[iEvent].nupdg = static_cast<int>(Matrix(iEvent,Variables::nupdg));
    dunemcSamples[iEvent].nupdgUnosc = static_cast<int>(Matrix(iEvent,Variables::nupdgUnosc));
    dunemcSamples[iEvent].rw_isCC = static_cast<int>(Matrix(iEvent,Variables::rw_isCC));
    dunemcSamples[iEvent].OscChannelIndex = Matrix(iEvent,Variables::OscChannelIndex);
    dunemcSamples[iEvent].rw_erec = Matrix(iEvent,Variables::rw_erec);
    dunemcSamples[iEvent].rw_etru = Matrix(iEvent,Variables::rw_etru);
    dunemcSamples[iEvent].flux_w = Matrix(iEvent,Variables::flux_w);
    dunemcSamples[iEvent].mode = Matrix(iEvent,Variables::mode);
    dunemcSamples[iEvent].rw_theta = Matrix(iEvent,Variables::rw_theta);
    dunemcSamples[iEvent].rw_truecz = Matrix(iEvent,Variables::rw_truecz);
  }

  return nEntries;
}

void SampleHandlerAtm::TransferToEigen(std::string FileName) {
  Eigen::MatrixXd Matrix = Eigen::MatrixXd(dunemcSamples.size(),Variables::nVariables);

  for (size_t iEvent=0;iEvent<dunemcSamples.size();++iEvent) {
    Matrix(iEvent,Variables::Target) = dunemcSamples[iEvent].Target;
    Matrix(iEvent,Variables::nupdg) = dunemcSamples[iEvent].nupdg;
    Matrix(iEvent,Variables::nupdgUnosc) = dunemcSamples[iEvent].nupdgUnosc;
    Matrix(iEvent,Variables::rw_isCC) = dunemcSamples[iEvent].rw_isCC;
    Matrix(iEvent,Variables::OscChannelIndex) = dunemcSamples[iEvent].OscChannelIndex;
    Matrix(iEvent,Variables::rw_erec) = dunemcSamples[iEvent].rw_erec;
    Matrix(iEvent,Variables::rw_etru) = dunemcSamples[iEvent].rw_etru;
    Matrix(iEvent,Variables::flux_w) = dunemcSamples[iEvent].flux_w;
    Matrix(iEvent,Variables::mode) = dunemcSamples[iEvent].mode;
    Matrix(iEvent,Variables::rw_theta) = dunemcSamples[iEvent].rw_theta;
    Matrix(iEvent,Variables::rw_truecz) = dunemcSamples[iEvent].rw_truecz;
  }

  MACH3LOG_INFO("Writing Eigen matrix to:{}",FileName);
  WriteBinaryToFile(FileName.c_str(),Matrix);
}

int SampleHandlerAtm::SetupExperimentMC() {
  if (EigenInputFile != "") {
    int nEntries = ReadFromEigen();
    return nEntries;
  }
  
  int CurrErrorLevel = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kFatal;
  
  caf::StandardRecord* sr = new caf::StandardRecord();

  TChain* Chain = new TChain("cafTree");
  for (size_t iSample=0;iSample<mc_files.size();iSample++) {
    Chain->Add(mc_files[iSample].c_str());
  }
  
  Chain->SetBranchStatus("*", 1);
  Chain->SetBranchAddress("rec", &sr);

  int nEntries = static_cast<int>(Chain->GetEntries());
  dunemcSamples.resize(nEntries);
 
  for (int iEvent=0;iEvent<nEntries;iEvent++) {
    Chain->GetEntry(iEvent);

    if ((iEvent % (nEntries/10))==0) {
      MACH3LOG_INFO("\tProcessing event: {}/{}",iEvent,nEntries);
    }
    
    if(sr->common.ixn.pandora.size() != 1) {
      MACH3LOG_WARN("Skipping event {}/{} -> Number of neutrino slices found in event: {}",iEvent,nEntries,sr->common.ixn.pandora.size());
      continue;
    }
    
    TVector3 RecoNuMomentumVector;
    double RecoENu;
    if (IsELike) {
      RecoENu = sr->common.ixn.pandora[0].Enu.e_calo;
      RecoNuMomentumVector = (TVector3(sr->common.ixn.pandora[0].dir.heshw.X(),sr->common.ixn.pandora[0].dir.heshw.Y(),sr->common.ixn.pandora[0].dir.heshw.Z())).Unit();
    } else {
      RecoENu = sr->common.ixn.pandora[0].Enu.lep_calo;
      RecoNuMomentumVector = (TVector3(sr->common.ixn.pandora[0].dir.lngtrk.X(),sr->common.ixn.pandora[0].dir.lngtrk.Y(),sr->common.ixn.pandora[0].dir.lngtrk.Z())).Unit();      
    }
    double RecoCZ = -RecoNuMomentumVector.Y(); // +Y in CAF files translates to +Z in typical CosZ
    if (std::isnan(RecoCZ)) {
      MACH3LOG_WARN("Skipping event {}/{} -> Reconstructed Cosine Z is NAN",iEvent,nEntries);
      continue;
    }
    if (std::isnan(RecoENu)) {
      MACH3LOG_WARN("Skipping event {}/{} -> Reconstructed Neutrino Energy is NAN",iEvent,nEntries);
      continue;
    }
    dunemcSamples[iEvent].rw_erec = RecoENu;
    dunemcSamples[iEvent].rw_theta = RecoCZ;
    
    std::string CurrFileName = Chain->GetCurrentFile()->GetName();
    dunemcSamples[iEvent].nupdgUnosc = GetInitPDGFromFileName(CurrFileName);
    dunemcSamples[iEvent].nupdg = GetFinalPDGFromFileName(CurrFileName);
    dunemcSamples[iEvent].OscChannelIndex = static_cast<double>(GetOscChannel(OscChannels, dunemcSamples[iEvent].nupdgUnosc, dunemcSamples[iEvent].nupdg));
    
    int M3Mode = Modes->GetModeFromGenerator(std::abs(sr->mc.nu[0].mode));
    if (!sr->mc.nu[0].iscc) M3Mode += 14; //Account for no ability to distinguish CC/NC
    if (M3Mode > 15) M3Mode -= 1; //Account for no NCSingleKaon
    dunemcSamples[iEvent].mode = M3Mode;
    
    dunemcSamples[iEvent].rw_isCC = sr->mc.nu[0].iscc;
    dunemcSamples[iEvent].Target = kTarget_Ar;
    
    dunemcSamples[iEvent].rw_etru = static_cast<double>(sr->mc.nu[0].E);

    TVector3 TrueNuMomentumVector = (TVector3(sr->mc.nu[0].momentum.X(),sr->mc.nu[0].momentum.Y(),sr->mc.nu[0].momentum.Z())).Unit();
    dunemcSamples[iEvent].rw_truecz = -TrueNuMomentumVector.Y(); // +Y in CAF files translates to +Z in typical CosZ

    dunemcSamples[iEvent].flux_w = sr->mc.nu[0].genweight;
  }

  delete Chain;
  gErrorIgnoreLevel = CurrErrorLevel;

  return nEntries;
}

void SampleHandlerAtm::SetupFDMC() {
  for(int iEvent = 0 ;iEvent < int(GetNEvents()) ; ++iEvent) {
    MCSamples[iEvent].rw_etru = &(dunemcSamples[iEvent].rw_etru);
    MCSamples[iEvent].mode = &(dunemcSamples[iEvent].mode);
    MCSamples[iEvent].Target = &(dunemcSamples[iEvent].Target);    
    MCSamples[iEvent].isNC = !dunemcSamples[iEvent].rw_isCC;
    MCSamples[iEvent].nupdg = &(dunemcSamples[iEvent].nupdg);
    MCSamples[iEvent].nupdgUnosc = &(dunemcSamples[iEvent].nupdgUnosc);

    MCSamples[iEvent].rw_truecz = &(dunemcSamples[iEvent].rw_truecz);
  }
}

const double* SampleHandlerAtm::GetPointerToKinematicParameter(KinematicTypes KinPar, int iEvent) {
  double* KinematicValue;

  switch (KinPar) {
  case kTrueNeutrinoEnergy:
    KinematicValue = &(dunemcSamples[iEvent].rw_etru);
    break;
  case kRecoNeutrinoEnergy:
    KinematicValue = &(dunemcSamples[iEvent].rw_erec);
    break;
  case kTrueCosZ:
    KinematicValue = &(dunemcSamples[iEvent].rw_truecz);
    break;
  case kRecoCosZ:
    KinematicValue = &(dunemcSamples[iEvent].rw_theta);
    break;
  case kOscChannel:
    KinematicValue = &(dunemcSamples[iEvent].OscChannelIndex);
    break;
  case kMode:
    KinematicValue = &(dunemcSamples[iEvent].mode);
    break;
  default:
    MACH3LOG_ERROR("Unknown KinPar: {}",static_cast<int>(KinPar));
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  
  return KinematicValue;
}

const double* SampleHandlerAtm::GetPointerToKinematicParameter(double KinematicVariable, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(KinematicVariable);
  return GetPointerToKinematicParameter(KinPar,iEvent);
}

const double* SampleHandlerAtm::GetPointerToKinematicParameter(std::string KinematicParameter, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(ReturnKinematicParameterFromString(KinematicParameter));
  return GetPointerToKinematicParameter(KinPar,iEvent);
}

double SampleHandlerAtm::ReturnKinematicParameter(int KinematicVariable, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(KinematicVariable);
  return *GetPointerToKinematicParameter(KinPar, iEvent);
}

double SampleHandlerAtm::ReturnKinematicParameter(std::string KinematicParameter, int iEvent) {
  return *GetPointerToKinematicParameter(KinematicParameter, iEvent);
}
