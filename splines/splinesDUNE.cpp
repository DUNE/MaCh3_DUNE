#include "splinesDUNE.h"
#include "samplePDFDUNE/StructsDUNE.h"

#include "TKey.h"
#include "TROOT.h"

splinesDUNE::splinesDUNE(covarianceXsec* xsec_cov) : splineFDBase(xsec_cov) {
  MACH3LOG_INFO("Created splinesDUNE object");
}

splinesDUNE::~splinesDUNE(){
  MACH3LOG_INFO("Deleting splineSKBase object");
}

//****************************************
void splinesDUNE::FillSampleArray(std::string SampleName, std::vector<std::string> OscChanFileNames)
// Performs two jobs
//  1. Fills indexing/each sample
//  2. Creates the big spline vector
//****************************************
{

  int iSample = getSampleIndex(SampleName);

  enum ModeState { kInConfig = 0, kInConfigAndInFile, kInFile };

  std::vector<std::map<std::string, ModeState>> ModeStatus;

  // Create map of Mode Status
  for (unsigned iSyst = 0; iSyst < SplineFileParPrefixNames[iSample].size(); iSyst++)
  {
    auto modes = SplineModeVecs[iSample][iSyst];
    ModeStatus.emplace_back();

    for (auto const &mode : modes)
    {
      // Add Modes from config
      ModeStatus.back()[MaCh3mode_ToDUNEString((MaCh3_Mode)mode).c_str()] = kInConfig;
    }
  }

  int nOscChannels = nOscChans[iSample];

  for (int iOscChan = 0; iOscChan < nOscChannels; iOscChan++)
  {
    std::cout << "Processing:" << OscChanFileNames[iOscChan] << std::endl;

    TFile *File = new TFile(OscChanFileNames[iOscChan].c_str(), "READ");
    if (!File || File->IsZombie())
    {
      std::cerr << "File " << OscChanFileNames[iOscChan] << " not found" << std::endl;
      throw;
    }

    // This is the MC specific part of the code
    // i.e. we always assume that the splines are just store in  single TDirectory and they're all
    // in there as single objects
    TIter Next(File->GetListOfKeys());
    TKey *Key;

    std::set<std::string> unique_spline_names;
    int nb_splines = 0;

    while ((Key = (TKey *)Next()))
    {
      TClass *Class = gROOT->GetClass(Key->GetClassName());

      // Skip the TGraphs also in the spline files
      if (!Class->InheritsFrom("TSpline3"))
      {
        continue;
      }
      const char *keyName = Key->GetName();

      char *SplineName = new char[strlen(keyName) + 1];
      strcpy(SplineName, keyName);

      nb_splines += 1;
      if (unique_spline_names.count(std::string(SplineName)) > 0)
      {
        if (std::string(SplineName).find("unknown") == std::string::npos)
        {
          // std::cout << "Repeated entry for spline named: " << std::string(SplineName) <<
          // std::endl;
          continue;
        }
      }
      unique_spline_names.insert(std::string(SplineName));

      char *Syst;
      char *Mode;
      int Var1Bin;
      int Var2Bin;
      int Var3Bin;

      char *Token = strtok(SplineName, "_");
      Token = strtok(NULL, "_");
      Syst = Token;

      int SystNum = -1;
      for (unsigned iSyst = 0; iSyst < SplineFileParPrefixNames[iSample].size(); iSyst++)
      {
        if (strcmp(Syst, SplineFileParPrefixNames[iSample][iSyst].c_str()) == 0)
        {
          SystNum = iSyst;
          break;
        }
      }

      // If the syst doesn't match any of the spline names then skip it
      if (SystNum == -1)
      {
        continue;
      }

      int ModeNum = -1;
      Mode = strtok(NULL, "_");
      for (unsigned int iMode = 0; iMode < SplineModeVecs[iSample][SystNum].size(); iMode++)
      {
        if (strcmp(Mode, MaCh3mode_ToDUNEString((MaCh3_Mode)SplineModeVecs[iSample][SystNum][iMode])
                             .c_str()) == 0)
        {
          ModeNum = iMode;
          break;
        }
      }

      // Check if mode has been registered already
      if (ModeStatus[SystNum].count(Mode))
      {
        // Chech if mode has been found in config
        if (ModeStatus[SystNum][Mode] != kInFile)
        {
          ModeStatus[SystNum][Mode] = kInConfigAndInFile;
        } else
        {
          continue; // Skip if mode has been found in file but not config
        }
      } else
      {
        ModeStatus[SystNum][Mode] =
            kInFile; // If mode hasn't been registered then skip because it's not in the config
        continue;
      }

      TSpline3 *Obj = (TSpline3 *)Key->ReadObj();
      TSpline3_red *Spline = new TSpline3_red(Obj);
      delete Obj;

      Token = strtok(NULL, "_"); // DB Needed to remove sp from spline name

      Var1Bin = atoi(strtok(NULL, "_"));
      Var2Bin = atoi(strtok(NULL, "_"));

      char *Var3Bin_Char = strtok(NULL, "_");
      if (Var3Bin_Char == NULL)
      {
        Var3Bin = 0;
      } else
      {
        Var3Bin = atoi(Var3Bin_Char);
      }

      if (isValidSplineIndex(SampleName, iOscChan, SystNum, ModeNum, Var1Bin, Var2Bin, Var3Bin))
      {
        // loop over all the spline knots and check their value
        // if the value is 1 then set the flat bool to false
        int nKnots = Spline->GetNp();
        bool isFlat = true;

        for (int iKnot = 0; iKnot < nKnots; iKnot++)
        {
          double x = -999;
          double y = -999;
          Spline->GetKnot(iKnot, x, y);

          if (x == -999 || y == -999)
          {
            std::cerr << "Something has gone wrong... knot position is at -999" << std::endl;
            std::cerr << "This error brought you by the folks at : " << __FILE__ << " : "
                      << __LINE__ << std::endl;
            throw;
          }

          double Eval = Spline->Eval(x);
          if (Eval < 0.99999 || Eval > 1.00001)
          {
            isFlat = false;
            break;
          }
        }

        // Rather than keeping a mega vector of splines then converting, this should just keep
        // everything nice in memory!
        indexvec[iSample][iOscChan][SystNum][ModeNum][Var1Bin][Var2Bin][Var3Bin] = MonolithIndex;

        coeffindexvec.push_back(CoeffIndex);
        // Should save memory rather saving [x_i_0 ,... x_i_maxknots] for every spline!
        if (isFlat)
        {
          splinevec_Monolith.push_back(NULL);
          delete Spline;
        } else
        {
          splinevec_Monolith.push_back(Spline);
          int np = Spline->GetNp();
          uniquecoeffindices.push_back(MonolithIndex); // So we can get the unique coefficients and
                                                       // skip flat splines later on!
          CoeffIndex += np;
        }

        MonolithIndex += 1;
      }
    } // End of loop over all TKeys in file
    // ETA - I have no idea why but this breaks in ROOT 6.24 :/
    std::cout << "Got " << nb_splines << " total splines with " << unique_spline_names.size()
              << " unique names." << std::endl;
    delete File;
  } // End of oscillation channel loop

  // Find all modes which have been found in the spline file but have not been specified in the
  // config
  for (unsigned iSyst = 0; iSyst < SplineFileParPrefixNames[iSample].size(); iSyst++)
  {
    std::vector<std::string> MissedModes;
    for (auto const &ModeUsage : ModeStatus[iSyst])
    {
      if (ModeUsage.second == kInFile)
      {
        MissedModes.push_back(ModeUsage.first);
      }
    }

    if (MissedModes.size() != 0)
    {
      MACH3LOG_INFO("Parameter {} has splines for {} modes which have not been read in!",
                    SplineFileParPrefixNames[iSample][iSyst].c_str(), MissedModes.size());
      std::cout << "Modes: ";
      for (auto const &Mode : MissedModes)
      {
        std::cout << Mode << " ";
      }
      std::cout << std::endl;
    }
  }

  return;
}

//****************************************
std::vector< std::vector<int> > splinesDUNE::GetEventSplines(std::string SampleName, int iOscChan, int EventMode, double Var1Val, double Var2Val, double Var3Val)
//****************************************
{
  std::vector<std::vector<int>> ReturnVec;
  int SampleIndex = -1;
  for (unsigned int iSample = 0; iSample < SampleNames.size(); iSample++) {
    if (SampleName == SampleNames[iSample]) {
      SampleIndex = iSample;
    }
  }

  if (SampleIndex == -1)
  {
    MACH3LOG_ERROR("Sample not found: {}", SampleName);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  int nSplineSysts = (int)indexvec[SampleIndex][iOscChan].size();
  //ETA- this is already a MaCh3 mode
  //int Mode = MaCh3Mode_to_SplineMode(MaCh3_Mode(EventMode));
  int Mode = MaCh3_Mode(EventMode);

  for (int iSyst = 0; iSyst < nSplineSysts; iSyst++) {
    std::vector<int> spline_modes = SplineModeVecs[SampleIndex][iSyst];
    int nSampleModes = (int)spline_modes.size();

    // ETA - look here at the length of spline_modes and what you're actually
    // comparing against
    for (int iMode = 0; iMode < nSampleModes; iMode++) {
      if (Mode == spline_modes[iMode]) {
        // HH: Moved this chunk to be after you do the mode checking
        int Var1Bin =
            SplineBinning[SampleIndex][iOscChan][0]->FindBin(Var1Val) - 1;
        if (Var1Bin < 0 ||
            Var1Bin >= SplineBinning[SampleIndex][iOscChan][0]->GetNbins()) {
          // HH: Removing push back and instead change Var1Bin to the last bin
          // std::cout << "======================================================"
          //           << std::endl;
          // std::cout << "I'm moving event with Var1Val: " << Var1Val
          //           << " to the last Var1Bin!" << std::endl;
          // std::cout << "Sample name: " << SampleName << std::endl;
          // std::cout << "Oscillation channel: " << iOscChan << std::endl;
          // std::cout << "MaCh3 mode: " << Mode << std::endl;
          // std::cout << "Var1Val: " << Var1Val << std::endl;
          // std::cout << "Var2Val: " << Var2Val << std::endl;
          // std::cout << "Var3Val: " << Var3Val << std::endl;
          // std::cout << "======================================================"
          //           << std::endl;
          Var1Bin = SplineBinning[SampleIndex][iOscChan][0]->GetNbins() - 1;
          // Explicitly push back with an empty vector
          // ReturnVec.push_back(std::vector<int>());
          // return ReturnVec;
        }

        int Var2Bin =
            SplineBinning[SampleIndex][iOscChan][1]->FindBin(Var2Val) - 1;
        if (Var2Bin < 0 ||
            Var2Bin >= SplineBinning[SampleIndex][iOscChan][1]->GetNbins()) {
          // HH: Removing push back and instead change Var2Bin to last bin
          // std::cout << "======================================================"
          //           << std::endl;
          // std::cout << "I'm moving event with Var2Val: " << Var2Val
          //           << " to the last Var2Bin!" << std::endl;
          // std::cout << "Sample name: " << SampleName << std::endl;
          // std::cout << "Oscillation channel: " << iOscChan << std::endl;
          // std::cout << "MaCh3 mode: " << Mode << std::endl;
          // std::cout << "Var1Val: " << Var1Val << std::endl;
          // std::cout << "Var2Val: " << Var2Val << std::endl;
          // std::cout << "Var3Val: " << Var3Val << std::endl;
          // std::cout << "======================================================"
          //           << std::endl;
          Var2Bin = SplineBinning[SampleIndex][iOscChan][1]->GetNbins() - 1;
          // Explicitly push back with an empty vector
          // ReturnVec.push_back(std::vector<int>());
          // return ReturnVec;
        }

        int Var3Bin =
            SplineBinning[SampleIndex][iOscChan][2]->FindBin(Var3Val) - 1;

        if (Var3Bin < 0 ||
            Var3Bin >= SplineBinning[SampleIndex][iOscChan][2]->GetNbins()) {
          // HH: Removing push back and instead change Var3Bin to last bin
          // std::cout << "======================================================"
          //           << std::endl;
          // std::cout << "I'm moving event with Var3Val: " << Var3Val
          //           << " to the last Var3Bin!" << std::endl;
          // std::cout << "Sample name: " << SampleName << std::endl;
          // std::cout << "Oscillation channel: " << iOscChan << std::endl;
          // std::cout << "MaCh3 mode: " << Mode << std::endl;
          // std::cout << "Var1Val: " << Var1Val << std::endl;
          // std::cout << "Var2Val: " << Var2Val << std::endl;
          // std::cout << "Var3Val: " << Var3Val << std::endl;
          // std::cout << "======================================================"
          //           << std::endl;
          Var3Bin = SplineBinning[SampleIndex][iOscChan][2]->GetNbins() - 1;
          // Explicitly push back with an empty vector
          // ReturnVec.push_back(std::vector<int>());
          // return ReturnVec;
        }

        std::vector<int> event_vec(7);
        event_vec[0] = SampleIndex;
        event_vec[1] = iOscChan;
        event_vec[2] = iSyst;
        event_vec[3] = iMode;
        event_vec[4] = Var1Bin;
        event_vec[5] = Var2Bin;
        event_vec[6] = Var3Bin;
        int splineID = indexvec[SampleIndex][iOscChan][iSyst][iMode][Var1Bin]
                               [Var2Bin][Var3Bin];
        if (!isflatarray[splineID]) {
          ReturnVec.push_back(event_vec);
        }
      }
    }
  }

  return ReturnVec;
}

// Converts MaCh3 modes to unique modes that splines apply to
// (for example if CCRES and CCCoherent are treated as one spline mode)
std::vector< std::vector<int> > splinesDUNE::StripDuplicatedModes(std::vector< std::vector<int> > InputVector) {
  
  //ETA - this is of size nPars from the xsec model 
  int InputVectorSize = InputVector.size();
  std::vector< std::vector<int> > ReturnVec(InputVectorSize);

  //ETA - loop over all systematics
  for (int iVec=0;iVec<InputVectorSize;iVec++) {
    std::vector<int> TmpVec;
	
	//Loop over the modes that we've listed in xsec cov
	for (unsigned int iMode = 0 ; iMode < InputVector[iVec].size() ; iMode++){
	  //Convert the MaCh3 mode to spline mode (could use MaCh3Mode_SplineMode_Map here)
	  int iSplineMode = MaCh3Mode_to_SplineMode(InputVector[iVec][iMode]);
	  bool IncludeMode = true;

	  //Now check to see if we've already included this unique spline mode
	  for(unsigned int entry_i = 0 ; entry_i < TmpVec.size() ; entry_i++){
	    //Check to see if we've already included this unique SK spline mode
	    if(iSplineMode == TmpVec[entry_i]){
	      IncludeMode = false;
	    }
	  }
	  
	  if(IncludeMode){
	    //Push back with the spline mode a systematic applies to
	    TmpVec.push_back(iSplineMode);
	  }
	  
	}//end of loop over modes for a syst
	
	ReturnVec[iVec] = TmpVec;
  }//end of loop over systs
  
  return ReturnVec;
}
