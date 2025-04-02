#include "splinesDUNE.h"
#include "TROOT.h"

splinesDUNE::splinesDUNE(covarianceXsec* xsec_cov, MaCh3Modes* Modes_) : splineFDBase(xsec_cov, Modes_) {
}

splinesDUNE::~splinesDUNE() {
}

std::vector<std::string> splinesDUNE::GetTokensFromSplineName(std::string FullSplineName) {
  std::vector<std::string> ReturnVec(TokenOrdering::kNTokens);

  TObjArray *tokens = TString(FullSplineName).Tokenize("_");
  
  ReturnVec[TokenOrdering::kSystToken] = ((TObjString *)(tokens->At(1)))->GetString();
  ReturnVec[TokenOrdering::kModeToken] = ((TObjString *)(tokens->At(2)))->GetString();
  // Skip 3 because it's "sp"
  ReturnVec[TokenOrdering::kVar1BinToken] = ((TObjString *)(tokens->At(4)))->GetString();
  ReturnVec[TokenOrdering::kVar2BinToken] = ((TObjString *)(tokens->At(5)))->GetString();
  ReturnVec[TokenOrdering::kVar3BinToken] = "0";
  
  if (tokens->GetEntries() == 7) {
    ReturnVec[TokenOrdering::kVar3BinToken] = ((TObjString *)(tokens->At(6)))->GetString();
  }

  return ReturnVec;
}

/*
///
DB: Warning to anyone looking at the below function and thinking why doesn't it exist in core?
    - It's because there are duplicate splines and special treatment for the unknown spline mode
      I have zero idea of why unknown spline modes even exist for one thing..
      It's such a mess - the DUNE splines need remaking to be renamed properly and duplicates removed
      Then the below can be removed
///
*/

void splinesDUNE::FillSampleArray(std::string SampleName, std::vector<std::string> OscChanFileNames) {
  int iSample = getSampleIndex(SampleName);
  int nOscChannels = nOscChans[iSample];
  
  for (int iOscChan = 0; iOscChan < nOscChannels; iOscChan++) {
    MACH3LOG_INFO("Processing: {}", OscChanFileNames[iOscChan]);
    TSpline3* mySpline = nullptr;
    TSpline3_red* Spline = nullptr;
    TString Syst, Mode;
    int nKnots, SystNum, ModeNum, Var1Bin, Var2Bin, Var3Bin = M3::_BAD_INT_;
    double x,y, Eval = M3::_BAD_DOUBLE_;
    bool isFlat = true;

    std::set<std::string> SplineFileNames;
    auto File = std::unique_ptr<TFile>(TFile::Open(OscChanFileNames[iOscChan].c_str()));

    if (!File || File->IsZombie()) {
      MACH3LOG_ERROR("File {} not found", OscChanFileNames[iOscChan]);
      throw MaCh3Exception(__FILE__, __LINE__);
    }

    //This is the MC specific part of the code
    //i.e. we always assume that the splines are just store in  single TDirectory and they're all in there as single objects   
    for (auto k : *File->GetListOfKeys()) {
      auto Key = static_cast<TKey*>(k);
      TClass *Class = gROOT->GetClass(Key->GetClassName(), false);
      if(!Class->InheritsFrom("TSpline3")) {
        continue;
      }

      std::string FullSplineName = std::string(Key->GetName());
      
      if (SplineFileNames.count(FullSplineName) > 0) {
	//if (FullSplineName.find("unknown") == std::string::npos) {
	  MACH3LOG_CRITICAL("Skipping spline - Found a spline whose name has already been encountered before: {}", FullSplineName); 
	  continue;
	  //}
      }
      SplineFileNames.insert(FullSplineName);

      std::vector<std::string> Tokens = GetTokensFromSplineName(FullSplineName);

      if (Tokens.size() != kNTokens) {
        MACH3LOG_ERROR("Invalid tokens from spline name - Expected {} tokens. Check implementation in GetTokensFromSplineName()", static_cast<int>(kNTokens));
        throw MaCh3Exception(__FILE__, __LINE__);
      }
      
      Syst = Tokens[kSystToken];
      Mode = Tokens[kModeToken];
      Var1Bin = std::stoi(Tokens[kVar1BinToken]);
      Var2Bin = std::stoi(Tokens[kVar2BinToken]);
      Var3Bin = std::stoi(Tokens[kVar3BinToken]);

      SystNum = -1;
      for (unsigned iSyst = 0; iSyst < SplineFileParPrefixNames[iSample].size(); iSyst++) {
        if (Syst == SplineFileParPrefixNames[iSample][iSyst]) {
          SystNum = iSyst;
          break;
        }
      }

      // If the syst doesn't match any of the spline names then skip it
      if (SystNum == -1){
        MACH3LOG_WARN("Couldn't match!!");
        MACH3LOG_DEBUG("Couldn't Match any systematic name in xsec yaml with spline name: {}" , FullSplineName);
        continue;
      }

      ModeNum = -1;
      for (unsigned int iMode = 0; iMode < SplineModeVecs[iSample][SystNum].size(); iMode++) {
        if (Mode == Modes->GetSplineSuffixFromMaCh3Mode(SplineModeVecs[iSample][SystNum][iMode])) {
          ModeNum = iMode;
          break;
        }
      }

      if (ModeNum == -1) {
	//DB - If you have splines in the root file that you don't want to use (e.g. removing a mode from a syst), this will cause a throw
	//     Therefore include as debug warning and continue instead
        MACH3LOG_DEBUG("Couldn't find mode for {} in {}. Problem Spline is : {} ", Mode, Syst, FullSplineName);
	continue;
      }

      mySpline = Key->ReadObject<TSpline3>();

      if (isValidSplineIndex(SampleName, iOscChan, SystNum, ModeNum, Var1Bin, Var2Bin, Var3Bin)) { // loop over all the spline knots and check their value
        MACH3LOG_DEBUG("Pushed back monolith for spline {}", FullSplineName);
        // if the value is 1 then set the flat bool to false
        nKnots = mySpline->GetNp();
        isFlat = true;
        for (int iKnot = 0; iKnot < nKnots; iKnot++) {
          mySpline->GetKnot(iKnot, x, y);

          Eval = mySpline->Eval(x);
          if (Eval < 0.99999 || Eval > 1.00001)
          {
            isFlat = false;
            break;
          }
        }

        //Rather than keeping a mega vector of splines then converting, this should just keep everything nice in memory!
        indexvec[iSample][iOscChan][SystNum][ModeNum][Var1Bin][Var2Bin][Var3Bin]=MonolithIndex;
        coeffindexvec.push_back(CoeffIndex);
        // Should save memory rather saving [x_i_0 ,... x_i_maxknots] for every spline!
        if (isFlat) {
          splinevec_Monolith.push_back(nullptr);
          delete mySpline;
        } else {
          Spline = new TSpline3_red(mySpline, SplineInterpolationTypes[iSample][SystNum]);
          delete mySpline;

          splinevec_Monolith.push_back(Spline);
          uniquecoeffindices.push_back(MonolithIndex); //So we can get the unique coefficients and skip flat splines later on!
          CoeffIndex+=nKnots;
        }
        //Incrementing MonolithIndex to keep track of number of valid spline indices
        MonolithIndex+=1;
      } else {
        //Potentially you are not a valid spline index
        delete mySpline;
      }
    }//End of loop over all TKeys in file

    //A bit of clean up
    File->Delete("*");
    File->Close();
  } //End of oscillation channel loop
}
