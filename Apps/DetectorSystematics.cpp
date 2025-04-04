#include <vector>

#include <TH1D.h>
#include <TH2D.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TLegend.h>

#include "samplePDFDUNE/MaCh3DUNEFactory.h"
#include "samplePDFDUNE/StructsDUNE.h"

struct KinematicCut {
  std::string Name;
  std::string VarString;
  std::vector<double> Range;
};

struct CategoryCut {
  std::string Name;
  std::string VarString;
  std::vector< std::vector<double> > Breakdown;
  std::vector<double> Colours;
  std::vector<std::string> CategoryNames;
};

struct DetectorSystematics1DVariable {
  std::string Name;
  std::string VarString;
  std::vector<double> BinEdges;
  
  std::vector<KinematicCut> KinematicCuts;
  std::vector<CategoryCut> CategoryCuts;
};

struct DetectorSystematics2DVariable {
  std::string XName;
  std::string YName;
  std::string XVarString;
  std::string YVarString;
  std::vector<double> XBinEdges;
  std::vector<double> YBinEdges;
  
  std::vector<KinematicCut> KinematicCuts;
  std::vector<CategoryCut> CategoryCuts;
};

std::string ReturnFormattedHistogramNameFromDetectorSystematic(DetectorSystematics1DVariable Proj) {
  std::string ReturnStr;

  for (size_t iKinematicCut=0;iKinematicCut<Proj.KinematicCuts.size();iKinematicCut++) {
    if (iKinematicCut > 0) {
      ReturnStr += " && ";
    }
    ReturnStr += Form("(%4.2f < %s < %4.2f)",Proj.KinematicCuts[iKinematicCut].Range[0],Proj.KinematicCuts[iKinematicCut].Name.c_str(),Proj.KinematicCuts[iKinematicCut].Range[1]);
  }

  ReturnStr += ";"+Proj.Name+";Events";
  return ReturnStr;
}

void PrintTH1Histogram(TH1* Hist, std::string OutputName) {
  TCanvas Canv = TCanvas();
  // Hist->SetStats(0);
  Hist->Draw("HIST");
  Canv.Print(OutputName.c_str());
}

void PrintTH2Histogram(TH2* Hist, std::string OutputName) {
    TCanvas Canv = TCanvas();
    Hist->Draw("COLZ");
    Canv.Print(OutputName.c_str());
  }

void PrintCategoryLegends(std::vector<DetectorSystematics1DVariable> DetectorSystematics) {
  TLegend Legend = TLegend(0.1,0.1,0.9,0.9);

  TCanvas Canv = TCanvas();

  std::vector<TH1D*> HistVec;
  
  for (size_t iProj=0;iProj<DetectorSystematics.size();iProj++) {
    for (size_t iCat=0;iCat<DetectorSystematics[iProj].CategoryCuts.size();iCat++) {
      CategoryCut Cat = DetectorSystematics[iProj].CategoryCuts[iCat];

      Legend.SetHeader(DetectorSystematics[iProj].CategoryCuts[iCat].Name.c_str());

      HistVec.resize(DetectorSystematics[iProj].CategoryCuts[iCat].Breakdown.size());
      for (size_t iBreak=0;iBreak<Cat.Breakdown.size();iBreak++) {
	HistVec[iBreak] = new TH1D(Form("DummyHist_%i",(int)iBreak),"",1,0,1);
      }
      
      for (size_t iBreak=0;iBreak<Cat.Breakdown.size();iBreak++) {
	HistVec[iBreak]->SetFillColor(Cat.Colours[iBreak]);
	Legend.AddEntry(HistVec[iBreak],Cat.CategoryNames[iBreak].c_str(),"f");
      }

      Legend.Draw();
      Canv.Print(("Legend_"+DetectorSystematics[iProj].CategoryCuts[iCat].Name+".png").c_str());
      Legend.Clear();

      for (size_t iBreak=0;iBreak<Cat.Breakdown.size();iBreak++) {
	delete HistVec[iBreak];
      }
    }
  }

}

void PrintTHStackHistogram(THStack* Hist, std::string OutputName) {
  TCanvas Canv = TCanvas();
  Hist->Draw("HIST");
  Canv.Print(OutputName.c_str());
}

int main(int argc, char * argv[]) {
    if(argc == 1){
      MACH3LOG_ERROR("Usage: bin/EventRatesDUNEBeam config.cfg");
      return 1;
    }
    auto fitMan = std::unique_ptr<manager>(new manager(argv[1]));
  
    int WeightStyle = 1;
    
    //###############################################################################################################################
    //Create samplePDFFD objects
    
    covarianceXsec* xsec = nullptr;
    covarianceOsc* osc = nullptr;
    
    std::vector<samplePDFFDBase*> DUNEPdfs;
    MakeMaCh3DuneInstance(fitMan.get(), DUNEPdfs, xsec, osc);
    
    //###############################################################################################################################
    //Perform reweight and print total integral for sanity check
  
    MACH3LOG_INFO("=================================================");
    std::vector<TH1D*> DUNEHists;
    for(auto Sample : DUNEPdfs){
      Sample->reweight();
      DUNEHists.push_back(Sample->get1DHist());
      
      std::string EventRateString = fmt::format("{:.2f}", Sample->get1DHist()->Integral());
      MACH3LOG_INFO("Event rate for {} : {:<5}", Sample->GetName(), EventRateString);
    }
  
    //###############################################################################################################################
    //Grab DetectorSystematics from the config
  
    std::vector<DetectorSystematics1DVariable> DetectorSystematics1D;
    std::vector<DetectorSystematics2DVariable> DetectorSystematics2D;
    
    for (auto &DetectorSystematicConfig: fitMan->raw()["DetectorSystematics"]) {
      if (DetectorSystematicConfig["Type"].as<std::string>() == "1D") {
        std::string VarName = DetectorSystematicConfig["Name"].as<std::string>();
        std::string VarString = DetectorSystematicConfig["VarString"].as<std::string>();
  
        std::vector<double> VarBinEdges = DetectorSystematicConfig["VarBins"].as< std::vector<double> >();
  
        std::vector<KinematicCut> KinematicCuts;
        std::vector<CategoryCut> CategoryCuts;
        
        for (auto &KinematicCutConfig: DetectorSystematicConfig["KinematicCuts"]) {
          std::string KinematicCutName = KinematicCutConfig["Name"].as<std::string>();
          std::string KinematicCutVarString = KinematicCutConfig["VarString"].as<std::string>();
          std::vector<double> KinematicCutRange = KinematicCutConfig["Range"].as< std::vector<double> >();
          
          KinematicCut Cut = KinematicCut{KinematicCutName,KinematicCutVarString,KinematicCutRange};
          KinematicCuts.emplace_back(Cut);
        }
  
        for (auto &CategoryCutConfig: DetectorSystematicConfig["CategoryCuts"]) {
          std::string CategoryCutName = CategoryCutConfig["Name"].as<std::string>();
          std::string CategoryCutVarString = CategoryCutConfig["VarString"].as<std::string>();
          std::vector< std::vector<double> > CategoryCutBreakdown = CategoryCutConfig["Breakdown"].as< std::vector< std::vector<double> > >();
  
          std::vector<double> CategoryCutColours;
          if (CategoryCutConfig["Colours"]) {
            CategoryCutColours = CategoryCutConfig["Colours"].as< std::vector<double> >();
          } else {
            CategoryCutColours.resize(CategoryCutBreakdown.size());
            int colour = 20.;
            for (size_t iColour=0;iColour<CategoryCutColours.size();iColour++) {
              CategoryCutColours[iColour] = colour;
              colour += 4;
              if (colour > 50) {colour -= 30;}
            }
          }
  
          std::vector<std::string> CategoryCutNames;
          if (CategoryCutConfig["Names"]) {
            CategoryCutNames = CategoryCutConfig["Names"].as< std::vector<std::string> >();
          } else {
            CategoryCutNames.resize(CategoryCutBreakdown.size());
            for (size_t i=0;i<CategoryCutBreakdown.size();i++) {
              CategoryCutNames[i] = Form("%i",(int)i);
            }
          }
  
          CategoryCut Cut = CategoryCut{CategoryCutName,CategoryCutVarString,CategoryCutBreakdown,CategoryCutColours,CategoryCutNames};
          CategoryCuts.emplace_back(Cut);
        }
        
        DetectorSystematics1DVariable Proj = DetectorSystematics1DVariable{VarName,VarString,VarBinEdges,KinematicCuts,CategoryCuts};
        DetectorSystematics1D.emplace_back(Proj);
      } else if (DetectorSystematicConfig["Type"].as<std::string>() == "2D") {
        std::string XName = DetectorSystematicConfig["XName"].as<std::string>();
        std::string YName = DetectorSystematicConfig["YName"].as<std::string>();
        std::string XVarString = DetectorSystematicConfig["XVarString"].as<std::string>();
        std::string YVarString = DetectorSystematicConfig["YVarString"].as<std::string>();
  
        std::vector<double> XBinEdges = DetectorSystematicConfig["XVarBins"].as< std::vector<double> >();
        std::vector<double> YBinEdges = DetectorSystematicConfig["YVarBins"].as< std::vector<double> >();
  
        std::vector<KinematicCut> KinematicCuts;
        std::vector<CategoryCut> CategoryCuts;
        
        for (auto &KinematicCutConfig: DetectorSystematicConfig["KinematicCuts"]) {
          std::string KinematicCutName = KinematicCutConfig["Name"].as<std::string>();
          std::string KinematicCutVarString = KinematicCutConfig["VarString"].as<std::string>();
          std::vector<double> KinematicCutRange = KinematicCutConfig["Range"].as< std::vector<double> >();
          
          KinematicCut Cut = KinematicCut{KinematicCutName,KinematicCutVarString,KinematicCutRange};
          KinematicCuts.emplace_back(Cut);
        }
  
        for (auto &CategoryCutConfig: DetectorSystematicConfig["CategoryCuts"]) {
          std::string CategoryCutName = CategoryCutConfig["Name"].as<std::string>();
          std::string CategoryCutVarString = CategoryCutConfig["VarString"].as<std::string>();
          std::vector< std::vector<double> > CategoryCutBreakdown = CategoryCutConfig["Breakdown"].as< std::vector< std::vector<double> > >();
  
          std::vector<double> CategoryCutColours;
          if (CategoryCutConfig["Colours"]) {
            CategoryCutColours = CategoryCutConfig["Colours"].as< std::vector<double> >();
          } else {
            CategoryCutColours.resize(CategoryCutBreakdown.size());
            int colour = 20.;
            for (size_t iColour=0;iColour<CategoryCutColours.size();iColour++) {
              CategoryCutColours[iColour] = colour;
              colour += 4;
              if (colour > 50) {colour -= 30;}
            }
          }
  
          std::vector<std::string> CategoryCutNames;
          if (CategoryCutConfig["Names"]) {
            CategoryCutNames = CategoryCutConfig["Names"].as< std::vector<std::string> >();
          } else {
            CategoryCutNames.resize(CategoryCutBreakdown.size());
            for (size_t i=0;i<CategoryCutBreakdown.size();i++) {
              CategoryCutNames[i] = Form("%i",(int)i);
            }
          }
  
          CategoryCut Cut = CategoryCut{CategoryCutName,CategoryCutVarString,CategoryCutBreakdown,CategoryCutColours,CategoryCutNames};
          CategoryCuts.emplace_back(Cut);
        }
        
        DetectorSystematics2DVariable Proj = DetectorSystematics2DVariable{XName, YName, XVarString, YVarString, XBinEdges, YBinEdges, KinematicCuts, CategoryCuts};
        DetectorSystematics2D.emplace_back(Proj);
      }
    }
  
    MACH3LOG_INFO("=================================================");
    MACH3LOG_INFO("DetectorSystematics pulled from Config..");
    MACH3LOG_INFO("================================");
    
    for (size_t iProj=0;iProj<DetectorSystematics1D.size();iProj++) {
      MACH3LOG_INFO("DetectorSystematic {:<2} - Name : {:<20} , VarString : {:<20}",iProj,DetectorSystematics1D[iProj].Name,DetectorSystematics1D[iProj].VarString);
      MACH3LOG_INFO("\t\tBinning: {}", fmt::join(DetectorSystematics1D[iProj].BinEdges, ", "));
  
      if (DetectorSystematics1D[iProj].KinematicCuts.size()>0) {
        MACH3LOG_INFO("\t\tKinematicCuts:");
        for (size_t iCut=0;iCut<DetectorSystematics1D[iProj].KinematicCuts.size();iCut++) {
          MACH3LOG_INFO("\t\t\tCut {:<2} - Name : {:<20} , Lower Bound : {:<10} , Upper Bound : {:<10}",iCut,DetectorSystematics1D[iProj].KinematicCuts[iCut].Name,DetectorSystematics1D[iProj].KinematicCuts[iCut].Range[0],DetectorSystematics1D[iProj].KinematicCuts[iCut].Range[1]);
        }
      }
  
      if (DetectorSystematics1D[iProj].CategoryCuts.size()>0) {
        MACH3LOG_INFO("\t\tCategoryCuts:");
        for (size_t iCut=0;iCut<DetectorSystematics1D[iProj].CategoryCuts.size();iCut++) {
  
          std::vector<std::string> BreakdownStrs(DetectorSystematics1D[iProj].CategoryCuts[iCut].Breakdown.size());
          for (size_t iBreak=0;iBreak<DetectorSystematics1D[iProj].CategoryCuts[iCut].Breakdown.size();iBreak++) {
            BreakdownStrs[iBreak] = fmt::format("{}",fmt::join(DetectorSystematics1D[iProj].CategoryCuts[iCut].Breakdown[iBreak], ", "));
          }
          MACH3LOG_INFO("\t\t\tCategory {:<2} - Name : {:<20} , Category Breakdown : {}",iCut,DetectorSystematics1D[iProj].CategoryCuts[iCut].Name,fmt::join(BreakdownStrs, ", "));
          
        }
      }
      MACH3LOG_INFO("================================");
    }
  
    PrintCategoryLegends(DetectorSystematics1D);
  
    //###############################################################################################################################
  
    MACH3LOG_INFO("=================================================");
    MACH3LOG_INFO("Building DetectorSystematics..");
  
    TH1* Hist;
    THStack* Stack;
    TH2D* Hist2D;
    
    for (size_t iProj=0;iProj<DetectorSystematics1D.size();iProj++) {
      MACH3LOG_INFO("================================");
      MACH3LOG_INFO("DetectorSystematic {}/{}",iProj,DetectorSystematics1D.size());
    
      std::string DetectorSystematicVar_Str = DetectorSystematics1D[iProj].VarString;
      TAxis Axis = TAxis(DetectorSystematics1D[iProj].BinEdges.size()-1,DetectorSystematics1D[iProj].BinEdges.data());
  
      for (auto Sample: DUNEPdfs) {
  
        std::vector< std::vector<double> > SelectionVector;
        for (size_t iCut=0;iCut<DetectorSystematics1D[iProj].KinematicCuts.size();iCut++) {
          std::vector<double> Selection(3);
          Selection[0] = Sample->ReturnKinematicParameterFromString(DetectorSystematics1D[iProj].KinematicCuts[iCut].VarString);
          Selection[1] = DetectorSystematics1D[iProj].KinematicCuts[iCut].Range[0];
          Selection[2] = DetectorSystematics1D[iProj].KinematicCuts[iCut].Range[1];
          
          SelectionVector.emplace_back(Selection);
        }
        Hist = Sample->get1DVarHist(DetectorSystematicVar_Str,SelectionVector,WeightStyle,&Axis);
        Hist->Scale(1.0,"Width");
        Hist->SetTitle(ReturnFormattedHistogramNameFromDetectorSystematic(DetectorSystematics1D[iProj]).c_str());
  
        MACH3LOG_INFO("\tSample: {:<20} - Integral: {:<10}",Sample->GetName(),Hist->Integral());
        PrintTH1Histogram(Hist,Sample->GetName()+"_"+DetectorSystematicVar_Str+".png");
        
        for (size_t iCat=0;iCat<DetectorSystematics1D[iProj].CategoryCuts.size();iCat++) {
          MACH3LOG_INFO("\t\tCategory: {:<10} - Name : {:<20}",iCat,DetectorSystematics1D[iProj].CategoryCuts[iCat].Name);
  
          Stack = new THStack(DetectorSystematics1D[iProj].CategoryCuts[iCat].Name.c_str(),ReturnFormattedHistogramNameFromDetectorSystematic(DetectorSystematics1D[iProj]).c_str());
  
          for (size_t iBreak=0;iBreak<DetectorSystematics1D[iProj].CategoryCuts[iCat].Breakdown.size();iBreak++) {
  
            TH1* BreakdownHist = nullptr;
  
            for (size_t iGroup=0;iGroup<DetectorSystematics1D[iProj].CategoryCuts[iCat].Breakdown[iBreak].size();iGroup++) {
              std::vector< std::vector<double> > SelectionVector_IncCategory = std::vector< std::vector<double> >(SelectionVector);
              
              std::vector<double> Selection(3);
              Selection[0] = Sample->ReturnKinematicParameterFromString(DetectorSystematics1D[iProj].CategoryCuts[iCat].VarString);
              Selection[1] = DetectorSystematics1D[iProj].CategoryCuts[iCat].Breakdown[iBreak][iGroup];
              Selection[2] = DetectorSystematics1D[iProj].CategoryCuts[iCat].Breakdown[iBreak][iGroup]+1;
              SelectionVector_IncCategory.emplace_back(Selection);
              
              Hist = Sample->get1DVarHist(DetectorSystematicVar_Str,SelectionVector_IncCategory,WeightStyle,&Axis);
              Hist->SetFillColor(DetectorSystematics1D[iProj].CategoryCuts[iCat].Colours[iBreak]);
              Hist->Scale(1.0,"Width");
  
              if (BreakdownHist == nullptr) {
                BreakdownHist = Hist;
              } else {
                BreakdownHist->Add(Hist);
              }
            }
  
            MACH3LOG_INFO("\t\t\tBreakdown: {:<10} - Integral: {:<10}",iBreak,BreakdownHist->Integral());
            Stack->Add(BreakdownHist);
          }
  
          PrintTHStackHistogram(Stack,Sample->GetName()+"_"+DetectorSystematicVar_Str+"_"+DetectorSystematics1D[iProj].CategoryCuts[iCat].Name+"_Stack.png");
        }
        
      }
    }
    
  // 2D Histograms
  for (size_t iProj = 0; iProj < DetectorSystematics2D.size(); iProj++) {
    MACH3LOG_INFO("================================");
    MACH3LOG_INFO("DetectorSystematic {}/{}", iProj, DetectorSystematics2D.size());
  
    std::string XVarString = DetectorSystematics2D[iProj].XVarString;
    std::string YVarString = DetectorSystematics2D[iProj].YVarString;
    TAxis XAxis = TAxis(DetectorSystematics2D[iProj].XBinEdges.size() - 1, DetectorSystematics2D[iProj].XBinEdges.data());
    TAxis YAxis = TAxis(DetectorSystematics2D[iProj].YBinEdges.size() - 1, DetectorSystematics2D[iProj].YBinEdges.data());
  
    for (auto Sample : DUNEPdfs) {
        std::vector<std::vector<double>> SelectionVector;
        for (size_t iCut = 0; iCut < DetectorSystematics2D[iProj].KinematicCuts.size(); iCut++) {
            std::vector<double> Selection(3);
            Selection[0] = Sample->ReturnKinematicParameterFromString(DetectorSystematics2D[iProj].KinematicCuts[iCut].VarString);
            Selection[1] = DetectorSystematics2D[iProj].KinematicCuts[iCut].Range[0];
            Selection[2] = DetectorSystematics2D[iProj].KinematicCuts[iCut].Range[1];
        
            SelectionVector.emplace_back(Selection);
        }
      
        // Initialize the 2D histogram with the binning specified in the YAML configuration
        TH2D* Hist2D = new TH2D("Hist2D", "2D Histogram",
                                DetectorSystematics2D[iProj].XBinEdges.size() - 1, DetectorSystematics2D[iProj].XBinEdges.data(),
                                DetectorSystematics2D[iProj].YBinEdges.size() - 1, DetectorSystematics2D[iProj].YBinEdges.data());
        
        // Fill the 2D histogram manually using the 1D histograms
        TH1D* XHist = dynamic_cast<TH1D*>(Sample->get1DVarHist(XVarString, SelectionVector, WeightStyle, &XAxis));
        TH1D* YHist = dynamic_cast<TH1D*>(Sample->get1DVarHist(YVarString, SelectionVector, WeightStyle, &YAxis));
        
        if (XHist == nullptr || YHist == nullptr) {
            std::cerr << "Failed to get 1D histograms" << std::endl;
            delete Hist2D;
            continue;
        }
      
        for (int xBin = 1; xBin <= XHist->GetNbinsX(); ++xBin) {
            for (int yBin = 1; yBin <= YHist->GetNbinsY(); ++yBin) {
                double x = XHist->GetXaxis()->GetBinCenter(xBin);
                double y = YHist->GetXaxis()->GetBinCenter(yBin);
                double content = XHist->GetBinContent(xBin) * YHist->GetBinContent(yBin);
                double error = std::sqrt(std::pow(XHist->GetBinError(xBin), 2) + std::pow(YHist->GetBinError(yBin), 2));
                Hist2D->Fill(x, y, content);
                int newBin = Hist2D->FindBin(x, y);
                Hist2D->SetBinError(newBin, error);
            }
        }
      
        // Log the number of entries in the final histogram
        MACH3LOG_INFO("Final Hist2D entries: {}", Hist2D->GetEntries());
      
        Hist2D->SetTitle((DetectorSystematics2D[iProj].XName + " vs " + DetectorSystematics2D[iProj].YName).c_str());
      
        MACH3LOG_INFO("\tSample: {:<20} - Integral: {:<10}", Sample->GetName(), Hist2D->Integral());
        PrintTH2Histogram(Hist2D, Sample->GetName() + "_" + XVarString + "_" + YVarString + ".png");
        delete Hist2D;
        delete XHist;
        delete YHist;
    }
  }
    MACH3LOG_INFO("=================================================");
  }