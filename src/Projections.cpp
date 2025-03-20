#include <vector>

#include <TCanvas.h>
#include <TH1D.h>
#include <THStack.h>
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
  std::vector<std::vector<double>> Breakdown;
  std::vector<double> Colours;
  std::vector<std::string> CategoryNames;
};

struct ProjectionVariable {
  std::string Name;
  std::string VarString;
  std::vector<double> BinEdges;

  std::vector<KinematicCut> KinematicCuts;
  std::vector<CategoryCut> CategoryCuts;
};

std::string ReturnFormattedHistogramNameFromProjection(ProjectionVariable Proj) {
  std::string ReturnStr;

  for (size_t iKinematicCut = 0; iKinematicCut < Proj.KinematicCuts.size(); iKinematicCut++) {
    if (iKinematicCut > 0) {
      ReturnStr += " && ";
    }
    ReturnStr += Form("(%4.2f < %s < %4.2f)", Proj.KinematicCuts[iKinematicCut].Range[0], Proj.KinematicCuts[iKinematicCut].Name.c_str(), Proj.KinematicCuts[iKinematicCut].Range[1]);
  }

  ReturnStr += ";" + Proj.Name + ";Events";
  return ReturnStr;
}

void PrintTH1Histogram(TH1 *Hist, std::string OutputName) {
  TCanvas Canv = TCanvas();
  Hist->Draw();
  Canv.Print(OutputName.c_str());
}

void PrintCategoryLegends(std::vector<ProjectionVariable> Projections) {
  TLegend Legend = TLegend(0.1, 0.1, 0.9, 0.9);

  TCanvas Canv = TCanvas();

  std::vector<TH1D *> HistVec;

  for (size_t iProj = 0; iProj < Projections.size(); iProj++) {
    for (size_t iCat = 0; iCat < Projections[iProj].CategoryCuts.size(); iCat++) {
      CategoryCut Cat = Projections[iProj].CategoryCuts[iCat];

      Legend.SetHeader(Projections[iProj].CategoryCuts[iCat].Name.c_str());

      HistVec.resize(Projections[iProj].CategoryCuts[iCat].Breakdown.size());
      for (size_t iBreak = 0; iBreak < Cat.Breakdown.size(); iBreak++) {
        HistVec[iBreak] = new TH1D(Form("DummyHist_%i", (int)iBreak), "", 1, 0, 1);
      }

      for (size_t iBreak = 0; iBreak < Cat.Breakdown.size(); iBreak++) {
        HistVec[iBreak]->SetFillColor(Cat.Colours[iBreak]);
        Legend.AddEntry(HistVec[iBreak], Cat.CategoryNames[iBreak].c_str(), "f");
      }

      Legend.Draw();
      Canv.Print(("Legend_" + Projections[iProj].CategoryCuts[iCat].Name + ".png").c_str());
      Legend.Clear();

      for (size_t iBreak = 0; iBreak < Cat.Breakdown.size(); iBreak++) {
        delete HistVec[iBreak];
      }
    }
  }
}

void PrintTHStackHistogram(THStack *Hist, std::string OutputName) {
  TCanvas Canv = TCanvas();
  Hist->Draw("HIST");
  Canv.Print(OutputName.c_str());
}

int main(int argc, char *argv[]) {
  if (argc == 1) {
    MACH3LOG_ERROR("Usage: bin/EventRatesDUNEBeam config.cfg");
    return 1;
  }
  auto fitMan = std::unique_ptr<manager>(new manager(argv[1]));

  int WeightStyle = 1;

  // ###############################################################################################################################
  // Create samplePDFFD objects

  covarianceXsec *xsec = nullptr;
  covarianceOsc *osc = nullptr;

  std::vector<samplePDFFDBase *> DUNEPdfs;
  MakeMaCh3DuneInstance(fitMan.get(), DUNEPdfs, xsec, osc);

  // ###############################################################################################################################
  // Perform reweight and print total integral for sanity check

  MACH3LOG_INFO("=================================================");
  std::vector<TH1D *> DUNEHists;
  for (auto Sample : DUNEPdfs) {
    Sample->reweight();
    DUNEHists.push_back(Sample->get1DHist());

    std::string EventRateString = fmt::format("{:.2f}", Sample->get1DHist()->Integral());
    MACH3LOG_INFO("Event rate for {} : {:<5}", Sample->GetName(), EventRateString);
  }

  // ###############################################################################################################################
  // Grab Projections from the config

  std::vector<ProjectionVariable> Projections;

  for (auto &ProjectionConfig : fitMan->raw()["Projections"]) {
    std::string VarName = ProjectionConfig["Name"].as<std::string>();
    std::string VarString = ProjectionConfig["VarString"].as<std::string>();

    // Could replace this with uniform [lbin, hbin, nbins] for example
    std::vector<double> VarBinEdges = ProjectionConfig["VarBins"].as<std::vector<double>>();

    std::vector<KinematicCut> KinematicCuts;
    std::vector<CategoryCut> CategoryCuts;

    for (auto &KinematicCutConfig : ProjectionConfig["KinematicCuts"]) {
      std::string KinematicCutName = KinematicCutConfig["Name"].as<std::string>();
      std::string KinematicCutVarString = KinematicCutConfig["VarString"].as<std::string>();
      std::vector<double> KinematicCutRange = KinematicCutConfig["Range"].as<std::vector<double>>();

      KinematicCut Cut = KinematicCut{KinematicCutName, KinematicCutVarString, KinematicCutRange};
      KinematicCuts.emplace_back(Cut);
    }

    for (auto &CategoryCutConfig : ProjectionConfig["CategoryCuts"]) {
      std::string CategoryCutName = CategoryCutConfig["Name"].as<std::string>();
      std::string CategoryCutVarString = CategoryCutConfig["VarString"].as<std::string>();
      std::vector<std::vector<double>> CategoryCutBreakdown = CategoryCutConfig["Breakdown"].as<std::vector<std::vector<double>>>();

      std::vector<double> CategoryCutColours;
      if (CategoryCutConfig["Colours"]) {
        CategoryCutColours = CategoryCutConfig["Colours"].as<std::vector<double>>();
      } else {
        CategoryCutColours.resize(CategoryCutBreakdown.size());
        int colour = 20.;
        for (size_t iColour = 0; iColour < CategoryCutColours.size(); iColour++) {
          CategoryCutColours[iColour] = colour;
          colour += 4;
          if (colour > 50) {
            colour -= 30;
          }
        }
      }

      std::vector<std::string> CategoryCutNames;
      if (CategoryCutConfig["Names"]) {
        CategoryCutNames = CategoryCutConfig["Names"].as<std::vector<std::string>>();
      } else {
        CategoryCutNames.resize(CategoryCutBreakdown.size());
        for (size_t i = 0; i < CategoryCutBreakdown.size(); i++) {
          CategoryCutNames[i] = Form("%i", (int)i);
        }
      }

      CategoryCut Cut = CategoryCut{CategoryCutName, CategoryCutVarString, CategoryCutBreakdown, CategoryCutColours, CategoryCutNames};
      CategoryCuts.emplace_back(Cut);
    }

    ProjectionVariable Proj = ProjectionVariable{VarName, VarString, VarBinEdges, KinematicCuts, CategoryCuts};
    Projections.emplace_back(Proj);
  }

  MACH3LOG_INFO("=================================================");
  MACH3LOG_INFO("Projections pulled from Config..");
  MACH3LOG_INFO("================================");

  for (size_t iProj = 0; iProj < Projections.size(); iProj++) {
    MACH3LOG_INFO("Projection {:<2} - Name : {:<20} , VarString : {:<20}", iProj, Projections[iProj].Name, Projections[iProj].VarString);
    MACH3LOG_INFO("\t\tBinning: {}", fmt::join(Projections[iProj].BinEdges, ", "));

    if (Projections[iProj].KinematicCuts.size() > 0) {
      MACH3LOG_INFO("\t\tKinematicCuts:");
      for (size_t iCut = 0; iCut < Projections[iProj].KinematicCuts.size(); iCut++) {
        MACH3LOG_INFO("\t\t\tCut {:<2} - Name : {:<20} , Lower Bound : {:<10} , Upper Bound : {:<10}", iCut, Projections[iProj].KinematicCuts[iCut].Name,
                      Projections[iProj].KinematicCuts[iCut].Range[0], Projections[iProj].KinematicCuts[iCut].Range[1]);
      }
    }

    if (Projections[iProj].CategoryCuts.size() > 0) {
      MACH3LOG_INFO("\t\tCategoryCuts:");
      for (size_t iCut = 0; iCut < Projections[iProj].CategoryCuts.size(); iCut++) {

        std::vector<std::string> BreakdownStrs(Projections[iProj].CategoryCuts[iCut].Breakdown.size());
        for (size_t iBreak = 0; iBreak < Projections[iProj].CategoryCuts[iCut].Breakdown.size(); iBreak++) {
          BreakdownStrs[iBreak] = fmt::format("{}", fmt::join(Projections[iProj].CategoryCuts[iCut].Breakdown[iBreak], ", "));
        }
        MACH3LOG_INFO("\t\t\tCategory {:<2} - Name : {:<20} , Category Breakdown : {}", iCut, Projections[iProj].CategoryCuts[iCut].Name, fmt::join(BreakdownStrs, ", "));
      }
    }
    MACH3LOG_INFO("================================");
  }

  PrintCategoryLegends(Projections);

  // ###############################################################################################################################
  // Make the plots..

  MACH3LOG_INFO("=================================================");
  MACH3LOG_INFO("Building Projections..");

  TH1 *Hist;
  THStack *Stack;

  for (size_t iProj = 0; iProj < Projections.size(); iProj++) {
    MACH3LOG_INFO("================================");
    MACH3LOG_INFO("Projection {}/{}", iProj, Projections.size());

    std::string ProjectionVar_Str = Projections[iProj].VarString;
    TAxis Axis = TAxis(Projections[iProj].BinEdges.size() - 1, Projections[iProj].BinEdges.data());

    for (auto Sample : DUNEPdfs) {

      std::vector<std::vector<double>> SelectionVector;
      for (size_t iCut = 0; iCut < Projections[iProj].KinematicCuts.size(); iCut++) {
        std::vector<double> Selection(3);
        Selection[0] = Sample->ReturnKinematicParameterFromString(Projections[iProj].KinematicCuts[iCut].VarString);
        Selection[1] = Projections[iProj].KinematicCuts[iCut].Range[0];
        Selection[2] = Projections[iProj].KinematicCuts[iCut].Range[1];

        SelectionVector.emplace_back(Selection);
      }
      Hist = Sample->get1DVarHist(ProjectionVar_Str, SelectionVector, WeightStyle, &Axis);
      Hist->Scale(1.0, "Width");
      Hist->SetTitle(ReturnFormattedHistogramNameFromProjection(Projections[iProj]).c_str());

      MACH3LOG_INFO("\tSample: {:<20} - Integral: {:<10}", Sample->GetName(), Hist->Integral());
      PrintTH1Histogram(Hist, Sample->GetName() + "_" + ProjectionVar_Str + ".png");

      for (size_t iCat = 0; iCat < Projections[iProj].CategoryCuts.size(); iCat++) {
        MACH3LOG_INFO("\t\tCategory: {:<10} - Name : {:<20}", iCat, Projections[iProj].CategoryCuts[iCat].Name);

        Stack = new THStack(Projections[iProj].CategoryCuts[iCat].Name.c_str(), ReturnFormattedHistogramNameFromProjection(Projections[iProj]).c_str());

        for (size_t iBreak = 0; iBreak < Projections[iProj].CategoryCuts[iCat].Breakdown.size(); iBreak++) {

          TH1 *BreakdownHist = nullptr;

          for (size_t iGroup = 0; iGroup < Projections[iProj].CategoryCuts[iCat].Breakdown[iBreak].size(); iGroup++) {
            std::vector<std::vector<double>> SelectionVector_IncCategory = std::vector<std::vector<double>>(SelectionVector);

            std::vector<double> Selection(3);
            Selection[0] = Sample->ReturnKinematicParameterFromString(Projections[iProj].CategoryCuts[iCat].VarString);
            Selection[1] = Projections[iProj].CategoryCuts[iCat].Breakdown[iBreak][iGroup];
            Selection[2] = Projections[iProj].CategoryCuts[iCat].Breakdown[iBreak][iGroup] + 1;
            SelectionVector_IncCategory.emplace_back(Selection);

            Hist = Sample->get1DVarHist(ProjectionVar_Str, SelectionVector_IncCategory, WeightStyle, &Axis);
            Hist->SetFillColor(Projections[iProj].CategoryCuts[iCat].Colours[iBreak]);
            Hist->Scale(1.0, "Width");

            if (BreakdownHist == nullptr) {
              BreakdownHist = Hist;
            } else {
              BreakdownHist->Add(Hist);
            }
          }

          MACH3LOG_INFO("\t\t\tBreakdown: {:<10} - Integral: {:<10}", iBreak, BreakdownHist->Integral());
          Stack->Add(BreakdownHist);
        }

        PrintTHStackHistogram(Stack, Sample->GetName() + "_" + ProjectionVar_Str + "_" + Projections[iProj].CategoryCuts[iCat].Name + "_Stack.png");
      }
    }
  }

  MACH3LOG_INFO("=================================================");
}
