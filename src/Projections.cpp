#include <TDirectory.h>
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
	std::vector<std::string> VarStrings;
	std::vector<std::vector<double>> BinEdges;

	std::vector<KinematicCut> KinematicCuts;
	std::vector<CategoryCut> CategoryCuts;
};

std::string ReturnFormattedHistogramNameFromProjection(ProjectionVariable Proj) {
	std::string ReturnStr;

	for (size_t iKinematicCut=0;iKinematicCut<Proj.KinematicCuts.size();iKinematicCut++) {
		if (iKinematicCut > 0) {
			ReturnStr += " && ";
		}
		ReturnStr += Form("(%4.2f < %s < %4.2f)",Proj.KinematicCuts[iKinematicCut].Range[0],Proj.KinematicCuts[iKinematicCut].Name.c_str(),Proj.KinematicCuts[iKinematicCut].Range[1]);
	}
	std::string y_axis_title;
	if (Proj.VarStrings.size()==1) {y_axis_title = "Events";}
	else {y_axis_title = Proj.VarStrings[1];}
	ReturnStr += Proj.Name+";"+Proj.VarStrings[0]+";"+y_axis_title+";";
	return ReturnStr;
}

void WriteTH1Histogram(TH1 *Hist, TDirectory *Dir, std::string Name) {
	Dir->cd();
	Hist->Write(Name.c_str());
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

void WriteTHStackHistogram(THStack *Hist, TDirectory *Dir, std::string Name) {
	Dir->cd();
	Hist->Write(Name.c_str());
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

	std::string OutputFileName = fitMan->raw()["General"]["OutputFile"].as<std::string>();
	TFile *File = TFile::Open(OutputFileName.c_str(), "RECREATE");

	for (auto &ProjectionConfig : fitMan->raw()["Projections"]) {
		std::string VarName = ProjectionConfig["Name"].as<std::string>();
		//JM now a vector of size 1 (for 1d hists) or 2 (for 2d hists)
		std::vector<std::string> VarStrings = ProjectionConfig["VarStrings"].as< std::vector<std::string> >();

		//Could replace this with uniform [lbin, hbin, nbins] for example
		//JM have included this option (as [nbins, lbin, hbin])
		std::vector<std::vector<double>> VarBinnings = ProjectionConfig["VarBins"].as<std::vector<std::vector<double>>>();

		std::vector<KinematicCut> KinematicCuts;
		std::vector<CategoryCut> CategoryCuts;

		if ((VarStrings.size()!=1 && VarStrings.size()!=2) || VarStrings.size() != VarBinnings.size()) {
			MACH3LOG_ERROR("Projections: {} VarStrings specified, {} VarBinnings specified. Specify 1 or 2 of both.", VarStrings.size(), VarBinnings.size());
			throw MaCh3Exception(__FILE__,__LINE__);
		}
		for (int iBinning=0; iBinning<VarBinnings.size(); iBinning++) {
			if (VarBinnings[iBinning].size() == 3) {
				double nbins = VarBinnings[iBinning][0];
				double xmin = VarBinnings[iBinning][1];
				double xmax = VarBinnings[iBinning][2];
				double step = (xmax-xmin)/(nbins-1);
				VarBinnings[iBinning] = {};
				std::cout << nbins << " " << xmin << " " << xmax << std::endl;
				for (double iBinEdge=xmin; iBinEdge<=xmax; iBinEdge+=step) {
					VarBinnings[iBinning].push_back(iBinEdge);
					//std::cout<<iBinEdge;
				}
			}
		}

		for (auto &KinematicCutConfig: ProjectionConfig["KinematicCuts"]) {
			std::string KinematicCutName = KinematicCutConfig["Name"].as<std::string>();
			std::string KinematicCutVarString = KinematicCutConfig["VarString"].as<std::string>();
			std::vector<double> KinematicCutRange = KinematicCutConfig["Range"].as< std::vector<double> >();

			KinematicCut Cut = KinematicCut{KinematicCutName,KinematicCutVarString,KinematicCutRange};
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

			CategoryCut Cut = CategoryCut{CategoryCutName,CategoryCutVarString,CategoryCutBreakdown,CategoryCutColours,CategoryCutNames};
			CategoryCuts.emplace_back(Cut);
		}

		ProjectionVariable Proj = ProjectionVariable{VarName,VarStrings,VarBinnings,KinematicCuts,CategoryCuts};
		Projections.emplace_back(Proj);
	}

	MACH3LOG_INFO("=================================================");
	MACH3LOG_INFO("Projections pulled from Config..");
	MACH3LOG_INFO("================================");

	for (size_t iProj=0;iProj<Projections.size();iProj++) {
		int histdim = Projections[iProj].VarStrings.size();

		if (histdim == 1) {
			MACH3LOG_INFO("Projection {:<2} - Name : {} , VarString : {} , Binning : {}, {}, {}"
					,iProj,Projections[iProj].Name,
					Projections[iProj].VarStrings[0],Projections[iProj].BinEdges[0].size(),Projections[iProj].BinEdges[0][0],Projections[iProj].BinEdges[0][Projections[iProj].BinEdges[0].size()-1]);
		} else {
			MACH3LOG_INFO("Projection {:<2} - Name : {} , VarString1 : {} , Binning : {} , {} , {} , VarString2 : {} , Binning : {}, {}, {}"
					,iProj,Projections[iProj].Name,
					Projections[iProj].VarStrings[0],Projections[iProj].BinEdges[0].size(),Projections[iProj].BinEdges[0][0],Projections[iProj].BinEdges[0][Projections[iProj].BinEdges[0].size()-1],
					Projections[iProj].VarStrings[1],Projections[iProj].BinEdges[1].size(),Projections[iProj].BinEdges[1][0],Projections[iProj].BinEdges[1][Projections[iProj].BinEdges[1].size()-1]);
		}

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

	TH1* Hist;
	THStack* Stack;

	for (size_t iProj=0;iProj<Projections.size();iProj++) {
		MACH3LOG_INFO("================================");
		MACH3LOG_INFO("Projection {}/{}",iProj+1,Projections.size());

		std::vector<std::string> ProjectionVar_Str = Projections[iProj].VarStrings;
		int histdim = ProjectionVar_Str.size();
		TAxis AxisX = TAxis(Projections[iProj].BinEdges[0].size()-1,Projections[iProj].BinEdges[0].data());
		TAxis AxisY;
		if (histdim == 2) {AxisY = TAxis(Projections[iProj].BinEdges[1].size()-1,Projections[iProj].BinEdges[1].data());}

		for (auto Sample : DUNEPdfs) {

			File->mkdir(Sample->GetName().c_str());
			TDirectory *dir = File->GetDirectory(Sample->GetName().c_str());
			std::vector< std::vector<double> > SelectionVector;
			for (size_t iCut=0;iCut<Projections[iProj].KinematicCuts.size();iCut++) {
				std::vector<double> Selection(3);
				Selection[0] = Sample->ReturnKinematicParameterFromString(Projections[iProj].KinematicCuts[iCut].VarString);
				Selection[1] = Projections[iProj].KinematicCuts[iCut].Range[0];
				Selection[2] = Projections[iProj].KinematicCuts[iCut].Range[1];

				SelectionVector.emplace_back(Selection);
			}

			std::string outputname;
			if (histdim==1) {
				Hist = Sample->get1DVarHist(ProjectionVar_Str[0],SelectionVector,WeightStyle,&AxisX);
				outputname = Sample->GetName()+"_"+ProjectionVar_Str[0]+".png";
			} 
			else {
				Hist = (TH1*)Sample->get2DVarHist(ProjectionVar_Str[0],ProjectionVar_Str[1],SelectionVector,WeightStyle,&AxisX,&AxisY);
				outputname = Sample->GetName()+"_"+ProjectionVar_Str[0]+"_"+ProjectionVar_Str[1]+".png";
			}
			Hist->Scale(1.0,"Width");
			Hist->SetTitle(ReturnFormattedHistogramNameFromProjection(Projections[iProj]).c_str());
			MACH3LOG_INFO("\tSample: {:<20} - Integral: {:<10}",Sample->GetName(),Hist->Integral());
			PrintTH1Histogram(Hist,outputname);
			WriteTH1Histogram(Hist, dir, Sample->GetName() + "_" + ProjectionVar_Str);

			for (size_t iCat=0;iCat<Projections[iProj].CategoryCuts.size();iCat++) {
				MACH3LOG_INFO("\t\tCategory: {:<10} - Name : {:<20}",iCat,Projections[iProj].CategoryCuts[iCat].Name);

				Stack = new THStack(Projections[iProj].CategoryCuts[iCat].Name.c_str(), ReturnFormattedHistogramNameFromProjection(Projections[iProj]).c_str());
				TLegend *leg = new TLegend(0.6, 0.7, 0.88, 0.88);
				leg->SetHeader(Projections[iProj].CategoryCuts[iCat].Name.c_str());

				for (size_t iBreak = 0; iBreak < Projections[iProj].CategoryCuts[iCat].Breakdown.size(); iBreak++) {

					TH1 *BreakdownHist = nullptr;

					for (size_t iGroup=0;iGroup<Projections[iProj].CategoryCuts[iCat].Breakdown[iBreak].size();iGroup++) {
						std::vector< std::vector<double> > SelectionVector_IncCategory = std::vector< std::vector<double> >(SelectionVector);

						std::vector<double> Selection(3);
						Selection[0] = Sample->ReturnKinematicParameterFromString(Projections[iProj].CategoryCuts[iCat].VarString);
						Selection[1] = Projections[iProj].CategoryCuts[iCat].Breakdown[iBreak][iGroup];
						Selection[2] = Projections[iProj].CategoryCuts[iCat].Breakdown[iBreak][iGroup]+1;
						SelectionVector_IncCategory.emplace_back(Selection);

						if (histdim==1) {
							Hist = Sample->get1DVarHist(ProjectionVar_Str[0],SelectionVector_IncCategory,WeightStyle,&AxisX);
						}	else {
							Hist = (TH1*)Sample->get2DVarHist(ProjectionVar_Str[0],ProjectionVar_Str[1],SelectionVector_IncCategory,WeightStyle,&AxisX,&AxisY);
						}
						Hist->SetFillColor(Projections[iProj].CategoryCuts[iCat].Colours[iBreak]);
						Hist->Scale(1.0,"Width");

						if (BreakdownHist == nullptr) {
							BreakdownHist = Hist;
						} else {
							BreakdownHist->Add(Hist);
						}
						leg->AddEntry(Hist, Projections[iProj].CategoryCuts[iCat].CategoryNames[iBreak].c_str(), "f");
					}

					MACH3LOG_INFO("\t\t\tBreakdown: {:<10} - Integral: {:<10}", iBreak, BreakdownHist->Integral());
					Stack->Add(BreakdownHist);
				}

				if (histdim==1) {
					PrintTHStackHistogram(Stack,Sample->GetName()+"_"+ProjectionVar_Str[0]+"_"+Projections[iProj].CategoryCuts[iCat].Name+"_Stack.png");
				}
				else {
					PrintTHStackHistogram(Stack,Sample->GetName()+"_"+ProjectionVar_Str[0]+"_"+ProjectionVar_Str[1]+"_"+Projections[iProj].CategoryCuts[iCat].Name+"_Stack.png");
				}
				// HH: Add the stack to the file
				TCanvas Canv = TCanvas();
				Stack->Draw("HIST");
				leg->Draw();
				Canv.Write((ProjectionVar_Str + "_" + Projections[iProj].CategoryCuts[iCat].Name + "_Stack_Canvas").c_str());
				WriteTHStackHistogram(Stack, dir, ProjectionVar_Str + "_" + Projections[iProj].CategoryCuts[iCat].Name + "_Stack");
			}
		}
	}
	File->Close();
	MACH3LOG_INFO("=================================================");
}
