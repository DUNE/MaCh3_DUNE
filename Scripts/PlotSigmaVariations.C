#include "TString.h"
#include "TFile.h"
#include "TH1D.h"

#include <memory>
#include <vector>

// HW: Temporary Hacky macro to plot Sigma Variations for DDAS!
void PlotSigmaVariations(TString file_name){
    auto file = std::make_unique<TFile>(file_name, "OPEN");

    std::vector<TString> labels = {
        "-3#sigma",
        "-1#sigma",
        "Nominal",
        "1#sigma",
        "3#sigma",
    };

    std::vector<int> colours = {
        kBlue,
        kCyan,
        kBlack,
        kMagenta,
        kRed
    };

    std::vector<int> line_styles {
        kDotted,
        kDashed,
        kSolid,
        kDashed,
        kDotted
    };

    TString base_name = "Variation";
    int n_variations = static_cast<int>(labels.size());
    int nominal_variation = 3;

    TString nominal_hist_name;
    nominal_hist_name.Form(base_name+"_%d", nominal_variation-1);

    if(!file || file->IsZombie()){
        std::cerr<<"Cannot find "<<file_name<<std::endl; 
        return;
    }

    auto c = std::make_unique<TCanvas>("c", "c", 1600, 1600);
    c->Print("sigmavar.pdf[");

    auto key_list = file->GetListOfKeys();

    for(int i =0; i<key_list->GetEntries(); ++i){
        TString name = key_list->At(i)->GetName();
        
        TString dir_path = name+"/BeamFD_FHC_Numu";

        auto hist_dir = file->Get<TDirectory>(dir_path);
        if(!hist_dir){
            std::cerr<<"Cannot find "<<dir_path<<" in "<<file_name<<std::endl;
            return;
        }
        
        // Now we can get the histograms
        THStack* s = new THStack(name,name+";E_{rec} (GeV);Variation");

        TLegend* leg = new TLegend(0.6, 0.6, 0.9, 0.9);

        TH1D* nominal = hist_dir->Get<TH1D>(nominal_hist_name);

        double y_min=10;
        double y_max = 0;

        for(int v=0; v<n_variations; ++v){
            TString hist_name;
            hist_name.Form(base_name+"_%d", v);
            // Now we get the histogram
            auto hist = hist_dir->Get<TH1D>(hist_name);
            
            hist->Divide(nominal);

            hist->SetLineColor(colours[v]);
            hist->SetLineStyle(line_styles[v]);
            hist->SetLineWidth(2);

            if(!hist){
                std::cerr<<"Cannot find "<<hist_name<<" in "<<dir_path<<std::endl;
                return;
            }
            for (int b = 1; b <= hist->GetNbinsX(); ++b) {
                double val = hist->GetBinContent(b);
                if (val != 0) {  // skip empty bins
                    y_min = std::min(y_min, val);
                    y_max = std::max(y_max, val);
                }
            }


            leg->AddEntry(hist, labels[v]);
            s->Add(hist);
        }
        
        double padding = (y_max - y_min) * 0.1;
        s->SetMinimum(y_min - padding);  // <-- use SetMinimum/SetMaximum
        s->SetMaximum(y_max + padding);  // <-- NOT GetYaxis()->SetRangeUser
        s->Draw("hist nostack");
        
        leg->Draw();
        c->Update();



        c->Print("sigmavar.pdf");
        delete hist_dir;
    }

    c->Print("sigmavar.pdf]");
    file->Close();
}