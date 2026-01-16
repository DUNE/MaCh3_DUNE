void Plot_SigmaVariations(const char* filename = "/scratch/abipeake/MaCh3_DUNE_Nov2025/MaCh3_DUNE/eventratetest_template.root") {

    TFile *f = TFile::Open(filename);
    if (!f || f->IsZombie()) { std::cout << "File not found!\n"; return; }
    gStyle->SetOptStat(0);

    // -----------------------------------------
    // Build Tol Incandescent palette
    // -----------------------------------------
    std::vector<int> TolIncandescent;
    std::vector<std::tuple<int,int,int>> incandescentRGB = {
        {255,221, 85},{255,187,34},{255,136,0},
        {255,85,0},{255,34,0},{221,0,0},{170,0,0}
    };

    for (auto &rgb : incandescentRGB) {
        int idx = TColor::GetFreeColorIndex();
        new TColor(idx,
                   std::get<0>(rgb)/255.0,
                   std::get<1>(rgb)/255.0,
                   std::get<2>(rgb)/255.0);
        TolIncandescent.push_back(idx);
    }

    std::vector<double> sigmas = {-3, -2, -1,  0, 1, 2 , 3};
    std::vector<int> sigmaToTol = {0,1, 2, 4, 5, 6};

    std::vector<int> lineColors;
    for (int i : sigmaToTol) lineColors.push_back(TolIncandescent[i]);

    // -----------------------------------------
    // Canvas with ratio split
    // -----------------------------------------
    TCanvas *c = new TCanvas("c","c",900,700);
    c->Divide(1,2);

    TPad *pad1 = (TPad*)c->cd(1);
    pad1->SetPad(0, 0.30, 1, 1.0);
    pad1->SetBottomMargin(0.02);

    TPad *pad2 = (TPad*)c->cd(2);
    pad2->SetPad(0, 0.00, 1, 0.30);
    pad2->SetTopMargin(0.05);
    pad2->SetBottomMargin(0.28);

    c->Print("SystematicsNuWro.pdf[");

    // -----------------------------------------
    // Loop directories
    // -----------------------------------------
    TIter parDirIter(f->GetListOfKeys());
    TKey *parKey;

    while ((parKey = (TKey*)parDirIter())) {

        TString parName = parKey->GetName();
        TDirectory *parDir = (TDirectory*)f->Get(parName);
        if (!parDir) continue;

        TIter sampleIter(parDir->GetListOfKeys());
        TKey *sampleKey;

        while ((sampleKey = (TKey*)sampleIter())) {

            TString sampleName = sampleKey->GetName();
            TDirectory *sampleDir = (TDirectory*)parDir->Get(sampleName);
            if (!sampleDir) continue;

            // --------------------------
            // Legend
            // --------------------------
            TLegend *leg = new TLegend(0.65,0.6,0.88,0.88);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);

            bool first = true;
            TH1 *h0 = nullptr;

            // --------------------------
            // Draw top pad
            // --------------------------
            pad1->cd();

            for (size_t i = 0; i < sigmas.size(); ++i) {

                TH1 *h = (TH1*)sampleDir->Get(Form("Variation_%d", (int)i));
                if (!h) continue;

                h->SetLineColor(lineColors[i]);
                h->SetLineWidth(3);
                h->SetMinimum(0);
                h->Scale(1.0/(h->Integral()),"width");
                double max = h->GetMaximum();
                h->SetMaximum(1.1*max);

                if (first) {
                    h0 = (TH1*)h->Clone("ref");
                    h0->SetTitle(Form("%s -- %s", parName.Data(), sampleName.Data()));
                    h0->Draw("hist");
                    first = false;
                } else {
                    h->Draw("hist same");
                }

                leg->AddEntry(h, Form("%+.0f#sigma", sigmas[i]), "l");
            }

            leg->Draw();

            // --------------------------
            // Ratio pad
            // --------------------------
            if (h0) {
                pad2->cd();

                for (size_t i = 0; i < sigmas.size(); ++i) {

                    TH1 *h = (TH1*)sampleDir->Get(Form("Variation_%d", (int)i));
                    if (!h) continue;

                    TH1 *r = (TH1*)h->Clone(Form("ratio_%d", (int)i));
                    r->Divide(h0);

                    r->SetLineColor(lineColors[i]);
                    r->SetLineWidth(2);

                    r->GetYaxis()->SetTitle("Var/Nom");
                    r->GetYaxis()->SetTitleSize(0.08);
                    r->GetYaxis()->SetTitleOffset(0.5);
                    r->GetYaxis()->SetLabelSize(0.07);

                    r->GetXaxis()->SetTitleSize(0.10);
                    r->GetXaxis()->SetLabelSize(0.09);

                    r->SetMinimum(0.7);
                    r->SetMaximum(1.1);

                    if (i == 2) { // nominal
                        r->Draw("hist");
                    } else {
                        r->Draw("hist same");
                    }
                }
            }

            c->Print("SystematicsNuWro.pdf");
            delete leg;
        }
    }

    c->Print("SystematicsNuwro.pdf]");
}
