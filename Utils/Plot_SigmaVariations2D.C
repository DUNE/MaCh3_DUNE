void Plot_SigmaVariations2D(
    const char* filename = "/scratch/abipeake/MaCh3_DUNE_Nov2025/MaCh3_DUNE/eventratetest.root")
{
    TFile *f = TFile::Open(filename);
    if (!f || f->IsZombie()) {
        std::cout << "File not found!\n";
        return;
    }

    gStyle->SetOptStat(0);
    gStyle->SetNumberContours(50);

    std::vector<double> sigmas = {-3, -1, 0, 1, 3};


    TCanvas *c = new TCanvas("c","c",900,900);
    c->Divide(1,2);

    TPad *pad1 = (TPad*)c->cd(1);
    pad1->SetPad(0, 0.35, 1, 1.0);
    pad1->SetBottomMargin(0.02);
    pad1->SetRightMargin(0.15);

    TPad *pad2 = (TPad*)c->cd(2);
    pad2->SetPad(0, 0.00, 1, 0.35);
    pad2->SetTopMargin(0.05);
    pad2->SetBottomMargin(0.28);
    pad2->SetRightMargin(0.15);

    c->Print("Systematics_2D.pdf[");

    
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

            
            // Nominal (sigma = 0)
            TH2 *hNom = (TH2*)sampleDir->Get("Variation_2");
            if (!hNom) continue;

            // +1 sigma only (index 3)
                const size_t i = 3;

                TH2 *hPlus = (TH2*)sampleDir->Get("Variation_3");
                if (!hPlus) continue;

                TH2 *hVar   = (TH2*)hPlus->Clone("hVar");
                TH2 *hRatio = (TH2*)hPlus->Clone("hRatio");
                hRatio->Divide(hNom);

                // Titles
                hVar->SetTitle(Form("%s — %s (+1#sigma)",
                                    parName.Data(),
                                    sampleName.Data()));

                hRatio->SetTitle("+1#sigma / Nominal");

                // Draw variation
                pad1->cd();
                hVar->GetZaxis()->SetTitle("Events");
                hVar->Draw("COLZ");

                // Draw ratio
                pad2->cd();
                hRatio->GetZaxis()->SetTitle("+1#sigma / Nom");
                hRatio->SetMinimum(0.5);
                hRatio->SetMaximum(0.996);
                hRatio->Draw("COLZ");

                c->Print("Systematics_2D.pdf");

                delete hVar;
                delete hRatio;


            
            // for (size_t i = 0; i < sigmas.size(); ++i) {

            //     TH2 *h = (TH2*)sampleDir->Get(Form("Variation_%zu", i));
            //     if (!h) continue;

            //     // Clone to be safe
            //     TH2 *hVar = (TH2*)h->Clone("hVar");
            //     TH2 *hRatio = (TH2*)hVar->Clone("hRatio");
            //     hRatio->Divide(hNom);

                
            //     hVar->SetTitle(Form("%s — %s (%+.0f#sigma)",
            //                         parName.Data(),
            //                         sampleName.Data(),
            //                         sigmas[i]));

            //     hRatio->SetTitle("Variation / Nominal");

               
            //     pad1->cd();
            //     hVar->GetZaxis()->SetTitle("Events");
            //     hVar->Draw("COLZ");

                
            //     pad2->cd();
            //     hRatio->GetZaxis()->SetTitle("Var / Nom");
            //     hRatio->SetMinimum(0.5);
            //     hRatio->SetMaximum(1.5);
            //     hRatio->Draw("COLZ");

            //     c->Print("Systematics_2D.pdf");

            //     delete hVar;
            //     delete hRatio;
            // }
        }
    }

    c->Print("Systematics_2D.pdf]");
}
