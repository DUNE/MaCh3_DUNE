void Plot_FluxSigma2D_PRISM(const char* filename =
    "/scratch/abipeake/MaCh3_DUNE_Nov2025/MaCh3_DUNE/eventratetest.root")
{
    TFile* f = TFile::Open(filename);
    if (!f || f->IsZombie()) {
        std::cout << "Could not open file: " << filename << "\n";
        return;
    }

    gStyle->SetOptStat(0);
    gStyle->SetNumberContours(50);

    TCanvas* c = new TCanvas("c", "Flux systematics", 900, 800);
    c->SetRightMargin(0.15);

    // Start multi-page PDF
    c->Print("FluxSystematics_OffAxisND.pdf[");

    // --- Loop over all top-level keys ---
    TIter parIter(f->GetListOfKeys());
    TKey* parKey;

    while ((parKey = (TKey*)parIter())) {

        TString parName = parKey->GetName();

        TDirectory* parDir = (TDirectory*)f->Get(parName);
        if (!parDir) continue;

        TDirectory* offAxisDir = (TDirectory*)parDir->Get("OffAxisND");
        if (!offAxisDir) continue;

        TH2D* hNom  = (TH2D*)offAxisDir->Get("Variation_3");
        TH2D* hPlus = (TH2D*)offAxisDir->Get("Variation_4");

        if (!hNom || !hPlus) {
            std::cout << "Skipping " << parName << ": missing Variation_3 or Variation_4\n";
            continue;
        }

        TH2D* hRatio = (TH2D*)hPlus->Clone(Form("%s_ratio", parName.Data()));
        hRatio->Divide(hNom);
        hRatio->SetTitle(Form("%s (+1#sigma / Nominal)", parName.Data()));
        hRatio->GetXaxis()->SetTitle("Neutrino Energy [GeV]");
        hRatio->GetYaxis()->SetTitle("Off-axis [m]");
        hRatio->SetMinimum(0.98);
        hRatio->SetMaximum(1.01);

        hRatio->Draw("COLZ");

        // Add this plot as a page in the multi-page PDF
        c->Print("FluxSystematics_OffAxisND.pdf");

        delete hRatio;
    }

    // Close the multi-page PDF
    c->Print("FluxSystematics_OffAxisND.pdf]");
    delete c;
    f->Close();
}
