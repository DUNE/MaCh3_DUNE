
void findfixedparams(const char* filename = "/scratch/abipeake/Off_axis_chains_new/LONGERWALLTIME2_q0q3_pTpzEnu_afteradaptive/LONGERWALLTIME2_q0q3_pTpzEnu_afteradaptive_chain_0_job_0.root",
                     const char* outname = "frozen_params.txt") {
    TFile *f = TFile::Open(filename);
    if (!f || f->IsZombie()) {
        std::cerr << "Error: could not open file " << filename << std::endl;
        return;
    }

    TTree *t = (TTree*)f->Get("posteriors");
    if (!t) {
        std::cerr << "Error: could not find TTree 'posteriors'" << std::endl;
        return;
    }

    TObjArray *branches = t->GetListOfBranches();
    std::vector<std::string> frozen;

    for (int i = 0; i < branches->GetEntries(); i++) {
        TString name = branches->At(i)->GetName();

        // check all parameters, not only xsec_*
        auto h = new TH1D("h","",100,-1,1);
        t->Draw(name+">>h","","goff");
        if (h->GetRMS() == 0) {
            frozen.push_back(name.Data());
        }
        delete h;
    }

    f->Close();

    // Write to text file
    std::ofstream out(outname);
    if (!out) {
        std::cerr << "Error: could not open output file " << outname << std::endl;
        return;
    }
    for (size_t i = 0; i < frozen.size(); i++) {
        out << frozen[i];
        if (i < frozen.size() - 1) out << ",";
    }
    out << std::endl;
    out.close();

    std::cout << "Wrote " << frozen.size() << " frozen parameters to " << outname << std::endl;
}
