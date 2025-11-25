// Get main tree
TTree *t = (TTree*)_file0->Get("posteriors");
Long64_t nentries = t->GetEntries();
Long64_t half = nentries / 2;

// Create output files
TFile *f1 = new TFile("Saturday_enubiasenu_newchain_part1.root", "RECREATE");
TFile *f2 = new TFile("Saturday_enubiasenu_newchain_part2.root", "RECREATE");

// --- Copy everything except the posterior tree ---
TIter nextkey(_file0->GetListOfKeys());
TKey *key;
while ((key = (TKey*)nextkey())) {
    TString name = key->GetName();
    if (name == "posteriors") continue; // skip, handled separately

    TObject *obj = key->ReadObj();

    // Write into both files
    f1->cd();
    obj->Write();
    f2->cd();
    obj->Write();

    delete obj;
}

// --- Now copy half the posteriors into each file ---
f1->cd();
TTree *t1 = t->CopyTree("", "", 0, half);
t1->Write();

f2->cd();
TTree *t2 = t->CopyTree("", "", half, nentries);
t2->Write();

// Save and close
f1->Close();
f2->Close();

cout << "✅ Split complete!" << endl;
cout << "Created:" << endl;
cout << "  • Saturday_enubiasenu_newchain_part1.root (first half)" << endl;
cout << "  • Saturday_enubiasenu_newchain_part2.root (second half)" << endl;
