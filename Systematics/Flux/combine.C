#include "TH2Jagged/build/Linux/include/TH2Jagged.h"

#pragma cling load("TH2Jagged/build/Linux/lib/libTH2Jagged.so")

void TH2JToDir(TH2Jagged<float> *thj, TDirectory *dir) {
  dir->WriteTObject(thj->fUniformAxis.Clone("OffAxisTAxis"), "OffAxisTAxis");
  for (size_t ui = 0; ui < thj->fUniformAxis.GetNbins(); ++ui) {
    auto thf = thj->NonUniformSlice(ui);
    dir->WriteTObject(
        thf, (std::string(thj->GetName()) + "_" + std::to_string(ui)).c_str());
    delete thf;
  }
}

// adapted from https://root.cern/doc/v632/copyFiles_8C.html
void CopyDir(TDirectory *source, TDirectory *dest, bool pca_limit) {

  // loop on all entries of this directory
  TKey *key;
  // Loop in reverse order to make sure that the order of cycles is
  // preserved.
  TIter nextkey(source->GetListOfKeys(), kIterBackward);
  while ((key = (TKey *)nextkey())) {

    const char *classname = key->GetClassName();
    TClass *cl = gROOT->GetClass(classname);

    if (!cl)
      continue;

    if (cl->InheritsFrom(TDirectory::Class())) {

      if (pca_limit) {
        if (std::string(key->GetName()).substr(0, 9) == "param_pca") {
          int pca_num = std::stoi(std::string(key->GetName()).substr(10));
          if (pca_num >= 20) {
            continue;
          }
        } else {
          continue;
        }
      }

      std::string dirname = key->GetName();
      if (dirname.substr(0, 6) == "param_") {
        dirname = dirname.substr(6);
      }
      if (dirname.substr(0, 3) == "pca") {
        dirname = "HadronProduction_" + dirname;
      }

      std::cout << "Copying directory: " << source->GetName() << "/" << dirname
                << " to " << dest->GetName() << std::endl;
      CopyDir(source->GetDirectory(key->GetName()),
              dest->mkdir(dirname.c_str()), pca_limit);
    } else {
      if (std::string(key->GetName()) == "param_names") {
        continue;
      }
      std::cout << "  writing object: " << key->GetName() << " of type "
                << cl->GetName() << " to " << dest->GetName() << std::endl;
      if (std::string(cl->GetName()) == "TH2Jagged<float>") {
        TH2JToDir(static_cast<TH2Jagged<float> *>(key->ReadObj()),
                  dest->mkdir(key->GetName()));
      } else {
        dest->WriteTObject(key->ReadObj(), key->GetName());
      }
    }
  }
}

void combine() {
  TFile fcomb("flux_variations_FD_and_PRISM_2023.root", "Recreate");
  TFile ffoc("flux_shifts_OffAxis2023.root", "READ");

  TFile fhp("flux_shifts_OffAxis.root", "READ");

  CopyDir(ffoc.GetDirectory("FluxParameters"),
          fcomb.mkdir("FluxParameters")->mkdir("Focussing"), false);

  CopyDir(fhp.GetDirectory("FluxParameters"),
          fcomb.GetDirectory("FluxParameters")->mkdir("HadronProduction"),
          true);

  fcomb.Write();
}