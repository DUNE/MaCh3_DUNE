

#include "Samples/MaCh3DUNEFactory.h"
#include "Samples/StructsDUNE.h"
#include "Fitters/FitterBase.h"
#include "Manager/Manager.h"
#include "Parameters/ParameterHandlerBase.h"

#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <filesystem>
#include <map>
#include "Samples/MaCh3DUNEFactory.h"
#include "Samples/StructsDUNE.h"
#include "Samples/SampleHandlerBase.h"
#include "Fitters/FitterBase.h"
#include "Manager/Manager.h"
#include "Parameters/ParameterHandlerBase.h"

#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <filesystem>
#include <map>

void MakeSpectaVariations2D(
    SampleHandlerFD* pdf,
    const std::string& varX,
    const std::string& varY,
    TFile* fout,
    const std::string& ND_or_FD,
    const std::string& pdfTitle,
    int p,
    int sample_idx = 0,
    std::map<std::string,
    std::vector<std::vector<std::vector<double>>>>* posteriorStore = nullptr)
{
    // ---------------------------------------------------------
    // Empty cuts
    // ---------------------------------------------------------

    std::vector<KinematicCut> emptyEventCuts;
    std::vector<KinematicCut> emptySubEventCuts;

    // ---------------------------------------------------------
    // Get histogram
    // ---------------------------------------------------------

    TH2* h = pdf->Get2DVarHist(
        sample_idx,
        varX,
        varY,
        emptyEventCuts,
        0,              // weighted
        nullptr,        // x-axis binning
        nullptr,        // y-axis binning
        emptySubEventCuts
    );

    if (!h) {
        std::cerr
            << "[WARN] Null TH2 for "
            << varX << " vs " << varY
            << " from " << pdfTitle
            << std::endl;
        return;
    }

    TH2D* h2d = dynamic_cast<TH2D*>(h);

    if (!h2d) {
        std::cerr
            << "[WARN] Could not cast TH2 to TH2D for "
            << varX << " vs " << varY
            << std::endl;

        delete h;
        return;
    }

    // ---------------------------------------------------------
    // VERY IMPORTANT
    // Prevent ROOT ownership issues
    // ---------------------------------------------------------

    h2d->SetDirectory(nullptr);

    // ---------------------------------------------------------
    // Clean histogram name
    // ---------------------------------------------------------

    std::string cleanTitle = pdfTitle;

    std::replace(
        cleanTitle.begin(),
        cleanTitle.end(),
        ' ',
        '_'
    );

    std::replace(
        cleanTitle.begin(),
        cleanTitle.end(),
        '/',
        '_'
    );

    // ---------------------------------------------------------
    // Output names
    // ---------------------------------------------------------

    std::string baseDir;
    std::string histName;

    if (p < 0) {

        baseDir =
            "Asimov/" +
            ND_or_FD +
            "/2D";

        histName = Form(
            "%s_%s_%s_vs_%s_Asimov",
            ND_or_FD.c_str(),
            cleanTitle.c_str(),
            varX.c_str(),
            varY.c_str()
        );
    }
    else {

        baseDir =
            ND_or_FD +
            "/2D";

        histName = Form(
            "%s_%s_%s_vs_%s_posterior_toy_%03d",
            ND_or_FD.c_str(),
            cleanTitle.c_str(),
            varX.c_str(),
            varY.c_str(),
            p
        );
    }

    // ---------------------------------------------------------
    // Create directory safely
    // ---------------------------------------------------------

    TDirectory* dir =
        fout->GetDirectory(baseDir.c_str());

    if (!dir)
        dir = fout->mkdir(baseDir.c_str());

    dir->cd();

    // ---------------------------------------------------------
    // Clone histogram safely
    // ---------------------------------------------------------

    TH2D* cloneH =
        static_cast<TH2D*>(
            h2d->Clone(histName.c_str())
        );

    cloneH->SetDirectory(dir);

    cloneH->Write();

    // ---------------------------------------------------------
    // Optional posterior storage
    // ---------------------------------------------------------

    if (p >= 0 && posteriorStore) {

        std::string key =
            ND_or_FD + "_" +
            cleanTitle + "_" +
            varX + "_vs_" + varY;

        int nX = cloneH->GetNbinsX();
        int nY = cloneH->GetNbinsY();

        auto& grid = (*posteriorStore)[key];

        // initialize once only
        if (grid.empty()) {

            grid.resize(nX);

            for (int ix = 0; ix < nX; ++ix)
                grid[ix].resize(nY);
        }

        for (int ix = 1; ix <= nX; ++ix) {

            for (int iy = 1; iy <= nY; ++iy) {

                double val =
                    cloneH->GetBinContent(ix, iy);

                grid[ix-1][iy-1]
                    .push_back(val);
            }
        }
    }

    fout->cd();

    // ---------------------------------------------------------
    // Cleanup
    // ---------------------------------------------------------

    delete cloneH;
    delete h2d;
}