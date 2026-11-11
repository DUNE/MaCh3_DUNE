// #include <iostream>
// #include <iomanip>
// #include <sstream>
// #include <vector>
// #include <memory>

// #include <TH1D.h>
// #include <TH2D.h>
// #include <THStack.h>
// #include <TStyle.h>
// #include <TCanvas.h>
// #include <TRint.h>
// #include <TLegend.h>
// #include <TColor.h>
// #include <TMath.h>
// #include <TFile.h>

// #include "Samples/MaCh3DUNEFactory.h"
// #include "Samples/StructsDUNE.h"
// #include "Fitters/MaCh3Factory.h"


// // ============================================================================
// // Write histograms to ROOT file
// // ============================================================================

// void Write1DHistogramsToFile(
//     const std::string& OutFileName,
//     const std::vector<TH1*>& Histograms)
// {
//     auto OutputFile = std::unique_ptr<TFile>(
//         TFile::Open(OutFileName.c_str(), "RECREATE")
//     );

//     OutputFile->cd();

//     for (auto Hist : Histograms) {
//         Hist->Write();
//     }

//     OutputFile->Close();
// }


// // ============================================================================
// // Write histograms to PDF
// // ============================================================================

// void Write1DHistogramsToPdf(
//     std::string OutFileName,
//     const std::vector<TH1*>& Histograms)
// {
//     // Replace .root with .pdf
//     OutFileName.erase(OutFileName.find('.'));
//     OutFileName += ".pdf";

//     auto c1 = std::unique_ptr<TCanvas>(
//         new TCanvas("c1", "c1", 800, 600)
//     );

//     c1->cd();

//     c1->Print((OutFileName + "[").c_str());

//     for (auto Hist : Histograms) {
//         Hist->Draw("HIST");
//         c1->Print(OutFileName.c_str());
//     }

//     c1->Print((OutFileName + "]").c_str());
// }



// void Print2DHistogramTable(TH2* hist)
// {
//     double totalEvents = hist->Integral();

//     std::string xVarName = hist->GetXaxis()->GetTitle();
//     std::string yVarName = hist->GetYaxis()->GetTitle();

//     if (xVarName.empty()) xVarName = "X Variable";
//     if (yVarName.empty()) yVarName = "Y Variable";

    
//     std::cout << "Histogram: " << hist->GetName() << "\n";
//     std::cout << "Total events: "
//               << std::fixed << std::setprecision(2)
//               << totalEvents << "\n";
//     std::cout << "Only showing template bins with > 0.01% contribution\n";
    

//     std::cout
//         << std::left
//         << std::setw(25) << (xVarName + " Bin Range")
//         << std::setw(25) << (yVarName + " Bin Range")
//         << std::setw(15) << "Events"
//         << std::setw(15) << "% of Total"
//         << "\n";

//     std::cout << std::string(80, '-') << "\n";

//     for (int ix = 1; ix <= hist->GetNbinsX(); ++ix)
//     {
//         double xLow  = hist->GetXaxis()->GetBinLowEdge(ix);
//         double xHigh = hist->GetXaxis()->GetBinUpEdge(ix);

//         for (int iy = 1; iy <= hist->GetNbinsY(); ++iy)
//         {
//             double yLow  = hist->GetYaxis()->GetBinLowEdge(iy);
//             double yHigh = hist->GetYaxis()->GetBinUpEdge(iy);

//             double content = hist->GetBinContent(ix, iy);

//             double percent = 0.0;
//             if (totalEvents > 0.0)
//             {
//                 percent = 100.0 * content / totalEvents;
//             }

//             // Only print bins contributing more than 1%
//             if (percent < 0.01) continue;

//             std::stringstream xRange;
//             xRange << "["
//                    << std::fixed << std::setprecision(2)
//                    << xLow << ", " << xHigh << "]";

//             std::stringstream yRange;
//             yRange << "["
//                    << std::fixed << std::setprecision(2)
//                    << yLow << ", " << yHigh << "]";

//             std::cout
//                 << std::left
//                 << std::setw(25) << xRange.str()
//                 << std::setw(25) << yRange.str()
//                 << std::setw(15) << std::fixed << std::setprecision(2)
//                 << content
//                 << std::setw(15) << std::fixed << std::setprecision(2)
//                 << percent
//                 << "\n";
//         }
//     }

//     std::cout << "========================================================================\n";
// }


// int main(int argc, char* argv[])
// {
//     MaCh3Utils::MaCh3Usage(argc, argv);
//     auto fitMan = MaCh3ManagerFactory(argc, argv);

    
//     ParameterHandlerGeneric* xsec = nullptr;

//     std::vector<SampleHandlerFD*> DUNEPdfs;
//     MakeMaCh3DuneInstance(fitMan, DUNEPdfs, xsec);


    
//     std::vector<TH1*> DUNEHists;

//     for (auto handler : DUNEPdfs) {

//         handler->Reweight();

//         for (int iSample = 0; iSample < handler->GetNsamples(); ++iSample) {

//             TH1* hist = handler->GetMCHist(iSample);

            
//             double totalEvents = hist->Integral(); // True total event count

            
//             TH1* plotHist = (TH1*)hist->Clone(
//                 Form("%s_binWidthNorm", hist->GetName())
//             );

//             plotHist->Scale(1.0, "width");
//             plotHist->GetYaxis()->SetTitle("Events / GeV");

           
//             hist->SetTitle(
//                 Form("%s (Raw Counts)", hist->GetTitle())
//             );

//             plotHist->SetTitle(
//                 Form("%s (Bin Width Normalized)", hist->GetTitle())
//             );

            
//             //DUNEHists.push_back(hist);
//             DUNEHists.push_back(plotHist);

//             std::string EventRateString =
//                 fmt::format("{:.2f}", totalEvents);

//             MACH3LOG_INFO(
//                 "Event rate for {} : {:<5}",
//                 handler->GetSampleTitle(iSample),
//                 EventRateString
//             );

//             handler->PrintIntegral(iSample);

            
//             std::vector<KinematicCut> SelectionVec;

            

//             TH2* Hist2D = handler->Get2DVarHist(
//                 iSample,
//                 "TrueNeutrinoEnergy",   // X variable
//                 "Enubias",   // Y variable
//                 SelectionVec
//             );

//             Print2DHistogramTable(Hist2D);

//             delete Hist2D;
//         }
//     }

//     std::string OutFileName =
//         GetFromManager<std::string>(
//             fitMan->raw()["General"]["OutputFile"],
//             "EventRatesOutput.root"
//         );

//     Write1DHistogramsToFile(OutFileName, DUNEHists);
//     Write1DHistogramsToPdf(OutFileName, DUNEHists);


//     // =========================================================================
//     // Oscillation Mode Breakdown
//     // =========================================================================

//     MACH3LOG_INFO("========================================================================");
//     MACH3LOG_INFO("Oscillation Mode Breakdown:");

//     for (auto handler : DUNEPdfs) {

//         for (int iSample = 0; iSample < handler->GetNsamples(); ++iSample) {

//             MACH3LOG_INFO("======================");

//             int nOscChannels = handler->GetNOscChannels(iSample);

//             for (int iOscChan = 0; iOscChan < nOscChannels; ++iOscChan) {

//                 std::vector<KinematicCut> SelectionVec;

//                 KinematicCut SelecChannel;
//                 SelecChannel.ParamToCutOnIt =
//                     handler->ReturnKinematicParameterFromString(
//                         "OscillationChannel"
//                     );

//                 SelecChannel.LowerBound = iOscChan;
//                 SelecChannel.UpperBound = iOscChan + 1;

//                 SelectionVec.push_back(SelecChannel);

//                 TH1* Hist = handler->Get1DVarHist(
//                     iSample,
//                     handler->GetXBinVarName(iSample),
//                     SelectionVec
//                 );

//                 MACH3LOG_INFO(
//                     "{:<20} : {:<20} : {:<20.2f}",
//                     handler->GetSampleTitle(iSample),
//                     handler->GetFlavourName(iSample, iOscChan),
//                     Hist->Integral()
//                 );
//             }
//         }
//     }


//     // =========================================================================
//     // Interaction Mode Breakdown
//     // =========================================================================

//     MACH3LOG_INFO("========================================================================");
//     MACH3LOG_INFO("Interaction Mode Breakdown:");

//     for (auto handler : DUNEPdfs) {

//         for (int iSample = 0; iSample < handler->GetNsamples(); ++iSample) {

//             MACH3LOG_INFO("======================");

//             MaCh3Modes* Modes = handler->GetMaCh3Modes();
//             int nModeChannels = Modes->GetNModes();

//             for (int iModeChan = 0; iModeChan < nModeChannels; ++iModeChan) {

//                 std::vector<KinematicCut> SelectionVec;

//                 KinematicCut SelecChannel;
//                 SelecChannel.ParamToCutOnIt =
//                     handler->ReturnKinematicParameterFromString("Mode");

//                 SelecChannel.LowerBound = iModeChan;
//                 SelecChannel.UpperBound = iModeChan + 1;

//                 SelectionVec.push_back(SelecChannel);

//                 TH1* Hist = handler->Get1DVarHist(
//                     iSample,
//                     handler->GetXBinVarName(iSample),
//                     SelectionVec
//                 );

//                 MACH3LOG_INFO(
//                     "{:<20} : {:<20} : {:<20.2f}",
//                     handler->GetSampleTitle(iSample),
//                     Modes->GetMaCh3ModeName(iModeChan),
//                     Hist->Integral()
//                 );
//             }
//         }
//     }

//     return 0;
// }


// #include <iostream>
// #include <chrono>
// #include <iomanip>
// #include <vector>
// #include <map>
// #include <algorithm>
// #include <limits>

// #include <TH1D.h>
// #include <THStack.h>
// #include <TStyle.h>
// #include <TCanvas.h>
// #include <TRint.h>
// #include <TLegend.h>
// #include <TColor.h>
// #include <TMath.h>
// #include <TFile.h>
// #include <TTree.h>
// #include <TBranch.h>
// #include <TGraph.h>
// #include <TMultiGraph.h>

// #include "Samples/MaCh3DUNEFactory.h"
// #include "Samples/StructsDUNE.h"
// #include "Fitters/MaCh3Factory.h"

// void WriteOverlayPdf(const std::string& OutFileName,
//                       std::map<std::string, std::map<double, TH1*>>& SysHists,
//                       const std::string& ParamName) {

//   std::string PdfName = OutFileName;
//   PdfName.erase(PdfName.find('.'));
//   PdfName += "_overlay.pdf";

//   gStyle->SetOptStat(0);

//   auto c1 = std::unique_ptr<TCanvas>(new TCanvas("c_overlay", "c_overlay", 900, 700));
//   c1->cd();
//   c1->Print(std::string(PdfName + "[").c_str());

//   for (auto& [title, valMap] : SysHists) {
//     // Clone so we don't mutate the histograms other code might still use
//     std::vector<std::pair<double, TH1*>> sortedVals(valMap.begin(), valMap.end());
//     std::sort(sortedVals.begin(), sortedVals.end(),
//               [](auto& a, auto& b){ return a.first < b.first; });

//     std::vector<TH1*> clones;
//     double maxY = 0;

//     int colours[] = {kBlack, kRed+1, kAzure+2, kGreen+2, kMagenta+1};
//     int ic = 0;

//     auto leg = std::make_unique<TLegend>(0.65, 0.7, 0.88, 0.88);
//     leg->SetBorderSize(0);
//     leg->SetFillStyle(0);

//     for (auto& [val, hOrig] : sortedVals) {
//       TH1* h = (TH1*)hOrig->Clone(Form("%s_clone_%.3f", hOrig->GetName(), val));
//       h->SetDirectory(nullptr);

//       // Divide out bin widths: content -> content / bin_width
//       h->Scale(1.0, "width");

//       h->SetLineColor(colours[ic % 5]);
//       h->SetLineWidth(2);
//       h->SetMarkerColor(colours[ic % 5]);
//       h->SetMarkerStyle(20);
//       h->SetTitle(Form("%s;%s;Events / bin width", title.c_str(), h->GetXaxis()->GetTitle()));

//       maxY = std::max(maxY, h->GetMaximum());
//       leg->AddEntry(h, Form("%s = %.2f", ParamName.c_str(), val), "l");
//       clones.push_back(h);
//       ic++;
//     }

//     for (size_t i = 0; i < clones.size(); ++i) {
//       clones[i]->SetMaximum(maxY * 1.3);
//       clones[i]->Draw(i == 0 ? "HIST" : "HIST SAME");
//     }
//     leg->Draw();
//     c1->Print(PdfName.c_str());

//     for (auto h : clones) delete h;
//   }

//   c1->Print(std::string(PdfName + "]").c_str());
//   MACH3LOG_INFO("Wrote overlay PDF: {}", PdfName);
// }

// void Write1DHistogramsToFile(std::string OutFileName, std::vector<TH1*> Histograms) {
//   auto OutputFile = std::unique_ptr<TFile>(TFile::Open(OutFileName.c_str(), "RECREATE"));
//   OutputFile->cd();
//   for(auto Hist : Histograms){
//     Hist->Write();
//   }
//   OutputFile->Close();
// }

// void Write1DHistogramsToPdf(std::string OutFileName, std::vector<TH1*> Histograms) {
//   //Remove root from end of file
//   OutFileName.erase(OutFileName.find('.'));
//   OutFileName+=".pdf";

//   auto c1 = std::unique_ptr<TCanvas>(new TCanvas("c1", "c1", 800, 600));
//   c1->cd();
//   c1->Print(std::string(OutFileName+"[").c_str());
//   for(auto Hist : Histograms){
//     Hist->Draw("HIST COLZ");
//     c1->Print(OutFileName.c_str());
//   }
//   c1->Print(std::string(OutFileName+"]").c_str());
// }

// //###############################################################################################################################
// // Posterior-chain event-rate / likelihood check
// //
// // Reads a standard MaCh3 MCMC output tree (one branch per parameter, plus a
// // LogL branch), and for a user-specified list of steps: pushes those
// // parameter values into the systematic handler, reweights every sample, and
// // records the predicted event rate per sample together with both the LogL
// // stored in the chain at that step and a freshly recomputed LogL from the
// // reweighted samples. Results are written out as a multi-page PDF.
// //###############################################################################################################################

// struct ChainStepResult {
//   Long64_t Step = 0;
//   double LogL_Chain = std::numeric_limits<double>::quiet_NaN();
//   double LogL_Recomputed = std::numeric_limits<double>::quiet_NaN();
//   std::map<std::string, double> SampleRates; // sample title -> predicted event rate
//   double TotalRate = 0.0;
// };

// // Try to find, for every parameter known to `xsec`, which branch in the
// // posterior chain corresponds to it.
// //
// // MaCh3 has used a couple of different naming conventions for these branches
// // across versions ("xsec_<i>" in older releases, "param_<i>" in newer ones),
// // and some setups additionally store the human-readable parameter name in
// // the branch *title* rather than the branch name. This function tries the
// // obvious candidates and falls back to matching on branch title, so in most
// // cases you shouldn't need to edit it — but if none of your parameters are
// // found, open the chain in `root -l` and check `Chain->Print()` to see what
// // your branches are actually called, then add the pattern below.
// std::map<int, std::string> MapParamsToChainBranches(TTree* Chain, ParameterHandlerGeneric* xsec) {
//   std::map<int, std::string> ParToBranch;

//   std::vector<std::string> BranchNames;
//   std::vector<std::string> BranchTitles;
//   for (auto* obj : *Chain->GetListOfBranches()) {
//     auto* br = (TBranch*)obj;
//     BranchNames.push_back(br->GetName());
//     BranchTitles.push_back(br->GetTitle());
//   }

//   int NPars = xsec->GetNumParams(); // <-- confirm this matches the real API
//   for (int i = 0; i < NPars; ++i) {
//     std::string ParName = xsec->GetParName(i); // <-- confirm this matches the real API

//     std::vector<std::string> Candidates = {
//       Form("param_%d", i),
//       Form("xsec_%d", i),
//       ParName,
//       Form("param_%s", ParName.c_str()),
//       Form("xsec_%s", ParName.c_str())
//     };

//     bool Found = false;
//     for (auto& cand : Candidates) {
//       if (std::find(BranchNames.begin(), BranchNames.end(), cand) != BranchNames.end()) {
//         ParToBranch[i] = cand;
//         Found = true;
//         break;
//       }
//     }

//     // Fall back to matching on branch title, which MaCh3 often sets to the
//     // human-readable parameter name even when the branch itself is indexed
//     if (!Found) {
//       for (size_t b = 0; b < BranchTitles.size(); ++b) {
//         if (BranchTitles[b] == ParName) {
//           ParToBranch[i] = BranchNames[b];
//           Found = true;
//           break;
//         }
//       }
//     }

//     if (!Found) {
//       MACH3LOG_WARN("Could not find a chain branch for parameter {} ({}) "
//                      "— it will be left at its current value for every step",
//                      i, ParName);
//     }
//   }

//   return ParToBranch;
// }

// void RunPosteriorChainCheck(const std::string& ChainFileName,
//                              const std::string& ChainTreeName,
//                              const std::string& LogLBranchName,
//                              const std::vector<Long64_t>& StepsToProcess,
//                              ParameterHandlerGeneric* xsec,
//                              std::vector<SampleHandlerFD*>& DUNEPdfs,
//                              const std::string& OutPdfName) {

//   MACH3LOG_INFO("========================================================================");
//   MACH3LOG_INFO("Posterior-chain event-rate / likelihood check");
//   MACH3LOG_INFO("Chain file: {}, tree: {}", ChainFileName, ChainTreeName);

//   auto ChainFile = std::unique_ptr<TFile>(TFile::Open(ChainFileName.c_str(), "READ"));
//   if (!ChainFile || ChainFile->IsZombie()) {
//     MACH3LOG_ERROR("Could not open chain file {}", ChainFileName);
//     return;
//   }
//   auto* Chain = (TTree*)ChainFile->Get(ChainTreeName.c_str());
//   if (!Chain) {
//     MACH3LOG_ERROR("Could not find tree {} in {}", ChainTreeName, ChainFileName);
//     return;
//   }

//   auto ParToBranch = MapParamsToChainBranches(Chain, xsec);
//   if (ParToBranch.empty()) {
//     MACH3LOG_ERROR("No parameter branches were matched in the chain, aborting check");
//     return;
//   }

//   int NPars = xsec->GetNumParams(); // <-- confirm this matches the real API
//   std::vector<double> ParVals(NPars, 0.0);
//   for (auto& [i, name] : ParToBranch) {
//     Chain->SetBranchAddress(name.c_str(), &ParVals[i]);
//   }

//   double LogL_Chain = std::numeric_limits<double>::quiet_NaN();
//   bool HaveLogLBranch = (Chain->GetBranch(LogLBranchName.c_str()) != nullptr);
//   if (HaveLogLBranch) {
//     Chain->SetBranchAddress(LogLBranchName.c_str(), &LogL_Chain);
//   } else {
//     MACH3LOG_WARN("Could not find LogL branch '{}' in chain, "
//                    "stored LogL will be reported as NaN — check the branch name "
//                    "(e.g. with root -l {} -e '_file0->Get(\"{}\")->Print()')",
//                    LogLBranchName, ChainFileName, ChainTreeName);
//   }

//   Long64_t NEntries = Chain->GetEntries();
//   std::vector<ChainStepResult> Results;

//   for (Long64_t step : StepsToProcess) {
//     if (step < 0 || step >= NEntries) {
//       MACH3LOG_WARN("Requested step {} is outside the chain (0..{}) — skipping", step, NEntries - 1);
//       continue;
//     }
//     Chain->GetEntry(step);

//     // Push the parameter values from this step into the systematic handler
//     for (auto& [i, name] : ParToBranch) {
//       xsec->SetParProp(i, ParVals[i]); // <-- confirm this matches the real API
//     }

//     ChainStepResult Res;
//     Res.Step = step;
//     Res.LogL_Chain = LogL_Chain;

//     double TotalRate = 0.0;
//     double RecomputedLogL = 0.0;

//     for (auto handler : DUNEPdfs) {
//       handler->Reweight();
//       RecomputedLogL += handler->GetLikelihood(); // <-- confirm this matches the real API
//       for (int iSample = 0; iSample < handler->GetNsamples(); iSample++) {
//         std::string title = handler->GetSampleTitle(iSample);
//         double Rate = handler->GetMCHist(iSample)->Integral();
//         Res.SampleRates[title] += Rate;
//         TotalRate += Rate;
//       }
//     }
//     // Add the systematic/prior penalty term so the recomputed LogL is on the
//     // same footing as what's stored in the chain (typically sample
//     // likelihoods + systematic priors, summed). If your fit has other
//     // ParameterHandlerBase objects beyond `xsec` (e.g. a separate
//     // oscillation-parameter covariance), add their GetLikelihood() here too,
//     // or the comparison to LogL_Chain won't fully close.
//     RecomputedLogL += xsec->GetLikelihood(); // <-- confirm this matches the real API

//     Res.TotalRate = TotalRate;
//     Res.LogL_Recomputed = RecomputedLogL;

//     MACH3LOG_INFO("------------------------------------------------------------------------");
//     MACH3LOG_INFO("Step {:<8} : TotalRate={:.3f}  LogL(chain)={:.3f}  LogL(recomputed)={:.3f}  diff={:.3f}",
//                   Res.Step, Res.TotalRate, Res.LogL_Chain, Res.LogL_Recomputed,
//                   Res.LogL_Chain - Res.LogL_Recomputed);
//     for (auto& [title, rate] : Res.SampleRates) {
//       MACH3LOG_INFO("    {:<25} : {:.3f}", title, rate);
//     }

//     Results.push_back(Res);
//   }

//   if (Results.empty()) {
//     MACH3LOG_ERROR("No valid steps were processed, not writing PDF");
//     return;
//   }

//   //###################################
//   // Make the PDF
//   //###################################
//   gStyle->SetOptStat(0);
//   auto c1 = std::unique_ptr<TCanvas>(new TCanvas("c_chain", "c_chain", 900, 700));
//   c1->cd();
//   c1->Print(std::string(OutPdfName + "[").c_str());

//   std::vector<std::string> SampleTitles;
//   for (auto& [title, rate] : Results.front().SampleRates) SampleTitles.push_back(title);

//   int NSteps = (int)Results.size();

//   // --- Page 1: total predicted event rate vs step ---
//   {
//     TGraph g(NSteps);
//     for (int i = 0; i < NSteps; ++i) g.SetPoint(i, (double)Results[i].Step, Results[i].TotalRate);
//     g.SetTitle("Total predicted event rate;Chain step;Events");
//     g.SetLineColor(kBlack);
//     g.SetLineWidth(2);
//     g.SetMarkerStyle(20);
//     g.SetMarkerColor(kBlack);
//     g.Draw("ALP");
//     c1->Print(OutPdfName.c_str());
//   }

//   // --- Page 2: per-sample predicted event rate vs step ---
//   {
//     TMultiGraph mg;
//     auto leg = std::make_unique<TLegend>(0.65, 0.6, 0.88, 0.88);
//     leg->SetBorderSize(0);
//     leg->SetFillStyle(0);
//     int colours[] = {kBlack, kRed+1, kAzure+2, kGreen+2, kMagenta+1, kOrange+1};
//     int ic = 0;
//     std::vector<std::unique_ptr<TGraph>> Graphs;
//     for (auto& title : SampleTitles) {
//       auto g = std::make_unique<TGraph>(NSteps);
//       for (int i = 0; i < NSteps; ++i) g->SetPoint(i, (double)Results[i].Step, Results[i].SampleRates.at(title));
//       g->SetLineColor(colours[ic % 6]);
//       g->SetMarkerColor(colours[ic % 6]);
//       g->SetMarkerStyle(20);
//       g->SetLineWidth(2);
//       leg->AddEntry(g.get(), title.c_str(), "l");
//       mg.Add(g.get(), "LP");
//       Graphs.push_back(std::move(g));
//       ic++;
//     }
//     mg.SetTitle("Predicted event rate per sample;Chain step;Events");
//     mg.Draw("A");
//     leg->Draw();
//     c1->Print(OutPdfName.c_str());
//   }

//   // --- Page 3: LogL stored in chain vs recomputed LogL ---
//   {
//     TGraph gChain(NSteps), gRecomp(NSteps);
//     for (int i = 0; i < NSteps; ++i) {
//       gChain.SetPoint(i, (double)Results[i].Step, Results[i].LogL_Chain);
//       gRecomp.SetPoint(i, (double)Results[i].Step, Results[i].LogL_Recomputed);
//     }
//     gChain.SetLineColor(kBlack);   gChain.SetMarkerColor(kBlack);   gChain.SetMarkerStyle(20); gChain.SetLineWidth(2);
//     gRecomp.SetLineColor(kRed+1);  gRecomp.SetMarkerColor(kRed+1);  gRecomp.SetMarkerStyle(21); gRecomp.SetLineWidth(2);

//     TMultiGraph mg;
//     mg.Add(&gChain, "LP");
//     mg.Add(&gRecomp, "LP");
//     mg.SetTitle("-log L at each step;Chain step;-log L");
//     mg.Draw("A");

//     auto leg = std::make_unique<TLegend>(0.6, 0.75, 0.88, 0.88);
//     leg->SetBorderSize(0);
//     leg->SetFillStyle(0);
//     leg->AddEntry(&gChain, "Stored in chain", "l");
//     leg->AddEntry(&gRecomp, "Recomputed after reweight", "l");
//     leg->Draw();
//     c1->Print(OutPdfName.c_str());
//   }

//   // --- Page 4: difference between stored and recomputed LogL ---
//   {
//     TGraph g(NSteps);
//     for (int i = 0; i < NSteps; ++i)
//       g.SetPoint(i, (double)Results[i].Step, Results[i].LogL_Chain - Results[i].LogL_Recomputed);
//     g.SetTitle("LogL(chain) - LogL(recomputed);Chain step;#Delta(-log L)");
//     g.SetLineColor(kAzure+2);
//     g.SetLineWidth(2);
//     g.SetMarkerColor(kAzure+2);
//     g.SetMarkerStyle(20);
//     g.Draw("ALP");
//     c1->Print(OutPdfName.c_str());
//   }

//   c1->Print(std::string(OutPdfName + "]").c_str());
//   MACH3LOG_INFO("Wrote posterior-chain check PDF: {}", OutPdfName);
// }

// int main(int argc, char * argv[]) {
//   MaCh3Utils::MaCh3Usage(argc, argv);
//   auto fitMan = MaCh3ManagerFactory(argc, argv);

//   //###############################################################################################################################
//   //Create SampleHandlerFD objects
  
//   ParameterHandlerGeneric* xsec = nullptr;
  
//   std::vector<SampleHandlerFD*> DUNEPdfs;
//   MakeMaCh3DuneInstance(fitMan, DUNEPdfs, xsec);

//   //###############################################################################################################################
//   //Perform reweight and print total integral

//   std::vector<TH1*> DUNEHists;
//   for(auto handler : DUNEPdfs){
//     handler->Reweight();
//     for (int iSample=0; iSample<handler->GetNsamples(); iSample++) {
//       DUNEHists.push_back(handler->GetMCHist(iSample));

//       std::string EventRateString = fmt::format("{:.2f}", handler->GetMCHist(iSample)->Integral());
//       MACH3LOG_INFO("Event rate for {} : {:<5}", handler->GetSampleTitle(iSample), EventRateString);
//       handler->PrintIntegral(iSample);
//     }
//   }

//   std::string OutFileName = GetFromManager<std::string>(fitMan->raw()["General"]["OutputFile"], "EventRatesOutput.root");
//   Write1DHistogramsToFile(OutFileName, DUNEHists); 
//   Write1DHistogramsToPdf(OutFileName, DUNEHists);


//    //###############################################################################################################################
//   // MissingProtonFD systematic check: set to 0.0 and 0.2, reweight, compare
//   MACH3LOG_INFO("========================================================================");
//   MACH3LOG_INFO("Testing MissingProtonFD functional parameter");

//   const std::string ParamName = "MissingProtonFD";
//   int ParIndex = xsec->GetParIndex(ParamName); // <-- confirm this matches the real API

//   if (ParIndex < 0) {
//     MACH3LOG_ERROR("Could not find parameter {} in ParameterHandler, check it's defined in your systematics YAML", ParamName);
//   } else {
//     std::vector<double> TestValues = {0.0, 0.2};
//     // sample-title -> {value -> histogram}
//     std::map<std::string, std::map<double, TH1*>> SysHists;

//     for (double val : TestValues) {
//       MACH3LOG_INFO("------------------------------------------------------------------------");
//       MACH3LOG_INFO("Setting {} = {}", ParamName, val);
//       xsec->SetParProp(ParIndex, val); // <-- confirm this matches the real API

//       for (auto handler : DUNEPdfs) {
//         handler->Reweight();
//         for (int iSample=0; iSample<handler->GetNsamples(); iSample++) {
//           std::string title = handler->GetSampleTitle(iSample);
//           TH1* Hist = (TH1*)handler->GetMCHist(iSample)->Clone(
//               Form("%s_%s_%.2f", title.c_str(), ParamName.c_str(), val));
//           Hist->SetDirectory(nullptr); // detach from any open TFile so it survives
//           SysHists[title][val] = Hist;

//           MACH3LOG_INFO("[{} = {:.2f}] Event rate for {:<20} : {:.2f}",
//                         ParamName, val, title, Hist->Integral());
//         }
//       }
//     }

//     // Print per-bin differences so you can see the shape effect, not just the integral
//     MACH3LOG_INFO("------------------------------------------------------------------------");
//     MACH3LOG_INFO("Bin-by-bin comparison ({} = 0.0 vs {} = 0.2):", ParamName, ParamName);
//     for (auto& [title, valMap] : SysHists) {
//       TH1* h0 = valMap.at(0.0);
//       TH1* h2 = valMap.at(0.2);
//       MACH3LOG_INFO("Sample: {}", title);
//       for (int b = 1; b <= h0->GetNbinsX(); ++b) {
//         double c0 = h0->GetBinContent(b);
//         double c2 = h2->GetBinContent(b);
//         double diff = c2 - c0;
//         MACH3LOG_INFO("  bin {:<3} : nominal(0.0)={:.3f}  shifted(0.2)={:.3f}  diff={:.3f}",
//                       b, c0, c2, diff);
//       }
//       if (std::abs(h0->Integral() - h2->Integral()) < 1e-9 &&
//           [&](){ for(int b=1;b<=h0->GetNbinsX();++b) if(std::abs(h0->GetBinContent(b)-h2->GetBinContent(b))>1e-9) return false; return true; }()) {
//         MACH3LOG_WARN("Histograms for {} are IDENTICAL between par=0.0 and par=0.2 — systematic is having NO effect!", title);
//       }
//     }

//     // Dump both sets of histograms to file/pdf for visual inspection
//     std::vector<TH1*> FlatHists;
//     for (auto& [title, valMap] : SysHists)
//       for (auto& [val, h] : valMap)
//         FlatHists.push_back(h);

//     Write1DHistogramsToFile("MissingProtonFD_test.root", FlatHists);
//     Write1DHistogramsToPdf("MissingProtonFD_test.root", FlatHists);
//     WriteOverlayPdf("MissingProtonFD_test.root", SysHists, ParamName);
//     // Reset parameter back to nominal before continuing with the rest of the app
//     xsec->SetParProp(ParIndex, 0.0);
//     for (auto handler : DUNEPdfs) handler->Reweight();
//   }

//   //###############################################################################################################################
//   // Posterior-chain event-rate / likelihood check
//   //
//   // Reads a posterior chain, and for a hand-picked list of MCMC steps,
//   // reweights the samples with the parameter values at that step and
//   // compares the resulting event rate and LogL to what's stored in the chain.
//   //
//   // Set General:PosteriorChain:FileName in your config YAML to point at the
//   // chain, e.g.:
//   //   General:
//   //     PosteriorChain:
//   //       FileName: "/path/to/mcmc_chain.root"
//   //       TreeName: "posteriors"     # optional, defaults shown
//   //       LogLBranch: "LogL"         # optional, defaults shown
//   //       OutputPdf: "PosteriorChainCheck.pdf"  # optional

//   std::string ChainFileName = GetFromManager<std::string>(fitMan->raw()["General"]["PosteriorChain"]["FileName"], "");
//   std::string ChainTreeName = GetFromManager<std::string>(fitMan->raw()["General"]["PosteriorChain"]["TreeName"], "posteriors");
//   std::string LogLBranchName = GetFromManager<std::string>(fitMan->raw()["General"]["PosteriorChain"]["LogLBranch"], "LogL");
//   std::string ChainPdfName = GetFromManager<std::string>(fitMan->raw()["General"]["PosteriorChain"]["OutputPdf"], "PosteriorChainCheck.pdf");

//   // <-- Edit this to whichever chain steps you want to check
//   std::vector<Long64_t> StepsToProcess = {0, 1000, 5000, 10000, 20000};

//   if (!ChainFileName.empty()) {
//     RunPosteriorChainCheck(ChainFileName, ChainTreeName, LogLBranchName, StepsToProcess, xsec, DUNEPdfs, ChainPdfName);
//     // Reset parameters back to nominal before continuing with the rest of the
//     // app (mirrors the reset done above after the MissingProtonFD check —
//     // adjust if your parameters' nominal values aren't all 0)
//     for (int i = 0; i < xsec->GetNumParams(); ++i) xsec->SetParProp(i, 0.0); // <-- confirm this matches the real API
//     for (auto handler : DUNEPdfs) handler->Reweight();
//   } else {
//     MACH3LOG_INFO("No General:PosteriorChain:FileName set in config, skipping posterior chain check");
//   }

//   //###############################################################################################################################
//   //Make oscillation channel breakdown

//   MACH3LOG_INFO("========================================================================");
//   MACH3LOG_INFO("========================================================================");
//   MACH3LOG_INFO("Oscillation Mode Breakdown:");
  
//   for(auto handler : DUNEPdfs) {
//     for (int iSample = 0; iSample < handler->GetNsamples(); iSample++) {
//       MACH3LOG_INFO("======================");
//       int nOscChannels = handler->GetNOscChannels(iSample);
//       for (int iOscChan=0;iOscChan<nOscChannels;iOscChan++) {
//         std::vector< KinematicCut > SelectionVec;

//         KinematicCut SelecChannel;
//         SelecChannel.ParamToCutOnIt = handler->ReturnKinematicParameterFromString("OscillationChannel");
//         SelecChannel.LowerBound = iOscChan;
//         SelecChannel.UpperBound = iOscChan+1;
//         SelectionVec.push_back(SelecChannel);
        
//         TH1* Hist = handler->Get1DVarHist(iSample, handler->GetXBinVarName(iSample),SelectionVec);
//         MACH3LOG_INFO("{:<20} : {:<20} : {:<20.2f}",handler->GetSampleTitle(iSample),handler->GetFlavourName(iSample, iOscChan),Hist->Integral());
//       }

//       TH1* Hist = handler->Get1DVarHist(iSample, handler->GetXBinVarName(iSample));
//       MACH3LOG_INFO("{:<20} : {:<20.2f}",handler->GetSampleTitle(iSample),Hist->Integral());
//     }
//   }

//   //###############################################################################################################################
//   //Make interaction channel breakdown

//   MACH3LOG_INFO("========================================================================");
//   MACH3LOG_INFO("========================================================================");
//   MACH3LOG_INFO("Interaction Mode Breakdown:");

//   for(auto handler : DUNEPdfs) {
//     for (int iSample = 0; iSample < handler->GetNsamples(); iSample++) {
//       MACH3LOG_INFO("======================");

//       MaCh3Modes* Modes = handler->GetMaCh3Modes();
//       int nModeChannels = Modes->GetNModes();
//       for (int iModeChan=0;iModeChan<nModeChannels;iModeChan++) {
//         std::vector< KinematicCut > SelectionVec;

//         KinematicCut SelecChannel;
//         SelecChannel.ParamToCutOnIt = handler->ReturnKinematicParameterFromString("Mode");
//         SelecChannel.LowerBound = iModeChan;
//         SelecChannel.UpperBound = iModeChan+1;
//         SelectionVec.push_back(SelecChannel);

//         TH1* Hist = handler->Get1DVarHist(iSample, handler->GetXBinVarName(iSample),SelectionVec);
//         MACH3LOG_INFO("{:<20} : {:<20} : {:<20.2f}",handler->GetSampleTitle(iSample),Modes->GetMaCh3ModeName(iModeChan),Hist->Integral());
//       }

//       TH1* Hist = handler->Get1DVarHist(iSample, handler->GetXBinVarName(iSample));
//       MACH3LOG_INFO("{:<20} : {:<20.2f}",handler->GetSampleTitle(iSample),Hist->Integral());
//     }
//   }

//   //###############################################################################################################################
// }


#include <iostream>
#include <chrono>
#include <iomanip>
#include <vector>
#include <cmath>

#include <TH1D.h>
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TRint.h>
#include <TLegend.h>
#include <TColor.h>
#include <TMath.h>
#include <TLatex.h>
#include <TLine.h>

#include "Samples/MaCh3DUNEFactory.h"
#include "Samples/StructsDUNE.h"
#include "Fitters/MaCh3Factory.h"

//=====================================================================================
// Global style: call once, before making any plots. Gives you a clean, uncluttered
// look that reproduces well in print (no default ROOT grey stats box, larger fonts,
// sensible margins, thin serif-free axis labels).
//=====================================================================================
void SetThesisStyle() {
  TStyle* thesisStyle = new TStyle("Thesis", "Thesis style");

  // Canvas / pad
  thesisStyle->SetCanvasColor(kWhite);
  thesisStyle->SetPadColor(kWhite);
  thesisStyle->SetFrameBorderMode(0);
  thesisStyle->SetCanvasBorderMode(0);
  thesisStyle->SetPadBorderMode(0);
  thesisStyle->SetPadTickX(1);
  thesisStyle->SetPadTickY(1);

  // Stats / title boxes off
  thesisStyle->SetOptStat(0);
  thesisStyle->SetOptTitle(0); // we'll draw our own title with TLatex for control over placement

  // Fonts: 42 = Helvetica, precision 2 (scales with pad, not fixed pixel size)
  const int font = 42;
  thesisStyle->SetTextFont(font);
  thesisStyle->SetLabelFont(font, "XYZ");
  thesisStyle->SetTitleFont(font, "XYZ");
  thesisStyle->SetLegendFont(font);
  thesisStyle->SetLabelSize(0.045, "XYZ");
  thesisStyle->SetTitleSize(0.050, "XYZ");
  thesisStyle->SetTitleOffset(1.05, "X");
  thesisStyle->SetTitleOffset(1.25, "Y");

  // Margins
  thesisStyle->SetPadTopMargin(0.08);
  thesisStyle->SetPadBottomMargin(0.14);
  thesisStyle->SetPadLeftMargin(0.14);
  thesisStyle->SetPadRightMargin(0.05);

  // Line/marker defaults
  thesisStyle->SetHistLineWidth(2);
  thesisStyle->SetEndErrorSize(0);

  thesisStyle->cd();
  gROOT->ForceStyle();
}

//=====================================================================================
// Sanitise a sample title into something safe to use as a filename: spaces and
// anything that isn't alphanumeric/underscore/hyphen become underscores.
//=====================================================================================
std::string SanitiseFileName(const std::string& in) {
  std::string out = in;
  for (auto& c : out) {
    if (!std::isalnum(static_cast<unsigned char>(c)) && c != '_' && c != '-') c = '_';
  }
  return out;
}

//=====================================================================================
// Overlay plot for a systematic scan: main panel with each parameter value overlaid,
// plus a ratio panel underneath (each variation / nominal, with stat error bars) so
// shape effects that are invisible in the integral are obvious at a glance.
//
// Writes ONE PDF (+ one PNG) PER SAMPLE, named <OutFileName>_overlay_<sample>.pdf,
// rather than a single multi-page file -- makes it much easier to grab just the one
// plot you want for a given chapter/section without paging through a combined PDF.
//
// Assumes the FIRST entry (sorted by value) is nominal for the ratio panel.
//=====================================================================================
void WriteOverlayPdf(const std::string& OutFileName,
                      std::map<std::string, std::map<double, TH1*>>& SysHists,
                      const std::string& ParamName) {

  std::string BaseName = OutFileName;
  BaseName.erase(BaseName.find('.'));

  const int colours[] = {kBlack, kRed + 1, kAzure + 2, kGreen + 2, kMagenta + 1};

  for (auto& [title, valMap] : SysHists) {
    std::string PdfName = BaseName + "_overlay_" + SanitiseFileName(title) + ".pdf";
    std::string PngName = BaseName + "_overlay_" + SanitiseFileName(title) + ".png";

    auto c1 = std::unique_ptr<TCanvas>(new TCanvas(
        Form("c_overlay_%s", SanitiseFileName(title).c_str()), "c_overlay", 1000, 800));

    // --- Split canvas: main plot on top (70%), ratio panel below (30%) ---
    TPad* padMain = new TPad("padMain", "padMain", 0.0, 0.30, 1.0, 1.0);
    TPad* padRatio = new TPad("padRatio", "padRatio", 0.0, 0.0, 1.0, 0.30);
    padMain->SetBottomMargin(0.02);   // touches the ratio pad
    padMain->SetTopMargin(0.08);
    padMain->SetLeftMargin(0.14);
    padMain->SetRightMargin(0.05);
    padRatio->SetTopMargin(0.03);
    padRatio->SetBottomMargin(0.35);  // room for the x-axis title, shared with main
    padRatio->SetLeftMargin(0.14);
    padRatio->SetRightMargin(0.05);
    padMain->Draw();
    padRatio->Draw();

    std::vector<std::pair<double, TH1*>> sortedVals(valMap.begin(), valMap.end());
    std::sort(sortedVals.begin(), sortedVals.end(),
              [](auto& a, auto& b) { return a.first < b.first; });

    std::vector<TH1*> clones;
    std::vector<TH1*> ratios;
    double maxY = 0;
    int ic = 0;

    padMain->cd();
    auto leg = std::make_unique<TLegend>(0.62, 0.68, 0.90, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.045);

    TH1* nominal = nullptr;

    for (auto& [val, hOrig] : sortedVals) {
      TH1* h = (TH1*)hOrig->Clone(Form("%s_clone_%.3f", hOrig->GetName(), val));
      h->SetDirectory(nullptr);

      // Divide out bin widths so differently-binned histograms are comparable
      h->Scale(1.0, "width");

      h->SetLineColor(colours[ic % 5]);
      h->SetLineWidth(3);
      h->SetMarkerColor(colours[ic % 5]);
      h->SetMarkerStyle(20);
      h->SetMarkerSize(0.6); // small marker so the stat error bar caps are visible without cluttering the line
      h->GetXaxis()->SetLabelSize(0); // hidden on main pad; ratio pad shows the axis
      h->GetYaxis()->SetTitle("Events / bin width");
      h->GetYaxis()->SetTitleSize(0.06);
      h->GetYaxis()->SetLabelSize(0.05);
      h->GetYaxis()->SetTitleOffset(1.05);

      maxY = std::max(maxY, h->GetMaximum());
      leg->AddEntry(h, Form("%s = %.2f", ParamName.c_str(), val), "lep");
      clones.push_back(h);

      if (ic == 0) nominal = h; // first (smallest) value treated as nominal/reference
      ic++;
    }

    for (size_t i = 0; i < clones.size(); ++i) {
      clones[i]->SetMaximum(maxY * 1.35);
      clones[i]->SetMinimum(0.0);
      // "HIST" draws the connecting line, "E1" overlays stat error bars with caps
      clones[i]->Draw(i == 0 ? "HIST E1" : "HIST E1 SAME");
    }
    leg->Draw();

    // Title drawn by hand so we control font/size consistently with the rest of the style
    TLatex titleLat;
    titleLat.SetNDC();
    titleLat.SetTextFont(42);
    titleLat.SetTextSize(0.055);
    titleLat.DrawLatex(0.14, 0.94, title.c_str());

    // --- Ratio panel ---
    padRatio->cd();
    ic = 0;
    for (auto h : clones) {
      TH1* r = (TH1*)h->Clone(Form("%s_ratio", h->GetName()));
      // Divide() propagates stat errors from both h and nominal into r's bin errors
      r->Divide(nominal);
      r->SetLineColor(colours[ic % 5]);
      r->SetLineWidth(3);
      r->SetMarkerColor(colours[ic % 5]);
      r->SetMarkerStyle(20);
      r->SetMarkerSize(0.6);
      r->GetYaxis()->SetTitle(Form("Ratio to %s=%.2f", ParamName.c_str(), sortedVals.front().first));
      r->GetYaxis()->SetTitleSize(0.11);
      r->GetYaxis()->SetLabelSize(0.10);
      r->GetYaxis()->SetTitleOffset(0.55);
      r->GetYaxis()->SetNdivisions(505);
      r->GetXaxis()->SetTitle(nominal->GetXaxis()->GetTitle());
      r->GetXaxis()->SetTitleSize(0.13);
      r->GetXaxis()->SetLabelSize(0.11);
      r->SetMinimum(0.80);
      r->SetMaximum(1.20);
      // "HIST E1" so the ratio line still reads cleanly at a glance, with error bar caps on top
      r->Draw(ic == 0 ? "HIST E1" : "HIST E1 SAME");
      ratios.push_back(r);
      ic++;
    }
    TLine unity(nominal->GetXaxis()->GetXmin(), 1.0, nominal->GetXaxis()->GetXmax(), 1.0);
    unity.SetLineStyle(2);
    unity.SetLineColor(kGray + 2);
    unity.Draw("SAME");

    c1->cd();
    c1->Print(PdfName.c_str());
    c1->Print(PngName.c_str()); // print-resolution PNG, for LaTeX \includegraphics

    for (auto h : clones) delete h;
    for (auto r : ratios) delete r;
    delete padMain;
    delete padRatio;

    MACH3LOG_INFO("Wrote overlay plot for {}: {}", title, PdfName);
  }
}

void Write1DHistogramsToFile(std::string OutFileName, std::vector<TH1*> Histograms) {
  auto OutputFile = std::unique_ptr<TFile>(TFile::Open(OutFileName.c_str(), "RECREATE"));
  OutputFile->cd();
  for (auto Hist : Histograms) {
    Hist->Write();
  }
  OutputFile->Close();
}

void Write1DHistogramsToPdf(std::string OutFileName, std::vector<TH1*> Histograms) {
  OutFileName.erase(OutFileName.find('.'));
  OutFileName += ".pdf";

  auto c1 = std::unique_ptr<TCanvas>(new TCanvas("c1", "c1", 800, 600));
  c1->cd();
  c1->Print(std::string(OutFileName + "[").c_str());
  for (auto Hist : Histograms) {
    Hist->Draw("HIST COLZ");
    c1->Print(OutFileName.c_str());
  }
  c1->Print(std::string(OutFileName + "]").c_str());
}

int main(int argc, char* argv[]) {
  MaCh3Utils::MaCh3Usage(argc, argv);
  auto fitMan = MaCh3ManagerFactory(argc, argv);

  SetThesisStyle();

  //###############################################################################################################################
  //Create SampleHandlerFD objects
  ParameterHandlerGeneric* xsec = nullptr;
  std::vector<SampleHandlerFD*> DUNEPdfs;
  MakeMaCh3DuneInstance(fitMan, DUNEPdfs, xsec);

  //###############################################################################################################################
  //Perform reweight and print total integral

  std::vector<TH1*> DUNEHists;
  for (auto handler : DUNEPdfs) {
    handler->Reweight();
    for (int iSample = 0; iSample < handler->GetNsamples(); iSample++) {
      DUNEHists.push_back(handler->GetMCHist(iSample));

      std::string EventRateString = fmt::format("{:.2f}", handler->GetMCHist(iSample)->Integral());
      MACH3LOG_INFO("Event rate for {} : {:<5}", handler->GetSampleTitle(iSample), EventRateString);
      handler->PrintIntegral(iSample);
    }
  }

  std::string OutFileName = GetFromManager<std::string>(fitMan->raw()["General"]["OutputFile"], "EventRatesOutput.root");
  Write1DHistogramsToFile(OutFileName, DUNEHists);
  Write1DHistogramsToPdf(OutFileName, DUNEHists);

  //###############################################################################################################################
  // MissingProtonFD systematic check: set to 0.0 and 0.2, reweight, compare
  MACH3LOG_INFO("========================================================================");
  MACH3LOG_INFO("Testing MissingProtonFD functional parameter");

  const std::string ParamName = "MissingProtonFD";
  int ParIndex = xsec->GetParIndex(ParamName); // <-- confirm this matches the real API

  if (ParIndex < 0) {
    MACH3LOG_ERROR("Could not find parameter {} in ParameterHandler, check it's defined in your systematics YAML", ParamName);
  } else {
    std::vector<double> TestValues = {0.0, 0.2};
    // sample-title -> {value -> histogram}
    std::map<std::string, std::map<double, TH1*>> SysHists;

    for (double val : TestValues) {
      MACH3LOG_INFO("------------------------------------------------------------------------");
      MACH3LOG_INFO("Setting {} = {}", ParamName, val);
      xsec->SetParProp(ParIndex, val); // <-- confirm this matches the real API

      for (auto handler : DUNEPdfs) {
        handler->Reweight();
        for (int iSample = 0; iSample < handler->GetNsamples(); iSample++) {
          std::string title = handler->GetSampleTitle(iSample);
          TH1* Hist = (TH1*)handler->GetMCHist(iSample)->Clone(
              Form("%s_%s_%.2f", title.c_str(), ParamName.c_str(), val));
          Hist->SetDirectory(nullptr); // detach from any open TFile so it survives
          SysHists[title][val] = Hist;

          MACH3LOG_INFO("[{} = {:.2f}] Event rate for {:<20} : {:.2f}",
                        ParamName, val, title, Hist->Integral());
        }
      }
    }

    // Print per-bin differences so you can see the shape effect, not just the integral
    MACH3LOG_INFO("------------------------------------------------------------------------");
    MACH3LOG_INFO("Bin-by-bin comparison ({} = 0.0 vs {} = 0.2):", ParamName, ParamName);
    for (auto& [title, valMap] : SysHists) {
      TH1* h0 = valMap.at(0.0);
      TH1* h2 = valMap.at(0.2);
      MACH3LOG_INFO("Sample: {}", title);
      for (int b = 1; b <= h0->GetNbinsX(); ++b) {
        double c0 = h0->GetBinContent(b);
        double c2 = h2->GetBinContent(b);
        double diff = c2 - c0;
        MACH3LOG_INFO("  bin {:<3} : nominal(0.0)={:.3f}  shifted(0.2)={:.3f}  diff={:.3f}",
                      b, c0, c2, diff);
      }
      if (std::abs(h0->Integral() - h2->Integral()) < 1e-9 &&
          [&]() { for (int b = 1; b <= h0->GetNbinsX(); ++b) if (std::abs(h0->GetBinContent(b) - h2->GetBinContent(b)) > 1e-9) return false; return true; }()) {
        MACH3LOG_WARN("Histograms for {} are IDENTICAL between par=0.0 and par=0.2 — systematic is having NO effect!", title);
      }
    }

    // Dump both sets of histograms to file/pdf for visual inspection
    std::vector<TH1*> FlatHists;
    for (auto& [title, valMap] : SysHists)
      for (auto& [val, h] : valMap)
        FlatHists.push_back(h);

    Write1DHistogramsToFile("MissingProtonFD_test.root", FlatHists);
    Write1DHistogramsToPdf("MissingProtonFD_test.root", FlatHists);
    WriteOverlayPdf("MissingProtonFD_test.root", SysHists, ParamName);

    // Reset parameter back to nominal before continuing with the rest of the app
    xsec->SetParProp(ParIndex, 0.0);
    for (auto handler : DUNEPdfs) handler->Reweight();
  }

  //###############################################################################################################################
  //Make oscillation channel breakdown

  MACH3LOG_INFO("========================================================================");
  MACH3LOG_INFO("========================================================================");
  MACH3LOG_INFO("Oscillation Mode Breakdown:");
  for (auto handler : DUNEPdfs) {
    for (int iSample = 0; iSample < handler->GetNsamples(); iSample++) {
      MACH3LOG_INFO("======================");
      int nOscChannels = handler->GetNOscChannels(iSample);
      for (int iOscChan = 0; iOscChan < nOscChannels; iOscChan++) {
        std::vector<KinematicCut> SelectionVec;

        KinematicCut SelecChannel;
        SelecChannel.ParamToCutOnIt = handler->ReturnKinematicParameterFromString("OscillationChannel");
        SelecChannel.LowerBound = iOscChan;
        SelecChannel.UpperBound = iOscChan + 1;
        SelectionVec.push_back(SelecChannel);
        TH1* Hist = handler->Get1DVarHist(iSample, handler->GetXBinVarName(iSample), SelectionVec);
        MACH3LOG_INFO("{:<20} : {:<20} : {:<20.2f}", handler->GetSampleTitle(iSample), handler->GetFlavourName(iSample, iOscChan), Hist->Integral());
      }

      TH1* Hist = handler->Get1DVarHist(iSample, handler->GetXBinVarName(iSample));
      MACH3LOG_INFO("{:<20} : {:<20.2f}", handler->GetSampleTitle(iSample), Hist->Integral());
    }
  }

  //###############################################################################################################################
  //Make interaction channel breakdown

  MACH3LOG_INFO("========================================================================");
  MACH3LOG_INFO("========================================================================");
  MACH3LOG_INFO("Interaction Mode Breakdown:");

  for (auto handler : DUNEPdfs) {
    for (int iSample = 0; iSample < handler->GetNsamples(); iSample++) {
      MACH3LOG_INFO("======================");

      MaCh3Modes* Modes = handler->GetMaCh3Modes();
      int nModeChannels = Modes->GetNModes();
      for (int iModeChan = 0; iModeChan < nModeChannels; iModeChan++) {
        std::vector<KinematicCut> SelectionVec;

        KinematicCut SelecChannel;
        SelecChannel.ParamToCutOnIt = handler->ReturnKinematicParameterFromString("Mode");
        SelecChannel.LowerBound = iModeChan;
        SelecChannel.UpperBound = iModeChan + 1;
        SelectionVec.push_back(SelecChannel);

        TH1* Hist = handler->Get1DVarHist(iSample, handler->GetXBinVarName(iSample), SelectionVec);
        MACH3LOG_INFO("{:<20} : {:<20} : {:<20.2f}", handler->GetSampleTitle(iSample), Modes->GetMaCh3ModeName(iModeChan), Hist->Integral());
      }

      TH1* Hist = handler->Get1DVarHist(iSample, handler->GetXBinVarName(iSample));
      MACH3LOG_INFO("{:<20} : {:<20.2f}", handler->GetSampleTitle(iSample), Hist->Integral());
    }
  }

  //###############################################################################################################################
}