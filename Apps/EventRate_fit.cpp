// #include <iostream>
// #include <chrono>
// #include <iomanip>
// #include <vector>

// #include <TH1D.h>
// #include <THStack.h>
// #include <TStyle.h>
// #include <TCanvas.h>
// #include <TRint.h>
// #include <TLegend.h>
// #include <TColor.h>
// #include <TMath.h>
// #include <map>
// #include <algorithm>
// #include <limits>
// #include <TFile.h>
// #include <TTree.h>
// #include <TBranch.h>
// #include <TGraph.h>
// #include <TMultiGraph.h>

// #include "Samples/MaCh3DUNEFactory.h"
// #include "Samples/StructsDUNE.h"
// #include "Fitters/MaCh3Factory.h"

// //###############################################################################################################################
// // Generic helpers
// //###############################################################################################################################

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
// // reweighted samples. Results are written out as a multi-page PDF, plus an
// // overlay-per-sample PDF and a one-histogram-per-page PDF of the spectrum
// // at each step.
// //###############################################################################################################################

// struct ChainStepResult {
//   Long64_t Step = 0;
//   double LogL_Chain = std::numeric_limits<double>::quiet_NaN();
//   double LogL_Recomputed = std::numeric_limits<double>::quiet_NaN();
//   std::map<std::string, double> SampleRates; // sample title -> predicted event rate
//   std::map<std::string, TH1*> SampleHists;   // sample title -> cloned spectrum at this step
//   double TotalRate = 0.0;
// };

// // Try to find, for every parameter known to xsec, which branch in the
// // posterior chain corresponds to it.
// //
// // MaCh3 has used a couple of different naming conventions for these branches
// // across versions ("xsec_<i>" in older releases, "param_<i>" in newer ones),
// // and some setups additionally store the human-readable parameter name in
// // the branch *title* rather than the branch name. This function tries the
// // obvious candidates and falls back to matching on branch title, so in most
// // cases you shouldn't need to edit it — but if none of your parameters are
// // found, open the chain in root -l and check Chain->Print() to see what
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

//   int NPars = xsec->GetNumParams(); // <-- confirm this matches 
//   for (int i = 0; i < NPars; ++i) {
//     std::string ParName = xsec->GetParName(i); 

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

// // One page per (sample, step), unoverlaid — grouped sample-by-sample so all
// // steps for a given sample are together in the PDF.
// void WriteSpectraPerStepPdf(const std::string& OutFileName,
//                              std::vector<ChainStepResult>& Results) {
//   std::string PdfName = OutFileName;
//   PdfName.erase(PdfName.find('.'));
//   PdfName += "_spectra_per_step.pdf";

//   gStyle->SetOptStat(0);

//   auto c1 = std::unique_ptr<TCanvas>(new TCanvas("c_spectra", "c_spectra", 900, 700));
//   c1->cd();
//   c1->Print(std::string(PdfName + "[").c_str());

//   // Collect sample titles from the first step's results, so pages are grouped by sample
//   std::vector<std::string> SampleTitles;
//   if (!Results.empty()) {
//     for (auto& [title, h] : Results.front().SampleHists) SampleTitles.push_back(title);
//   }

//   for (auto& title : SampleTitles) {
//     for (auto& Res : Results) {
//       auto it = Res.SampleHists.find(title);
//       if (it == Res.SampleHists.end()) continue;

//       TH1* h = (TH1*)it->second->Clone(Form("%s_step%lld_page", title.c_str(), (long long)Res.Step));
//       h->SetDirectory(nullptr);
//       h->SetLineColor(kAzure+2);
//       h->SetLineWidth(2);
//       h->SetFillColorAlpha(kAzure+2, 0.25);
//       h->SetTitle(Form("%s : Step %lld;%s;Events",
//                         title.c_str(), (long long)Res.Step, h->GetXaxis()->GetTitle()));

//       h->Draw("HIST");
//       c1->Print(PdfName.c_str());
//       delete h;
//     }
//   }

//   c1->Print(std::string(PdfName + "]").c_str());
//   MACH3LOG_INFO("Wrote per-step spectra PDF: {}", PdfName);
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

//   int NPars = xsec->GetNumParams();
//   std::vector<double> ParVals(NPars, 0.0);
//   for (auto& [i, name] : ParToBranch) {
//     Chain->SetBranchAddress(name.c_str(), &ParVals[i]);
//   }

//   double LogL_Chain = std::numeric_limits<double>::quiet_NaN();
//   bool HaveLogLBranch = (Chain->GetBranch(LogLBranchName.c_str()) != nullptr);
//   if (HaveLogLBranch) {
//     Chain->SetBranchAddress(LogLBranchName.c_str(), &LogL_Chain);
//   } else {
//     MACH3LOG_WARN("Could not find LogL branch '{}' in chain, stored LogL will be reported as NaN", LogLBranchName);
//   }

//   Long64_t NEntries = Chain->GetEntries();
//   std::vector<ChainStepResult> Results;

//   for (Long64_t step : StepsToProcess) {
//     if (step < 0 || step >= NEntries) {
//       MACH3LOG_WARN("Requested step {} is outside the chain (0..{}) — skipping", step, NEntries - 1);
//       continue;
//     }
//     Chain->GetEntry(step);
//     std::cout << "\n=====================================================\n";
//     std::cout << "Processing chain step " << step << std::endl;
//     std::cout << "=====================================================\n";

//     for (auto& [i, name] : ParToBranch) {
//       xsec->SetParProp(i, ParVals[i]);
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
//         TH1* h = handler->GetMCHist(iSample);
//         double Rate = h->Integral();
//         Res.SampleRates[title] += Rate;
//         TotalRate += Rate;

//         // Clone and stash the spectrum for this step so we can plot it later
//         TH1* hClone = (TH1*)h->Clone(Form("%s_step%lld", title.c_str(), (long long)step));
//         hClone->SetDirectory(nullptr);
//         Res.SampleHists[title] = hClone;
//       }
//     }
 
//     RecomputedLogL += xsec->GetLikelihood();

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
//   // Make the summary PDF (rates + likelihoods vs step)
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

//   //###################################
//   // Overlay the event-rate spectrum at each step, one page per sample
//   //###################################
//   {
//     std::map<std::string, std::map<double, TH1*>> SpectraByStep;
//     for (auto& Res : Results) {
//       for (auto& [title, h] : Res.SampleHists) {
//         SpectraByStep[title][(double)Res.Step] = h;
//       }
//     }
//     WriteOverlayPdf(OutPdfName, SpectraByStep, "Step");
//   }

//   //###################################
//   // Separate page per (sample, step) — no overlay, one histogram per page
//   //###################################
//   WriteSpectraPerStepPdf(OutPdfName, Results);

//   // Clean up the cloned histograms we stashed in Res.SampleHists
//   for (auto& Res : Results) {
//     for (auto& [title, h] : Res.SampleHists) {
//       delete h;
//     }
//   }
// }

// //###############################################################################################################################
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


//   //###############################################################################################################################
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

//       //double current = xsec->GetParProp(i);

      

//       for (auto handler : DUNEPdfs) {
//         double before = 0.0;

// for(int i=0;i<handler->GetNsamples();i++)
//     before += handler->GetMCHist(i)->Integral();

//     std::cout << "Total integral before Reweight = "
//               << before
//               << std::endl;

//     handler->Reweight();

//     double after = 0.0;

//     for(int i=0;i<handler->GetNsamples();i++)
//         after += handler->GetMCHist(i)->Integral();

//     std::cout << "Total integral after Reweight  = "
//               << after
//               << std::endl;
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
//   // Set General:PosteriorChain:FileName in your config YAML to point at the
//   // chain, e.g.:
//   //   General:
//   //     PosteriorChain:
//   //       FileName: "/path/to/mcmc_chain.root"
//   //       TreeName: "posteriors"     # optional, defaults shown
//   //       LogLBranch: "LogL"         # optional, defaults shown
//   //       OutputPdf: "PosteriorChainCheck.pdf"  # optional
//   //
//   // This produces three PDFs:
//   //   <OutputPdf>                          - rate/LogL summary vs step (4 pages)
//   //   <OutputPdf minus ext>_overlay.pdf     - one page per sample, all steps overlaid
//   //   <OutputPdf minus ext>_spectra_per_step.pdf - one page per (sample, step)

//   std::string ChainFileName = GetFromManager<std::string>(fitMan->raw()["General"]["PosteriorChain"]["FileName"], "");
//   std::string ChainTreeName = GetFromManager<std::string>(fitMan->raw()["General"]["PosteriorChain"]["TreeName"], "posteriors");
//   std::string LogLBranchName = GetFromManager<std::string>(fitMan->raw()["General"]["PosteriorChain"]["LogLBranch"], "LogL");
//   std::string ChainPdfName = GetFromManager<std::string>(fitMan->raw()["General"]["PosteriorChain"]["OutputPdf"], "PosteriorChainCheck.pdf");

//   // <-- Edit this to whichever chain steps you want to check
//   std::vector<Long64_t> StepsToProcess = {0, 5, 10, 70, 100 , 5000, 10000, 50000, 60000};

//   if (!ChainFileName.empty()) {
//     RunPosteriorChainCheck(ChainFileName, ChainTreeName, LogLBranchName, StepsToProcess, xsec, DUNEPdfs, ChainPdfName);
//     for (int i = 0; i < xsec->GetNumParams(); ++i) xsec->SetParProp(i, 0.0);
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

#include <TH1D.h>
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRint.h>
#include <TLegend.h>
#include <TColor.h>
#include <TMath.h>
#include <map>
#include <algorithm>
#include <limits>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TGraph.h>
#include <TMultiGraph.h>

#include "Samples/MaCh3DUNEFactory.h"
#include "Samples/StructsDUNE.h"
#include "Fitters/MaCh3Factory.h"

//###############################################################################################################################
// Generic helpers
//###############################################################################################################################

// Pulls the currently-attached data histogram (set via handler->AddData(...))
// for every sample across every handler, keyed by sample title, and returns
// owned clones (detached from any TFile) that the caller is responsible for
// deleting.
//
// Uses SampleHandlerBase::GetDataHist(iSample), which returns a const TH1*
// pointing at whatever was last passed to AddData() for that sample.
std::map<std::string, TH1*> CollectDataHists(std::vector<SampleHandlerFD*>& DUNEPdfs) {
  std::map<std::string, TH1*> DataHistBySample;

  for (auto handler : DUNEPdfs) {
    for (int iSample = 0; iSample < handler->GetNsamples(); iSample++) {
      std::string title = handler->GetSampleTitle(iSample);
      const TH1* hData = handler->GetDataHist(iSample);
      if (!hData) {
        MACH3LOG_WARN("No data histogram attached for sample {}, it will be skipped in data comparisons", title);
        continue;
      }
      TH1* hClone = (TH1*)hData->Clone(Form("%s_data_clone", title.c_str()));
      hClone->SetDirectory(nullptr);
      DataHistBySample[title] = hClone;
    }
  }

  return DataHistBySample;
}

void WriteOverlayPdf(const std::string& OutFileName,
                      std::map<std::string, std::map<double, TH1*>>& SysHists,
                      const std::string& ParamName,
                      std::map<std::string, TH1*>* DataHists = nullptr) {

  std::string PdfName = OutFileName;
  PdfName.erase(PdfName.find('.'));
  PdfName += "_overlay.pdf";

  gStyle->SetOptStat(0);

  auto c1 = std::unique_ptr<TCanvas>(new TCanvas("c_overlay", "c_overlay", 900, 700));
  c1->cd();
  c1->Print(std::string(PdfName + "[").c_str());

  for (auto& [title, valMap] : SysHists) {
    // Clone so we don't mutate the histograms other code might still use
    std::vector<std::pair<double, TH1*>> sortedVals(valMap.begin(), valMap.end());
    std::sort(sortedVals.begin(), sortedVals.end(),
              [](auto& a, auto& b){ return a.first < b.first; });

    std::vector<TH1*> clones;
    double maxY = 0;

    int colours[] = {kBlack, kRed+1, kAzure+2, kGreen+2, kMagenta+1};
    int ic = 0;

    auto leg = std::make_unique<TLegend>(0.65, 0.65, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    for (auto& [val, hOrig] : sortedVals) {
      TH1* h = (TH1*)hOrig->Clone(Form("%s_clone_%.3f", hOrig->GetName(), val));
      h->SetDirectory(nullptr);

      // Divide out bin widths: content -> content / bin_width
      h->Scale(1.0, "width");

      h->SetLineColor(colours[ic % 5]);
      h->SetLineWidth(2);
      h->SetMarkerColor(colours[ic % 5]);
      h->SetMarkerStyle(20);
      h->SetTitle(Form("%s;%s;Events / bin width", title.c_str(), h->GetXaxis()->GetTitle()));

      maxY = std::max(maxY, h->GetMaximum());
      leg->AddEntry(h, Form("%s = %.2f", ParamName.c_str(), val), "l");
      clones.push_back(h);
      ic++;
    }

    // Data overlay, if provided and present for this sample
    TH1* hData = nullptr;
    if (DataHists && DataHists->count(title)) {
      hData = (TH1*)DataHists->at(title)->Clone(Form("%s_data_overlay_clone", title.c_str()));
      hData->SetDirectory(nullptr);
      hData->Scale(1.0, "width");
      hData->SetMarkerStyle(20);
      hData->SetMarkerColor(kBlack);
      hData->SetLineColor(kBlack);
      maxY = std::max(maxY, hData->GetMaximum());
      leg->AddEntry(hData, "Data", "p");
    }

    for (size_t i = 0; i < clones.size(); ++i) {
      clones[i]->SetMaximum(maxY * 1.3);
      clones[i]->Draw(i == 0 ? "HIST" : "HIST SAME");
    }
    if (hData) hData->Draw("P SAME");
    leg->Draw();
    c1->Print(PdfName.c_str());

    for (auto h : clones) delete h;
    if (hData) delete hData;
  }

  c1->Print(std::string(PdfName + "]").c_str());
  MACH3LOG_INFO("Wrote overlay PDF: {}", PdfName);
}

void Write1DHistogramsToFile(std::string OutFileName, std::vector<TH1*> Histograms) {
  auto OutputFile = std::unique_ptr<TFile>(TFile::Open(OutFileName.c_str(), "RECREATE"));
  OutputFile->cd();
  for(auto Hist : Histograms){
    Hist->Write();
  }
  OutputFile->Close();
}

void Write1DHistogramsToPdf(std::string OutFileName, std::vector<TH1*> Histograms) {
  //Remove root from end of file
  OutFileName.erase(OutFileName.find('.'));
  OutFileName+=".pdf";

  auto c1 = std::unique_ptr<TCanvas>(new TCanvas("c1", "c1", 800, 600));
  c1->cd();
  c1->Print(std::string(OutFileName+"[").c_str());
  for(auto Hist : Histograms){
    Hist->Draw("HIST COLZ");
    c1->Print(OutFileName.c_str());
  }
  c1->Print(std::string(OutFileName+"]").c_str());
}

//###############################################################################################################################
// Posterior-chain event-rate / likelihood check
//
// Reads a standard MaCh3 MCMC output tree (one branch per parameter, plus a
// LogL branch), and for a user-specified list of steps: pushes those
// parameter values into the systematic handler, reweights every sample, and
// records the predicted event rate per sample together with both the LogL
// stored in the chain at that step and a freshly recomputed LogL from the
// reweighted samples. The (fixed) data rate attached via AddData() in main()
// is also recorded for comparison. Results are written out as a multi-page
// PDF, plus an overlay-per-sample PDF (with data overlaid) and a
// one-histogram-per-page PDF of the spectrum at each step.
//###############################################################################################################################

struct ChainStepResult {
  Long64_t Step = 0;
  double LogL_Chain = std::numeric_limits<double>::quiet_NaN();
  double LogL_Recomputed = std::numeric_limits<double>::quiet_NaN();
  std::map<std::string, double> SampleRates;     // sample title -> predicted event rate
  std::map<std::string, double> SampleDataRates; // sample title -> data event rate (fixed across steps)
  std::map<std::string, TH1*> SampleHists;       // sample title -> cloned spectrum at this step
  double TotalRate = 0.0;
  double TotalDataRate = 0.0;
};

// Try to find, for every parameter known to xsec, which branch in the
// posterior chain corresponds to it.
//
// MaCh3 has used a couple of different naming conventions for these branches
// across versions ("xsec_<i>" in older releases, "param_<i>" in newer ones),
// and some setups additionally store the human-readable parameter name in
// the branch *title* rather than the branch name. This function tries the
// obvious candidates and falls back to matching on branch title, so in most
// cases you shouldn't need to edit it — but if none of your parameters are
// found, open the chain in root -l and check Chain->Print() to see what
// your branches are actually called, then add the pattern below.
std::map<int, std::string> MapParamsToChainBranches(TTree* Chain, ParameterHandlerGeneric* xsec) {
  std::map<int, std::string> ParToBranch;

  std::vector<std::string> BranchNames;
  std::vector<std::string> BranchTitles;
  for (auto* obj : *Chain->GetListOfBranches()) {
    auto* br = (TBranch*)obj;
    BranchNames.push_back(br->GetName());
    BranchTitles.push_back(br->GetTitle());
  }

  int NPars = xsec->GetNumParams(); // <-- confirm this matches
  for (int i = 0; i < NPars; ++i) {
    std::string ParName = xsec->GetParName(i);

    std::vector<std::string> Candidates = {
      Form("param_%d", i),
      Form("xsec_%d", i),
      ParName,
      Form("param_%s", ParName.c_str()),
      Form("xsec_%s", ParName.c_str())
    };

    bool Found = false;
    for (auto& cand : Candidates) {
      if (std::find(BranchNames.begin(), BranchNames.end(), cand) != BranchNames.end()) {
        ParToBranch[i] = cand;
        Found = true;
        break;
      }
    }

    // Fall back to matching on branch title, which MaCh3 often sets to the
    // human-readable parameter name even when the branch itself is indexed
    if (!Found) {
      for (size_t b = 0; b < BranchTitles.size(); ++b) {
        if (BranchTitles[b] == ParName) {
          ParToBranch[i] = BranchNames[b];
          Found = true;
          break;
        }
      }
    }

    if (!Found) {
      MACH3LOG_WARN("Could not find a chain branch for parameter {} ({}) "
                     "— it will be left at its current value for every step",
                     i, ParName);
    }
  }

  return ParToBranch;
}

// One page per (sample, step), unoverlaid — grouped sample-by-sample so all
// steps for a given sample are together in the PDF.
void WriteSpectraPerStepPdf(const std::string& OutFileName,
                             std::vector<ChainStepResult>& Results) {
  std::string PdfName = OutFileName;
  PdfName.erase(PdfName.find('.'));
  PdfName += "_spectra_per_step.pdf";

  gStyle->SetOptStat(0);

  auto c1 = std::unique_ptr<TCanvas>(new TCanvas("c_spectra", "c_spectra", 900, 700));
  c1->cd();
  c1->Print(std::string(PdfName + "[").c_str());

  // Collect sample titles from the first step's results, so pages are grouped by sample
  std::vector<std::string> SampleTitles;
  if (!Results.empty()) {
    for (auto& [title, h] : Results.front().SampleHists) SampleTitles.push_back(title);
  }

  for (auto& title : SampleTitles) {
    for (auto& Res : Results) {
      auto it = Res.SampleHists.find(title);
      if (it == Res.SampleHists.end()) continue;

      TH1* h = (TH1*)it->second->Clone(Form("%s_step%lld_page", title.c_str(), (long long)Res.Step));
      h->SetDirectory(nullptr);
      h->SetLineColor(kAzure+2);
      h->SetLineWidth(2);
      h->SetFillColorAlpha(kAzure+2, 0.25);
      h->SetTitle(Form("%s : Step %lld;%s;Events",
                        title.c_str(), (long long)Res.Step, h->GetXaxis()->GetTitle()));

      h->Draw("HIST");
      c1->Print(PdfName.c_str());
      delete h;
    }
  }

  c1->Print(std::string(PdfName + "]").c_str());
  MACH3LOG_INFO("Wrote per-step spectra PDF: {}", PdfName);
}

void RunPosteriorChainCheck(const std::string& ChainFileName,
                             const std::string& ChainTreeName,
                             const std::string& LogLBranchName,
                             const std::vector<Long64_t>& StepsToProcess,
                             ParameterHandlerGeneric* xsec,
                             std::vector<SampleHandlerFD*>& DUNEPdfs,
                             const std::string& OutPdfName) {

  MACH3LOG_INFO("========================================================================");
  MACH3LOG_INFO("Posterior-chain event-rate / likelihood check");
  MACH3LOG_INFO("Chain file: {}, tree: {}", ChainFileName, ChainTreeName);

  auto ChainFile = std::unique_ptr<TFile>(TFile::Open(ChainFileName.c_str(), "READ"));
  if (!ChainFile || ChainFile->IsZombie()) {
    MACH3LOG_ERROR("Could not open chain file {}", ChainFileName);
    return;
  }
  auto* Chain = (TTree*)ChainFile->Get(ChainTreeName.c_str());
  if (!Chain) {
    MACH3LOG_ERROR("Could not find tree {} in {}", ChainTreeName, ChainFileName);
    return;
  }

  auto ParToBranch = MapParamsToChainBranches(Chain, xsec);
  if (ParToBranch.empty()) {
    MACH3LOG_ERROR("No parameter branches were matched in the chain, aborting check");
    return;
  }

  int NPars = xsec->GetNumParams();
  std::vector<double> ParVals(NPars, 0.0);
  for (auto& [i, name] : ParToBranch) {
    Chain->SetBranchAddress(name.c_str(), &ParVals[i]);
  }

  double LogL_Chain = std::numeric_limits<double>::quiet_NaN();
  bool HaveLogLBranch = (Chain->GetBranch(LogLBranchName.c_str()) != nullptr);
  if (HaveLogLBranch) {
    Chain->SetBranchAddress(LogLBranchName.c_str(), &LogL_Chain);
  } else {
    MACH3LOG_WARN("Could not find LogL branch '{}' in chain, stored LogL will be reported as NaN", LogLBranchName);
  }

  // Data was attached via AddData() in main() before this function was
  // called, so pull the (fixed) data histogram/rate per sample once here —
  // it doesn't change step to step.
  std::map<std::string, TH1*> DataHistBySample = CollectDataHists(DUNEPdfs);
  std::map<std::string, double> DataRateBySample;
  double TotalDataRate = 0.0;
  for (auto& [title, h] : DataHistBySample) {
    double rate = h->Integral();
    DataRateBySample[title] = rate;
    TotalDataRate += rate;
  }

  Long64_t NEntries = Chain->GetEntries();
  std::vector<ChainStepResult> Results;

  for (Long64_t step : StepsToProcess) {
    if (step < 0 || step >= NEntries) {
      MACH3LOG_WARN("Requested step {} is outside the chain (0..{}) — skipping", step, NEntries - 1);
      continue;
    }
    Chain->GetEntry(step);
    std::cout << "\n=====================================================\n";
    std::cout << "Processing chain step " << step << std::endl;
    std::cout << "=====================================================\n";

    for (auto& [i, name] : ParToBranch) {
      xsec->SetParProp(i, ParVals[i]);
    }

    ChainStepResult Res;
    Res.Step = step;
    Res.LogL_Chain = LogL_Chain;
    Res.SampleDataRates = DataRateBySample;
    Res.TotalDataRate = TotalDataRate;

    double TotalRate = 0.0;
    double RecomputedLogL = 0.0;

    for (auto handler : DUNEPdfs) {
      handler->Reweight();
      RecomputedLogL += handler->GetLikelihood(); // <-- confirm this matches the real API
      for (int iSample = 0; iSample < handler->GetNsamples(); iSample++) {
        std::string title = handler->GetSampleTitle(iSample);
        TH1* h = handler->GetMCHist(iSample);
        double Rate = h->Integral();
        Res.SampleRates[title] += Rate;
        TotalRate += Rate;

        // Clone and stash the spectrum for this step so we can plot it later
        TH1* hClone = (TH1*)h->Clone(Form("%s_step%lld", title.c_str(), (long long)step));
        hClone->SetDirectory(nullptr);
        Res.SampleHists[title] = hClone;
      }
    }

    RecomputedLogL += xsec->GetLikelihood();

    Res.TotalRate = TotalRate;
    Res.LogL_Recomputed = RecomputedLogL;

    MACH3LOG_INFO("------------------------------------------------------------------------");
    MACH3LOG_INFO("Step {:<8} : TotalRate={:.3f}  TotalData={:.3f}  diff={:.3f}  LogL(chain)={:.3f}  LogL(recomputed)={:.3f}  diff={:.3f}",
                  Res.Step, Res.TotalRate, Res.TotalDataRate, Res.TotalRate - Res.TotalDataRate,
                  Res.LogL_Chain, Res.LogL_Recomputed, Res.LogL_Chain - Res.LogL_Recomputed);
    for (auto& [title, rate] : Res.SampleRates) {
      double dataRate = Res.SampleDataRates.count(title) ? Res.SampleDataRates.at(title) : 0.0;
      MACH3LOG_INFO("    {:<25} : MC={:.3f}  Data={:.3f}  diff={:.3f}", title, rate, dataRate, rate - dataRate);
    }

    Results.push_back(Res);
  }

  if (Results.empty()) {
    MACH3LOG_ERROR("No valid steps were processed, not writing PDF");
    for (auto& [title, h] : DataHistBySample) delete h;
    return;
  }

  //###################################
  // Make the summary PDF (rates + likelihoods vs step)
  //###################################
  gStyle->SetOptStat(0);
  auto c1 = std::unique_ptr<TCanvas>(new TCanvas("c_chain", "c_chain", 900, 700));
  c1->cd();
  c1->Print(std::string(OutPdfName + "[").c_str());

  std::vector<std::string> SampleTitles;
  for (auto& [title, rate] : Results.front().SampleRates) SampleTitles.push_back(title);

  int NSteps = (int)Results.size();

  // --- Page 1: total predicted event rate vs step (+ flat data reference line) ---
  {
    TGraph g(NSteps);
    for (int i = 0; i < NSteps; ++i) g.SetPoint(i, (double)Results[i].Step, Results[i].TotalRate);
    g.SetTitle("Total predicted event rate;Chain step;Events");
    g.SetLineColor(kBlack);
    g.SetLineWidth(2);
    g.SetMarkerStyle(20);
    g.SetMarkerColor(kBlack);

    TGraph gData(2);
    gData.SetPoint(0, (double)Results.front().Step, Results.front().TotalDataRate);
    gData.SetPoint(1, (double)Results.back().Step, Results.back().TotalDataRate);
    gData.SetLineColor(kRed+1);
    gData.SetLineStyle(2);
    gData.SetLineWidth(2);

    TMultiGraph mg;
    mg.Add(&g, "LP");
    mg.Add(&gData, "L");
    mg.SetTitle("Total predicted event rate (dashed red = data);Chain step;Events");
    mg.Draw("A");

    auto leg = std::make_unique<TLegend>(0.6, 0.75, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(&g, "MC prediction", "l");
    leg->AddEntry(&gData, "Data", "l");
    leg->Draw();
    c1->Print(OutPdfName.c_str());
  }

  // --- Page 2: per-sample predicted event rate vs step (+ flat data reference line each) ---
  {
    TMultiGraph mg;
    auto leg = std::make_unique<TLegend>(0.6, 0.55, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    int colours[] = {kBlack, kRed+1, kAzure+2, kGreen+2, kMagenta+1, kOrange+1};
    int ic = 0;
    std::vector<std::unique_ptr<TGraph>> Graphs;
    for (auto& title : SampleTitles) {
      auto g = std::make_unique<TGraph>(NSteps);
      for (int i = 0; i < NSteps; ++i) g->SetPoint(i, (double)Results[i].Step, Results[i].SampleRates.at(title));
      g->SetLineColor(colours[ic % 6]);
      g->SetMarkerColor(colours[ic % 6]);
      g->SetMarkerStyle(20);
      g->SetLineWidth(2);
      leg->AddEntry(g.get(), title.c_str(), "l");
      mg.Add(g.get(), "LP");
      Graphs.push_back(std::move(g));

      // Flat dashed line at the data rate for this sample
      double dataRate = Results.front().SampleDataRates.count(title) ? Results.front().SampleDataRates.at(title) : 0.0;
      auto gData = std::make_unique<TGraph>(2);
      gData->SetPoint(0, (double)Results.front().Step, dataRate);
      gData->SetPoint(1, (double)Results.back().Step, dataRate);
      gData->SetLineColor(colours[ic % 6]);
      gData->SetLineStyle(2);
      gData->SetLineWidth(1);
      mg.Add(gData.get(), "L");
      Graphs.push_back(std::move(gData));

      ic++;
    }
    mg.SetTitle("Predicted event rate per sample (dashed = data);Chain step;Events");
    mg.Draw("A");
    leg->Draw();
    c1->Print(OutPdfName.c_str());
  }

  // --- Page 3: LogL stored in chain vs recomputed LogL ---
  {
    TGraph gChain(NSteps), gRecomp(NSteps);
    for (int i = 0; i < NSteps; ++i) {
      gChain.SetPoint(i, (double)Results[i].Step, Results[i].LogL_Chain);
      gRecomp.SetPoint(i, (double)Results[i].Step, Results[i].LogL_Recomputed);
    }
    gChain.SetLineColor(kBlack);   gChain.SetMarkerColor(kBlack);   gChain.SetMarkerStyle(20); gChain.SetLineWidth(2);
    gRecomp.SetLineColor(kRed+1);  gRecomp.SetMarkerColor(kRed+1);  gRecomp.SetMarkerStyle(21); gRecomp.SetLineWidth(2);

    TMultiGraph mg;
    mg.Add(&gChain, "LP");
    mg.Add(&gRecomp, "LP");
    mg.SetTitle("-log L at each step;Chain step;-log L");
    mg.Draw("A");

    auto leg = std::make_unique<TLegend>(0.6, 0.75, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(&gChain, "Stored in chain", "l");
    leg->AddEntry(&gRecomp, "Recomputed after reweight", "l");
    leg->Draw();
    c1->Print(OutPdfName.c_str());
  }

  // --- Page 4: difference between stored and recomputed LogL ---
  {
    TGraph g(NSteps);
    for (int i = 0; i < NSteps; ++i)
      g.SetPoint(i, (double)Results[i].Step, Results[i].LogL_Chain - Results[i].LogL_Recomputed);
    g.SetTitle("LogL(chain) - LogL(recomputed);Chain step;#Delta(-log L)");
    g.SetLineColor(kAzure+2);
    g.SetLineWidth(2);
    g.SetMarkerColor(kAzure+2);
    g.SetMarkerStyle(20);
    g.Draw("ALP");
    c1->Print(OutPdfName.c_str());
  }

  c1->Print(std::string(OutPdfName + "]").c_str());
  MACH3LOG_INFO("Wrote posterior-chain check PDF: {}", OutPdfName);

  //###################################
  // Overlay the event-rate spectrum at each step, one page per sample,
  // with the fixed data histogram overlaid on top
  //###################################
  {
    std::map<std::string, std::map<double, TH1*>> SpectraByStep;
    for (auto& Res : Results) {
      for (auto& [title, h] : Res.SampleHists) {
        SpectraByStep[title][(double)Res.Step] = h;
      }
    }
    WriteOverlayPdf(OutPdfName, SpectraByStep, "Step", &DataHistBySample);
  }

  //###################################
  // Separate page per (sample, step) — no overlay, one histogram per page
  //###################################
  WriteSpectraPerStepPdf(OutPdfName, Results);

  // Clean up the cloned histograms we stashed in Res.SampleHists
  for (auto& Res : Results) {
    for (auto& [title, h] : Res.SampleHists) {
      delete h;
    }
  }
  for (auto& [title, h] : DataHistBySample) {
    delete h;
  }
}

//###############################################################################################################################
int main(int argc, char * argv[]) {
  MaCh3Utils::MaCh3Usage(argc, argv);
  auto fitMan = MaCh3ManagerFactory(argc, argv);

  //###############################################################################################################################
  //Create SampleHandlerFD objects

  ParameterHandlerGeneric* xsec = nullptr;

  std::vector<SampleHandlerFD*> DUNEPdfs;
  MakeMaCh3DuneInstance(fitMan, DUNEPdfs, xsec);

  //###############################################################################################################################
  //Perform reweight and print total integral

  std::vector<TH1*> DUNEHists;
  for(auto handler : DUNEPdfs){
    handler->Reweight();
    for (int iSample=0; iSample<handler->GetNsamples(); iSample++) {
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
  // Build the fake ("Asimov") data exactly the way the production fit app
  // does, so that every LogL evaluation and data-comparison plot from this
  // point on (the MissingProtonFD systematic check and the posterior-chain
  // check) is measured against the SAME dataset the chain was actually fit
  // to - not a generic nominal MC prediction.
  //
  // The fit app (EventRate_fit's companion) generates its fake data with
  // MissingProtonFD = 0.2, then fixes MissingProtonFD at 0.0 for the actual
  // MCMC - i.e. the chain can never directly correct for this offset via
  // that parameter and has to compensate through other, correlated
  // parameters. Reproducing a plain nominal (all-params-at-0) Asimov set
  // here would silently compare the chain against the wrong target,
  // especially for steps that have moved away from the starting position to
  // compensate for the baked-in offset - exactly the "further from step 0
  // fits worse" pattern being investigated.
  MACH3LOG_INFO("========================================================================");
  MACH3LOG_INFO("Setting fake ('Asimov') data with MissingProtonFD = 0.2, matching the fit app");

  int MissingProtonFDIdx = xsec->GetParIndex("MissingProtonFD");
  if (MissingProtonFDIdx < 0) {
    MACH3LOG_ERROR("Could not find parameter MissingProtonFD in ParameterHandler - "
                    "cannot reproduce the fit's fake-data recipe, data comparisons below will be wrong");
  } else {
    xsec->SetSingleParameter(MissingProtonFDIdx, 0.2);
    for (auto handler : DUNEPdfs) handler->Reweight();
  }

  // Keep these alive for the lifetime of the run: AddData is assumed to
  // store the pointer rather than copy it (mirroring the pattern used in
  // the accompanying fit app), so the histograms must not be deleted or go
  // out of scope while DUNEPdfs are in use below.
  std::vector<TH1*> DataHistograms;
  for (auto handler : DUNEPdfs) {
    for (unsigned iSample = 0; iSample < handler->GetNsamples(); ++iSample) {
      std::string name = handler->GetSampleTitle(iSample);
      TString NameTString = TString(name.c_str());

      TH1* DataHist = static_cast<TH1*>(handler->GetMCHist(iSample)->Clone(NameTString + "_DataHist"));
      DataHist->SetDirectory(nullptr);
      DataHistograms.push_back(DataHist);

      if (handler->GetNDim(iSample) == 1) {
        handler->AddData(iSample, static_cast<TH1D*>(DataHist));
      } else if (handler->GetNDim(iSample) == 2) {
        handler->AddData(iSample, static_cast<TH2D*>(DataHist));
      } else {
        MACH3LOG_ERROR("Unsupported number of dimensions > 2 - Quitting");
        throw MaCh3Exception(__FILE__, __LINE__);
      }

      MACH3LOG_INFO("Fake data integral for {} : {:.2f}", name, DataHist->Integral());
    }
  }

  // Reset MissingProtonFD back to 0.0, matching the value it's fixed at for
  // the actual MCMC fit, so every subsequent reweight in this script
  // (the systematic check and the posterior-chain check) predicts under the
  // same fixed condition the chain was run under.
  if (MissingProtonFDIdx >= 0) {
    xsec->SetSingleParameter(MissingProtonFDIdx, 0.0);
    for (auto handler : DUNEPdfs) handler->Reweight();
  }

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

    // Data is fixed regardless of the systematic value being tested, so
    // grab it once, up front.
    std::map<std::string, TH1*> DataHistBySample = CollectDataHists(DUNEPdfs);

    for (double val : TestValues) {
      MACH3LOG_INFO("------------------------------------------------------------------------");
      MACH3LOG_INFO("Setting {} = {}", ParamName, val);
      xsec->SetParProp(ParIndex, val); // <-- confirm this matches the real API

      for (auto handler : DUNEPdfs) {
        double before = 0.0;
        for (int i = 0; i < handler->GetNsamples(); i++)
          before += handler->GetMCHist(i)->Integral();

        std::cout << "Total integral before Reweight = " << before << std::endl;

        handler->Reweight();

        double after = 0.0;
        for (int i = 0; i < handler->GetNsamples(); i++)
          after += handler->GetMCHist(i)->Integral();

        std::cout << "Total integral after Reweight  = " << after << std::endl;

        for (int iSample=0; iSample<handler->GetNsamples(); iSample++) {
          std::string title = handler->GetSampleTitle(iSample);
          TH1* Hist = (TH1*)handler->GetMCHist(iSample)->Clone(
              Form("%s_%s_%.2f", title.c_str(), ParamName.c_str(), val));
          Hist->SetDirectory(nullptr); // detach from any open TFile so it survives
          SysHists[title][val] = Hist;

          double dataRate = DataHistBySample.count(title) ? DataHistBySample.at(title)->Integral() : 0.0;
          MACH3LOG_INFO("[{} = {:.2f}] Event rate for {:<20} : MC={:.2f}  Data={:.2f}  diff={:.2f}",
                        ParamName, val, title, Hist->Integral(), dataRate, Hist->Integral() - dataRate);
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
          [&](){ for(int b=1;b<=h0->GetNbinsX();++b) if(std::abs(h0->GetBinContent(b)-h2->GetBinContent(b))>1e-9) return false; return true; }()) {
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
    WriteOverlayPdf("MissingProtonFD_test.root", SysHists, ParamName, &DataHistBySample);

    for (auto& [title, h] : DataHistBySample) delete h;

    // Reset parameter back to nominal before continuing with the rest of the app
    xsec->SetParProp(ParIndex, 0.0);
    for (auto handler : DUNEPdfs) handler->Reweight();
  }

  //###############################################################################################################################
  // Posterior-chain event-rate / likelihood check
  //
  // Set General:PosteriorChain:FileName in your config YAML to point at the
  // chain, e.g.:
  //   General:
  //     PosteriorChain:
  //       FileName: "/path/to/mcmc_chain.root"
  //       TreeName: "posteriors"     # optional, defaults shown
  //       LogLBranch: "LogL"         # optional, defaults shown
  //       OutputPdf: "PosteriorChainCheck.pdf"  # optional
  //
  // This produces three PDFs:
  //   <OutputPdf>                          - rate/data/LogL summary vs step (4 pages)
  //   <OutputPdf minus ext>_overlay.pdf     - one page per sample, all steps overlaid + data
  //   <OutputPdf minus ext>_spectra_per_step.pdf - one page per (sample, step)

  std::string ChainFileName = GetFromManager<std::string>(fitMan->raw()["General"]["PosteriorChain"]["FileName"], "");
  std::string ChainTreeName = GetFromManager<std::string>(fitMan->raw()["General"]["PosteriorChain"]["TreeName"], "posteriors");
  std::string LogLBranchName = GetFromManager<std::string>(fitMan->raw()["General"]["PosteriorChain"]["LogLBranch"], "LogL");
  std::string ChainPdfName = GetFromManager<std::string>(fitMan->raw()["General"]["PosteriorChain"]["OutputPdf"], "PosteriorChainCheck.pdf");

  // <-- Edit this to whichever chain steps you want to check
  std::vector<Long64_t> StepsToProcess = {0, 155229};

  if (!ChainFileName.empty()) {
    RunPosteriorChainCheck(ChainFileName, ChainTreeName, LogLBranchName, StepsToProcess, xsec, DUNEPdfs, ChainPdfName);
    for (int i = 0; i < xsec->GetNumParams(); ++i) xsec->SetParProp(i, 0.0);
    for (auto handler : DUNEPdfs) handler->Reweight();
  } else {
    MACH3LOG_INFO("No General:PosteriorChain:FileName set in config, skipping posterior chain check");
  }

  //###############################################################################################################################
  //Make oscillation channel breakdown

  MACH3LOG_INFO("========================================================================");
  MACH3LOG_INFO("========================================================================");
  MACH3LOG_INFO("Oscillation Mode Breakdown:");

  for(auto handler : DUNEPdfs) {
    for (int iSample = 0; iSample < handler->GetNsamples(); iSample++) {
      MACH3LOG_INFO("======================");
      int nOscChannels = handler->GetNOscChannels(iSample);
      for (int iOscChan=0;iOscChan<nOscChannels;iOscChan++) {
        std::vector< KinematicCut > SelectionVec;

        KinematicCut SelecChannel;
        SelecChannel.ParamToCutOnIt = handler->ReturnKinematicParameterFromString("OscillationChannel");
        SelecChannel.LowerBound = iOscChan;
        SelecChannel.UpperBound = iOscChan+1;
        SelectionVec.push_back(SelecChannel);

        TH1* Hist = handler->Get1DVarHist(iSample, handler->GetXBinVarName(iSample),SelectionVec);
        MACH3LOG_INFO("{:<20} : {:<20} : {:<20.2f}",handler->GetSampleTitle(iSample),handler->GetFlavourName(iSample, iOscChan),Hist->Integral());
      }

      TH1* Hist = handler->Get1DVarHist(iSample, handler->GetXBinVarName(iSample));
      MACH3LOG_INFO("{:<20} : {:<20.2f}",handler->GetSampleTitle(iSample),Hist->Integral());
    }
  }

  //###############################################################################################################################
  //Make interaction channel breakdown

  MACH3LOG_INFO("========================================================================");
  MACH3LOG_INFO("========================================================================");
  MACH3LOG_INFO("Interaction Mode Breakdown:");

  for(auto handler : DUNEPdfs) {
    for (int iSample = 0; iSample < handler->GetNsamples(); iSample++) {
      MACH3LOG_INFO("======================");

      MaCh3Modes* Modes = handler->GetMaCh3Modes();
      int nModeChannels = Modes->GetNModes();
      for (int iModeChan=0;iModeChan<nModeChannels;iModeChan++) {
        std::vector< KinematicCut > SelectionVec;

        KinematicCut SelecChannel;
        SelecChannel.ParamToCutOnIt = handler->ReturnKinematicParameterFromString("Mode");
        SelecChannel.LowerBound = iModeChan;
        SelecChannel.UpperBound = iModeChan+1;
        SelectionVec.push_back(SelecChannel);

        TH1* Hist = handler->Get1DVarHist(iSample, handler->GetXBinVarName(iSample),SelectionVec);
        MACH3LOG_INFO("{:<20} : {:<20} : {:<20.2f}",handler->GetSampleTitle(iSample),Modes->GetMaCh3ModeName(iModeChan),Hist->Integral());
      }

      TH1* Hist = handler->Get1DVarHist(iSample, handler->GetXBinVarName(iSample));
      MACH3LOG_INFO("{:<20} : {:<20.2f}",handler->GetSampleTitle(iSample),Hist->Integral());
    }
  }

  //###############################################################################################################################
}
