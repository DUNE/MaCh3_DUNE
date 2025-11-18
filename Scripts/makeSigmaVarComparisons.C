#include "TH1.h"
#include "TList.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"


#include <iostream>

TString getXsecParName(TString keyname)
{
  if (keyname.EndsWith("_par_0"))
    keyname.ReplaceAll("_par_0","_0 MaQE");
  else if (keyname.EndsWith("_par_1"))
    keyname.ReplaceAll("_par_1","_1 VecFFCCQEshape");
  else if (keyname.EndsWith("_par_2"))
    keyname.ReplaceAll("_par_2","_2 CCQEPauliSupViaKF");
  else if (keyname.EndsWith("_par_3"))
    keyname.ReplaceAll("_par_3","_3 MaCCRES");
  else if (keyname.EndsWith("_par_4"))
    keyname.ReplaceAll("_par_4","_4 MvCCRES");
  else if (keyname.EndsWith("_par_5"))
    keyname.ReplaceAll("_par_5","_5 MaNCRES");
  else if (keyname.EndsWith("_par_6"))
    keyname.ReplaceAll("_par_6","_6 MvNCRES");
  else if (keyname.EndsWith("_par_7"))
    keyname.ReplaceAll("_par_7","_7 Theta_Delta2Npi");
  else if (keyname.EndsWith("_par_8"))
    keyname.ReplaceAll("_par_8","_8 AhtBY");
  else if (keyname.EndsWith("_par_9"))
    keyname.ReplaceAll("_par_9","_9 BhtBY");
  else if (keyname.EndsWith("_par_10"))
    keyname.ReplaceAll("_par_10","_10 CV1uBY");
  else if (keyname.EndsWith("_par_11"))
    keyname.ReplaceAll("_par_11","_11 CV2uBY");
  else if (keyname.EndsWith("_par_12"))
    keyname.ReplaceAll("_par_12","_12 FrCEx_pi");
  else if (keyname.EndsWith("_par_13"))
    keyname.ReplaceAll("_par_13","_13 FrElas_pi");
  else if (keyname.EndsWith("_par_14"))
    keyname.ReplaceAll("_par_14","_14 FrInel_pi");
  else if (keyname.EndsWith("_par_15"))
    keyname.ReplaceAll("_par_15","_15 FrAbs_pi");
  else if (keyname.EndsWith("_par_16"))
    keyname.ReplaceAll("_par_16","_16 FrPiProd_pi");
  else if (keyname.EndsWith("_par_17"))
    keyname.ReplaceAll("_par_17","_17 FrCEx_N");
  else if (keyname.EndsWith("_par_18"))
    keyname.ReplaceAll("_par_18","_18 FrElas_N");
  else if (keyname.EndsWith("_par_19"))
    keyname.ReplaceAll("_par_19","_19 FrInel_N");
  else if (keyname.EndsWith("_par_20"))
    keyname.ReplaceAll("_par_20","_20 FrAbs_N");
  else if (keyname.EndsWith("_par_21"))
    keyname.ReplaceAll("_par_21","_21 FrPiProd_N");
  else if (keyname.EndsWith("_par_22"))
    keyname.ReplaceAll("_par_22","_22 E2p2h_A_nu");
  else if (keyname.EndsWith("_par_23"))
    keyname.ReplaceAll("_par_23","_23 E2p2h_B_nu");
  else if (keyname.EndsWith("_par_24"))
    keyname.ReplaceAll("_par_24","_24 E2p2h_A_nubar");
  else if (keyname.EndsWith("_par_25"))
    keyname.ReplaceAll("_par_25","_25 E2p2h_B_nubar");
  else if (keyname.EndsWith("_par_26"))
    keyname.ReplaceAll("_par_26","_26 C12ToAr40_2p2hScaling_nu");
  else if (keyname.EndsWith("_par_27"))
    keyname.ReplaceAll("_par_27","_27 C12ToAr40_2p2hScaling_nubar");
  else if (keyname.EndsWith("_par_28"))
    keyname.ReplaceAll("_par_28","_28 NR_nu_n_CC_2pi");
  else if (keyname.EndsWith("_par_29"))
    keyname.ReplaceAll("_par_29","_29 NR_nu_n_CC_3pi");
  else if (keyname.EndsWith("_par_30"))
    keyname.ReplaceAll("_par_30","_30 NR_nu_p_CC_2pi");
  else if (keyname.EndsWith("_par_31"))
    keyname.ReplaceAll("_par_31","_31 NR_nu_p_CC_3pi");
  else if (keyname.EndsWith("_par_32"))
    keyname.ReplaceAll("_par_32","_32 NR_nu_np_CC_1pi");
  else if (keyname.EndsWith("_par_33"))
    keyname.ReplaceAll("_par_33","_33 NR_nu_n_NC_1pi");
  else if (keyname.EndsWith("_par_34"))
    keyname.ReplaceAll("_par_34","_34 NR_nu_n_NC_2pi");
  else if (keyname.EndsWith("_par_35"))
    keyname.ReplaceAll("_par_35","_35 NR_nu_n_NC_3pi");
  else if (keyname.EndsWith("_par_36"))
    keyname.ReplaceAll("_par_36","_36 NR_nu_p_NC_1pi");
  else if (keyname.EndsWith("_par_37"))
    keyname.ReplaceAll("_par_37","_37 NR_nu_p_NC_2pi");
  else if (keyname.EndsWith("_par_38"))
    keyname.ReplaceAll("_par_38","_38 NR_nu_p_NC_3pi");
  else if (keyname.EndsWith("_par_39"))
    keyname.ReplaceAll("_par_39","_39 NR_nubar_n_CC_1pi");
  else if (keyname.EndsWith("_par_40"))
    keyname.ReplaceAll("_par_40","_40 NR_nubar_n_CC_2pi");
  else if (keyname.EndsWith("_par_41"))
    keyname.ReplaceAll("_par_41","_41 NR_nubar_n_CC_3pi");
  else if (keyname.EndsWith("_par_42"))
    keyname.ReplaceAll("_par_42","_42 NR_nubar_p_CC_1pi");
  else if (keyname.EndsWith("_par_43"))
    keyname.ReplaceAll("_par_43","_43 NR_nubar_p_CC_2pi");
  else if (keyname.EndsWith("_par_44"))
    keyname.ReplaceAll("_par_44","_44 NR_nubar_p_CC_3pi");
  else if (keyname.EndsWith("_par_45"))
    keyname.ReplaceAll("_par_45","_45 NR_nubar_n_NC_1pi");
  else if (keyname.EndsWith("_par_46"))
    keyname.ReplaceAll("_par_46","_46 NR_nubar_n_NC_2pi");
  else if (keyname.EndsWith("_par_47"))
    keyname.ReplaceAll("_par_47","_47 NR_nubar_n_NC_3pi");
  else if (keyname.EndsWith("_par_48"))
    keyname.ReplaceAll("_par_48","_48 NR_nubar_p_NC_1pi");
  else if (keyname.EndsWith("_par_49"))
    keyname.ReplaceAll("_par_49","_49 NR_nubar_p_NC_2pi");
  else if (keyname.EndsWith("_par_50"))
    keyname.ReplaceAll("_par_50","_50 NR_nubar_p_NC_3pi");
  else if (keyname.EndsWith("_par_51"))
    keyname.ReplaceAll("_par_51","_51 BeRPA_A");
  else if (keyname.EndsWith("_par_52"))
    keyname.ReplaceAll("_par_52","_52 BeRPA_B");
  else if (keyname.EndsWith("_par_53"))
    keyname.ReplaceAll("_par_53","_53 BeRPA_D");
  else if (keyname.EndsWith("_par_54"))
    keyname.ReplaceAll("_par_54","_54 nuenuebar_xsec_ratio");
  else if (keyname.EndsWith("_par_55"))
    keyname.ReplaceAll("_par_55","_55 nuenumu_xsec_ratio");

  return keyname;
}

TString FindXsecValName(TString keyname, TString keyname2)
{
  if (keyname.Contains("_par_0_"))
    keyname2 += "_MaCCQE_-3";
  else if (keyname.Contains("_par_1_"))
    keyname2 += "_VecFFCCQEshape_-3";
  else if (keyname.Contains("_par_2_"))
    keyname2 += "_CCQEPauliSupViaKF_-3";
  else if (keyname.Contains("_par_3_"))
    keyname2 += "_MaCCRES_-3";
  else if (keyname.Contains("_par_4_"))
    keyname2 += "_MvCCRES_-3";
  else if (keyname.Contains("_par_5_"))
    keyname2 += "_MaNCRES_-3";
  else if (keyname.Contains("_par_6_"))
    keyname2 += "_MvNCRES_-3";
  else if (keyname.Contains("_par_7_"))
    keyname2 += "_Theta_Delta2Npi_-3";
  else if (keyname.Contains("_par_8_"))
    keyname2 += "_AhtBY_-3";
  else if (keyname.Contains("_par_9_"))
    keyname2 += "_BhtBY_-3";
  else if (keyname.Contains("_par_10"))
    keyname2 += "_CV1uBY_-3";
  else if (keyname.Contains("_par_11"))
    keyname2 += "_CV2uBY_-3";
  else if (keyname.Contains("_par_12"))
    keyname2 += "_FrCEx_pi_-3";
  else if (keyname.Contains("_par_13"))
    keyname2 += "_FrElas_pi_-3";
  else if (keyname.Contains("_par_14"))
    keyname2 += "_FrInel_pi_-3";
  else if (keyname.Contains("_par_15"))
    keyname2 += "_FrAbs_pi_-3";
  else if (keyname.Contains("_par_16"))
    keyname2 += "_FrPiProd_pi_-3";
  else if (keyname.Contains("_par_17"))
    keyname2 += "_FrCEx_N_-3";
  else if (keyname.Contains("_par_18"))
    keyname2 += "_FrElas_N_-3";
  else if (keyname.Contains("_par_19"))
    keyname2 += "_FrInel_N_-3";
  else if (keyname.Contains("_par_20"))
    keyname2 += "_FrAbs_N_-3";
  else if (keyname.Contains("_par_21"))
    keyname2 += "_FrPiProd_N_-3";
  else if (keyname.Contains("_par_22"))
    keyname2 += "_E2p2h_A_nu_-3";
  else if (keyname.Contains("_par_23"))
    keyname2 += "_E2p2h_B_nu_-3";
  else if (keyname.Contains("_par_24"))
    keyname2 += "_E2p2h_A_nubar_-3";
  else if (keyname.Contains("_par_25"))
    keyname2 += "_E2p2h_B_nubar_-3";
  else if (keyname.Contains("_par_26"))
    keyname2 += "_C12ToAr40_2p2hScaling_nu_-3";
  else if (keyname.Contains("_par_27"))
    keyname2 += "_C12ToAr40_2p2hScaling_nubar_-3";
  else if (keyname.Contains("_par_28"))
    keyname2 += "_NR_nu_n_CC_2Pi_-3";
  else if (keyname.Contains("_par_29"))
    keyname2 += "_NR_nu_n_CC_3Pi_-3";
  else if (keyname.Contains("_par_30"))
    keyname2 += "_NR_nu_p_CC_2Pi_-3";
  else if (keyname.Contains("_par_31"))
    keyname2 += "_NR_nu_p_CC_3Pi_-3";
  else if (keyname.Contains("_par_32"))
    keyname2 += "_NR_nu_np_CC_1Pi_-3";
  else if (keyname.Contains("_par_33"))
    keyname2 += "_NR_nu_n_NC_1Pi_-3";
  else if (keyname.Contains("_par_34"))
    keyname2 += "_NR_nu_n_NC_2Pi_-3";
  else if (keyname.Contains("_par_35"))
    keyname2 += "_NR_nu_n_NC_3Pi_-3";
  else if (keyname.Contains("_par_36"))
    keyname2 += "_NR_nu_p_NC_1Pi_-3";
  else if (keyname.Contains("_par_37"))
    keyname2 += "_NR_nu_p_NC_2Pi_-3";
  else if (keyname.Contains("_par_38"))
    keyname2 += "_NR_nu_p_NC_3Pi_-3";
  else if (keyname.Contains("_par_39"))
    keyname2 += "_NR_nubar_n_CC_1Pi_-3";
  else if (keyname.Contains("_par_40"))
    keyname2 += "_NR_nubar_n_CC_2Pi_-3";
  else if (keyname.Contains("_par_41"))
    keyname2 += "_NR_nubar_n_CC_3Pi_-3";
  else if (keyname.Contains("_par_42"))
    keyname2 += "_NR_nubar_p_CC_1Pi_-3";
  else if (keyname.Contains("_par_43"))
    keyname2 += "_NR_nubar_p_CC_2Pi_-3";
  else if (keyname.Contains("_par_44"))
    keyname2 += "_NR_nubar_p_CC_3Pi_-3";
  else if (keyname.Contains("_par_45"))
    keyname2 += "_NR_nubar_n_NC_1Pi_-3";
  else if (keyname.Contains("_par_46"))
    keyname2 += "_NR_nubar_n_NC_2Pi_-3";
  else if (keyname.Contains("_par_47"))
    keyname2 += "_NR_nubar_n_NC_3Pi_-3";
  else if (keyname.Contains("_par_48"))
    keyname2 += "_NR_nubar_p_NC_1Pi_-3";
  else if (keyname.Contains("_par_49"))
    keyname2 += "_NR_nubar_p_NC_2Pi_-3";
  else if (keyname.Contains("_par_50"))
    keyname2 += "_NR_nubar_p_NC_3Pi_-3";
  else if (keyname.Contains("_par_51"))
    keyname2 += "_BeRPA_A_-3";
  else if (keyname.Contains("_par_52"))
    keyname2 += "_BeRPA_B_-3";
  else if (keyname.Contains("_par_53"))
    keyname2 += "_BeRPA_D_-3";
  else if (keyname.Contains("_par_54"))
    keyname2 += "_nuenuebar_xsec_ratio_-3";
  else if (keyname.Contains("_par_55"))
    keyname2 += "_nuenumu_xsec_ratio_-3";

  return keyname2;
}


TString FindLIBAnaName(TString keyname, TString keyname3)
{
  if (keyname.Contains("_par_0_"))
    keyname3 += "wgt_MaCCQE";
  else if (keyname.Contains("_par_1_"))
    keyname3 += "wgt_VecFFCCQEshape";
  else if (keyname.Contains("_par_2_"))
    keyname3 += "wgt_CCQEPauliSupViaKF";
  else if (keyname.Contains("_par_3_"))
    keyname3 += "wgt_MaCCRES";
  else if (keyname.Contains("_par_4_"))
    keyname3 += "wgt_MvCCRES";
  else if (keyname.Contains("_par_5_"))
    keyname3 += "wgt_MaNCRES";
  else if (keyname.Contains("_par_6_"))
    keyname3 += "wgt_MvNCRES";
  else if (keyname.Contains("_par_7_"))
    keyname3 += "wgt_Theta_Delta2Npi";
  else if (keyname.Contains("_par_8_"))
    keyname3 += "wgt_AhtBY";
  else if (keyname.Contains("_par_9_"))
    keyname3 += "wgt_BhtBY";
  else if (keyname.Contains("_par_10"))
    keyname3 += "wgt_CV1uBY";
  else if (keyname.Contains("_par_11"))
    keyname3 += "wgt_CV2uBY";
  else if (keyname.Contains("_par_12"))
    keyname3 += "wgt_FrCEx_pi";
  else if (keyname.Contains("_par_13"))
    keyname3 += "wgt_FrElas_pi";
  else if (keyname.Contains("_par_14"))
    keyname3 += "wgt_FrInel_pi";
  else if (keyname.Contains("_par_15"))
    keyname3 += "wgt_FrAbs_pi";
  else if (keyname.Contains("_par_16"))
    keyname3 += "wgt_FrPiProd_pi";
  else if (keyname.Contains("_par_17"))
    keyname3 += "wgt_FrCEx_N";
  else if (keyname.Contains("_par_18"))
    keyname3 += "wgt_FrElas_N";
  else if (keyname.Contains("_par_19"))
    keyname3 += "wgt_FrInel_N";
  else if (keyname.Contains("_par_20"))
    keyname3 += "wgt_FrAbs_N";
  else if (keyname.Contains("_par_21"))
    keyname3 += "wgt_FrPiProd_N";
  else if (keyname.Contains("_par_22"))
    keyname3 += "wgt_E2p2h_A_nu";
  else if (keyname.Contains("_par_23"))
    keyname3 += "wgt_E2p2h_B_nu";
  else if (keyname.Contains("_par_24"))
    keyname3 += "wgt_E2p2h_A_nubar";
  else if (keyname.Contains("_par_25"))
    keyname3 += "wgt_E2p2h_B_nubar";
  else if (keyname.Contains("_par_26"))
    keyname3 += "wgt_C12ToAr40_2p2hScaling_nu";
  else if (keyname.Contains("_par_27"))
    keyname3 += "wgt_C12ToAr40_2p2hScaling_nubar";
  else if (keyname.Contains("_par_28"))
    keyname3 += "wgt_NR_nu_n_CC_2Pi";
  else if (keyname.Contains("_par_29"))
    keyname3 += "wgt_NR_nu_n_CC_3Pi";
  else if (keyname.Contains("_par_30"))
    keyname3 += "wgt_NR_nu_p_CC_2Pi";
  else if (keyname.Contains("_par_31"))
    keyname3 += "wgt_NR_nu_p_CC_3Pi";
  else if (keyname.Contains("_par_32"))
    keyname3 += "wgt_NR_nu_np_CC_1Pi";
  else if (keyname.Contains("_par_33"))
    keyname3 += "wgt_NR_nu_n_NC_1Pi";
  else if (keyname.Contains("_par_34"))
    keyname3 += "wgt_NR_nu_n_NC_2Pi";
  else if (keyname.Contains("_par_35"))
    keyname3 += "wgt_NR_nu_n_NC_3Pi";
  else if (keyname.Contains("_par_36"))
    keyname3 += "wgt_NR_nu_p_NC_1Pi";
  else if (keyname.Contains("_par_37"))
    keyname3 += "wgt_NR_nu_p_NC_2Pi";
  else if (keyname.Contains("_par_38"))
    keyname3 += "wgt_NR_nu_p_NC_3Pi";
  else if (keyname.Contains("_par_39"))
    keyname3 += "wgt_NR_nubar_n_CC_1Pi";
  else if (keyname.Contains("_par_40"))
    keyname3 += "wgt_NR_nubar_n_CC_2Pi";
  else if (keyname.Contains("_par_41"))
    keyname3 += "wgt_NR_nubar_n_CC_3Pi";
  else if (keyname.Contains("_par_42"))
    keyname3 += "wgt_NR_nubar_p_CC_1Pi";
  else if (keyname.Contains("_par_43"))
    keyname3 += "wgt_NR_nubar_p_CC_2Pi";
  else if (keyname.Contains("_par_44"))
    keyname3 += "wgt_NR_nubar_p_CC_3Pi";
  else if (keyname.Contains("_par_45"))
    keyname3 += "wgt_NR_nubar_n_NC_1Pi";
  else if (keyname.Contains("_par_46"))
    keyname3 += "wgt_NR_nubar_n_NC_2Pi";
  else if (keyname.Contains("_par_47"))
    keyname3 += "wgt_NR_nubar_n_NC_3Pi";
  else if (keyname.Contains("_par_48"))
    keyname3 += "wgt_NR_nubar_p_NC_1Pi";
  else if (keyname.Contains("_par_49"))
    keyname3 += "wgt_NR_nubar_p_NC_2Pi";
  else if (keyname.Contains("_par_50"))
    keyname3 += "wgt_NR_nubar_p_NC_3Pi";
  else if (keyname.Contains("_par_51"))
    keyname3 += "wgt_BeRPA_A";
  else if (keyname.Contains("_par_52"))
    keyname3 += "wgt_BeRPA_B";
  else if (keyname.Contains("_par_53"))
    keyname3 += "wgt_BeRPA_D";
  else if (keyname.Contains("_par_54"))
    keyname3 += "wgt_nuenuebar_xsec_ratio";
  else if (keyname.Contains("_par_55"))
    keyname3 += "wgt_nuenumu_xsec_ratio";

  return keyname3;
}

void makeSigmaVarComparisons(TString inputfile, TString valfile)
{
  std::cout << "honk" << std::endl;
  gStyle->SetOptStat(0);
   TFile* file = new TFile(inputfile);
   TList* list = file->GetListOfKeys();
   
   TFile* file2 = new TFile(valfile);
   //TFile* file3 = new TFile("LIBAna_splines_1sig_test.root");
   
   TH1D* fhc_nue_nom = (TH1D*)file->Get("fhc_nue_osc");
   fhc_nue_nom->SetLineColor(kBlack);
   fhc_nue_nom->GetXaxis()->SetTitle("E_{reconstructed} (GeV)");
   TH1D* fhc_numu_nom = (TH1D*)file->Get("fhc_numu_osc");
   fhc_numu_nom->SetLineColor(kBlack);
   fhc_numu_nom->GetXaxis()->SetTitle("E_{reconstructed} (GeV)");
   TH1D* rhc_nue_nom = (TH1D*)file->Get("rhc_nue_osc");
   rhc_nue_nom->SetLineColor(kBlack);
   rhc_nue_nom->GetXaxis()->SetTitle("E_{reconstructed} (GeV)");
   TH1D* rhc_numu_nom = (TH1D*)file->Get("rhc_numu_osc");
   rhc_numu_nom->SetLineColor(kBlack);
   rhc_numu_nom->GetXaxis()->SetTitle("E_{reconstructed} (GeV)");
   
   TH1D* fhc_nue_nom2 = (TH1D*)file2->Get("FD_FHC_nue/FD_FHC_Nue_unosc_total_flux_Nov17_0_0");
   fhc_nue_nom2->SetLineColor(kBlack);
   fhc_nue_nom2->SetLineStyle(2);
   TH1D* fhc_numu_nom2 = (TH1D*)file2->Get("FD_FHC_numu/FD_FHC_Numu_unosc_total_flux_Nov17_0_0");
   fhc_numu_nom2->SetLineColor(kBlack);
   fhc_numu_nom2->SetLineStyle(2);
   TH1D* rhc_nue_nom2 = (TH1D*)file2->Get("FD_RHC_nue/FD_RHC_Nue_unosc_total_flux_Nov17_0_0");
   rhc_nue_nom2->SetLineColor(kBlack);
   rhc_nue_nom2->SetLineStyle(2);
   TH1D* rhc_numu_nom2 = (TH1D*)file2->Get("FD_RHC_numu/FD_RHC_Numu_unosc_total_flux_Nov17_0_0");
   rhc_numu_nom2->SetLineColor(kBlack);
   rhc_numu_nom2->SetLineStyle(2);

   std::cout << "honk2" << std::endl;
   std::cout << list->GetEntries() << "  number of Hists" << std::endl;

   TCanvas* c0 = new TCanvas("c0","c0",0,0,700,900);
   c0->Divide(1,2);
   c0->Print("sigmavar_all.pdf[");
   for(int i=0; i<list->GetEntries(); i++)
   {
     //if(i == 2120 ) {break;}
     TString keyname = list->At(i)->GetName();
     TString keyname2;
     TString keyname3;
      if(keyname=="fhc_nue_osc" || keyname=="fhc_numu_osc" 
	 || keyname=="rhc_nue_osc" || keyname=="rhc_numu_osc" || !keyname.Contains("n3"))
	 continue;
      TH1D* divisor;
      TH1D* divisor2;
      if(keyname.Contains("fhc_numu")) {
	divisor=fhc_numu_nom;
	divisor2=fhc_numu_nom2; 
        keyname2 = "FD_FHC_numu/FD_FHC_Numu_unosc_total"; }
      if(keyname.Contains("fhc_nue")) {
	divisor=fhc_nue_nom;
	divisor2=fhc_nue_nom2;
        keyname2 = "FD_FHC_nue/FD_FHC_Nue_unosc_total"; }
      if(keyname.Contains("rhc_numu")) {
	divisor=rhc_numu_nom;
	divisor2=rhc_numu_nom2;
        keyname2 = "FD_RHC_numu/FD_RHC_Numu_unosc_total"; }
      if(keyname.Contains("rhc_nue")) {
	divisor=rhc_nue_nom;
	divisor2=rhc_nue_nom2;
        keyname2 = "FD_RHC_nue/FD_RHC_Nue_unosc_total"; }
      
      keyname2 = FindXsecValName(keyname, keyname2);
      keyname3 = FindLIBAnaName(keyname, keyname3);



      std::cout << keyname << " || " << keyname2 << " || " << keyname3 << std::endl; 

      //TH1D* lib = (TH1D*)file3->Get(keyname3);
      
      TH1D* hVar[4];
      TH1D* hVar2[4];
      hVar[0]=(TH1D*)file->Get(keyname);
      hVar2[0]=(TH1D*)file2->Get(keyname2);
      keyname.ReplaceAll("n3","n1");
      keyname2.ReplaceAll("_-3","_-1");
      hVar[1]=(TH1D*)file->Get(keyname);
      hVar2[1]=(TH1D*)file2->Get(keyname2);
      keyname.ReplaceAll("n1","p1");
      keyname2.ReplaceAll("_-1","_1");
      hVar[2]=(TH1D*)file->Get(keyname);
      hVar2[2]=(TH1D*)file2->Get(keyname2);
      keyname.ReplaceAll("p1","p3");
      keyname2.Chop();
      keyname2+="3";
      hVar[3]=(TH1D*)file->Get(keyname);
      hVar2[3]=(TH1D*)file2->Get(keyname2);

      hVar[0]->SetLineColor(kRed);
      hVar[1]->SetLineColor(kMagenta);
      hVar[2]->SetLineColor(kCyan);
      hVar[3]->SetLineColor(kBlue);
      //lib->SetLineColor(kYellow);
      
      hVar2[0]->SetLineColor(kRed);
      hVar2[0]->SetLineStyle(2);
      hVar2[1]->SetLineColor(kMagenta);
      hVar2[1]->SetLineStyle(2);
      hVar2[2]->SetLineColor(kCyan);
      hVar2[2]->SetLineStyle(2);
      hVar2[3]->SetLineColor(kBlue);
      hVar2[3]->SetLineStyle(2);


      keyname.ReplaceAll("_sig_p3","");
      
      // if xsec, find the name of the parameter
      if(keyname.Contains("xsec"))
	keyname = getXsecParName(keyname);
      divisor->SetTitle(keyname);
      
      c0->cd(1);
      
      // Set the axis range to max value of all variations
      std::vector<double> maxlist;
      maxlist.push_back(divisor->GetMaximum());
      for(int j=0; j<4; j++) {
         maxlist.push_back(hVar[j]->GetMaximum());
         maxlist.push_back(hVar2[j]->GetMaximum()); }
      double maxmax = *std::max_element(maxlist.begin(), maxlist.end());
      divisor->GetYaxis()->SetRangeUser(0,1.05*maxmax);
      

     divisor->DrawCopy("HIST");
     divisor2->DrawCopy("HIST SAME");
     
     if(keyname.Contains("E2p2h")) {
       for(int j=1; j<3; j++) { 
	   hVar[j]->DrawCopy("HIST SAME");
	   hVar2[j]->DrawCopy("HIST SAME"); }
    } 
    else {
       for(int j=0; j<4; j++) { 
	   hVar[j]->DrawCopy("HIST SAME");
	   hVar2[j]->DrawCopy("HIST SAME"); }
    }
     
      c0->Update();

      TH1D* divcopy = (TH1D*)divisor->Clone("divcopy");
      divcopy->SetTitle("Response (Sigma/Nominal)");
      for(int j=0; j<4; j++)
      {
         hVar[j]->Divide(divcopy);
         hVar2[j]->Divide(divisor2);
      }
 
      // Set the axis range to max value of all variations
      std::vector<double> maxlist2;
      std::vector<double> minlist2;
      

     if(keyname.Contains("E2p2h")) {
       for(int j=1; j<3; j++) {
         maxlist2.push_back(hVar[j]->GetMaximum()); 
         maxlist2.push_back(hVar2[j]->GetMaximum()); 
         //maxlist2.push_back(lib->GetMaximum()); 
         minlist2.push_back(hVar[j]->GetMinimum()); 
         minlist2.push_back(hVar2[j]->GetMinimum());} 
         //minlist2.push_back(lib->GetMinimum()); }
      double maxmax2 = *std::max_element(maxlist2.begin(), maxlist2.end());
      double minmin2 = *std::min_element(minlist2.begin(), minlist2.end());
      divcopy->GetYaxis()->SetRangeUser(1.0001*minmin2,1.0001*maxmax2);}
     else {
       for(int j=0; j<4; j++) {
         maxlist2.push_back(hVar[j]->GetMaximum()); 
         maxlist2.push_back(hVar2[j]->GetMaximum()); 
         //maxlist2.push_back(lib->GetMaximum()); 
         minlist2.push_back(hVar[j]->GetMinimum()); 
         minlist2.push_back(hVar2[j]->GetMinimum());} 
         //minlist2.push_back(lib->GetMinimum()); }
      double maxmax2 = *std::max_element(maxlist2.begin(), maxlist2.end());
      double minmin2 = *std::min_element(minlist2.begin(), minlist2.end());
      divcopy->GetYaxis()->SetRangeUser(1.0001*minmin2,1.0001*maxmax2);}
      
      c0->cd(2);
      divcopy->Draw("HIST");
     if(keyname.Contains("E2p2h")) {
       for(int j=1; j<3; j++) { 
         //lib->DrawCopy("HIST SAME");}
	   hVar[j]->DrawCopy("HIST SAME");
	   hVar2[j]->DrawCopy("HIST SAME"); }
    } 
    else {
       for(int j=0; j<4; j++) { 
         //lib->DrawCopy("HIST SAME");}
	   hVar[j]->DrawCopy("HIST SAME");
	   hVar2[j]->DrawCopy("HIST SAME"); }
    }
      
     c0->Update();

      TLegend* leg = new TLegend(0.65,0.55,0.8,0.8);
      leg->AddEntry(divisor,"MaCh3 nominal","L");
      leg->AddEntry(divisor2,"CAFANA nominal","L");
      leg->AddEntry(hVar[0],"-3#sigma MaCh3","L");
      leg->AddEntry(hVar2[0],"-3#sigma CAFAna","L");
      leg->AddEntry(hVar[1],"-1#sigma MaCh3","L");
      leg->AddEntry(hVar2[1],"-1#sigma CAFAna","L");
      leg->AddEntry(hVar[2],"+1#sigma MaCh3","L");
      leg->AddEntry(hVar2[2],"+1#sigma CAFAna","L");
      leg->AddEntry(hVar[3],"+3#sigma MaCh3","L");
      leg->AddEntry(hVar2[3],"+3#sigma CAFAna","L");
      leg->SetFillColor(0);
      c0->cd(1);
      leg->Draw("SAME");
      c0->Update();

      TLegend* leg2 = new TLegend(0.65,0.55,0.8,0.8);
      leg2->AddEntry(hVar[0],"-3#sigma MaCh3","L");
      leg2->AddEntry(hVar2[0],"-3#sigma CAFAna","L");
      leg2->AddEntry(hVar[1],"-1#sigma MaCh3","L");
      leg2->AddEntry(hVar2[1],"-1#sigma CAFAna","L");
      leg2->AddEntry(hVar[2],"+1#sigma MaCh3","L");
      leg2->AddEntry(hVar2[2],"+1#sigma CAFAna","L");
      leg2->AddEntry(hVar[3],"+3#sigma MaCh3","L");
      leg2->AddEntry(hVar2[3],"+3#sigma CAFAna","L");
      //leg2->AddEntry(lib,"+1#sigma LIBAna","L");
      leg2->SetFillColor(0);
      c0->cd(2);
      leg2->Draw("SAME");
      
      c0->Update();
//      getchar();
      c0->Print("sigmavar_all.pdf");
      std::cout << " DONE." << std::endl; 

   }

   c0->Print("sigmavar_all.pdf]");

 }

