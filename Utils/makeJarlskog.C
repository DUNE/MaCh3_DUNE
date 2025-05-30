{

  TFile* outFile = new TFile("/vols/dune/ljw20/fit_results/DUNE_Perl_NuPhys_NDDet/HaddedChains/Jarlskog_DUNE_NuPhys_withRC.root","RECREATE");

  TChain* post = new TChain("posteriors");
  post->Add("/vols/dune/ljw20/fit_results/DUNE_Perl_NuPhys_NDDet/HaddedChains/DUNE_NuPhys_wNDDet.root");

  unsigned int nSteps = post->GetEntries();

  std::cout << "number of mcmc steps = " << nSteps << std::endl;

  double s2th13, s2th23, s2th12, dcp, dm2, rc;
  int step;
  post->SetBranchAddress("sin2th_13",&s2th13);
  post->SetBranchAddress("sin2th_23",&s2th23);
  post->SetBranchAddress("sin2th_12",&s2th12);
  post->SetBranchAddress("delta_cp",&dcp);
  post->SetBranchAddress("delm2_23",&dm2);

  rc = 1.0;
  double RCreweight = 1.0;
  post->SetBranchAddress("step",&step);

  double prior_s2th13, prior_s2th23, prior_s2th12, prior_dcp, prior_sindcp, prior_dm32, prior_dm21;
  double prior_wRC_s2th13;

  double s13, s23, s12, sdcp, c13, c23, c12, j;
  double prior_s13, prior_s23, prior_s12, prior_sdcp, prior_c13, prior_c23, prior_c12, prior_j;
  double prior_wRC_s13, prior_wRC_c13, prior_wRC_j, prior_wRC_flatsindcp_j;

  // Variables needed for RC reweight in sin^2 theta_13 (NuFit 4.0)
  double new_central = 0.02241;
  double new_error = 0.00065;
  double  new_chi, new_prior;

  TH1D* jarl = new TH1D("jarl","jarl",1000,-0.05,0.05);
  TH2D* jarl_th23 = new TH2D("jarl_th23","jarl_th23",500,-0.05,0.05,500,0.3,0.7);
  TH2D* jarl_dcp = new TH2D("jarl_dcp","jarl_dcp",500,-0.05,0.05,500,-1.*TMath::Pi(),TMath::Pi());
  jarl->SetTitle("Jarlskog Invariant;J #equiv s_{13}c_{13}^{2}s_{12}c_{12}s_{23}c_{23}sin#delta");
  jarl_th23->SetTitle("Jarlskog Invariant;J #equiv s_{13}c_{13}^{2}s_{12}c_{12}s_{23}c_{23}sin#delta");

  TH1D* jarl_IH = (TH1D*)jarl->Clone();
  TH1D* jarl_NH = (TH1D*)jarl->Clone();
  TH2D* jarl_th23_IH = (TH2D*)jarl_th23->Clone();
  TH2D* jarl_th23_NH = (TH2D*)jarl_th23->Clone();
  TH2D* jarl_dcp_IH = (TH2D*)jarl_dcp->Clone();
  TH2D* jarl_dcp_NH = (TH2D*)jarl_dcp->Clone();

  TH1D* jarl_flatsindcp = (TH1D*)jarl->Clone();
  TH1D* jarl_RC = (TH1D*)jarl->Clone();
  TH1D* jarl_IH_flatsindcp = (TH1D*)jarl->Clone();
  TH1D* jarl_NH_flatsindcp = (TH1D*)jarl->Clone();

  TH2D* jarl_th23_flatsindcp = (TH2D*)jarl_th23->Clone();
  TH2D* jarl_th23_IH_flatsindcp = (TH2D*)jarl_th23->Clone();
  TH2D* jarl_th23_NH_flatsindcp = (TH2D*)jarl_th23->Clone();

  TCanvas* c = new TCanvas("c","c",600,600);
  c->Draw();

  // to apply a prior that is flat in sin(dcp) intead of dcp
  TF1 *prior3 = new TF1("prior3","TMath::Abs(TMath::Cos(x))");

  // T2K prior is flat (and uncorrelated) in dcp, sin^2(th13), sin^2(th23)
  TRandom3* randGen = new TRandom3();

  for(unsigned int i = 0;i<nSteps; i++) {
  //for(unsigned int i = 0;i<2000000; i++) {

    if(i%1000000==0) std::cout << "step # " << i <<std::endl;
    post->GetEntry(i);
    if(step<80000) continue; // burn-in cut
    s13 = TMath::Sqrt(s2th13);
    s23 = TMath::Sqrt(s2th23);
    s12 = TMath::Sqrt(s2th12);
    sdcp = TMath::Sin(dcp);
    c13 = TMath::Sqrt(1.-s2th13);
    c12 = TMath::Sqrt(1.-s2th12);
    c23 = TMath::Sqrt(1.-s2th23);

	new_chi = (s2th13 - new_central)/new_error;
	new_prior = std::exp(-0.5 * new_chi * new_chi);
	RCreweight = new_prior;

    j = s13*c13*c13*s12*c12*s23*c23*sdcp;

    double prior_weight = prior3->Eval(dcp);

    jarl->Fill(j,rc);
	jarl_RC->Fill(j, RCreweight);
    jarl_th23->Fill(j,s2th23,rc);
    jarl_dcp->Fill(j,dcp,rc);

    jarl_flatsindcp->Fill(j,prior_weight);
    jarl_th23_flatsindcp->Fill(j,s2th23,prior_weight*rc);

    if(dm2 > 0.) {
      jarl_NH->Fill(j,rc);
      jarl_th23_NH->Fill(j,s2th23,rc);
      jarl_dcp_NH->Fill(j,dcp,rc);
      jarl_NH_flatsindcp->Fill(j,prior_weight*rc);
      jarl_th23_NH_flatsindcp->Fill(j,s2th23,prior_weight*rc);
    }
    else if(dm2 < 0.) {
      jarl_IH->Fill(j,rc);
      jarl_th23_IH->Fill(j,s2th23,rc);
      jarl_dcp_IH->Fill(j,dcp,rc);
      jarl_IH_flatsindcp->Fill(j,prior_weight*rc);
      jarl_th23_IH_flatsindcp->Fill(j,s2th23,prior_weight*rc);
    }
  }

  outFile->cd();
  jarl->Write("jarlskog_both");
  jarl_RC->Write("jarlskog_both_RC");
  jarl_NH->Write("jarlskog_NH");
  jarl_IH->Write("jarlskog_IH");
  jarl_th23->Write("jarlskog_th23_both");
  jarl_th23_NH->Write("jarlskog_th23_NH");
  jarl_th23_IH->Write("jarlskog_th23_IH");

  jarl_dcp->Write("jarlskog_dcp_both");
  jarl_dcp_NH->Write("jarlskog_dcp_NH");
  jarl_dcp_IH->Write("jarlskog_dcp_IH");


  jarl_flatsindcp->Write("jarlskog_both_flatsindcp");
  jarl_NH_flatsindcp->Write("jarlskog_NH_flatsindcp");
  jarl_IH_flatsindcp->Write("jarlskog_IH_flatsindcp");
  jarl_th23_flatsindcp->Write("jarlskog_th23_both_flatsindcp");
  jarl_th23_NH_flatsindcp->Write("jarlskog_th23_NH_flatsindcp");
  jarl_th23_IH_flatsindcp->Write("jarlskog_th23_IH_flatsindcp");

}
