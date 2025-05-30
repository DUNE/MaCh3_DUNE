void makePiePlot(bool wRC=true, int H=0) //H=0 for both, H=1 for NH, H=-1 for IH
{
//   double bds[6] = {2.80,-2.95,-2.39,-1.13,-0.5,0.13};//3 sig, 2 sig, 1 sig, 1 sig, 2 sig, 3 sig 2019
//   double bestf=-1.74; //best fit 2019
   TString plotname;
   double bds[6];
   double bestf;
   TString bfst;
   if(wRC && H==0)
   {
      bds[0]=2.545; //3 sig
      bds[1]=-TMath::Pi(); //2 sig (MUST BE A NEGATIVE NUMBER)
      bds[2]=-2.545; //1 sig
      bds[3]=-1.037; //1 sig
      bds[4]=-0.346; //2 sig
      bds[5]=0.346; //3 sig
      bestf=-1.97; //best fit 2019 wRC
      bfst ="-1.97";
      plotname="circleplotwRC_both";
   }
   else if(!wRC && H==0)
   {
      bds[0]=2.549; //3 sig
      bds[1]=-3.7338; //2 sig (MUST BE A NEGATIVE NUMBER)
      bds[2]=-2.796; //1 sig
      bds[3]=-0.723; //1 sig
      bds[4]=0.220; //2 sig
      bds[5]=0.220; //3 sig
      bestf=-2.09; //best fit 2019 woRC
      bfst ="-2.09";
      plotname="circleplotwoRC_both";
   }
   else if(wRC && H==1)
   {
      bds[0]=2.51327; //3 sig
      bds[1]=-3.2044; //2 sig (MUST BE A NEGATIVE NUMBER)
      bds[2]=-2.670; //1 sig
      bds[3]=-1.162; //1 sig
      bds[4]=-0.3456; //2 sig
      bds[5]=0.4084; //3 sig
      bestf=-2.09; //best fit 2019 woRC
      bfst ="-2.09";
      plotname="circleplotwRC_NH";
   }
   else if(!wRC && H==1)
   {
      bds[0]=2.325; //3 sig
      bds[1]=-3.958; //2 sig (MUST BE A NEGATIVE NUMBER)
      bds[2]=-TMath::Pi(); //1 sig
      bds[3]=-1.13097; //1 sig
      bds[4]=0.2199; //2 sig
      bds[5]=0.2199; //3 sig
      bestf=-2.35; //best fit 2019 woRC
      bfst ="-2.35";
      plotname="circleplotwoRC_NH";
   }
   else if(wRC && H==-1)
   {
      bds[0]=3.14159; //3 sig
      bds[1]=-2.545; //2 sig (MUST BE A NEGATIVE NUMBER)
      bds[2]=-1.9792; //1 sig
      bds[3]=-0.8796; //1 sig
      bds[4]=-0.37699; //2 sig
      bds[5]=0.125664; //3 sig
      bestf=-1.45958; //best fit 2019 woRC
      bfst ="-1.46";
      plotname="circleplotwRC_IH";
   }
   else if(!wRC && H==-1)
   {
      bds[0]=3.3929; //3 sig
      bds[1]=-2.89027; //2 sig (MUST BE A NEGATIVE NUMBER)
      bds[2]=-2.042; //1 sig
      bds[3]=-0.503; //1 sig
      bds[4]=0.1885; //2 sig
      bds[5]=0.1885; //3 sig
      bestf=-1.21; //best fit 2019 woRC
      bfst ="-1.21";
      plotname="circleplotwoRC_IH";
   }


   double radius=0.4;

   TEllipse* onesig = new TEllipse(0.5,0.5,radius,radius,bds[2]*180/3.14159,bds[4]*180/3.14159);
   TEllipse* twosigA = new TEllipse(0.5,0.5,radius,radius,bds[1]*180/3.14159,bds[2]*180/3.14159);
   TEllipse* twosigB = new TEllipse(0.5,0.5,radius,radius,bds[3]*180/3.14159,bds[4]*180/3.14159);
   TEllipse* threesigA = new TEllipse(0.5,0.5,radius,radius,bds[0]*180/3.14159,(bds[1]+2*3.14159)*180/3.14159);
   TEllipse* threesigB = new TEllipse(0.5,0.5,radius,radius,bds[4]*180/3.14159,bds[5]*180/3.14159);
   TEllipse* rest = new TEllipse(0.5,0.5,radius,radius,bds[5]*180/3.14159,bds[0]*180/3.14159);

   onesig->SetFillColor(13);
   twosigA->SetFillColor(12);
   twosigB->SetFillColor(12);
   threesigA->SetFillColor(11);
   threesigB->SetFillColor(11);
   TLine* l = new TLine(0.5-radius,0.5,0.5+radius,0.5);
   l->SetLineWidth(3);
   
   TLine* l2 = new TLine(0.5,0.5-radius,0.5,0.5+radius);
   l2->SetLineWidth(3);

   TArrow* bf = new TArrow(0.5,0.5,0.5+radius*cos(bestf),0.5+radius*sin(bestf),0.04,"|>");
   bf->SetLineWidth(3);
   bf->SetLineColor(kRed);
   bf->SetFillColor(kRed);
      
   TCanvas* c = new TCanvas("c","",0,0,1000,1000);
   onesig->Draw();
   twosigA->Draw();
   twosigB->Draw();
   threesigA->Draw();
   threesigB->Draw();
   rest->Draw();
   l->Draw();
   l2->Draw();
   bf->Draw();
      

   TLegend* leg = new TLegend(0.0,0.8,0.23,0.95);
   leg->AddEntry(bf,"Best Fit","L");
   leg->AddEntry(onesig,"1#sigma","F");
   leg->AddEntry(twosigA,"2#sigma","F");
   if(wRC)
      leg->AddEntry(threesigA,"3#sigma","F");
   leg->Draw();

   TText *t0 = new TText(0.5+radius+0.02,0.5,"0");
   t0->SetTextAlign(22);
   t0->SetTextColor(kBlack);
   t0->SetTextFont(43);
   t0->SetTextSize(40);
   t0->SetTextAngle(0);
   t0->Draw();

   TLatex *tp = new TLatex(0.5-radius-0.02,0.5,"#pi");
   tp->SetTextAlign(22);
   tp->SetTextColor(kBlack);
   tp->SetTextFont(43);
   tp->SetTextSize(40);
   tp->SetTextAngle(0);
   tp->Draw();

   TLatex *tp2 = new TLatex(0.5,0.5+radius+0.04,"#frac{#pi}{2}");
   tp2->SetTextAlign(22);
   tp2->SetTextColor(kBlack);
   tp2->SetTextFont(43);
   tp2->SetTextSize(40);
   tp2->SetTextAngle(0);
   tp2->Draw();

   TLatex *tmp2 = new TLatex(0.5,0.5-radius-0.04,"-#frac{#pi}{2}");
   tmp2->SetTextAlign(22);
   tmp2->SetTextColor(kBlack);
   tmp2->SetTextFont(43);
   tmp2->SetTextSize(40);
   tmp2->SetTextAngle(0);
   tmp2->Draw();



   TLatex *tbf = new TLatex(0.5+(radius+0.02)*cos(bestf),0.5+(radius+0.02)*sin(bestf),bfst);
   tbf->SetTextAlign(22);
   tbf->SetTextColor(kRed);
   tbf->SetTextFont(43);
   tbf->SetTextSize(40);
   tbf->SetTextAngle(0);
   tbf->Draw();

   c->SaveAs(plotname+".pdf");
   c->SaveAs(plotname+".png");


}
