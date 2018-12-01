{
  cout<<""<<endl;
  gROOT->ProcessLine(".L ../BetaSpec.C+");

  BetaSpec Kr85 ("Kr85", 36,85, -1);
  BetaSpec Bi210("Bi210",83,210,-1);
  BetaSpec C11  ("C11",  6, 11,  1);
  BetaSpec C14  ("C14",  6, 14, -1);
  BetaSpec C10  ("C10",  6, 10,  1);
  BetaSpec Tl208("Tl208",81,208,-1);
  BetaSpec Be11 ("Be11", 4, 11, -1);

  TCanvas c("Background", "Background spectra", 750, 750);

  c.Divide(3,3);
  TVirtualPad * pPad;
  
  pPad = c.cd(1);
  pPad->SetLeftMargin(0.15);
  Kr85.Draw();

  pPad = c.cd(2);
  pPad->SetLeftMargin(0.15);
  Bi210.Draw();

  /// About C11 decay, what is its Q-value? wiki (1 MeV) or nndc (2 MeV)?
  pPad = c.cd(3);
  pPad->SetLeftMargin(0.15);
  C11.Draw();

  pPad = c.cd(4);
  pPad->SetLeftMargin(0.15);
  C14.Draw();

  pPad = c.cd(5);
  pPad->SetLeftMargin(0.15);
  C10.Draw();

  pPad = c.cd(6);
  pPad->SetLeftMargin(0.15);
  Tl208.Draw();

  pPad = c.cd(7);
  pPad->SetLeftMargin(0.15);
  Be11.Draw();

  gROOT->ProcessLine(".L ../../ExtGamma/ExtGammaSpec.C+");
  ExtGammaSpec ExtTl208("ExtTl208", 0.4);
  pPad = c.cd(8);
  pPad->SetLeftMargin(0.15);

  ExtTl208.hVisTotal->SetLineWidth(2);
  ExtTl208.hVisTotal->SetLineColor(kBlack);
  ExtTl208.hVisTotal->GetXaxis()->SetTitle("E (MeV)");
  ExtTl208.hVisTotal->GetXaxis()->SetRangeUser(0,2.8);
  ExtTl208.hVisTotal->Draw();

  double xMax = ExtTl208.hVisTotal->GetXaxis()->GetXmax();
  double yMax = ExtTl208.hVisTotal->GetMaximum();

  TLatex latex;
  latex.SetTextSize(0.08);
  latex.SetTextAlign(11);  //align at top                                                                                                    
  latex.DrawLatex( xMax*0.05, yMax*0.6, "Ext-Tl208");
  
  c.Print("BkgSpec.ps");
}
