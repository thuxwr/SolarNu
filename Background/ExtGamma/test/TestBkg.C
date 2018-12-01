{
  cout<<""<<endl;
  gROOT->ProcessLine(".L ../ExtGammaSpec.C+");

  ExtGammaSpec ExtTl208("ExtTl208", 0.4);

  TCanvas c("Background", "Background spectra", 600, 600);
  c->SetLeftMargin(0.20);
  c->SetTopMargin(0.20);

  ExtTl208.hVisTotal->Draw();
  ExtTl208.hVisTotal->SetLineColor(kBlue);
  ExtTl208.hVisTotal->SetLineWidth(2);
  ExtTl208.hVisTotal->GetXaxis()->SetTitle("Visible energy [MeV]");

  ExtTl208.hVisTotal->GetYaxis()->SetTitle("Strength [Arbitary unit]");
  ExtTl208.hVisTotal->GetYaxis()->SetTitleOffset(2.1);
  
  c.Print("ExtGamSpec.ps");
}
