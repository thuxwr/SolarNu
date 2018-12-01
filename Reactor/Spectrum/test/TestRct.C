{
  cout<<"Test the standard reactor neutrino spectrum"<<endl;
  gROOT->ProcessLine(".L ../RctSpec.C+");

  RctSpec *Rct = new RctSpec;
  
  TCanvas* c = new TCanvas;
  Rct->Get()->Draw();
  cout<<Rct->Get()->Integral(0,20)<<endl;

  TCanvas* c2 = new TCanvas;
  Rct->GetSpec()->Draw();
}
