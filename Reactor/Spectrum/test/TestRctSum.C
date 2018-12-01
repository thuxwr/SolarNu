#include "../RctSum.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TApplication.h"

int main()
{
  cout<<"Test reactor neutrino background spectrum at each site"<<endl;
  //gROOT->ProcessLine(".include ${EIGEN}"); 
  //gROOT->ProcessLine(".L ../RctSum.C+");

  TApplication* pApp =  new TApplication("appKey",0,0);

  RctSum *RctBkg = new RctSum;

  // All raw neutrino spectra
  TCanvas *c1 = new TCanvas("c1","c1");
  c1->Divide(3,2);
  for( int Exp=0; Exp<NExp_RctSum; Exp++ )  {
    c1->cd(Exp+1);
    RctBkg->NueSpec[Exp]->Draw();
  }

  // Visiable energy through IBD process
  cout<<"IBD events:"<<endl;
  TCanvas *c2 = new TCanvas("c2","c2");
  c2->Divide(3,2);
  for( int Exp=0; Exp<NExp_RctSum; Exp++ )  {
    c2->cd(Exp+1);
    RctBkg->GetVisableIBD(Exp)->Draw();
    cout<<"IBD events of "<<RctBkg->ExpName[Exp]<<" "
	<<RctBkg->GetVisableIBD(Exp)->Integral() * RctBkg->GetVisableIBD(Exp)->GetBinWidth(1)<<endl;
  }
  c2->Print("IBD.ps");

  // Test that with Daya Bay assumption
  TCanvas *c3 = new TCanvas("c3","c3");
  // The above output has the unit of #/proton/year
  TH1D* hDYB = new TH1D( *(RctBkg->GetVisableIBD(3)) );
  // Daya Bay 7.169*10^25 hydrogens/kg * 20*10^3 kg, 1 day
  hDYB->Scale( 7.169e25 * 20e3 * 24*60*60 );
  hDYB->Draw();
  hDYB->GetYaxis()->SetTitle("N IBD/MeV/day");
  hDYB->GetXaxis()->SetTitle("MeV");

  cout<<"Daya Bay can detect "<<hDYB->Integral()*hDYB->GetBinWidth(1)<<" IBD event per day."<<endl;;

  // Visible kinetic energy spectrum
  cout<<"Elastic scattering events"<<endl;
  TCanvas *c4 = new TCanvas("c4","c4");
  c4->Divide(3,2);
  for( int Exp=0; Exp<NExp_RctSum; Exp++ )  {
    c4->cd(Exp+1);
    RctBkg->GetVisableElas(Exp)->Draw();
    cout<<"Elastic events of "<<RctBkg->ExpName[Exp]<<" "
	<<RctBkg->GetVisableElas(Exp)->Integral() * RctBkg->GetVisableElas(Exp)->GetBinWidth(1)<<endl;
  }
  c4->Print("Elas.ps");

  // Save to file for other uses
  TFile f("RctBkgSum.root", "recreate");
  for( int Exp=0; Exp<NExp_RctSum; Exp++ )  {
    RctBkg->NueSpec[Exp]->Write();    
    RctBkg->GetVisableIBD(Exp)->Write();    
    RctBkg->GetVisableElas(Exp)->Write();
  }
  f.Write();

  /* Just waiting for a key hit */
  TCanvas *empty = new TCanvas;
  empty->WaitPrimitive();

  return 1;
}
