{
	TFile* file = TFile::Open("SK4.root");

	TH1D* spec = (TH1D*)file->Get("spectrum");
	spec->SetTitle("");
	spec->GetXaxis()->SetTitle("E_{kin} [MeV]");
	spec->GetXaxis()->SetTitleSize(0.045);
	spec->GetYaxis()->SetTitleSize(0.045);
	spec->SetLineWidth(2);

	TCanvas* c1 = new TCanvas;
	c1->DrawFrame(3.5,0.2,20,0.75);
	//c1->DrawFrame(4.5,0,20,0.85);

	TLatex latex;
	latex.SetTextSize(0.08);
	latex.SetTextAlign(12);
	latex.DrawLatex(4.4,0.63, "SK-IV");

	spec->Draw("same");

}
