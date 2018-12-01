void Draw()
{
	TFile* file = TFile::Open("ClAbsorpCS.root");

	TGraph* g = (TGraph*)file->Get("CrossSection");

	g->GetXaxis()->SetRangeUser(0,20);
	g->GetXaxis()->SetTitle("Energy [MeV]");
	g->GetXaxis()->SetTitleSize(0.045);

	g->GetYaxis()->SetTitle("Cross Section [#times 10^{-46} cm^{2}]");
	g->GetYaxis()->SetTitleSize(0.045);

	g->SetTitle("");

	g->SetLineWidth(2);

	g->Draw("al");

}
