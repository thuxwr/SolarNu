{
	TFile* f1 = new TFile("./null/up1.root");
	TH2D* up1 = (TH2D*)f1->Get("JPcontour");
	TFile* f2 = new TFile("./null/up2.root");
	TH2D* up2 = (TH2D*)f2->Get("JPcontour");

	TH2D* total = (TH2D*)up1->Clone();

	for(int xbin=1; xbin<=up1->GetXaxis()->GetNbins(); xbin++) for(int ybin=1; ybin<=up1->GetYaxis()->GetNbins(); ybin++)
	{
		double normal = up1->GetBinContent(xbin, ybin);
		double invert = up2->GetBinContent(xbin, ybin);

		/* Choose a minimum. */
		double min = normal;
		if(normal>invert) min = invert;
		total->SetBinContent(xbin, ybin, min);
	}
	total->SetTitle("");
	total->GetXaxis()->SetTitle("#varepsilon_{D}^{u}");
	total->GetXaxis()->SetTitleSize(0.045);
	total->GetYaxis()->SetTitle("#varepsilon_{N}^{u}");
	total->GetYaxis()->SetTitleSize(0.045);
	total->Draw("colz");
}
