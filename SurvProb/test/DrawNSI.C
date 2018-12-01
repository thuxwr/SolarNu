{
	TFile* file = TFile::Open("NSI.root");

	upNSI->SetLineColor(kViolet);
	upNSI->SetLineWidth(2);
	downNSI->SetLineColor(kOrange);
	downNSI->SetLineWidth(2);
	MSW->SetLineColor(kBlue);
	MSW->SetLineWidth(2);

	upNSI->GetXaxis()->SetRangeUser(0.1,14);
	upNSI->SetTitle("");
	upNSI->GetXaxis()->SetTitle("E_{#nu} [MeV]");
	upNSI->GetXaxis()->SetTitleSize(0.045);
	upNSI->GetYaxis()->SetTitle("P_{ee}");
	upNSI->GetYaxis()->SetTitleSize(0.045);

	TLegend* leg = new TLegend(0.1,0.1,0.4,0.4);
	leg->SetBorderSize(0);
	leg->AddEntry(MSW,"Standard","l");
	leg->AddEntry(upNSI,"NSI-up","l");
	leg->AddEntry(downNSI,"NSI-dw","l");

	upNSI->Draw("al");
	downNSI->Draw("lsame");
	MSW->Draw("lsame");
	leg->Draw("same");
}
