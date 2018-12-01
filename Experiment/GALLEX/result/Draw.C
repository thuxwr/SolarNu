{
	TFile *_file0 = TFile::Open("GS98GaContour.root");
	int palette[3];
	palette[0] = 38;
	palette[1] = kBlack;
	palette[2] = 10;
	gStyle->SetPalette(3, palette);

	GS98Ga->SetContour(3);
	GS98Ga->SetContourLevel(0,0);
	GS98Ga->SetContourLevel(1,5.99146);
	GS98Ga->SetContourLevel(2, 6.5);

	GS98Ga->Draw("cont4");

	c1->SaveAs("ga.pdf");

}
