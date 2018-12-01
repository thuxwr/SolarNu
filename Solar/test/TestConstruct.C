/*
	 A very trivial test.

*/
#include <iostream>
#include "TFile.h"
#include "TGraph.h"

using namespace std;

int main()
{
	TFile* file = new TFile("./SurvProb.root","READ");
	TGraph* graph = (TGraph*)file->Get("B8daye");
	TGraph* hehe = new TGraph(*graph);
//	TGraph* hehe2 = &hehe;
	file->Close();
	hehe->SaveAs("ttt.root");
	return 1;
}


