/*
	 Draw Fig.10 in PHYSICAL REVIEW C 88, 025501 (2013).
	 This program can only be run offline, using ROOT version >6, since no SetFillColorAlpha can be used in ROOT 5.

	 Weiran, Apr. 1, 2018.

*/

#include "TGraph.h"
#include "Eigen/Dense"
#include <iostream>
#include "TMath.h"
#include <complex>
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TAttFill.h"

using namespace Eigen;

int main()
{
	int npoints = 1000;
	double x[npoints], y[npoints], ex[npoints], ey[npoints];

	Matrix3d R; //SNO correlation matrix for survival probability parameterization.
	/* Pee_day = c0 + c1(Enu-10MeV) + c2(Enu-10MeV)^2 */

							/*   c0        c1        c2    */
	R << /* c0 */  1.000,   -0.299,   -0.366,
			 /* c1 */ -0.299,    1.000,   -0.206,
			 /* c2 */ -0.366,   -0.206,    1.000;

	/* Fit result and uncertainty. */
	double c[3] = {0.317, 0.0039, -0.0010};
	double errc[3] = {sqrt(0.016*0.016+0.009*0.009), sqrt(0.0066*0.0066+0.0045*0.0045), sqrt(0.0029*0.0029+0.0015*0.0015)};

	/* Covariance matrix. */
	Matrix3d Cov = Matrix3d::Zero();
	for(int i=0; i<3; i++) for(int j=0; j<3; j++)
	{
		Cov(i,j) = R(i,j) * errc[i] * errc[j];
	}

	for(int i=0; i<npoints; i++)
	{
		Vector3d v;
		x[i] = 4 + 11./npoints * i; //4~15 MeV
		v(0) = 1; v(1) = x[i] - 10; v(2) = v(1) * v(1); 
		ex[i] = 11. / npoints /2.;
		y[i] = c[0]*v(0) + c[1]*v(1) + c[2]*v(2);
		ey[i] = sqrt(v.transpose() * Cov * v);
	}

	TGraph* graph = new TGraph(npoints, x, y);
	TGraphErrors* err = new TGraphErrors(npoints, x, y, ex, ey);

	TCanvas* band = new TCanvas();
	graph->GetYaxis()->SetRangeUser(0,0.6);
	graph->GetXaxis()->SetRangeUser(4,15);
	graph->SetLineColor(2);
	graph->Draw("aplsame");
	err->SetFillColorAlpha(6,0.5);
	err->SetFillStyle(3005);
	err->Draw("3");
	band->SaveAs("this.png");
	return 0;
}


