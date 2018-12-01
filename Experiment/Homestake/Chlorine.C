/*
	 Test constraint from Homestake.
	 Use model prediction from both GS and AGS.

	 Weiran, Mar. 19, 2018.

*/

#include "../../Target/ClCapture.h"
#include "../../Solar/SolarNu.h"
#include "TH2D.h"
#include "TMath.h"
#include <iostream>
#include <string>

using namespace std;

/* See Astrophys.J. 496 (1998) 505-526 */
double ClCenter = 2.56, ClError = 0.23;

int main()
{
	ClCapture capt;
	SolarNu solar;
	int nbinstheta = 200;
	int nbinsmass = 200;
	for(int model=1; model<=2; model++)
	{
		string histname = solar.GetModelName(model) + "Cl";
		TH2D* ClContour = new TH2D(histname.c_str(), histname.c_str(), nbinstheta, -4, 1, nbinsmass, -12, -3);
		double flux[9], fluxerr[9], hflux[9];
		for(int comp=1; comp<=9; comp++) 
		{
			flux[comp-1] = solar.GetFlux(model, comp);
			fluxerr[comp-1] = solar.GetError(model, comp);
			hflux[comp-1] = flux[comp-1] + fluxerr[comp-1];
		}
		for(int binx=1; binx<=nbinstheta; binx++) for(int biny=1; biny<=nbinsmass; biny++)
		{
			double log_tan_2_theta = ClContour->GetXaxis()->GetBinCenter(binx);
			double tan_2_theta = pow(10, log_tan_2_theta);
			double sin_2_theta = tan_2_theta / (1 + tan_2_theta);
			double log_ms = ClContour->GetYaxis()->GetBinCenter(biny);
			double ms = pow(10, log_ms);

			double ClCSErr, FluxErr, TotErr;
			double ClCap;
			capt.SetCapture(sin_2_theta, ms);
			ClCap = capt.GetFlux(flux);
			ClCSErr = capt.GetError(flux);

			capt.SetCapture(sin_2_theta, ms);
			FluxErr = capt.GetFlux(hflux)-ClCap;

			TotErr = sqrt(ClCSErr*ClCSErr + FluxErr*FluxErr);

			/* Calculate chi2. */
			double sigma = sqrt(TotErr*TotErr + ClError*ClError);
			double chi2 = pow((ClCap-ClCenter)/sigma,2);
			ClContour->SetBinContent(binx, biny, chi2);
		}
		ClContour->GetXaxis()->SetTitle("log_{10}tan^{2}#theta_{12}");
		ClContour->GetYaxis()->SetTitle("log_{10}#Delta m^{2}_{12} [eV^{2}]");
		ClContour->SetStats(kFALSE);
		string filename = "./result/" + histname + "Contour.root";
		ClContour->SaveAs(filename.c_str());
	}

	return 0;
}








