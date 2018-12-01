/*
	 Test constraint from all gallium experiments including SAGE/GALLEX and GNO.
	 Use model prediction from both GS and AGS.

	 Weiran, Mar. 19, 2018.

*/

#include "../../Target/GaCapture.h"
#include "../../Solar/SolarNu.h"
#include "TH2D.h"
#include "TMath.h"
#include <iostream>
#include <string>

using namespace std;

/* Combined result from SAGE, GALLEX and GNO gives 66.1 \pm 3.1 SNU, Phys. Rev., C80:015807, 2009. */
double GaCenter = 66.1, GaError = 3.1;

int main()
{
	GaCapture capt;
	SolarNu solar;
	int nbinstheta = 200;
	int nbinsmass = 200;
	for(int model=1; model<=2; model++)
	{
		string histname = solar.GetModelName(model) + "Ga";
		TH2D* GaContour = new TH2D(histname.c_str(), histname.c_str(), nbinstheta, -4, 1, nbinsmass, -12, -3);
		double flux[9], fluxerr[9], hflux[9];
		for(int comp=1; comp<=9; comp++) 
		{
			flux[comp-1] = solar.GetFlux(model, comp);
			fluxerr[comp-1] = solar.GetError(model, comp);
			hflux[comp-1] = flux[comp-1] + fluxerr[comp-1];
		}
		for(int binx=1; binx<=nbinstheta; binx++) for(int biny=1; biny<=nbinsmass; biny++)
		{
			double log_tan_2_theta = GaContour->GetXaxis()->GetBinCenter(binx);
			double tan_2_theta = pow(10, log_tan_2_theta);
			double sin_2_theta = tan_2_theta / (1 + tan_2_theta);
			double log_ms = GaContour->GetYaxis()->GetBinCenter(biny);
			double ms = pow(10, log_ms);

			double GaCSErr, FluxErr, TotErr;
			double GaCap;
			capt.SetCapture(sin_2_theta, ms);
			GaCap = capt.GetFlux(flux);
			GaCSErr = capt.GetError(flux);

			capt.SetCapture(sin_2_theta, ms);
			FluxErr = capt.GetFlux(hflux)-GaCap;

			TotErr = sqrt(GaCSErr*GaCSErr + FluxErr*FluxErr);

			/* Calculate chi2. */
			double sigma = sqrt(TotErr*TotErr + GaError*GaError);
			double chi2 = pow((GaCap-GaCenter)/sigma,2);
			GaContour->SetBinContent(binx, biny, chi2);
		}
		GaContour->GetXaxis()->SetTitle("log_{10}tan^{2}#theta_{12}");
		GaContour->GetYaxis()->SetTitle("log_{10}#Delta m^{2}_{12} [eV^{2}]");
		GaContour->SetStats(kFALSE);
		string filename = "./result/" + histname + "Contour.root";
		GaContour->SaveAs(filename.c_str());
	}

	return 0;
}








