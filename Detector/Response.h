/*
	 Give detector response function for SK I~IV, and convert true spectra to SK's detected ones.
	 All can be found on SK's website and arXiv.

	 SNO's will be added in the future.
	 Weiran, Mar. 19, 2018.

*/

#ifndef RESPONSE_H
#define RESPONSE_H

#include <string>
#include "TH1D.h"

using namespace std;

class Response
{
	public:
		Response();
		~Response(){};

		/* Experiment names are: SKI, SKII, SKIII, SKIV, (SNO). */
		void Convert(TH1D*& DetS, TH1D* TrueS, string ExpName);

		double Resolution(double Energy, string ExpName);

	private:
		const static int nGauss = 1000000;
		double gauss_cdf(double x);
		double gauss_array[nGauss];

};

#endif

