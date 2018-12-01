/*
	 Neutrino flux changes if distance between sun and earth changes.
	 This header calculates seasonal change of total neutrino flux.

	 Weiran, Mar. 7, 2018.

*/
#ifndef SEASON_H
#define SEASON_H

using namespace std;

class Season
{
	public:
		Season()
		{
			SunEarthDistance = 7.58122e11;
			DistPeriAp = 0.01672;
		}
		~Season(){};

		/* I divide the whole flux into four parts according to their season. The average of these four seasons should be the same as that indicated by SSM. */
		double GetFlux(double totalflux, int season) /* Season: 1-spring, 2-summer, 3-autumn, 4-winter*/
		{





	private:
		double SunEarthDistance; //The same as in SurvProb.h, unit:MeV/eV^2
		double DistPeriAp;

