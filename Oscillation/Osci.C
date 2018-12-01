#include "Osci.h"
#include <TMath.h>
#include <iostream>

using namespace Eigen;

Oscillation::Oscillation(int generation)
{
	gene = generation;
	IsNSI = 0;
	Fermion = -1;
}

void Oscillation::SetupParameters(double sin_2_t12, double sin_2_t13, double sin_2_t23, double sin_2_t14, double ms12, double ms13, double ms14, double cp_phase)
{
	s12 = sqrt(sin_2_t12);
	c12 = sqrt(1-sin_2_t12);
	s13 = sqrt(sin_2_t13);
	c13 = sqrt(1-sin_2_t13);
	s23 = sqrt(sin_2_t23);
	c23 = sqrt(1-sin_2_t23);
	s14 = sqrt(sin_2_t14);
	c14 = sqrt(1-sin_2_t14);

	/* Get Vacuum mixing matrix. */
	Uv = MatrixXcd::Zero(gene, gene);
	Mv = MatrixXcd::Zero(gene, gene);
	if(gene==2) 
	{
		Uv <<       c12 ,  s12,
						 		-s12,  c12;
		Mv(1,1) += ms12;
	}

	if(gene==3)
	{
		Uv <<       c12*c13                                     ,  s12*c13                                     ,  std::polar(s13, -cp_phase),
						 		-s12*c23 - std::polar(c12*s13*s23, cp_phase),  c12*c23 - std::polar(s12*s13*s23, cp_phase) ,  c13*s23                   ,
								s12*s23 - std::polar(c12*s13*c23, cp_phase) ,  -c12*s23 - std::polar(s12*s13*c23, cp_phase),  c13*c23                   ;
		Mv(1,1) += ms12;
		Mv(2,2) += ms13;
	}

	if(gene==4)
	{
		Uv <<       c12*c13                                     ,  s12*c13                                     ,  std::polar(s13, -cp_phase) ,  0,
						 		-s12*c23 - std::polar(c12*s13*s23, cp_phase),  c12*c23 - std::polar(s12*s13*s23, cp_phase) ,  c13*s23                    ,  0,
								s12*s23 - std::polar(c12*s13*c23, cp_phase) ,  -c12*s23 - std::polar(s12*s13*c23, cp_phase),  c13*c23                    ,  0,
								0                                           ,  0                                           ,  0                          ,  1;
		MatrixXcd Rotate = MatrixXcd::Identity(4,4);
		/* Assuming no 14, 24 and 34 CP violation, since they are not detectable till now. */
		Rotate << c14 , 0 , 0 , s14 ,
					    0   , 1 , 0 , 0   ,
							0   , 0 , 1 , 0   ,
							-s14, 0 , 0 , c14 ;
		Uv = Uv * Rotate;
		Mv(1,1) += ms12;
		Mv(2,2) += ms13;
		Mv(3,3) += ms14;
	}

	A = Uv * Mv * Uv.transpose().conjugate();
}

void Oscillation::SetupNSI(double eeel, double eeer, double etel, double eter, double eeu, double eed, double euu, double eud, double etu, double etd, double uuu, double uud, double utu, double utd, double ttu, double ttd)
{
	IsNSI = 1;
	/* Should first setup vacuum matrix, using results from KamLAND and others, without all solar nu experiments. */
	SetupParameters(0.325, 0.0241, 0.5, 0, 7.54e-5, 2.5e-3, 0, 0); //KamLAND-only.
	if(gene==4) 
	{
		std::cout << "Sterile degenerates with NSI. Please decrease dof and try again." << std::endl;
		exit(0);
	}
	if(gene==2)
	{
		/* To be finished. */
		exit(0);
	}

	/* All parameters. UUe, UTe, TTe do not affect oscillation. EUe has been constrained to almost 0. */
	EEe = eeel + eeer;	EEu = eeu;	EEd = eed;	EUu = euu;	EUd = eud;	ETe = etel + eter;
	ETu = etu;	ETd = etd;	UUu = uuu;	UUd = uud;	UTu = utu;	UTd = utd;	TTu = ttu;	TTd = ttd;
}

void Oscillation::SetupNSI(double sin_2_t12, double ms12, int fermion, double ed, double en)
{
	if(gene!=2)
	{
		std::cout << "Use parameterization for 2nu. Change number of generations." << std::endl;
		exit(0);
	}
	IsNSI = 2;
	SetupParameters(sin_2_t12, 0.022, 0.5, 0, ms12);
	Fermion = fermion;
	Ed = ed; En = en;
}


MatrixXcd Oscillation::GetVacuumMixingMatrix()
{
	return Uv;
}

MatrixXcd Oscillation::GetMatterMixingMatrix(double Energy, double ne, double nn, double np)
{
	MatrixXcd Ham = Hamiltonian(Energy, ne, nn, np);
	MatrixXcd B = Ham * 2 * Energy;
	solver.compute(B);
	MatrixXcd Um = solver.eigenvectors();
	return Um;
}

double Oscillation::AdiabaticPropagation(double Energy, double ne, double nn, double np)
{
	MatrixXcd Um = GetMatterMixingMatrix(Energy, ne, nn, np);
//	return (init.transpose() * Um * Uv.transpose().conjugate()).transpose();
	std::complex<double> AM[3];
	for(int mass=0; mass<gene; mass++) AM[mass] = Um(0, mass);
	double prob = 0;
	for(int mass=0; mass<gene; mass++) prob += pow(abs(AM[mass]),2) * pow(abs(Uv(0,mass)),2);
	if(IsNSI==2) { //P_ee = cos^4\theta_13*P_{2\nu} + sin^4\theta_13
		double sin_2_t13 = 0.0241;
		double P2nu = prob;
		prob = pow(1-sin_2_t13, 2) * P2nu + pow(sin_2_t13, 2);
	}
	//std::cout << Um << std::endl;
	//std::cout << Uv << std::endl;
	//return Uv * Um.transpose().conjugate() * init;
	return prob;
}


MatrixXcd Oscillation::Hamiltonian(double Energy, double ne, double nn, double np)
{
	double Acc = 2.5348e-31 /*2sqrt(2)G_{F}*/ * ne * Energy;
	if(IsNSI==2) Acc *= (1-0.022);
	double Anc = -0.5 * 2.5348e-31 * nn * Energy;
	MatrixXcd Ham = A;
	Ham(0,0) = Ham(0,0) + Acc + Anc;
	Ham(1,1) += Anc;
	if(gene>2) Ham(2,2) += Anc;
	Ham /= (2*Energy);
	if(IsNSI==1)
	{
		double nu = 2 * np + nn;
		double nd = 2 * nn + np; 
		MatrixXcd C = MatrixXcd::Zero(3,3);
		double Const = 1.2674e-31; //sqrt(2)*G_F
		C(0,0) = ne*EEe + nu*EEu + nd*EEd;
		C(0,1) =          nu*EUu + nd*EUd;
		C(0,2) = ne*ETe + nu*ETu + nd*ETd;
		C(1,1) =          nu*UUu + nd*UUd;
		C(1,2) =          nu*UTu + nd*UTd;
		C(2,2) =          nu*TTu + nd*TTd;
		C(1,0) = C(0,1); C(2,0) = C(0,2); C(2,1) = C(1,2);
		C = Const * C;
		Ham = Ham + C; //Add NSI contribution.
	}
	if(IsNSI==2)
	{
		double nu = 2 * np + nn;
		double nd = 2 * nn + np;
		MatrixXcd C = MatrixXcd::Zero(2,2);
		double Const = 1.2674e-31;
		double N = 0;
		if(Fermion==0) N = ne;
		else if(Fermion==1) N = nu;
		else if(Fermion==2) N = nd;
		else {
			std::cout << "An error occured due to unsupported NSI scheme. " << std::endl;
			exit(0);
		}

		C(0,0) = N*(-Ed); C(0,1) = N*En; C(1,0) = N*En; C(1,1) = N*Ed;
		C = Const * C;
		Ham = Ham + C;
	}

	return Ham;
}

double Oscillation::SurvProb(double Energy, double ne, double nn)
{
	MatrixXcd Um = GetMatterMixingMatrix(Energy, ne, nn);
	double prob = 0;
	for(int i=0; i<3; i++)
	{
		prob += pow(std::abs(Um(0,i)),2) * pow(std::abs(Uv(0,i)),2);
	}
	return prob;
}

double Oscillation::GetTheta12(double Energy, double ne, double nn, double np)
{
	MatrixXcd Mm = GetMatterMixingMatrix(Energy, ne, nn, np);
	if(gene==2) return asin(abs(Mm(0,1)));
	if(gene==3 || gene==4) return atan(abs(Mm(0,1))/abs(Mm(0,0))); //Note: if gene=4 this is not precise. However, for those points where 1 and 4 do not mix maximally, this is enough.
	else
	{
		std::cout << "Generation is incorrect. Only support 2, 3 or 4 generations. " << std::endl;
		return 0;
	}
}

double Oscillation::GetTheta13(double Energy, double ne, double nn, double np)
{
	MatrixXcd Mm = GetMatterMixingMatrix(Energy, ne, nn, np);
	if(gene==2)
	{
		std::cout << "Theta13 is unavailable for 2 flavor neutrinos. " << std::endl;
		return 0;
	}
	if(gene==3 || gene==4) return asin(abs(Mm(0,2)));
	else
	{
		std::cout << "Generation is incorrect. Only support 2, 3 or 4 generations. " << std::endl;
		return 0;
	}
}

double Oscillation::GetTheta14(double Energy, double ne, double nn)
{
	MatrixXcd Mm = GetMatterMixingMatrix(Energy, ne, nn);
	if(gene==2 || gene==3)
	{
		std::cout << "Theta14 is unavailable for no sterile neutrinos. " << std::endl;
		return 0;
	}
	if(gene==4) return asin(abs(Mm(3,0)));
	else
	{
		std::cout << "Generation is incorrect. Only support 2, 3 or 4 generations. " << std::endl;
		return 0;
	}
}

double Oscillation::GetMass12(double Energy, double ne, double nn, double np)
{
	MatrixXcd Ham = Hamiltonian(Energy, ne, nn, np);
	MatrixXcd B = Ham * 2 * Energy;
	solver.compute(B);
	return solver.eigenvalues()[1].real() - solver.eigenvalues()[0].real();
}

double Oscillation::GetMass13(double Energy, double ne, double nn, double np)
{
	MatrixXcd Ham = Hamiltonian(Energy, ne, nn, np);
	MatrixXcd B = Ham * 2 * Energy;
	solver.compute(B);
	return solver.eigenvalues()[2].real() - solver.eigenvalues()[0].real();
}
