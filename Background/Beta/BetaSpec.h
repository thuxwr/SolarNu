//
// Extracted from Daya Bay simulation package.
// Gang Li and Zhe Wang verified and corrected. Jun. 4, 2014
//
#ifndef BETASPEC_H
#define BETASPEC_H

#include <string>
#include "TString.h"
#include "TComplex.h"

#define MAXBRANCH 200
#define MAXGAMMA  40

class TH1D;
class TGraph;

class BetaSpec {
public:
    BetaSpec();
    BetaSpec(std::string name, int z, int a, int charge=-1, bool useHuber=false);
    virtual ~BetaSpec();
    
    void SetIsotope(std::string name, int z, int a, int charge);
    void LoadSpectra(std::string filename);
    void MakeBetaSpectra();
    void MakeTree(std::string filename="", int nEvents=1000000);
    
    // correction functions
    double Fermi(int z, double Te, int charge=-1);
    double FiniteSizeEM(int z, double Te, int charge=-1);
    double FiniteSizeWI(int z, double Te, double Q0, int charge=-1);
    double Screening(int z, double Te, int charge=-1);
    double WeakMagnetism(double Te, int charge=-1);
    bool UseForbiddenCorrection(int deltaSpin, int deltaParity);
    double Forbiddenness(int z, double Te, double Q0, int deltaSpin, int deltaParity, int charge=-1);
    double ForbiddennessHuber(double Te, double Q0, int deltaSpin);
    
    
    void Draw();
    TGraph* Draw(TString option);
    void CompareSpectra();
    void CompareSize();
    void Glance();
    void Align(TGraph* g);
    
    static TComplex CGamma(TComplex z);
    static int WeightedChoice(double* weights, int size);
    
    std::string isotope;
    int Z;
    int A;
    int C;  // charge (-1 for electron, 1 for positron, 0 EC) (EC is not implemented)
    int Zdaughter;
    double Qmax;    // maximum energy level difference
    double Emax;    // maximum (energy level difference + all gamma energies)
    double R;       // nucleus charge radius
    
    double me;
    double alpha;
    TGraph* gScreeningV;
    bool useForbiddennessHuber;
    std::string baseDir;
    
    int nBranch;    // number of branches
    double Q[MAXBRANCH];    // Q-value of branches, size nBranch
    double BR[MAXBRANCH];   // branching ratio.
    int dSpin[MAXBRANCH];   // delta spin
    int dParity[MAXBRANCH]; // delta parity
    int nGamma[MAXBRANCH];  // number of correlated gammas
    double eGamma[MAXBRANCH][MAXGAMMA]; //energy of the correlated gammas, size (nBranch x nGamma[br_index])
    TH1D* hBeta[MAXBRANCH]; // the corrected beta spectra, size nBranch, binwidth 1keV, Integral()==1
    TH1D* hBetaNoCorrection[MAXBRANCH];  // the un-corrected beta spectra (no Fermi/Screening/etc.) 
    
    
    // Visable energy: Sum beta and related gammas. If beta+, add 2*0.511.
    // Unfortunately, beta+ and EC cannot be distinguished here.
    // All gammas and betas energy are added linearly. No nonlineaity considered.
    // binwidth 1 keV, Integral() == BR
    TH1D* hVis[MAXBRANCH];
    // Sum all branches together. Add 2*0.511 if it is beta+ decay.
    // binwidth 1keV, Integral() == 1
    TH1D* hVisTotal;
};

#endif
