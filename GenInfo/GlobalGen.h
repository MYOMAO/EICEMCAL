#include "TVector3.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"




float GammaPX1 = 0;
float GammaPY1 = 0;
float GammaPZ1 = 0;


float GammaPX2 = 0;
float GammaPY2 = 0;
float GammaPZ2 = 0;


float GammaPX[2]; 
float GammaPY[2];
float GammaPZ[2];
float GammaE[2];

int iGamma;
int nElectrons;

double TotalENERGY = 0; 
int NSUCKPAR = 0;

int NMissMass = 0;
double MissMassEnergy = 0;
double etotPre = 999;



float OutSideE;

int LayerID;
int BlackHole;
int IsActive;
int EMCAL;
int BHLAYER;


TCanvas * c = new TCanvas("c","c",600,600);

TH1D * ZPos = new TH1D("ZPos","",5000,-1000,1000);
TH2D * XYPos = new TH2D("XYPos","",100,-1000,1000,100,-1000,1000);


TFile * finDebug; 
TTree * DebugTree; 

