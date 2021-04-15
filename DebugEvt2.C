#include "TROOT.h"
#include "TH1.h"
#include "TTree.h"
#include "TH2.h"
#include "TAxis.h"
#include "TF1.h"
#include "TFile.h"
#include "TMath.h"
#include "TSystem.h"
#include "TVector2.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom.h"
#include <iostream>
#include <fstream>


using namespace std;

using std::cout;
using std::endl;



void DebugEvt2(int Index){
	
	
	gStyle->SetOptStat(0);

	float Ratio;
	float Energy = 1;

	TString infile = "DebugEvt/out.root";
	
	float BHBLedep;
	float BHFLedep;
	float FEMCALedep;
	int evtid;
	float evtidfloat;

	TFile * fin = new TFile(infile.Data());
	fin->cd();
	TTree * ShowerInfo = (TTree *) fin->Get("ShowerInfo");
	ShowerInfo->SetBranchAddress("BHBLedep",&BHBLedep);
	ShowerInfo->SetBranchAddress("BHFLedep",&BHFLedep);
	ShowerInfo->SetBranchAddress("FEMCALedep",&FEMCALedep);
	ShowerInfo->SetBranchAddress("evtid",&evtidfloat);

	int NEntries = ShowerInfo->GetEntries();

	float TotalEnergy = 0;

	for(int i = 0; i < NEntries; i++){
		
		evtid = int(evtidfloat);

		ShowerInfo->GetEntry(i);
		
		TotalEnergy = TotalEnergy + FEMCALedep + BHBLedep + BHFLedep;

	

	}

	Ratio = TotalEnergy/Energy;

    std::ofstream DebugBro2("DatFiles/DebugEvt2.dat",ios::app);
	DebugBro2  << "Index = " << Index << "    FUCKIN TotalEnergy = " << TotalEnergy << "   Energy = " << Energy << "   Ratio = " << Ratio << endl;




}
