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



void DebugEvtAna(int Energy, int BeamOpt, int ScinOpt){
	
	


	gStyle->SetOptStat(0);
	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();

	TString BeamName;


	if(BeamOpt == 0) BeamName = "pi0";
	if(BeamOpt == 1) BeamName = "gamma";	
	if(BeamOpt == 2) BeamName = "e";
	if(BeamOpt == 3) BeamName = "mu";

	TString BeamNameFull;

	if(BeamOpt == 0) BeamNameFull = "#pi^{0}";
	if(BeamOpt == 1) BeamNameFull = "#gamma";	
	if(BeamOpt == 2) BeamNameFull = "e";
	if(BeamOpt == 3) BeamNameFull = "#mu";


	float BeamMass;

	if(BeamOpt == 0) BeamMass = 0;
	if(BeamOpt == 1) BeamMass = 0;	
	if(BeamOpt == 2) BeamMass = 0.00051;
	if(BeamOpt == 3) BeamMass = 0.1057;



	TLatex *lat = new TLatex();
	lat->SetNDC();
	lat->SetTextSize(0.035);

	float Ratio;
	float IncEnergy = Energy;

	TString infile = Form("DebugEvt/ShowerInfo_%s.root",BeamName.Data());

	if(ScinOpt == 1) infile = Form("DebugEvt/ShowerInfo_%s_New.root",BeamName.Data());

	
	float BHBLedep;
	float BHFLedep;
	float FEMCALedep;

	int evtid = 0;
	int evtidPre = 0;


	int NBins = 100; 
	float EWidthDown = 0.2;
	float EWidthUp = 0.1;
	


	TFile * fin = new TFile(infile.Data());
	fin->cd();
	TTree * ShowerInfo = (TTree *) fin->Get("ShowerInfo");
	ShowerInfo->SetBranchAddress("BHBLedep",&BHBLedep);
	ShowerInfo->SetBranchAddress("BHFLedep",&BHFLedep);
	ShowerInfo->SetBranchAddress("FEMCALedep",&FEMCALedep);

	ShowerInfo->SetBranchAddress("evtid",&evtid);

	int NEntries = ShowerInfo->GetEntries();

	float TotalEnergyBefore = 0;

	float TotalEnergyAfter = 0;

	TH1D * EnergyHisBefore = new TH1D("EnergyHisBefore","",NBins,Energy-EWidthDown,Energy+EWidthUp);
	EnergyHisBefore->GetXaxis()->SetTitle("Total Energy Deposited to the EMCAL (GeV)");
	EnergyHisBefore->GetYaxis()->SetTitle("Number of Events");
	EnergyHisBefore->GetXaxis()->CenterTitle();
	EnergyHisBefore->GetYaxis()->CenterTitle();
	EnergyHisBefore->GetYaxis()->SetTitleOffset(1.4);
	
	TH1D * EnergyHisAfter = new TH1D("EnergyHisAfter","",NBins,Energy-EWidthDown,Energy+EWidthUp);
	EnergyHisAfter->GetXaxis()->SetTitle("Total Energy After Addition (GeV)");
	EnergyHisAfter->GetYaxis()->SetTitle("Number of Events");
	EnergyHisAfter->GetXaxis()->CenterTitle();
	EnergyHisAfter->GetYaxis()->CenterTitle();
	EnergyHisAfter->GetYaxis()->SetTitleOffset(1.4);
		

	TotalEnergyBefore = BeamMass;
	for(int i = 0; i < NEntries; i++){
		
		ShowerInfo->GetEntry(i);



		
		//if(evtid > 10) break;
		if(evtidPre != evtid){

			
			
			EnergyHisBefore->Fill(TotalEnergyBefore);
			EnergyHisAfter->Fill(TotalEnergyAfter);

			cout << "evtidPre = " << evtidPre << "    TotalEnergyAfter = " << TotalEnergyAfter << endl;

			TotalEnergyBefore = BeamMass;
			TotalEnergyAfter =  BeamMass;
	
		
		}
		evtidPre = evtid;

			TotalEnergyBefore = FEMCALedep  + TotalEnergyBefore;
			TotalEnergyAfter = BHFLedep + FEMCALedep + BHBLedep + TotalEnergyAfter;
			
	//		cout << "TotalEnergyAfter = " << TotalEnergyAfter << endl;
	}

	//Ratio = TotalEnergy/Energy;
	
	//cout << "TotalEnergyAfter = " << TotalEnergyAfter << endl;

//    std::ofstream DebugBro2("DatFiles/DebugEvt2.dat",ios::app);
//	DebugBro2  << "Index = " << Index << "    FUCKIN TotalEnergy = " << TotalEnergy << "   Energy = " << Energy << "   Ratio = " << Ratio << endl;
	
	if(ScinOpt == 0){
	EnergyHisBefore->Draw();
	lat->DrawLatex(0.30,0.80,Form("Incident %s Eneregy = %.1f GeV",BeamNameFull.Data(),IncEnergy));	
	lat->DrawLatex(0.30,0.75,"Standard PHENIX Pb Shashlik EMCAL");	

	c->SaveAs(Form("DebugEvtAnaPlots/EnergyHisBefore_%d_%s.png",Energy,BeamName.Data()));
	
	EnergyHisAfter->Draw();
	lat->DrawLatex(0.30,0.80,Form("Incident %s Eneregy = %.1f GeV",BeamNameFull.Data(),IncEnergy));	
	lat->DrawLatex(0.30,0.75,"Standard PHENIX Pb Shashlik EMCAL");	

	c->SaveAs(Form("DebugEvtAnaPlots/EnergyHisAfter_%d_%s.png",Energy,BeamName.Data()));

	}

	if(ScinOpt == 1){
	EnergyHisBefore->Draw();
	lat->DrawLatex(0.18,0.80,Form("Incident %s Eneregy = %.1f GeV",BeamNameFull.Data(),IncEnergy));	
	lat->DrawLatex(0.18,0.75,"Pb Shashlik EMCAL 1.5 mm Pb and 1.0 mm Scintillator");	

	c->SaveAs(Form("DebugEvtAnaPlots/EnergyHisBefore_%d_%s_New.png",Energy,BeamName.Data()));
	
	EnergyHisAfter->Draw();
	lat->DrawLatex(0.18,0.80,Form("Incident %s Eneregy = %.1f GeV",BeamNameFull.Data(),IncEnergy));	
	lat->DrawLatex(0.18,0.75,"Pb Shashlik EMCAL 1.5 mm Pb and 1.0 mm Scintillator");	

	c->SaveAs(Form("DebugEvtAnaPlots/EnergyHisAfter_%d_%s_New.png",Energy,BeamName.Data()));

	}

}
