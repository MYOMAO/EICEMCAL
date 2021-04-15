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



void PlotEDepLoop(){

	gStyle->SetOptStat(0);




	TLatex *lat = new TLatex();
	lat->SetNDC();
	lat->SetTextSize(0.035);

	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();

	int evtid;
	float X1;
	float X2;
	float Y1;
	float Y2;
	float Z1;
	float Z2;
	float T1;
	float T2;
	float BHBLedep;
	float BHFLedep;
	float FEMCALedep;
	int hitid;
	int trkid;
	int showerid;
	int layer;	

	int XBin;

	float Energy;
	float EStep = 0.5;
	float EInit = 0.5;


	TH1D * EnergyLeakHis = new TH1D("EnergyLeakHis","",110,0,11); 
	EnergyLeakHis->GetXaxis()->SetTitle("Incident Electron Energy (GeV)");
	EnergyLeakHis->GetYaxis()->SetTitle("Average Back Energy Leakage (MeV)");
	EnergyLeakHis->GetXaxis()->CenterTitle();
	EnergyLeakHis->GetYaxis()->CenterTitle();
	EnergyLeakHis->GetYaxis()->SetTitleOffset(1.4);
		
	EnergyLeakHis->SetMarkerSize(1.5);
	EnergyLeakHis->SetMarkerStyle(20);
	EnergyLeakHis->SetMarkerColor(1);



	float MeanEnergy;
	float RMSEnergy;

	int NFiles = 16;

	for(int q = 0; q < NFiles; q++){

		double TotalBHEnergy = 0;
		double TotalFEMCALEnergy = 0;

	
		Energy = EStep * q + EInit;
		




		TString infile = Form("EnergyScanNew/ShowerInfo_%d.root",q);
		TFile * fin = new TFile(infile.Data());
		TTree * ShowerInfo = (TTree *) fin->Get("ShowerInfo");

		ShowerInfo->SetBranchAddress("evtid",&evtid);

		ShowerInfo->SetBranchAddress("hitid",&hitid);
		ShowerInfo->SetBranchAddress("trkid",&trkid);
		ShowerInfo->SetBranchAddress("showerid",&showerid);

		ShowerInfo->SetBranchAddress("X1",&X1);
		ShowerInfo->SetBranchAddress("X2",&X2);
		ShowerInfo->SetBranchAddress("Y1",&Y1);
		ShowerInfo->SetBranchAddress("Y2",&Y2);
		ShowerInfo->SetBranchAddress("Z1",&Z1);
		ShowerInfo->SetBranchAddress("Z2",&Z2);
		ShowerInfo->SetBranchAddress("T1",&T1);
		ShowerInfo->SetBranchAddress("BHBLedep",&BHBLedep);
		ShowerInfo->SetBranchAddress("BHFLedep",&BHFLedep);

		ShowerInfo->SetBranchAddress("FEMCALedep",&FEMCALedep);

		TH1D * BHedepHis = new TH1D("BHedepHis","",100,10,100);
		BHedepHis->GetXaxis()->SetTitle("Total Back Leakage Energy Per Event (MeV)");
		BHedepHis->GetYaxis()->SetTitle("Number of Events");
		BHedepHis->GetXaxis()->CenterTitle();
		BHedepHis->GetYaxis()->CenterTitle();
		BHedepHis->GetYaxis()->SetTitleOffset(1.4);





		TH1D * EnergyDepHis = new TH1D("EnergyDepHis","",100,10,5000); 
		EnergyDepHis->GetXaxis()->SetTitle("Total Energy Deposited to FEMCAL Per Event (MeV)");
		EnergyDepHis->GetYaxis()->SetTitle("Number of Events");
		EnergyDepHis->GetXaxis()->CenterTitle();
		EnergyDepHis->GetYaxis()->CenterTitle();
		EnergyDepHis->GetYaxis()->SetTitleOffset(1.4);




		TH1D * SamFracHis = new TH1D("SamFracHis","",100,0.1,0.5); 
		SamFracHis->GetXaxis()->SetTitle("Sampling Fraction");
		SamFracHis->GetYaxis()->SetTitle("Number of Events");
		SamFracHis->GetXaxis()->CenterTitle();
		SamFracHis->GetYaxis()->CenterTitle();
		SamFracHis->GetYaxis()->SetTitleOffset(1.4);

		int NEntries = ShowerInfo->GetEntries();
		int evtidPre = 0;



		for(int i = 0; i < NEntries; i++){

			ShowerInfo->GetEntry(i);


			//cout << "BHedep = " << BHBLedep << endl;

			if(evtid != evtidPre){


				BHedepHis->Fill(TotalBHEnergy);
				EnergyDepHis->Fill(TotalFEMCALEnergy);
				SamFracHis->Fill(TotalFEMCALEnergy/(Energy*1000));

				TotalBHEnergy = 0;
				TotalFEMCALEnergy = 0;


			}


			TotalBHEnergy = BHBLedep * 1000 + TotalBHEnergy;
			TotalFEMCALEnergy = FEMCALedep * 1000 + TotalFEMCALEnergy;


			evtidPre = evtid;

		}


		BHedepHis->Draw();



		lat->DrawLatex(0.30,0.80,Form("Incident Electron Energy = %.1f GeV",Energy));	

		c->SaveAs(Form("ELeakAnaLoop/BlackHole/BackLeakage_%d.png",q));

	


		EnergyDepHis->Draw();

		lat->DrawLatex(0.30,0.80,Form("Incident Electron Energy = %.1f GeV",Energy));	

		c->SaveAs(Form("ELeakAnaLoop/FEMCAL/EnergyDepHis_%d.png",q));

		SamFracHis->Draw();
		lat->DrawLatex(0.30,0.80,Form("Incident Electron Energy = %.1f GeV",Energy));	
		
		c->SaveAs(Form("ELeakAnaLoop/FEMCAL/SamFracHis_%d.png",q));
	
		MeanEnergy = BHedepHis->GetMean();
		RMSEnergy = BHedepHis->GetRMS();

		XBin = EnergyLeakHis->GetXaxis()->FindBin(Energy); 

		EnergyLeakHis->SetBinContent(XBin,MeanEnergy);
		EnergyLeakHis->SetBinError(XBin,RMSEnergy);

	}

	EnergyLeakHis->Draw("ep");
	c->SaveAs("ELeakAnaLoop/EnergyScan.png");

	

}
