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


void QuickAna(){

	gStyle->SetOptStat(0);

	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();

	const int NFiles = 29;
	int FileMin = 0;
	float EStep = 0.001;
	float EMin = FileMin * EStep;	
	float EMax = NFiles * EStep;

	const int NDecayMax = 5; 

	int Event = 0;
	int EventPre = 0;

	float EventFloat;
	float clusterIDFloat;

	int clusterID;

	float RECOEnergy;
	float E[NDecayMax];


	for(int i = 0; i < NDecayMax; i++){

		E[i] = 0;

	}

	TString infile;

	TH1D * EVSEThres = new  TH1D("EVSEThres","",NFiles+1,EMin-EStep/2,EMax+EStep/2);

	EVSEThres->GetXaxis()->SetTitle("Minimum Tower Energy Threshold (GeV)");
	EVSEThres->GetYaxis()->SetTitle("Reconstructed Cluster Energy (GeV)");
	EVSEThres->GetXaxis()->CenterTitle();
	EVSEThres->GetYaxis()->CenterTitle();
	EVSEThres->GetYaxis()->SetTitleOffset(1.4);
	EVSEThres->SetMarkerStyle(20);
	EVSEThres->SetMarkerSize(1);
	EVSEThres->SetMarkerColor(1);

	EVSEThres->SetMinimum(0.84);
	EVSEThres->SetMaximum(0.89);


	for(int i = FileMin; i < NFiles; i++){

		infile = Form("EThresholdScan/G4EICDetector.root_g4femc_eval_%d.root",i);

		TFile * fin = new TFile(infile.Data());
		fin->cd();

		TH1D * RECOEHis = new TH1D("RECOEHis","",100,0,1);
		RECOEHis->GetXaxis()->SetTitle("Reconstructed Cluster Energy (GeV)");
		RECOEHis->GetYaxis()->SetTitle("Counts");
		RECOEHis->GetXaxis()->CenterTitle();
		RECOEHis->GetYaxis()->CenterTitle();
		RECOEHis->GetYaxis()->SetTitleOffset(1.4);


		RECOEHis->SetMarkerStyle(20);
		RECOEHis->SetMarkerSize(1);
		RECOEHis->SetMarkerColor(1);


		TTree * ntp_cluster = (TTree * ) fin->Get("ntp_cluster");
		ntp_cluster->SetBranchAddress("event",&EventFloat);		
		ntp_cluster->SetBranchAddress("e",&RECOEnergy);
		ntp_cluster->SetBranchAddress("clusterID",&clusterIDFloat);		


		int NEvents = ntp_cluster->GetEntries();

		for(int j = 0; j < NEvents; j++){

			ntp_cluster->GetEntry(j);

			Event = int(EventFloat);
			clusterID = int(clusterIDFloat);

			if(EventPre != Event){

				if(E[0] > 0 && E[1] > 0 && E[2] == 0 && E[3] == 0 && E[4] == 0){

					RECOEHis->Fill(E[0] + E[1]);

				}

				for(int i = 0; i < NDecayMax; i++){

					E[i] = 0;

				}

			}

			E[clusterID] = RECOEnergy;


			EventPre = Event;

		}


		TF1 * func = new TF1("func","gaus",0.75,1);

		RECOEHis->Fit(func,"R");
		RECOEHis->Draw("ep");


		c->SaveAs(Form("EThresholdScan/Fits/RECOEHisFits_%d.png",i));


		float Mean = func->GetParameter(1);
		float MeanErr = func->GetParError(1);

		float RMS = func->GetParameter(2);

		EVSEThres->SetBinContent(i+1,Mean);
		//EVSEThres->SetBinError(i+1,RMS);
		EVSEThres->SetBinError(i+1,MeanErr);

	}

	EVSEThres->Draw("ep");

	c->SaveAs("EThresholdScan/EThresScan.png");


}

