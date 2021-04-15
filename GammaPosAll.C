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
#include "PlotEDep.h"
using namespace std;

using std::cout;
using std::endl;


void GammaPosAll(){

	gStyle->SetOptStat(0);

	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();

	
	const int NRealEvents = 5000;
	const int NFiles = 27;


	float PVx = 25;
	float PVy = 25;
	float PVz = 100;
	float EnterZ = 190.275;
	float RealZ;
	float Phi = 0;

	float etaMin = 1.5;
	float etaStep = 0.25;
	float etaMax = etaMin + etaStep * NFiles;
	float eta;

	int NEtaBin = NFiles + 1;
	float EtaHisMin = etaMin  - etaStep/2;
	float EtaHisMax = etaMax  + etaStep/2;


	TF1 * Xfunc = new TF1("XFunc","(16 + 190.275) * 2/(TMath::Exp(x) - TMath::Exp(-x))",etaMin,etaMax);
	Xfunc->SetLineColor(kRed);

	TH1D * ClusZvsEta = new TH1D("ClusZvsEta","",NEtaBin,EtaHisMin,EtaHisMax);
	ClusZvsEta->GetXaxis()->SetTitle("#gamma Beam #eta");
	ClusZvsEta->GetYaxis()->SetTitle("Clustered Z Value (cm)");
	ClusZvsEta->GetXaxis()->CenterTitle();
	ClusZvsEta->GetYaxis()->CenterTitle();
	ClusZvsEta->GetYaxis()->SetTitleOffset(1.4);
	ClusZvsEta->SetTitle("Cluster Reconstructed Z vs #gamma Beam #eta");

	ClusZvsEta->SetMarkerColor(kBlack);
	ClusZvsEta->SetLineColor(kBlack);	
	ClusZvsEta->SetMarkerStyle(20);
	ClusZvsEta->SetMarkerSize(1);

	TH1D * ClusXvsEta = new TH1D("ClusXvsEta","",NEtaBin,EtaHisMin,EtaHisMax);
	ClusXvsEta->GetXaxis()->SetTitle("#gamma Beam #eta");
	ClusXvsEta->GetYaxis()->SetTitle("Clustered X Value (cm)");
	ClusXvsEta->GetXaxis()->CenterTitle();
	ClusXvsEta->GetYaxis()->CenterTitle();
	ClusXvsEta->GetYaxis()->SetTitleOffset(1.4);
	ClusXvsEta->SetTitle("Cluster Reconstructed X vs #gamma Beam #eta");


	ClusXvsEta->SetMarkerColor(kBlack);
	ClusXvsEta->SetLineColor(kBlack);	
	ClusXvsEta->SetMarkerStyle(20);
	ClusXvsEta->SetMarkerSize(1);



	for(int q = 1; q < NFiles; q++){

		eta = etaMin + etaStep * q;

		float XPos[NRealEvents];
		float YPos[NRealEvents];
		float ZPos[NRealEvents];


		for(int i = 0 ; i < NRealEvents; i ++){




			XPos[i] = 0;
			YPos[i] = 0;
			ZPos[i] = 0;	

		}

		float clusterIDFloat;
		float EventFloat;
		float RECOEnergy;
		float x;
		float y;
		float z;

		int clusterID;
		int Event = 0;
		int EventPre = 0;
		int evtid = 0;
		int evtidPre = 0;


		float MeanX = 0;
		float SumX = 0;
		float SumSquareX = 0;
		float RMSX = 0;

		float MeanZ = 0;
		float SumZ = 0;
		float SumSquareZ = 0;
		float RMSZ = 0;
	
		int NClus = 0;


		int N = 0;

		TString infile = Form("GammaPosAna/G4EICDetector.root_g4femc_eval_%d.root",q);
		TFile * fin = new TFile(infile.Data());

		TTree * ntp_cluster = (TTree * ) fin->Get("ntp_cluster");
		ntp_cluster->SetBranchAddress("event",&EventFloat);
		ntp_cluster->SetBranchAddress("clusterID",&clusterIDFloat);
		ntp_cluster->SetBranchAddress("x",&x);
		ntp_cluster->SetBranchAddress("z",&z);

		int NEvents = ntp_cluster->GetEntries();

		for(int i = 0; i < NEvents; i++){

			ntp_cluster->GetEntry(i);

			Event = int(EventFloat);
			clusterID = int(clusterIDFloat);	

			//cout << "EventPre = " << EventPre << "    Event = " << Event << endl;

			if(EventPre != Event){

				//cout << "NClus = " << NClus << endl;
				if(NClus == 1){

					XPos[EventPre] = x - PVx;
					YPos[EventPre] = y - PVy;
					ZPos[EventPre] = z - EnterZ - PVz;

					SumX = SumX + XPos[EventPre];
					SumZ = SumZ + ZPos[EventPre];

					N = N + 1;
				}

				NClus = 0;

			}

			NClus = NClus + 1;

			EventPre = Event;
		}

		MeanX = SumX/N;
		MeanZ = SumZ/N;



		for(int i = 0; i < NEvents; i++){

			ntp_cluster->GetEntry(i);

			Event = int(EventFloat);
			clusterID = int(clusterIDFloat);	

			//cout << "EventPre = " << EventPre << "    Event = " << Event << endl;

			if(EventPre != Event){

				//cout << "NClus = " << NClus << endl;
				if(NClus == 1){

			
					SumSquareX = SumSquareX + (XPos[EventPre] - MeanX) * (XPos[EventPre] - MeanX);
					SumSquareZ = SumSquareZ + (ZPos[EventPre] - MeanZ) * (ZPos[EventPre] - MeanZ);

					N = N + 1;
				}

				NClus = 0;

			}

			NClus = NClus + 1;

			EventPre = Event;
		}


		RMSX = sqrt(SumSquareX)/N;
		RMSZ = sqrt(SumSquareZ)/N;


		ClusZvsEta->SetBinContent(q,MeanZ);
		ClusZvsEta->SetBinError(q,RMSZ);


		ClusXvsEta->SetBinContent(q,MeanX);
		ClusXvsEta->SetBinError(q,RMSX);


	}

	ClusZvsEta->Draw("ep");
	c->SaveAs("GammaPosAnaPlots/All/ClusZvsEta.png");

	ClusXvsEta->Draw("ep");
	Xfunc->Draw("SAME");

	TLegend* leg = new TLegend(0.50,0.60,0.80,0.80,NULL,"brNDC");
	leg->SetBorderSize(0);
	leg->SetTextSize(0.032);
	leg->SetTextFont(42);
	leg->SetFillStyle(0);

	leg->AddEntry(ClusXvsEta,"Simulation","LP");
	leg->AddEntry(Xfunc,"Expectation","L");
	
	leg->Draw("SAME");

	c->SaveAs("GammaPosAnaPlots/All/ClusXvsEta.png");


}
