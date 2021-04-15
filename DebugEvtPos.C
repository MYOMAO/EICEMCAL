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



void DebugEvtPos(int Index){

	gStyle->SetOptStat(0);
	

	float RadLength = 0.56;

	float PVx = 75;
	float PVy = 75;
	float PVz = 100;
	float EnterZ = 290.275;
	
	float PeakZ = 15.90;

	float x;
	float y;
	float z;
	float RECOEnergy;


	float xReal;
	float yReal;
	float zReal;


	float eventFloat;
	float clusterIDFloat;
	int event;
	int cluster;

	int NBinsX = 200;
	int NBinsY = 200;

	float XMin = -20;
	float XMax = 20;
	float YMin = -20;
	float YMax = 20;

	int NBinsZ = 100;
	float ZMin = 14.0;
	float ZMax = 17.0;

	int NBinsE = 100;
	float EMin = 0.5;
	float EMax = 1.1;


	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();

	TString infile = "PosEvtDebug/G4EICDetector.root_g4femc_eval.root";

	TFile * fin = new TFile(infile.Data());
	fin->cd();

	
	TTree * ntp_cluster = (TTree * ) fin->Get("ntp_cluster");

	ntp_cluster->SetBranchAddress("event",&eventFloat);
	ntp_cluster->SetBranchAddress("clusterID",&clusterIDFloat);
	ntp_cluster->SetBranchAddress("e",&RECOEnergy);
	ntp_cluster->SetBranchAddress("x",&x);
	ntp_cluster->SetBranchAddress("y",&y);
	ntp_cluster->SetBranchAddress("z",&z);



	TH2D * SinglePhotonXY = new TH2D("SinglePhotonXY","",NBinsX,XMin,XMax,NBinsY,YMin,YMax);
	SinglePhotonXY->GetXaxis()->SetTitle("Cluster Reconstructed X (cm)");
	SinglePhotonXY->GetYaxis()->SetTitle("Cluster Reconstructed Y (cm)");
	SinglePhotonXY->SetTitle("EMCAL Cluster Reconstructed XY Distribution with Respect to the #pi^{0} Decay Vertex");	
	SinglePhotonXY->GetXaxis()->CenterTitle();
	SinglePhotonXY->GetYaxis()->CenterTitle();
	SinglePhotonXY->GetYaxis()->SetTitleOffset(1.4);




	TH1D * SinglePhotonZ = new TH1D("SinglePhotonZ","",NBinsZ,ZMin,ZMax);
	SinglePhotonZ->GetXaxis()->SetTitle("Cluster Reconstructed Z (cm)");
	SinglePhotonZ->GetYaxis()->SetTitle("Number of Events");
	SinglePhotonZ->SetTitle("EMCAL Cluster Reconstructed Z Distribution with Respect to the EMCAL Front Face");	
	SinglePhotonZ->GetXaxis()->CenterTitle();
	SinglePhotonZ->GetYaxis()->CenterTitle();
	SinglePhotonZ->GetYaxis()->SetTitleOffset(1.4);



	TH1D * ClusE = new TH1D("ClusE","",NBinsE,EMin,EMax);
	ClusE->GetXaxis()->SetTitle("Cluster Reconstructed Energy (GeV)");
	ClusE->GetYaxis()->SetTitle("Number of Events");
	ClusE->SetTitle("EMCAL Cluster Reconstructed Energy Distribution");	
	ClusE->GetXaxis()->CenterTitle();
	ClusE->GetYaxis()->CenterTitle();
	ClusE->GetYaxis()->SetTitleOffset(1.4);



	int NEntries = ntp_cluster->GetEntries();

	std::ofstream DebugBro2("DatFilesPos/DebugEvt2.dat",ios::app);
	
	for(int i = 0; i < NEntries; i++){
		
		ntp_cluster->GetEntry(i);
		
		xReal = x - PVx;
		yReal = y - PVy;
		zReal = z - EnterZ;

		SinglePhotonXY->Fill(xReal,yReal);
		SinglePhotonZ->Fill(zReal);
		ClusE->Fill(RECOEnergy);
		DebugBro2  << "Index = " << Index << "   xReal = " << xReal << "   yReal = " << yReal  << endl;

	}
	


	SinglePhotonXY->Draw("COLZ");

	c->SaveAs("PosDebug/SinglePhotonXY.png");

	SinglePhotonZ->Draw("COLZ");
	
	TLine *lPeak = new TLine(PeakZ,0,PeakZ,SinglePhotonZ->GetMaximum());
	lPeak->SetLineStyle(2);
	lPeak->SetLineWidth(2);
	lPeak->SetLineColor(2);
	lPeak->Draw("SAME");

	c->SaveAs("PosDebug/SinglePhotonZ.png");

	ClusE->Draw();
	c->SaveAs("PosDebug/ClusE.png");


}
