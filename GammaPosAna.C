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


void GammaPosAna(int AngleBin){

	int DoCorr = 1;

	gStyle->SetOptStat(0);

	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();

	TLatex *lat = new TLatex();
	lat->SetNDC();
	lat->SetTextSize(0.035);

	const int NRealEvents = 5000;


	float PVx = 25;
	float PVy = 25;
	float PVz = 100;
	float EnterZ = 190.275;
	float RealZ;
	float Phi = 0;



	TString infile = Form("GammaPosAna/G4EICDetector.root_g4femc_eval_%d.root",AngleBin);
	TString infileTruth  = Form("GammaPosAna/ShowerInfo_%d.root",AngleBin);

	TFile * fin = new TFile(infile.Data());
	TFile * finTruth = new TFile(infileTruth.Data());

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
	float EnterX1[NRealEvents];
	float EnterY1[NRealEvents];
	float XPos[NRealEvents];
	float YPos[NRealEvents];
	float ZPos[NRealEvents];

	float XPosCorr[NRealEvents];
	float YPosCorr[NRealEvents];
	float ZPosCorr[NRealEvents];

	float XOverZ = 0;


	float etaMin = 1.25;
	float etaStep = 0.25;
	
	float eta = AngleBin * etaStep + etaMin;

	TVector3 * TruthVec = new TVector3;
	TVector3 * ClusVec = new TVector3;
	TVector3 * ClusCorrVec = new TVector3;



	float CosAngle;
	float CosAngleCorr;

	int NBinsX = 100;
	int ClusXMin = -5;
	int ClusXMax = 50;

	int NBinsXY = 100;
	int XMin = -5;
	int YMin = -5;
	int XMax = 5;
	int YMax= 5;
	
	int NBinsZ = 100;
	int ZMin = 12;
	int ZMax = 18;


	int NAngle = 104;
	float MinAngle = -1.02;
	float MaxAngle = 1.02;
	

	int NBinsXOverZ = 100;
	float XOverZMin = 0.0;
	float XOverZMax = 3.0;


	TH2D * ClusXY = new TH2D("ClusXY","",NBinsXY,XMin,XMax,NBinsXY,YMin,YMax);
	ClusXY->GetXaxis()->SetTitle("Reco Cluster: X (cm)");
	ClusXY->GetYaxis()->SetTitle("Reco Cluster: Y (cm)");
	ClusXY->GetXaxis()->CenterTitle();
	ClusXY->GetYaxis()->CenterTitle();
	ClusXY->GetYaxis()->SetTitleOffset(1.4);
	ClusXY->SetTitle("Cluster Reconstructed XY Distribution");

	TH2D * TruthXY = new TH2D("TruthXY","",NBinsXY,XMin,XMax,NBinsXY,YMin,YMax);
	TruthXY->GetXaxis()->SetTitle("EMCAL Front Face Truth X (cm)");
	TruthXY->GetYaxis()->SetTitle("EMCAL Front Face Truth Y (cm)");
	TruthXY->SetTitle("EMCAL Front Face Photon XY Distribution");	
	TruthXY->GetXaxis()->CenterTitle();
	TruthXY->GetYaxis()->CenterTitle();
	TruthXY->GetYaxis()->SetTitleOffset(1.4);


	TH1D * ClusZ = new TH1D("ClusXZ","",NBinsZ,ZMin,ZMax);
	ClusZ->GetXaxis()->SetTitle("Reco Cluster: Z (cm) - Distance from the Front Face");
	ClusZ->GetYaxis()->SetTitle("Number of Events");
	ClusZ->GetXaxis()->CenterTitle();
	ClusZ->GetYaxis()->CenterTitle();
	ClusZ->GetYaxis()->SetTitleOffset(1.4);
	ClusZ->SetTitle("Cluster Reconstructed Z Distribution");

	TH1D * CosAngleHis = new TH1D("CosAngleHis","",NAngle,MinAngle,MaxAngle);
	CosAngleHis->GetXaxis()->SetTitle("cos(#alpha)");
	CosAngleHis->GetYaxis()->SetTitle("Number of Events");
	CosAngleHis->GetXaxis()->CenterTitle();
	CosAngleHis->GetYaxis()->CenterTitle();
	CosAngleHis->GetYaxis()->SetTitleOffset(1.4);
	CosAngleHis->SetTitle("cos(#alpha) Distribution");

	TH1D * CosAngleHisCorr = new TH1D("CosAngleHisCorr","",NAngle,MinAngle,MaxAngle);
	CosAngleHisCorr->GetXaxis()->SetTitle("cos(#alpha) Corrected");
	CosAngleHisCorr->GetYaxis()->SetTitle("Number of Events");
	CosAngleHisCorr->GetXaxis()->CenterTitle();
	CosAngleHisCorr->GetYaxis()->CenterTitle();
	CosAngleHisCorr->GetYaxis()->SetTitleOffset(1.4);
	CosAngleHisCorr->SetTitle("After Correction");
	

	TH1D * ClusX = new TH1D("ClusX","",NBinsX,ClusXMin,ClusXMax);
	ClusX->GetXaxis()->SetTitle("Reco Cluster: X (cm) - Distance from the Front Face");
	ClusX->GetYaxis()->SetTitle("Number of Events");
	ClusX->GetXaxis()->CenterTitle();
	ClusX->GetYaxis()->CenterTitle();
	ClusX->GetYaxis()->SetTitleOffset(1.4);
	ClusX->SetTitle("Cluster Reconstructed X Distribution");



	TH1D * ClusXOverZ = new TH1D("ClusXOverZ","",NBinsXOverZ,XOverZMin,XOverZMax);
	ClusXOverZ->GetXaxis()->SetTitle("Reco Cluster: X (cm) - Distance from the Front Face");
	ClusXOverZ->GetYaxis()->SetTitle("Number of Events");
	ClusXOverZ->GetXaxis()->CenterTitle();
	ClusXOverZ->GetYaxis()->CenterTitle();
	ClusXOverZ->GetYaxis()->SetTitleOffset(1.4);
	ClusXOverZ->SetTitle("Cluster Reconstructed X Distribution");




	for(int i = 0 ; i < NRealEvents; i ++){

		EnterX1[i] = 0;
		EnterY1[i] = 0;
		
		XPosCorr[i] = 0;
		YPosCorr[i] = 0;
		ZPosCorr[i] = 0;	


		XPos[i] = 0;
		YPos[i] = 0;
		ZPos[i] = 0;	
		
	}
	

	int NClus = 0;
	
	int Outlier = 0;

	float GPX1;
	float GPX2;
	float GPY1;
	float GPY2;
	float GPZ1;
	float GPZ2;

	TTree * ShowerInfo = (TTree *) finTruth->Get("ShowerInfo");

	ShowerInfo->SetBranchAddress("evtid",&evtid);
	ShowerInfo->SetBranchAddress("GPX1",&GPX1);
	ShowerInfo->SetBranchAddress("GPY1",&GPY1);
	ShowerInfo->SetBranchAddress("GPZ1",&GPZ1);
	ShowerInfo->SetBranchAddress("GPX2",&GPX2);
	ShowerInfo->SetBranchAddress("GPY2",&GPY2);
	ShowerInfo->SetBranchAddress("GPZ2",&GPZ2);


	TTree * ntp_cluster = (TTree * ) fin->Get("ntp_cluster");

	ntp_cluster->SetBranchAddress("event",&EventFloat);
	ntp_cluster->SetBranchAddress("clusterID",&clusterIDFloat);
	ntp_cluster->SetBranchAddress("e",&RECOEnergy);
	ntp_cluster->SetBranchAddress("x",&x);
	ntp_cluster->SetBranchAddress("y",&y);
	ntp_cluster->SetBranchAddress("z",&z);

	
	int NEventsTruth = ShowerInfo->GetEntries();
	int NEvents  = ntp_cluster->GetEntries();
	
	//cout << "NEvents = " << NEvents << endl;

	for(int i = 0; i < NEventsTruth; i++){
	
		ShowerInfo->GetEntry(i);
		

		if(evtidPre != evtid){

			EnterX1[evtidPre] = GPX1/GPZ1 * EnterZ;
			EnterY1[evtidPre] = GPY1/GPZ1 * EnterZ;
	
			
		}

		evtidPre = evtid;

	}

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
					
				XPosCorr[EventPre] = x - PVx - ((z - PVz) * 2/(TMath::Exp(eta) - TMath::Exp(-eta)));
				YPosCorr[EventPre] = YPos[EventPre];	
				ZPosCorr[EventPre] = ZPos[EventPre];

			

				if(XPosCorr[EventPre] > 5 || YPosCorr[EventPre] > 5) Outlier  = Outlier + 1;

			}
		
			NClus = 0;
			
		}
	
		NClus = NClus + 1;
		
		EventPre = Event;
	}



	for(int i = 0; i < NRealEvents; i++){
	
//		cout  << "i = " << i << "   EnterX1[i] =  "  << EnterX1[i] << "   EnterY1[i] =  "  <<  EnterY1[i]   <<    "   XPos[i] =  "    <<   XPos[i]  <<    "   YPos[i] =  "  << YPos[i] << endl;

		if(XPos[i] * YPos[i] != 0 ){
			
			TruthXY->Fill(EnterX1[i],EnterY1[i]);	
			ClusXY->Fill(XPosCorr[i],YPosCorr[i]);
			ClusZ->Fill(ZPos[i]);
			ClusX->Fill(XPos[i]);

			TruthVec->SetXYZ(EnterX1[i],EnterY1[i],EnterZ);
			ClusVec->SetXYZ(XPos[i],YPos[i],ZPos[i]);

			
			CosAngle = ClusVec->Unit().Dot(TruthVec->Unit());

			

			CosAngleHis->Fill(CosAngle);


			XOverZ = XPos[i]/ZPos[i];
			
			ClusXOverZ->Fill(XOverZ);

		}

	}
	

	cout << "XY Center: " << "  X = " << EnterX1[1] << "   Y =  " << EnterY1[1] << endl;

	//TruthXY->Draw("COLZ");
	//c->SaveAs(Form("GammaPosAnaPlots/%d/TruthXY.png",AngleBin));

	ClusXY->Draw("COLZ");
	lat->DrawLatex(0.20,0.85,"Incident Single #gamma Eneregy = 1 GeV");		
	lat->DrawLatex(0.20,0.80,Form("Incident Single #gamma #eta = %.1f (#phi = 0)",eta));	
	lat->DrawLatex(0.20,0.75,"Standard PHENIX Pb Shashlik EMCAL");	
	
	c->SaveAs(Form("GammaPosAnaPlots/%d/ClusXY.png",AngleBin));



	ClusXOverZ->Draw();
	lat->DrawLatex(0.20,0.85,"Incident Single #gamma Eneregy = 1 GeV");		
	lat->DrawLatex(0.20,0.80,Form("Incident Single #gamma #eta = %.1f (#phi = 0)",eta));	
	lat->DrawLatex(0.20,0.75,"Standard PHENIX Pb Shashlik EMCAL");	

	c->SaveAs(Form("GammaPosAnaPlots/%d/ClusXOverZ.png",AngleBin));



	ClusZ->Draw();
	lat->DrawLatex(0.20,0.85,"Incident Single #gamma Eneregy = 1 GeV");		
	lat->DrawLatex(0.20,0.80,Form("Incident Single #gamma #eta = %.1f (#phi = 0)",eta));	
	lat->DrawLatex(0.20,0.75,"Standard PHENIX Pb Shashlik EMCAL");	

	c->SaveAs(Form("GammaPosAnaPlots/%d/ClusZ.png",AngleBin));

	
	ClusX->Draw();
	lat->DrawLatex(0.20,0.85,"Incident Single #gamma Eneregy = 1 GeV");		
	lat->DrawLatex(0.20,0.80,Form("Incident Single #gamma #eta = %.1f (#phi = 0)",eta));	
	lat->DrawLatex(0.20,0.75,"Standard PHENIX Pb Shashlik EMCAL");	

	c->SaveAs(Form("GammaPosAnaPlots/%d/ClusX.png",AngleBin));


	c->SetLogy();

	CosAngleHis->Draw();
	lat->DrawLatex(0.20,0.85,"Incident Single #gamma Eneregy = 1 GeV");		
	lat->DrawLatex(0.20,0.80,Form("Incident Single #gamma #eta = %.1f (#phi = 0)",eta));	
	lat->DrawLatex(0.20,0.75,"Standard PHENIX Pb Shashlik EMCAL");	
	c->SaveAs(Form("GammaPosAnaPlots/%d/CosAngleHis.png",AngleBin));





	ClusXOverZ->Draw();
	lat->DrawLatex(0.20,0.85,"Incident Single #gamma Eneregy = 1 GeV");		
	lat->DrawLatex(0.20,0.80,Form("Incident Single #gamma #eta = %.1f (#phi = 0)",eta));	
	lat->DrawLatex(0.20,0.75,"Standard PHENIX Pb Shashlik EMCAL");	
	c->SaveAs(Form("GammaPosAnaPlots/%d/ClusXOverZLog.png",AngleBin));



	cout << "Outlier = " << Outlier << endl;




}
