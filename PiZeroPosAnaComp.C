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
#include "TLorentzVector.h"

using namespace std;

using std::cout;
using std::endl;


void PiZeroPosAnaComp(int doCut){
	
	TString Name;

	if(doCut == 0) Name = "NoCut";
	if(doCut == 1) Name = "DoCut";


	TString Comments;



	if(doCut == 0) Comments = "No Selection Applied";
	if(doCut == 1) Comments = "Tower Near Boundary Required to Have 0 Energy";


	gStyle->SetOptStat(0);




	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();
	c->SetLogy();

	TLatex *lat = new TLatex();
	lat->SetNDC();
	lat->SetTextSize(0.035);

	float PVx = 75;
	float PVy = 75;
	float PVz = 100;

	float IncEnergy = 1.0;

	float EnterZ = 190.27500;

	//TFile * fin = new TFile("PosAna/G4EICDetector.root_g4femc_eval.root");

	TFile * fin = new TFile(Form("PosAnaComp/G4EICDetector.root_g4femc_eval_%d.root",doCut));


	float x;
	float y;
	float z;
	int event;
	int clusterID;
	float eventFloat;
	float clusterIDFloat;
	float RECOEnergy;

	const int NEvents = 1000;



	float ClusX1[NEvents];
	float ClusX2[NEvents];

	float ClusY1[NEvents];
	float ClusY2[NEvents];

	float ClusZ1[NEvents];
	float ClusZ2[NEvents];

	float ClusR1[NEvents];
	float ClusR2[NEvents];


	float ClusE1[NEvents];
	float ClusE2[NEvents];
	

	float ClusPX1[NEvents];
	float ClusPX2[NEvents];

	float ClusPY1[NEvents];
	float ClusPY2[NEvents];

	float ClusPZ1[NEvents];
	float ClusPZ2[NEvents];



	TTree * ntp_cluster = (TTree * ) fin->Get("ntp_cluster");

	ntp_cluster->SetBranchAddress("event",&eventFloat);
	ntp_cluster->SetBranchAddress("clusterID",&clusterIDFloat);
	ntp_cluster->SetBranchAddress("e",&RECOEnergy);
	ntp_cluster->SetBranchAddress("x",&x);
	ntp_cluster->SetBranchAddress("y",&y);
	ntp_cluster->SetBranchAddress("z",&z);



	int NEntries = ntp_cluster->GetEntries(); 

	int NBinsXY = 100;
	int XMin = -200;
	int YMin = -200;
	int XMax = 200;
	int YMax= 200;

	int NBinsZ = 150;
	int ZMin = 185;
	int ZMax = 320;

	int NBinsR = 100;
	int RMin = 0;
	int RMax = 250;
	
	int NBinsE = 110;
	float EMin = 0;
	float EMax = 1.10;

	int NBinsPiMass = 100;
	float PiMassMin = 0;
	float PiMassMax = 1.0;

	float EAsym;


	TH2D * TruthXY = new TH2D("TruthXY","",NBinsXY,XMin,XMax,NBinsXY,YMin,YMax);
	TruthXY->GetXaxis()->SetTitle("EMCAL Front Face Truth X (cm)");
	TruthXY->GetYaxis()->SetTitle("EMCAL Front Face Truth Y (cm)");
	TruthXY->SetTitle("EMCAL Front Face Photon XY Distribution");	
	TruthXY->GetXaxis()->CenterTitle();
	TruthXY->GetYaxis()->CenterTitle();
	TruthXY->GetYaxis()->SetTitleOffset(1.4);


	TH2D * ClusXY = new TH2D("ClusXY","",NBinsXY,XMin,XMax,NBinsXY,YMin,YMax);
	ClusXY->GetXaxis()->SetTitle("Cluster Reconstructed X (cm)");
	ClusXY->GetYaxis()->SetTitle("Cluster Reconstructed Y (cm)");
	ClusXY->SetTitle("EMCAL Cluster Reconstructed XY Distribution");	
	ClusXY->GetXaxis()->CenterTitle();
	ClusXY->GetYaxis()->CenterTitle();
	ClusXY->GetYaxis()->SetTitleOffset(1.4);



	TH1D * ClusZ = new TH1D("ClusZ","",NBinsZ,ZMin,ZMax);
	ClusZ->GetXaxis()->SetTitle("Cluster Reconstructed Z (cm)");
	ClusZ->GetYaxis()->SetTitle("Number of Events");
	ClusZ->SetTitle("EMCAL Cluster Reconstructed Z Distribution");	
	ClusZ->GetXaxis()->CenterTitle();
	ClusZ->GetYaxis()->CenterTitle();
	ClusZ->GetYaxis()->SetTitleOffset(1.4);


	TH1D * TruthR = new TH1D("TruthR","",NBinsR,RMin,RMax);
	TruthR->GetXaxis()->SetTitle("Distance Between Two Photons on EMCAL Front Face");
	TruthR->GetYaxis()->SetTitle("Number of Events");
	TruthR->GetXaxis()->CenterTitle();
	TruthR->GetYaxis()->CenterTitle();
	TruthR->GetYaxis()->SetTitleOffset(1.4);



	TH1D * ClusR = new TH1D("ClusR","",NBinsR,RMin,RMax);
	ClusR->GetXaxis()->SetTitle("Distance Between Two Clusters Centers");
	ClusR->GetYaxis()->SetTitle("Number of Events");
	ClusR->GetXaxis()->CenterTitle();
	ClusR->GetYaxis()->CenterTitle();
	ClusR->GetYaxis()->SetTitleOffset(1.4);




	TH1D * ClusE = new TH1D("ClusE","",NBinsE,EMin,EMax);
	ClusE->GetXaxis()->SetTitle("Total Energy of the Two Photon Clusters");
	ClusE->GetYaxis()->SetTitle("Number of Events");
	ClusE->GetXaxis()->CenterTitle();
	ClusE->GetYaxis()->CenterTitle();
	ClusE->GetYaxis()->SetTitleOffset(1.4);


	TH1D * ClusEAsym = new TH1D("ClusEAsym","",200,-1,1);
	ClusEAsym->GetXaxis()->SetTitle("(E_{1} - E_{2})/(E_{1} + E_{2})");
	ClusEAsym->GetYaxis()->SetTitle("Number of Events");
	ClusEAsym->GetXaxis()->CenterTitle();
	ClusEAsym->GetYaxis()->CenterTitle();
	ClusEAsym->GetYaxis()->SetTitleOffset(1.4);



	TH1D * PiMassRECO = new TH1D("PiMassRECO","",NBinsPiMass,PiMassMin,PiMassMax);
	PiMassRECO->GetXaxis()->SetTitle("#pi Invariant Mass: #pi #rightarrow #gamma #gamma: (GeV/c^{3})");
	PiMassRECO->GetYaxis()->SetTitle("Number of Events");
	PiMassRECO->GetXaxis()->CenterTitle();
	PiMassRECO->GetYaxis()->CenterTitle();
	PiMassRECO->GetYaxis()->SetTitleOffset(1.4);


	TFile * finGen = new TFile(Form("PosAnaComp/ShowerInfo_%d.root",doCut));

	float X1;
	float Y1;
	float Z1;

	float X2;
	float Y2;
	float Z2;


	float GPX1;
	float GPY1;
	float GPZ1;

	float GPX2;
	float GPY2;
	float GPZ2;
	int evtid;
	int layer;


	float PiMass;

	TTree * ShowerInfo = (TTree * ) finGen->Get("ShowerInfo");

	ShowerInfo->SetBranchAddress("evtid",&evtid);
	ShowerInfo->SetBranchAddress("layer",&layer);

	ShowerInfo->SetBranchAddress("GPX1",&GPX1);
	ShowerInfo->SetBranchAddress("GPY1",&GPY1);
	ShowerInfo->SetBranchAddress("GPZ1",&GPZ1);
	ShowerInfo->SetBranchAddress("GPX2",&GPX2);
	ShowerInfo->SetBranchAddress("GPY2",&GPY2);
	ShowerInfo->SetBranchAddress("GPZ2",&GPZ2);


	float EnterX1[NEvents];
	float EnterX2[NEvents];

	float EnterY1[NEvents];
	float EnterY2[NEvents];

	float EnterR;

	float ClusterR;


	float TotalClusE;
	 
	TLorentzVector * Gamma1 = new TLorentzVector;
	TLorentzVector * Gamma2 = new TLorentzVector;

	TLorentzVector * PiZero = new TLorentzVector; 


	for(int i = 0; i < NEvents; i++){

		EnterX1[i] = 0;
		EnterX2[i] = 0;

		EnterY1[i] = 0;
		EnterY2[i] = 0;
		ClusX1[i] = 0;
		ClusX2[i] = 0;

		ClusY1[i] = 0;
		ClusY2[i] = 0;

		ClusZ1[i] = 0;
		ClusZ2[i] = 0;

		ClusR1[i] = 0;
		ClusR2[i] = 0;

	
		ClusPX1[i] = 0;
		ClusPX2[i] = 0;

		ClusPY1[i] = 0;
		ClusPY2[i] = 0;


		ClusPZ1[i] = 0;
		ClusPZ2[i] = 0;


	}

	

	int NEntriesGen = ShowerInfo->GetEntries(); 

	for(int i = 0; i < NEntries; i++){

		ntp_cluster->GetEntry(i);



		clusterID = int(clusterIDFloat);
		event = int(eventFloat);
		//cout << "event = " << event << endl;

		if(clusterID == 0){
			ClusX1[event] = x - PVx;
			ClusY1[event] = y - PVy;
			ClusZ1[event] = z - PVz;
			
			ClusR1[event] = sqrt(ClusX1[event] * ClusX1[event] + ClusY1[event] * ClusY1[event] + ClusZ1[event] * ClusZ1[event] );
			
			ClusE1[event] = RECOEnergy;
		

			ClusPX1[i] = ClusE1[event]  * ClusX1[event]/ClusR1[event];
			ClusPY1[i] = ClusE1[event]  * ClusY1[event]/ClusR1[event];
			ClusPZ1[i] = ClusE1[event]  * ClusZ1[event]/ClusR1[event];

		}

		if(clusterID == 1){
			ClusX2[event] = x - PVx;
			ClusY2[event] = y - PVy;
			ClusZ2[event] = z - PVz;
			ClusR2[event] = sqrt(ClusX2[event] * ClusX2[event] + ClusY2[event] * ClusY2[event] + ClusZ2[event] * ClusZ2[event] );
			
			ClusE2[event] = RECOEnergy;

			ClusPX2[i] = ClusE2[event]  * ClusX2[event]/ClusR2[event];
			ClusPY2[i] = ClusE2[event]  * ClusY2[event]/ClusR2[event];
			ClusPZ2[i] = ClusE2[event]  * ClusZ2[event]/ClusR2[event];


		}

	}



	for(int i = 0; i < NEntriesGen; i++){

		ShowerInfo->GetEntry(i);
		
		//cout << "evtid = " << evtid << endl;

		EnterX1[evtid] = GPX1/GPZ1 * EnterZ;
		EnterX2[evtid] = GPX2/GPZ2 * EnterZ;
		EnterY1[evtid] = GPY1/GPZ1 * EnterZ;
		EnterY2[evtid] = GPY2/GPZ2 * EnterZ;


	}

	TVector3 * PosClus1 = new TVector3;
	TVector3 * PosTruth1  = new TVector3;
	TVector3 * PosClus2  = new TVector3;
	TVector3 * PosTruth2  = new TVector3;
	TVector3 * Dis1  = new TVector3;
	TVector3 * Dis2  = new TVector3;

	double CosAngle;

	TH1D * AngleDis = new TH1D("AngleDis","",100,-1.1,1.1);
	AngleDis->GetXaxis()->SetTitle("cos(#theta)");
	AngleDis->GetYaxis()->SetTitle("Number of Events");
	AngleDis->GetXaxis()->CenterTitle();
	AngleDis->GetYaxis()->CenterTitle();
	AngleDis->GetYaxis()->SetTitleOffset(1.4);

	for(int i = 0; i < NEvents; i++){
	
		TruthXY->Fill(EnterX1[i],EnterY1[i]);
		TruthXY->Fill(EnterX2[i],EnterY2[i]);

		ClusXY->Fill(ClusX1[i],ClusY1[i]);
		ClusXY->Fill(ClusX2[i],ClusY2[i]);
	
		ClusZ->Fill(ClusZ1[i] + PVz);
		ClusZ->Fill(ClusZ2[i] + PVz);
	
		EnterR = sqrt((EnterX1[i] - EnterX2[i]) * (EnterX1[i] - EnterX2[i]) + (EnterY1[i] - EnterY2[i]) * (EnterY1[i] - EnterY2[i])) ;
		ClusterR = sqrt((ClusX1[i] - ClusX2[i]) * (ClusX1[i] - ClusX2[i]) + (ClusY1[i] - ClusY2[i]) * (ClusY1[i] -ClusY2[i]) + (ClusZ1[i] - ClusZ2[i]) * (ClusZ1[i] - ClusZ2[i]));
	
		PosClus1->SetXYZ(ClusX1[i],ClusY1[i],ClusZ1[i]);
		PosTruth1->SetXYZ(EnterX1[i],EnterY1[i],EnterZ);
		PosClus2->SetXYZ(ClusX2[i],ClusY2[i],ClusZ2[i]);
		PosTruth2->SetXYZ(EnterX2[i],EnterY2[i],EnterZ);

		* Dis1 = * PosClus1 - *PosTruth1;
		* Dis2 = * PosClus2 - *PosTruth1;
	
		if(Dis1->Mag() < Dis2->Mag())  	CosAngle = PosClus1->Unit().Dot(PosTruth1->Unit());
		if(Dis1->Mag() > Dis2->Mag())  	CosAngle = PosClus1->Unit().Dot(PosTruth2->Unit());


		Gamma1->SetPxPyPzE(ClusPX1[i],ClusPY1[i],ClusPZ1[i],ClusE1[i]);
		Gamma2->SetPxPyPzE(ClusPX2[i],ClusPY2[i],ClusPZ2[i],ClusE2[i]);

		* PiZero = * Gamma1 + *Gamma2;
		

		PiMass = PiZero->M();

		TotalClusE = PiZero->E();
		EAsym = (ClusE1[i] - ClusE2[i])/(ClusE2[i] + ClusE2[i]);
	
			
		AngleDis->Fill(CosAngle);
		TruthR->Fill(EnterR);
		ClusR->Fill(ClusterR);
		ClusE->Fill(TotalClusE);
		ClusEAsym->Fill(EAsym);
		PiMassRECO->Fill(PiMass);

	}





	AngleDis->Draw();
	lat->DrawLatex(0.30,0.80,Form("Incident #pi^{0} Eneregy = %.1f GeV",IncEnergy));	
	lat->DrawLatex(0.30,0.75,"Standard PHENIX Pb Shashlik EMCAL");	
	lat->DrawLatex(0.30,0.70,Form("%s",Comments.Data()));	
	
	c->SaveAs("PosAnaPlots/AngleDis.png");

	
	TCanvas * c2 = new TCanvas("c2","c2",600,600);
	c2->cd();

	TruthXY->Draw("COLZ");
	c2->SaveAs(Form("PosAnaPlots/%s/TruthXY.png",Name.Data()));


	ClusXY->Draw("COLZ");
	c2->SaveAs(Form("PosAnaPlots/%s/ClusXY.png",Name.Data()));

	ClusZ->Draw();

	TLine *line = new TLine(EnterZ,0,EnterZ,ClusZ->GetMaximum());
	line->SetLineStyle(2);
	line->SetLineWidth(2);
	line->SetLineColor(2);
	line->Draw("SAME");
	c2->SaveAs(Form("PosAnaPlots/%s/ClusZ.png",Name.Data()));

	TruthR->Draw();
	c2->SaveAs(Form("PosAnaPlots/%s/TruthR.png",Name.Data()));

	ClusR->Draw();
	c2->SaveAs(Form("PosAnaPlots/%s/ClusR.png",Name.Data()));

	ClusE->Draw();
	lat->DrawLatex(0.30,0.80,Form("Incident #pi^{0} Eneregy = %.1f GeV",IncEnergy));	
	lat->DrawLatex(0.30,0.75,"Standard PHENIX Pb Shashlik EMCAL");	
	lat->DrawLatex(0.30,0.70,Form("%s",Comments.Data()));	
	
	c2->SaveAs(Form("PosAnaPlots/%s/ClusE.png",Name.Data()));

	ClusEAsym->Draw();
	lat->DrawLatex(0.30,0.80,Form("Incident #pi^{0} Eneregy = %.1f GeV",IncEnergy));	
	lat->DrawLatex(0.30,0.75,"Standard PHENIX Pb Shashlik EMCAL");	
	lat->DrawLatex(0.30,0.70,Form("%s",Comments.Data()));	
	c2->SaveAs(Form("PosAnaPlots/%s/ClusEAsym.png",Name.Data()));


	PiMassRECO->Draw();
	lat->DrawLatex(0.30,0.80,Form("Incident #pi^{0} Eneregy = %.1f GeV",IncEnergy));	
	lat->DrawLatex(0.30,0.75,"Standard PHENIX Pb Shashlik EMCAL");	
	lat->DrawLatex(0.30,0.70,Form("%s",Comments.Data()));	
	c2->SaveAs(Form("PosAnaPlots/%s/PiMassRECO.png",Name.Data()));


	TFile * fout = new TFile(Form("PosAna_%s.root",Name.Data()),"RECREATE");
	fout->cd();
	AngleDis->Write();
	TruthXY->Write();
	ClusXY->Write();
	ClusZ->Write();
	TruthR->Write();
	ClusR->Write();
	ClusE->Write();

	fout->Close();
}
