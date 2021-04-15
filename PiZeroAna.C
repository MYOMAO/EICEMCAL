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


void PiZeroAna(int EnergyBin){




	gStyle->SetOptStat(0);


	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();

	TLatex *lat = new TLatex();
	lat->SetNDC();
	lat->SetTextSize(0.035);

	TString infile = Form("PiZeroAnaFile/G4EICDetector.root_g4femc_eval_%d.root",EnergyBin);

	TFile * fin = new TFile(infile.Data());
	fin->cd();

	float EnergyStep = 0.1;
	float IncEnergy = EnergyBin * EnergyStep;
	float PVx = 75;
	float PVy = 75;
	float PVz = 100;

	//float PiZeroMassPDG = 0.13957;
	float MassUpRatio = 1.4;
	float MassLowRatio = 0.0;

	float PiZeroMassPDG = 0.1349766;


	float MassUpRatioCorr = 1.3;
	float MassLowRatioCorr = 0.3;


	float EnergyUpRatio = 1.3;
	float EnergyLowRatio = 0.5;

	int NBins = 100;

	const int NPiEvents = 5000;

	float PiMassUp = PiZeroMassPDG * MassUpRatio;
	float PiMassLow = PiZeroMassPDG * MassLowRatio;


	float PiMassUpCorr = PiZeroMassPDG * MassUpRatioCorr;
	float PiMassLowCorr = PiZeroMassPDG * MassLowRatioCorr;

	float PiEUp = IncEnergy * EnergyUpRatio;
	float PiEDown = IncEnergy * EnergyLowRatio;



	const int NDecay = 2;
	const int NDecayMax = 10;

	float clusterIDFloat;
	float EventFloat;


	int EventPre = 0;
	int Event = 0;



	int clusterID;

	float RECOEnergy;

	float XPos;
	float YPos;
	float ZPos;


	float E[NDecayMax];
	float x[NDecayMax];
	float y[NDecayMax];
	float z[NDecayMax];
	float r[NDecayMax];

	for(int i = 0; i < NDecayMax; i++){

		E[i] = 0;

	}


	float xCorr[NDecay][NPiEvents];
	float yCorr[NDecay][NPiEvents];
	float zCorr[NDecay][NPiEvents];
	float rCorr[NDecay][NPiEvents];

	for(int i = 0; i < NDecay; i ++){


		for(int j = 0; j < NPiEvents; j++){

			xCorr[i][j] = 0;
			yCorr[i][j] = 0;
			zCorr[i][j] = 0;
			rCorr[i][j] = 0;
		}


	}


	E[0] = 0;
	E[1] = 0;

	float Px[NDecayMax];
	float Py[NDecayMax];
	float Pz[NDecayMax];

	float E1;
	float E2;


	float Asymmetry;
	float PiZeroMass;
	float PiZeroEnergy;
	float Eff;

	TLorentzVector * Gamma1 = new TLorentzVector;
	TLorentzVector * Gamma2 = new TLorentzVector;

	TLorentzVector * PiZero = new TLorentzVector;


	float PiZeroMassCorr;

	float PxCorr[NDecayMax];
	float PyCorr[NDecayMax];
	float PzCorr[NDecayMax];

	int NClusterBins = 6;
	float NClusMin = 0;
	float NClusMax = 5;


	TLorentzVector * Gamma1Corr = new TLorentzVector;
	TLorentzVector * Gamma2Corr = new TLorentzVector;

	TLorentzVector * PiZeroCorr = new TLorentzVector;




	TTree * ntp_cluster = (TTree * ) fin->Get("ntp_cluster");

	ntp_cluster->SetBranchAddress("event",&EventFloat);
	ntp_cluster->SetBranchAddress("clusterID",&clusterIDFloat);
	ntp_cluster->SetBranchAddress("e",&RECOEnergy);
	ntp_cluster->SetBranchAddress("x",&XPos);
	ntp_cluster->SetBranchAddress("y",&YPos);
	ntp_cluster->SetBranchAddress("z",&ZPos);


	int NEvents = ntp_cluster->GetEntries();
	float Total = 0;
	float Pass = 0;


	int NClus = 0;

	//Histograms for RECO//

	TH1D * PiMassHis = new TH1D("PiMassHis","",NBins,PiMassLow,PiMassUp);
	PiMassHis->GetXaxis()->SetTitle("Reconstructed m(#pi^{0} -> #gamma #gamma) (GeV/c^{2})");
	PiMassHis->GetYaxis()->SetTitle("Counts");
	PiMassHis->GetXaxis()->CenterTitle();
	PiMassHis->GetYaxis()->CenterTitle();
	PiMassHis->GetYaxis()->SetTitleOffset(1.4);



	TH1D * PiMassHisCorr = new TH1D("PiMassHisCorr","",NBins,PiMassLowCorr,PiMassUpCorr);
	PiMassHisCorr->GetXaxis()->SetTitle("Reconstructed m(#pi^{0} -> #gamma #gamma) (GeV/c^{2}) - Vertex Corrected");
	PiMassHisCorr->GetYaxis()->SetTitle("Counts");
	PiMassHisCorr->GetXaxis()->CenterTitle();
	PiMassHisCorr->GetYaxis()->CenterTitle();
	PiMassHisCorr->GetYaxis()->SetTitleOffset(1.4);


	TH1D * PiEnergyHis = new TH1D("PiEnergyHis","",NBins,PiEDown,PiEUp);
	PiEnergyHis->GetXaxis()->SetTitle("Reconstructed #pi^{0} Energy (GeV)");
	PiEnergyHis->GetYaxis()->SetTitle("Counts");
	PiEnergyHis->GetXaxis()->CenterTitle();
	PiEnergyHis->GetYaxis()->CenterTitle();
	PiEnergyHis->GetYaxis()->SetTitleOffset(1.4);

	TH1D * AsymmetryHis = new TH1D("AsymmetryHis","",NBins,-1,1);
	AsymmetryHis->GetXaxis()->SetTitle("(E_{1} - E_{2})/(E_{1} + E_{2})");
	AsymmetryHis->GetYaxis()->SetTitle("Counts");
	AsymmetryHis->GetXaxis()->CenterTitle();
	AsymmetryHis->GetYaxis()->CenterTitle();
	AsymmetryHis->GetYaxis()->SetTitleOffset(1.4);


	TH1D * NClusterHis = new TH1D("NClusterHis","",NClusterBins,NClusMin - 0.5,NClusMax + 0.5);
	NClusterHis->GetXaxis()->SetTitle("Number of Cluster In One Event");
	NClusterHis->GetYaxis()->SetTitle("Number of Events");
	NClusterHis->GetXaxis()->CenterTitle();
	NClusterHis->GetYaxis()->CenterTitle();
	NClusterHis->GetYaxis()->SetTitleOffset(1.4);


	//	int NBinsXY = 200;
	float WidthXY = 75;

	int NBinsXY = 100;
	int XMin = -200;
	int YMin = -200;
	int XMax = 200;
	int YMax= 200;

	int NBinsR = 100;
	int RMin = 0;
	int RMax = 250;

	TH2D * EXYHis = new TH2D("EXYHis","",NBinsXY,PVx-WidthXY,PVx+WidthXY,NBinsXY,PVy-WidthXY,PVy+WidthXY);
	EXYHis->GetXaxis()->SetTitle("Reco Cluster: X (cm)");
	EXYHis->GetYaxis()->SetTitle("Reco Cluster: Y (cm)");
	EXYHis->GetXaxis()->CenterTitle();
	EXYHis->GetYaxis()->CenterTitle();
	EXYHis->GetYaxis()->SetTitleOffset(1.4);
	EXYHis->SetTitle("EMCAL Shower Energy vs XY Distribution");

	TH2D * ClusXY = new TH2D("EXYHis","",NBinsXY,XMin,XMax,NBinsXY,YMin,YMax);
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




	TH1D * TruthR = new TH1D("TruthR","",NBinsR,RMin,RMax);
	TruthR->GetXaxis()->SetTitle("Distance Between Two Photons on EMCAL Front Face");
	TruthR->GetYaxis()->SetTitle("Number of Events");
	TruthR->GetXaxis()->CenterTitle();
	TruthR->GetYaxis()->CenterTitle();
	TruthR->GetYaxis()->SetTitleOffset(1.4);


	TH1D * TruthE = new TH1D("TruthE","",NBins,PiEDown,PiEDown);
	TruthE->GetXaxis()->SetTitle("Total Energy Deposited to the EMCAL");
	TruthE->GetYaxis()->SetTitle("Number of Events");
	TruthE->GetXaxis()->CenterTitle();
	TruthE->GetYaxis()->CenterTitle();
	TruthE->GetYaxis()->SetTitleOffset(1.4);



	TH1D * ClusR = new TH1D("ClusR","",NBinsR,RMin,RMax);
	ClusR->GetXaxis()->SetTitle("Distance Between Two Clusters Centers");
	ClusR->GetYaxis()->SetTitle("Number of Events");
	ClusR->GetXaxis()->CenterTitle();
	ClusR->GetYaxis()->CenterTitle();
	ClusR->GetYaxis()->SetTitleOffset(1.4);





	//Vertex Corrction// 

	cout << "Now Correcting the Fuckin Vertex" << endl;

	TString infile2 = Form("PiZeroAnaFile/ShowerInfo_%d.root",EnergyBin);

	TFile * fin2 = new TFile(infile2.Data());
	fin2->cd();
	TTree * ShowerInfo = (TTree * ) fin2->Get("ShowerInfo");
	int NEvents2 = ShowerInfo->GetEntries();



	float X1;
	float Y1;
	float Z1;

	float X2;
	float Y2;
	float Z2;

	float XPosCorr1 = 999;
	float YPosCorr1 = 999;
	float ZPosCorr1 = 999;

	float XPosCorr2 = 999;
	float YPosCorr2 = 999;
	float ZPosCorr2 = 999;


	float GPX1;
	float GPX2;
	float GPY1;
	float GPY2;
	float GPZ1;
	float GPZ2;

	int evtid;
	int evtidPre;
	int layer;
	float FEMCALedep;

	float ZMin = 999;

	ShowerInfo->SetBranchAddress("X1",&X1);
	ShowerInfo->SetBranchAddress("X2",&X2);
	ShowerInfo->SetBranchAddress("Y1",&Y1);
	ShowerInfo->SetBranchAddress("Y2",&Y2);
	ShowerInfo->SetBranchAddress("Z1",&Z1);
	ShowerInfo->SetBranchAddress("Z2",&Z2);
	ShowerInfo->SetBranchAddress("evtid",&evtid);
	ShowerInfo->SetBranchAddress("layer",&layer);
	ShowerInfo->SetBranchAddress("FEMCALedep",&FEMCALedep);



	ShowerInfo->SetBranchAddress("GPX1",&GPX1);
	ShowerInfo->SetBranchAddress("GPY1",&GPY1);
	ShowerInfo->SetBranchAddress("GPZ1",&GPZ1);
	ShowerInfo->SetBranchAddress("GPX2",&GPX2);
	ShowerInfo->SetBranchAddress("GPY2",&GPY2);
	ShowerInfo->SetBranchAddress("GPZ2",&GPZ2);

	//Truth Info//

	float EnterZ = 190.27500;
	float TotalTruthE = 0;
	float EnterX1[NPiEvents];
	float EnterY1[NPiEvents];

	float EnterX2[NPiEvents];
	float EnterY2[NPiEvents];

	float EnterR[NPiEvents];

	float ClusDis;


	for(int i = 0; i < NEvents2; i++){


		ShowerInfo->GetEntry(i);


		EnterX1[evtid] = GPX1/GPZ1 * EnterZ;
		EnterX2[evtid] = GPX2/GPZ2 * EnterZ;
		EnterY1[evtid] = GPY1/GPZ1 * EnterZ;
		EnterY2[evtid] = GPY2/GPZ2 * EnterZ;

		EnterR[evtid] = sqrt((EnterX1[evtid] - EnterX2[evtid]) * (EnterX1[evtid] - EnterX2[evtid])  + (EnterY1[evtid] - EnterY2[evtid])  * (EnterY1[evtid] - EnterY2[evtid]) );

		if(evtidPre != evtid){

			TruthE->Fill(TotalTruthE);


			//Real Corr//
			
			XPosCorr1 = EnterX1[evtidPre];
			XPosCorr2 = EnterX2[evtidPre];
			YPosCorr1 = EnterY1[evtidPre];
			YPosCorr2 = EnterY2[evtidPre];
			ZPosCorr1 = EnterZ;
			ZPosCorr2 = EnterZ;

			//Done Modification//

			xCorr[0][evtidPre] = XPosCorr1;
			yCorr[0][evtidPre] = YPosCorr1;
			zCorr[0][evtidPre] = ZPosCorr1;
			rCorr[0][evtidPre] = sqrt(xCorr[0][evtidPre] * xCorr[0][evtidPre] + yCorr[0][evtidPre] * yCorr[0][evtidPre] + zCorr[0][evtidPre] * zCorr[0][evtidPre]);

			xCorr[1][evtidPre] = XPosCorr2;
			yCorr[1][evtidPre] = YPosCorr2;
			zCorr[1][evtidPre] = ZPosCorr2;
			rCorr[1][evtidPre] = sqrt(xCorr[1][evtidPre] * xCorr[1][evtidPre] + yCorr[1][evtidPre] * yCorr[1][evtidPre] + zCorr[1][evtidPre] * zCorr[1][evtidPre]);


			TruthXY->Fill(EnterX1[evtidPre],EnterY1[evtidPre]);
			TruthXY->Fill(EnterX2[evtidPre],EnterY2[evtidPre]);

			TruthR->Fill(EnterR[evtidPre]);

			ZMin = 999;

			TotalTruthE = 0;

		}


		if(ZMin > (Z1+Z2)/2 && layer == -1){
			ZMin = (Z1+Z2)/2;

			if((X1+X2)/2 < PVx){
				XPosCorr1 = (X1+X2)/2 - PVx;
				YPosCorr1 = (Y1+Y2)/2 - PVy;
				ZPosCorr1 = (Z1+Z2)/2 - PVz;				
			}

			if((X1+X2)/2 > PVx){
				XPosCorr2 = (X1+X2)/2 - PVx;
				YPosCorr2 = (Y1+Y2)/2 - PVy;
				ZPosCorr2 = (Z1+Z2)/2- PVz;				
			}
		}

		TotalTruthE = TotalTruthE + FEMCALedep;


		evtidPre = evtid;
	}

	cout << "DONE CORRECTING THE VERTEX" << endl;

	for(int i = 0; i < NDecay; i ++){


		for(int j = 0; j < NPiEvents; j++){

			//cout << "xCorr[i][j] = " << xCorr[i][j]   << "    yCorr[i][j] = " << yCorr[i][j]  << "   zCorr[i][j] = " << zCorr[i][j] << endl;

		}

	}

	//DONE Vertex Correction//



	for(int i = 0; i < NEvents; i++){

//		cout << "Now Working on FUCKIN Event = " << i << "   clusterID = " << clusterID << endl;

		ntp_cluster->GetEntry(i);

		Event = int(EventFloat);
		clusterID = int(clusterIDFloat);


		//RECOMSTRUCTION AND FILL HISTOGRAM//	
		if(EventPre !=  Event){
			Total = Total + 1;	
			//		if(E[0] == 0 || E[1] == 0) continue; 
			if(E[0] > 0 && E[1] > 0 && E[2] == 0 && E[3] == 0 && E[4] == 0){
				Pass = Pass + 1;

				Gamma1->SetPxPyPzE(Px[0],Py[0],Pz[0],E[0]);
				Gamma2->SetPxPyPzE(Px[1],Py[1],Pz[1],E[1]);

				* PiZero = * Gamma1 + * Gamma2;

				PiZeroMass = PiZero->M();
				PiZeroEnergy = PiZero->E();

				if(PiZeroMass < 0.06) cout << "FUCKED EVENT: " << EventPre << endl;

				PiMassHis->Fill(PiZeroMass);
				PiEnergyHis->Fill(PiZeroEnergy);


				EXYHis->Fill(x[0] + PVx,y[0] + PVy,E[0]);
				EXYHis->Fill(x[1] + PVx,y[1] + PVy,E[1]);

				ClusXY->Fill(x[0],y[0]);
				ClusXY->Fill(x[1],y[1]);

				ClusDis = sqrt((x[0] - x[1]) * (x[0] - x[1]) + (y[0] - y[1]) * (y[0] - y[1]));

				ClusR->Fill(ClusDis);

				E1 = Gamma1->E();
				E2 = Gamma2->E();
				Asymmetry = (E1-E2)/(E1+E2);


				AsymmetryHis->Fill(Asymmetry);
				//cout << "PiZeroMass = " << PiZeroMass << "   PiZeroEnergy =  " << PiZeroEnergy << "    Asymmetry = " << Asymmetry << endl;

				//Corrected Inv Mass//

				Gamma1Corr->SetPxPyPzE(PxCorr[0],PyCorr[0],PzCorr[0],E[0]);
				Gamma2Corr->SetPxPyPzE(PxCorr[1],PyCorr[1],PzCorr[1],E[1]);

				* PiZeroCorr = * Gamma1Corr + * Gamma2Corr;

				PiZeroMassCorr = PiZeroCorr->M();
				//cout << "PiZeroMassCorr = " << PiZeroMassCorr << endl;
				PiMassHisCorr->Fill(PiZeroMassCorr);




			}



			for(int i = 0; i < NDecayMax; i++){

				E[i] = 0;

			}

			NClusterHis->Fill(NClus);
			NClus = 0;
		}


		E[clusterID] = RECOEnergy;
		x[clusterID] = XPos - PVx;
		y[clusterID] = YPos - PVy;
		z[clusterID] = ZPos - PVz;
		z[clusterID] = EnterZ - 5;  //PHENIX Pb EMCAL Front Face z coordinate
		r[clusterID] = sqrt(x[clusterID]*x[clusterID] + y[clusterID]*y[clusterID] + z[clusterID] * z[clusterID]);


		Px[clusterID] = E[clusterID]/r[clusterID]*x[clusterID];
		Py[clusterID] = E[clusterID]/r[clusterID]*y[clusterID];
		Pz[clusterID] = E[clusterID]/r[clusterID]*z[clusterID];



		if(clusterID < 2){
			PxCorr[clusterID] = E[clusterID]/rCorr[clusterID][Event]*xCorr[clusterID][Event];
			PyCorr[clusterID] = E[clusterID]/rCorr[clusterID][Event]*yCorr[clusterID][Event];
			PzCorr[clusterID] = E[clusterID]/rCorr[clusterID][Event]*zCorr[clusterID][Event];

			//	cout << "xRECO =  " <<  x[clusterID] << "   xTruth = " << xCorr[clusterID][Event] << endl;
			//	cout << "yRECO =  " <<  y[clusterID] << "   yTruth = " << yCorr[clusterID][Event] << endl;
			//	cout << "zRECO =  " <<  z[clusterID] << "   zTruth = " << zCorr[clusterID][Event] << endl;

		}

		NClus = NClus + 1;

		EventPre = Event;
	}

	Eff = Pass/Total;

	cout << "Pion Reconstruction Efficiency = " << Eff << endl;


	//Draw Shits//


	PiMassHis->Draw("hist");
	lat->DrawLatex(0.30,0.80,Form("Incident #pi^{0} Eneregy = %.1f GeV",IncEnergy));	
	lat->DrawLatex(0.30,0.75,"Standard PHENIX Pb Shashlik EMCAL");	
	




	TLine *line = new TLine(PiZeroMassPDG,0,PiZeroMassPDG,PiMassHis->GetMaximum());
	line->SetLineStyle(2);
	line->SetLineWidth(2);
	line->SetLineColor(2);
	line->Draw("SAME");

	c->SaveAs(Form("PiZeroPlots/%d/Mass.png",EnergyBin));


	PiMassHisCorr->Draw("hist");
	lat->DrawLatex(0.30,0.80,Form("Incident #pi^{0} Eneregy = %.1f GeV",IncEnergy));	
	lat->DrawLatex(0.30,0.75,"Standard PHENIX Pb Shashlik EMCAL");	

	TLine *line2 = new TLine(PiZeroMassPDG,0,PiZeroMassPDG,PiMassHisCorr->GetMaximum());
	line2->SetLineStyle(2);
	line2->SetLineWidth(2);
	line2->SetLineColor(2);
	line2->Draw("SAME");

	
	c->SaveAs(Form("PiZeroPlots/%d/MassCorr.png",EnergyBin));

	PiEnergyHis->Draw("hist");	
	lat->DrawLatex(0.30,0.80,Form("Incident #pi^{0} Eneregy = %.1f GeV",IncEnergy));	
	lat->DrawLatex(0.30,0.75,"Standard PHENIX Pb Shashlik EMCAL");	

	c->SaveAs(Form("PiZeroPlots/%d/Energy.png",EnergyBin));

	AsymmetryHis->Draw("hist");	
	lat->DrawLatex(0.30,0.80,Form("Incident #pi^{0} Eneregy = %.1f GeV",IncEnergy));	
	lat->DrawLatex(0.30,0.75,"Standard PHENIX Pb Shashlik EMCAL");	

	c->SaveAs(Form("PiZeroPlots/%d/Asymmetry.png",EnergyBin));


	EXYHis->Draw("COLZ");
	c->SaveAs(Form("PiZeroPlots/%d/EXYHis.png",EnergyBin));


	ClusXY->Draw("COLZ");
	c->SaveAs(Form("PiZeroPlots/%d/ClusXY.png",EnergyBin));

	TruthXY->Draw("COLZ");
	c->SaveAs(Form("PiZeroPlots/%d/TruthXY.png",EnergyBin));


	TruthR->Draw();
	lat->DrawLatex(0.30,0.80,Form("Incident #pi^{0} Eneregy = %.1f GeV",IncEnergy));	
	lat->DrawLatex(0.30,0.75,"Standard PHENIX Pb Shashlik EMCAL");	
	c->SaveAs(Form("PiZeroPlots/%d/TruthR.png",EnergyBin));

	TruthE->Draw();
	lat->DrawLatex(0.30,0.80,Form("Incident #pi^{0} Eneregy = %.1f GeV",IncEnergy));	
	lat->DrawLatex(0.30,0.75,"Standard PHENIX Pb Shashlik EMCAL");	
	c->SaveAs(Form("PiZeroPlots/%d/TruthE.png",EnergyBin));

	ClusR->Draw();
	lat->DrawLatex(0.30,0.80,Form("Incident #pi^{0} Eneregy = %.1f GeV",IncEnergy));	
	lat->DrawLatex(0.30,0.75,"Standard PHENIX Pb Shashlik EMCAL");	
	c->SaveAs(Form("PiZeroPlots/%d/ClusR.png",EnergyBin));


	NClusterHis->Draw();
	lat->DrawLatex(0.30,0.80,Form("Incident #pi^{0} Eneregy = %.1f GeV",IncEnergy));	
	lat->DrawLatex(0.30,0.75,"Standard PHENIX Pb Shashlik EMCAL");	
	c->SaveAs(Form("PiZeroPlots/%d/NClusterHis.png",EnergyBin));


	float PiMassBad = 0.06;
	float ClusRBad = 48;	
	int XBinMax;
	int NBadEvents;

	XBinMax = PiMassHis->GetXaxis()->FindBin(PiMassBad);
	NBadEvents = PiMassHis->Integral(0,XBinMax);
	cout << "PiMass: NBadEvents = " << NBadEvents << endl;

	XBinMax = ClusR->GetXaxis()->FindBin(ClusRBad);
	NBadEvents = ClusR->Integral(0,XBinMax);
	cout << "ClusR: NBadEvents = " << NBadEvents << endl;

	


}
