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


void PiZeroAnaGen(int EnergyBin){




	gStyle->SetOptStat(0);


	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();

	TLatex *lat = new TLatex();
	lat->SetNDC();
	lat->SetTextSize(0.035);
	
	float eMass = 0.000511;

	float EnergyStep = 0.1;
	float IncEnergy = EnergyBin * EnergyStep;
	float PVx = 75;
	float PVy = 75;
	float PVz = 100;
	
//	float shift = 120;
	float EnterZ = 190.27500;

	float PiZeroMassPDG = 0.13957;
	float MassUpRatio = 1.3;
	float MassLowRatio = 0.3;

	int NumOfE;
	
	float MassUpRatioCorr = 1.3;
	float MassLowRatioCorr = 0.3;


	float EnergyUpRatio = 1.02;
	float EnergyLowRatio = 0.98;

	int NBins = 100;

	const int NPiEvents = 3000;

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

	float XPos;PiZeroAnaGen.C
	float YPos;
	float ZPos;


	float E[NDecayMax];
	float x[NDecayMax];
	float y[NDecayMax];
	float z[NDecayMax];
	float r[NDecayMax];

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

	TLorentzVector * Gamma1Corr = new TLorentzVector;
	TLorentzVector * Gamma2Corr = new TLorentzVector;

	TLorentzVector * PiZeroCorr = new TLorentzVector;

	//

	TLorentzVector * Gamma1Truth = new TLorentzVector;
	TLorentzVector * Gamma2Truth = new TLorentzVector;

	TLorentzVector * PiZeroTruth = new TLorentzVector;
	

	float PiZeroTruthMass;


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


	TH1D * PiMassHisTruth = new TH1D("PiMassHisTruth","",NBins,PiMassLowCorr,PiMassUpCorr);
	PiMassHisTruth->GetXaxis()->SetTitle("Total #pi^{0} Mass (MeV/c^{2})");
	PiMassHisTruth->GetYaxis()->SetTitle("Counts");
	PiMassHisTruth->GetXaxis()->CenterTitle();
	PiMassHisTruth->GetYaxis()->CenterTitle();
	PiMassHisTruth->GetYaxis()->SetTitleOffset(1.4);


	TH1D * PiEnergyHisCorr = new TH1D("PiEnergyHisCorr","",NBins,0.95,1.01);
	PiEnergyHisCorr->GetXaxis()->SetTitle("Total #pi^{0} Energy (GeV)");
	PiEnergyHisCorr->GetYaxis()->SetTitle("Counts");
	PiEnergyHisCorr->GetXaxis()->CenterTitle();
	PiEnergyHisCorr->GetYaxis()->CenterTitle();
	PiEnergyHisCorr->GetYaxis()->SetTitleOffset(1.4);


	TH1D * PiEnergyHisBeforeCorr = new TH1D("PiEnergyHisBeforeCorr","",NBins,0,1);
	PiEnergyHisBeforeCorr->GetXaxis()->SetTitle("Total #pi^{0} Energy (GeV)");
	PiEnergyHisBeforeCorr->GetYaxis()->SetTitle("Counts");
	PiEnergyHisBeforeCorr->GetXaxis()->CenterTitle();
	PiEnergyHisBeforeCorr->GetYaxis()->CenterTitle();
	PiEnergyHisBeforeCorr->GetYaxis()->SetTitleOffset(1.4);






	TH1D * PiEnergyHis = new TH1D("PiEnergyHis","",NBins,PiEDown,PiEUp);
	PiEnergyHis->GetXaxis()->SetTitle("Reconstructed #pi^{0} Energy (GeV)");
	PiEnergyHis->GetYaxis()->SetTitle("Counts");
	PiEnergyHis->GetXaxis()->CenterTitle();
	PiEnergyHis->GetYaxis()->CenterTitle();
	PiEnergyHis->GetYaxis()->SetTitleOffset(1.4);



	TH1D * GammaEHis = new TH1D("GammaEHis","",NBins,0.90,1.10);
	GammaEHis->GetXaxis()->SetTitle("Potons Energy Sum: E_{1} + E_{2}");
	GammaEHis->GetYaxis()->SetTitle("Number of Events");
	GammaEHis->GetXaxis()->CenterTitle();
	GammaEHis->GetYaxis()->CenterTitle();
	GammaEHis->GetYaxis()->SetTitleOffset(1.4);

	TH1D * AsymmetryHis = new TH1D("AsymmetryHis","",NBins,-1,1);
	AsymmetryHis->GetXaxis()->SetTitle("(E_{1} - E_{2})/(E_{1} + E_{2})");
	AsymmetryHis->GetYaxis()->SetTitle("Counts");
	AsymmetryHis->GetXaxis()->CenterTitle();
	AsymmetryHis->GetYaxis()->CenterTitle();
	AsymmetryHis->GetYaxis()->SetTitleOffset(1.4);


	TH1D * TotalNEDistribution = new TH1D("TotalNEDistribution","",1000,0,1000);
	TotalNEDistribution->GetXaxis()->SetTitle("Number of Electrons Per Events");
	TotalNEDistribution->GetYaxis()->SetTitle("Number of Events");
	TotalNEDistribution->GetXaxis()->CenterTitle();
	TotalNEDistribution->GetYaxis()->CenterTitle();
	TotalNEDistribution->GetYaxis()->SetTitleOffset(1.4);


	int NBinsXY = 50;
	float WidthXY = 30;

	int NBinsR = 200;
	float WidthR = 20;
	float MeanR = 25.6;


	TH2D * EXYHis = new TH2D("EXYHis","",NBinsXY,PVx-WidthXY,PVx+WidthXY,NBinsXY,PVy-WidthXY,PVy+WidthXY);
	EXYHis->GetXaxis()->SetTitle("X (cm)");
	EXYHis->GetYaxis()->SetTitle("Y (cm)");
	EXYHis->GetXaxis()->CenterTitle();
	EXYHis->GetYaxis()->CenterTitle();
	EXYHis->GetYaxis()->SetTitleOffset(1.4);
	EXYHis->SetTitle("EMCAL Shower Energy vs XY Distribution");


	TH2D * XYHisDis = new TH2D("XYHisDis","",NBinsXY,-WidthXY,WidthXY,NBinsXY,-WidthXY,+WidthXY);
	XYHisDis->GetXaxis()->SetTitle("X (cm)");
	XYHisDis->GetYaxis()->SetTitle("Y (cm)");
	XYHisDis->GetXaxis()->CenterTitle();
	XYHisDis->GetYaxis()->CenterTitle();
	XYHisDis->GetYaxis()->SetTitleOffset(1.4);
	XYHisDis->SetTitle("EMCAL Front Face Phton XY Distribution For Each Event");



	TH1D * RHisDis = new TH1D("RHisDis","",NBinsR,0,100);
	RHisDis->GetXaxis()->SetTitle("Distance Between the Two Photon on the EMCAL Front Face (cm)");
	RHisDis->GetYaxis()->SetTitle("Number of Events");
	RHisDis->GetXaxis()->CenterTitle();
	RHisDis->GetYaxis()->CenterTitle();
	RHisDis->GetYaxis()->SetTitleOffset(1.4);




	TH1D * RHisDisCorr = new TH1D("RHisDisCorr","",NBinsR,0,100);
	RHisDisCorr->GetXaxis()->SetTitle("Distance Between the Two Photon on the EMCAL Front Face - From Shower (cm)");
	RHisDisCorr->GetYaxis()->SetTitle("Number of Events");
	RHisDisCorr->GetXaxis()->CenterTitle();
	RHisDisCorr->GetYaxis()->CenterTitle();
	RHisDisCorr->GetYaxis()->SetTitleOffset(1.4);



	TH1D * ClusEDev = new TH1D("ClusEDev","",100,-1,1); 
	ClusEDev->GetXaxis()->SetTitle("(E_{1} - E_{#gamma 1})/E_{#gamma 1} Or (E_{2} - E_{#gamma 2})/E_{#gamma 2}");
	ClusEDev->GetYaxis()->SetTitle("Number of Events");
	ClusEDev->GetXaxis()->CenterTitle();
	ClusEDev->GetYaxis()->CenterTitle();
	ClusEDev->GetYaxis()->SetTitleOffset(1.4);


	//Vertex Corrction// 
	
	cout << "Now Correcting the Fuckin Vertex" << endl;

	TString infile2 = Form("PiZeroAnaFileGen/ShowerInfo_%d.root",EnergyBin);

	TFile * fin2 = new TFile(infile2.Data());
	fin2->cd();PiZeroAnaGen.C
	TTree * ShowerInfo = (TTree * ) fin2->Get("ShowerInfo");
	int NEvents2 = ShowerInfo->GetEntries();


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


	float E1FracDev;
	float E2FracDev;

	float EnterX1;
	float EnterY1;

	float EnterX2;
	float EnterY2;
	float EnterR;

	float TotalR1;
	float TotalR2;

	float XPosCorr1 = 999;
	float YPosCorr1 = 999;
	float ZPosCorr1 = 999;
	
	float XPosCorr2 = 999;
	float YPosCorr2 = 999;
	float ZPosCorr2 = 999;
	

	float RCorr;

	int evtid = 0;
	int evtidPre = 0;
	int layer;
		
	
	float BHBLedep;
	float BHFLedep;
	
	float PiZeroMassTruth;

	float GE1;
	float GE2;


	float GE1Ana[NPiEvents];
	float GE2Ana[NPiEvents];

	float GPX1Ana[NPiEvents];
	float GPY1Ana[NPiEvents];
	float GPZ1Ana[NPiEvents];

	float GPX2Ana[NPiEvents];
	float GPY2Ana[NPiEvents];
	float GPZ2Ana[NPiEvents];


	float ZMin = 999;

	ShowerInfo->SetBranchAddress("X1",&X1);
	ShowerInfo->SetBranchAddress("X2",&X2);
	ShowerInfo->SetBranchAddress("Y1",&Y1);
	ShowerInfo->SetBranchAddress("Y2",&Y2);
	ShowerInfo->SetBranchAddress("Z1",&Z1);
	ShowerInfo->SetBranchAddress("Z2",&Z2);
	ShowerInfo->SetBranchAddress("evtid",&evtid);
	ShowerInfo->SetBranchAddress("layer",&layer);

	ShowerInfo->SetBranchAddress("GPX1",&GPX1);
	ShowerInfo->SetBranchAddress("GPY1",&GPY1);
	ShowerInfo->SetBranchAddress("GPZ1",&GPZ1);
	ShowerInfo->SetBranchAddress("GPX2",&GPX2);
	ShowerInfo->SetBranchAddress("GPY2",&GPY2);
	ShowerInfo->SetBranchAddress("GPZ2",&GPZ2);

	ShowerInfo->SetBranchAddress("GE1",&GE1);
	ShowerInfo->SetBranchAddress("GE2",&GE2);

	ShowerInfo->SetBranchAddress("NumOfE",&NumOfE);


	ShowerInfo->SetBranchAddress("FEMCALedep",&FEMCALedep);
//	ShowerInfo->SetBranchAddress("FEMCALLYDep",&FEMCALLYDep);
	ShowerInfo->SetBranchAddress("BHBLedep",&BHBLedep);
	ShowerInfo->SetBranchAddress("BHFLedep",&BHFLedep);
	
	float TotalEnergy  = 0;
	

	float TotalX1 = 0;
	float TotalY1 = 0;
	float TotalZ1 = 0;

	float TotalX2 = 0;
	float TotalY2 = 0;
	float TotalZ2 = 0;


	float MeanX1 = 0;
	float MeanY1 = 0;
	float MeanZ1 = 0;

	float MeanX2 = 0;
	float MeanY2 = 0;
	float MeanZ2 = 0;
	

	float TotalN1 = 0;
	float TotalN2 = 0;


	float ToalELeak = 0;


	for(int i = 0; i < NEvents2; i++){
		
		ShowerInfo->GetEntry(i);
		GE1Ana[evtid] = GE1;
		GE2Ana[evtid] = GE2;
		GPX1Ana[evtid] = GPX1;
		GPX2Ana[evtid] = GPX2;
		GPY1Ana[evtid] = GPY1;
		GPY2Ana[evtid] = GPY2;
		GPZ1Ana[evtid] = GPZ1;
		GPZ2Ana[evtid] = GPZ2;
		
	}
	evtid  = 0;

	for(int i = 0; i < NEvents2; i++){
				
		ShowerInfo->GetEntry(i);
	

		
		if(evtidPre != evtid ){
			
			TotalNEDistribution->Fill(NumOfE);
			PiEnergyHisBeforeCorr->Fill(TotalEnergy);
		
			//cout << "TotalEnergy = " << TotalEnergy << endl;
			PiEnergyHisCorr->Fill(TotalEnergy + NumOfE * eMass);
			//cout << "E[0] = " << E[0]  << endl;
			//cout << "E[1] = " << E[1]  << endl;
			

			//cout << "TotalN1 = " <<  TotalN1 << "   TotalN2 = " <<  TotalN2 << endl;
			MeanX1 = TotalX1/TotalN1;
			MeanY1 = TotalY1/TotalN1;
			MeanZ1 = TotalZ1/TotalN1;
	
			MeanX2 = TotalX2/TotalN2;
			MeanY2 = TotalY2/TotalN2;
			MeanZ2 = TotalZ2/TotalN2;
			

			x[0] = MeanX1;
			y[0] = MeanY1;
			z[0] = EnterZ;
			r[0] = sqrt(x[0] * x[0] +y[0]  * y[0] + z[0] * z[0]);

			x[1] = MeanX2;
			y[1] = MeanY2;
			z[1]  = EnterZ;
			r[1] = sqrt(x[1] * x[1] +y[1]  * y[1] + z[1] * z[1]);
				
			Px[0] = E[0]/r[0]*x[0];
			Py[0] = E[0]/r[0]*y[0];
			Pz[0] = E[0]/r[0]*z[0];
			
			Px[1] = E[1]/r[1]*x[1];
			Py[1] = E[1]/r[1]*y[1];
			Pz[1] = E[1]/r[1]*z[1];


			Gamma1->SetPxPyPzE(Px[0],Py[0],Pz[0],E[0]);
			Gamma2->SetPxPyPzE(Px[1],Py[1],Pz[1],E[1]);

			* PiZero = * Gamma1 + * Gamma2;
			
			PiZeroMass = PiZero->M();
			PiMassHis->Fill(PiZeroMass);
		

			EXYHis->Fill(MeanX1 ,MeanY1,E[0]);
			EXYHis->Fill(MeanX2 ,MeanY2,E[1]);
	

			EnterX1 = GPX1/GPZ1 * EnterZ;
			EnterY1 = GPY1/GPZ1 * EnterZ;


			EnterX2 = GPX2/GPZ2 * EnterZ;
			EnterY2 = GPY2/GPZ2 * EnterZ;
	
		

			//cout << "EnterX1 = " << EnterX1 << "   EnterY1 = " << EnterY1 << "  EnterX2 = " << EnterX2 << "  EnterY2 = " << EnterY2 << endl;

			EnterR = sqrt((EnterX1 -EnterX2) * (EnterX1 -EnterX2) + (EnterY1 -EnterY2) * (EnterY1 -EnterY2));

			TotalR1 = sqrt(EnterX1 * EnterX1 + EnterY1 * EnterY1 + EnterZ * EnterZ);
			TotalR2 = sqrt(EnterX2 * EnterX2 + EnterY2 * EnterY2 + EnterZ * EnterZ);

			XYHisDis->Fill(EnterX1,EnterY1);
			XYHisDis->Fill(EnterX2,EnterY2);
			
			RHisDis->Fill(EnterR);
			GammaEHis->Fill(GE1Ana[evtid-1] + GE2Ana[evtid-1]);

		//	Gamma1Truth->SetPxPyPzE(E[0]* EnterX1/TotalR1,E[0]* EnterY1/TotalR1,E[0]* EnterZ/TotalR1,E[0]);

		//	Gamma2Truth->SetPxPyPzE(E[1]* EnterX1/TotalR2,E[1]* EnterY1/TotalR2,E[1]* EnterZ/TotalR2,E[1]);
		




			if(TMath::Abs(E[0] - GE1Ana[evtid-1]) < TMath::Abs(E[1] - GE1Ana[evtid-1]) ){
				

				E1 = E[0];
				E2 = E[1];
				
			}
			if(TMath::Abs(E[0] - GE1Ana[evtid-1]) > TMath::Abs(E[1] - GE1Ana[evtid-1]) ){

				E1 = E[1];
				E2 = E[0];
			}

			cout << "Check: E1 " <<  E1 << "   GE1 = " << GE1Ana[evtid-1] <<  "   E2 = " << E2 << "   GE2 = " << GE2Ana[evtid-1] << endl;
			
		//	Gamma1Truth->SetPxPyPzE(GPX1Ana[evtid-1]/1000,GPY1Ana[evtid-1]/1000,GPZ1Ana[evtid-1]/1000,GE1Ana[evtid-1]);
		//	Gamma2Truth->SetPxPyPzE(GPX2Ana[evtid-1]/1000,GPY2Ana[evtid-1]/1000,GPZ2Ana[evtid-1]/1000,GE2Ana[evtid-1] );
			Gamma1Truth->SetPxPyPzE(GPX1Ana[evtid-1]/1000,GPY1Ana[evtid-1]/1000,GPZ1Ana[evtid-1]/1000,E1);
			Gamma2Truth->SetPxPyPzE(GPX2Ana[evtid-1]/1000,GPY2Ana[evtid-1]/1000,GPZ2Ana[evtid-1]/1000,E2);

		
			*PiZeroTruth = *Gamma1Truth + *Gamma2Truth;
	
			PiZeroTruthMass = PiZeroTruth->M();
				
			PiMassHisTruth->Fill(PiZeroTruthMass);
			
		
			E1FracDev = (E1 - GE1Ana[evtid-1])/GE1Ana[evtid-1];
			E2FracDev = (E2 - GE2Ana[evtid-1])/GE2Ana[evtid-1];

	
			ClusEDev->Fill(E1FracDev);
			ClusEDev->Fill(E2FracDev);
			
			
			
			
			
			
			
			//Min Corr
		

			RCorr = sqrt( (XPosCorr1 - XPosCorr2) * (XPosCorr1 - XPosCorr2) + (YPosCorr1 - YPosCorr2) * (YPosCorr1 - YPosCorr2));
			RHisDisCorr->Fill(RCorr);

			xCorr[0][evtidPre] = XPosCorr1;
			yCorr[0][evtidPre] = YPosCorr1;
			zCorr[0][evtidPre] = ZPosCorr1;
			rCorr[0][evtidPre] = sqrt(xCorr[0][evtidPre] * xCorr[0][evtidPre] + yCorr[0][evtidPre] * yCorr[0][evtidPre] + zCorr[0][evtidPre] * zCorr[0][evtidPre]);

			xCorr[1][evtidPre] = XPosCorr2;
			yCorr[1][evtidPre] = YPosCorr2;
			zCorr[1][evtidPre] = ZPosCorr2;
			rCorr[1][evtidPre] = sqrt(xCorr[1][evtidPre] * xCorr[1][evtidPre] + yCorr[1][evtidPre] * yCorr[1][evtidPre] + zCorr[1][evtidPre] * zCorr[1][evtidPre]);
	
			PxCorr[0] = E[0]/rCorr[0][evtidPre]*xCorr[0][evtidPre];
			PyCorr[0] = E[0]/rCorr[0][evtidPre]*yCorr[0][evtidPre];
			PzCorr[0] = E[0]/rCorr[0][evtidPre]*zCorr[0][evtidPre];

	
			PxCorr[1] = E[1]/rCorr[1][evtidPre]*xCorr[1][evtidPre];
			PyCorr[1] = E[1]/rCorr[1][evtidPre]*yCorr[1][evtidPre];
			PzCorr[1] = E[1]/rCorr[1][evtidPre]*zCorr[1][evtidPre];

			
			Gamma1Corr->SetPxPyPzE(PxCorr[0],PyCorr[0],PzCorr[0],E[0]);
			Gamma2Corr->SetPxPyPzE(PxCorr[1],PyCorr[1],PzCorr[1],E[1]);

			* PiZeroCorr = * Gamma1Corr + * Gamma2Corr;
			
			PiZeroMassCorr = PiZeroCorr->M();
			//cout << "PiZeroMassCorr = " << PiZeroMassCorr << endl;
			PiMassHisCorr->Fill(PiZeroMassCorr);
	
			
			PiZeroEnergy = PiZeroCorr->E();
			PiEnergyHis->Fill(PiZeroEnergy);

			//cout << "PiZeroEnergy = " << PiZeroEnergy << endl;

			E1 = Gamma1Corr->E();
			E2 = Gamma2Corr->E();
			Asymmetry = (E1-E2)/(E1+E2);
			AsymmetryHis->Fill(Asymmetry);
			
			ZMin = 999;
			TotalEnergy  = 0;
			E[0] = 0;
			E[1] = 0;
			
			TotalX1 = 0;
			TotalY1 = 0;
			TotalZ1 = 0;
			TotalX2 = 0;
			TotalY2 = 0;
			TotalZ2 = 0;
			
			TotalN1 = 0;
			TotalN2 = 0;
			ToalELeak = 0;

			evtidPre = evtid;

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

			
		if((X1+X2)/2 < PVx&& layer == -1 ){

			TotalX1 = (X1+X2)/2 + TotalX1;
			TotalY1 = (Y1+Y2)/2 + TotalY1;
			TotalZ1 = (Z1+Z2)/2 + TotalZ1;
			TotalN1 = TotalN1 + 1.0;
			//E[0] = FEMCALedep + BHBLedep + BHFLedep + E[0];
			E[0] = FEMCALedep + E[0];

		}
		TotalEnergy =  FEMCALedep + BHBLedep + BHFLedep + TotalEnergy;

	
		if((X1+X2)/2 > PVx && layer == -1){

			TotalX2 = (X1+X2)/2 + TotalX2;
			TotalY2 = (Y1+Y2)/2 + TotalY2;
			TotalZ2 = (Z1+Z2)/2 + TotalZ2;
			TotalN2 = TotalN2 + 1.0;
	//		E[1] = FEMCALedep + BHBLedep + BHFLedep + E[1];				
			E[1] = FEMCALedep + E[1];
		
		}


		
		
	}



	PiMassHis->Draw("hist");
	lat->DrawLatex(0.30,0.80,Form("Incident #pi^{0} Eneregy = %.1f GeV",IncEnergy));	
	lat->DrawLatex(0.30,0.75,"Standard PHENIX Pb Shashlik EMCAL");	

	c->SaveAs(Form("PiZeroPlotsGen/%d/Mass.png",EnergyBin));

	PiMassHisCorr->Draw("hist");
	lat->DrawLatex(0.30,0.80,Form("Incident #pi^{0} Eneregy = %.1f GeV",IncEnergy));	
	lat->DrawLatex(0.30,0.75,"Standard PHENIX Pb Shashlik EMCAL");	

	c->SaveAs(Form("PiZeroPlotsGen/%d/MassCorr.png",EnergyBin));

	PiEnergyHis->Draw("hist");	
	lat->DrawLatex(0.30,0.80,Form("Incident #pi^{0} Eneregy = %.1f GeV",IncEnergy));	
	lat->DrawLatex(0.30,0.75,"Standard PHENIX Pb Shashlik EMCAL");	

	c->SaveAs(Form("PiZeroPlotsGen/%d/Energy.png",EnergyBin));

	PiEnergyHisBeforeCorr->Draw("hist");	
	lat->DrawLatex(0.30,0.80,Form("Incident #pi^{0} Eneregy = %.1f GeV",IncEnergy));	
	lat->DrawLatex(0.30,0.75,"Standard PHENIX Pb Shashlik EMCAL");	

	c->SaveAs(Form("PiZeroPlotsGen/%d/EnergyBeforeCorr.png",EnergyBin));


	PiEnergyHisCorr->Draw("hist");	
	lat->DrawLatex(0.30,0.80,Form("Incident #pi^{0} Eneregy = %.1f GeV",IncEnergy));	
	lat->DrawLatex(0.30,0.75,"Standard PHENIX Pb Shashlik EMCAL");	

	c->SaveAs(Form("PiZeroPlotsGen/%d/EnergyCorr.png",EnergyBin));

	AsymmetryHis->Draw("hist");	
	lat->DrawLatex(0.30,0.80,Form("Incident #pi^{0} Eneregy = %.1f GeV",IncEnergy));	
	lat->DrawLatex(0.30,0.75,"Standard PHENIX Pb Shashlik EMCAL");	

	c->SaveAs(Form("PiZeroPlotsGen/%d/Asymmetry.png",EnergyBin));
	

	EXYHis->Draw("COLZ");
	c->SaveAs(Form("PiZeroPlotsGen/%d/EXYHis.png",EnergyBin));


	XYHisDis->Draw("COLZ");
	c->SaveAs(Form("PiZeroPlotsGen/%d/XYHisDis.png",EnergyBin));
	

	RHisDis->Draw();
	lat->DrawLatex(0.30,0.80,Form("Incident #pi^{0} Eneregy = %.1f GeV",IncEnergy));	
	lat->DrawLatex(0.30,0.75,"Standard PHENIX Pb Shashlik EMCAL");	
	c->SaveAs(Form("PiZeroPlotsGen/%d/RHisDis.png",EnergyBin));


	RHisDisCorr->Draw();
	lat->DrawLatex(0.30,0.80,Form("Incident #pi^{0} Eneregy = %.1f GeV",IncEnergy));	
	lat->DrawLatex(0.30,0.75,"Standard PHENIX Pb Shashlik EMCAL");	
	c->SaveAs(Form("PiZeroPlotsGen/%d/RHisDisCorr.png",EnergyBin));


	GammaEHis->Draw();
	lat->DrawLatex(0.30,0.80,Form("Incident #pi^{0} Eneregy = %.1f GeV",IncEnergy));	
	lat->DrawLatex(0.30,0.75,"Standard PHENIX Pb Shashlik EMCAL");	
	c->SaveAs(Form("PiZeroPlotsGen/%d/GammaEHis.png",EnergyBin));
	

	PiMassHisTruth->Draw();
	lat->DrawLatex(0.30,0.80,Form("Incident #pi^{0} Eneregy = %.1f GeV",IncEnergy));	
	lat->DrawLatex(0.30,0.75,"Standard PHENIX Pb Shashlik EMCAL");	
	c->SaveAs(Form("PiZeroPlotsGen/%d/PiMassHisTruth.png",EnergyBin));


	

	TotalNEDistribution->Draw();
	lat->DrawLatex(0.30,0.80,Form("Incident #pi^{0} Eneregy = %.1f GeV",IncEnergy));	
	lat->DrawLatex(0.30,0.75,"Standard PHENIX Pb Shashlik EMCAL");	
	c->SaveAs(Form("PiZeroPlotsGen/%d/TotalNEDistribution.png",EnergyBin));


	ClusEDev->Draw();
	lat->DrawLatex(0.30,0.80,Form("Incident #pi^{0} Eneregy = %.1f GeV",IncEnergy));	
	lat->DrawLatex(0.30,0.75,"Standard PHENIX Pb Shashlik EMCAL");	
	c->SaveAs(Form("PiZeroPlotsGen/%d/ClusEDev.png",EnergyBin));

}
