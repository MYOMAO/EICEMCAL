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




void PlotEDep(int Energy, int BeamOpt){
	
	TString BeamName;
	
	if(BeamOpt == 0) BeamName = "e";
	if(BeamOpt == 1) BeamName = "mu";

	TString BeamNameFull;
	if(BeamOpt == 0) BeamNameFull = "Electron";
	if(BeamOpt == 1) BeamNameFull = "Muon";


	if(BeamOpt == 0){
		
		TotalEMin = 0.99;
		TotalEMax = 1.01;

	}


	if(BeamOpt == 1){
		
		TotalEMin = 0.89;
		TotalEMax = 0.91;

	}


	const int NEvents = 1000;

	double MinZ = 999;

	TLatex *lat = new TLatex();
	lat->SetNDC();
	lat->SetTextSize(0.035);


	gStyle->SetOptStat(0);

	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();

	double TotalBHEnergy = 0;
	double TotalFEMCALEnergy = 0;
	double TotalBHEnergyBF = 0;
	double Factor = 4;


	int NBins = 200;
	int LeakHisMin = 0;


	double Total;
	double Cumulative;
	double Percentage;
	double PercentDepth;
	double Percent20X0;
	int XBinLine;

	
	double WCharge = 74;
	double MoliereRadius = 0.0265 * (WCharge + 1.2);

	double value;
	double width;
	double nomalized;

	TString infile = Form("ShowerInfo_%d_%s.root",Energy,BeamName.Data());
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
	ShowerInfo->SetBranchAddress("layer",&layer);

	ShowerInfo->SetBranchAddress("FEMCALedep",&FEMCALedep);

	TH1D * BHedepHis = new TH1D("BHedepHis","",200,0,1000);
	BHedepHis->GetXaxis()->SetTitle("Total Back Leakage Energy Per Event (MeV)");
	BHedepHis->GetYaxis()->SetTitle("Counts");
	BHedepHis->GetXaxis()->CenterTitle();
	BHedepHis->GetYaxis()->CenterTitle();
	BHedepHis->GetYaxis()->SetTitleOffset(1.4);


	TH1D * BHedepTotalHis = new TH1D("BHedepTotalHis","",200,0,1000);
	BHedepTotalHis->GetXaxis()->SetTitle("Total Back + Front Leakage Energy Per Event (MeV)");
	BHedepTotalHis->GetYaxis()->SetTitle("Counts");
	BHedepTotalHis->GetXaxis()->CenterTitle();
	BHedepTotalHis->GetYaxis()->CenterTitle();
	BHedepTotalHis->GetYaxis()->SetTitleOffset(1.4);

	int NEntries = ShowerInfo->GetEntries();
	int evtidPre = 0;



	
	TH1D * EvsZ = new TH1D("EvsZ","",100,0,50);
	EvsZ->GetXaxis()->SetTitle("z (mm)");
	EvsZ->GetYaxis()->SetTitle("Energy (GeV)");
	EvsZ->GetXaxis()->CenterTitle();
	EvsZ->GetYaxis()->CenterTitle();
	EvsZ->GetYaxis()->SetTitleOffset(1.4);

	
	TH1D * ETotalHis = new TH1D("ETotalHis","",100,TotalEMin,TotalEMax);
	ETotalHis->GetXaxis()->SetTitle(Form("(Scintillator + Absorber Energy + Back + Front Leakage)/Incident %s Energy",BeamNameFull.Data()));
	ETotalHis->GetYaxis()->SetTitle("Counts");
	ETotalHis->GetXaxis()->CenterTitle();
	ETotalHis->GetYaxis()->CenterTitle();
	ETotalHis->GetYaxis()->SetTitleOffset(1.4);
	ETotalHis->SetTitle("Energy Conservation Check");



	TH1D * EDepLongTungs1GeV = new TH1D("EDepLongTungs1GeV","",NPoints1GeV,EdepCurve1GeVX);
	EDepLongTungs1GeV->GetXaxis()->SetTitle("Calorimeter Depth (X0)");
	EDepLongTungs1GeV->GetYaxis()->SetTitle("Longitudinal Energy Profile (% of E_{inc})");
	EDepLongTungs1GeV->SetTitle("Incident e^{-} Longitudinal Energy Deposition");
	EDepLongTungs1GeV->GetXaxis()->CenterTitle();
	EDepLongTungs1GeV->GetYaxis()->CenterTitle();
	EDepLongTungs1GeV->GetYaxis()->SetTitleOffset(1.2);
	EDepLongTungs1GeV->SetLineColor(kRed);



	TH1D * EDepLongTungs1GeVCum = new TH1D("EDepLongTungs1GeVCum","",NPoints1GeV,EdepCurve1GeVX);
	EDepLongTungs1GeVCum->GetXaxis()->SetTitle("Calorimeter Depth (X0)");
	EDepLongTungs1GeVCum->GetYaxis()->SetTitle("Cumulative Longitudinal Energy Profile (% of E_{inc})");
	EDepLongTungs1GeVCum->SetTitle("Incident e^{-} Longitudinal Energy Deposition");
	EDepLongTungs1GeVCum->GetXaxis()->CenterTitle();
	EDepLongTungs1GeVCum->GetYaxis()->CenterTitle();
	EDepLongTungs1GeVCum->GetYaxis()->SetTitleOffset(1.2);
	EDepLongTungs1GeVCum->SetLineColor(kRed);



	TH1D * EDepLongTungs10GeV = new TH1D("EDepLongTungs10GeV","",NPoints10GeV,EdepCurve10GeVX);
	EDepLongTungs10GeV->GetXaxis()->SetTitle("Calorimeter Depth (X0)");
	EDepLongTungs10GeV->GetYaxis()->SetTitle("Longitudinal Energy Profile (% of E_{inc})");
	EDepLongTungs10GeV->SetTitle("Incident e^{-} Longitudinal Energy Deposition");
	EDepLongTungs10GeV->GetXaxis()->CenterTitle();
	EDepLongTungs10GeV->GetYaxis()->CenterTitle();
	EDepLongTungs10GeV->GetYaxis()->SetTitleOffset(1.2);
	EDepLongTungs10GeV->SetLineColor(kGreen);



	TH1D * EDepLongTungs10GeVCum = new TH1D("EDepLongTungs10GeVCum","",NPoints1GeV,EdepCurve1GeVX);
	EDepLongTungs10GeVCum->GetXaxis()->SetTitle("Calorimeter Depth (X0)");
	EDepLongTungs10GeVCum->GetYaxis()->SetTitle("Cumulative Longitudinal Energy Profile (% of E_{inc})");
	EDepLongTungs10GeVCum->SetTitle("Incident e^{-} Longitudinal Energy Deposition");
	EDepLongTungs10GeVCum->GetXaxis()->CenterTitle();
	EDepLongTungs10GeVCum->GetYaxis()->CenterTitle();
	EDepLongTungs10GeVCum->GetYaxis()->SetTitleOffset(1.2);
	EDepLongTungs10GeVCum->SetLineColor(kRed);


	TH1D * EDepLongW = new TH1D("EDepLongW","",NPoints1GeV,EdepCurve1GeVX);
	EDepLongW->GetXaxis()->SetTitle("Calorimeter Depth (X0)");
	EDepLongW->GetYaxis()->SetTitle("Longitudinal Energy Profile (% of E_{inc})");
	EDepLongW->SetTitle("Incident e^{-} Longitudinal Energy Deposition");
	EDepLongW->GetXaxis()->CenterTitle();
	EDepLongW->GetYaxis()->CenterTitle();
	EDepLongW->GetYaxis()->SetTitleOffset(1.2);
	EDepLongW->SetLineColor(kBlue);



	TH1D * EDepLongWCum = new TH1D("EDepLongWCum","",NPoints1GeV,EdepCurve1GeVX);
	EDepLongWCum->GetXaxis()->SetTitle("Calorimeter Depth (X0)");
	EDepLongWCum->GetYaxis()->SetTitle("Cumulative Longitudinal Energy Profile (% of E_{inc})");
	EDepLongWCum->SetTitle("Incident e^{-} Longitudinal Energy Deposition");
	EDepLongWCum->GetXaxis()->CenterTitle();
	EDepLongWCum->GetYaxis()->CenterTitle();
	EDepLongWCum->GetYaxis()->SetTitleOffset(1.2);
	EDepLongWCum->SetLineColor(kBlue);


	TH1D * EDepLongWFT = new TH1D("EDepLongWFT","",NPoints1GeV,EdepCurve1GeVX);
	EDepLongWFT->GetXaxis()->SetTitle("Calorimeter Depth (X0)");
	EDepLongWFT->GetYaxis()->SetTitle("Transverse Energy Profile (% of E_{inc})");
	EDepLongWFT->SetTitle("Incident e^{-} Longitudinal Energy Deposition");
	EDepLongWFT->GetXaxis()->CenterTitle();
	EDepLongWFT->GetYaxis()->CenterTitle();
	EDepLongWFT->GetYaxis()->SetTitleOffset(1.2);
	EDepLongWFT->SetLineColor(kGreen);



	//Transverse

	TH1D * EDepTransTungs1GeVCum = new TH1D("EDepTransTungs1GeVCum","",NPoints1GeVTrans,EdepCurve1GeVXTransCum);
	EDepTransTungs1GeVCum->GetXaxis()->SetTitle("Calorimeter Radius (X0)");
	EDepTransTungs1GeVCum->GetYaxis()->SetTitle("Cumulative Transverse Energy Profile (% of E_{inc})");
	EDepTransTungs1GeVCum->SetTitle("Incident e^{-} Transverse Energy Deposition");
	EDepTransTungs1GeVCum->GetXaxis()->CenterTitle();
	EDepTransTungs1GeVCum->GetYaxis()->CenterTitle();
	EDepTransTungs1GeVCum->GetYaxis()->SetTitleOffset(1.2);
	EDepTransTungs1GeVCum->SetLineColor(kRed);


	TH1D * EDepTransW = new TH1D("EDepTransW","",NPoints1GeVTrans,EdepCurve1GeVXTransCum);
	EDepTransW->GetXaxis()->SetTitle("Calorimeter Radius (X0)");
	EDepTransW->GetYaxis()->SetTitle("Transverse Energy Profile (% of E_{inc})");
	EDepTransW->SetTitle("Incident e^{-} Transverse Energy Deposition");
	EDepTransW->GetXaxis()->CenterTitle();
	EDepTransW->GetYaxis()->CenterTitle();
	EDepTransW->GetYaxis()->SetTitleOffset(1.2);
	EDepTransW->SetLineColor(kBlack);



	TH1D * EDepTransWCorr = new TH1D("EDepTransWCorr","",NPoints1GeVTrans,EdepCurve1GeVXTransCum);
	EDepTransWCorr->GetXaxis()->SetTitle("Calorimeter Radius (X0)");
	EDepTransWCorr->GetYaxis()->SetTitle("Transverse Energy Profile (% of E_{inc})");
	EDepTransWCorr->SetTitle("Incident e^{-} Transverse Energy Deposition");
	EDepTransWCorr->GetXaxis()->CenterTitle();
	EDepTransWCorr->GetYaxis()->CenterTitle();
	EDepTransWCorr->GetYaxis()->SetTitleOffset(1.2);
	EDepTransWCorr->SetLineColor(kBlue);



	TH1D * EDepTransWCorrFT = new TH1D("EDepTransWCorrFT","",NPoints1GeVTrans,EdepCurve1GeVXTransCum);
	EDepTransWCorrFT->GetXaxis()->SetTitle("Calorimeter Radius (X0)");
	EDepTransWCorrFT->GetYaxis()->SetTitle("Transverse Energy Profile (% of E_{inc})");
	EDepTransWCorrFT->SetTitle("Incident e^{-} Transverse Energy Deposition");
	EDepTransWCorrFT->GetXaxis()->CenterTitle();
	EDepTransWCorrFT->GetYaxis()->CenterTitle();
	EDepTransWCorrFT->GetYaxis()->SetTitleOffset(1.2);
	EDepTransWCorrFT->SetLineColor(kGreen);


	TH1D * EDepTransWCum = new TH1D("EDepTransWCum","",NPoints1GeVTrans,EdepCurve1GeVXTransCum);
	EDepTransWCum->GetXaxis()->SetTitle("Calorimeter Radius (X0)");
	EDepTransWCum->GetYaxis()->SetTitle("Cumulative Transverse Energy Profile (% of E_{inc})");
	EDepTransWCum->SetTitle("Incident e^{-} Transverse Energy Deposition");
	EDepTransWCum->GetXaxis()->CenterTitle();
	EDepTransWCum->GetYaxis()->CenterTitle();
	EDepTransWCum->GetYaxis()->SetTitleOffset(1.2);
	EDepTransWCum->SetLineColor(kBlack);
	//	EDepTransWCum->SetLineColor(kBlue);

	TH1D * EDepTransWCumCorr = new TH1D("EDepTransWCumCorr","",NPoints1GeVTrans,EdepCurve1GeVXTransCum);
	EDepTransWCumCorr->GetXaxis()->SetTitle("Calorimeter Radius (X0)");
	EDepTransWCumCorr->GetYaxis()->SetTitle("Cumulative Transverse Energy Profile (% of E_{inc})");
	EDepTransWCumCorr->SetTitle("Incident e^{-} Transverse Energy Deposition");
	EDepTransWCumCorr->GetXaxis()->CenterTitle();
	EDepTransWCumCorr->GetYaxis()->CenterTitle();
	EDepTransWCumCorr->GetYaxis()->SetTitleOffset(1.2);
	EDepTransWCumCorr->SetLineColor(kBlue);

	TH1D * EDepTransWCumCorrFT = new TH1D("EDepTransWCumCorrFT","",NPoints1GeVTrans,EdepCurve1GeVXTransCum);
	EDepTransWCumCorrFT->GetXaxis()->SetTitle("Calorimeter Radius (X0)");
	EDepTransWCumCorrFT->GetYaxis()->SetTitle("Cumulative Transverse Energy Profile (% of E_{inc})");
	EDepTransWCumCorrFT->SetTitle("Incident e^{-} Transverse Energy Deposition");
	EDepTransWCumCorrFT->GetXaxis()->CenterTitle();
	EDepTransWCumCorrFT->GetYaxis()->CenterTitle();
	EDepTransWCumCorrFT->GetYaxis()->SetTitleOffset(1.2);
	EDepTransWCumCorrFT->SetLineColor(kGreen);


	TH2D * EDepXY = new TH2D("EDepXY","",50,0,100,50,0,100);
	EDepXY->GetXaxis()->SetTitle("X (cm)");
	EDepXY->GetYaxis()->SetTitle("Y (cm)");
	EDepXY->SetTitle("Incident e^{-} Energy Deposition vs XY (GeV)");
	EDepXY->GetXaxis()->CenterTitle();
	EDepXY->GetYaxis()->CenterTitle();
	EDepXY->GetYaxis()->SetTitleOffset(1.2);

	TH2D * EDepXYCorr = new TH2D("EDepXYCorr","",50,-5,5,50,-5,5);
	EDepXYCorr->GetXaxis()->SetTitle("Corrected X (cm)");
	EDepXYCorr->GetYaxis()->SetTitle("Corrected Y (cm)");
	EDepXYCorr->SetTitle("Incident e^{-} Energy Deposition vs XY (GeV)");
	EDepXYCorr->GetXaxis()->CenterTitle();
	EDepXYCorr->GetYaxis()->CenterTitle();
	EDepXYCorr->GetYaxis()->SetTitleOffset(1.2);



	TH2D * EDepXYCorrFT = new TH2D("EDepXYCorrFT","",50,-5,5,50,-5,5);
	EDepXYCorrFT->GetXaxis()->SetTitle("Corrected + Fine Toned X (cm)");
	EDepXYCorrFT->GetYaxis()->SetTitle("Corrected + Fine Toned Y (cm)");
	EDepXYCorrFT->SetTitle("Incident e^{-} Energy Deposition vs XY (GeV)");
	EDepXYCorrFT->GetXaxis()->CenterTitle();
	EDepXYCorrFT->GetYaxis()->CenterTitle();
	EDepXYCorrFT->GetYaxis()->SetTitleOffset(1.2);


	double ZPos;
	double ZPosFT;


	double ZBin;

	double XPos;
	double YPos;
	double RPos;
	double XBin;
	double YBin;

	double XPosCorr;
	double YPosCorr;
	double RPosCorr;




	double XPosCorrFT;
	double YPosCorrFT;
	double RPosCorrFT;


	double RPosRECO;
	scalefactor = (WRadLength + ScinThick/SamFrac);
	
	//New Parameters
//	scalefactor = (WRadLength + ScinThick);
//	MinDepthZ = 305;

	
	/*
	//Setup RECO Shower Position//
	float x;
	float y;


	double XPosRECO[NEvents];
	double YPosRECO[NEvents];



	TString infile2 = "G4EICDetector.root_g4femc_eval.root";
	TFile * fin2 = new TFile(infile2.Data());
	TTree * ntp_gshower = (TTree * ) fin2->Get("ntp_gshower");

	ntp_gshower->SetBranchAddress("x",&x);
	ntp_gshower->SetBranchAddress("y",&y);


	for(int i = 0; i <NEvents; i++){

	ntp_gshower->GetEntry(i);

	XPosRECO[i] = x;
	YPosRECO[i] = y;


	}



	//DONE
	*/

	MinDepthXY = 2 * MinDepthZ/(TMath::Exp(eta) - TMath::Exp(-1 * eta)); 
	//MinDepthXYFT = 2 * MinDepthZFT/(TMath::Exp(eta) - TMath::Exp(-1 * eta)); 

	cout << "MinDepthZ = " << MinDepthZ << "   MinDepthXY = " << MinDepthXY << endl;

	//	double WRadLength = 0.56;

	double VtxX = MinDepthX;
	double VtxY = MinDepthY;
	double VtxZ = MinDepthZ;

	TVector3 Vtx(VtxX,VtxY,VtxZ);
	TVector3 Pos;
	TVector3 RelPos;
	TVector3 RelLong;
	TVector3 RelTrans;


	for(int i = 0; i < NPoints1GeV; i++){

		EDepLongTungs1GeV->SetBinContent(i+1,EdepCurve1GeVY[i+1]);

	}

	for(int i = 0; i < NPoints1GeVTrans+1; i++){

		EDepTransTungs1GeVCum->SetBinContent(i+1,EdepCurve1GeVYTransCum[i+1]);

	}



	double TotalEnergyCounted = 0;
	double TotalEnergyCountedCum = 0;


	for(int i = 0; i < NEntries; i++){

		ShowerInfo->GetEntry(i);


		ZPos = ((Z1 + Z2)/2 - MinDepthZ )/scalefactor;
		//ZPosFT = (Z1 + Z2 - MinDepthZFT * 2);



		XPos = (X1 + X2)/2;
		YPos = (Y1 + Y2)/2;
		if(ZPos < (Z1 + Z2)/2 && layer == -1)  MinZ = (Z1 + Z2)/2;

		Pos.SetXYZ(XPos,YPos,(Z1+Z2)/2);
		
	    RelPos = Pos - Vtx;

		RPos = (sqrt(X1 * X1 + Y1 * Y1) - MinDepthXY)/scalefactor; 

		//XPosCorr = XPos - MinDepthXY * cos(phi);
		//YPosCorr = YPos - MinDepthXY * sin(phi);
		XPosCorr = XPos - MinDepthX;
		YPosCorr = YPos - MinDepthY;

		RPosCorr = sqrt(XPosCorr * XPosCorr + YPosCorr * YPosCorr)/scalefactor;



//		XPosCorrFT = XPos - MinDepthXYFT * cos(phi);
//		YPosCorrFT = YPos - MinDepthXYFT * sin(phi);
//		RPosCorrFT = sqrt(XPosCorrFT * XPosCorrFT + YPosCorrFT * YPosCorrFT);

		ZPosFT = RelPos.Dot(Vtx)/Vtx.Mag() * 2;
		RPosCorrFT =  sqrt(RelPos.Mag() * RelPos.Mag() - RelPos.Dot(Vtx)/Vtx.Mag() * RelPos.Dot(Vtx)/Vtx.Mag()); 
	
		RelLong = Vtx * (1/(Vtx.Mag() * RelPos.Dot(Vtx)/Vtx.Mag()) ); 
		RelTrans = RelPos - RelLong;
		
		XPosCorrFT = RelTrans.X();
		YPosCorrFT = RelTrans.Y();
		

		//	RPosRECO = sqrt( (XPos -  XPosRECO[evtid]) * (XPos -  XPosRECO[evtid]) + (YPos -  YPosRECO[evtid]) * (YPos -  YPosRECO[evtid]));

		//cout << "ZPos = " << ZPos << endl;

		if(layer == -1){
			EvsZ->Fill(ZPos,FEMCALedep);
		
			EDepLongW->Fill(ZPos,FEMCALedep * 1000);
			EDepLongWFT->Fill(ZPosFT,FEMCALedep * 1000);
		
			EDepXY->Fill(XPos,YPos,FEMCALedep);
			EDepXYCorr->Fill(XPosCorr,YPosCorr,FEMCALedep);
			EDepXYCorrFT->Fill(XPosCorrFT,YPosCorrFT,FEMCALedep);

			EDepTransW->Fill(RPos,FEMCALedep * 1000);
			EDepTransWCorr->Fill(RPosCorr,FEMCALedep * 1000);
			EDepTransWCorrFT->Fill(RPosCorrFT,FEMCALedep * 1000);
		}
		//cout << "RPos = " << RPos << endl;

		//cout << "BHedep = " << BHBLedep << endl;

		if(evtid != evtidPre){

			//cout << "Been here: " << "evtid = " << evtid << "    evtidPre = " << evtidPre << "  TotalEnergyCounted = " << TotalEnergyCounted << endl;

			BHedepHis->Fill(TotalBHEnergy);
			BHedepTotalHis->Fill(TotalBHEnergyBF);
			ETotalHis->Fill(TotalEnergyCounted/Energy);
		

			TotalBHEnergy = 0;
			TotalFEMCALEnergy = 0;
			TotalBHEnergyBF = 0;
			TotalEnergyCounted = 0;
		}


		TotalBHEnergy = BHBLedep * 1000 + TotalBHEnergy;
		TotalFEMCALEnergy = FEMCALedep * 1000 + TotalFEMCALEnergy;
		TotalBHEnergyBF =  BHBLedep * 1000 + BHFLedep * 1000 + TotalBHEnergyBF;

		evtidPre = evtid;

		TotalEnergyCounted = FEMCALedep + BHBLedep + BHFLedep + TotalEnergyCounted;
		TotalEnergyCountedCum = FEMCALedep + BHBLedep + BHFLedep + TotalEnergyCountedCum;
	}



	


	BHedepHis->Draw();

	c->SaveAs("ELeakAna/BackLeakage.png");


	TCanvas *c2 = new TCanvas("c2","c2",600,600);
	c2->SetLogy();
	c2->cd();

	BHedepTotalHis->Draw();

	c2->SaveAs("ELeakAna/TotalLeakage.png");

	c->cd();


	cout << "Total Percentage Sum = " << EDepLongTungs1GeV->Integral() << "    Total Percentage Width =  " << EDepLongTungs1GeV->Integral("width") << endl;

	for(int i = 0; i < EDepLongW->GetNbinsX(); i++){

		value = EDepLongW->GetBinContent(i+1);
		width = EDepLongW->GetBinWidth(i+1);
		nomalized = value/width;

		EDepLongW->SetBinContent(i+1, nomalized);

	}
	EDepLongW->Scale(100/EDepLongW->Integral("width"));


	for(int i = 0; i < EDepLongWFT->GetNbinsX(); i++){

		value = EDepLongWFT->GetBinContent(i+1);
		width = EDepLongWFT->GetBinWidth(i+1);
		nomalized = value/width;

		EDepLongWFT->SetBinContent(i+1, nomalized);

	}
	EDepLongWFT->Scale(100/EDepLongWFT->Integral("width"));


	//EDepLongW->Scale(EDepLongTungs1GeV->Integral()/EDepLongW->Integral());
	EDepLongW->Draw("hist");



	//Longitudinal 

	Total = EDepLongW->Integral();

	for(int i = 0; i < NPoints1GeV + 1;i++){
	
		

		Cumulative = EDepLongW->Integral(1,i+1);
		Percentage = Cumulative/Total;
		EDepLongWCum->SetBinContent(i+1,Percentage);


	}

	Total = EDepLongTungs1GeV->Integral();

	for(int i = 0; i < NPoints1GeV + 1; i++){

		Cumulative = EDepLongTungs1GeV->Integral(1,i+1);
		Percentage = Cumulative/Total;
		EDepLongTungs1GeVCum->SetBinContent(i+1,Percentage);

	}



	TLegend * leg  = new TLegend(0.12,0.35,0.55,0.55);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.037);
	leg->SetTextFont(42);
	leg->SetFillStyle(0);
	leg->SetLineWidth(3);





	TLegend * legTrans  = new TLegend(0.12,0.35,0.55,0.55);
	legTrans->SetBorderSize(0);
	legTrans->SetTextSize(0.037);
	legTrans->SetTextFont(42);
	legTrans->SetFillStyle(0);
	legTrans->SetLineWidth(3);


	TLegend * legLong  = new TLegend(0.12,0.35,0.55,0.55);
	legLong->SetBorderSize(0);
	legLong->SetTextSize(0.037);
	legLong->SetTextFont(42);
	legLong->SetFillStyle(0);
	legLong->SetLineWidth(3);




	Total = EDepTransW->Integral();

	


	for(int i = 0; i < NPoints1GeVTrans+1; i++){

		Cumulative = EDepTransW->Integral(1,i+1);

		Percentage = Cumulative/(NEvents * 1000);
		EDepTransWCum->SetBinContent(i+1,Percentage);

	}
	
	Total = EDepTransWCorr->Integral();



	for(int i = 0; i < NPoints1GeVTrans+1; i++){

		Cumulative = EDepTransWCorr->Integral(1,i+1);
		//Percentage = Cumulative/(NEvents * 1000);
		Percentage = Cumulative/Total;
		EDepTransWCumCorr->SetBinContent(i+1,Percentage);

	}
	cout << " ----------------------------------------------------------------------------------"  << endl;

	//	cout << "Total Cumulative Fraction = " << TotalEnergyCountedCum/NEvents << endl;

	cout << "MinZ = " << MinZ << endl;
	cout << " ---------------------------------------------------------------------------------- "  << endl;


	for(int i = 0; i < NPoints1GeVTrans+1; i++){

		Cumulative = EDepTransWCorrFT->Integral(1,i+1);
		Percentage = Cumulative/Total;
		EDepTransWCumCorrFT->SetBinContent(i+1,Percentage);

	}





	TLine *line90Percent = new TLine(0,90,18,90);
	line90Percent->SetLineStyle(2);
	line90Percent->SetLineWidth(2);
	line90Percent->SetLineColor(1);

	TLine *MoliereR = new TLine(MoliereRadius,0,MoliereRadius,100);
	MoliereR->SetLineStyle(2);
	MoliereR->SetLineWidth(2);
	MoliereR->SetLineColor(1);

	if(Energy == 1){
		lat->DrawLatex(0.30,0.80,Form("Incident %s Eneregy = 1 GeV",BeamNameFull.Data()));	

		EDepLongTungs1GeV->Draw("histSAME");
	//	EDepLongWFT->Draw("histSAME");
		leg->AddEntry(EDepLongTungs1GeV,"W Plates (John)","l");
		leg->AddEntry(EDepLongW,"W-Shashlik (Mine, Position Corrected)","l");	
//		leg->AddEntry(EDepLongWFT,"W-Shashlik (Mine, Position Corrected + Angle Rotated)","l");	
	
		leg->Draw("SAME");
		c->SaveAs(Form("NewStudy/%s/Longitudinal/EDepthLongStudy1GeV.png",BeamName.Data()));


		EDepLongWCum->Scale(100);
		EDepLongTungs1GeVCum->Scale(100);


		//EDepTransWCum->Draw("hist");
		EDepLongTungs1GeVCum->SetMinimum(0);
		EDepLongTungs1GeVCum->Draw("hist");
		EDepLongWCum->Draw("histSAME");
//		EDepTransWCumCorrFT->Draw("histSAME");	

		lat->DrawLatex(0.30,0.80,Form("Incident %s Eneregy = 1 GeV",BeamName.Data()));	
		lat->DrawLatex(0.30,0.75,"Sampling Thickness: 1/4 X0");	
		lat->DrawLatex(0.30,0.70,"Total RadLength: 20 X0");
		lat->DrawLatex(0.30,0.65,"Scintillator: 1.4 mm Polystyrene");

		legLong->AddEntry(EDepLongTungs1GeVCum,"W Plates (John)","l");
//		legTrans->AddEntry(EDepTransWCum,"W-Shashlik (Mine, No Position Correction)","l");
		legLong->AddEntry(EDepLongWCum,"W-Shashlik (Mine, Position Corrected)","l");
//		legTrans->AddEntry(EDepTransWCumCorrFT,"W-Shashlik (Mine, Position Corrected + Angle Rotated)","l");

		legLong->Draw("SAME");

		c->SaveAs(Form("NewStudy/%s/Longitudinal/EDepthLingStudy1GeVCum.png",BeamName.Data()));


		EvsZ->Draw();
		c->SaveAs(Form("NewStudy/%s/Debug/EvsZ.png",BeamName.Data()));


		//Transverse Profile

		EDepTransWCum->Scale(100);
		EDepTransWCumCorr->Scale(100);
		EDepTransWCumCorrFT->Scale(100);

		//EDepTransWCum->Draw("hist");
		EDepTransTungs1GeVCum->SetMinimum(0);
		EDepTransTungs1GeVCum->Draw("hist");
		EDepTransWCumCorr->Draw("histSAME");
//		EDepTransWCumCorrFT->Draw("histSAME");	
	
		lat->DrawLatex(0.30,0.80,"Incident %s Eneregy = 1 GeV");	
		lat->DrawLatex(0.30,0.75,"Sampling Thickness: 1/4 X0");	
		lat->DrawLatex(0.30,0.70,"Total RadLength: 20 X0");
		lat->DrawLatex(0.30,0.65,"Scintillator: 1.4 mm Polystyrene");

		legTrans->AddEntry(EDepTransTungs1GeVCum,"W Plates (John)","l");
//		legTrans->AddEntry(EDepTransWCum,"W-Shashlik (Mine, No Position Correction)","l");
		legTrans->AddEntry(EDepTransWCumCorr,"W-Shashlik (Mine, Position Corrected)","l");
//		legTrans->AddEntry(EDepTransWCumCorrFT,"W-Shashlik (Mine, Position Corrected + Angle Rotated)","l");

		legTrans->Draw("SAME");

		line90Percent->Draw("SAME");
		MoliereR->Draw("SAME");

		c->SaveAs(Form("NewStudy/%s/Transverse/EDepthTransStudy1GeV.png",BeamName.Data()));

		EDepTransWCumCorr->Draw("hist");
		c->SaveAs(Form("NewStudy/%s/Transverse/EDepTransWCumCorr.png",BeamName.Data()));


		EDepXY->Draw("COLZ");

		lat->DrawLatex(0.30,0.80,"Incident %s Eneregy = 1 GeV");	
		lat->DrawLatex(0.30,0.75,"Sampling Thickness: 1/4 X0");	
		lat->DrawLatex(0.30,0.70,"Total RadLength: 20 X0");
		lat->DrawLatex(0.30,0.65,"Scintillator: 1.4 mm Polystyrene");

		c->SaveAs(Form("NewStudy/%s/Transverse/EDepXY.png",BeamName.Data()));



		TCanvas * c3 = new TCanvas("c3","c3",600,600);
		c3->cd();
		c3->SetLogz();

		EDepXYCorr->Draw("COLZ");
		lat->DrawLatex(0.30,0.80,"Incident %s Eneregy = 1 GeV");	
		lat->DrawLatex(0.30,0.75,"Sampling Thickness: 1/4 X0");	
		lat->DrawLatex(0.30,0.70,"Total RadLength: 20 X0");
		lat->DrawLatex(0.30,0.65,"Scintillator: 1.4 mm Polystyrene");
		c3->SaveAs(Form("NewStudy/%s/Transverse/EDepXYCorr.png",BeamName.Data()));


	//	EDepXYCorrFT->Draw("COLZ");
		lat->DrawLatex(0.30,0.80,"Incident %s Eneregy = 1 GeV");	
		lat->DrawLatex(0.30,0.75,"Sampling Thickness: 1/4 X0");	
		lat->DrawLatex(0.30,0.70,"Total RadLength: 20 X0");
		lat->DrawLatex(0.30,0.65,"Scintillator: 1.4 mm Polystyrene");
		c3->SaveAs(Form("NewStudy/%s/Transverse/EDepXYCorrFT.png",BeamName.Data()));

		c->cd();


	}

	ETotalHis->Draw();
	c->SaveAs(Form("NewStudy/%s/Debug/EnergyConservationCheck.png",BeamName.Data()));



	if(Energy == 10){
		EDepLongTungs10GeV->Draw("histSAME");
		lat->DrawLatex(0.40,0.80,Form("Incident %s Eneregy = 10 GeV",BeamNameFull.Data()));		
		leg->AddEntry(EDepLongTungs1GeV,"W Plates (John)","l");
		leg->AddEntry(EDepLongW,"W-Shashlik (Mine)","l");	
		leg->Draw("SAME");
		c->SaveAs(Form("NewStudy/%s/Comparison/EDepthLongStudy10GeV.png",BeamName.Data()));
	}



	Total = EDepLongTungs1GeV->Integral();

	for(int i = 0; i < NPoints1GeV; i++){

		Cumulative = EDepLongTungs1GeV->Integral(1,i+1);

		Percentage = Cumulative/Total;


		EDepLongTungs1GeVCum->SetBinContent(i+1,Percentage);

	}




	EDepLongTungs1GeVCum->SetMinimum(0);
	EDepLongTungs1GeVCum->SetMaximum(1.0);

	EDepLongTungs1GeVCum->Draw("hist");

	XBinLine = EDepLongTungs1GeVCum->GetXaxis()->FindBin(4.0);
	PercentDepth = EDepLongTungs1GeVCum->GetBinContent(XBinLine) *100;

	lat->DrawLatex(0.20,0.75,Form("Inciden %s Beam Energy = 1 GeV",BeamNameFull.Data()));
	lat->DrawLatex(0.20,0.70,Form("4X0 Cumulative Energy Deposition = %.1f %%",PercentDepth));

	XBinLine = EDepLongTungs1GeVCum->GetXaxis()->FindBin(20.0);
	PercentDepth = EDepLongTungs1GeVCum->GetBinContent(XBinLine) * 100;

	lat->DrawLatex(0.20,0.65,Form("20X0 Cumulative Energy Deposition = %.1f %%",PercentDepth));


	TLine *line4X0 = new TLine(4,0,4,1);
	line4X0->SetLineStyle(2);
	line4X0->SetLineWidth(2);
	line4X0->SetLineColor(4);
	line4X0->Draw("SAME");

	TLine *line20X0 = new TLine(20,0,20,1);
	line20X0->SetLineStyle(2);
	line20X0->SetLineWidth(2);
	line20X0->SetLineColor(4);
	line20X0->Draw("SAME");



	c->SaveAs("ELeakAna/Cumulative/Cumulative1GeV.png");


	Total = EDepLongTungs10GeV->Integral();


	for(int i = 0; i < NPoints10GeV; i++){

		Cumulative = EDepLongTungs10GeV->Integral(1,i+1);

		Percentage = Cumulative/Total;


		EDepLongTungs10GeVCum->SetBinContent(i+1,Percentage);

	}


	EDepLongTungs10GeVCum->Draw("hist");


	XBinLine = EDepLongTungs10GeVCum->GetXaxis()->FindBin(4.0);
	PercentDepth = EDepLongTungs10GeVCum->GetBinContent(XBinLine) * 100;

	lat->DrawLatex(0.20,0.75,Form("Inciden %s Beam Energy = 10 GeV" ,BeamNameFull.Data()));
	lat->DrawLatex(0.20,0.70,Form("4X0 Cumulative Energy Deposition = %.1f %%",PercentDepth ));

	XBinLine = EDepLongTungs10GeVCum->GetXaxis()->FindBin(20.0);
	PercentDepth = EDepLongTungs10GeVCum->GetBinContent(XBinLine) * 100;

	lat->DrawLatex(0.20,0.65,Form("20X0 Cumulative Energy Deposition = %.1f %%",PercentDepth));


	line4X0->Draw("SAME");
	line20X0->Draw("SAME");

	c->SaveAs("ELeakAna/Cumulative/Cumulative10GeV.png");


}
