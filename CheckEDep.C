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
#include "CheckEDep.h"

using namespace std;

using std::cout;
using std::endl;


void CheckEDep(int Energy, int BeamOpt){
	
	TString BeamName;
	
	float NominalEnergy;

	if(BeamOpt == 0) BeamName = "e";
	if(BeamOpt == 1) BeamName = "mu";

	TString BeamNameFull;
	if(BeamOpt == 0) BeamNameFull = "Electron";
	if(BeamOpt == 1) BeamNameFull = "Muon";


	if(Energy == 0){
		MeanAbs = 157.79;
		MeanScin = 23.0554;
	}


	if(Energy == 1){
		MeanAbs = 172.795;
		MeanScin = 24.313856;
	}

	if(BeamOpt == 0){
		
		TotalEMin = 0.99;
		TotalEMax = 1.01;

	}

	if(BeamOpt == 0){
		
		RatioUp = 0.995;
		RatioDown = 1.001;

	}

	if(BeamOpt == 1){
		
		TotalEMin = 0.89;
		TotalEMax = 0.91;

	}

	if(BeamOpt == 1){
		
		RatioUp = 0.892;
		RatioDown = 0.896;

	}


	if(Energy == 0){
		
		NominalEnergy = 0.4;

	}
	if(Energy == 1){
		
		NominalEnergy = 1.0;

	}

	const int NEvents = 100;



	TLatex *lat = new TLatex();
	lat->SetNDC();
	lat->SetTextSize(0.035);


	gStyle->SetOptStat(0);

	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();


	double TotalBackLeak = 0;
	double TotalFrontLeak = 0;
	
	double TotalScinEnergy = 0;
	double TotalAbsoEnergy = 0;
	double TotalAirEnergy = 0;
	
	double TotalEnergy = 0;




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
	ShowerInfo->SetBranchAddress("FEMCALLYDep",&FEMCALLYDep);
	

	TH1D * AbsEnergyDis = new TH1D("AbsEnergyDis","",NHisBins,MeanAbs-WidthAbsDown, MeanAbs+WidthAbsUp);
	AbsEnergyDis->GetXaxis()->SetTitle(Form("Total Incident %s Energy Deposited to the W Absorber (MeV)",BeamNameFull.Data()));
	AbsEnergyDis->GetYaxis()->SetTitle("Number of Events");
	AbsEnergyDis->GetXaxis()->CenterTitle();
	AbsEnergyDis->GetYaxis()->CenterTitle();
	AbsEnergyDis->GetYaxis()->SetTitleOffset(1.2);


	TH1D * ScinEnergyDis = new TH1D("ScinEnergyDis","",NHisBins,MeanScin-WidthScinDown,MeanScin+WidthScinUp);
	ScinEnergyDis->GetXaxis()->SetTitle(Form("Total Incident %s Energy Deposited to the Polystyrene Scintillator (MeV)",BeamNameFull.Data()));
	ScinEnergyDis->GetYaxis()->SetTitle("Number of Events");
	ScinEnergyDis->GetXaxis()->CenterTitle();
	ScinEnergyDis->GetYaxis()->CenterTitle();
	ScinEnergyDis->GetYaxis()->SetTitleOffset(1.2);




	TH1D * ForwardLeakDis = new TH1D("ForwardLeakDis","",NHisBins,MinFLeak,MaxFLeak);
	ForwardLeakDis->GetXaxis()->SetTitle(Form("Total Incident %s Forward Energy Leakage (MeV)",BeamNameFull.Data()));
	ForwardLeakDis->GetYaxis()->SetTitle("Number of Events");
	ForwardLeakDis->GetXaxis()->CenterTitle();
	ForwardLeakDis->GetYaxis()->CenterTitle();
	ForwardLeakDis->GetYaxis()->SetTitleOffset(1.2);




	TH1D * BackLeakDis = new TH1D("BackLeakDis","",NHisBins,MinBLeak,MaxBLeak);
	BackLeakDis->GetXaxis()->SetTitle(Form("Total Incident %s Backward Energy Leakage (MeV)",BeamNameFull.Data()));
	BackLeakDis->GetYaxis()->SetTitle("Number of Events");
	BackLeakDis->GetXaxis()->CenterTitle();
	BackLeakDis->GetYaxis()->CenterTitle();
	BackLeakDis->GetYaxis()->SetTitleOffset(1.2);


	TH1D * AirDis = new TH1D("AirDis","",NHisBins,MinAir,MaxAir);
	AirDis->GetXaxis()->SetTitle(Form("Total %s Energy Deposited to the Air (MeV)",BeamNameFull.Data()));
	AirDis->GetYaxis()->SetTitle("Number of Events");
	AirDis->GetXaxis()->CenterTitle();
	AirDis->GetYaxis()->CenterTitle();
	AirDis->GetYaxis()->SetTitleOffset(1.2);





	TH1D * TotalEnergyHis = new TH1D("TotalEnergyHis","",NHisBins,Energy * 1000 * RatioUp,Energy * 1000 * RatioDown);
	TotalEnergyHis->GetXaxis()->SetTitle(Form("Total Sum of %s Energy Per Event (MeV)",BeamNameFull.Data()));
	TotalEnergyHis->GetYaxis()->SetTitle("Number of Events");
	TotalEnergyHis->GetXaxis()->CenterTitle();
	TotalEnergyHis->GetYaxis()->CenterTitle();
	TotalEnergyHis->GetYaxis()->SetTitleOffset(1.2);



	int NEntries = ShowerInfo->GetEntries();
	int evtidPre = 0;

	for(int i = 0; i < NEntries; i++){

		ShowerInfo->GetEntry(i);


		if(evtid != evtidPre)
		{

			AbsEnergyDis->Fill(TotalAbsoEnergy);
			ScinEnergyDis->Fill(TotalScinEnergy);
			ForwardLeakDis->Fill(TotalFrontLeak);
			BackLeakDis->Fill(TotalBackLeak);
			TotalEnergyHis->Fill(TotalEnergy);
			AirDis->Fill(TotalAirEnergy);

			TotalBackLeak = 0;
			TotalFrontLeak = 0;
			TotalScinEnergy = 0;
			TotalAbsoEnergy = 0;
			TotalAirEnergy = 0;
			TotalEnergy = 0;

		}
	

		if(FEMCALLYDep == 2) TotalAbsoEnergy  = TotalAbsoEnergy + FEMCALedep * 1000;
		if(FEMCALLYDep == 3) TotalAirEnergy = TotalAirEnergy + FEMCALedep * 1000;
		if(FEMCALLYDep == 4) TotalScinEnergy  = TotalScinEnergy + FEMCALedep * 1000;

		TotalFrontLeak = TotalFrontLeak + BHFLedep * 1000;
		TotalBackLeak = TotalBackLeak + BHBLedep * 1000;

		TotalEnergy = TotalEnergy + FEMCALedep * 1000 + BHFLedep * 1000 + BHBLedep * 1000;
		
		evtidPre = evtid;
	}

	MinHisAbs = AbsEnergyDis->GetMinimum();
	MaxHisAbs = AbsEnergyDis->GetMaximum() * 1.4;
 	AbsEnergyDis->SetMaximum(MaxHisAbs);
	HisMean = AbsEnergyDis->GetMean();

	AbsEnergyDis->Draw();

	lat->DrawLatex(0.30,0.85,Form("Incident %s Eneregy = %.1f GeV",BeamName.Data(),NominalEnergy));	
	lat->DrawLatex(0.30,0.80,"Sampling Thickness: 1/4 X0");	
	lat->DrawLatex(0.30,0.75,"Total RadLength: 20 X0");
	lat->DrawLatex(0.30,0.70,"Scintillator: 1.4 mm Polystyrene");
	lat->DrawLatex(0.30,0.65,Form("Mean = %.1f MeV", HisMean));
	lat->DrawLatex(0.30,0.60,Form("Expected Energy Deposited = %.1f MeV",MeanAbs));



	TLine *AbsLine = new TLine(MeanAbs,MinHisAbs,MeanAbs,MaxHisAbs);
	AbsLine->SetLineStyle(2);
	AbsLine->SetLineWidth(2);
	AbsLine->SetLineColor(2);
	AbsLine->Draw("SAME");

	c->SaveAs(Form("CheckDepPlots/%s/AbsEnergyDis_%d.png",BeamName.Data(),Energy));


	MinHisScin = ScinEnergyDis->GetMinimum();
	MaxHisScin = ScinEnergyDis->GetMaximum() * 1.4;
 	ScinEnergyDis->SetMaximum(MaxHisScin);
	HisMean = ScinEnergyDis->GetMean();


	ScinEnergyDis->Draw();

	lat->DrawLatex(0.30,0.85,Form("Incident %s Eneregy = %.1f GeV",BeamName.Data(),NominalEnergy));	
	lat->DrawLatex(0.30,0.80,"Sampling Thickness: 1/4 X0");	
	lat->DrawLatex(0.30,0.75,"Total RadLength: 20 X0");
	lat->DrawLatex(0.30,0.70,"Scintillator: 1.4 mm Polystyrene");
	lat->DrawLatex(0.30,0.65,Form("Mean = %.1f MeV", HisMean));
	lat->DrawLatex(0.30,0.60,Form("Expected Energy Deposited = %.1f MeV",MeanScin));



	TLine *ScinLine = new TLine(MeanScin,MinHisScin,MeanScin,MaxHisScin);
	ScinLine->SetLineStyle(2);
	ScinLine->SetLineWidth(2);
	ScinLine->SetLineColor(2);
	ScinLine->Draw("SAME");

	c->SaveAs(Form("CheckDepPlots/%s/ScinEnergyDis_%d.png",BeamName.Data(),Energy));



	ForwardLeakDis->Draw();
	lat->DrawLatex(0.30,0.85,Form("Incident %s Eneregy = %.1f GeV",BeamName.Data(),NominalEnergy));	
	lat->DrawLatex(0.30,0.80,"Sampling Thickness: 1/4 X0");	
	lat->DrawLatex(0.30,0.75,"Total RadLength: 20 X0");
	lat->DrawLatex(0.30,0.70,"Scintillator: 1.4 mm Polystyrene");
	c->SaveAs(Form("CheckDepPlots/%s/ForwardLeakDis_%d.png",BeamName.Data(),Energy));



	BackLeakDis->Draw();
	lat->DrawLatex(0.30,0.85,Form("Incident %s Eneregy = %.1f GeV",BeamName.Data(),NominalEnergy));	
	lat->DrawLatex(0.30,0.80,"Sampling Thickness: 1/4 X0");	
	lat->DrawLatex(0.30,0.75,"Total RadLength: 20 X0");
	lat->DrawLatex(0.30,0.70,"Scintillator: 1.4 mm Polystyrene");
	c->SaveAs(Form("CheckDepPlots/%s/BackLeakDis_%d.png",BeamName.Data(),Energy));

	
	TotalEnergyHis->Draw();
	lat->DrawLatex(0.30,0.85,Form("Incident %s Eneregy = %.1f GeV",BeamName.Data(),NominalEnergy));	
	lat->DrawLatex(0.30,0.80,"Sampling Thickness: 1/4 X0");	
	lat->DrawLatex(0.30,0.75,"Total RadLength: 20 X0");
	lat->DrawLatex(0.30,0.70,"Scintillator: 1.4 mm Polystyrene");
	c->SaveAs(Form("CheckDepPlots/%s/TotalEnergyHis_%d.png",BeamName.Data(),Energy));


	AirDis->Draw();
	lat->DrawLatex(0.30,0.85,Form("Incident %s Eneregy = %.1f GeV",BeamName.Data(),NominalEnergy));	
	lat->DrawLatex(0.30,0.80,"Sampling Thickness: 1/4 X0");	
	lat->DrawLatex(0.30,0.75,"Total RadLength: 20 X0");
	lat->DrawLatex(0.30,0.70,"Scintillator: 1.4 mm Polystyrene");
	c->SaveAs(Form("CheckDepPlots/%s/AirEnergyHis_%d.png",BeamName.Data(),Energy));



}
