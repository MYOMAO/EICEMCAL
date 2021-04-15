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
#include "PlotEDepComp.h"

using namespace std;

using std::cout;
using std::endl;




void PlotEDepComp(int Energy){



	TLatex *lat = new TLatex();
	lat->SetNDC();
	lat->SetTextSize(0.035);


	gStyle->SetOptStat(0);

	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();

	double TotalBHEnergy = 0;
	double TotalFEMCALEnergy = 0;
	double TotalBHEnergyBF = 0;


	int NBins = 200;
	int LeakHisMin = 0;




	TString infile = Form("ShowerInfo_%d.root",Energy);
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







	TH1D * EDepLongTungs1GeV = new TH1D("EDepLongTungs1GeV","",NPoints1GeV,EdepCurve1GeVX);
	EDepLongTungs1GeV->GetXaxis()->SetTitle("Calorimeter Depth (X0)");
	EDepLongTungs1GeV->GetYaxis()->SetTitle("Longitudinal Energy Profile (% of E_{inc})");
	EDepLongTungs1GeV->SetTitle("Incident e^{-} Longitudinal Energy Deposition");
	EDepLongTungs1GeV->GetXaxis()->CenterTitle();
	EDepLongTungs1GeV->GetYaxis()->CenterTitle();
	EDepLongTungs1GeV->GetYaxis()->SetTitleOffset(1.2);
	EDepLongTungs1GeV->SetLineColor(kRed);

	TH1D * EDepLongTungs10GeV = new TH1D("EDepLongTungs1GeV","",NPoints10GeV,EdepCurve10GeVX);
	EDepLongTungs10GeV->GetXaxis()->SetTitle("Calorimeter Depth (X0)");
	EDepLongTungs10GeV->GetYaxis()->SetTitle("Longitudinal Energy Profile (% of E_{inc})");
	EDepLongTungs10GeV->SetTitle("Incident e^{-} Longitudinal Energy Deposition");
	EDepLongTungs10GeV->GetXaxis()->CenterTitle();
	EDepLongTungs10GeV->GetYaxis()->CenterTitle();
	EDepLongTungs10GeV->GetYaxis()->SetTitleOffset(1.2);
	EDepLongTungs10GeV->SetLineColor(kGreen);





	TH1D * EDepTotalW = new TH1D("EDepTotalW","",100,10,float(Energy)/2.5 * 100);
	EDepTotalW->GetXaxis()->SetTitle("Total Visible Energy (MeV)");
	EDepTotalW->GetYaxis()->SetTitle("Fraction");
	EDepTotalW->GetXaxis()->CenterTitle();
	EDepTotalW->GetYaxis()->CenterTitle();
	EDepTotalW->GetYaxis()->SetTitleOffset(1.4);
	EDepTotalW->SetLineColor(kBlack);





	TH1D * EDepTotalWRef = new TH1D("EDepTotalWRef","",100,10,float(Energy)/2.5 * 100);
	EDepTotalWRef->GetXaxis()->SetTitle("Total Visible Energy (MeV)");
	EDepTotalWRef->GetYaxis()->SetTitle("Fraction");
	EDepTotalWRef->GetXaxis()->CenterTitle();
	EDepTotalWRef->GetYaxis()->CenterTitle();
	EDepTotalWRef->GetYaxis()->SetTitleOffset(1.4);
	EDepTotalWRef->SetLineColor(kBlue);




	TH1D * EDepLongW = new TH1D("EDepLongW","",NPoints1GeV,EdepCurve1GeVX);
	EDepLongW->GetXaxis()->SetTitle("Calorimeter Depth (X0)");
	EDepLongW->GetYaxis()->SetTitle("Longitudinal Energy Profile (% of E_{inc})");
	EDepLongW->SetTitle("Incident e^{-} Longitudinal Energy Deposition");
	EDepLongW->GetXaxis()->CenterTitle();
	EDepLongW->GetYaxis()->CenterTitle();
	EDepLongW->GetYaxis()->SetTitleOffset(1.2);
	EDepLongW->SetLineColor(kBlack);




	double ZPos;
	double ZBin;

	double WRadLength = 0.56;

	for(int i = 0; i < EVis8GeV; i++){

		EDepTotalWRef->Fill(EVis8GeVX[i],EVis8GeVY[i]);

	}


	for(int i = 0; i < NPoints1GeV; i++){

		EDepLongTungs1GeV->SetBinContent(i+1,EdepCurve1GeVY[i+1]);

	}

	for(int i = 0; i < NPoints10GeV; i++){

		EDepLongTungs10GeV->SetBinContent(i+1,EdepCurve10GeVY[i+1]);

	}




	for(int i = 0; i < NEntries; i++){

		ShowerInfo->GetEntry(i);

		ZPos = (Z1 + Z2 - MinDepth * 2) * WRadLength/2;

		//cout << "ZPos = " << ZPos << endl;

		EDepLongW->Fill(ZPos,FEMCALedep * 100);

		//cout << "BHedep = " << BHBLedep << endl;

		if(evtid != evtidPre){


			BHedepHis->Fill(TotalBHEnergy);
			BHedepTotalHis->Fill(TotalBHEnergyBF);
			EDepTotalW->Fill(TotalFEMCALEnergy);


			TotalBHEnergy = 0;
			TotalFEMCALEnergy = 0;
			TotalBHEnergyBF = 0;

		}


		TotalBHEnergy = BHBLedep * 1000 + TotalBHEnergy;
		TotalFEMCALEnergy = FEMCALedep * 1000 + TotalFEMCALEnergy;
		TotalBHEnergyBF =  BHBLedep * 1000 + BHFLedep * 1000 + TotalBHEnergyBF;

		evtidPre = evtid;

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

	EDepLongW->Scale(EDepLongTungs1GeV->Integral()/EDepLongW->Integral());
	EDepLongW->Draw("hist");



	TLegend * leg  = new TLegend(0.40,0.55,0.75,0.75);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.040);
	leg->SetTextFont(42);
	leg->SetFillStyle(0);
	leg->SetLineWidth(3);


	if(Energy == 1){
		lat->DrawLatex(0.40,0.80,"Incident Electron Eneregy = 1 GeV");	
		EDepLongTungs1GeV->Draw("histSAME");
		leg->AddEntry(EDepLongTungs1GeV,"W Plates (John)","l");
		leg->AddEntry(EDepLongW,"W-Shashlik (Mine)","l");	
		leg->Draw("SAME");
		c->SaveAs("ELeakAna/Comparison/EDepthLongStudy1GeV.png");
	}




	if(Energy == 10){
		EDepLongTungs10GeV->Draw("histSAME");
		lat->DrawLatex(0.40,0.80,"Incident Electron Eneregy = 10 GeV");		
		leg->AddEntry(EDepLongTungs1GeV,"W Plates (John)","l");
		leg->AddEntry(EDepLongW,"W-Shashlik (Mine)","l");	
		leg->Draw("SAME");
		c->SaveAs("ELeakAna/Comparison/EDepthLongStudy10GeV.png");
	}



	if(Energy == 8){

		leg  = new TLegend(0.10,0.60,0.60,0.75);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.040);
		leg->SetTextFont(42);
		leg->SetFillStyle(0);
		leg->SetLineWidth(3);


		EDepTotalW->Scale(1.0/EDepTotalW->Integral());
		EDepTotalWRef->Scale(1.0/EDepTotalWRef->Integral());
		EDepTotalW->SetMaximum(0.22);
		EDepTotalW->Draw("hist");
		EDepTotalWRef->Draw("histSAME");
		lat->DrawLatex(0.12,0.85,"W Plate/Shashlik Calorimeter");		
		lat->DrawLatex(0.12,0.80,"20 Plates/Layer and 1X0 Width Layer Thickness (3.5mm)");		
		lat->DrawLatex(0.12,0.75,"Incident Electron Eneregy = 10 GeV");		
		leg->AddEntry(EDepTotalWRef,"Cheung-Haggerty","l");
		leg->AddEntry(EDepTotalW,"Zhaozhong (W : Scintillator Width = 3.5 : 1)","l");	
		leg->Draw("SAME");
		c->SaveAs("ELeakAna/Comparison/VisEnergy.png");
	}

}
