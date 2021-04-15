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




void PlotEDepCompLay(int Energy, int Option){

	TString Thickness;

	if(Option == 0) Thickness = "Quarter";
	if(Option == 1) Thickness = "Half";
	if(Option == 2) Thickness = "Unity";


	int NLayers;

	if(Option == 0) NLayers = 16;
	if(Option == 1) NLayers = 20;
	if(Option == 2) NLayers = 20;

	int NBins = 40;
	double EMin = 100;
	double EMax = 300;


	TLatex *lat = new TLatex();
	lat->SetNDC();
	lat->SetTextSize(0.035);


	gStyle->SetOptStat(0);

	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();

	double TotalBHEnergy = 0;
	double TotalFEMCALEnergy = 0;
	double TotalBHEnergyBF = 0;


	int LeakHisMin = 0;



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



	TH1D * EDepTotalWRef1 = new TH1D("EDepTotalWRef1","",NBins,EMin,EMax);
	EDepTotalWRef1->GetXaxis()->SetTitle("Total Visible Energy (MeV)");
	EDepTotalWRef1->GetYaxis()->SetTitle("Fraction");
	EDepTotalWRef1->GetXaxis()->CenterTitle();
	EDepTotalWRef1->GetYaxis()->CenterTitle();
	EDepTotalWRef1->GetYaxis()->SetTitleOffset(1.4);
	EDepTotalWRef1->SetLineColor(kBlue);
	EDepTotalWRef1->SetLineWidth(2);




	TH1D * EDepTotalWRef2 = new TH1D("EDepTotalWRef2","",NBins,EMin,EMax);
	EDepTotalWRef2->GetXaxis()->SetTitle("Total Visible Energy (MeV)");
	EDepTotalWRef2->GetYaxis()->SetTitle("Fraction");
	EDepTotalWRef2->GetXaxis()->CenterTitle();
	EDepTotalWRef2->GetYaxis()->CenterTitle();
	EDepTotalWRef2->GetYaxis()->SetTitleOffset(1.4);
	EDepTotalWRef2->SetLineColor(kBlue);
	EDepTotalWRef2->SetLineWidth(2);





	TH1D * EDepTotalWRef4 = new TH1D("EDepTotalWRef2","",NBins,EMin,EMax);
	EDepTotalWRef4->GetXaxis()->SetTitle("Total Visible Energy (MeV)");
	EDepTotalWRef4->GetYaxis()->SetTitle("Fraction");
	EDepTotalWRef4->GetXaxis()->CenterTitle();
	EDepTotalWRef4->GetYaxis()->CenterTitle();
	EDepTotalWRef4->GetYaxis()->SetTitleOffset(1.4);
	EDepTotalWRef4->SetLineColor(kBlue);
	EDepTotalWRef4->SetLineWidth(2);

	TH1D * EDepTotalWRef8 = new TH1D("EDepTotalWRef2","",NBins,EMin,EMax);
	EDepTotalWRef8->GetXaxis()->SetTitle("Total Visible Energy (MeV)");
	EDepTotalWRef8->GetYaxis()->SetTitle("Fraction");
	EDepTotalWRef8->GetXaxis()->CenterTitle();
	EDepTotalWRef8->GetYaxis()->CenterTitle();
	EDepTotalWRef8->GetYaxis()->SetTitleOffset(1.4);
	EDepTotalWRef8->SetLineColor(kBlue);
	EDepTotalWRef8->SetLineWidth(2);



	for(int i = 0; i < EVis1GeV; i++){

		EDepTotalWRef1->Fill(EVis1GeVX[i],EVis1GeVY[i]);

	}

	for(int i = 0; i < EVis2GeV; i++){

		EDepTotalWRef2->Fill(EVis2GeVX[i],EVis2GeVY[i]);

	}


	for(int i = 0; i < EVis4GeV; i++){

		EDepTotalWRef4->Fill(EVis4GeVX[i],EVis4GeVY[i]);

	}

	for(int i = 0; i < EVis8GeV; i++){

		EDepTotalWRef8->Fill(EVis8GeVX[i],EVis8GeVY[i]);

	}








	for(int i = 0; i < NPoints1GeV; i++){

		EDepLongTungs1GeV->SetBinContent(i+1,EdepCurve1GeVY[i+1]);

	}

	for(int i = 0; i < NPoints10GeV; i++){

		EDepLongTungs10GeV->SetBinContent(i+1,EdepCurve10GeVY[i+1]);

	}


	const int NFiles = 3;

	TString infile[NFiles]={ Form("CompLay/%s/%d/ShowerInfo_35to10.root",Thickness.Data(),Energy), Form("CompLay/%s/%d/ShowerInfo_35to14.root",Thickness.Data(),Energy), Form("CompLay/%s/%d/ShowerInfo_35to15.root",Thickness.Data(),Energy)};



	int color[NFiles] = {1,2,3};


	TString Length[NFiles] = {"1.0 mm" , "1.4 mm" , "1.5 mm"};


	TH1D * EDepTotalW[NFiles]; 


	for(int q = 0; q < NFiles;q++){

		TFile * fin = new TFile(infile[q].Data());
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

		int NEntries = ShowerInfo->GetEntries();
		int evtidPre = 0;






		EDepTotalW[q]	= new TH1D(Form("EDepTotalW_q"),"",NBins,EMin,EMax);
		EDepTotalW[q]->GetXaxis()->SetTitle("Total Visible Energy (MeV)");
		EDepTotalW[q]->GetYaxis()->SetTitle("Fraction");
		EDepTotalW[q]->GetXaxis()->CenterTitle();
		EDepTotalW[q]->GetYaxis()->CenterTitle();
		EDepTotalW[q]->GetYaxis()->SetTitleOffset(1.4);
		EDepTotalW[q]->SetLineColor(color[q]);
		EDepTotalW[q]->SetLineWidth(2);





		double ZPos;
		double ZBin;

		double WRadLength = 0.56;





		for(int i = 0; i < NEntries; i++){

			ShowerInfo->GetEntry(i);

			ZPos = (Z1 + Z2 - MinDepth * 2) * WRadLength/2;



			if(evtid != evtidPre){


				EDepTotalW[q]->Fill(TotalFEMCALEnergy);


				TotalBHEnergy = 0;
				TotalFEMCALEnergy = 0;
				TotalBHEnergyBF = 0;

			}


			TotalBHEnergy = BHBLedep * 1000 + TotalBHEnergy;
			TotalFEMCALEnergy = FEMCALedep * 1000 + TotalFEMCALEnergy;
			TotalBHEnergyBF =  BHBLedep * 1000 + BHFLedep * 1000 + TotalBHEnergyBF;

			evtidPre = evtid;

		}



	}


	/*

	   TLegend * leg  = new TLegend(0.40,0.42,0.75,0.74);
	   leg->SetBorderSize(0);
	   leg->SetTextSize(0.040);
	   leg->SetTextFont(42);
	   leg->SetFillStyle(0);
	   leg->SetLineWidth(3);

*/

	TLegend * leg  = new TLegend(0.10,0.485,0.60,0.735);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.040);
	leg->SetTextFont(42);
	leg->SetFillStyle(0);
	leg->SetLineWidth(3);


	if(Option == 2){


		lat->DrawLatex(0.12,0.85,"W Plate/Shashlik Calorimeter");		
		lat->DrawLatex(0.12,0.80,Form("%d Plates/Layer and 1X0 Width Layer Thickness (3.5mm)",NLayers));		
		lat->DrawLatex(0.12,0.75,Form("Incident Electron Eneregy = %d GeV",Energy));	


		if(Energy == 1){


			EDepTotalWRef1->Scale(1.0/EDepTotalWRef1->Integral());
			EDepTotalWRef1->SetMaximum(0.32);
			EDepTotalWRef1->Draw("hist");

			leg->AddEntry(EDepTotalWRef1,"Cheung-Haggerty","l");

			for(int q = 0; q < NFiles; q++){
				EDepTotalW[q]->Scale(1.0/EDepTotalW[q]->Integral());
				EDepTotalW[q]->Draw("histSAME");
				leg->AddEntry(EDepTotalW[q],Form("Zhaozhong (%s Polystyrene)",Length[q].Data()),"l");	
			}


			leg->Draw("SAME");
		}



		if(Energy == 2){



			EDepTotalWRef2->Scale(1.0/EDepTotalWRef2->Integral());
			EDepTotalWRef2->SetMaximum(0.32);
			EDepTotalWRef2->Draw("hist");
			lat->DrawLatex(0.12,0.85,"W Plate/Shashlik Calorimeter");		
			lat->DrawLatex(0.12,0.80,"20 Plates/Layer and 1X0 Width Layer Thickness (3.5mm)");		
			lat->DrawLatex(0.12,0.75,"Incident Electron Eneregy = 8 GeV");		
			leg->AddEntry(EDepTotalWRef2,"Cheung-Haggerty","l");

			for(int q = 0; q < NFiles; q++){
				EDepTotalW[q]->Scale(1.0/EDepTotalW[q]->Integral());
				EDepTotalW[q]->Draw("histSAME");
				leg->AddEntry(EDepTotalW[q],Form("Zhaozhong (%s Polystyrene)",Length[q].Data()),"l");	
			}


			leg->Draw("SAME");
		}



		if(Energy == 4){


			EDepTotalWRef4->Scale(1.0/EDepTotalWRef4->Integral());
			EDepTotalWRef4->SetMaximum(0.32);
			EDepTotalWRef4->Draw("hist");
			lat->DrawLatex(0.12,0.85,"W Plate/Shashlik Calorimeter");		
			lat->DrawLatex(0.12,0.80,"20 Plates/Layer and 1X0 Width Layer Thickness (3.5mm)");		
			lat->DrawLatex(0.12,0.75,"Incident Electron Eneregy = 8 GeV");		
			leg->AddEntry(EDepTotalWRef4,"Cheung-Haggerty","l");

			for(int q = 0; q < NFiles; q++){
				EDepTotalW[q]->Scale(1.0/EDepTotalW[q]->Integral());
				EDepTotalW[q]->Draw("histSAME");
				leg->AddEntry(EDepTotalW[q],Form("Zhaozhong (%s Polystyrene)",Length[q].Data()),"l");	
			}


			leg->Draw("SAME");
		}



		if(Energy == 8){

			EDepTotalWRef8->Scale(1.0/EDepTotalWRef8->Integral());
			EDepTotalWRef8->SetMaximum(0.32);
			EDepTotalWRef8->Draw("hist");
			lat->DrawLatex(0.12,0.85,"W Plate/Shashlik Calorimeter");		
			lat->DrawLatex(0.12,0.80,"20 Plates/Layer and 1X0 Width Layer Thickness (3.5mm)");		
			lat->DrawLatex(0.12,0.75,"Incident Electron Eneregy = 8 GeV");		
			leg->AddEntry(EDepTotalWRef8,"Cheung-Haggerty","l");

			for(int q = 0; q < NFiles; q++){
				EDepTotalW[q]->Scale(1.0/EDepTotalW[q]->Integral());
				EDepTotalW[q]->Draw("histSAME");
				leg->AddEntry(EDepTotalW[q],Form("Zhaozhong (%s Polystyrene)",Length[q].Data()),"l");	
			}


			leg->Draw("SAME");
		}

	}


	if(Option == 0){



		lat->DrawLatex(0.12,0.85,"W Plate/Shashlik Calorimeter");		
		lat->DrawLatex(0.12,0.80,Form("%d Plates/Layer and 1/4 X0 Width Layer Thickness (3.5mm)",NLayers));		
		lat->DrawLatex(0.12,0.75,"Incident Electron Eneregy = 8 GeV");		

		if(Energy == 8){



			EDepTotalWRef8->Scale(1.0/EDepTotalWRef8->Integral());
			EDepTotalWRef8->SetMaximum(0.32);
			EDepTotalWRef8->Draw("hist");

			leg->AddEntry(EDepTotalWRef8,"Cheung-Haggerty","l");

			for(int q = 0; q < NFiles; q++){
				EDepTotalW[q]->Scale(1.0/EDepTotalW[q]->Integral());
				EDepTotalW[q]->Draw("histSAME");
				leg->AddEntry(EDepTotalW[q],Form("Zhaozhong (%s Polystyrene)",Length[q].Data()),"l");	
			}


			leg->Draw("SAME");
		}


	}


	c->SaveAs(Form("ELeakAna/Comparison/%s/VisEnergyComp%d.png",Thickness.Data(),Energy));



}
