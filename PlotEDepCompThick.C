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




void PlotEDepCompThick(int Option){

	int Mode = 0;




	const int NRefEnergy = 4;
	const int NRefThick = 2;
	double WRadLength = 0.35;

	double MaxHeight[NRefThick] = {0.60,0.30};


	TString Thickness[NRefThick] = {"Unity","Quarter"};
	TString ThickName[NRefThick] = {"1", "1/4"};
	int RadLength[NRefThick] = {20,4};

	int NLayers [NRefThick] = {20, 16};



	int Energy[NRefEnergy] ={1,2,4,8};

	TString infile[NRefEnergy];

	for(int i = 0; i < NRefEnergy; i ++){

		infile[i] =  Form("EnergyScanLay/%s/ShowerInfo_%d.root",Thickness[Option].Data(),Energy[i]);

	}

	int NBins = 100;
	double EMin = 10;
	double EMax = 260;


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

	double MeanConFrac;
	int color[NRefEnergy] = {1,2,3,4};
	int colorZZ[NRefEnergy] = {5,6,7,8};


	TH1D * EDepTotalW[NRefEnergy];


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



	TH1D * EDepTotalWRef[NRefEnergy][NRefThick];

	for(int i = 0; i < NRefEnergy; i++){

		for(int j = 0; j < NRefThick; j++){

			EDepTotalWRef[i][j] = new TH1D(Form("EDepTotalWRef_%d_%d",i,j),"",NBins,EMin,EMax);
			EDepTotalWRef[i][j]->GetXaxis()->SetTitle("Total Visible Energy (MeV)");
			EDepTotalWRef[i][j]->GetYaxis()->SetTitle("Fraction");
			EDepTotalWRef[i][j]->GetXaxis()->CenterTitle();
			EDepTotalWRef[i][j]->GetYaxis()->CenterTitle();
			EDepTotalWRef[i][j]->GetYaxis()->SetTitleOffset(1.4);
			EDepTotalWRef[i][j]->SetLineWidth(2);
			EDepTotalWRef[i][j]->SetLineColor(color[i]);


		}
	}




	for(int i = 0; i < EVis1GeV; i++){

		EDepTotalWRef[0][0]->Fill(EVis1GeVX[i],EVis1GeVY[i]);

	}



	for(int i = 0; i < EVis2GeV; i++){

		EDepTotalWRef[1][0]->Fill(EVis2GeVX[i],EVis2GeVY[i]);

	}



	for(int i = 0; i < EVis4GeV; i++){

		EDepTotalWRef[2][0]->Fill(EVis4GeVX[i],EVis4GeVY[i]);

	}



	for(int i = 0; i < EVis8GeV; i++){

		EDepTotalWRef[3][0]->Fill(EVis8GeVX[i],EVis8GeVY[i]);

	}


	//Quarter


	for(int i = 0; i < EVis1GeVQuarter; i++){

		EDepTotalWRef[0][1]->Fill(EVis1GeVXQuarter[i],EVis1GeVYQuarter[i]);

	}



	for(int i = 0; i < EVis2GeVQuarter; i++){

		EDepTotalWRef[1][1]->Fill(EVis2GeVXQuarter[i],EVis2GeVYQuarter[i]);

	}



	for(int i = 0; i < EVis4GeVQuarter; i++){

		EDepTotalWRef[2][1]->Fill(EVis4GeVXQuarter[i],EVis4GeVYQuarter[i]);

	}



	for(int i = 0; i < EVis8GeVQuarter; i++){

		EDepTotalWRef[3][1]->Fill(EVis8GeVXQuarter[i],EVis8GeVYQuarter[i]);

	}






	for(int i = 0; i < NPoints1GeV; i++){

		EDepLongTungs1GeV->SetBinContent(i+1,EdepCurve1GeVY[i+1]);

	}

	for(int i = 0; i < NPoints10GeV; i++){

		EDepLongTungs10GeV->SetBinContent(i+1,EdepCurve10GeVY[i+1]);

	}


	//Extra Debug//


	TH1D * EDepAll[NRefEnergy]; EnergyThickScan/Debug/


	for(int i = 0; i < NRefEnergy; i++){

		EDepAll[i] = new TH1D(Form("EDepAll_%d",i),"",NBins,0,0.95);
		EDepAll[i]->GetXaxis()->SetTitle("Electron Energy Containment Fraction");
		EDepAll[i]->GetYaxis()->SetTitle("Counts");
		EDepAll[i]->GetXaxis()->CenterTitle();
		EDepAll[i]->GetYaxis()->CenterTitle();
		EDepAll[i]->GetYaxis()->SetTitleOffset(1.4);
		EDepAll[i]->SetLineWidth(2);
		EDepAll[i]->SetLineColor(color[i] + 1);

	}


	double TotalDepFrac;

	for(int q = 0; q < NRefEnergy;q++){

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






		EDepTotalW[q]	= new TH1D(Form("EDepTotalW_%d",q),"",NBins,EMin,EMax);
		EDepTotalW[q]->GetXaxis()->SetTitle("Total Visible Energy (MeV)");
		EDepTotalW[q]->GetYaxis()->SetTitle("Fraction");
		EDepTotalW[q]->GetXaxis()->CenterTitle();
		EDepTotalW[q]->GetYaxis()->CenterTitle();
		EDepTotalW[q]->GetYaxis()->SetTitleOffset(1.4);
		EDepTotalW[q]->SetLineColor(colorZZ[q]);
		EDepTotalW[q]->SetLineWidth(2);





		double ZPos;
		double ZBin;


		for(int i = 0; i < NEntries; i++){

			ShowerInfo->GetEntry(i);

			ZPos = (Z1 + Z2 - MinDepth * 2) * WRadLength/2;



			if(evtid != evtidPre){


				EDepTotalW[q]->Fill(TotalFEMCALEnergy);
				TotalDepFrac = (Energy[q] - TotalBHEnergyBF/1000)/Energy[q];
				EDepAll[q]->Fill(TotalDepFrac);


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





	TLegend * leg  = new TLegend(0.10,0.42,0.60,0.80);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.035);
	leg->SetTextFont(42);
	leg->SetFillStyle(0);
	leg->SetLineWidth(3);



	if(Mode == 1){

		for(int i = 0; i < NRefEnergy; i++){

			TString infileRECO = Form("EnergyScanLay/%s/G4EICDetector.root_g4femc_eval_%d.root",Thickness[Option].Data(),Energy[i]);

			TFile * finRECO = new TFile(infileRECO.Data());

			TTree * ntp_gshower = (TTree *) finRECO->Get("ntp_gshower");


			EDepTotalW[i]	= new TH1D(Form("EDepTotalW_%d",i),"",NBins,EMin,EMax);
			EDepTotalW[i]->GetXaxis()->SetTitle("Total Visible Energy (MeV)");
			EDepTotalW[i]->GetYaxis()->SetTitle("Fraction");
			EDepTotalW[i]->GetXaxis()->CenterTitle();
			EDepTotalW[i]->GetYaxis()->CenterTitle();
			EDepTotalW[i]->GetYaxis()->SetTitleOffset(1.4);
			EDepTotalW[i]->SetLineColor(colorZZ[i]);
			EDepTotalW[i]->SetLineWidth(2);

			ntp_gshower->Project(Form("EDepTotalW_%d",i),"e * 1000");


		}
	}

	for(int q = 0; q < NRefEnergy; q++){

		MeanConFrac = EDepAll[q]->GetMean() * 100;
		EDepAll[q]->SetMaximum(EDepAll[q]->GetMaximum()*1.3);
		EDepAll[q]->SetTitle(Form("Total EMCAL Radiation Length = %d X0",RadLength[Option]));
		EDepAll[q]->Draw("hist");
		lat->DrawLatex(0.12,0.85,Form("Incident Electron Beam Energy = %d GeV",Energy[q]));		
		lat->DrawLatex(0.12,0.80,Form("Mean Containment Fraction = %.1f %%",MeanConFrac));		
		lat->DrawLatex(0.12,0.75,"John's Note = 39.2 %");	

		TLine *line4X0 = new TLine(0.392,0,0.392,EDepAll[q]->GetMaximum());
		line4X0->SetLineStyle(2);
		line4X0->SetLineWidth(2);
		line4X0->SetLineColor(4);
		line4X0->Draw("SAME");

		c->SaveAs(Form("NewComp/EMCALContainment_%d_%s.png",Energy[q],Thickness[Option].Data()));

	}





	for(int q = 0; q < NRefEnergy; q++){

		EDepTotalWRef[q][Option]->Scale(1.0/EDepTotalWRef[q][Option]->Integral());

		if(q == 0){ 
			EDepTotalWRef[q][Option]->SetMaximum(MaxHeight[Option]);
			EDepTotalWRef[q][Option]->SetTitle(Form("Total EMCAL Radiation Length = %d X0",RadLength[Option]));
			EDepTotalWRef[q][Option]->Draw("hist");
		}
		if(q != 0) EDepTotalWRef[q][Option]->Draw("histSAME");

		leg->AddEntry(EDepTotalWRef[q][Option],Form("Cheung-Haggerty %d GeV",Energy[q]),"l");


		EDepTotalW[q]->Scale(1.0/EDepTotalW[q]->Integral());
		EDepTotalW[q]->Draw("histSAME");
		leg->AddEntry(EDepTotalW[q],Form("Zhaozhong (1.4 mm Polystyrene) - %d GeV",Energy[q]),"l");	
	}


	lat->DrawLatex(0.12,0.85,"W Plate/Shashlik Calorimeter");		
	lat->DrawLatex(0.12,0.80,Form("%d Plates/Layer and %s X0 Width Layer Thickness (3.5mm)",NLayers[Option],ThickName[Option].Data()));		

	leg->Draw("SAME");


	if(Mode == 0) c->SaveAs(Form("NewComp/VisComp_%s.png",Thickness[Option].Data()));
	if(Mode == 1) c->SaveAs(Form("NewComp/VisComp_RECO_%s.png",Thickness[Option].Data()));



	TLegend * Refleg  = new TLegend(0.10,0.52,0.60,0.85);
	Refleg->SetBorderSize(0);
	Refleg->SetTextSize(0.035);
	Refleg->SetTextFont(42);
	Refleg->SetFillStyle(0);
	Refleg->SetLineWidth(3);


	for(int r = 0; r < NRefThick; r++){
		for(int q = 0; q < NRefEnergy; q++){
			EDepTotalWRef[q][r]->Scale(1.0/EDepTotalWRef[q][r]->Integral());
			if(r == 0) EDepTotalWRef[q][r]->SetLineColor(color[q]);
			if(r == 1) EDepTotalWRef[q][r]->SetLineColor(colorZZ[q]);		
			if(q == 0 && r == 0){ 
				EDepTotalWRef[q][r]->SetMaximum(MaxHeight[Option]);
				EDepTotalWRef[q][r]->SetTitle("John's Note Reference");
				EDepTotalWRef[q][r]->Draw("hist");
				
			}
			else{
				EDepTotalWRef[q][r]->Draw("histSAME");
			}
			
			Refleg->AddEntry(EDepTotalWRef[q][r],Form("%d GeV - %s X0 Sampling Thickness %d X0 Rad Length",Energy[q],ThickName[r].Data(),RadLength[r]),"l");


		}
		
	}

	Refleg->Draw("SAME");
	c->SaveAs("NewComp/JohnRefComp.png");

}
