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




void EnergyThickScan(){
	
	gStyle->SetOptStat(0);

	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();

	const int NThick = 5;
	const int NEnergy = 9;

	int evtid;
	float FEMCALedep;
	int XBin;

	TString ThickName[NThick] = {"Quarter","Half","1","2","4"}; 
	TString ThickName2[NThick] = {"1/4","1/2","1","2","4"}; 

	int Energy[NEnergy] = {1,2,3,4,5,6,7,8,9};

	TString infile;

	double MeanVisEnergy;
	double RMSVisEnergy;
	
	double TotalEnergy = 0;


	TH1D * ThickEnergyScan[NThick];


	//Reference//
	
	double x1[NThick] = {3.3246384297520657, 3.230544077134986,2.978822314049587,4.38752582644628,6.046229338842975};
	double x2[NThick] = {5.100809228650138, 6.07433712121212, 6.7535296143250685,6.649492079889806,8.569344008264464};
	double y1[NThick] = {0.34318181818181814, 0.17499999999999993, 0.084090909090909,0.06363636363636371,0.04318181818181799};
	double y2[NThick] = {0.525, 0.32727272727272727, 0.1886363636363636,0.09545454545454546,0.06363636363636371};

	double a[NThick];
	double b[NThick];

	TF1 * f[NThick];


	for(int i  = 0; i < NThick; i++){

		a[i] = (y2[i] - y1[i])/(x2[i] - x1[i]);
		b[i] = y1[i] - (y2[i] - y1[i])/(x2[i] - x1[i]) * x1[i];
		f[i] = new TF1("f1",Form("%f * x + %f",a[i],b[i]),0,10);
		f[i]->SetLineColor(i+1);
		f[i]->SetLineWidth(2);
	}


	



	for(int i = 0; i < NThick; i++){

		int NBins = 100;
		int EMin = 0;
		int EMax = 10;

		ThickEnergyScan[i] = new TH1D(Form("ThickEnergyScan_%d",i),"",NBins,EMin,EMax);
		ThickEnergyScan[i]->GetXaxis()->SetTitle("Incident Electron Energy (GeV)");
		ThickEnergyScan[i]->GetYaxis()->SetTitle("Visible Energy (GeV)");
		ThickEnergyScan[i]->GetYaxis()->SetTitleOffset(1.3);
	
		ThickEnergyScan[i]->GetXaxis()->CenterTitle();
		ThickEnergyScan[i]->GetYaxis()->CenterTitle();
		ThickEnergyScan[i]->SetMarkerColor(i+1);
		ThickEnergyScan[i]->SetMarkerStyle(20);
		ThickEnergyScan[i]->SetMarkerSize(1);

		TotalEnergy = 0;

		for(int j = 0; j < NEnergy; j++){

			infile = Form("EnergyScanThick/%s/ShowerInfo_%d.root",ThickName[i].Data(),Energy[j]);
			TFile * fin = new TFile(infile.Data());
			fin->cd();

			TTree * ShowerInfo = (TTree *) fin->Get("ShowerInfo");

			ShowerInfo->SetBranchAddress("evtid",&evtid);
			ShowerInfo->SetBranchAddress("FEMCALedep",&FEMCALedep);

			int NEntries = ShowerInfo->GetEntries();
			int evtidPre = 0;
			float TotalFEMCALEnergy = 0;
	
			TH1D * TotalEnergy = new TH1D("TotalEnergy","",100,0.01,1.0);

			for(int k = 0; k < NEntries; k++){

				ShowerInfo->GetEntry(k);



				if(evtid != evtidPre){


	
					TotalEnergy->Fill(TotalFEMCALEnergy);



					TotalFEMCALEnergy = 0;


				}


				TotalFEMCALEnergy = FEMCALedep + TotalFEMCALEnergy;

				evtidPre = evtid;

			}

		//	TotalEnergy->Draw();
		//	c->SaveAs(Form("EnergyThickScan/Debug/EnergyDis%s_%d.png",ThickName[i].Data(),Energy[j]));
			
			MeanVisEnergy = TotalEnergy->GetMean();
			RMSVisEnergy = TotalEnergy->GetRMS();

			XBin = ThickEnergyScan[i]->GetXaxis()->FindBin(Energy[j]);
			ThickEnergyScan[i]->SetBinContent(XBin,MeanVisEnergy);
			ThickEnergyScan[i]->SetBinError(XBin,RMSVisEnergy);

		}


	
		

	}

	
	


	TLegend * leg  = new TLegend(0.10,0.50,0.50,0.88);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.030);
	leg->SetTextFont(42);
	leg->SetFillStyle(0);
	leg->SetLineWidth(3);






	for(int i = 0; i < NThick; i++){
		
		if(i == 0){	
			ThickEnergyScan[i]->SetMaximum(1.0);
			ThickEnergyScan[i]->Draw("ep");
	
		}
			if(i > 0)	ThickEnergyScan[i]->Draw("epSAME");

		leg->AddEntry(ThickEnergyScan[i],Form("Zhaozhong Simulation: %s X0 W-Shashlik EMCAL",ThickName2[i].Data()),"pl");	
		f[i]->Draw("SAME");
		leg->AddEntry(f[i],Form("John's Note: %s X0 W EMCAL",ThickName2[i].Data()),"l");	

	}


	leg->Draw("SAME");

	


	c->SaveAs("EnergyThickScan/20X0EScanComp.png");




}


