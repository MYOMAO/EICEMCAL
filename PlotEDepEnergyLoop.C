#include "TROOT.h"
#include "TH1.h"
#include "TTree.h"
#include "TH2.h"
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
using namespace std;

using std::cout;
using std::endl;



void PlotEDepEnergyLoop(){


	gStyle->SetOptStat(0);

	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();

	c->SetLogy();



	TLatex *lat = new TLatex();
	lat->SetNDC();
	lat->SetTextSize(0.035);


	double Energy;
	double Energymin = 0.5;
	double Energystep = 0.5;

	double ELeakMean;
	double ELeakErr;

	int NFiles = 20;

	double Energymax =  0.5 + Energystep * NFiles;


	TH1D * EnergyLeak = new TH1D("EnergyLeak","",NFiles,Energymin-Energystep/2, Energymax-Energystep/2);

	EnergyLeak->GetXaxis()->SetTitle("Incident Electron Energy (GeV)");
	EnergyLeak->GetYaxis()->SetTitle("Average Energy Back Leakage (MeV)");
	EnergyLeak->GetXaxis()->CenterTitle();
	EnergyLeak->GetYaxis()->CenterTitle();

	EnergyLeak->SetMarkerSize(1.2);
	EnergyLeak->SetMarkerStyle(20);
	EnergyLeak->SetMarkerColor(1);

	TString infile;



	for(int q = 0; q < NFiles; q++){

		cout << "Now Working on File = " << q << endl;

		Energy = Energymin + Energystep * q;

		infile	= Form("EnergyScan/G4EICDetector.root_DSTReader_%d.root",q);


		TFile * fin = new TFile(infile.Data());

		TTree * T = (TTree *) fin->Get("T");


		TH1D * FEMALEdep = new TH1D("FEMALEdep","",100,0,6);

		TH1D * BHForwardEdep = new TH1D("BHForwardEdep","",100,0,6);
		TH1D * BHBackwardEdep = new TH1D("BHBackwardEdep","",100,0,6);

		FEMALEdep->GetXaxis()->SetTitle("Energy Deposited (MeV) Per Hit");
		FEMALEdep->GetYaxis()->SetTitle("Number of Hits");
		FEMALEdep->SetTitle("G4 Electron Energy Deposition In One Event");

		FEMALEdep->GetXaxis()->CenterTitle();
		FEMALEdep->GetYaxis()->CenterTitle();

		FEMALEdep->GetXaxis()->SetTitleOffset(1.2);
		FEMALEdep->GetYaxis()->SetTitleOffset(1.4);

		FEMALEdep->SetMinimum(0.1);




		T->Project("FEMALEdep","G4HIT_FEMC.edep * 1000","","",1,0);

		T->Project("BHForwardEdep","G4HIT_BH_FORWARD_PLUS.edep * 1000","","",1,100);


		T->Project("BHBackwardEdep","G4HIT_BH_FORWARD_NEG.edep * 1000","","",1,100);




		FEMALEdep->SetLineWidth(2);
		BHForwardEdep->SetLineWidth(2);
		BHBackwardEdep->SetLineWidth(2);


		FEMALEdep->SetLineColor(kGreen);
		BHForwardEdep->SetLineColor(kRed);
		BHBackwardEdep->SetLineColor(kBlue);

		FEMALEdep->Draw();


		BHForwardEdep->Draw("SAME");
		BHBackwardEdep->Draw("SAME");



		TLegend * leg  = new TLegend(0.40,0.60,0.85,0.85);
		leg->SetTextSize(0.035);
		leg->SetTextFont(42);
		leg->AddEntry(FEMALEdep,"Forward EMCAL","l");
		leg->AddEntry(BHForwardEdep,"Forward Leakage","l");	
		leg->AddEntry(BHBackwardEdep,"Backward Leakage","l");
		leg->Draw("SAME");


		c->SaveAs("ELeakAnaEnergy/SingleEvent.png");


		int NEvents = T->GetEntries();

		/*

		   cout << "NEvents =  " << NEvents << endl;

		   int NFEMCALHits;
		   int NPLUSBHHits;
		   int NNEGBHHits;







		   T->SetBranchAddress("n_G4HIT_FEMC",&NFEMCALHits);
		   T->SetBranchAddress("G4HIT_BH_FORWARD_PLUS",&NPLUSBHHits);
		   T->SetBranchAddress("G4HIT_BH_FORWARD_NEG",&NNEGBHHits);



		   const int NHits = 5000;

		   float G4HIT_FEMC_edep[NHits];
		   float G4HIT_BH_FORWARD_PLUS_edep[NHits];
		   float G4HIT_BH_FORWARD_NEG_edep[NHits];


		   T->SetBranchAddress("G4HIT_FEMC.edep",G4HIT_FEMC_edep);
		   T->SetBranchAddress("G4HIT_BH_FORWARD_PLUS.edep",G4HIT_BH_FORWARD_PLUS_edep);
		   T->SetBranchAddress("G4HIT_BH_FORWARD_NEG.edep",G4HIT_BH_FORWARD_NEG_edep);


		   float TotalFEMC;
		   float TotalBackwardBH;
		   float TotalForwardBH;



		   TH1D * TotalFEMCHis = new TH1D("TotalFEMCHis","",100,0,1000);
		   TH1D * TotalForwardBHHis = new TH1D("TotalForwardBHHis","",100,0,1000);
		   TH1D * TotalBackwardBHHis = new TH1D("TotalBackwardBHHis","",100,0,1000);


		   TotalFEMCHis->GetXaxis()->SetTitle("Total Deposited Energy Per Event (MeV)");
		   TotalFEMCHis->GetYaxis()->SetTitle("Number of Events");

		   TotalFEMCHis->GetXaxis()->CenterTitle();
		   TotalFEMCHis->GetYaxis()->CenterTitle();

		   TotalFEMCHis->GetXaxis()->SetTitleOffset(1.2);
		   TotalFEMCHis->GetYaxis()->SetTitleOffset(1.4);

*/
	
		float TotalEnergyForward;
		float MeanEnergyForward;


		float TotalFEMC;
		float TotalBackwardBH;
		float TotalForwardBH;

		TH1D * FEMALEdepTest = new TH1D("FEMALEdepTest","",5000,0,2000);
		TH1D * BHForwardEdepTest = new TH1D("BHForwardEdepTest","",5000,0,2000);
		TH1D * BHBackwardEdepTest = new TH1D("BHBackwardEdepTest","",5000,0,2000);



		TH1D * TotalFEMCHis = new TH1D("TotalFEMCHis","",25,1.6,2);
		TH1D * TotalForwardBHHis = new TH1D("TotalForwardBHHis","",5000,0,2000);
		TH1D * TotalBackwardBHHis = new TH1D("TotalBackwardBHHis","",25,0,0.02);


		TotalFEMCHis->GetXaxis()->SetTitle("Total Deposited Energy Per Event (MeV)");
		TotalFEMCHis->GetYaxis()->SetTitle("Number of Events");

		TotalFEMCHis->GetXaxis()->CenterTitle();
		TotalFEMCHis->GetYaxis()->CenterTitle();

		TotalFEMCHis->GetXaxis()->SetTitleOffset(1.2);
		TotalFEMCHis->GetYaxis()->SetTitleOffset(1.4);



		TotalForwardBHHis->GetXaxis()->SetTitle("Total Deposited Energy Per Event (MeV)");
		TotalForwardBHHis->GetYaxis()->SetTitle("Number of Events");

		TotalForwardBHHis->GetXaxis()->CenterTitle();
		TotalForwardBHHis->GetYaxis()->CenterTitle();

		TotalForwardBHHis->GetXaxis()->SetTitleOffset(1.2);
		TotalForwardBHHis->GetYaxis()->SetTitleOffset(1.4);


		TotalBackwardBHHis->GetXaxis()->SetTitle("Total Deposited Energy Per Event (MeV)");
		TotalBackwardBHHis->GetYaxis()->SetTitle("Number of Events");

		TotalBackwardBHHis->GetXaxis()->CenterTitle();
		TotalBackwardBHHis->GetYaxis()->CenterTitle();

		TotalBackwardBHHis->GetXaxis()->SetTitleOffset(1.2);
		TotalBackwardBHHis->GetYaxis()->SetTitleOffset(1.4);





		cout << "Pass here" << endl;

		for(int i = 0; i < NEvents; i++){

			//	cout << "i = " << i << endl;

			/*
			   T->GetEntry(i);

			   cout << "Got Entry " << endl;

			   TotalFEMC = 0;
			   TotalForwardBH = 0;
			   TotalBackwardBH = 0;



			   for(int j = 0; j < NFEMCALHits; j++){

			   TotalFEMC = G4HIT_FEMC_edep[j] + TotalFEMC;

			   }

			   for(int j = 0; j < NPLUSBHHits; j++){

			   TotalForwardBH = G4HIT_BH_FORWARD_PLUS_edep[j] + TotalBackwardBH;

			   }

			   for(int j = 0; j < NNEGBHHits; j++){

			   TotalBackwardBH = G4HIT_BH_FORWARD_NEG_edep[j] + TotalForwardBH;

			   }

*/




			//cout << "Project " << endl;
			T->Project("FEMALEdepTest","G4HIT_FEMC.edep * 1000","","",1,i);
			T->Project("BHForwardEdepTest","G4HIT_BH_FORWARD_PLUS.edep * 1000","","",1,i);
			T->Project("BHBackwardEdepTest","G4HIT_BH_FORWARD_NEG.edep * 1000","","",1,i);

			TotalEnergyForward = 0;
			MeanEnergyForward = 0;


			for(int j = 0; j < BHForwardEdepTest->GetNbinsX(); j++){
	
				TotalEnergyForward = BHForwardEdepTest->GetBinContent(j+ 1) * BHForwardEdepTest->GetBinCenter(j+1)  + TotalEnergyForward;

			}

			TotalFEMC =	FEMALEdepTest->Integral();
			TotalForwardBH = BHForwardEdepTest->Integral();
			TotalBackwardBH = BHBackwardEdepTest->Integral();
			
			//return;
			//		cout << " TotalForwardBH = " << TotalForwardBH << endl;
			//		cout << "TotalBackwardBH = " << TotalBackwardBH << endl;
			//		cout << "TotalFEMC = " << TotalFEMC << endl;

	
			//MeanEnergyForward = TotalEnergyForward/TotalForwardBH;

			//if(i%100 == 0) cout << "MeanEnergyForward = " << MeanEnergyForward << endl;


			TotalFEMCHis->Fill(TotalFEMC);
			TotalForwardBHHis->Fill(TotalEnergyForward);
			TotalBackwardBHHis->Fill(TotalBackwardBH);

		}
	


		cout << "TotalForwardBHHis->GetNbinsX() = " << TotalForwardBHHis->GetNbinsX() << endl;


		TotalFEMCHis->SetLineWidth(2);
		TotalForwardBHHis->SetLineWidth(2);
		TotalBackwardBHHis->SetLineWidth(2);


		TotalFEMCHis->SetLineColor(kGreen);
		TotalForwardBHHis->SetLineColor(kRed);
		TotalBackwardBHHis->SetLineColor(kBlue);


		TCanvas * c2 = new TCanvas("c2","c2",600,600);
		c2->cd();
		//	c2->SetLogx();


		TotalForwardBHHis->Draw();	
		lat->DrawLatex(0.30,0.80,"Forward Black Hole Energy Response");	

		c2->SaveAs(Form("ELeakAnaEnergy/LeakSumDisForward_%d.png",q));

		TotalBackwardBHHis->Draw();
		c2->SaveAs(Form("ELeakAnaEnergy/LeakSumDisBackward_%d.png",q));

		TotalFEMCHis->Draw();
		c2->SaveAs(Form("ELeakAnaEnergy/LeakSumDisFEMCAL_%d.png",q));


		ELeakMean = TotalForwardBHHis->GetMean();
		ELeakErr = TotalForwardBHHis->GetRMS();


		EnergyLeak->SetBinContent(q+1,ELeakMean);
		EnergyLeak->SetBinError(q+1,ELeakMean/sqrt(1000));


	}

	TCanvas * c3 = new TCanvas("c3","c3",600,600);
	c3->cd();

	EnergyLeak->Draw("ep");

	TF1 * f = new TF1("f","[0] + [1] * x + [2] * x * x",Energymin-Energystep/2, Energymax-Energystep/2);
	EnergyLeak->Fit(f,"R");
	double p0 = f->GetParameter(0);
	double p1 = f->GetParameter(1);
	double p2 = f->GetParameter(2);


	lat->DrawLatex(0.15,0.80,Form("Fit Function: y = %.1f + %.1f * x +  %.1f * x^{2}",p0,p1,p2));
	lat->DrawLatex(0.15,0.75,"Number of Layers = 66");

	c3->SaveAs("ELeakAnaEnergy/All/EnergyLeak.png");



}
