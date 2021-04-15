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
using namespace std;

using std::cout;
using std::endl;



void EnergyDepthAna(){


	double AverageE;
	int NHits;
	double TotalE;
	double AverageEerr;
	double NowE;
	int NowCounts;



	double MinDepth = 290;
	double MaxDepth = 330;
	double Width = 0.56;

	double DepthWidth = MaxDepth - MinDepth;

	int NBins = (MaxDepth - MinDepth)/Width + 1; 


	gStyle->SetOptStat(0);

	const int NFile = 20;

	double Energy;
	double EnergyMin = 0.5;
	double EnergyStep = 0.5;

	const int NFilePlot = 5;  

	int Step = NFile/NFilePlot;

	int iFile;

	TString infile;

	TH1D * EnergyDepth[NFile]; 



	TCanvas * c4 = new TCanvas("c4","c4",600,600);
	c4->cd();




	TLatex *lat = new TLatex();
	lat->SetNDC();
	lat->SetTextSize(0.035);
	lat->DrawLatex(0.30,0.80,"Forward Black Hole Energy Response");	
	



	for(int q = 0; q < NFile; q++){

		Energy  = q * EnergyStep + EnergyMin;


		cout << "Now Working on File index = " << q << endl;


		infile = Form("EnergyScan/G4EICDetector.root_DSTReader_%d.root",q);

		TFile * fin = new TFile(infile.Data());
		TTree * T = (TTree *) fin->Get("T");




		TH2D * EdepDepth = new TH2D("EdepDepth","",NBins,MinDepth,MaxDepth,10000,0,0.020);
		T->Project("EdepDepth","G4HIT_FEMC.edep:G4HIT_FEMC.z");


		EdepDepth->Draw();

		c4->SaveAs("LongStudy/EdepDepth.png");


		EnergyDepth[q] = new TH1D(Form("EnergyDepth_%d",q),Form("EnergyDepth_%d",q),NBins,0,DepthWidth);
		EnergyDepth[q]->GetXaxis()->SetTitle("Depth: t (cm)");
		EnergyDepth[q]->GetYaxis()->SetTitle("Average Energy Depositon: dE/dt (MeV/cm)");
		EnergyDepth[q]->GetYaxis()->SetTitleOffset(1.4);

		EnergyDepth[q]->GetXaxis()->CenterTitle();
		EnergyDepth[q]->GetYaxis()->CenterTitle();


		for(int i = 0; i < NBins; i++){

			AverageE = 0;
			AverageEerr = 0;
			TotalE = 0;
			NHits = 0;
			for(int j = 0 ; j < 10000; j++){

				NowCounts = EdepDepth->GetBinContent(i+1,j+1);


				if(NowCounts > 0){
					NowE = EdepDepth->GetYaxis()->GetBinCenter(j+1) * 1000;

					TotalE = TotalE + NowE;

					NHits = NHits + 1;
					//		cout << "NowE = " << NowE << endl;
				}
			}

			if(NHits > 0){
				AverageE = TotalE/NHits/Width;
				AverageEerr = AverageE/sqrt(NHits);
			}


			EnergyDepth[q]->SetBinContent(i+1,AverageE);
			EnergyDepth[q]->SetBinError(i+1,AverageEerr);

		}

		EnergyDepth[q]->SetMarkerStyle(20);
		EnergyDepth[q]->SetMarkerColor(1);
		EnergyDepth[q]->SetMarkerSize(1);


		//	EnergyDepth[q]->Draw("ep");


		TF1 * f1 = new TF1("f1","[0] * x**[1] * TMath::Exp(-[2]*x)",3,36);
		f1->SetParLimits(0,0,10);
		f1->SetParLimits(1,0,1);
		f1->SetParLimits(2,0,5);

		//f1->SetLineColor(q+1);
		f1->SetLineColor(kMagenta);

		EnergyDepth[q]->Fit(f1,"R");

		float p0 = f1->GetParameter(0);
		float p1 = f1->GetParameter(1);
		float p2 = f1->GetParameter(2);


		EnergyDepth[q]->Draw("ep");
		f1->Draw("SAME");


		lat->DrawLatex(0.15,0.30,Form("Incident Electron Energy = %.1f GeV",Energy));	
		lat->DrawLatex(0.15,0.25,Form("Fit Function: dE/dt = %.2f * t ^ %.2f *  e^ (-%.2f * t)",p0,p1,p2));	

		c4->SaveAs(Form("LongStudy/EnergyScan/LongEnergyProfile_%d.png",q));

	
	}


	TLegend * leg  = new TLegend(0.60,0.70,0.90,0.90);
	leg->SetTextSize(0.035);
	leg->SetTextFont(42);



	for(int q = 0; q < NFilePlot; q++){

		iFile = q * Step;

		Energy = iFile * EnergyStep + EnergyMin;
	
		EnergyDepth[iFile]->SetMarkerColor(q+1);
		EnergyDepth[iFile]->SetLineColor(q+1);
		if(q == 0){
			EnergyDepth[iFile]->SetTitle("");
			EnergyDepth[iFile]->SetMaximum(EnergyDepth[iFile]->GetMaximum() * 2);
			EnergyDepth[iFile]->Draw();
		}
		if(q > 0)	EnergyDepth[iFile]->Draw("SAME");
	
		

		leg->AddEntry(EnergyDepth[iFile],Form("E_{e^{-}} = %.1f GeV",Energy),"lp");

	}


	leg->Draw("SAME");



	c4->SaveAs("LongStudy/Inclusive/Comparison.png");


}
