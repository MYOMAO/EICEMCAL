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


void ClusterAnalyzer(int BeamOpt){

	gStyle->SetOptStat(0);

	TString BeamName;

	if(BeamOpt == 0) BeamName = "e";
	if(BeamOpt == 1) BeamName = "pi0";

	TString BeamNameFull;

	if(BeamOpt == 0) BeamNameFull = "e^{-}";
	if(BeamOpt == 1) BeamNameFull = "pi^{0}";


	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();
	
	const int NFiles = 11;
	double RatioUp = 1;
	double RatioDown = 0.6;
	int NBins = 100;

	double EStep = 1.0;
	double Energy;
	


	double ERECO;
	double ERECOError;
	int XBinERECO;



	TH1D * Linearity = new TH1D("Linearity","",100,0,25);
	TH1D * LinearityRatio = new TH1D("LinearityRatio","",100,0,25);


	for(int q = 1; q < NFiles; q++){

		Energy = (q+1) * EStep;

		TString infile = Form("ERECOFiles/%s/G4EICDetector.root_g4femc_eval_%d.root",BeamName.Data(),q);

		TFile * fin = new TFile(infile.Data());
		fin->cd();
		TTree * ntp_gshower = (TTree * ) fin->Get("ntp_gshower");
	
		
		TH1D * ERECODis = new TH1D("ERECODis","",NBins,Energy*RatioDown,Energy*RatioUp);
		
		ntp_gshower->Project("ERECODis","e");
		
		ERECODis->GetXaxis()->SetTitle(Form("Reconstructed %s Energy",BeamNameFull.Data()));
		ERECODis->GetYaxis()->SetTitle("Counts");
		ERECODis->GetXaxis()->CenterTitle();
		ERECODis->GetYaxis()->CenterTitle();
		ERECODis->GetYaxis()->SetTitleOffset(1.4);

		ERECODis->SetMarkerStyle(20);
		ERECODis->SetMarkerSize(1);
		ERECODis->SetMarkerColor(kBlack);
		ERECO = ERECODis->GetMean();
		ERECOError = ERECODis->GetRMS();
			

		XBinERECO = Linearity->GetXaxis()->FindBin(Energy);

		Linearity->SetBinContent(XBinERECO,ERECO);
		Linearity->SetBinError(XBinERECO,ERECOError);
		
		LinearityRatio->SetBinContent(XBinERECO,ERECO/Energy);
		LinearityRatio->SetBinError(XBinERECO,ERECOError/Energy);


		ERECODis->Draw("ep");
		c->SaveAs(Form("RECOPlots/%s/RECOEnergyDis/ERECODis_%d.png",BeamName.Data(),q+1));

	}
	
	Linearity->GetXaxis()->SetTitle(Form("Incident %s Energy",BeamNameFull.Data()));
	Linearity->GetYaxis()->SetTitle(Form("Reconstructed %s Energy",BeamNameFull.Data()));
	Linearity->GetXaxis()->CenterTitle();
	Linearity->GetYaxis()->CenterTitle();
	Linearity->GetYaxis()->SetTitleOffset(1.4);


	Linearity->SetMarkerStyle(20);
	Linearity->SetMarkerSize(1);
	Linearity->SetMarkerColor(kBlack);

	Linearity->Draw("ep");

	c->SaveAs(Form("RECOPlots/%s/Summary/EnergyLinearity.png",BeamName.Data()));

	
	
	LinearityRatio->GetXaxis()->SetTitle(Form("Incident %s Energy",BeamNameFull.Data()));
	LinearityRatio->GetYaxis()->SetTitle(Form("Reconstructed %s Energy/Incident %s Energy",BeamNameFull.Data(),BeamNameFull.Data()));
	LinearityRatio->GetXaxis()->CenterTitle();
	LinearityRatio->GetYaxis()->CenterTitle();
	LinearityRatio->GetYaxis()->SetTitleOffset(1.4);


	LinearityRatio->SetMarkerStyle(20);
	LinearityRatio->SetMarkerSize(1);
	LinearityRatio->SetMarkerColor(kBlack);

	LinearityRatio->Draw("ep");

	c->SaveAs(Form("RECOPlots/%s/Summary/EnergyLinearityRatio.png",BeamName.Data()));

}
