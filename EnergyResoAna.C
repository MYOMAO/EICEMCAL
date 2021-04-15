#include <cmath>
#include <TFile.h>
#include <TString.h>
#include <TLine.h>
#include <TTree.h>
#include <TLatex.h>
#include <TGraphErrors.h>
#include <cassert>
#include <iostream>
#include <fstream>
#include "TROOT.h"
#include "TH1.h"
#include "TTree.h"
#include "TH2.h"
#include "TF1.h"
#include "TFile.h"
#include "TMath.h"

using namespace std;



using std::cout;
using std::endl;



void EnergyResoAna(double EtaCutDown, 	double EtaCutUp, int NLayers){


	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();

	double ScaleFactor = 0.90; 
	double Width = 0.2;

	double ScaleDownFactor = 0.23;

	double Energy = 5;

	double EnergyCent = Energy * ScaleFactor;

	double UpperE = Energy * ScaleFactor + Energy * Width;
	double LowerE = Energy * ScaleFactor - Energy * Width;


	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);


	//TFile * fin = new TFile("G4EICDetector.root_g4cemc_eval.root");
	TFile * fin = new TFile("G4EICDetector.root_g4femc_eval.root");
	fin->cd();


	TLatex *latE = new TLatex();
	latE->SetNDC();
	latE->SetTextSize(0.03);
	latE->SetTextColor(kBlack);


	TTree * t = (TTree * ) fin->Get("ntp_gshower");



	TH1D * RealEnergy = new TH1D("RealEnergy","",100,4.9,5.1);


	t->Project("RealEnergy","ge");

	RealEnergy->GetXaxis()->SetTitle("Generated Electron (Real) Energy (GeV)");
	RealEnergy->GetYaxis()->SetTitle("Counts");
	RealEnergy->GetYaxis()->SetTitleOffset(1.3);

	RealEnergy->GetXaxis()->CenterTitle();
	RealEnergy->GetYaxis()->CenterTitle();


	RealEnergy->Draw();


	c->SaveAs(Form("AnaPlots/RealEnergy_%0.f.png",Energy));



	TH1D * RecoEnergy = new TH1D("RecoEnergy","",60,3.6,4.6);


	t->Project("RecoEnergy","e");

	RecoEnergy->GetXaxis()->SetTitle("Reconstructed Electron Energy (GeV)");
	RecoEnergy->GetYaxis()->SetTitle("Counts");
	RecoEnergy->GetYaxis()->SetTitleOffset(1.3);

	RecoEnergy->GetXaxis()->CenterTitle();
	RecoEnergy->GetYaxis()->CenterTitle();


	RecoEnergy->Draw();


	c->SaveAs(Form("AnaPlots/RecoEnergy_%0.f.png",Energy));



	//Linearity


	TH1D * EnergyRatio = new TH1D("EnergyRatio","",60,0.6,1.2);


	t->Project("EnergyRatio","e/ge");

	EnergyRatio->GetXaxis()->SetTitle("Reco Energy (From Clusters)/Real Energy");
	EnergyRatio->GetYaxis()->SetTitle("Counts");
	EnergyRatio->GetYaxis()->SetTitleOffset(1.3);

	EnergyRatio->GetXaxis()->CenterTitle();
	EnergyRatio->GetYaxis()->CenterTitle();


	EnergyRatio->SetMarkerStyle(20);
	EnergyRatio->SetMarkerSize(1);
	EnergyRatio->SetMarkerColor(kBlack);

	EnergyRatio->SetMaximum(EnergyRatio->GetMaximum()*1.5);

	EnergyRatio->Draw("ep");


	TF1 * fLinearity = new TF1("fLinearity","gaus",0.6,1.2);
	fLinearity->SetLineColor(kRed);

	EnergyRatio->Fit(fLinearity,"R");



	double LinearityMean = fLinearity->GetParameter(2);
	double LinearityMeanErr = fLinearity->GetParError(2);
	double LinearityWidth = fLinearity->GetParameter(1);
	double LinearityWidthErr = fLinearity->GetParError(1);


	latE->DrawLatex(0.20,0.85,Form("Number of Layers: %d",NLayers));
	latE->DrawLatex(0.20,0.80,Form("Rapidity Range: %.1f < |#eta| < %.1f",EtaCutDown,EtaCutUp));
	latE->DrawLatex(0.20,0.75,Form("Mean Real/Reco Energy = %.3f #pm %.3f",LinearityMean,LinearityMeanErr));
	latE->DrawLatex(0.20,0.70,Form("Width Real/Reco Energy = %.3f #pm %.3f",LinearityWidth,LinearityWidthErr));



	c->SaveAs(Form("AnaPlots/EnergyRatio/EnergyRatio_%0.f_%d.png",Energy,NLayers));


	//2D Linearity//



	TH2D * EnergyLinearity = new TH2D("EnergyLinearity","",30,3.8,5.5,30,4.6,5.4);


	t->Project("EnergyLinearity","ge:e");

	EnergyLinearity->GetXaxis()->SetTitle("Reconstructed Electron Energy (GeV)");
	EnergyLinearity->GetYaxis()->SetTitle("Generated Electron (Real) Energy (GeV)");
	EnergyLinearity->GetYaxis()->SetTitleOffset(1.3);

	EnergyLinearity->GetXaxis()->CenterTitle();
	EnergyLinearity->GetYaxis()->CenterTitle();


	EnergyLinearity->Draw("COLZ");


	c->SaveAs(Form("AnaPlots/EnergyLinearity_%0.f.png",Energy));




	TH2D * EnergyReulation = new TH2D("EnergyReulation","",30,3.0,4.6,30,0.93,1.3);


	t->Project("EnergyReulation","(ge/e):e");

	EnergyReulation->GetXaxis()->SetTitle("Reconstructed Electron Energy (GeV)");
	EnergyReulation->GetYaxis()->SetTitle("Reco/Real Energy");
	EnergyReulation->GetYaxis()->SetTitleOffset(1.3);

	EnergyReulation->GetXaxis()->CenterTitle();
	EnergyReulation->GetYaxis()->CenterTitle();



	EnergyReulation->Draw("COLZ");


	c->SaveAs(Form("AnaPlots/ErrorvsEnergy_%0.f.png",Energy));



//	LowerE = 0;
//	UpperE = 60;

	TH1D * EnergyDis = new TH1D("EnergyDis","",70,LowerE,UpperE);


	t->Project("EnergyDis","e",Form("abs(eta) > %f && abs(eta) < %f",EtaCutDown,EtaCutUp));


	EnergyDis->GetXaxis()->SetTitle("Reconstructed Electron Energy (GeV)");
	EnergyDis->GetYaxis()->SetTitle("Counts");
	EnergyDis->GetYaxis()->SetTitleOffset(1.3);

	EnergyDis->GetXaxis()->CenterTitle();
	EnergyDis->GetYaxis()->CenterTitle();

	EnergyDis->SetMarkerStyle(20);
	EnergyDis->SetMarkerSize(1.0);
	EnergyDis->SetMarkerColor(kBlack);


	TF1 * f = new TF1("f","gaus",LowerE,UpperE);
	f->SetLineColor(kRed);

	/*	
		TF1 * f = new TF1("f","[0] * TMath::Exp(-(x-[1])*(x-[1])/(2 * [2]*[2]))",LowerE,UpperE);
		f->SetParLimits(1,EnergyCent-0.1,EnergyCent+0.1);
		f->SetParLimits(0,30,60);


		EnergyDis->Fit(f);


*/
	EnergyDis->SetMaximum(EnergyDis->GetMaximum()* 1.6);
	EnergyDis->Draw("ep");

	EnergyDis->Fit("f");



	TLegend* leg = new TLegend(0.30,0.80,0.75,0.90,NULL,"brNDC");
	leg->SetBorderSize(0);
	leg->SetTextSize(0.030);
	leg->SetTextFont(42);
	leg->SetFillStyle(0);

	leg->AddEntry(EnergyDis,"EIC Forward EMCAL Reconstructed Energy","EP");
	leg->AddEntry(f,"Gaussian Fit Function","l");

	leg->Draw("SAME");




	double EnergyWidth = f->GetParameter(2);
	double EnergyWidthError = f->GetParError(2);
	double EnergyMean = f->GetParameter(1);
	double EnergyMeanError = f->GetParError(1);

	double EnergyResolution = EnergyWidth/EnergyMean;

	double EnergyResolutionError = EnergyResolution *  sqrt(EnergyWidthError/EnergyWidth * EnergyWidthError/EnergyWidth + EnergyMeanError/EnergyMean *  EnergyMeanError/EnergyMean );




	cout << "Energy = " << Energy << "      Energy Resolution = " << EnergyResolution << "    Energy Resolution Error =   " << EnergyResolutionError << endl;


	latE->DrawLatex(0.30,0.75,Form("Rapidity Range: %.1f < |#eta| < %.1f",EtaCutDown,EtaCutUp));
	latE->DrawLatex(0.30,0.70,Form("Energy Resolution = %.3f #pm %.3f",EnergyResolution,EnergyResolutionError));

	c->SaveAs(Form("AnaPlots/EnergyReso_%.0f_%.0f_%.0f.png",Energy,EtaCutDown,EtaCutUp));


	//Energy Dep//

	TH1D * EnergyDisDep = new TH1D("EnergyDisDep","",70,LowerE*ScaleDownFactor,UpperE*ScaleDownFactor);


	t->Project("EnergyDisDep","gedep",Form("abs(eta) > %f && abs(eta) < %f",EtaCutDown,EtaCutUp));



	EnergyDisDep->GetXaxis()->SetTitle("Deposited Electron Energy All Scintillators (GeV)");
	EnergyDisDep->GetYaxis()->SetTitle("Counts");
	EnergyDisDep->GetYaxis()->SetTitleOffset(1.3);

	EnergyDisDep->GetXaxis()->CenterTitle();
	EnergyDisDep->GetYaxis()->CenterTitle();

	TF1 * fdep = new TF1("fdep","gaus",LowerE*ScaleDownFactor,UpperE*ScaleDownFactor);
	fdep->SetLineColor(kRed);

	EnergyDisDep->SetMarkerStyle(20);
	EnergyDisDep->SetMarkerSize(1.0);
	EnergyDisDep->SetMarkerColor(kBlack);

	EnergyDisDep->SetMaximum(EnergyDisDep->GetMaximum()* 1.6);

	EnergyDisDep->Draw("ep");
	EnergyDisDep->Fit("fdep");


	double EnergyWidthDep = fdep->GetParameter(2);
	double EnergyWidthDepError = fdep->GetParError(2);
	double EnergyMeanDep  = fdep->GetParameter(1);
	double EnergyMeanDepError = fdep->GetParError(1);

	double EnergyResolutionDep  = EnergyWidthDep/EnergyMeanDep;

	double EnergyResolutionDepError = EnergyResolutionDep *  sqrt(EnergyWidthDepError/EnergyWidthDep * EnergyWidthDepError/EnergyWidthDep + EnergyMeanDepError/EnergyMeanDep *  EnergyMeanDepError/EnergyMeanDep );



	TLegend* legdep = new TLegend(0.30,0.80,0.75,0.90,NULL,"brNDC");
	legdep->SetBorderSize(0);
	legdep->SetTextSize(0.030);
	legdep->SetTextFont(42);
	legdep->SetFillStyle(0);

	legdep->AddEntry(EnergyDisDep,"EIC Forward EMCAL Deposited Energy","EP");
	legdep->AddEntry(fdep,"Gaussian Fit Function","l");

	legdep->Draw("SAME");

	latE->DrawLatex(0.30,0.75,Form("Rapidity Range: %.1f < |#eta| < %.1f",EtaCutDown,EtaCutUp));
	latE->DrawLatex(0.30,0.70,Form("Energy Resolution = %.3f #pm %.3f",EnergyResolutionDep,EnergyResolutionDepError));


	c->SaveAs(Form("AnaPlots/EnergyResoDep_%.0f_%.0f_%.0f.png",Energy,EtaCutDown,EtaCutUp));



	cout << "Energy = " << Energy << "      Energy Resolution Deposited = " << EnergyResolutionDep << "    Energy Resolution Deposited Error =   " << EnergyResolutionDepError << endl;

}
