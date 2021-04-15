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
#include "TCanvas.h"
#include "TStyle.h"


#include "sPhenixStyle.C"
#include "sPhenixStyle.h"


using namespace std;



using std::cout;
using std::endl;



void EnergyResoAnaLoop(){

	gStyle->SetPadLeftMargin(0.15);
	gStyle->SetPadBottomMargin(0.13);
	gStyle->SetOptStat(0);


	SetsPhenixStyle();
	

	double ScaleDownFactor = 0.13;


	double EtaMax = 1.24;
	double EtaMin = 4.00;

//	const int NetaCut = 5;
	const int NetaCut = 1;

	double EtaStep = (EtaMax - EtaMin)/NetaCut;

	double EtaCutUp[NetaCut];
	double EtaCutDown[NetaCut];

	for(int i = 0; i < NetaCut; i++){
		EtaCutDown[i] = EtaMin + EtaStep * i;
		EtaCutUp[i] = EtaMin + EtaStep * (i+1);

	}





	double ScaleFactor = 0.90; 
	double Width = 0.2;


	int InitFile = 10;
	int NFiles = 250;

	int Index = 0;

	int StepSize = 5;
	
	int FileIndex = 2222715;

//	int TotalFiles = (NFiles - InitFile)/StepSize;
	int TotalFiles = (NFiles - InitFile);

	double EnergyStep = 0.10; 
	double Emin = EnergyStep * InitFile;

	double Emax = EnergyStep * NFiles;

//	gStyle->SetOptTitle(0);

	double Energy;

	TString InfileName;
	TFile * fin;


	TH1D * EnergyResolutionHis = new TH1D("EnergyResolutionHis","",TotalFiles ,Emin,Emax);
	EnergyResolutionHis->GetXaxis()->SetTitle("Geant4 truth electron energy (GeV)");
	EnergyResolutionHis->GetYaxis()->SetTitle("#Delta E/E");
//	EnergyResolutionHis->GetYaxis()->SetTitleOffset(1.8);
//	EnergyResolutionHis->GetXaxis()->SetTitleOffset(1.5);
	EnergyResolutionHis->GetYaxis()->SetTitleSize(0.05);
	EnergyResolutionHis->GetXaxis()->SetTitleSize(0.05);
	EnergyResolutionHis->GetXaxis()->SetTitleOffset(1.0);
	EnergyResolutionHis->GetYaxis()->SetTitleOffset(1.7);


	EnergyResolutionHis->GetXaxis()->CenterTitle();
	EnergyResolutionHis->GetYaxis()->CenterTitle();

	EnergyResolutionHis->SetMarkerStyle(20);
	EnergyResolutionHis->SetMarkerColor(kBlack);
	EnergyResolutionHis->SetMarkerSize(0.7);



	TH1D * EnergyResolutionDepHis = new TH1D("EnergyResolutionDepHis","",TotalFiles ,Emin,Emax);
	EnergyResolutionDepHis->GetXaxis()->SetTitle("Geant4 truth electron energy (GeV)");
	EnergyResolutionDepHis->GetYaxis()->SetTitle("#Delta Edep/Edep");
//	EnergyResolutionDepHis->GetYaxis()->SetTitleOffset(1.8);
//	EnergyResolutionDepHis->GetXaxis()->SetTitleOffset(1.6);
	EnergyResolutionDepHis->GetYaxis()->SetTitleSize(0.05);
	EnergyResolutionDepHis->GetXaxis()->SetTitleSize(0.05);	
	EnergyResolutionDepHis->GetXaxis()->SetTitleOffset(1.0);
	EnergyResolutionDepHis->GetYaxis()->SetTitleOffset(1.7);

	EnergyResolutionDepHis->GetXaxis()->CenterTitle();
	EnergyResolutionDepHis->GetYaxis()->CenterTitle();

	EnergyResolutionDepHis->SetMarkerStyle(20);
	EnergyResolutionDepHis->SetMarkerColor(kBlack);
	EnergyResolutionDepHis->SetMarkerSize(0.7);



	float eta;
	float e;
	float ge;


	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();


	TH1D * StatTermEta = new TH1D("StatTermEta","",5,0,0.1);
	StatTermEta->GetXaxis()->SetTitle("|#eta|");
	StatTermEta->GetYaxis()->SetTitle("Statistical Term Term");
	StatTermEta->GetYaxis()->SetTitleOffset(1.5);
	StatTermEta->GetXaxis()->CenterTitle();
	StatTermEta->GetYaxis()->CenterTitle();


	StatTermEta->SetMarkerStyle(20);
	StatTermEta->SetMarkerColor(kBlack);
	StatTermEta->SetMarkerSize(0.7);


	TH1D * ConstTermEta = new TH1D("ConstTermEta","",5,0,0.1);
	ConstTermEta->GetXaxis()->SetTitle("|#eta|");
	ConstTermEta->GetYaxis()->SetTitle("Constant Term Term");
	ConstTermEta->GetYaxis()->SetTitleOffset(1.5);
	ConstTermEta->GetXaxis()->CenterTitle();
	ConstTermEta->GetYaxis()->CenterTitle();

	ConstTermEta->SetMarkerStyle(20);
	ConstTermEta->SetMarkerColor(kBlack);
	ConstTermEta->SetMarkerSize(0.7);


	TH2D * Linearity[NetaCut];

	for(int i = 0; i < NetaCut; i++){

		Linearity[i] = new TH2D(Form("Linearity%d",i),"",TotalFiles,Emin,Emax,100,0.5,1.2);
		Linearity[i]->GetXaxis()->SetTitle("Gen E (GeV)");
		Linearity[i]->GetYaxis()->SetTitle("Vis E/Gen E");
		Linearity[i]->GetYaxis()->SetTitleOffset(1.5);
		Linearity[i]->GetXaxis()->CenterTitle();
		Linearity[i]->GetYaxis()->CenterTitle();

	}

	cout << "Now Making Linearity" << endl;


	for(int i = InitFile; i < NFiles; i = i + StepSize){


		Energy = EnergyStep * i + 0.1;


		InfileName = Form("ToMerge/G4EICDetector.root_g4cemc_eval_%d_%.1f.root",FileIndex,Energy);


		fin = new TFile(InfileName.Data());
		fin->cd();


		TTree * t = (TTree * ) fin->Get("ntp_gshower");


		t->SetBranchAddress("e",&e);	
		t->SetBranchAddress("ge",&ge);	
		t->SetBranchAddress("eta",&eta);	

		int NEvents = 	t->GetEntries();

		for(int k = 0; k < NEvents; k++){

			t->GetEntry(k);

			//cout << "e = " << e << "   ge = " <<    ge  << "  eta = " << eta << endl;

			for(int l= 0; l < NetaCut; l++){

				if(eta > EtaCutDown[l] && eta <  EtaCutUp[l] )Linearity[l]->Fill(ge,e/ge);

			}
		}
	

			FileIndex = FileIndex + StepSize * 10;


	}
	
	FileIndex = 2222715;

	TLine *unity = new TLine(Emin,1,Emax,1);
	unity->SetLineStyle(2);
	unity->SetLineWidth(3);
	unity->SetLineColor(2);


	TLatex *lat4 = new TLatex();
	lat4->SetNDC();
	lat4->SetTextSize(0.04);
	lat4->SetTextColor(kBlack);


	for(int j= 0; j < NetaCut; j++){

		Linearity[j]->Draw("COLZ");
		unity->Draw("SAME");
		lat4->DrawLatex(0.40,0.2,Form("%.2f < |#eta| < %.2f",EtaCutDown[j],EtaCutUp[j]));

		c->SaveAs(Form("AnaPlots/Linearity/EnergyLinearity_%.2f_%.2f.png",EtaCutDown[j],EtaCutUp[j]));

	}

	cout << "DONE Making Linearity" << endl;

	for(int j = 0; j < NetaCut; j++){

		cout << "Now Working on Eta Range = " << EtaCutDown[j] <<  " - "  << EtaCutUp[j] << endl;
	//	Index = 0;

		for(int i = InitFile; i < NFiles; i = i + StepSize){

			Energy = EnergyStep * i + 0.1;

	

			double EnergyCent = Energy * ScaleFactor;

			double UpperE = Energy * ScaleFactor + Energy * Width;
			double LowerE = Energy * ScaleFactor - Energy * Width;


		//	InfileName = Form("Merged/G4EICDetector.root_g4cemc_eval_%.1f.root",Energy);
	
		//	InfileName = Form("ToMerge/G4EICDetector.root_g4cemc_eval_%.1f.root",Energy);


			InfileName = Form("ToMerge/G4EICDetector.root_g4cemc_eval_%d_%.1f.root",FileIndex,Energy);


			cout << "Now Working on i = " << i << endl;
			cout << "InfileName = " << InfileName.Data() << endl;



			fin = new TFile(InfileName.Data());
			fin->cd();


			TTree * t = (TTree * ) fin->Get("ntp_gshower");


			TH2D * EnergyReulation = new TH2D("EnergyReulation","",30,3.8,5.5,30,0.93,1.3);


			cout << "Pass 1 " << endl;
			/*
			   t->Project("EnergyReulation","(ge/e):e");

			   cout << "Pass 2 " << endl;

			   EnergyReulation->GetXaxis()->SetTitle("Reconstructed Electron Energy (GeV)");
			   EnergyReulation->GetYaxis()->SetTitle("Reco/Real Energy");
			   EnergyReulation->GetYaxis()->SetTitleOffset(1.3);

			   EnergyReulation->GetXaxis()->CenterTitle();
			   EnergyReulation->GetYaxis()->CenterTitle();



			   EnergyReulation->Draw("COLZ");



			   c->SaveAs(Form("AnaPlots/ErrorvsEnergy_%0.f.png",Energy));


*/

			TH1D * EnergyDis = new TH1D("EnergyDis","",70,LowerE,UpperE);


			t->Project("EnergyDis","e",Form("abs(eta) > %f && abs(eta) < %f",EtaCutDown[j],EtaCutUp[j]));


			cout << "Pass 3" << endl;

			EnergyDis->GetXaxis()->SetTitle("Reconstructed Electron Energy (GeV)");
			EnergyDis->GetYaxis()->SetTitle("Counts");
			EnergyDis->GetYaxis()->SetTitleOffset(1.3);

			EnergyDis->GetXaxis()->CenterTitle();
			EnergyDis->GetYaxis()->CenterTitle();
			TF1 * f = new TF1("f","gaus",LowerE,UpperE);




			/*	
				TF1 * f = new TF1("f","[0] * TMath::Exp(-(x-[1])*(x-[1])/(2 * [2]*[2]))",LowerE,UpperE);
				f->SetParLimits(1,EnergyCent-0.1,EnergyCent+0.1);
				f->SetParLimits(0,30,60);


				EnergyDis->Fit(f);
				*/
			EnergyDis->Fit("f");

			EnergyDis->Draw();

			TLatex *latE = new TLatex();
			latE->SetNDC();
			latE->SetTextSize(0.04);
			latE->SetTextColor(kBlack);
			latE->DrawLatex(0.57,0.8,Form("Energy = %.1f GeV",Energy));





			c->SaveAs(Form("AnaPlots/EnergyReso/EnergyReso_%d.png",i));


			TH2D * TransProfHis = new TH2D("TransProfHis","",200,-100,100,200,-100,100);
			TransProfHis->GetXaxis()->SetTitle("y (cm)");
			TransProfHis->GetYaxis()->SetTitle("x (cm)");
			TransProfHis->GetYaxis()->SetTitleOffset(1.4);
			TransProfHis->GetXaxis()->CenterTitle();
			TransProfHis->GetYaxis()->CenterTitle();

			t->Project("TransProfHis","y:x");

			TransProfHis->Draw("COLZ");

			latE->DrawLatex(0.27,0.6,Form("Energy = %.1f",Energy));

			c->SaveAs(Form("AnaPlots/TransProfile/TransProfile_%d.png",i));





			TH1D * LongProfHis = new TH1D("LongProfHis","",100,-10,10);
			LongProfHis->GetXaxis()->SetTitle("z (cm)");
			LongProfHis->GetYaxis()->SetTitle("Counts");
			LongProfHis->GetYaxis()->SetTitleOffset(1.4);
			LongProfHis->GetXaxis()->CenterTitle();
			LongProfHis->GetYaxis()->CenterTitle();

			t->Project("LongProfHis","z");

			LongProfHis->Draw("COLZ");

			latE->DrawLatex(0.27,0.6,Form("Energy = %.1f",Energy));

			c->SaveAs(Form("AnaPlots/LongProfile/LongProfile_%d.png",i));





			TH1D * EtaDisHis = new TH1D("EtaDisHis","",50,-0.15,0.15);
			EtaDisHis->GetXaxis()->SetTitle("#eta");
			EtaDisHis->GetYaxis()->SetTitle("Counts");
			EtaDisHis->GetYaxis()->SetTitleOffset(1.4);
			EtaDisHis->GetXaxis()->CenterTitle();
			EtaDisHis->GetYaxis()->CenterTitle();

			t->Project("EtaDisHis","eta");
			EtaDisHis->Draw("COLZ");

			latE->DrawLatex(0.27,0.85,Form("Energy = %.1f",Energy));

			c->SaveAs(Form("AnaPlots/EtaDisHis/EtaDisHis_%d.png",i));



			double EnergyWidth = f->GetParameter(2);
			double EnergyWidthError = f->GetParError(2);
			double EnergyMean = f->GetParameter(1);
			double EnergyMeanError = f->GetParError(1);

			double EnergyResolution = EnergyWidth/EnergyMean;

			double EnergyResolutionError = EnergyResolution *  sqrt(EnergyWidthError/EnergyWidth * EnergyWidthError/EnergyWidth + EnergyMeanError/EnergyMean *  EnergyMeanError/EnergyMean );







			cout << "Energy = " << Energy << "      Energy Resolution = " << EnergyResolution << "    Energy Resolution Error =   " << EnergyResolutionError << endl;

			int EnergyXBin = EnergyResolutionHis->GetXaxis()->FindBin(Energy);

			EnergyResolutionHis->SetBinContent(EnergyXBin,EnergyResolution);
			EnergyResolutionHis->SetBinError(EnergyXBin,EnergyResolutionError);


			TH1D * EnergyDisDep = new TH1D("EnergyDisDep","",70,LowerE*ScaleDownFactor,UpperE*ScaleDownFactor);


			t->Project("EnergyDisDep","gedep",Form("abs(eta) > %f && abs(eta) < %f",EtaCutDown[j],EtaCutUp[j]));



			EnergyDisDep->GetXaxis()->SetTitle("Deposited Electron Energy All Scintillators (GeV)");
			EnergyDisDep->GetYaxis()->SetTitle("Counts");
			EnergyDisDep->GetYaxis()->SetTitleOffset(1.3);

			EnergyDisDep->GetXaxis()->CenterTitle();
			EnergyDisDep->GetYaxis()->CenterTitle();
			
			TF1 * fdep = new TF1("fdep","gaus",LowerE*ScaleDownFactor,UpperE*ScaleDownFactor);
			EnergyDisDep->Fit("fdep");
			latE->DrawLatex(0.67,0.7,Form("Energy = %.1f",Energy));


			c->SaveAs(Form("AnaPlots/EnergyResoDep/EnergyReso_%d.png",i));

			double EnergyWidthDep = fdep->GetParameter(2);
			double EnergyWidthDepError = fdep->GetParError(2);
			double EnergyMeanDep = fdep->GetParameter(1);
			double EnergyMeanDepError = fdep->GetParError(1);

			double EnergyResolutionDep = EnergyWidthDep/EnergyMeanDep;

			double EnergyResolutionDepError = EnergyResolutionDep *  sqrt(EnergyWidthDepError/EnergyWidthDep * EnergyWidthDepError/EnergyWidthDep + EnergyMeanDepError/EnergyMeanDep *  EnergyMeanDepError/EnergyMeanDep );

			EnergyResolutionDepHis->SetBinContent(EnergyXBin,EnergyResolutionDep);
			EnergyResolutionDepHis->SetBinError(EnergyXBin,EnergyResolutionDepError);


			FileIndex = FileIndex + StepSize * 10;

		//	Index = Index + 1;
		}

		

		//	TF1 * f2 = new TF1("f2","[0] + [1]/sqrt(x) + [2]/x",2,Emax);
		TF1 * f2 = new TF1("f2","[0] + [1]/sqrt(x)",Emin,Emax);
		f2->SetLineColor(kRed);

		EnergyResolutionHis->Draw("ep");
		EnergyResolutionHis->Fit(f2,"R");


		double ConstTerm = f2->GetParameter(0);
		double ConstTermError = f2->GetParError(0);

		double StatTerm = f2->GetParameter(1);
		double StatTermError = f2->GetParError(1);


		//	double ExtraTerm = f2->GetParameter(2);
		//	double ExtraTermError = f2->GetParError(2);




		//	cout << "ConstTerm = " << ConstTerm << "    ConstTerm Error = " << ConstTermError << "  StatTerm =  " << StatTerm << "  StatTermError = " << StatTermError <<   "  ExtraTerm = "  << ExtraTerm << "  ExtraTermError =  "  <<  ExtraTermError  << endl;

		cout << "ConstTerm = " << ConstTerm << "    ConstTerm Error = " << ConstTermError << "  StatTerm =  " << StatTerm << "  StatTermError = " << StatTermError << endl;
		//	TString FuncName = Form("#Delta E/E = %.2f + %.2f/#sqrt{E} + %.4f/E",ConstTerm,StatTerm,ExtraTerm);
		TString FuncName = Form("#Delta E/E = %.3f + %.3f/#sqrt{E}",ConstTerm,StatTerm);


		TLatex *lat = new TLatex();
		lat->SetNDC();
		lat->SetTextSize(0.04);
		lat->SetTextColor(kRed);


		TF1 * f3 = new TF1("f3","sqrt([0] * [0] + [1] * [1]/x)",Emin,Emax);
		f3->SetLineColor(kBlue);

		EnergyResolutionHis->SetTitle(Form("Reconstructed Energy Resolution in |#eta| < %.2f",EtaCutUp[j]));
		EnergyResolutionHis->SetMaximum(0.10);

		EnergyResolutionHis->Fit(f3,"R");
		EnergyResolutionHis->Draw("ep");
		f2->Draw("SAME");
		f3->Draw("SAME");
		double ConstTerm2 = f3->GetParameter(0);
		double ConstTermError2 = f3->GetParError(0);

		double StatTerm2 = f3->GetParameter(1);
		double StatTermError2 = f3->GetParError(1);

		cout << "ConstTerm = " << ConstTerm2 << "    ConstTerm Error = " << ConstTermError2 << "  StatTerm =  " << StatTerm2 << "  StatTermError = " << StatTermError2 << endl;
		//	TString FuncName = Form("#Delta E/E = %.2f + %.2f/#sqrt{E} + %.4f/E",ConstTerm,StatTerm,ExtraTerm);
	//	TString FuncName2 = Form("#Delta E/E = #sqrt{(%.3f)^{2} + (%.3f/#sqrt{E})^{2}}",ConstTerm2,StatTerm2);

		TString FuncName2 = Form("#Delta E/E = %.3f/#sqrt{E} #oplus %.3f",StatTerm2,ConstTerm2);


		//lat->DrawLatex(0.27,0.7,FuncName.Data());

		TLatex *lat2 = new TLatex();
		lat2->SetNDC();
		lat2->SetTextSize(0.04);
		lat2->SetTextColor(kBlue);
//		lat2->DrawLatex(0.35,0.6,FuncName2.Data());
	
		

		
		TLatex *lat3 = new TLatex();
		lat3->SetNDC();
		lat3->SetTextSize(0.04);
		lat3->SetTextColor(kBlack);
	//	lat3->DrawLatex(0.27,0.8,Form("%.2f < |#eta| < %.2f",EtaCutDown[j],EtaCutUp[j]));

		TLegend* leg = new TLegend(0.30,0.65,0.75,0.90,NULL,"brNDC");
		leg->SetBorderSize(0);
		leg->SetTextSize(0.038);
		leg->SetTextFont(42);
		leg->SetFillStyle(0);

	    leg->AddEntry("","#it{#bf{Fun4All-EIC}} Simulation","");
//		leg->AddEntry("","eRD1 W-scintillator Shashlik EMCal","");
//		leg->AddEntry("",Form("Rapidity Range: |#eta| < %.1f",EtaCutUp[j]),"");
//		leg->AddEntry("","Indenting angle < 0.1 rad","");
		leg->AddEntry("","Geant4, digitization and clustering",""); 
		leg->AddEntry(EnergyResolutionHis,"Electron energy resolution","EP");
    	leg->AddEntry(f3, Form("#bf{#Delta E/E = %.3f/#sqrt{E} #oplus %.3f}",StatTerm2,ConstTerm2),"l");
		leg->Draw("SAME");



		c->SaveAs(Form("AnaPlots/Final/EnergyResolution_%.2f_%.2f.png",EtaCutDown[j],EtaCutUp[j]));

		StatTermEta->SetBinContent(j+1,StatTerm);
		StatTermEta->SetBinError(j+1,StatTermError);

		ConstTermEta->SetBinContent(j+1,ConstTerm);
		ConstTermEta->SetBinError(j+1,ConstTermError);

		cout << "EnergyResolutionHis Before = " << EnergyResolutionHis->Integral() << endl;

		EnergyResolutionHis->Reset();


		//Energy Dep//


		TF1 * f2dep = new TF1("f2dep","[0] + [1]/sqrt(x)",Emin,Emax);
		f2dep->SetLineColor(kRed);

		EnergyResolutionDepHis->Draw("ep");
		EnergyResolutionDepHis->Fit(f2dep,"R");


		double ConstTermDep = f2dep->GetParameter(0);
		double ConstTermDepError = f2dep->GetParError(0);

		double StatTermDep = f2dep->GetParameter(1);
		double StatTermDepError = f2dep->GetParError(1);


		//	double ExtraTerm = f2->GetParameter(2);
		//	double ExtraTermError = f2->GetParError(2);




		//	cout << "ConstTerm = " << ConstTerm << "    ConstTerm Error = " << ConstTermError << "  StatTerm =  " << StatTerm << "  StatTermError = " << StatTermError <<   "  ExtraTerm = "  << ExtraTerm << "  ExtraTermError =  "  <<  ExtraTermError  << endl;

		cout << "ConstTerm = " << ConstTermDep << "    ConstTerm Error = " << ConstTermDepError << "  StatTerm =  " << StatTermDep << "  StatTermError = " << StatTermDepError << endl;
		//	TString FuncName = Form("#Delta E/E = %.2f + %.2f/#sqrt{E} + %.4f/E",ConstTerm,StatTerm,ExtraTerm);
		TString FuncNameDep = Form("#Delta E/E = %.3f + %.3f/#sqrt{E}",ConstTermDep,StatTermDep);


		TLatex *latdep = new TLatex();
		latdep->SetNDC();
		latdep->SetTextSize(0.04);
		latdep->SetTextColor(kRed);


		TF1 * f3dep = new TF1("f3dep","sqrt([0] * [0] + [1] * [1]/x)",Emin,Emax);
		f3dep->SetLineColor(kBlue);


		EnergyResolutionDepHis->SetTitle(Form("Total Deposited Energy Resolution in |#eta| < %.2f",EtaCutUp[j]));
		EnergyResolutionDepHis->SetMaximum(0.10);

		EnergyResolutionDepHis->Fit(f3dep,"R");
		EnergyResolutionDepHis->Draw("ep");
//		f2dep->Draw("SAME");
		f3dep->Draw("SAME");
		double ConstTerm2Dep = f3dep->GetParameter(0);
		double ConstTermDepError2 = f3dep->GetParError(0);

		double StatTerm2Dep = f3dep->GetParameter(1);
		double StatTermDepError2 = f3dep->GetParError(1);

		cout << "ConstTerm = " << ConstTerm2Dep << "    ConstTerm Error = " << ConstTermDepError2 << "  StatTerm =  " << StatTerm2Dep << "  StatTermError = " << StatTermDepError2 << endl;
		//	TString FuncName = Form("#Delta E/E = %.2f + %.2f/#sqrt{E} + %.4f/E",ConstTerm,StatTerm,ExtraTerm);
//		TString FuncName2Dep = Form("#Delta E/E = #sqrt{(%.3f)^{2} + (%.3f/#sqrt{E})^{2}}",ConstTerm2Dep,StatTerm2Dep);
		TString FuncName2Dep = Form("#Delta E/E = %.3f/#sqrt{E} #oplus  %.3f",StatTerm2Dep,ConstTerm2Dep);


//		latdep->DrawLatex(0.27,0.7,FuncNameDep.Data());

		TLatex *lat2dep = new TLatex();
		lat2dep->SetNDC();
		lat2dep->SetTextSize(0.04);
		lat2dep->SetTextColor(kBlue);
//		lat2dep->DrawLatex(0.35,0.6,FuncName2Dep.Data());


		TLatex *lat3dep = new TLatex();
		lat3dep->SetNDC();
		lat3dep->SetTextSize(0.04);
		lat3dep->SetTextColor(kBlack);
//		lat3dep->DrawLatex(0.27,0.8,Form("%.2f < |#eta| < %.2f",EtaCutDown[j],EtaCutUp[j]));



		TLegend* legdep = new TLegend(0.30,0.65,0.75,0.90,NULL,"brNDC");
		legdep->SetBorderSize(0);
		legdep->SetTextSize(0.038);
		legdep->SetTextFont(42);
		legdep->SetFillStyle(0);

	    legdep->AddEntry("","#it{#bf{Fun4All-EIC}} Simulation","");
//    legdep->AddEntry("","eRD1 W-scintillator Shashlik EMCal","");
//		legdep->AddEntry("",Form("Rapidity Range: |#eta| < %.1f",EtaCutUp[j]),"");
//	    legdep->AddEntry("","Indenting angle < 0.1 rad","");
		legdep->AddEntry("","Geant4, truth energy deposition",""); 
		legdep->AddEntry(EnergyResolutionHis,"Electron energy resolution","EP");
    	legdep->AddEntry(f3, Form("#bf{#Delta E/E = %.3f/#sqrt{E} #oplus %.3f}",StatTerm2Dep,ConstTerm2Dep),"l");
		legdep->Draw("SAME");



		c->SaveAs(Form("AnaPlots/FinalDep/EnergyResolutionDep_%.2f_%.2f.png",EtaCutDown[j],EtaCutUp[j]));



		StatTermEta->SetBinContent(j+1,StatTerm);
		StatTermEta->SetBinError(j+1,StatTermError);

		ConstTermEta->SetBinContent(j+1,ConstTerm);
		ConstTermEta->SetBinError(j+1,ConstTermError);

		cout << "EnergyResolutionHis Before = " << EnergyResolutionHis->Integral() << endl;

		EnergyResolutionHis->Reset();

		cout << "EnergyResolutionHis After = " << EnergyResolutionHis->Integral() << endl;






	}


	StatTermEta->Draw("ep");
	c->SaveAs("AnaPlots/Final/EtaStatPlot.png");


	ConstTermEta->Draw("ep");
	c->SaveAs("AnaPlots/Final/EtaConstPlot.png");
	

	

}
