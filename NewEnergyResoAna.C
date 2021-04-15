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



void NewEnergyResoAna(int ConfigType, int ThickType, int DoWhat){



	TString FolderName;



	double MaxHis = 0.10;

	if(DoWhat == 0) MaxHis = 0.10;
	if(DoWhat == 1) MaxHis = 0.20;
	if(DoWhat == 2) MaxHis = 0.17;



	int Color[5] = {1,2,3,4,6};
	TString ConfigName[5] = {"AQuater","AHalf","3Quarters","Unity","AFifth"};


	int Type = ConfigType + ThickType * 10;
	
	int TotalRadLength;
	
	if(ThickType == 0) TotalRadLength = 20;
	if(ThickType == 1) TotalRadLength = 18;
	if(ThickType == 2) TotalRadLength = 22;
	if(ThickType == 3) TotalRadLength = 24;


	if(Type == 0) FolderName =  "AQuarter/20/";
	if(Type == 1) FolderName =  "AHalf/20/";
	if(Type == 2) FolderName =  "3Quarters/20/";
	if(Type == 3) FolderName =  "Unity/20/";
	if(Type == 4) FolderName =  "AFifth/20/";


	if(Type == 10) FolderName =  "AQuarter/18/";
	if(Type == 11) FolderName =  "AHalf/18/";
	if(Type == 12) FolderName =  "3Quarters/18/";
	if(Type == 13) FolderName =  "Unity/18/";
	


	if(Type == 20) FolderName =  "AQuarter/22/";
	if(Type == 21) FolderName =  "AHalf/22/";
	if(Type == 22) FolderName =  "3Quarters/22/";
	if(Type == 23) FolderName =  "Unity/22/";


	if(Type == 30) FolderName =  "AQuarter/24/";
	if(Type == 31) FolderName =  "AHalf/24/";
	if(Type == 32) FolderName =  "3Quarters/24/";
	if(Type == 33) FolderName =  "Unity/24/";



	const int NGood = 36;
	double ListofGoodEnergies[NGood] = {0.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0,10.5,11.0,11.5,12.0,12.5,13.0,13.5,14.0,14.5,15.0,15.5,16.5,17.0,17.5,18.5,20.5,22.5,25.0};
	int Good = 0;

	gStyle->SetPadLeftMargin(0.15);
	gStyle->SetPadBottomMargin(0.13);
	gStyle->SetOptStat(0);

	

	SetsPhenixStyle();
	

	double ScaleDownFactor = 0.13;


	double EtaMax = 0.10;
	double EtaMin = 0.0;

//	const int NetaCut = 5;
	const int NetaCut = 1;

	double EtaStep = (EtaMax - EtaMin)/NetaCut;

	double EtaCutUp[NetaCut];
	double EtaCutDown[NetaCut];

	for(int i = 0; i < NetaCut; i++){
		EtaCutDown[i] = EtaMin + EtaStep * i;
		EtaCutUp[i] = EtaMin + EtaStep * (i+1);

	}




	const int NType = 5;


	double ScaleFactor[NType]; 
	double ScaleFactorConst[NType] = {1.2,0.8,0.5,0.4,1.4};
	double Width[NType] = {0.2,0.2,0.2,0.2,0.2};
	double Coeff[NType] = {0.1,0.1,0.1,0.1,0.1};
	double Rate[NType] = {0,1,1,1,0.1};

	int InitFile = 1;
	int NFiles = 51;

	int Index = 0;

	int StepSize = 1;
	
	//int FileIndex = 2222715;

//	int TotalFiles = (NFiles - InitFile)/StepSize;
	int TotalFiles = (NFiles - InitFile);

	double EnergyStep = 0.50; 
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
	//	cout << "Now Working on i = " << i << endl;

	//	Energy = EnergyStep * i + 0.1;

		Energy = EnergyStep * i + 0.0;

		if(Type == 4 && (Energy == 15 || Energy == 15.5 || Energy == 17)) continue; 	
		if(Type == 10 && Energy == 19.0) continue; 
	//	if(Type == 11 && Energy > 3) continue; 
	
		/*
		if(Type == 30){
				Good = 0;
				for(int s = 0;  s < NGood; s++){

					if(Energy == ListofGoodEnergies[s])  Good = 1;
				}

				if(Good == 0) continue;

		}
		*/	
		cout << "Now Working on Energy = " << Energy << endl;


//		InfileName = Form("ToMerge/G4EICDetector.root_g4cemc_eval_%d_%.1f.root",FileIndex,Energy);

		InfileName = Form("Merged/%s/G4EICDetector.root_g4cemc_eval_%.1f.root",FolderName.Data(),Energy);

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
	

	//		FileIndex = FileIndex + StepSize * 10;


	}
	
	//FileIndex = 2222715;

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

		c->SaveAs(Form("AnaPlots/%s/Linearity/EnergyLinearity_%.2f_%.2f.png",FolderName.Data(),EtaCutDown[j],EtaCutUp[j]));

	}

	cout << "DONE Making Linearity" << endl;

	for(int j = 0; j < NetaCut; j++){

		cout << "Now Working on Eta Range = " << EtaCutDown[j] <<  " - "  << EtaCutUp[j] << endl;
	//	Index = 0;

		for(int i = InitFile; i < NFiles; i = i + StepSize){
	
	
	//		Energy = EnergyStep * i + 0.1;

			Energy = EnergyStep * i + 0.0;

			if(Type == 4 && (Energy == 15 || Energy == 15.5 || Energy == 17)) continue; 	
			if(Type == 10 && Energy == 19.0) continue; 
			//if(Type == 11 && Energy > 3) continue; 
			/*
			if(Type == 30){
				Good = 0;
				for(int s = 0;  s < NGood; s++){

					if(Energy == ListofGoodEnergies[s])  Good = 1;
				}

				if(Good == 0) continue;

			}
			*/
			ScaleFactor[ConfigType] =  ScaleFactorConst[ConfigType] - Rate[ConfigType] * Energy * 0.01;


			double EnergyCent = Energy * ScaleFactor[ConfigType];

			Width[ConfigType] = 0.2 + Coeff[ConfigType]/sqrt(Energy);



			double UpperE = Energy * ScaleFactor[ConfigType] + Energy * Width[ConfigType];
			double LowerE = Energy * ScaleFactor[ConfigType] - Energy * Width[ConfigType];

			InfileName = Form("Merged/%s/G4EICDetector.root_g4cemc_eval_%.1f.root",FolderName.Data(),Energy);

		//	InfileName = Form("ToMerge/G4EICDetector.root_g4cemc_eval_%.1f.root",Energy);


	//		InfileName = Form("ToMerge/G4EICDetector.root_g4cemc_eval_%d_%.1f.root",FileIndex,Energy);


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



			   c->SaveAs(Form("AnaPlots/%s/ErrorvsEnergy_%0.f.png",Energy));


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
			f->SetLineColor(kRed);



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





			c->SaveAs(Form("AnaPlots/%s/EnergyReso/EnergyReso_%d.png",FolderName.Data(),i));
			if(Energy == 20.5) 		c->SaveAs(Form("CollectAllResults/EnergyReso_%s_%d_%d.png",ConfigName[ConfigType].Data(),ConfigType,ThickType));

			TH2D * TransProfHis = new TH2D("TransProfHis","",200,-100,100,200,-100,100);
			TransProfHis->GetXaxis()->SetTitle("y (cm)");
			TransProfHis->GetYaxis()->SetTitle("x (cm)");
			TransProfHis->GetYaxis()->SetTitleOffset(1.4);
			TransProfHis->GetXaxis()->CenterTitle();
			TransProfHis->GetYaxis()->CenterTitle();

			t->Project("TransProfHis","y:x");

			TransProfHis->Draw("COLZ");

			latE->DrawLatex(0.27,0.6,Form("Energy = %.1f",Energy));

			c->SaveAs(Form("AnaPlots/%s/TransProfile/TransProfile_%d.png",FolderName.Data(),i));





			TH1D * LongProfHis = new TH1D("LongProfHis","",100,-10,10);
			LongProfHis->GetXaxis()->SetTitle("z (cm)");
			LongProfHis->GetYaxis()->SetTitle("Counts");
			LongProfHis->GetYaxis()->SetTitleOffset(1.4);
			LongProfHis->GetXaxis()->CenterTitle();
			LongProfHis->GetYaxis()->CenterTitle();

			t->Project("LongProfHis","z");

			LongProfHis->Draw("COLZ");

			latE->DrawLatex(0.27,0.6,Form("Energy = %.1f",Energy));

			c->SaveAs(Form("AnaPlots/%s/LongProfile/LongProfile_%d.png",FolderName.Data(),i));





			TH1D * EtaDisHis = new TH1D("EtaDisHis","",50,-0.15,0.15);
			EtaDisHis->GetXaxis()->SetTitle("#eta");
			EtaDisHis->GetYaxis()->SetTitle("Counts");
			EtaDisHis->GetYaxis()->SetTitleOffset(1.4);
			EtaDisHis->GetXaxis()->CenterTitle();
			EtaDisHis->GetYaxis()->CenterTitle();

			t->Project("EtaDisHis","eta");
			EtaDisHis->Draw("COLZ");

			latE->DrawLatex(0.27,0.85,Form("Energy = %.1f",Energy));

			c->SaveAs(Form("AnaPlots/%s/EtaDisHis/EtaDisHis_%d.png",FolderName.Data(),i));
		



			double EnergyWidth = f->GetParameter(2);
			double EnergyWidthError = f->GetParError(2);
			double EnergyMean = f->GetParameter(1);
			double EnergyMeanError = f->GetParError(1);

			double EnergyResolution = EnergyWidth/EnergyMean;

			double EnergyResolutionError = EnergyResolution *  sqrt(EnergyWidthError/EnergyWidth * EnergyWidthError/EnergyWidth + EnergyMeanError/EnergyMean *  EnergyMeanError/EnergyMean );








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
			fdep->SetLineColor(kRed);
			EnergyDisDep->Fit(fdep,"R");
			EnergyDisDep->Draw();
			latE->DrawLatex(0.67,0.7,Form("Energy = %.1f",Energy));
			fdep->Draw("SAME");

			c->SaveAs(Form("AnaPlots/%s/EnergyResoDep/EnergyReso_%d.png",FolderName.Data(),i));
			
			if(Energy == 20.5) 		c->SaveAs(Form("CollectAllResults/EnergyResoDep_%s_%d_%d.png",ConfigName[ConfigType].Data(),ConfigType,ThickType));


			double EnergyWidthDep = fdep->GetParameter(2);
			double EnergyWidthDepError = fdep->GetParError(2);
			double EnergyMeanDep = fdep->GetParameter(1);
			double EnergyMeanDepError = fdep->GetParError(1);

			double EnergyResolutionDep = EnergyWidthDep/EnergyMeanDep;

			double EnergyResolutionDepError = EnergyResolutionDep *  sqrt(EnergyWidthDepError/EnergyWidthDep * EnergyWidthDepError/EnergyWidthDep + EnergyMeanDepError/EnergyMeanDep *  EnergyMeanDepError/EnergyMeanDep );


			cout << "Energy = " << Energy << "      Energy Resolution = " << EnergyResolution << "  Deposited Energy Resolution =   " << EnergyResolutionDep << endl;
	

			EnergyResolutionDepHis->SetBinContent(EnergyXBin,EnergyResolutionDep);
			EnergyResolutionDepHis->SetBinError(EnergyXBin,EnergyResolutionDepError);




		//	FileIndex = FileIndex + StepSize * 10;

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
		if(DoWhat == 0) f3->SetLineColor(kBlue);
		if(DoWhat == 1) f3->SetLineColor(Color[ConfigType]);	
		if(DoWhat == 2) f3->SetLineColor(Color[ThickType]);	

		EnergyResolutionHis->SetTitle(Form("Reconstructed Energy Resolution in |#eta| < %.2f",EtaCutUp[j]));
		EnergyResolutionHis->SetMaximum(MaxHis);

		EnergyResolutionHis->Fit(f3,"R");
		EnergyResolutionHis->Draw("ep");
	//	f2->Draw("SAME");
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

		TLegend* leg = new TLegend(0.30,0.55,0.75,0.90,NULL,"brNDC");
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
		leg->AddEntry("",Form("Total Radiation Length = %d",TotalRadLength),""); 	
    	leg->AddEntry(f3, Form("#bf{#Delta E/E = %.4f/#sqrt{E} #oplus %.4f}",StatTerm2,ConstTerm2),"l");
		leg->Draw("SAME");



		c->SaveAs(Form("AnaPlots/%s/FinalPlots/EnergyResolution_%.2f_%.2f.png",FolderName.Data(),EtaCutDown[j],EtaCutUp[j]));
		c->SaveAs(Form("CollectAllResults/EnergyResolution_%s_%d_%d.png",ConfigName[ConfigType].Data(),ConfigType,ThickType));

		StatTermEta->SetBinContent(j+1,StatTerm);
		StatTermEta->SetBinError(j+1,StatTermError);

		ConstTermEta->SetBinContent(j+1,ConstTerm);
		ConstTermEta->SetBinError(j+1,ConstTermError);

		cout << "EnergyResolutionHis Before = " << EnergyResolutionHis->Integral() << endl;

		EnergyResolutionHis->Reset();


		//Energy Dep//

	
		TF1 * f2dep = new TF1("f2dep","[0] + [1]/sqrt(x)",Emin,Emax);
		f2dep->SetLineColor(kRed);

		if(DoWhat == 2) EnergyResolutionDepHis->SetMinimum(0);
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
		if(DoWhat == 0) f3dep->SetLineColor(kBlue);
		if(DoWhat == 1)  f3dep->SetLineColor(Color[ConfigType]);	
		if(DoWhat == 2)  f3dep->SetLineColor(Color[ThickType]);	

		EnergyResolutionDepHis->SetTitle(Form("Total Deposited Energy Resolution in |#eta| < %.2f",EtaCutUp[j]));
		EnergyResolutionDepHis->SetMaximum(MaxHis);

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
		TString FuncName2Dep = Form("#Delta E/E = %.4f/#sqrt{E} #oplus  %.4f",StatTerm2Dep,ConstTerm2Dep);


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



		TLegend* legdep = new TLegend(0.30,0.55,0.75,0.90,NULL,"brNDC");
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
		legdep->AddEntry("",Form("Total Radiation Length = %d",TotalRadLength),""); 	
    	legdep->AddEntry(f3dep, Form("#bf{#Delta E/E = %.4f/#sqrt{E} #oplus %.4f}",StatTerm2Dep,ConstTerm2Dep),"l");
		legdep->Draw("SAME");



		c->SaveAs(Form("AnaPlots/%s/FinalPlots/EnergyResolutionDep_%.2f_%.2f.png",FolderName.Data(),EtaCutDown[j],EtaCutUp[j]));
		c->SaveAs(Form("CollectAllResults/EnergyResolutionDep_%s_%d_%d.png",ConfigName[ConfigType].Data(),ConfigType,ThickType));



		StatTermEta->SetBinContent(j+1,StatTerm);
		StatTermEta->SetBinError(j+1,StatTermError);

		ConstTermEta->SetBinContent(j+1,ConstTerm);
		ConstTermEta->SetBinError(j+1,ConstTermError);

		cout << "EnergyResolutionHis Before = " << EnergyResolutionHis->Integral() << endl;

		EnergyResolutionHis->Reset();

		cout << "EnergyResolutionHis After = " << EnergyResolutionHis->Integral() << endl;
	





	}


	StatTermEta->Draw("ep");
//	c->SaveAs(Form("AnaPlots/%s/Final/EtaStatPlot.png",FolderName.Data()));


	ConstTermEta->Draw("ep");
//	c->SaveAs(Form("AnaPlots/%s/Final/EtaConstPlot.png",FolderName.Data()));
	
	TString outfilename = Form("OutFiles/EnResoRad_%d_%d.root",ConfigType,ThickType);
	
	TFile * fout = new TFile(outfilename.Data(),"RECREATE");
	fout->cd();

	EnergyResolutionDepHis->Write();

	fout->Close();





}
