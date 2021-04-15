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



void PlotStudy(int AbsOpt, int RadLengthOpt){



	int RadLength;

	if(RadLengthOpt == 0) RadLength = 18;
	if(RadLengthOpt == 1) RadLength = 20;
	if(RadLengthOpt == 2) RadLength = 22;

	double Norm;
	gStyle->SetOptStat(0);

	TString BeamName = "Electron";

	TString AbsMaterial;

	//	if(AbsOpt == 0) AbsMaterial = "W";
	if(AbsOpt == 0) AbsMaterial = "WCu";
	if(AbsOpt == 1) AbsMaterial = "Pb";
	if(AbsOpt == 2) AbsMaterial = "Li";
	if(AbsOpt == 3) AbsMaterial = "Al";

	float MoliereR;
	if(AbsOpt == 0) MoliereR = 2.5;
	if(AbsOpt == 1) MoliereR = 4.5;

	int ColorFunc;

	if(AbsOpt == 0 && RadLengthOpt == 0) ColorFunc = 4;
	if(AbsOpt == 0 && RadLengthOpt == 1) ColorFunc = 1;
	if(AbsOpt == 0 && RadLengthOpt == 2) ColorFunc = 6;
	if(AbsOpt == 1 && RadLengthOpt == 0) ColorFunc = 3;
	if(AbsOpt == 1 && RadLengthOpt == 1) ColorFunc = 2;


	double MinZ= 999;

	if(RadLengthOpt == 0) MinDepthZ = 301.1;
	if(RadLengthOpt == 1 && AbsOpt == 0) MinDepthZ = 300;
	if(RadLengthOpt == 2) MinDepthZ = 299;

	if(RadLengthOpt == 1 && AbsOpt == 1) MinDepthZ = 298.55;
	if(RadLengthOpt == 1 && AbsOpt == 2) MinDepthZ = 298.55 + 1550/2;
	if(RadLengthOpt == 1 && AbsOpt == 3) MinDepthZ = 306;


	TCanvas * cLog = new TCanvas("cLog","cLog",600,600);
	cLog->SetLogy();


	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();



	TCanvas * cBoth = new TCanvas("cBoth","cBoth",1200,600);
	cBoth->Divide(2,1);



	TCanvas * c4 = new TCanvas("c4","c4",600,600);
	c4->SetLogy();






	TCanvas * c4Long = new TCanvas("c4Long","c4Long",600,600);



	double ERecoMin = 1;
	double ERecoMax= 23;


	ScinThick = 0.15;



	TLatex *lat = new TLatex();
	lat->SetNDC();
	lat->SetTextSize(0.035);

	double SamplingFrac = 0.10;
	double SamplingFracDown = 0.06;	
	double SamplingFracUp = 0.14;	

	double EFitMinRatio = 0.06;	
	double EFitMaxRatio = 0.14;	

	double TotalEvents;

	double Emin;
	double Emax;

	double VtxX = MinDepthX;
	double VtxY = MinDepthY;

	double Energy;
	double EnergyStep = 1.0;



	double EFitMin;
	double EFitMax;
	

	 

	if(AbsOpt == 1) WRadLength = 0.41;
	if(AbsOpt == 1) WRadLength = 0.53;
	
	const int NFiles = 23;

	TString infile;
	TFile * fin;
	TTree * ShowerInfo;


	TH1D * ScinEDepHis[NFiles];


	TH1D * TotalEAllHis[NFiles]; 


	TH1D * EDepLongWAbsScale[NFiles]; 


	TH1D * EDepTransWAbsScale[NFiles]; 
	TH1D * EDepTransCum[NFiles]; 

	TH2D * EDepXY[NFiles]; 

	TH1D * EDepTransDiff[NFiles]; 
	TH1D * EDepTransDiffCM[NFiles]; 

	TH1D * EDepTransInc[NFiles]; 


	TH1D * EDepTransWAbsScaleNorm[NFiles]; 


	TH1D * EnergyReso = new TH1D("EnergyReso","",100,0,NFiles);
	EnergyReso->GetXaxis()->SetTitle("Incident Electron Energy (GeV)");
	EnergyReso->GetYaxis()->SetTitle("Energy Resolution");
	EnergyReso->GetXaxis()->CenterTitle();
	EnergyReso->GetYaxis()->CenterTitle();
	EnergyReso->GetYaxis()->SetTitleOffset(1.4);
	EnergyReso->SetLineColor(kBlack);
	EnergyReso->SetMarkerStyle(20);
	EnergyReso->SetMarkerSize(1.0);
	EnergyReso->SetMarkerColor(kBlack);



	TH1D * LeakEDepHis = new TH1D("LeakEDepHis","",ERecoMax-ERecoMin+1,ERecoMin-0.5,ERecoMax+0.5);
	LeakEDepHis->GetXaxis()->SetTitle("Incident Electron Energy (GeV)");
	LeakEDepHis->GetYaxis()->SetTitle("Energy Leakage Fraction (%)");
	LeakEDepHis->GetXaxis()->CenterTitle();
	LeakEDepHis->GetYaxis()->CenterTitle();
	LeakEDepHis->GetYaxis()->SetTitleOffset(1.4);
	LeakEDepHis->SetLineColor(kBlack);

	double StatTerm;
	double ConstTerm;


	double TotalEScindep;
	double TotalLeakEnergy;
	double LeakPercent = 0;

	double GausMean;
	double GausRMS;

	double GausMeanErr;
	double GausRMSErr;

	double EReso;
	double EResoErr;


	float FEMCALLYDep;
	int XResoBin;

	const int NSpecialE = 4;
	int SpecialE[NSpecialE] = {20,10,5,1};
	int SpecialColor[NSpecialE] = {1,2,3,4};
	


	//PHENIX Pb Shashlik EMCAL Configuration 
	if(RadLengthOpt == 0 && AbsOpt == 1){
		ScinThick = 0.4;
		MinDepthZ = 289.84;
	    SamplingFracDown = 0.15;	
		SamplingFracUp = 0.25;	
		EFitMinRatio = 0.15;
		EFitMaxRatio = 0.25;

	}
	


	//Declare Histograms//

	for(int i = 0; i < NFiles; i++){
	


		Energy = (i + 1) * EnergyStep;

		Emin = Energy * SamplingFracDown;
		Emax = Energy * SamplingFracUp;

		EFitMin = Energy * EFitMinRatio;
		EFitMax = Energy * EFitMaxRatio;


		EDepLongWAbsScale[i] = new TH1D(Form("EDepLongWAbsScale_%d",i),"",50,0,30);
		EDepLongWAbsScale[i]->GetXaxis()->SetTitle("EMCAL Depth (cm)");
		EDepLongWAbsScale[i]->GetYaxis()->SetTitle("dE/dz (GeV/cm)");
		EDepLongWAbsScale[i]->SetTitle("Incident e^{-} Longitudinal Energy Deposition");
		EDepLongWAbsScale[i]->GetXaxis()->CenterTitle();
		EDepLongWAbsScale[i]->GetYaxis()->CenterTitle();
		EDepLongWAbsScale[i]->GetYaxis()->SetTitleOffset(1.4);
		EDepLongWAbsScale[i]->SetLineColor(kBlack);



		EDepTransWAbsScale[i] = new TH1D(Form("EDepTransWAbsScale_%d",i),"",100,0,10);
		EDepTransWAbsScale[i]->GetXaxis()->SetTitle("EMCAL Radius (cm)");
		EDepTransWAbsScale[i]->GetYaxis()->SetTitle("dE/dR (GeV/cm)");
		EDepTransWAbsScale[i]->SetTitle("Incident e^{-} Transverse Energy Deposition");
		EDepTransWAbsScale[i]->GetXaxis()->CenterTitle();
		EDepTransWAbsScale[i]->GetYaxis()->CenterTitle();
		EDepTransWAbsScale[i]->GetYaxis()->SetTitleOffset(1.4);
		EDepTransWAbsScale[i]->SetLineColor(kBlack);



		EDepTransCum[i] = new TH1D(Form("EDepTransCum_%d",i),"",100,0,20);
		EDepTransCum[i]->GetXaxis()->SetTitle("EMCAL Radius (cm)");
		EDepTransCum[i]->GetYaxis()->SetTitle("Cumulative Transverse Energy Profile (% of E_{inc})");
		EDepTransCum[i]->SetTitle("Incident e^{-} Transverse Energy Deposition");
		EDepTransCum[i]->GetXaxis()->CenterTitle();
		EDepTransCum[i]->GetYaxis()->CenterTitle();
		EDepTransCum[i]->GetYaxis()->SetTitleOffset(1.4);
		EDepTransCum[i]->SetLineColor(kBlack);



		EDepTransDiff[i] = new TH1D(Form("EDepTransDiff_%d",i),"",100,0,10);
		EDepTransDiff[i]->GetXaxis()->SetTitle("EMCAL Radius/R_{M}");
		EDepTransDiff[i]->GetYaxis()->SetTitle("Differential Transverse Energy Profile (% of E_{inc})");
		EDepTransDiff[i]->SetTitle("Incident e^{-} Transverse Energy Deposition");
		EDepTransDiff[i]->GetXaxis()->CenterTitle();
		EDepTransDiff[i]->GetYaxis()->CenterTitle();
		EDepTransDiff[i]->GetYaxis()->SetTitleOffset(1.35);
		EDepTransDiff[i]->SetLineColor(kBlack);



		EDepTransDiffCM[i] = new TH1D(Form("EDepTransDiffCM_%d",i),"",100,0,20);
		EDepTransDiffCM[i]->GetXaxis()->SetTitle("EMCAL Radius (cm)");
		EDepTransDiffCM[i]->GetYaxis()->SetTitle("Differential Transverse Energy Profile (% of E_{inc})");
		EDepTransDiffCM[i]->SetTitle("Incident e^{-} Transverse Energy Deposition");
		EDepTransDiffCM[i]->GetXaxis()->CenterTitle();
		EDepTransDiffCM[i]->GetYaxis()->CenterTitle();
		EDepTransDiffCM[i]->GetYaxis()->SetTitleOffset(1.35);
		EDepTransDiffCM[i]->SetLineColor(kBlack);




		EDepTransInc[i] = new TH1D(Form("EDepTransInc_%d",i),"",100,0,10);
		EDepTransInc[i]->GetXaxis()->SetTitle("EMCAL Radius/R_{M}");
		EDepTransInc[i]->GetYaxis()->SetTitle("Cumulative Transverse Energy Profile (% of E_{inc})");
		EDepTransInc[i]->SetTitle("Incident e^{-} Transverse Energy Deposition");
		EDepTransInc[i]->GetXaxis()->CenterTitle();
		EDepTransInc[i]->GetYaxis()->CenterTitle();
		EDepTransInc[i]->GetYaxis()->SetTitleOffset(1.4);
		EDepTransInc[i]->SetLineColor(kBlack);


		EDepXY[i] = new TH2D(Form("EDepXY_%d",i),"",100,-2,2,100,-2,2);
		EDepXY[i]->GetXaxis()->SetTitle("XCorr (cm)");
		EDepXY[i]->GetYaxis()->SetTitle("YCorr (cm)");
		EDepXY[i]->SetTitle("Incident e^{-} Transverse Energy Deposition 2D Map");
		EDepXY[i]->GetXaxis()->CenterTitle();
		EDepXY[i]->GetYaxis()->CenterTitle();
		EDepXY[i]->GetYaxis()->SetTitleOffset(1.4);




		EDepTransWAbsScaleNorm[i] = new TH1D(Form("EDepTransWAbsScaleNorm_%d",i),"",200,0,20);
		EDepTransWAbsScaleNorm[i]->GetXaxis()->SetTitle("EMCAL Radius (cm)");
		EDepTransWAbsScaleNorm[i]->GetYaxis()->SetTitle("dE/dR Normalized to Unity");
		EDepTransWAbsScaleNorm[i]->SetTitle("Incident e^{-} Transverse Energy Deposition");
		EDepTransWAbsScaleNorm[i]->GetXaxis()->CenterTitle();
		EDepTransWAbsScaleNorm[i]->GetYaxis()->CenterTitle();
		EDepTransWAbsScaleNorm[i]->GetYaxis()->SetTitleOffset(1.4);
		EDepTransWAbsScaleNorm[i]->SetLineColor(kBlack);


		ScinEDepHis[i] = new TH1D(Form("ScinEDepHis_%d",i),"",100,Emin,Emax);
		ScinEDepHis[i]->GetXaxis()->SetTitle("Visible Energy (GeV)");
		ScinEDepHis[i]->GetYaxis()->SetTitle("Counts");
		ScinEDepHis[i]->SetTitle("e^{-} Transverse Energy Deposition in Scintillator");
		ScinEDepHis[i]->GetXaxis()->CenterTitle();
		ScinEDepHis[i]->GetYaxis()->CenterTitle();
		ScinEDepHis[i]->GetYaxis()->SetTitleOffset(1.4);
		ScinEDepHis[i]->SetMarkerStyle(20);
		ScinEDepHis[i]->SetMarkerSize(1.0);
		ScinEDepHis[i]->SetMarkerColor(kBlack);


		TotalEAllHis[i] = new TH1D(Form("TotalEAllHis_%d",i),"",100,Energy-0.1,Energy+0.1);
		TotalEAllHis[i]->GetXaxis()->SetTitle("Total Energy from EMCAL + Black Hole (GeV)");
		TotalEAllHis[i]->GetYaxis()->SetTitle("Number of Events");
		TotalEAllHis[i]->GetXaxis()->CenterTitle();
		TotalEAllHis[i]->GetYaxis()->CenterTitle();
		TotalEAllHis[i]->GetYaxis()->SetTitleOffset(1.4);




	}
	double AbsWidth = EDepTransWAbsScale[0]->GetBinWidth(1);


	for(int q = 1; q < NFiles+1; q++){

		Energy = q * EnergyStep;

		Emin = Energy * SamplingFracDown;
		Emax = Energy * SamplingFracUp;

		EFitMin = Energy * EFitMinRatio;
		EFitMax = Energy * EFitMaxRatio;




		infile	= Form("/sphenix/user/zshi/SimulationWork/EICShashlik/macros/macros/g4simulations/JobSub/SimulatedFiles/%s/%d/%d/ShowerInfo.root",AbsMaterial.Data(),RadLength,q);
	    //infile	= Form("/sphenix/user/zshi/SimulationWork/EICShashlik/macros/macros/g4simulations/JobSub/SimulatedFiles/%s/%dOld/%d/ShowerInfo.root",AbsMaterial.Data(),RadLength,q);

		fin = new TFile(infile.Data());
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
		ShowerInfo->SetBranchAddress("layer",&layer);
		ShowerInfo->SetBranchAddress("FEMCALedep",&FEMCALedep);
		ShowerInfo->SetBranchAddress("FEMCALLYDep",&FEMCALLYDep);


		int NEntries = ShowerInfo->GetEntries();
		int evtidPre = 0;


		double TotalE = 0;

		double ZPos;
		double ZPosCorr;
		
		double TotalEAll = 0;

		double XPos;
		double YPos;
		double RPos;
		double XBin;
		double YBin;

		double XPosCorr;
		double YPosCorr;
		double RPosCorr;

		TotalEScindep = 0;

		TotalEvents = 0;
		
		TotalLeakEnergy = 0;
		LeakPercent = 0;
	


		for(int i = 0; i < NEntries; i++){ 
		
			

			ShowerInfo->GetEntry(i);

			XPos = (X1 + X2)/2;
			YPos = (Y1 + Y2)/2;
			XPosCorr = XPos - MinDepthX;
			YPosCorr = YPos - MinDepthY;
			RPosCorr = sqrt(XPosCorr * XPosCorr + YPosCorr * YPosCorr);
			EDepTransWAbsScale[q-1]->Fill(RPosCorr,FEMCALedep/AbsWidth);
			EDepTransWAbsScaleNorm[q-1]->Fill(RPosCorr,FEMCALedep);

			EDepTransDiff[q-1]->Fill(RPosCorr/MoliereR,FEMCALedep);
			EDepTransDiffCM[q-1]->Fill(RPosCorr,FEMCALedep);

			
			EDepXY[q-1]->Fill(XPosCorr,YPosCorr,FEMCALedep);

			ZPos = (Z1 + Z2)/2;
			
			if(FEMCALLYDep == 4 && MinZ > ZPos) MinZ = ZPos; 
			
			ZPosCorr = ZPos - MinDepthZ;
			

			EDepLongWAbsScale[q-1]->Fill(ZPosCorr,FEMCALedep);
				
			
			if(evtid != evtidPre){
				//cout << "TotalEScindep = " << TotalEScindep << endl;
				ScinEDepHis[q-1]->Fill(TotalEScindep);
				TotalEAllHis[q-1]->Fill(TotalEAll);
				

				TotalEScindep = 0;

				TotalEAll = 0;

				TotalLeakEnergy = 0;

				TotalEvents = TotalEvents + 1.0;

			}

			if(FEMCALLYDep == 4) TotalEScindep = TotalEScindep + FEMCALedep; 
			TotalLeakEnergy =  BHBLedep + BHFLedep + TotalLeakEnergy;
			
			TotalEAll = FEMCALedep + BHBLedep + BHFLedep + TotalEAll;
			TotalE = FEMCALedep + TotalE;

			evtidPre = evtid;

		}

	
	//	LeakPercent = TotalLeakEnergy/(TotalEvents * Energy) * 100;

		if(RadLengthOpt == 0 && AbsOpt == 1) LeakPercent = (1 - TotalE/(TotalEvents * Energy)) * 100;

		LeakEDepHis->Fill(Energy,LeakPercent);



		TF1 * FitFunc = new TF1("FitFunc","gaus",EFitMin,EFitMax);
		ScinEDepHis[q-1]->Fit(FitFunc,"R");

		GausMean = FitFunc->GetParameter(1);
		GausRMS = FitFunc->GetParameter(2);
		GausMeanErr = FitFunc->GetParError(1);
		GausRMSErr = FitFunc->GetParError(2);

		EReso = GausRMS/GausMean;
		EResoErr = EReso * sqrt((GausMeanErr/GausMean)*(GausMeanErr/GausMean) + (GausRMSErr/GausRMS)*(GausRMSErr/GausRMS));
	
		XResoBin = EnergyReso->GetXaxis()->FindBin(Energy);

		EnergyReso->SetBinContent(XResoBin,EReso);
		EnergyReso->SetBinError(XResoBin,EResoErr);



		c->cd();



		TotalEAllHis[q-1]->Draw("hist");
		lat->DrawLatex(0.30,0.80,Form("Incident %s Eneregy = %d GeV",BeamName.Data(),q));	
		lat->DrawLatex(0.30,0.75,Form("Absorber Material: %s",AbsMaterial.Data()));			
		lat->DrawLatex(0.30,0.70,"Sampling Thickness: 1/4 X0");	
		lat->DrawLatex(0.30,0.65,Form("Total RadLength: %d X0",RadLength));		
		lat->DrawLatex(0.30,0.60,Form("Scintillator: %.1f mm Polystyrene",ScinThick*10));
		c->SaveAs(Form("EMCALComp/%s/%d/TotalEnergyHis/TotalEAllHis_%d.png",AbsMaterial.Data(),RadLength,q));

		EDepXY[q-1]->Scale(1.0/EDepXY[q-1]->Integral());
		EDepXY[q-1]->SetMinimum(0);
		EDepXY[q-1]->Draw("COLZ");
		lat->DrawLatex(0.30,0.80,Form("Incident %s Eneregy = %d GeV",BeamName.Data(),q));	
		lat->DrawLatex(0.30,0.75,Form("Absorber Material: %s",AbsMaterial.Data()));			
		lat->DrawLatex(0.30,0.70,"Sampling Thickness: 1/4 X0");	
		lat->DrawLatex(0.30,0.65,Form("Total RadLength: %d X0",RadLength));		
		lat->DrawLatex(0.30,0.60,Form("Scintillator: %.1f mm Polystyrene",ScinThick*10));

		c->SaveAs(Form("EMCALComp/%s/%d/XYShower/EDepXY_%d.png",AbsMaterial.Data(),RadLength,q));
	

		
		EDepTransDiff[q-1]->Scale(1.0/TotalEvents);
		EDepTransDiff[q-1]->Draw("hist");
		lat->DrawLatex(0.30,0.80,Form("Incident %s Eneregy = %d GeV",BeamName.Data(),q));	
		lat->DrawLatex(0.30,0.75,Form("Absorber Material: %s",AbsMaterial.Data()));			
		lat->DrawLatex(0.30,0.70,"Sampling Thickness: 1/4 X0");	
		lat->DrawLatex(0.30,0.65,Form("Total RadLength: %d X0",RadLength));		
		lat->DrawLatex(0.30,0.60,Form("Scintillator: %.1f mm Polystyrene",ScinThick*10));

		c->SaveAs(Form("EMCALComp/%s/%d/ShowWidthNew/EDepTransDiff_%d.png",AbsMaterial.Data(),RadLength,q));

		float PercentValue;
		float TotalValue = EDepTransDiff[q-1]->Integral();

		for(int i = 1; i < EDepTransInc[q-1]->GetNbinsX() + 1; i++){

			PercentValue = EDepTransDiff[q-1]->Integral(1,i)/TotalValue * 100;
			EDepTransInc[q-1]->SetBinContent(i,PercentValue);
			
		}

		EDepTransInc[q-1]->Scale(1.0);
		EDepTransInc[q-1]->Draw("hist");
		lat->DrawLatex(0.30,0.80,Form("Incident %s Eneregy = %d GeV",BeamName.Data(),q));	
		lat->DrawLatex(0.30,0.75,Form("Absorber Material: %s",AbsMaterial.Data()));			
		lat->DrawLatex(0.30,0.70,"Sampling Thickness: 1/4 X0");	
		lat->DrawLatex(0.30,0.65,Form("Total RadLength: %d X0",RadLength));		
		lat->DrawLatex(0.30,0.60,Form("Scintillator: %.1f mm Polystyrene",ScinThick*10));

		c->SaveAs(Form("EMCALComp/%s/%d/ShowWidthNew/EDepTransInc_%d.png",AbsMaterial.Data(),RadLength,q));





		EDepTransWAbsScale[q-1]->Scale(1.0/TotalEvents);
		EDepTransWAbsScale[q-1]->Draw();
		lat->DrawLatex(0.30,0.80,Form("Incident %s Eneregy = %d GeV",BeamName.Data(),q));	
		lat->DrawLatex(0.30,0.75,Form("Absorber Material: %s",AbsMaterial.Data()));			
		lat->DrawLatex(0.30,0.70,"Sampling Thickness: 1/4 X0");	
		lat->DrawLatex(0.30,0.65,Form("Total RadLength: %d X0",RadLength));		
		lat->DrawLatex(0.30,0.60,Form("Scintillator: %.1f mm Polystyrene",ScinThick*10));

		c->SaveAs(Form("EMCALComp/%s/%d/ShowWidth/EDepTransWAbsScale_%d.png",AbsMaterial.Data(),RadLength,q));

		float PercentValue2;
		float TotalValue2 = EDepTransDiffCM[q-1]->Integral();

		

		for(int i = 1; i < EDepTransCum[q-1]->GetNbinsX() + 1; i++){

			PercentValue2 = EDepTransDiffCM[q-1]->Integral(1,i)/TotalValue2 * 100;
			//PercentValue2 = EDepTransWAbsScaleNorm[q-1]->Integral(1,i);

			EDepTransCum[q-1]->SetBinContent(i,PercentValue2);
			
		}



		EDepTransCum[q-1]->Scale(1.0);
		EDepTransCum[q-1]->Draw("hist");
		lat->DrawLatex(0.30,0.80,Form("Incident %s Eneregy = %d GeV",BeamName.Data(),q));	
		lat->DrawLatex(0.30,0.75,Form("Absorber Material: %s",AbsMaterial.Data()));			
		lat->DrawLatex(0.30,0.70,"Sampling Thickness: 1/4 X0");	
		lat->DrawLatex(0.30,0.65,Form("Total RadLength: %d X0",RadLength));		
		lat->DrawLatex(0.30,0.60,Form("Scintillator: %.1f mm Polystyrene",ScinThick*10));

		c->SaveAs(Form("EMCALComp/%s/%d/ShowWidthNew/EDepTransCum_%d.png",AbsMaterial.Data(),RadLength,q));


		EDepTransWAbsScaleNorm[q-1]->Scale(1.0/EDepTransWAbsScaleNorm[q-1]->Integral());
		EDepTransWAbsScaleNorm[q-1]->Draw();
		lat->DrawLatex(0.30,0.80,Form("Incident %s Eneregy = %d GeV",BeamName.Data(),q));	
		lat->DrawLatex(0.30,0.75,Form("Absorber Material: %s",AbsMaterial.Data()));			
		lat->DrawLatex(0.30,0.70,"Sampling Thickness: 1/4 X0");	
		lat->DrawLatex(0.30,0.65,Form("Total RadLength: %d X0",RadLength));		
		lat->DrawLatex(0.30,0.60,Form("Scintillator: %.1f mm Polystyrene",ScinThick*10));

		c->SaveAs(Form("EMCALComp/%s/%d/ShowWidthNorm/EDepTransWAbsScale_%d.png",AbsMaterial.Data(),RadLength,q));


		EDepLongWAbsScale[q-1]->Scale(1.0/TotalEvents);
		EDepLongWAbsScale[q-1]->Draw("hist");
		lat->DrawLatex(0.30,0.80,Form("Incident %s Eneregy = %d GeV",BeamName.Data(),q));	
		lat->DrawLatex(0.30,0.75,Form("Absorber Material: %s",AbsMaterial.Data()));			
		lat->DrawLatex(0.30,0.70,"Sampling Thickness: 1/4 X0");	
		lat->DrawLatex(0.30,0.65,Form("Total RadLength: %d X0",RadLength));		
		lat->DrawLatex(0.30,0.60,Form("Scintillator: %.1f mm Polystyrene",ScinThick*10));

		c->SaveAs(Form("EMCALComp/%s/%d/Longitudinal/EDepLongWAbsScale_%d.png",AbsMaterial.Data(),RadLength,q));

		ScinEDepHis[q-1]->Draw("ep");
		lat->DrawLatex(0.30,0.80,Form("Incident %s Eneregy = %d GeV",BeamName.Data(),q));	
		lat->DrawLatex(0.30,0.75,Form("Absorber Material: %s",AbsMaterial.Data()));			
		lat->DrawLatex(0.30,0.70,"Sampling Thickness: 1/4 X0");	
		lat->DrawLatex(0.30,0.65,Form("Total RadLength: %d X0",RadLength));				
		lat->DrawLatex(0.30,0.60,Form("Scintillator: %.1f mm Polystyrene",ScinThick*10));
		c->SaveAs(Form("EMCALComp/%s/%d/EDepScin/ScinEDepHis_%d.png",AbsMaterial.Data(),RadLength,q));


		cLog->cd();
		EDepTransWAbsScale[q-1]->Draw();
		lat->DrawLatex(0.30,0.80,Form("Incident %s Eneregy = %d GeV",BeamName.Data(),q));	
		lat->DrawLatex(0.30,0.75,Form("Absorber Material: %s",AbsMaterial.Data()));			
		lat->DrawLatex(0.30,0.70,"Sampling Thickness: 1/4 X0");	
		lat->DrawLatex(0.30,0.45,Form("Total RadLength: %d X0",RadLength));
		lat->DrawLatex(0.30,0.60,Form("Scintillator: %.1f mm Polystyrene",ScinThick*10));	
		cLog->SaveAs(Form("EMCALComp/%s/%d/ShowWidth/EDepTransWAbsLogScale_%d.png",AbsMaterial.Data(),RadLength,q));

		EDepTransWAbsScaleNorm[q-1]->Draw();
		lat->DrawLatex(0.30,0.80,Form("Incident %s Eneregy = %d GeV",BeamName.Data(),q));	
		lat->DrawLatex(0.30,0.75,Form("Absorber Material: %s",AbsMaterial.Data()));			
		lat->DrawLatex(0.30,0.70,"Sampling Thickness: 1/4 X0");	
		lat->DrawLatex(0.30,0.65,Form("Total RadLength: %d X0",RadLength));		
		lat->DrawLatex(0.30,0.60,Form("Scintillator: %.1f mm Polystyrene",ScinThick*10));

		cLog->SaveAs(Form("EMCALComp/%s/%d/ShowWidthNorm/EDepTransWAbsLogScale_%d.png",AbsMaterial.Data(),RadLength,q));

		cBoth->cd(1);
		ScinEDepHis[q-1]->Draw("ep");
		lat->DrawLatex(0.30,0.80,Form("Incident %s Eneregy = %d GeV",BeamName.Data(),q));	
		lat->DrawLatex(0.30,0.75,Form("Absorber Material: %s",AbsMaterial.Data()));			
		lat->DrawLatex(0.30,0.70,"Sampling Thickness: 1/4 X0");	
		lat->DrawLatex(0.30,0.65,Form("Total RadLength: %d X0",RadLength));		
		lat->DrawLatex(0.30,0.60,Form("Scintillator: %.1f mm Polystyrene",ScinThick*10));
		cBoth->cd(2);
		EDepTransWAbsScale[q-1]->Draw();


		//Save Files//
		
	}



	c->cd();

	TF1 * FitReso = new TF1("FitReso","sqrt([0] * [0] + [1] * [1]/x)",ERecoMin,ERecoMax);
	FitReso->SetLineColor(ColorFunc);	

	EnergyReso->Fit(FitReso,"R");

	ConstTerm = FitReso->GetParameter(0) * 100;
	StatTerm = FitReso->GetParameter(1) * 100;
	EnergyReso->Draw("ep");

	TString FuncName = Form("#Delta E/E = %.2f %%/#sqrt{E} #oplus %.2f %%",StatTerm,ConstTerm);


	TLegend* leg = new TLegend(0.30,0.65,0.75,0.85,NULL,"brNDC");
	leg->SetBorderSize(0);
	leg->SetTextSize(0.038);
	leg->SetTextFont(42);
	leg->SetFillStyle(0);

	leg->AddEntry("","#it{#bf{Fun4All-EIC}} Simulation","");
	leg->AddEntry("","Geant4, Energy Deposition in Scintillator",""); 
	leg->AddEntry(EnergyReso,"Electron energy resolution","EP");	
	leg->AddEntry(FitReso,FuncName.Data(),"L");	
	leg->Draw("SAME");

	lat->DrawLatex(0.38,0.55,Form("Absorber Material: %s",AbsMaterial.Data()));			
	lat->DrawLatex(0.38,0.50,"Sampling Thickness: 1/4 X0");	
	lat->DrawLatex(0.38,0.45,Form("Total RadLength: %d X0",RadLength));
	lat->DrawLatex(0.38,0.40,Form("Scintillator: %.1f mm Polystyrene",ScinThick*10));	

	c->SaveAs(Form("EMCALComp/%s/Inclusive/EnergyReso_%d.png",AbsMaterial.Data(),RadLength));

	
	LeakEDepHis->Draw("hist");
	
	lat->DrawLatex(0.38,0.55,Form("Absorber Material: %s",AbsMaterial.Data()));			
	lat->DrawLatex(0.38,0.50,"Sampling Thickness: 1/4 X0");	
	lat->DrawLatex(0.38,0.45,Form("Total RadLength: %d X0",RadLength));
	lat->DrawLatex(0.38,0.40,Form("Scintillator: %.1f mm Polystyrene",ScinThick*10));	

	c->SaveAs(Form("EMCALComp/%s/Inclusive/ELeakFrac_%d.png",AbsMaterial.Data(),RadLength));


	

	c4->cd();



	//Not Norm

	TLegend* legComp = new TLegend(0.50,0.65,0.85,0.90,NULL,"brNDC");
	legComp->SetBorderSize(0);
	legComp->SetTextSize(0.038);
	legComp->SetTextFont(42);
	legComp->SetFillStyle(0);


	for(int i = 0; i < NSpecialE; i++){
	
		EDepTransWAbsScale[SpecialE[i]]->SetLineColor(SpecialColor[i]);
		if(i == 0) EDepTransWAbsScale[SpecialE[i]]->Draw("hist");
		if(i > 0) EDepTransWAbsScale[SpecialE[i]]->Draw("histSAME");

		legComp->AddEntry(EDepTransWAbsScale[SpecialE[NSpecialE-1-i]],Form("E = %d GeV",SpecialE[NSpecialE-1-i]),"L");

	}

	legComp->Draw("SAME");


	lat->DrawLatex(0.44,0.55,Form("Absorber Material: %s",AbsMaterial.Data()));			
	lat->DrawLatex(0.44,0.50,"Sampling Thickness: 1/4 X0");	
	lat->DrawLatex(0.44,0.45,Form("Total RadLength: %d X0",RadLength));
	lat->DrawLatex(0.44,0.40,Form("Scintillator: %.1f mm Polystyrene",ScinThick*10));	


	c4->SaveAs(Form("EMCALComp/%s/Inclusive/dE_dR_Comp_%d.png",AbsMaterial.Data(),RadLength));

	//Norm


	TLegend* legCompNorm = new TLegend(0.50,0.65,0.85,0.90,NULL,"brNDC");
	legCompNorm->SetBorderSize(0);
	legCompNorm->SetTextSize(0.038);
	legCompNorm->SetTextFont(42);
	legCompNorm->SetFillStyle(0);


	for(int i = 0; i < NSpecialE; i++){
	
		EDepTransWAbsScaleNorm[SpecialE[i]]->SetLineColor(SpecialColor[i]);
		if(i == 0) EDepTransWAbsScaleNorm[SpecialE[i]]->Draw("hist");
		if(i > 0) EDepTransWAbsScaleNorm[SpecialE[i]]->Draw("histSAME");

		legCompNorm->AddEntry(EDepTransWAbsScaleNorm[SpecialE[NSpecialE-1-i]],Form("E = %d GeV",SpecialE[NSpecialE-1-i]),"L");

	}

	legCompNorm->Draw("SAME");


	lat->DrawLatex(0.44,0.55,Form("Absorber Material: %s",AbsMaterial.Data()));			
	lat->DrawLatex(0.44,0.50,"Sampling Thickness: 1/4 X0");	
	lat->DrawLatex(0.44,0.45,Form("Total RadLength: %d X0",RadLength));
	lat->DrawLatex(0.44,0.40,Form("Scintillator: %.1f mm Polystyrene",ScinThick*10));	


	c4->SaveAs(Form("EMCALComp/%s/Inclusive/dE_dR_Comp_Norm_%d.png",AbsMaterial.Data(),RadLength));


	c4Long->cd();



	TLegend* legLongComp = new TLegend(0.50,0.65,0.85,0.90,NULL,"brNDC");
	legLongComp->SetBorderSize(0);
	legLongComp->SetTextSize(0.038);
	legLongComp->SetTextFont(42);
	legLongComp->SetFillStyle(0);


	for(int i = 0; i < NSpecialE; i++){
	
		EDepLongWAbsScale[SpecialE[i]]->SetLineColor(SpecialColor[i]);
		if(i == 0) EDepLongWAbsScale[SpecialE[i]]->Draw("hist");
		if(i > 0) EDepLongWAbsScale[SpecialE[i]]->Draw("histSAME");

		legLongComp->AddEntry(EDepLongWAbsScale[SpecialE[NSpecialE-1-i]],Form("E = %d GeV",SpecialE[NSpecialE-1-i]),"L");

	}

	legLongComp->Draw("SAME");


	lat->DrawLatex(0.44,0.55,Form("Absorber Material: %s",AbsMaterial.Data()));			
	lat->DrawLatex(0.44,0.50,"Sampling Thickness: 1/4 X0");	
	lat->DrawLatex(0.44,0.45,Form("Total RadLength: %d X0",RadLength));
	lat->DrawLatex(0.44,0.40,Form("Scintillator: %.1f mm Polystyrene",ScinThick*10));	


	c4Long->SaveAs(Form("EMCALComp/%s/Inclusive/dE_dz_Comp_%d.png",AbsMaterial.Data(),RadLength));




	//Write File

	TString outFileName = Form("outFile/PlotStudy_%s_%d.root",AbsMaterial.Data(),RadLength);
	
	TFile * fout = new TFile(outFileName.Data(),"RECREATE");
	fout->cd();



	for(int i = 0; i < NFiles; i++){

		EDepTransWAbsScale[i]->Write();
		ScinEDepHis[i]->Write();
		EDepTransWAbsScaleNorm[i]->Write();
		EDepLongWAbsScale[i]->Write();
		EDepTransInc[i]->Write();
		EDepTransCum[i]->Write();
	}

	EnergyReso->Write();
	LeakEDepHis->Write();
	fout->Close();
	
	cout << "Final: MinZ = " << MinZ << endl;

}
