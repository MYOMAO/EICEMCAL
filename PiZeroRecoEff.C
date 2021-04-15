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
#include "CalibFunc.h"
using namespace std;

using std::cout;
using std::endl;


void PiZeroRecoEff(int AbsOpt){

	TString AbsName;
	if(AbsOpt==0) AbsName = "PHENIX";
	if(AbsOpt==1) AbsName = "WCu";
	if(AbsOpt==2) AbsName = "PHENIXHalf";
	if(AbsOpt==3) AbsName = "PHENIXThird";
	if(AbsOpt==4) AbsName = "PHENIXSixth";


	TString EMCALName;
	if(AbsOpt==0) EMCALName = "PHENIX Pb Shashlik";
	if(AbsOpt==1) EMCALName = "WCu Shashlik 2.5 cm #times 2.5 cm";
	if(AbsOpt==2) EMCALName = "PHENIX Pb Shashlik 2x Granularity";
	if(AbsOpt==3) EMCALName = "PHENIX Pb Shashlik 3x Granularity";
	if(AbsOpt==4) EMCALName = "PHENIX Pb Shashlik 6x Granularity";


	/*
	   TString ThresName;

	   if(ThresOpt== 0) ThresName = "NoThres"; 
	   if(ThresOpt== 1) ThresName = "SomeThres"; 
	   if(ThresOpt== 2) ThresName = "AllThres"; 
	   */

	int RejNonAcc = 1;

	gStyle->SetOptStat(0);


	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();


	TLatex *lat = new TLatex();
	lat->SetNDC();
	lat->SetTextSize(0.035);


	int NBins= 100;
	float EnterZ ;
	if(AbsOpt == 0)	EnterZ = 290.275;
	if(AbsOpt == 1)	EnterZ = 299.93130;
	if(AbsOpt == 2)	EnterZ = 290.275;
	if(AbsOpt == 3)	EnterZ = 290.275;
	if(AbsOpt == 4)	EnterZ = 290.275;

	float OuterRadius = 182.655;


	const int NFiles = 78;
	float EStep = 1;
	float EMin = 0;
	float EMax = EStep * (NFiles + 1) + EMin;
	float EMinHis = EMin - EStep/2;
	float EMaxHis = EMax + EStep/2;

	const int NDecayMax = 100;


	const int NBinsSasha = 18 + 2;
	float PiESasha[NBinsSasha] = {0,10.12142,12.42188,15.87359,17.49123,19.58949,22.65711,24.10811,27.03246,28.97830,31.74008,35.18264,39.28364,41.91129,45.36097,47.16760,50.29416,55.06782,60.00609,70};
	float PiEffSasha[NBinsSasha] = {0,0.0092592592592593,0.021604938271604812,0.03703703703703698,0.12345679012345667,0.25,0.4320987654320988,0.5246913580246914,0.6419753086419754,0.7314814814814815,0.8425925925925927,0.8858024691358025,0.9290123456790123,0.9475308641975309,0.9691358024691359,0.9814814814814815,0.9845679012345678,0.9845679012345678,0.9845679012345678,1};

	/*
	   const int NBinsSasha = 17;
	   float PiESasha[NBinsSasha] = {10.12142,15.87359,17.49123,19.58949,22.65711,24.10811,27.03246,28.97830,31.74008,35.18264,39.28364,41.91129,45.36097,47.16760,50.29416,55.06782,60.00609};
	   float PiEffSasha[NBinsSasha] = {0.0092592592592593,0.03703703703703698,0.12345679012345667,0.25,0.4320987654320988,0.5246913580246914,0.6419753086419754,0.7314814814814815,0.8425925925925927,0.8858024691358025,0.9290123456790123,0.9475308641975309,0.9691358024691359,0.9814814814814815,0.9845679012345678,0.9845679012345678,0.9845679012345678};
	   */




	TGraph * SashaEff = new TGraph(NBinsSasha,PiESasha,PiEffSasha);
	SashaEff->SetMarkerColor(kBlue);
	SashaEff->SetMarkerStyle(20);
	SashaEff->SetMarkerSize(1);
	//SashaEff->SetName("#pi^{0} Merging Probability vs #pi^{0} Energy");
	SashaEff->SetName("SashaEff");
	SashaEff->SetLineColor(kBlue);
	SashaEff->SetLineWidth(2);



	SashaEff->GetXaxis()->SetTitle("#pi^{0} Energy (GeV)");
	SashaEff->GetYaxis()->SetTitle("#pi^{0} Merging Probability");
	SashaEff->GetXaxis()->CenterTitle();
	SashaEff->GetYaxis()->CenterTitle();
	SashaEff->GetYaxis()->SetTitleOffset(1.4);
	SashaEff->SetTitle("#pi^{0} Merging Probability vs #pi^{0} Energy");

	TH1D * PiZeroMassPeak = new TH1D("PiZeroEnergyPeak","",NFiles + 1, EMinHis, EMaxHis);
	PiZeroMassPeak->GetXaxis()->SetTitle("Incident #pi^{0} Energy (GeV)");
	PiZeroMassPeak->GetYaxis()->SetTitle("#pi^{0} Mass Peak (GeV/c^{2})");
	PiZeroMassPeak->GetXaxis()->CenterTitle();
	PiZeroMassPeak->GetYaxis()->CenterTitle();
	PiZeroMassPeak->GetYaxis()->SetTitleOffset(1.5);
	PiZeroMassPeak->SetTitle("#pi^{0} Mass Peak vs Incident #pi^{0} Energy");

	PiZeroMassPeak->SetMarkerStyle(20);
	PiZeroMassPeak->SetMarkerSize(1);
	PiZeroMassPeak->SetMarkerColor(1);
	PiZeroMassPeak->SetLineColor(1);



	TH1D * EventStat = new TH1D("EventStat","",NFiles + 1, EMinHis, EMaxHis);
	EventStat->GetXaxis()->SetTitle("Incident #pi^{0} Energy (GeV)");
	EventStat->GetYaxis()->SetTitle("Total Number of Events");
	EventStat->GetXaxis()->CenterTitle();
	EventStat->GetYaxis()->CenterTitle();
	EventStat->GetYaxis()->SetTitleOffset(1.4);
	EventStat->SetTitle("Event Statistics Summary");




	TH1D * PiZeroMassWidth = new TH1D("PiZeroMassWidth","",NFiles + 1, EMinHis, EMaxHis);
	PiZeroMassWidth->GetXaxis()->SetTitle("Incident #pi^{0} Energy (GeV)");
	PiZeroMassWidth->GetYaxis()->SetTitle("#pi^{0} Mass Width");
	PiZeroMassWidth->GetXaxis()->CenterTitle();
	PiZeroMassWidth->GetYaxis()->CenterTitle();
	PiZeroMassWidth->GetYaxis()->SetTitleOffset(1.4);
	PiZeroMassWidth->SetTitle("#pi^{0} Mass Width vs Incident #pi^{0} Energy");

	PiZeroMassWidth->SetMarkerStyle(20);
	PiZeroMassWidth->SetMarkerSize(1);
	PiZeroMassWidth->SetMarkerColor(1);
	PiZeroMassWidth->SetLineColor(1);


	TH1D * PiZeroMassReso = new TH1D("PiZeroMassReso","",NFiles + 1, EMinHis, EMaxHis);
	PiZeroMassReso->GetXaxis()->SetTitle("Incident #pi^{0} Energy (GeV)");
	PiZeroMassReso->GetYaxis()->SetTitle("#pi^{0} Mass Resolution");
	PiZeroMassReso->GetXaxis()->CenterTitle();
	PiZeroMassReso->GetYaxis()->CenterTitle();
	PiZeroMassReso->GetYaxis()->SetTitleOffset(1.4);
	PiZeroMassReso->SetTitle("#pi^{0} Mass Resolution vs Incident #pi^{0} Energy");

	PiZeroMassReso->SetMarkerStyle(20);
	PiZeroMassReso->SetMarkerSize(1);
	PiZeroMassReso->SetMarkerColor(1);
	PiZeroMassReso->SetLineColor(1);



	TH1D * PiZeroEnergyPeak = new TH1D("PiZeroEnergyPeak","",NFiles + 1, EMinHis, EMaxHis);
	PiZeroEnergyPeak->GetXaxis()->SetTitle("Incident #pi^{0} Energy (GeV)");
	PiZeroEnergyPeak->GetYaxis()->SetTitle("Clustered + Calibrated #pi^{0} Energy");
	PiZeroEnergyPeak->GetXaxis()->CenterTitle();
	PiZeroEnergyPeak->GetYaxis()->CenterTitle();
	PiZeroEnergyPeak->GetYaxis()->SetTitleOffset(1.4);
	PiZeroEnergyPeak->SetTitle("#pi^{0} Peak Energy vs Incident #pi^{0} Energy");

	PiZeroEnergyPeak->SetMarkerStyle(20);
	PiZeroEnergyPeak->SetMarkerSize(1);
	PiZeroEnergyPeak->SetMarkerColor(1);
	PiZeroEnergyPeak->SetLineColor(1);

	int NBinsE = 200;
	float EffMinHis = 0;
	float EffMaxHis = 70;

	int XEBin;



	float Acc;
	float AccErr;


	TH1D * PiZeroRECOEffHis = new TH1D("PiZeroRECOEffHis","", NFiles + 1, EMinHis, EMaxHis);
	PiZeroRECOEffHis->GetXaxis()->SetTitle("#pi^{0} Energy (GeV)");
	PiZeroRECOEffHis->GetYaxis()->SetTitle("#pi^{0} Merging Probability");
	PiZeroRECOEffHis->GetXaxis()->CenterTitle();
	PiZeroRECOEffHis->GetYaxis()->CenterTitle();
	PiZeroRECOEffHis->GetYaxis()->SetTitleOffset(1.4);
	PiZeroRECOEffHis->SetTitle("#pi^{0} Merging Probability vs #pi^{0} Energy");

	PiZeroRECOEffHis->SetMarkerStyle(20);
	PiZeroRECOEffHis->SetMarkerSize(1);
	PiZeroRECOEffHis->SetMarkerColor(1);
	PiZeroRECOEffHis->SetLineColor(1);





	TH1D * PiZeroAcc = new TH1D("PiZeroAcc","", NFiles + 1, EMinHis, EMaxHis);
	PiZeroAcc->GetXaxis()->SetTitle("#pi^{0} Energy (GeV)");
	PiZeroAcc->GetYaxis()->SetTitle("Fraction of #pi^{0} Acceptance by the EMCAL");
	PiZeroAcc->GetXaxis()->CenterTitle();
	PiZeroAcc->GetYaxis()->CenterTitle();
	PiZeroAcc->GetYaxis()->SetTitleOffset(1.4);
	PiZeroAcc->SetTitle("#pi^{0} Acceptance vs #pi^{0} Energy");

	PiZeroAcc->SetMarkerStyle(20);
	PiZeroAcc->SetMarkerSize(1);
	PiZeroAcc->SetMarkerColor(1);
	PiZeroAcc->SetLineColor(1);



	TH2D * PiEffEmpty = new TH2D("PiEffEmpty","", NBinsE, EffMinHis, EffMaxHis,100,0,1.5);
	PiEffEmpty->GetXaxis()->SetTitle("#pi^{0} Energy (GeV)");
	PiEffEmpty->GetYaxis()->SetTitle("#pi^{0} Merging Probability");
	PiEffEmpty->GetXaxis()->CenterTitle();
	PiEffEmpty->GetYaxis()->CenterTitle();
	PiEffEmpty->GetYaxis()->SetTitleOffset(1.4);
	PiEffEmpty->SetTitle("#pi^{0} Merging Probability vs #pi^{0} Energy");

	PiEffEmpty->SetMarkerStyle(20);
	PiEffEmpty->SetMarkerSize(1);
	PiEffEmpty->SetMarkerColor(1);
	PiEffEmpty->SetLineColor(1);






	float ClusDis;




	float E[NDecayMax];
	float x[NDecayMax];
	float y[NDecayMax];
	float z[NDecayMax];
	float r[NDecayMax];
	float RecoR[NDecayMax];
	float Px[NDecayMax];
	float Py[NDecayMax];
	float Pz[NDecayMax];




	TF2 * ThetaCalib; 
	TF2 * TruthECalib; 
	TF2 * TruthELowCalib  = new TF2("TruthELowCalib",TruthEFuncFormLow.Data());
 

	if(AbsOpt == 0 || AbsOpt == 2 || AbsOpt == 3 || AbsOpt == 4){

		ThetaCalib = new TF2("ThetaCalib",ThetaFuncForm.Data());
		TruthECalib = new TF2("TruthECalib",TruthEFuncForm.Data());
		
	}

	if(AbsOpt == 1){
		ThetaCalib = new TF2("WCuThetaCalib",WCuThetaFuncForm.Data());
		TruthECalib = new TF2("WCuTruthECalib",WCuTruthEFuncForm.Data());

	}


	TLorentzVector * Gamma1 = new TLorentzVector;
	TLorentzVector * Gamma2 = new TLorentzVector;
	TLorentzVector * PiZero = new TLorentzVector;


	float PiZeroMass;
	float PiZeroEnergy;




	float PiZeroMassPDG = 0.1349766;
	float PiMassFraction = 0.1;	
	float PiMassFitLow = PiZeroMassPDG * (1 - PiMassFraction);
	float PiMassFitHigh = PiZeroMassPDG * (1 + PiMassFraction);


	float PiEnergyHigh;
	float PiEnergyLow;

	float PiEnergyReso = 0.2;



	float MassUpRatio = 1.3;
	float MassLowRatio = 0.7;

	float PiMassUp = PiZeroMassPDG * MassUpRatio;
	float PiMassLow = PiZeroMassPDG * MassLowRatio;
	//	NFiles = 2;

	//PiMassFitLow = PiMassLow;
	//PiMassFitHigh = PiMassUp;


	const int NClusBin = 21;
	float NClusMin = -0.5;
	float NClusMax = 20.5;

	float Theta;
	float PVx = 0;
	float PVy = 0;
	float PVz = 0;


	float NPiEvents; 

	float ClusDisMin;

	int TotalEvents = 0;


	float Theta1Truth;
	float Theta2Truth;

	float Theta1Reco;
	float Theta2Reco;

	float Theta1Diff;
	float Theta2Diff;

	int NThetaDiff = 100;
	float ThetaDiffMin = -0.5;
	float ThetaDiffMax = 0.5;


	float Phi1Truth;
	float Phi2Truth;

	float Phi1Reco;
	float Phi2Reco;

	float Phi1Diff;
	float Phi2Diff;

	int NPhiDiff = 100;
	float PhiDiffMin = -3;
	float PhiDiffMax = 3;

	float AngleTruthMax = 0.8;

	int NClusE = 100;

	int FileName[21] = {0,3075652,3076307,3076330,3351899,3076363,3351938,3351958,3351984,3076450,3076472,2579748,3352071,3352079,3352107,3352126,3352142,3352159,3352181,3352208,3352228};

	float ClusEMin = 0.10; 
	float NSigma = 3;
	float EResoWidth = 0.06 * sqrt(ClusEMin);
	float ClusEMinTruth = ClusEMin + EResoWidth * NSigma;

	ClusEMinTruth = 0.20;

	for(int q = 1; q < NFiles; q++){ 

		cout << "Now Working at File = " << q << endl;

		float Energy = (q + 0) * EStep + EMin;

		float PiEDown = Energy * 0.7;
		float PiEUp = Energy * 1.3;

		float ClusDisMin = 2 * PiZeroMassPDG/sqrt(Energy * Energy - PiZeroMassPDG * PiZeroMassPDG) * EnterZ * 0.90;

		float ClusEHisMin = 0;
		float ClusEHisMax = Energy * 2;

		float N1Clus = 0;


		float Eff = 0;
		float EffErr = 0;

		float EventFloat;
		float clusterIDFloat;

		int clusterID;
		int NClus = 0;

		int Event = 0;
		int EventPre = 0;

		float XPos;
		float YPos;
		float ZPos;

		float RECOEnergy;


		int NETruth = 100;
		float ETruthMin = 0;
		float ETruthMax = Energy/2.0 + 0.1;


		float E1Truth;
		float E2Truth;
		float AngleTruth; 

		std::vector<float> E1TruthVec;	
		std::vector<float> E2TruthVec;
		std::vector<float> AngleTruthVec;


		TString infile; 
	
		if(AbsOpt == 0){
			infile = Form("/sphenix/user/zshi/SimulationWork/EICShashlik/macros/macros/g4simulations/JobSub4/ERECO2/%d/G4EICDetector.root_g4femc_eval_%d.root",q,q);
			if(q == 1) infile = Form("/sphenix/user/zshi/SimulationWork/EICShashlik/macros/macros/g4simulations/JobSub4/ERECO4/G4EICDetector.root_g4femc_eval_%d.root",q);
		}
		if(AbsOpt == 1){
			infile = Form("/sphenix/user/zshi/SimulationWork/EICShashlik/macros/macros/g4simulations/JobSub4/ERECOWCu/%d/G4EICDetector.root_g4femc_eval_%d.root",q,q);
		//	infile = Form("/sphenix/user/zshi/SimulationWork/EICShashlik/macros/macros/g4simulations/JobSub4/ERECOHalf/%d/G4EICDetector.root_g4femc_eval_%d.root",q,q);

		}
		if(AbsOpt == 2){
		//	infile = Form("/sphenix/user/zshi/SimulationWork/EICShashlik/macros/macros/g4simulations/JobSub4/ERECOWCu/%d/G4EICDetector.root_g4femc_eval_%d.root",q,q);
			infile = Form("/sphenix/user/zshi/SimulationWork/EICShashlik/macros/macros/g4simulations/JobSub4/ERECOHalf/%d/G4EICDetector.root_g4femc_eval_%d.root",q,q);

		}

		if(AbsOpt == 3){
		//	infile = Form("/sphenix/user/zshi/SimulationWork/EICShashlik/macros/macros/g4simulations/JobSub4/ERECOWCu/%d/G4EICDetector.root_g4femc_eval_%d.root",q,q);
			infile = Form("/sphenix/user/zshi/SimulationWork/EICShashlik/macros/macros/g4simulations/JobSub4/ERECOThird/%d/G4EICDetector.root_g4femc_eval_%d.root",q,q);

		}
		if(AbsOpt == 4){
		//	infile = Form("/sphenix/user/zshi/SimulationWork/EICShashlik/macros/macros/g4simulations/JobSub4/ERECOWCu/%d/G4EICDetector.root_g4femc_eval_%d.root",q,q);
			infile = Form("/sphenix/user/zshi/SimulationWork/EICShashlik/macros/macros/g4simulations/JobSub4/ERECOSixth/%d/G4EICDetector.root_g4femc_eval_%d.root",q,q);

		}

		//	TString infile = Form("/sphenix/user/zshi/SimulationWork/EICShashlik/macros/macros/g4simulations/JobSub4/ERECO2/%d/G4EICDetector.root_g4femc_eval_%d.root",q,FileName[q]);

		//	if(q == 1) infile = Form("/sphenix/user/zshi/SimulationWork/EICShashlik/macros/macros/g4simulations/JobSub4/ERECO2/%d/G4EICDetector.root_g4femc_eval_%d.root",q,FileName[q]);


		//	if(q == 1) infile = Form("/sphenix/user/zshi/SimulationWork/EICShashlik/macros/macros/g4simulations/JobSub4/ERECO3/G4EICDetector.root_g4femc_eval_%d.root",q);


		//	TString infile = Form("/sphenix/user/zshi/SimulationWork/EICShashlik/macros/macros/g4simulations/JobSub4/ERECO4/G4EICDetector.root_g4femc_eval_%d.root",q);

		//TString infile = Form("/sphenix/user/zshi/SimulationWork/EICShashlik/macros/macros/g4simulations/JobSub3/ERECOSucked/G4EICDetector.root_g4femc_eval_%d.root",q);



		TFile * fin = new TFile(infile.Data());
		fin->cd();



		TTree * ntp_gshower = (TTree * ) fin->Get("ntp_gshower");
		NPiEvents = ntp_gshower->GetEntries();

		TTree * ntp_cluster = (TTree * ) fin->Get("ntp_cluster");

		ntp_cluster->SetBranchAddress("event",&EventFloat);
		ntp_cluster->SetBranchAddress("clusterID",&clusterIDFloat);
		ntp_cluster->SetBranchAddress("e",&RECOEnergy);
		ntp_cluster->SetBranchAddress("x",&XPos);
		ntp_cluster->SetBranchAddress("y",&YPos);
		ntp_cluster->SetBranchAddress("z",&ZPos);



		TH1D * PiMassHis = new TH1D("PiMassHis","",NBins,PiMassLow,PiMassUp);
		PiMassHis->GetXaxis()->SetTitle(Form("Reconstructed m(#pi^{0} -> #gamma #gamma) (GeV/c^{2})"));
		PiMassHis->GetYaxis()->SetTitle("Counts");
		PiMassHis->GetXaxis()->CenterTitle();
		PiMassHis->GetYaxis()->CenterTitle();
		PiMassHis->GetYaxis()->SetTitleOffset(1.4);
		PiMassHis->SetTitle(Form("Incident #pi Energy = %.1f GeV",Energy));
		PiMassHis->SetTitle("#pi Invariant Mass Distribution");


		PiMassHis->SetMarkerStyle(20);
		PiMassHis->SetMarkerSize(1);
		PiMassHis->SetMarkerColor(1);
		PiMassHis->SetLineColor(1);



		TH1D * PiEnergyHis = new TH1D("PiEnergyHis","",NBins,PiEDown,PiEUp);
		PiEnergyHis->GetXaxis()->SetTitle("Reconstructed #pi^{0} Energy (GeV)");
		PiEnergyHis->GetYaxis()->SetTitle("Counts");
		PiEnergyHis->GetXaxis()->CenterTitle();
		PiEnergyHis->GetYaxis()->CenterTitle();
		PiEnergyHis->GetYaxis()->SetTitleOffset(1.4);
		PiEnergyHis->SetTitle(Form("Incident #pi Energy = %.1f GeV",Energy));



		PiEnergyHis->SetMarkerStyle(20);
		PiEnergyHis->SetMarkerSize(1);
		PiEnergyHis->SetMarkerColor(1);
		PiEnergyHis->SetLineColor(1);



		TH1D *  NClusHis = new TH1D("NClusHis","",NClusBin,NClusMin,NClusMax);
		NClusHis->GetXaxis()->SetTitle("Number of Clusters");
		NClusHis->GetYaxis()->SetTitle("Fraction of Events");
		NClusHis->GetXaxis()->CenterTitle();
		NClusHis->GetYaxis()->CenterTitle();
		NClusHis->GetYaxis()->SetTitleOffset(1.4);
		NClusHis->SetTitle(Form("Incident #pi Energy = %.1f GeV",Energy));





		TH1D *  ClusEHis = new TH1D("ClusEHis","",NClusE,ClusEHisMin,ClusEHisMax);
		ClusEHis->GetXaxis()->SetTitle("Clustered Energy (GeV)");
		ClusEHis->GetYaxis()->SetTitle("Number of Clusters");
		ClusEHis->GetXaxis()->CenterTitle();
		ClusEHis->GetYaxis()->CenterTitle();
		ClusEHis->GetYaxis()->SetTitleOffset(1.4);
		ClusEHis->SetTitle(Form("Incident #pi Energy = %.1f GeV",Energy));





		TH1D *  ThetaDiff = new TH1D("ThetaDiff","",NThetaDiff,ThetaDiffMin,ThetaDiffMax);
		ThetaDiff->GetXaxis()->SetTitle("Truth #theta - Reco #theta ");
		ThetaDiff->GetYaxis()->SetTitle("Fraction of Events");
		ThetaDiff->GetXaxis()->CenterTitle();
		ThetaDiff->GetYaxis()->CenterTitle();
		ThetaDiff->GetYaxis()->SetTitleOffset(1.4);
		ThetaDiff->SetTitle(Form("Incident #pi Energy = %.1f GeV",Energy));




		TH1D *  PhiDiff = new TH1D("PhiDiff","",NPhiDiff,PhiDiffMin,PhiDiffMax);
		PhiDiff->GetXaxis()->SetTitle("Truth #phi - Reco #phi ");
		PhiDiff->GetYaxis()->SetTitle("Fraction of Events");
		PhiDiff->GetXaxis()->CenterTitle();
		PhiDiff->GetYaxis()->CenterTitle();
		PhiDiff->GetYaxis()->SetTitleOffset(1.4);
		PhiDiff->SetTitle(Form("Incident #pi Energy = %.1f GeV",Energy));


		TH1D *  TruthE = new TH1D("TruthE","",NETruth,ETruthMin,ETruthMax);
		TruthE->GetXaxis()->SetTitle("Truth Energy of the Lower Energy Photon (GeV)");
		TruthE->GetYaxis()->SetTitle("Number Events");
		TruthE->GetXaxis()->CenterTitle();
		TruthE->GetYaxis()->CenterTitle();
		TruthE->GetYaxis()->SetTitleOffset(1.4);
		TruthE->SetTitle(Form("Incident #pi Energy = %.1f GeV",Energy));




		int NEvents = ntp_cluster->GetEntries();

		TotalEvents = TotalEvents + NEvents;

		//NPiEvents = NEvents;

		float TotalGoodEvents = 0;
		/*
		   if(ThresOpt == 0){

		   ClusEMin = 0;
		   }


		   if(ThresOpt == 1){
		   if(q < 6)	ClusEMin = 0;
		   if(q > 5)	ClusEMin = 0.1;
		   }



		   if(ThresOpt == 2){

		   ClusEMin = 0.1;
		   }
		   */


		//Prepare to Reject Non Accepted Events//

		float GPX1;
		float GPY1;
		float GPZ1;

		float GPX2;
		float GPY2;
		float GPZ2;

		TString infileGammaInfo;

		if(AbsOpt == 0){
			infileGammaInfo = Form("/sphenix/user/zshi/SimulationWork/EICShashlik/macros/macros/g4simulations/JobSub4/ERECO2/%d/PiZeroGammaInfo_%d.root",q,q);
			if(q == 1) infileGammaInfo = Form("/sphenix/user/zshi/SimulationWork/EICShashlik/macros/macros/g4simulations/JobSub4/ERECO4/PiZeroGammaInfo_%d.root",q);
		}
		if(AbsOpt == 1){
			infileGammaInfo = Form("/sphenix/user/zshi/SimulationWork/EICShashlik/macros/macros/g4simulations/JobSub4/ERECOWCu/%d/PiZeroGammaInfo_%d.root",q,q);
	//		infileGammaInfo = Form("/sphenix/user/zshi/SimulationWork/EICShashlik/macros/macros/g4simulations/JobSub4/ERECOHalf/%d/PiZeroGammaInfo_%d.root",q,q);

		}
		if(AbsOpt == 2){
	//		infileGammaInfo = Form("/sphenix/user/zshi/SimulationWork/EICShashlik/macros/macros/g4simulations/JobSub4/ERECOWCu/%d/PiZeroGammaInfo_%d.root",q,q);
			infileGammaInfo = Form("/sphenix/user/zshi/SimulationWork/EICShashlik/macros/macros/g4simulations/JobSub4/ERECOHalf/%d/PiZeroGammaInfo_%d.root",q,q);
		}
		if(AbsOpt == 3){
	//		infileGammaInfo = Form("/sphenix/user/zshi/SimulationWork/EICShashlik/macros/macros/g4simulations/JobSub4/ERECOWCu/%d/PiZeroGammaInfo_%d.root",q,q);
			infileGammaInfo = Form("/sphenix/user/zshi/SimulationWork/EICShashlik/macros/macros/g4simulations/JobSub4/ERECOThird/%d/PiZeroGammaInfo_%d.root",q,q);
		}
		if(AbsOpt == 4){
	//		infileGammaInfo = Form("/sphenix/user/zshi/SimulationWork/EICShashlik/macros/macros/g4simulations/JobSub4/ERECOWCu/%d/PiZeroGammaInfo_%d.root",q,q);
			infileGammaInfo = Form("/sphenix/user/zshi/SimulationWork/EICShashlik/macros/macros/g4simulations/JobSub4/ERECOSixth/%d/PiZeroGammaInfo_%d.root",q,q);
		}


		//	TString infileGammaInfo = Form("/sphenix/user/zshi/SimulationWork/EICShashlik/macros/macros/g4simulations/JobSub4/ERECO2/%d/PiZeroGammaInfo_%d.root",q,FileName[q]);
		//	TString infileGammaInfo = Form("/sphenix/user/zshi/SimulationWork/EICShashlik/macros/macros/g4simulations/JobSub4/ERECO4/PiZeroGammaInfo_%d.root",q);

		//	if(q == 1)  infileGammaInfo = Form("/sphenix/user/zshi/SimulationWork/EICShashlik/macros/macros/g4simulations/JobSub4/ERECO2/%d/PiZeroGammaInfo_%d.root",q,FileName[q]);
		//	if(q == 1) infileGammaInfo = Form("/sphenix/user/zshi/SimulationWork/EICShashlik/macros/macros/g4simulations/JobSub4/ERECO3/PiZeroGammaInfo_%d.root",q);




		TFile * finGammaInfo = new TFile(infileGammaInfo.Data());
		finGammaInfo->cd();

		TTree * EvtInfo = (TTree *) finGammaInfo->Get("EvtInfo");

		EvtInfo->SetBranchAddress("GPX1",&GPX1);
		EvtInfo->SetBranchAddress("GPY1",&GPY1);
		EvtInfo->SetBranchAddress("GPZ1",&GPZ1);
		EvtInfo->SetBranchAddress("GPX2",&GPX2);
		EvtInfo->SetBranchAddress("GPY2",&GPY2);
		EvtInfo->SetBranchAddress("GPZ2",&GPZ2);

		int NEventGamma = EvtInfo->GetEntries();

		std::vector<int> EvtAccVec;
		float CalR1;
		float CalR2;

		int IsAcc;

		int NonAccEvt = 0;
		std::vector<float> TruthVec1;
		std::vector<float> TruthVec2;


		if(RejNonAcc == 1){

			for(int i = 0; i < NEventGamma; i ++){

				EvtInfo->GetEntry(i);

				CalR1 = sqrt(GPX1 * GPX1 + GPY1 * GPY1)/GPZ1 * EnterZ;
				CalR2 = sqrt(GPX2 * GPX2 + GPY2 * GPY2)/GPZ2 * EnterZ;
				IsAcc = 0;

				Theta1Truth = TMath::ATan2(sqrt(GPX1 * GPX1 + GPY1 * GPY1),GPZ1);
				Theta2Truth = TMath::ATan2(sqrt(GPX2 * GPX2 + GPY2 * GPY2),GPZ2);

				E1Truth = 0.5* Energy + 0.5 * sqrt(Energy * Energy - 2 * PiZeroMassPDG * PiZeroMassPDG/(1 - GPX1 * GPX2 - GPY1 * GPY2 - GPZ1 * GPZ2));
				E2Truth = 0.5* Energy - 0.5 * sqrt(Energy * Energy - 2 * PiZeroMassPDG * PiZeroMassPDG/(1 - GPX1 * GPX2 - GPY1 * GPY2 - GPZ1 * GPZ2));
				AngleTruth = abs(E1Truth - E2Truth)/(E1Truth + E2Truth);

				TruthVec1.push_back(Theta1Truth);
				TruthVec2.push_back(Theta2Truth);

				E1TruthVec.push_back(E1Truth);				
				E2TruthVec.push_back(E2Truth);
				AngleTruthVec.push_back(AngleTruth);

				TruthE->Fill(E2Truth);




				if(GPZ1 > 0 && GPZ2 > 0 && CalR1 < OuterRadius  && CalR2 < OuterRadius ) IsAcc = 1;

				EvtAccVec.push_back(IsAcc);

				if(IsAcc == 0) NonAccEvt  = NonAccEvt + 1;
			}

			cout << "Check: " << EvtAccVec.size() <<  "   E2TruthVec.size() = " << E2TruthVec.size() << "   NEventGamma =   " << NEventGamma <<  "   NPiEvents =  " << NPiEvents << endl;

		}
		//Done Preparation Bro//

		cout << "NonAccEvt = " << NonAccEvt << endl;



		int Index = 0;

		for(int i = 0; i < NEvents; i++){

			ntp_cluster->GetEntry(i);

			Event = int(EventFloat);
			clusterID = int(clusterIDFloat);



			if(EventPre != Event && EvtAccVec[Index] == 1 && E2TruthVec[Index] > ClusEMinTruth){
		//	if(EventPre != Event && EvtAccVec[Index] == 1 && AngleTruthVec[Index] < AngleTruthMax){

				TotalGoodEvents = TotalGoodEvents + 1;

				ClusDis = sqrt((x[0] - x[1]) * (x[0] - x[1]) + (y[0] - y[1]) * (y[0] - y[1]));


				if(NClus == 1) N1Clus = N1Clus + 1;

				if(NClus == 2 && ClusDis > ClusDisMin){
					Gamma1->SetPxPyPzE(Px[0],Py[0],Pz[0],E[0]);
					Gamma2->SetPxPyPzE(Px[1],Py[1],Pz[1],E[1]);


					Theta1Reco = TMath::ATan2(sqrt(Px[0] * Px[0] + Py[0] * Py[0]),Pz[0]);
					Theta2Reco = TMath::ATan2(sqrt(Px[1] * Px[1] + Py[1] * Py[1]),Pz[1]);

					Theta1Diff = TruthVec1[Index] - Theta1Reco;
					Theta2Diff = TruthVec2[Index] - Theta2Reco;

					ThetaDiff->Fill(Theta1Diff + Theta2Diff);

					* PiZero = * Gamma1 + * Gamma2;

					PiZeroMass = PiZero->M();
					PiZeroEnergy = PiZero->E();

					PiMassHis->Fill(PiZeroMass);
					PiEnergyHis->Fill(PiZeroEnergy);




				}

				NClusHis->Fill(NClus);
				NClus = 0;



			}
			if(EventPre != Event) Index = Index + 1;


			E[clusterID] = RECOEnergy;
			x[clusterID] = XPos - PVx;
			y[clusterID] = YPos - PVy;
			z[clusterID] = ZPos - PVz;		
			//z[clusterID] = EnterZ - 5;  //PHENIX Pb EMCAL Front Face z coordinate

			r[clusterID] = sqrt(x[clusterID]*x[clusterID] + y[clusterID]*y[clusterID] + z[clusterID] * z[clusterID]);
			RecoR[clusterID] = sqrt(x[clusterID]*x[clusterID] + y[clusterID]*y[clusterID]);

			Theta = ThetaCalib->Eval(RecoR[clusterID],E[clusterID]);
			E[clusterID] = TruthECalib->Eval(RecoR[clusterID],E[clusterID]);

			//		if(E[clusterID] > 5)E[clusterID] = TruthECalib->Eval(RecoR[clusterID],E[clusterID]);
			//		if(E[clusterID] < 5) E[clusterID] = TruthELowCalib->Eval(RecoR[clusterID],E[clusterID]);



			Px[clusterID] = E[clusterID] * sin(Theta) * x[clusterID]/RecoR[clusterID];
			Py[clusterID] = E[clusterID] * sin(Theta) * y[clusterID]/RecoR[clusterID];
			Pz[clusterID] = E[clusterID] * cos(Theta);



			if(RECOEnergy > ClusEMin && EvtAccVec[Index] == 1 && E2TruthVec[Index] > ClusEMinTruth) NClus = NClus + 1;
			//if(RECOEnergy > ClusEMin && EvtAccVec[Index] == 1 && AngleTruthVec[Index] < AngleTruthMax) NClus = NClus + 1;	
			EventPre = Event;

			ClusEHis->Fill(E[clusterID]);






		}






		ClusEHis->Draw();
		c->SaveAs(Form("PiZeroRECOEffPlots/%s/ClusEHis/ClusEHis_%d.png",AbsName.Data(),q));

		cout << "TotalGoodEvents = " << TotalGoodEvents << "   NonAccEvt =   " << NonAccEvt <<   " TotalGoodEvents + NonAccEvt " << TotalGoodEvents + NonAccEvt << "   NPiEvents = " << NPiEvents << endl  << "   Index = " << Index << endl;

		NClusHis->Scale(1.0/NPiEvents);
		NClusHis->Draw("hist");

		cout << "NClusHis->Integral() = " << NClusHis->Integral() << endl;

		c->SaveAs(Form("PiZeroRECOEffPlots/%s/NClus/NClus_%d.png",AbsName.Data(),q));



		Eff = N1Clus/TotalGoodEvents;

		EffErr = 1/sqrt(N1Clus) * Eff;
		
		EffErr = 0.000001;

		cout << "Eff = " << Eff << "   EffErr = " << EffErr << endl;

		XEBin = PiZeroRECOEffHis->GetXaxis()->FindBin(Energy);

		PiZeroRECOEffHis->SetBinContent(XEBin,Eff);
		PiZeroRECOEffHis->SetBinError(XEBin,EffErr);

		EventStat->SetBinContent(XEBin,NEvents);


		Acc = TotalGoodEvents/NPiEvents;
		AccErr = 1/sqrt(TotalGoodEvents) * Acc;

		PiZeroAcc->SetBinContent(XEBin,Acc);
		PiZeroAcc->SetBinError(XEBin,AccErr);


		//Plotting 

		PiEnergyHigh = Energy * (1 + PiEnergyReso + 0.12);
		PiEnergyLow = Energy * (1 - PiEnergyReso);




		TF1 * PiMassFit = new TF1("PiMassFit","gaus",PiMassFitLow,PiMassFitHigh);
		PiMassFit->SetLineColor(kBlue);


		PiMassHis->Draw("ep");
		PiMassHis->SetMaximum(PiMassHis->GetMaximum() * 1.5);
		PiMassHis->Fit(PiMassFit,"R");
		PiMassFit->Draw("SAME");

		float PiMassPeak = PiMassFit->GetParameter(1);
		float PiMassPeakError = PiMassFit->GetParError(1);
		float PiMassWidth = PiMassFit->GetParameter(2);
		float PiMassWidthError = PiMassFit->GetParError(2);
		float PiMassReso = PiMassWidth/PiMassPeak;
		float PiMassResoError = PiMassReso * sqrt((PiMassPeakError/PiMassPeak) * (PiMassPeakError/PiMassPeak) + (PiMassWidthError/PiMassWidth) * (PiMassWidthError/PiMassWidth));

		lat->DrawLatex(0.15,0.85,Form("Incident #pi^{0} Eneregy = %.1f GeV",Energy));	
		lat->DrawLatex(0.15,0.80,"Standard PHENIX Pb Shashlik EMCAL");	
		lat->DrawLatex(0.15,0.75,Form("#pi^{0} Mass Peak: m = %.1f MeV/c^{2}",PiMassPeak * 1000));	
		lat->DrawLatex(0.15,0.70,Form("#pi^{0} Mass Width: #Delta m = %.1f MeV/c^{2}",PiMassWidth * 1000));	
		lat->DrawLatex(0.15,0.65,Form("#pi^{0} Mass Resolution: #Delta m/m = %.1f%%",PiMassReso * 100));	

		cout << "Pass 3" << endl;


		c->SaveAs(Form("PiZeroRECOEffPlots/%s/PiMass/PiMass_%d.png",AbsName.Data(),q));






		TF1 * PiEFit = new TF1("PiEFit","gaus",PiEnergyLow,PiEnergyHigh);
		PiEFit->SetLineColor(kBlue);

		PiEnergyHis->Draw("ep");
		PiEnergyHis->Fit(PiEFit,"R");
		PiEnergyHis->Draw("SAME");




		float PiEnergyPeak = PiEFit->GetParameter(1);
		float PiEnergyPeakError = PiEFit->GetParError(1);
		float PiEnergyWidth = PiEFit->GetParameter(2);
		float PiEnergyWidthError = PiEFit->GetParError(2);
		float PiEnergyReso = PiEnergyWidth/PiEnergyPeak;
		float PiEnergyResoError = PiEnergyReso * sqrt((PiEnergyPeakError/PiEnergyPeak) * (PiEnergyPeakError/PiEnergyPeak) + (PiEnergyWidthError/PiEnergyWidth) * (PiEnergyWidthError/PiEnergyWidth));


		c->SaveAs(Form("PiZeroRECOEffPlots/%s/PiEnergy/PiEnergy_%d.png",AbsName.Data(),q));


		PiZeroMassPeak->SetBinContent(q+1,PiMassPeak);
		PiZeroMassPeak->SetBinError(q+1,PiMassPeakError);

		PiZeroMassWidth->SetBinContent(q+1,PiMassWidth);
		PiZeroMassWidth->SetBinError(q+1,PiMassWidthError);

		PiZeroMassReso->SetBinContent(q+1,PiMassReso);
		PiZeroMassReso->SetBinError(q+1,PiMassResoError);


		PiZeroEnergyPeak->SetBinContent(q+1,PiEnergyPeak);
		PiZeroEnergyPeak->SetBinError(q+1,PiEnergyPeakError);


		cout << "NCluster = " << NEvents << "   NPiEvents = " << NPiEvents << endl; 

		EvtAccVec.clear();
		TruthVec1.clear();
		TruthVec2.clear();
		E1TruthVec.clear();
		E2TruthVec.clear();


		ThetaDiff->Draw("hist");
		c->SaveAs(Form("PiZeroRECOEffPlots/%s/ThetaCheck/ThetaDiff_%d.png",AbsName.Data(),q));


		TruthE->Draw("hist");
		c->SaveAs(Form("PiZeroRECOEffPlots/%s/TruthE/TruthE_%d.png",AbsName.Data(),q));


		}




		EventStat->Draw("hist");
		c->SaveAs(Form("PiZeroRECOEffPlots/%s/General/EventStat.png",AbsName.Data()));


		//Sanity Check//

		PiZeroMassPeak->SetMinimum(0);
		PiZeroMassPeak->SetMaximum(0.20);

		PiZeroMassPeak->Draw("ep");
		TLine *line = new TLine(EMin,PiZeroMassPDG,EMax,PiZeroMassPDG);
		PiZeroMassPeak->GetXaxis()->SetRangeUser(0,20);
		line = new TLine(0,PiZeroMassPDG,20,PiZeroMassPDG);

		line->SetLineStyle(2);
		line->SetLineWidth(2);
		line->SetLineColor(2);
		line->Draw("SAME");
		c->SaveAs(Form("PiZeroRECOEffPlots/%s/SanityCheck/PiZeroMassPeak.png",AbsName.Data()));

		PiZeroMassWidth->Draw("ep");
		c->SaveAs(Form("PiZeroRECOEffPlots/%s/SanityCheck/PiZeroMassWidth.png",AbsName.Data()));

		PiZeroMassReso->Draw("ep");
		c->SaveAs(Form("PiZeroRECOEffPlots/%s/SanityCheck/PiZeroMassReso.png",AbsName.Data()));

		PiZeroEnergyPeak->Draw("ep");
		TF1 * EFunc = new TF1("EFunc","x",EMin,EMax);
		EFunc->SetLineStyle(2);
		EFunc->SetLineWidth(2);
		EFunc->SetLineColor(2);
		EFunc->Draw("SAME");
		c->SaveAs(Form("PiZeroRECOEffPlots/%s/SanityCheck/PiZeroEnergyPeak.png",AbsName.Data()));


		//cout << "PiZeroRECOEffHis->Integral() = "  <<  PiZeroRECOEffHis->Integral() << endl;
		PiEffEmpty->Draw();

		SashaEff->Draw("ACSAME");
		SashaEff->SetMaximum(1.5);
		PiZeroRECOEffHis->Draw("epSAME");	



		TLegend *leg = new TLegend(0.12,0.75,0.65,0.90,NULL,"brNDC");
		leg->SetBorderSize(0);
		leg->SetTextSize(0.040);
		leg->SetTextFont(42);
		leg->SetFillStyle(0);
		leg->SetLineWidth(3);

		leg->AddEntry(PiZeroRECOEffHis,"PHENIX EMCAL at z = 2.9 m","LP");
		leg->AddEntry(SashaEff,"2.5 #times 2.5 cm^{2} at z = 3 m (A. Bazilevsky)","L");

		leg->Draw("SAME");
		lat->DrawLatex(0.40,0.50,Form("%s",EMCALName.Data()));				
		lat->DrawLatex(0.40,0.45,Form("Clustered Energy > %.0f MeV",ClusEMin * 1000));	
		lat->DrawLatex(0.40,0.40,Form("cos(#alpha) = |E_{1} - E_{2}|/|E_{1} + E_{2}| < %.1f",AngleTruthMax));	
	
		
		c->SaveAs(Form("PiZeroRECOEffPlots/%s/FinalResults/PiZeroRECOEffComp.png",AbsName.Data()));

		PiZeroRECOEffHis->Draw("ep");
		lat->DrawLatex(0.35,0.50,Form("Truth Photon Eneryg > %.0f MeV",ClusEMinTruth * 1000));	
		lat->DrawLatex(0.35,0.50,Form("%s",EMCALName.Data()));		
		lat->DrawLatex(0.35,0.45,Form("Clustered Energy > %.0f MeV",ClusEMin * 1000));	
		//lat->DrawLatex(0.35,0.40,Form("cos(#alpha) = |E_{1} - E_{2}|/|E_{1} + E_{2}| < %.1f",AngleTruthMax));	

		c->SaveAs(Form("PiZeroRECOEffPlots/%s/FinalResults/PiZeroRECOEffHis.png",AbsName.Data()));


		c->SaveAs(Form("PiZeroRECOEffPlots/%s/FinalResults/PiZeroMergeProbComparison.png",AbsName.Data()));


		PiZeroAcc->Draw("ep");
		c->SaveAs(Form("PiZeroRECOEffPlots/%s/FinalResults/PiZeroAcc.png",AbsName.Data()));

		cout << "TotalEvents = " << TotalEvents << endl;


	    TFile * fout = new TFile(Form("PiZeroMerge/%s.root",AbsName.Data()),"RECREATE");
	    fout->cd();

		PiZeroRECOEffHis->Write();
		SashaEff->Write();
		fout->Close();
		   
	
	}
