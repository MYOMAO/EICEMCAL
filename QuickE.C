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



void QuickE(){

	
	gStyle->SetOptStat(0);
	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();


	TLatex *lat = new TLatex();
	lat->SetNDC();
	lat->SetTextSize(0.035);

	TFile * fin = new TFile("G4EICDetector.root_g4femc_eval.root");
	fin->cd();

	TTree * ntp_gshower = (TTree *) fin->Get("ntp_gshower");

	TH1D * ERECOHis = new TH1D("ERECOHis","",100,0.5,1.3);

	ntp_gshower->Project("ERECOHis","e");


	ERECOHis->GetXaxis()->SetTitle("Reconstructed #gamma Energy (GeV)");
	ERECOHis->GetYaxis()->SetTitle("Number of Events");
	ERECOHis->GetXaxis()->CenterTitle();
	ERECOHis->GetYaxis()->CenterTitle();
	ERECOHis->GetYaxis()->SetTitleOffset(1.4);

	ERECOHis->Draw();
	
	
	lat->DrawLatex(0.30,0.80,"Incident #gamma Eneregy = 1 GeV");	
	lat->DrawLatex(0.30,0.75,"Standard PHENIX Pb Shashlik EMCAL");	
	
	c->SaveAs("RECOESinglePhoton.png");

}
