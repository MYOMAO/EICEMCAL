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


void PlotAll(){

	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();
	c->SetLogy();

	int RadLength = 20;
	double ScinThick = 0.15;

	const int NFiles = 5;
	const int NEnergy = 4;

	int FileUsed = 3;
	const int NSel = 2;

	int index;


	TLatex *lat = new TLatex();
	lat->SetNDC();
	lat->SetTextSize(0.035);


	TString infile[NFiles] = {"outFile/PlotStudy_WCu_20.root","outFile/PlotStudy_Pb_20.root","outFile/PlotStudy_Pb_18.root","outFile/PlotStudy_WCu_18.root","outFile/PlotStudy_WCu_22.root"};
	//TString infile[NFiles] = {"outFile/PlotStudy_WCu_20.root","outFile/PlotStudy_Pb_20.root"}; 




	int RadLengthAll[NFiles] = {20,20,18,18,22};
	int Color[NFiles] = {1,2,3,4,6};


	TString AbsMaterialName[NFiles] = {"WCu","Pb","PHENIX Pb","WCu","WCu"};


	int Energy[NEnergy] = {1,5, 10, 20};


	TLegend* leg = new TLegend(0.30,0.35,0.75,0.85,NULL,"brNDC");
	leg->SetBorderSize(0);
	leg->SetTextSize(0.038);
	leg->SetTextFont(42);
	leg->SetFillStyle(0);


	double TotalEnergy;
	double CumulativeEnergy;
	double MolireRadius;
	int done;


	for(int i = 0; i < NEnergy; i++){

		for(int j = 0; j < FileUsed; j++){

			done = 0;
			index = j * NEnergy + i + 1;

			TFile * fin = new TFile(infile[j].Data());
			TString TransHisName = Form("EDepTransWAbsScale_%d",Energy[i]);

			TH1D * EDepTransWAbsScale =  (TH1D *) fin->Get(TransHisName.Data());
			EDepTransWAbsScale->SetLineColor(index);

			TotalEnergy = EDepTransWAbsScale->Integral();

			for(int q = 0; q < EDepTransWAbsScale->GetNbinsX(); q++){
				if(done == 0){	
					CumulativeEnergy = EDepTransWAbsScale->Integral(1,q+1);

					if(CumulativeEnergy > TotalEnergy * 0.90){

						MolireRadius = EDepTransWAbsScale->GetXaxis()->GetBinCenter(q+1);
						done = 1;
						cout << "AbsMaterialName = " << AbsMaterialName[j].Data()  << "    MolireRadius = "  << MolireRadius  << endl;
					}
				}
			}


			if(j == 0 && i == 0){
				EDepTransWAbsScale->Draw("hist");
				EDepTransWAbsScale->SetMaximum(20);
			}
			else EDepTransWAbsScale->Draw("histSAME");

			leg->AddEntry(EDepTransWAbsScale,Form("%s EMCAL: %d GeV",AbsMaterialName[j].Data(),Energy[i]),"L");	


		}

	}

	//lat->DrawLatex(0.38,0.55,Form("Absorber Material: WCu and Pb"));			
	lat->DrawLatex(0.38,0.30,"Sampling Thickness: 1/4 X0");	
	lat->DrawLatex(0.38,0.25,Form("Total RadLength: %d X0",RadLength));
//	lat->DrawLatex(0.38,0.20,Form("Scintillator: %.1f mm Polystyrene",ScinThick*10));	

	leg->Draw("SAME");

	c->SaveAs("PlotAll/dE_dR.png");

	c->cd();


	TLegend* legNorm2 = new TLegend(0.25,0.55,0.75,0.85,NULL,"brNDC");
	legNorm2->SetBorderSize(0);
	legNorm2->SetTextSize(0.038);
	legNorm2->SetTextFont(42);
	legNorm2->SetFillStyle(0);


	for(int i = 1; i < 2; i++){

		for(int j = 0; j < NFiles; j++){
			


			index = Color[j];

			TFile * fin = new TFile(infile[j].Data());
			TString TransHisNameNorm = Form("EDepTransWAbsScaleNorm_%d",Energy[i]);

			TH1D * EDepTransWAbsScaleNorm =  (TH1D *) fin->Get(TransHisNameNorm.Data());
			EDepTransWAbsScaleNorm->SetLineColor(index);



			if(j == 0 && i == 1){
				EDepTransWAbsScaleNorm->Draw("hist");
				EDepTransWAbsScaleNorm->SetMaximum(0.4);
			}
			else EDepTransWAbsScaleNorm->Draw("histSAME");

			legNorm2->AddEntry(EDepTransWAbsScaleNorm,Form("%s EMCAL: RadLength %d X0",AbsMaterialName[j].Data(),RadLengthAll[j]),"L");	


		}

	}

	//lat->DrawLatex(0.38,0.55,Form("Absorber Material: WCu and Pb"));			
	lat->DrawLatex(0.24,0.30,"Incident Electron Energy: 5 GeV");	
	lat->DrawLatex(0.14,0.25,"Sampling Thickness: 1/4 X0");	
	
//	lat->DrawLatex(0.38,0.55,Form("Total RadLength: %d X0",RadLength));
//	lat->DrawLatex(0.38,0.50,Form("Scintillator: %.1f mm Polystyrene",ScinThick*10));	

	legNorm2->Draw("SAME");

	c->SaveAs("PlotAll/dE_dR_Norm_All.png");
	
	


	//Inclusive in Unit of Moliere R//
	TCanvas * cTransInc = new TCanvas("cTransInc","cTransInc",600,600);

	cTransInc->cd();




	TLegend* legInc = new TLegend(0.18,0.77,0.75,0.90,NULL,"brNDC");
	legInc->SetBorderSize(0);
	legInc->SetTextSize(0.035);
	legInc->SetTextFont(42);
	legInc->SetFillStyle(0);



	int FileSel[NSel] =  {0,2};

	for(int i = 1; i < 2; i++){

		for(int j = 0; j < NSel; j++){
			


			index = Color[j];

			TFile * fin = new TFile(infile[FileSel[j]].Data());
			TString TransHisNameInc = Form("EDepTransInc_%d",Energy[i]);

			TH1D * EDepTransInc =  (TH1D *) fin->Get(TransHisNameInc.Data());
			EDepTransInc->SetLineColor(index);



			if(j == 0 && i == 1){
				EDepTransInc->GetYaxis()->SetTitleOffset(1.25);
				EDepTransInc->Draw("hist");
				EDepTransInc->SetMaximum(120);
				EDepTransInc->SetMinimum(0);
			
			}
			else EDepTransInc->Draw("histSAME");

			legInc->AddEntry(EDepTransInc,Form("%s EMCAL: RadLength %d X0",AbsMaterialName[FileSel[j]].Data(),RadLengthAll[FileSel[j]]),"L");	


		}

	}

	//lat->DrawLatex(0.38,0.55,Form("Absorber Material: WCu and Pb"));			
	lat->DrawLatex(0.14,0.30,"Incident Electron Energy: 5 GeV");	
	lat->DrawLatex(0.14,0.25,"Sampling Thickness: 1/4 X0");	
	
//	lat->DrawLatex(0.38,0.55,Form("Total RadLength: %d X0",RadLength));
//	lat->DrawLatex(0.38,0.50,Form("Scintillator: %.1f mm Polystyrene",ScinThick*10));	

	legInc->Draw("SAME");

	cTransInc->SaveAs("PlotAll/CumEdepTransSelected.png");






	TLegend* legInc2 = new TLegend(0.18,0.65,0.75,0.90,NULL,"brNDC");
	legInc2->SetBorderSize(0);
	legInc2->SetTextSize(0.035);
	legInc2->SetTextFont(42);
	legInc2->SetFillStyle(0);


	for(int i = 1; i < 2; i++){

		for(int j = 0; j < NFiles; j++){
			


			index = Color[j];

			TFile * fin = new TFile(infile[j].Data());
			TString TransHisNameInc = Form("EDepTransInc_%d",Energy[i]);

			TH1D * EDepTransInc =  (TH1D *) fin->Get(TransHisNameInc.Data());
			EDepTransInc->SetLineColor(index);



			if(j == 0 && i == 1){
				EDepTransInc->GetYaxis()->SetTitleOffset(1.25);
				EDepTransInc->Draw("hist");
				EDepTransInc->SetMaximum(140);
			}
			else EDepTransInc->Draw("histSAME");

			legInc2->AddEntry(EDepTransInc,Form("%s EMCAL: RadLength %d X0",AbsMaterialName[j].Data(),RadLengthAll[j]),"L");	


		}

	}

	//lat->DrawLatex(0.38,0.55,Form("Absorber Material: WCu and Pb"));			
	lat->DrawLatex(0.14,0.30,"Incident Electron Energy: 5 GeV");	
	lat->DrawLatex(0.14,0.25,"Sampling Thickness: 1/4 X0");	
	
//	lat->DrawLatex(0.38,0.55,Form("Total RadLength: %d X0",RadLength));
//	lat->DrawLatex(0.38,0.50,Form("Scintillator: %.1f mm Polystyrene",ScinThick*10));	

	legInc2->Draw("SAME");

	cTransInc->SaveAs("PlotAll/CumEdepTransAll.png");
	
	
	//Cumalative All in cm//


	TLegend* legCum2 = new TLegend(0.18,0.65,0.75,0.90,NULL,"brNDC");
	legCum2->SetBorderSize(0);
	legCum2->SetTextSize(0.035);
	legCum2->SetTextFont(42);
	legCum2->SetFillStyle(0);


	for(int i = 1; i < 2; i++){

		for(int j = 0; j < NFiles; j++){
			


			index = Color[j];

			TFile * fin = new TFile(infile[j].Data());
			TString TransHisNameCum = Form("EDepTransCum_%d",Energy[i]);

			TH1D * EDepTransCum =  (TH1D *) fin->Get(TransHisNameCum.Data());
			EDepTransCum->SetLineColor(index);



			if(j == 0 && i == 1){
				EDepTransCum->GetYaxis()->SetTitleOffset(1.25);
				EDepTransCum->Draw("hist");
				EDepTransCum->SetMaximum(140);
			}
			else EDepTransCum->Draw("histSAME");

			legCum2->AddEntry(EDepTransCum,Form("%s EMCAL: RadLength %d X0",AbsMaterialName[j].Data(),RadLengthAll[j]),"L");	


		}

	}

	//lat->DrawLatex(0.38,0.55,Form("Absorber Material: WCu and Pb"));			
	lat->DrawLatex(0.24,0.30,"Incident Electron Energy: 5 GeV");	
	lat->DrawLatex(0.24,0.25,"Sampling Thickness: 1/4 X0");	
	
//	lat->DrawLatex(0.38,0.55,Form("Total RadLength: %d X0",RadLength));
//	lat->DrawLatex(0.38,0.50,Form("Scintillator: %.1f mm Polystyrene",ScinThick*10));	

	legCum2->Draw("SAME");

	cTransInc->SaveAs("PlotAll/CumEdepTransInCM.png");





	//Longitudinal//

	TLegend* legNorm3 = new TLegend(0.25,0.55,0.75,0.85,NULL,"brNDC");
	legNorm3->SetBorderSize(0);
	legNorm3->SetTextSize(0.038);
	legNorm3->SetTextFont(42);
	legNorm3->SetFillStyle(0);
	
	
	TCanvas * c12 = new TCanvas("c12","c12",600,600);
	c12->cd();

	for(int i = 1; i < 2; i++){

		for(int j = 0; j < NFiles; j++){
			


			index = Color[j];

			TFile * fin = new TFile(infile[j].Data());
			TString TransHisNameNorm = Form("EDepLongWAbsScale_%d",Energy[i]);

			TH1D * EDepTransWAbsScaleNorm =  (TH1D *) fin->Get(TransHisNameNorm.Data());
			EDepTransWAbsScaleNorm->Scale(1.0/EDepTransWAbsScaleNorm->Integral());
			EDepTransWAbsScaleNorm->SetLineColor(index);
			EDepTransWAbsScaleNorm->SetLineWidth(2);


			if(j == 0 && i == 1){
				EDepTransWAbsScaleNorm->Draw("hist");
				EDepTransWAbsScaleNorm->SetMaximum(0.14);
			}
			else EDepTransWAbsScaleNorm->Draw("histSAME");

			legNorm3->AddEntry(EDepTransWAbsScaleNorm,Form("%s EMCAL: RadLength %d X0",AbsMaterialName[j].Data(),RadLengthAll[j]),"L");	


		}

	}

	//lat->DrawLatex(0.38,0.55,Form("Absorber Material: WCu and Pb"));			
	lat->DrawLatex(0.36,0.47,"Incident Electron Energy: 5 GeV");	
	lat->DrawLatex(0.36,0.42,"Sampling Thickness: 1/4 X0");	
	
//	lat->DrawLatex(0.38,0.55,Form("Total RadLength: %d X0",RadLength));
//	lat->DrawLatex(0.38,0.50,Form("Scintillator: %.1f mm Polystyrene",ScinThick*10));	


	legNorm3->Draw("SAME");


	c12->SaveAs("PlotAll/dE_dz_Norm_All.png");

	
	c->cd();

	TLegend* legNorm = new TLegend(0.30,0.75,0.75,0.85,NULL,"brNDC");
	legNorm->SetBorderSize(0);
	legNorm->SetTextSize(0.038);
	legNorm->SetTextFont(42);
	legNorm->SetFillStyle(0);


	for(int i = 1; i < 2; i++){

		for(int j = 0; j < FileUsed; j++){

			index = j * NEnergy + i + 1;

			TFile * fin = new TFile(infile[j].Data());
			TString TransHisNameNorm = Form("EDepTransWAbsScaleNorm_%d",Energy[i]);

			TH1D * EDepTransWAbsScaleNorm =  (TH1D *) fin->Get(TransHisNameNorm.Data());
			EDepTransWAbsScaleNorm->SetLineColor(index);



			if(j == 0 && i == 1){
				EDepTransWAbsScaleNorm->Draw("hist");
				EDepTransWAbsScaleNorm->SetMaximum(0.4);
			}
			else EDepTransWAbsScaleNorm->Draw("histSAME");

			legNorm->AddEntry(EDepTransWAbsScaleNorm,Form("%s EMCAL: %d GeV",AbsMaterialName[j].Data(),Energy[i]),"L");	


		}

	}

	//lat->DrawLatex(0.38,0.55,Form("Absorber Material: WCu and Pb"));			
	lat->DrawLatex(0.38,0.60,"Sampling Thickness: 1/4 X0");	
	lat->DrawLatex(0.38,0.55,Form("Total RadLength: %d X0",RadLength));
//	lat->DrawLatex(0.38,0.50,Form("Scintillator: %.1f mm Polystyrene",ScinThick*10));	

	legNorm->Draw("SAME");

	c->SaveAs("PlotAll/dE_dR_Norm_18_20_22.png");




	TCanvas * c2 = new TCanvas("c2","c2",600,600);
	c2->cd();


	TLegend* leg2 = new TLegend(0.20,0.55,0.65,0.90,NULL,"brNDC");
	leg2->SetBorderSize(0);
	leg2->SetTextSize(0.038);
	leg2->SetTextFont(42);
	leg2->SetFillStyle(0);
	leg2->AddEntry("","#it{#bf{Fun4All-EIC}} Simulation","");
	leg2->AddEntry("","Geant4, Energy Deposition in Scintillator",""); 


	for(int j = 0; j < NFiles; j++){


		TFile * fin = new TFile(infile[j].Data());

		TH1D * EnergyReso =  (TH1D *) fin->Get("EnergyReso");
		EnergyReso->SetMarkerColor(Color[j]);
		EnergyReso->SetMarkerStyle(20+j);
		EnergyReso->SetLineColor(Color[j]);

		if(j == 0) EnergyReso->Draw("ep");
		else EnergyReso->Draw("epSAME");

		leg2->AddEntry(EnergyReso,Form("%s EMCAL %d X0",AbsMaterialName[j].Data(),RadLengthAll[j]),"EP");	


	}
	leg2->Draw("SAME");

	lat->DrawLatex(0.38,0.50,"Sampling Thickness: 1/4 X0");	
	//lat->DrawLatex(0.38,0.45,Form("Total RadLength: %d X0",RadLength));
//	lat->DrawLatex(0.38,0.40,Form("Scintillator: %.1f mm Polystyrene",ScinThick*10));	

	c2->SaveAs("PlotAll/EnergyReso.png");


	//ELeak

	TLegend* leg3 = new TLegend(0.10,0.62,0.55,0.90,NULL,"brNDC");
	leg3->SetBorderSize(0);
	leg3->SetTextSize(0.038);
	leg3->SetTextFont(42);
	leg3->SetFillStyle(0);


	for(int j = 0; j < NFiles; j++){


		TFile * fin = new TFile(infile[j].Data());

		TH1D * LeakEDepHis =  (TH1D *) fin->Get("LeakEDepHis");
		LeakEDepHis->SetMarkerColor(Color[j]);
		LeakEDepHis->SetMarkerStyle(20+j);
		LeakEDepHis->SetLineColor(Color[j]);
		LeakEDepHis->SetLineWidth(2);

		if(j == 0){
			LeakEDepHis->SetMinimum(0);
			LeakEDepHis->SetMaximum(6);
			LeakEDepHis->Draw("hist");
		}
			else LeakEDepHis->Draw("histSAME");

		leg3->AddEntry(LeakEDepHis,Form("%s EMCAL %d X0",AbsMaterialName[j].Data(),RadLengthAll[j]),"L");	


	}
	leg3->Draw("SAME");

	lat->DrawLatex(0.13,0.60,"Sampling Thickness: 1/4 X0");	
	//lat->DrawLatex(0.38,0.45,Form("Total RadLength: %d X0",RadLength));
	//lat->DrawLatex(0.10,0.60,Form("Scintillator: %.1f mm Polystyrene",ScinThick*10));	

	c2->SaveAs("PlotAll/ELeak.png");









}

