#include "L1Ntuple.h"
#include "hist.C"
#include "Style.C"
#include "TEfficiency.h"
#include "TFile.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TString.h"
#include <sstream>
#include <fstream>
#include <vector>
#include <iostream>
#include <algorithm>
#include <math.h>

// --------------------------------------------------------------------
//                       UpgradeAnalysis_12 macro definition
// --------------------------------------------------------------------

struct order_gt : public std::binary_function<double, double, bool> {
  bool operator()(const double& x, const double& y) {
    return ( x > y ) ;
  }
};


class UpgradeAnalysis_12 : public L1Ntuple
{
  public :

    //constructor    
    UpgradeAnalysis_12(std::string filename) : L1Ntuple(filename) {}
    UpgradeAnalysis_12() {}
    ~UpgradeAnalysis_12() {}
    
    void run(Long64_t nevents, int preScale, int lsStart, int lsFin);
    void DoRateCalc(TH1D* h1, TH1D *h2, int preScale, int nLumis);
    void makeRateVsPuPlots();
    void setChrisStyle();
    

  private : 
	
	void BookHistos();
	
	TH1D *met_hist;
	TH1D *met_rate;
	TH1D *ett_hist;
	TH1D *ett_rate;
	TH1D *mht_hist;
	TH1D *mht_rate;
	TH1D *htt_hist;
	TH1D *htt_rate;
	TH1D *tau_hist[2];
	TH1D *tau_rate[2];
	TH1D *isoEG_hist[4];
	TH1D *isoEG_rate[4];
	TH1D *muon_hist[4];
	TH1D *muon_rate[4];
	TH1D *muon_hist_hi[4];
	TH1D *muon_rate_hi[4];
	TH1D *combEG_hist[4];
	TH1D *combEG_rate[4];
	TH1D *combJets_hist[4];
	TH1D *combJets_rate[4];
	TH1D *allJets_hist;
	
	TH1I *lumi_hist;
	TH1I *lumi_distro;
	TH1I *muonQual_hist;
	
	TGraph *lumiPU;
	TGraph *lumiIntL;
	TGraph *lumiInstL;
	
	TH1D *BXnoZero;
	TH1D *BXZero;
	
	TGraphErrors *puSinJ;
	TGraphErrors *puQuadJ;
	TGraphErrors *puSinEG;
	TGraphErrors *puTripEG;
	TGraphErrors *puETM;
	TGraphErrors *puETT;
	TGraphErrors *puHTM;
	TGraphErrors *puHTT;
	TGraphErrors *puSinMu;
	TGraphErrors *puTripMu;

};

void UpgradeAnalysis_12::setChrisStyle()
{
	//gStyle->SetOptStats("oeRMen");
}	


// --------------------------------------------------------------------
//                         bookhistos function 
// --------------------------------------------------------------------
void UpgradeAnalysis_12::BookHistos()
{

   met_hist 	= 	new TH1D("met_distro", "met_distro", 500, 0., 500.);
   met_rate 	= 	new TH1D("ETM_rate", "ETM", 500, 0., 500.);
   ett_hist 	= 	new TH1D("ett_distro", "ett_distro", 500, 0., 500.);
   ett_rate 	= 	new TH1D("ETT_rate", "ETT", 500, 0., 500.);
   mht_hist 	= 	new TH1D("mht_distro", "mht_distro", 500, 0., 500.);
   mht_rate 	= 	new TH1D("HTM_rate", "HTM", 500, 0., 500.);
   htt_hist 	= 	new TH1D("ht_distro", "ht_distro", 500, 0., 500.);
   htt_rate 	= 	new TH1D("HTT_rate", "HTT", 500, 0., 500.);
   
   tau_hist[0]	= 	new TH1D("tau_1_distro", "tau_1_distro", 255, 0., 255.);
   tau_hist[1]	= 	new TH1D("tau_2_distro", "tau_2_distro", 255, 0., 255.);
   tau_rate[0]	=	new TH1D("SingleTau_rate", "Single Tau", 255, 0., 255.);
   tau_rate[1]	=	new TH1D("DiTau_rate", "Double Tau", 255, 0., 255.);
   
   isoEG_hist[0]	=	new TH1D("isoEG_1_distro", "isoEG_1_distro", 63, 0., 63.);
   isoEG_hist[1]	=	new TH1D("isoEG_2_distro", "isoEG_2_distro", 63, 0., 63.);
   isoEG_hist[2]	=	new TH1D("isoEG_3_distro", "isoEG_3_distro", 63, 0., 63.);
   isoEG_hist[3]	=	new TH1D("isoEG_4_distro", "isoEG_4_distro", 63, 0., 63.);
   isoEG_rate[0]	=	new TH1D("SingleIsoEG_rate", "Single IsoEG", 63, 0., 63.);
   isoEG_rate[1]	=	new TH1D("DoubleIsoEG_rate", "Double IsoEG", 63, 0., 63.);
   isoEG_rate[2]	=	new TH1D("TripleIsoEG_rate", "Triple IsoEG", 63, 0., 63.);
   isoEG_rate[3]	=	new TH1D("QuadIsoEG_rate", "Quad IsoEG", 63, 0., 63.);
   
   muon_hist[0]	=	new TH1D("muon_1_distro", "muon_1_distro", 140, 0., 140.);
   muon_hist[1]	=	new TH1D("muon_2_distro", "muon_2_distro", 140, 0., 140.);
   muon_hist[2]	=	new TH1D("muon_3_distro", "muon_3_distro", 140, 0., 140.);
   muon_hist[3]	=	new TH1D("muon_4_distro", "muon_4_distro", 140, 0., 140.);
   muon_rate[0]	=	new TH1D("SingleMu_rate", "Single Muon", 140, 0., 140.);
   muon_rate[1]	=	new TH1D("DoubleMu_rate", "Double Muon", 140, 0., 140.);
   muon_rate[2]	=	new TH1D("TripleMu_rate", "Triple Muon", 140, 0., 140.);
   muon_rate[3]	=	new TH1D("QuadMu_rate", "Quad Muon", 140, 0., 140.);
   
   muon_hist_hi[0]	=	new TH1D("muon_1_distro_hi", "muon_1_distro_hi", 140, 0., 140.);
   muon_hist_hi[1]	=	new TH1D("muon_2_distro_hi", "muon_2_distro_hi", 140, 0., 140.);
   muon_hist_hi[2]	=	new TH1D("muon_3_distro_hi", "muon_3_distro_hi", 140, 0., 140.);
   muon_hist_hi[3]	=	new TH1D("muon_4_distro_hi", "muon_4_distro_hi", 140, 0., 140.);
   muon_rate_hi[0]	=	new TH1D("SingleMu_rate_hi", "Single Muon hi", 140, 0., 140.);
   muon_rate_hi[1]	=	new TH1D("DoubleMu_rate_hi", "Double Muon hi", 140, 0., 140.);
   muon_rate_hi[2]	=	new TH1D("TripleMu_rate_hi", "Triple Muon hi", 140, 0., 140.);
   muon_rate_hi[3]	=	new TH1D("QuadMu_rate_hi", "Quad Muon hi", 140, 0., 140.);


   combEG_hist[0]	=	new TH1D("combEG_1_distro", "combEG_1_distro", 63, 0., 63.);
   combEG_hist[1]	=	new TH1D("combEG_2_distro", "combEG_2_distro", 63, 0., 63.);
   combEG_hist[2]	=	new TH1D("combEG_3_distro", "combEG_3_distro", 63, 0., 63.);
   combEG_hist[3]	=	new TH1D("combEG_4_distro", "combEG_4_distro", 63, 0., 63.);
   combEG_rate[0]	=	new TH1D("SingleEG_rate", "Single EG", 63, 0., 63.);
   combEG_rate[1]	=	new TH1D("DoubleEG_rate", "Double EG", 63, 0., 63.);
   combEG_rate[2]	=	new TH1D("TripleEG_rate", "Triple EG", 63, 0., 63.);
   combEG_rate[3]	=	new TH1D("QuadEG_rate", "Quad EG", 63, 0., 63.);

   combJets_hist[0]	=	new TH1D("combJets_1_distro", "combJets_1_distro", 255, 0., 255.);
   combJets_hist[1]	=	new TH1D("combJets_2_distro", "combJets_2_distro", 255, 0., 255.);
   combJets_hist[2]	=	new TH1D("combJets_3_distro", "combJets_3_distro", 255, 0., 255.);
   combJets_hist[3]	=	new TH1D("combJets_4_distro", "combJets_4_distro", 255, 0., 255.);
   combJets_rate[0]	=	new TH1D("SingleJet_rate", "Single Jet", 255, 0., 255.);
   combJets_rate[1]	=	new TH1D("DoubleJet_rate", "Double Jet", 255, 0., 255.);
   combJets_rate[2]	=	new TH1D("TripleJet_rate", "Triple Jet", 255, 0., 255.);
   combJets_rate[3]	=	new TH1D("QuadJet_rate", "Quad Jet", 255, 0., 255.);
   
   allJets_hist		=	new TH1D("All Jets", "All Jets", 255, 0., 255.);
   
   lumi_hist		=	new TH1I("lumi", "lumis", 500, 0., 500.);
   lumi_distro		=	new TH1I("lumi_event_distro", "lumi_event_distro", 500, 0., 500.);
   
   muonQual_hist	=	new TH1I("muonQual_hist", "muonQual_hist", 500, 0., 500.);
   
   BXnoZero			=	new TH1D("BX before Zero Trigger", "BX before Zero Trigger", 4000, 0., 4000.);
   BXZero			=	new TH1D("BX after Zero Trigger", "BX after Zero Trigger", 4000, 0., 4000.);
   
}


// --------------------------------------------------------------------
//                          ratecalc function
// --------------------------------------------------------------------
void UpgradeAnalysis_12::DoRateCalc(TH1D* h1, TH1D *h2, int preScale, int nLumis) //change defintion to include prescale + lumis
{

    Int_t nbins  = h1->GetNbinsX();
    
    Double_t tLumi = 23.3570304;
   
    Double_t lumiAdjust = 2000./4.7;  //translation factor for 2e34 inst lumi
    //Double_t lumiAdjust = 1; 
     
    Double_t normalization = (1.*preScale*lumiAdjust)/(nLumis*tLumi); //check!! add in tLumi (~23.45secs)

    h2->GetXaxis()->SetTitle("Threshold (GeV)");
    h2->GetYaxis()->SetTitle("Rate (Hz)");
    //h2->Sumw2();
    for(Int_t i=1; i < nbins; ++i)
    {
      float rate = h1->Integral(i, nbins+1);
      h2->SetBinContent(i, rate*normalization);
      //if (i==1) std::cout << h2->GetTitle() << " - " << rate*normalization << std::endl;
      //std::cout << i << " " << h2->GetBinError(i) << std::endl;
    }
    
}


// --------------------------------------------------------------------
//                      RatevsPU plotting function
// --------------------------------------------------------------------
void UpgradeAnalysis_12::makeRateVsPuPlots()
{
	TCanvas c2;
	double sinJ[128], quadJ[128], sinEG[128], tripEG[128], ETM[128], ETT[128], HTT[128], HTM[128], sinMu[128], tripMu[128];
	double sinJerr[128], quadJerr[128], sinEGerr[128], tripEGerr[128], ETMerr[128], ETTerr[128], HTTerr[128], HTMerr[128], sinMuerr[128], tripMuerr[128], avgPUerr[128], avgInstLerr[128];
	double LS[450], IntL[450], PU[450], InstL[450], avgInstL[450], avgPU[450];
	double normInstL = 20*pow(10,33);
	int counter=0, nGra=0, lsStart, lsFin;
	
	for(int i=0; i<450; i++){ 
		avgPU[i]=0;
		avgPUerr[i]=0;
		avgInstL[i]=0;
		avgInstLerr[i]=0;
	}
	
	//----------------------------------------
	// get run data
	ifstream ifs( "/users/cl7359/trigger_studies/CMSSW_5_0_1/src/RecoLuminosity/LumiDB/getLumi_out_179828_v5_Corr.txt" );
	
	while(ifs){
		ifs >> LS[counter];
		ifs >> IntL[counter];
		ifs >> PU[counter];
		counter++;
	}
	for(int i=0; i<counter-1; i++) InstL[i] = (IntL[i]*pow(10,30))/23.357; //calculate instant. lumi from int. lumi for each LS
	
	//----------------------------------------
	//----------------------------------------
	
	//start of loop i
	for(int i=0; i<=24; i++){
		double instLFact = 0, instLFacterr=0;
		
		// define the various LS bins
		if(i==0){ lsStart = 2; lsFin = 22; }
		if(i>0 && i<7){ lsStart = (20*i) + 3; lsFin = lsStart+19; }
		if(i==7) { lsStart = 143; lsFin = 165; }
		if(i==8) { lsStart = 172; lsFin = 177; }
		if(i==9) { lsStart = 179; lsFin = 189; }
		if(i==10) { lsStart = 191; lsFin = 201; }
		if(i==11) { lsStart = 204; lsFin = 214; }
		if(i==12) { lsStart = 216; lsFin = 226; }
		if(i==13) { lsStart = 230; lsFin = 239; }
		if(i==14) { lsStart = 241; lsFin = 244; }
		if(i==15) { lsStart = 248; lsFin = 250; }
		if(i==16) { lsStart = 254; lsFin = 273; }
		if(i==17) { lsStart = 274; lsFin = 293; }
		if(i==18) { lsStart = 294; lsFin = 313; }
		if(i==19) { lsStart = 314; lsFin = 320; }
		if(i==20) { lsStart = 323; lsFin = 333; }
		if(i==21) { lsStart = 337; lsFin = 352; }
		if(i==22) { lsStart = 353; lsFin = 365; }
		if(i==23) { lsStart = 367; lsFin = 388; }
		if(i==24) { lsStart = 389; lsFin = 401; }
		
		int nLS;
		nLS = lsFin - lsStart;
		
		// skip certain lumiSection bins
		if(i==8) continue;
		if(i==13) continue;
		//if((i>10) && (i<16)) continue;
		if(i>18) continue;
		
		//calculate averages
		for(int j=lsStart; j<lsFin; j++){
			avgPU[i] += PU[j];
			avgInstL[i] += InstL[j];
		}
		avgPU[i] /= nLS;
		avgInstL[i] /= nLS;
		
		// calculate std-devs for errors
		for(int k=lsStart; k<lsFin; k++){
			avgPUerr[i] += pow( (avgPU[i]-PU[k]), 2 );
			avgInstLerr[i] += pow( (avgInstL[i]-InstL[k]), 2 );
		}
		avgPUerr[i] = sqrt( avgPUerr[i]/nLS );
		avgInstLerr[i] = sqrt( avgInstLerr[i]/nLS );
		
		//std::cout << i << "\t" << avgPU[i] << "\t" << avgPUerr[i] << std::endl;
		
		//calculate normalisation factor
		instLFact = normInstL/avgInstL[i];
		// calculate instLumi error
		instLFacterr = instLFact * (avgInstLerr[i]/avgInstL[i]);
		
		// open each input root file
		std::ostringstream inFileName;
		inFileName << "out_files/output_" << lsStart << "-" << lsFin << ".root";
		TFile *inFile = TFile::Open(inFileName.str().c_str()); //open the file
		
		
		//get rate plots from root file
		inFile->GetObject("SingleJet_rate;1", combJets_rate[0]);
		inFile->GetObject("QuadJet_rate;1", combJets_rate[3]);
		inFile->GetObject("SingleEG_rate;1", combEG_rate[0]);
		inFile->GetObject("TripleEG_rate;1", combEG_rate[2]);
		inFile->GetObject("ETM_rate;1", met_rate);
		inFile->GetObject("ETT_rate;1", ett_rate);
		inFile->GetObject("HTM_rate;1", mht_rate);
		inFile->GetObject("HTT_rate;1", htt_rate);
		inFile->GetObject("SingleMu_rate;1", muon_rate[0]);
		inFile->GetObject("TripleMu_rate;1", muon_rate[2]);
		
		
		// get values from rate plots
		sinJ[i] 		= combJets_rate[0]->GetBinContent(92);
		sinJerr[i]		= combJets_rate[0]->GetBinError(92);
		quadJ[i] 		= combJets_rate[3]->GetBinContent(36);
		quadJerr[i] 	= combJets_rate[3]->GetBinError(36);
		sinEG[i]		= combEG_rate[0]->GetBinContent(24);
		sinEGerr[i]		= combEG_rate[0]->GetBinError(24);
		tripEG[i]		= combEG_rate[2]->GetBinContent(7);
		tripEGerr[i]	= combEG_rate[2]->GetBinError(7);
		ETM[i]			= met_rate->GetBinContent(50);
		ETMerr[i]		= met_rate->GetBinError(50);
		ETT[i]			= ett_rate->GetBinContent(300);
		ETTerr[i]		= ett_rate->GetBinError(300);
		HTM[i]			= mht_rate->GetBinContent(50);
		HTMerr[i]		= mht_rate->GetBinError(50);
		HTT[i]			= htt_rate->GetBinContent(150);
		HTTerr[i]		= htt_rate->GetBinError(150);
		sinMu[i]		= muon_rate[0]->GetBinContent(16);
		sinMuerr[i]		= muon_rate[0]->GetBinError(16);
		tripMu[i]		= muon_rate[2]->GetBinContent(5);
		tripMuerr[i]	= muon_rate[2]->GetBinError(5);
		
		//std::cout << instLFact << "\t" << instLFacterr << std::endl;
		
		// calculate rate-errors
		sinJerr[i]		= sinJ[i]*instLFact*sqrt( pow( (sinJerr[i]/sinJ[i]), 2 ) + pow( (instLFacterr/instLFact), 2) );
		quadJerr[i]		= quadJ[i]*instLFact*sqrt( pow( quadJerr[i]/quadJ[i], 2 ) + pow( instLFacterr/instLFact, 2) );
		sinEGerr[i]		= sinEG[i]*instLFact*sqrt( pow( sinEGerr[i]/sinEG[i], 2 ) + pow( instLFacterr/instLFact, 2) );
		tripEGerr[i]	= tripEG[i]*instLFact*sqrt( pow( tripEGerr[i]/tripEG[i], 2 ) + pow( instLFacterr/instLFact, 2) );
		ETMerr[i]		= ETM[i]*instLFact*sqrt( pow( ETMerr[i]/ETM[i], 2 ) + pow( instLFacterr/instLFact, 2) );
		ETTerr[i]		= ETT[i]*instLFact*sqrt( pow( ETTerr[i]/ETT[i], 2 ) + pow( instLFacterr/instLFact, 2) );
		HTMerr[i]		= HTM[i]*instLFact*sqrt( pow( HTMerr[i]/HTM[i], 2 ) + pow( instLFacterr/instLFact, 2) );
		HTTerr[i]		= HTT[i]*instLFact*sqrt( pow( HTTerr[i]/HTT[i], 2 ) + pow( instLFacterr/instLFact, 2) );
		sinMuerr[i]		= sinMu[i]*instLFact*sqrt( pow( sinMuerr[i]/sinMu[i], 2 ) + pow( instLFacterr/instLFact, 2) );
		tripMuerr[i]	= tripMu[i]*instLFact*sqrt( pow( tripMuerr[i]/tripMu[i], 2 ) + pow( instLFacterr/instLFact, 2) );
		
		std::cout << i << "\t Before:" << sinJerr[i] << std::endl;
		
		// apply lumi-normalisation factors
		sinJ[i]		*= instLFact;
		quadJ[i]	*= instLFact;
		sinEG[i]	*= instLFact;
		tripEG[i]	*= instLFact;
		ETM[i]		*= instLFact;
		ETT[i]		*= instLFact;
		HTM[i]		*= instLFact;
		HTT[i]		*= instLFact;
		sinMu[i]	*= instLFact;
		tripMu[i]	*= instLFact;
		
		std::cout << i << "\t After:" << sinJerr[i] << std::endl;
		
		nGra++;
		
		
		inFile->Close();
	}   
	//end of loop i
	
	nGra = 40; //force into a re-sizeable scale
	puSinJ 		= new TGraphErrors(nGra, avgPU, sinJ, avgPUerr, sinJerr);
	puQuadJ		= new TGraphErrors(nGra, avgPU, quadJ, avgPUerr, quadJerr);
	puSinEG		= new TGraphErrors(nGra, avgPU, sinEG, avgPUerr, sinEGerr);
	puTripEG	= new TGraphErrors(nGra, avgPU, tripEG, avgPUerr, tripEGerr);
	puETM		= new TGraphErrors(nGra, avgPU, ETM, avgPUerr, ETMerr);
	puETT		= new TGraphErrors(nGra, avgPU, ETT, avgPUerr, ETTerr);
	puHTM		= new TGraphErrors(nGra, avgPU, HTM, avgPUerr, HTMerr);
	puHTT		= new TGraphErrors(nGra, avgPU, HTT, avgPUerr, HTTerr);
	puSinMu		= new TGraphErrors(nGra, avgPU, sinMu, avgPUerr, sinMuerr);
	puTripMu 	= new TGraphErrors(nGra, avgPU, tripMu, avgPUerr, tripMuerr);
	
	
	TMultiGraph *mg1 = new TMultiGraph();
	TMultiGraph *mg2 = new TMultiGraph();
	TMultiGraph *mg3 = new TMultiGraph();
	TMultiGraph *mg4 = new TMultiGraph();
	TMultiGraph *mg5 = new TMultiGraph();
	mg1->SetTitle("Rate vs PU;Average Pileup; Rate (Hz)");
	mg2->SetTitle("Rate vs PU;Average Pileup; Rate (Hz)");
	mg3->SetTitle("Rate vs PU;Average Pileup; Rate (Hz)");
	mg4->SetTitle("Rate vs PU;Average Pileup; Rate (Hz)");
	mg5->SetTitle("Rate vs PU;Average Pileup; Rate (Hz)");
	
	puSinJ	->SetMarkerColor(1);
	puSinJ	->SetMarkerStyle(8);
	puSinJ	->SetMarkerSize(0.65);
	mg1		->Add(puSinJ);
	puQuadJ	->SetMarkerColor(2);
	puQuadJ	->SetMarkerStyle(8);
	puQuadJ	->SetMarkerSize(0.65);
	mg1		->Add(puQuadJ);
	puSinEG	->SetMarkerColor(3);
	puSinEG	->SetMarkerStyle(8);
	puSinEG	->SetMarkerSize(0.65);
	mg2		->Add(puSinEG);
	puTripEG->SetMarkerColor(4);
	puTripEG->SetMarkerStyle(8);
	puTripEG->SetMarkerSize(0.65);
	mg2		->Add(puTripEG);
	puETM	->SetMarkerColor(6);
	puETM	->SetMarkerStyle(8);
	puETM	->SetMarkerSize(0.65);
	mg3		->Add(puETM);
	puETT	->SetMarkerColor(3);
	puETT	->SetMarkerStyle(24);
	puETT	->SetMarkerSize(0.65);
	mg3		->Add(puETT);
	puHTM	->SetMarkerColor(2);
	puHTM	->SetMarkerStyle(24);
	puHTM	->SetMarkerSize(0.65);
	mg4		->Add(puHTM);
	puHTT	->SetMarkerColor(1);
	puHTT	->SetMarkerStyle(24);
	puHTT	->SetMarkerSize(0.65);
	mg4		->Add(puHTT);
	/*puSinMu->SetMarkerColor(4);
	puSinMu	->SetMarkerStyle(24);
	mg5		->Add(puSinMu);
	puTripMu	->SetMarkerColor(5);
	puTripMu	->SetMarkerStyle(24);
	mg5		->Add(puTripMu);*/
	mg5		->Add(puHTT);
	mg5		->Add(puETT);
	
	
	
	mg1->Draw("ap");
	mg1->GetXaxis()->SetRangeUser(7,35);
	mg1->GetYaxis()->SetRangeUser(0,800000);
	
	//TLegend *leg1 = new TLegend(0.2, 0.7, 0.4, 0.9);
	TLegend *leg1 = new TLegend(0.2, 0.78, 0.55, 0.88);
	leg1->AddEntry(puSinJ, "Single Jet - 92 GeV", "P");
	leg1->AddEntry(puQuadJ, "Quad Jet - 36 GeV", "P");
	
	leg1->SetFillColor(0);
	leg1->SetLineColor(0);
	leg1->Draw();  
	
	c2.Print("threshRate_output.pdf(");
	c2.Print("puJetRates.png");
	c2.Clear();
	
	mg2->Draw("ap");
	mg2->GetXaxis()->SetRangeUser(7,35);
	mg2->GetYaxis()->SetRangeUser(0,80000);
	
	leg1->Clear();
	leg1->AddEntry(puSinEG, "Single EG - 24 GeV", "P");
	leg1->AddEntry(puTripEG, "Triple EG - 7 GeV", "P");
	leg1->Draw();
	
	c2.Print("threshRate_output.pdf");
	c2.Print("puEGRates.png");
	c2.Clear();
	
	mg3->Draw("ap");
	mg3->GetXaxis()->SetRangeUser(7,35);
	mg3->GetYaxis()->SetRangeUser(0,10000);
	
	leg1->Clear();
	leg1->AddEntry(puETM, "ETM - 50 GeV", "P");
	leg1->AddEntry(puETT, "ETT - 300 GeV", "P");
	leg1->Draw();
	
	c2.Print("threshRate_output.pdf");
	c2.Print("puETRates.png");
	c2.Clear();
	
	mg4->Draw("ap");
	mg4->GetXaxis()->SetRangeUser(7,35);
	mg4->GetYaxis()->SetRangeUser(0,100000);
	
	leg1->Clear();
	leg1->AddEntry(puHTM, "HTM - 50 GeV", "P");
	leg1->AddEntry(puHTT, "HTT - 150 GeV", "P");
	leg1->Draw();
	
	c2.Print("threshRate_output.pdf");
	c2.Print("puHTRates.png");
	c2.Clear();
	
	mg5->Draw("ap");
	mg5->GetXaxis()->SetRangeUser(7,35);
	mg5->GetYaxis()->SetRangeUser(0,60000);
	
	leg1->Clear();
	leg1->AddEntry(puHTT, "HTT - 150 GeV", "P");
	leg1->AddEntry(puETT, "ETT - 150 GeV", "P");
	leg1->Draw();
	
	c2.Print("threshRate_output.pdf)");
	c2.Print("puEnSumRates.png");
	c2.Clear();
	
	/*leg1->AddEntry(puSinMu, "Single Mu - 16 GeV", "P");
	leg1->AddEntry(puTripMu, "Triple Mu - 5 GeV", "P");*/


}


// --------------------------------------------------------------------
//                             run function 
// --------------------------------------------------------------------
void UpgradeAnalysis_12::run(Long64_t nevents, int preScale, int lsStart, int lsFin)//, TString outFileName)
{
  
	//load TDR style
	setTDRStyle();
	gStyle->SetOptStat(0);
	
	std::vector<int>	lumis;
	std::vector<string> instLumis;

	int counter=0;
	
	
	std::ostringstream outFileName;
	
	outFileName << "out_files/output_" << lsStart << "-" << lsFin << ".root";
	
	TFile *outFile = TFile::Open(outFileName.str().c_str(), "RECREATE"); //CHANGE
	outFile->cd();
	
	BookHistos();
	
	int ntot = 0;
	bool anyGoodLumi = false;
	
	
	if (nevents==-1 || nevents>GetEntries()) nevents=GetEntries();
	std::cout << nevents << " to process ..." << std::endl;
	
	//loop over the events to fill all object vectors
	for (Long64_t i=0; i<nevents; i++)
	{
		//load the i-th event 
		Long64_t ientry = LoadTree(i); if (ientry < 0) break;
		GetEntry(i);
		BXnoZero->Fill(event_->bx);
		//process progress
		if(i!=0 && (i%1000)==0) {std::cout << "- processing event " << i << "\r" << std::flush;}
		if ((gt_->tw1[2] & 0x0001) > 0){ // check that zero bias trigger fired
			
			BXZero->Fill(event_->bx);
			
			int thisLumi = event_->lumi; //get this event's lumi
			
			if ((thisLumi>=lsStart) && (thisLumi<=lsFin)){
				anyGoodLumi = true;
				
				// create vectors for all jets and em's
				std::vector<double> 	combJets;
				std::vector<double>		combEm;
				
				ntot++; //counter for events ran over
				//avgPU += PU[thisLumi];
				
				
				   //------- Start Filling Objects --------//
				
				met_hist-> 		Fill(l1extra_->met);
				ett_hist->		Fill(l1extra_->et);
				mht_hist->		Fill(l1extra_->mht);
				htt_hist->		Fill(l1extra_->ht);
				
				//if((l1extra_->met) > 200.) std::cout << "MET OUTLIER: " << i << "\t" << event_->event << "\t" << thisLumi << "\t" << l1extra_->met << std::endl;
				//if((l1extra_->mht) > 200.) std::cout << "MHT OUTLIER: " << i << "\t" << event_->event << "\t" << thisLumi << "\t" << l1extra_->mht << std::endl;

				
				bool newLumi = true;
				lumi_distro->Fill(thisLumi);
				
				for(unsigned int i=0; i<lumis.size(); ++i){
					if(lumis.at(i)==thisLumi) newLumi = false;
				}
				if(newLumi) lumis.push_back(thisLumi); //enter unique lumis to vector
				
				
				//------------------muon--------------------
				for(unsigned int i=0; i<l1extra_->nMuons; i++){
					if(l1extra_->muonEt.at(i) > 0.){
						if(i<4){
							muon_hist[i]->Fill(l1extra_->muonEt.at(i));
						}
					}
				}
				
				for(int i=0; i<gmt_->N; i++){
					if(gmt_->Pt.at(i) > 0.){
						if(i<4){
							if((gmt_->CandBx.at(i))==0){ //check is the same BX as trigger
								if((gmt_->Qual.at(i)) >= 5.) muon_hist_hi[i]->Fill(gmt_->Pt.at(i)); //quality check
							}
						}
					}
				}
				
				//---------------combined EG----------------
				// isoEG
				for(unsigned int i=0; i<l1extra_->nIsoEm; i++){
					if(l1extra_->isoEmEt.at(i) > 0.){
						if((l1extra_->isoEmBx.at(i))!=0) break;
						double et = l1extra_->isoEmEt.at(i); 
						combEm.push_back(et);
						if (i<4) isoEG_hist[i]->Fill(et);
					}
				}
				
				// nonisoEG
				for(unsigned int i=0; i<l1extra_->nNonIsoEm; i++){
					if(l1extra_->nonIsoEmEt.at(i) > 0.){
						if((l1extra_->nonIsoEmBx.at(i))!=0) break;
						double et = l1extra_->nonIsoEmEt.at(i); 
						combEm.push_back(et);
					}
				}
				
				sort(combEm.begin(), combEm.end(), order_gt()); //sort them into descending order
				
				// fill distros for combined EG objects
				for(unsigned int i=0; i<4; i++){
					if (combEm.size()>i) combEG_hist[i]->Fill(combEm.at(i));
				}
				
				//---------------combined Jet----------------
				
				//tau-jets
				for(unsigned int i=0; i<l1extra_->nTauJets; i++){ 
					if(l1extra_->tauJetEt.at(i) > 0.){
						if((l1extra_->tauJetBx.at(i))!=0) break;
						double et = l1extra_->tauJetEt.at(i); 
						combJets.push_back(et);
						if (i<2) tau_hist[i]->Fill(et);
					}
				}
				
				//fwd-jets
				for(unsigned int i=0; i<l1extra_->nFwdJets; i++){
					if(l1extra_->fwdJetEt.at(i) > 0.){
						if((l1extra_->fwdJetBx.at(i))!=0) break;
						double et = l1extra_->fwdJetEt.at(i);
						combJets.push_back(et);
					}
				}
				
				//cen-jets
				for(unsigned int i=0; i<l1extra_->nCenJets; i++){
					if(l1extra_->cenJetEt.at(i) > 0.){
						if((l1extra_->cenJetBx.at(i))!=0) break;
						double et = l1extra_->cenJetEt.at(i);
						combJets.push_back(et);
					}
				}
				
				sort(combJets.begin(), combJets.end(), order_gt());
				
				for(unsigned int i=0; i<combJets.size(); i++) allJets_hist->Fill(combJets.at(i));
				
				for(unsigned int i=0; i<4; i++){
					if (combJets.size() > i) combJets_hist[i]->Fill(combJets.at(i));
				}
			}
		} //zerobias trigger fired
	} //nevents
	
	if (anyGoodLumi){
	
		//avgPU /= ntot; //calculate average PU for selected events
		
		//lumiPU	=	new TGraph(counter, LS, PU);
		//lumiIntL	=	new TGraph(counter, LS, IntL);
		//lumiInstL	=	new TGraph(counter, LS, InstL);
		
		
		//------- Start Rate Calculations --------//
		
		int nLumis = lumis.size(); //get number of lumis from vector
				TCanvas c1("c1");
		
		
		DoRateCalc(met_hist, met_rate, preScale, nLumis);
		
		met_rate->Draw();
		met_rate->GetYaxis()->SetRangeUser(100,100000);
		met_rate->SetLineColor(kBlack);
		
		DoRateCalc(ett_hist, ett_rate, preScale, nLumis);
		ett_rate->Draw("SAME");
		ett_rate->SetLineColor(kRed);
		
		DoRateCalc(mht_hist, mht_rate, preScale, nLumis);
		mht_rate->Draw("SAME");
		mht_rate->SetLineColor(kBlue);
		
		DoRateCalc(htt_hist, htt_rate, preScale, nLumis);
		htt_rate->Draw("SAME");
		htt_rate->SetLineColor(kGreen);
		
		TLegend *leg1 = new TLegend(0.6, 0.75, 0.8, 0.92);
		leg1->AddEntry(met_rate, met_rate->GetTitle(), "L");
		leg1->AddEntry(ett_rate, ett_rate->GetTitle(), "L");
		leg1->AddEntry(mht_rate, mht_rate->GetTitle(), "L");
		leg1->AddEntry(htt_rate, htt_rate->GetTitle(), "L");
		leg1->SetFillColor(0);
		leg1->Draw();
		
		c1.SetLogy(1);
		//c1.Print("enSumRate.png");		
		c1.Print("rate_output.pdf(");
		//c1.Print("en_lhcrates.png");
		
		for(int i=0;i<2;i++){ 
			DoRateCalc(tau_hist[i],tau_rate[i], preScale, nLumis);
		}
				tau_rate[0]->Draw();
		tau_rate[0]->GetYaxis()->SetRangeUser(100,100000);
		tau_rate[0]->SetLineColor(kBlack);
		tau_rate[1]->Draw("SAME");
		tau_rate[1]->SetLineColor(kRed);
		
		TLegend *leg2 = new TLegend(0.6, 0.82, 0.95, 0.92);
		leg2->AddEntry(tau_rate[0], tau_rate[0]->GetTitle(), "L");
		leg2->AddEntry(tau_rate[1], tau_rate[1]->GetTitle(), "L");
		leg2->SetFillColor(0);
		leg2->Draw();
				
		c1.Print("rate_output.pdf");
		//c1.Print("tau_lhcrates.png");
		
		for(int i=0;i<4;i++){ 
			DoRateCalc(isoEG_hist[i], isoEG_rate[i], preScale, nLumis);
		}
		isoEG_rate[0]->Draw();
		isoEG_rate[0]->GetYaxis()->SetRangeUser(100,100000);
		isoEG_rate[0]->SetLineColor(kBlack);
		isoEG_rate[1]->Draw("SAME");
		isoEG_rate[1]->SetLineColor(kRed);
		isoEG_rate[2]->Draw("SAME");
		isoEG_rate[2]->SetLineColor(kBlue);
		isoEG_rate[3]->Draw("SAME");
		isoEG_rate[3]->SetLineColor(kGreen+2);
		
		TLegend *leg3 = new TLegend(0.6, 0.75, 0.95, 0.92);
		leg3->Clear();
		leg3->AddEntry(isoEG_rate[0], isoEG_rate[0]->GetTitle(), "L");
		leg3->AddEntry(isoEG_rate[1], isoEG_rate[1]->GetTitle(), "L");
		leg3->AddEntry(isoEG_rate[2], isoEG_rate[2]->GetTitle(), "L");
		leg3->AddEntry(isoEG_rate[3], isoEG_rate[3]->GetTitle(), "L");
		leg3->SetFillColor(0);
		leg3->Draw();
				
		c1.Print("rate_output.pdf");
		//c1.Print("isoEG_lhcrates.png");
		
		
		for(int i=0;i<4;i++){
			DoRateCalc(combEG_hist[i], combEG_rate[i], preScale, nLumis);	
		}
		combEG_rate[0]->Draw();
		combEG_rate[0]->GetYaxis()->SetRangeUser(100,100000);
		combEG_rate[0]->SetLineColor(kBlack);
		combEG_rate[1]->Draw("SAME");
		combEG_rate[1]->SetLineColor(kRed);
		combEG_rate[2]->Draw("SAME");
		combEG_rate[2]->SetLineColor(kBlue);
		combEG_rate[3]->Draw("SAME");
		combEG_rate[3]->SetLineColor(kGreen+2);
		
		leg3->Clear();
		leg3->AddEntry(combEG_rate[0], combEG_rate[0]->GetTitle(), "L");
		leg3->AddEntry(combEG_rate[1], combEG_rate[1]->GetTitle(), "L");
		leg3->AddEntry(combEG_rate[2], combEG_rate[2]->GetTitle(), "L");
		leg3->AddEntry(combEG_rate[3], combEG_rate[3]->GetTitle(), "L");
		leg3->Draw();
				
		c1.Print("rate_output.pdf");
		//c1.Print("EG_lhcrates.png");
		
		
		for(int i=0;i<4;i++){
			DoRateCalc(muon_hist[i], muon_rate[i], preScale, nLumis);
		}
		muon_rate[0]->Draw();
		muon_rate[0]->GetYaxis()->SetRangeUser(100,100000);
		muon_rate[0]->SetLineColor(kBlack);
		muon_rate[1]->Draw("SAME");
		muon_rate[1]->SetLineColor(kRed);
		muon_rate[2]->Draw("SAME");
		muon_rate[2]->SetLineColor(kBlue);
		muon_rate[3]->Draw("SAME");
		muon_rate[3]->SetLineColor(kGreen+2);
		
		TLegend *leg4 = new TLegend(0.55, 0.25, 0.75, 0.42);
		leg4->Clear();
		leg4->AddEntry(muon_rate[0], muon_rate[0]->GetTitle(), "L");
		leg4->AddEntry(muon_rate[1], muon_rate[1]->GetTitle(), "L");
		leg4->AddEntry(muon_rate[2], muon_rate[2]->GetTitle(), "L");
		leg4->AddEntry(muon_rate[3], muon_rate[3]->GetTitle(), "L");
		leg4->SetFillColor(0);
		leg4->Draw();
				
		c1.Print("rate_output.pdf");
		c1.Print("muon_rates.png");
		
		for(int i=0;i<4;i++){
			DoRateCalc(muon_hist_hi[i], muon_rate_hi[i], preScale, nLumis);
		}

		muon_rate_hi[0]->Draw();
		muon_rate_hi[0]->GetYaxis()->SetRangeUser(100,100000);
		muon_rate_hi[0]->SetLineColor(kBlack);
		muon_rate_hi[1]->Draw("SAME");
		muon_rate_hi[1]->SetLineColor(kRed);
		muon_rate_hi[2]->Draw("SAME");
		muon_rate_hi[2]->SetLineColor(kBlue);
		muon_rate_hi[3]->Draw("SAME");
		muon_rate_hi[3]->SetLineColor(kGreen+2);

		leg4->Clear();
		leg4->AddEntry(muon_rate_hi[0], muon_rate_hi[0]->GetTitle(), "L");
		leg4->AddEntry(muon_rate_hi[1], muon_rate_hi[1]->GetTitle(), "L");
		leg4->AddEntry(muon_rate_hi[2], muon_rate_hi[2]->GetTitle(), "L");
		leg4->AddEntry(muon_rate_hi[3], muon_rate_hi[3]->GetTitle(), "L");
		leg4->SetFillColor(0);
		leg4->Draw();

		c1.Print("rate_output.pdf");
		
		for(int i=0;i<4;i++){
			DoRateCalc(combJets_hist[i], combJets_rate[i], preScale, nLumis);	
		}
		combJets_rate[0]->Draw();
		combJets_rate[0]->GetYaxis()->SetRangeUser(100,100000);
		combJets_rate[0]->SetLineColor(kBlack);
		combJets_rate[1]->Draw("SAME");
		combJets_rate[1]->SetLineColor(kRed);
		combJets_rate[2]->Draw("SAME");
		combJets_rate[2]->SetLineColor(kBlue);
		combJets_rate[3]->Draw("SAME");
		combJets_rate[3]->SetLineColor(kGreen+2);
		
		leg3->Clear();
		leg3->AddEntry(combJets_rate[0], combJets_rate[0]->GetTitle(), "L");
		leg3->AddEntry(combJets_rate[1], combJets_rate[1]->GetTitle(), "L");
		leg3->AddEntry(combJets_rate[2], combJets_rate[2]->GetTitle(), "L");
		leg3->AddEntry(combJets_rate[3], combJets_rate[3]->GetTitle(), "L");
		leg3->Draw();
				
		c1.Print("rate_output.pdf)");
		//c1.Print("jet_lhcrates.png");
		
		//		lumi_distro->Draw();
		//c1.Print("rate_output.pdf)");
		
		std::cout << std::endl << "Total number of events: " << ntot << std::endl;
		
		for(unsigned int i=0; i<lumis.size(); ++i) lumi_hist->Fill(lumis.at(i));
		
		
		//------- Find Specific Rates --------//
		
		ofstream rateFile;
		rateFile.open ("/users/cl7359/trigger_studies/CMSSW_4_4_2_patch9/src/UserCode/L1TriggerDPG/macros/upgrade/threshRates_179828.txt", fstream::app); //open file in append mode
		
		rateFile << std::endl;
		rateFile << combJets_rate[0]->GetBinContent(92) 	<< std::endl;
		rateFile << combJets_rate[3]->GetBinContent(36) 	<< std::endl;
		rateFile << combEG_rate[0]	->GetBinContent(24) 	<< std::endl;
		rateFile << combEG_rate[2]	->GetBinContent(7) 		<< std::endl;
		rateFile << met_rate		->GetBinContent(50) 	<< std::endl;
		rateFile << ett_rate		->GetBinContent(300) 	<< std::endl;
		rateFile << htt_rate		->GetBinContent(150) 	<< std::endl;
		rateFile << mht_rate		->GetBinContent(50) 	<< std::endl;
		rateFile << muon_rate[0]	->GetBinContent(16) 	<< std::endl;
		rateFile << muon_rate[2]	->GetBinContent(5) 		<< std::endl;
		
		rateFile.close();
		
		/*lumiInstL->GetXaxis()->SetTitle("LumiSection");
		lumiInstL->GetYaxis()->SetTitle("Instantaneous Lumi (cm^-2 s^-1)");
		lumiInstL->Draw("A*");
		c1.Print("rate_output.pdf(");
		
		lumiPU->GetXaxis()->SetTitle("LumiSection");
		lumiPU->GetYaxis()->SetTitle("PileUP");
		lumiPU->Draw("A*");
		c1.Print("rate_output.pdf)");
		
		lumiIntL->GetXaxis()->SetTitle("LumiSection");
		lumiIntL->GetYaxis()->SetTitle("Integrated Lumi (/ub)");
		//lumiIntL->Draw("A*");
		
		c1.Print("rate_output.pdf)");
		
		c1.cd();
		c1.Write();*/
		
		outFile->Write();
		outFile->Close();
	}
} 	
