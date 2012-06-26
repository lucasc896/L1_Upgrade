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
    
    void GetTriggerDistros(Long64_t nevents, TString outFileName);
    void setChrisStyle();
    

  private : 
	
	void BookHistos();
	
	TH1D *met_hist;
	TH1D *ett_hist;
	TH1D *mht_hist;
	TH1D *htt_hist;
	TH1D *tau_hist[2];
	TH1D *isoEG_hist[4];
	TH1D *muon_hist[4];
	TH1D *muon_hist_hi[4];
	TH1D *combEG_hist[4];
	TH1D *combJets_hist[4];
	TH1D *allJets_hist;
	
	TH1I *lumi_hist;
	TH1I *lumi_distro;
	TH1I *muonQual_hist;
	
	TGraph *lumiPU;
	TGraph *lumiIntL;
	TGraph *lumiInstL;
	
	TH1D *BXnoZero;
	TH1D *BXZero;
	
	TH1D *runHist;
	TH2D *nVtx2D;
	TH1D *nVtx1D;
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
   ett_hist 	= 	new TH1D("ett_distro", "ett_distro", 500, 0., 500.);
   mht_hist 	= 	new TH1D("mht_distro", "mht_distro", 500, 0., 500.);
   htt_hist 	= 	new TH1D("ht_distro", "ht_distro", 500, 0., 500.);
   
   tau_hist[0]	= 	new TH1D("tau_1_distro", "tau_1_distro", 255, 0., 255.);
   tau_hist[1]	= 	new TH1D("tau_2_distro", "tau_2_distro", 255, 0., 255.);
   
   isoEG_hist[0]	=	new TH1D("isoEG_1_distro", "isoEG_1_distro", 63, 0., 63.);
   isoEG_hist[1]	=	new TH1D("isoEG_2_distro", "isoEG_2_distro", 63, 0., 63.);
   isoEG_hist[2]	=	new TH1D("isoEG_3_distro", "isoEG_3_distro", 63, 0., 63.);
   isoEG_hist[3]	=	new TH1D("isoEG_4_distro", "isoEG_4_distro", 63, 0., 63.);
   
   muon_hist[0]	=	new TH1D("muon_1_distro", "muon_1_distro", 140, 0., 140.);
   muon_hist[1]	=	new TH1D("muon_2_distro", "muon_2_distro", 140, 0., 140.);
   muon_hist[2]	=	new TH1D("muon_3_distro", "muon_3_distro", 140, 0., 140.);
   muon_hist[3]	=	new TH1D("muon_4_distro", "muon_4_distro", 140, 0., 140.);
   
   muon_hist_hi[0]	=	new TH1D("muon_1_distro_hi", "muon_1_distro_hi", 140, 0., 140.);
   muon_hist_hi[1]	=	new TH1D("muon_2_distro_hi", "muon_2_distro_hi", 140, 0., 140.);
   muon_hist_hi[2]	=	new TH1D("muon_3_distro_hi", "muon_3_distro_hi", 140, 0., 140.);
   muon_hist_hi[3]	=	new TH1D("muon_4_distro_hi", "muon_4_distro_hi", 140, 0., 140.);

   combEG_hist[0]	=	new TH1D("combEG_1_distro", "combEG_1_distro", 63, 0., 63.);
   combEG_hist[1]	=	new TH1D("combEG_2_distro", "combEG_2_distro", 63, 0., 63.);
   combEG_hist[2]	=	new TH1D("combEG_3_distro", "combEG_3_distro", 63, 0., 63.);
   combEG_hist[3]	=	new TH1D("combEG_4_distro", "combEG_4_distro", 63, 0., 63.);

   combJets_hist[0]	=	new TH1D("combJets_1_distro", "combJets_1_distro", 255, 0., 255.);
   combJets_hist[1]	=	new TH1D("combJets_2_distro", "combJets_2_distro", 255, 0., 255.);
   combJets_hist[2]	=	new TH1D("combJets_3_distro", "combJets_3_distro", 255, 0., 255.);
   combJets_hist[3]	=	new TH1D("combJets_4_distro", "combJets_4_distro", 255, 0., 255.);
   
   allJets_hist		=	new TH1D("All Jets", "All Jets", 255, 0., 255.);
   
   lumi_hist		=	new TH1I("lumi", "lumis", 10000, 0., 10000.);
   lumi_distro		=	new TH1I("lumi_event_distro", "lumi_event_distro", 500, 0., 500.);
   
   muonQual_hist	=	new TH1I("muonQual_hist", "muonQual_hist", 500, 0., 500.);
   
   BXnoZero			=	new TH1D("BX before Zero Trigger", "BX before Zero Trigger", 4000, 0., 4000.);
   BXZero			=	new TH1D("BX after Zero Trigger", "BX after Zero Trigger", 4000, 0., 4000.);
   
   runHist			= 	new TH1D("Run Distro", "Run Distro", 5000, 190000., 195000.);
   nVtx2D			=	new TH2D("nVtx2D", "nVtx2D", 3291, 190450., 193741., 100, 0., 100.);
   nVtx1D			=	new TH1D("nVtx1D", "nVtx1D", 100, 0., 100.);
   
}


// --------------------------------------------------------------------
//                             run function 
// --------------------------------------------------------------------
void UpgradeAnalysis_12::GetTriggerDistros(Long64_t nevents, TString outFileName)
{
  
	//load TDR style
	setTDRStyle();
	gStyle->SetOptStat(0);
	
	std::vector<int>	lumis;
	std::vector<string> instLumis;

	int counter=0;
	
	//load the output rootfile
	TString outFileName1 = "output.root";
	TFile *outFile = TFile::Open(outFileName1, "RECREATE");
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
		//process progress
		if(i!=0 && (i%1000)==0) {std::cout << "- processing event " << i << "\r" << std::flush;}

		if (event_->run == 190949) continue;

		//if ((gt_->tw1[2] & 0x0001) > 0){ // check that zero bias trigger fired
				
			// create vectors for all jets and em's
			std::vector<double> 	combJets;
			std::vector<double>		combEm;

			ntot++; //counter for events ran over

			int thisLumi = event_->lumi; //get this event's lumi
			lumi_hist->Fill(thisLumi);

			nVtx2D->Fill(event_->run, recoVertex_->nVtx);
			nVtx1D->Fill(recoVertex_->nVtx);

			//------- Start Filling Objects --------//

			met_hist-> 		Fill(l1extra_->met);
			ett_hist->		Fill(l1extra_->et);
			mht_hist->		Fill(l1extra_->mht);
			htt_hist->		Fill(l1extra_->ht);

			bool newLumi = true;
			lumi_distro->Fill(thisLumi);

			//finds if there is a unique lumi being ran
			for(unsigned int i=0; i<lumis.size(); ++i){
				if(lumis.at(i)==thisLumi) newLumi = false;
			}
			if(newLumi) lumis.push_back(thisLumi); //enter unique lumis to vector


			//------------------muon--------------------//
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

			//---------------combined EG----------------//
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

			//---------------combined Jet----------------//

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
		//} //zerobias trigger fired
	} //nevents
	
	std::cout << ntot << " events run over" << std::endl;

	outFile->Write();
	outFile->Close();
}
