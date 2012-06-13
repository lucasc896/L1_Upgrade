{
	gROOT->ProcessLine(".X ~/trigger_studies/CMSSW_4_4_2_patch9/src/UserCode/L1TriggerDPG/macros/initL1Analysis.C");
	gROOT->ProcessLine(".X ~/trigger_studies/CMSSW_4_4_2_patch9/src/UserCode/L1TriggerDPG/macros/upgrade/style_macro.C");
	cout << "Loading Upgrade Analysis Macro" << endl;
	gROOT->ProcessLine(".L ~/trigger_studies/CMSSW_4_4_2_patch9/src/UserCode/L1TriggerDPG/macros/upgrade/UpgradeAnalysis_12.C+");
	//gROOT->ProcessLine(".! rm threshRates_179828.txt");
	//gROOT->ProcessLine(".! touch threshRates_179828.txt");
	UpgradeAnalysis_12 a;
	gROOT->ProcessLine("a->makeRateVsPuPlots()");
	//delete a;
	//gROOT->ProcessLine(".q");
}
