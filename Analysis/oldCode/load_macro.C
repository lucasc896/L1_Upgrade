{
	cout << "Initialising Trigger Analysis code..." << endl;
	gROOT->ProcessLine(".X ~/trigger_studies/CMSSW_4_4_2_patch9/src/UserCode/L1TriggerDPG/macros/initL1Analysis.C");
	gROOT->ProcessLine(".X ~/trigger_studies/CMSSW_4_4_2_patch9/src/UserCode/L1TriggerDPG/macros/upgrade/style_macro.C");
	cout << "Loading Upgrade Analysis Macro" << endl;
	gROOT->ProcessLine(".L ~/trigger_studies/CMSSW_4_4_2_patch9/src/UserCode/L1TriggerDPG/macros/upgrade/UpgradeAnalysis_12.C+");
	//gROOT->ProcessLine(".! rm threshRates_179828.txt");
	//gROOT->ProcessLine(".! touch threshRates_179828.txt");
	UpgradeAnalysis_12 a;
	a.OpenWithList("inputFiles_179828_v5.txt");
	gROOT->ProcessLine("a->run(-1, 444, 2, 22)");
	gROOT->ProcessLine("delete a");
	UpgradeAnalysis_12 b;
	b.OpenWithList("inputFiles_179828_v5.txt");
	gROOT->ProcessLine("b->run(-1, 444, 23, 42)");
	gROOT->ProcessLine("delete b");
	UpgradeAnalysis_12 c;
	c.OpenWithList("inputFiles_179828_v5.txt");
	gROOT->ProcessLine("c->run(-1, 444, 43, 62)");
	gROOT->ProcessLine("delete c");
	UpgradeAnalysis_12 d;
	d.OpenWithList("inputFiles_179828_v5.txt");
	gROOT->ProcessLine("d->run(-1, 444, 63, 82)");
	gROOT->ProcessLine("delete d");
	UpgradeAnalysis_12 e;
	e.OpenWithList("inputFiles_179828_v5.txt");
	gROOT->ProcessLine("e->run(-1, 444, 83, 102)");
	gROOT->ProcessLine("delete e");
	UpgradeAnalysis_12 f;
	f.OpenWithList("inputFiles_179828_v5.txt");
	gROOT->ProcessLine("f->run(-1, 444, 103, 122)");
	gROOT->ProcessLine("delete f");
	UpgradeAnalysis_12 g;
	g.OpenWithList("inputFiles_179828_v5.txt");
	gROOT->ProcessLine("g->run(-1, 444, 123, 142)");
	gROOT->ProcessLine("delete g");
	UpgradeAnalysis_12 h;
	h.OpenWithList("inputFiles_179828_v5.txt");
	gROOT->ProcessLine("h->run(-1, 444, 143, 165)");
	gROOT->ProcessLine("delete h");
	

//	gROOT->ProcessLine("a->makeRateVsPuPlots('/users/cl7359/trigger_studies/CMSSW_4_4_2_patch9/src/UserCode/L1TriggerDPG/macros/upgrade/threshRates_179828.txt')");
//	gROOT->ProcessLine("a->makeRateVsPuPlots()");
//	gROOT->ProcessLine("TBrowser asd");
}
