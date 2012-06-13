{
        cout << "Initialising Trigger Analysis code..." << endl;
        gROOT->ProcessLine(".X ~/trigger_studies/CMSSW_5_2_5/src/UserCode/L1TriggerDPG/macros/initL1Analysis.C");
        gROOT->ProcessLine(".X ~/trigger_studies/CMSSW_5_2_5/src/UserCode/L1TriggerDPG/macros/upgrade/style_macro.C");
        cout << "Loading Upgrade Analysis Macro" << endl;
        gROOT->ProcessLine(".L ~/trigger_studies/CMSSW_5_2_5/src/UserCode/L1TriggerDPG/macros/upgrade/UpgradeAnalysis_12_8TeV.C+");
        //gROOT->ProcessLine(".! rm threshRates_179828.txt");
        //gROOT->ProcessLine(".! touch threshRates_179828.txt");
        UpgradeAnalysis_12 a;
        a.OpenWithList("inputFiles_Min2012A.txt");
        gROOT->ProcessLine("a->run(100000, 1, 0, 420)");
//        gROOT->ProcessLine("delete a");
}
