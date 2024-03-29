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
	gROOT->ProcessLine("a->run(-1, 444, 172, 177)");
	gROOT->ProcessLine("delete a");
	UpgradeAnalysis_12 b;
	b.OpenWithList("inputFiles_179828_v5.txt");
	gROOT->ProcessLine("b->run(-1, 444, 179, 189)");
	gROOT->ProcessLine("delete b");
	UpgradeAnalysis_12 c;
	c.OpenWithList("inputFiles_179828_v5.txt");
	gROOT->ProcessLine("c->run(-1, 444, 191, 201)");
	gROOT->ProcessLine("delete c");
	UpgradeAnalysis_12 d;
	d.OpenWithList("inputFiles_179828_v5.txt");
	gROOT->ProcessLine("d->run(-1, 444, 204, 214)");
	gROOT->ProcessLine("delete d");
	UpgradeAnalysis_12 e;
	e.OpenWithList("inputFiles_179828_v5.txt");
	gROOT->ProcessLine("e->run(-1, 444, 216, 226)");
	gROOT->ProcessLine("delete e");
	UpgradeAnalysis_12 f;
	f.OpenWithList("inputFiles_179828_v5.txt");
	gROOT->ProcessLine("f->run(-1, 444, 230, 239)");
	gROOT->ProcessLine("delete f");
	UpgradeAnalysis_12 g;
	g.OpenWithList("inputFiles_179828_v5.txt");
	gROOT->ProcessLine("g->run(-1, 444, 241, 244)");
	gROOT->ProcessLine("delete g");
	UpgradeAnalysis_12 h;
	h.OpenWithList("inputFiles_179828_v5.txt");
	gROOT->ProcessLine("h->run(-1, 444, 248, 250)");
	gROOT->ProcessLine("delete h");
	UpgradeAnalysis_12 i;
	i.OpenWithList("inputFiles_179828_v5.txt");
	gROOT->ProcessLine("i->run(-1, 444, 254, 273)");
	gROOT->ProcessLine("delete i");
	UpgradeAnalysis_12 j;
	j.OpenWithList("inputFiles_179828_v5.txt");
	gROOT->ProcessLine("j->run(-1, 444, 274, 293)");
	gROOT->ProcessLine("delete j");
	UpgradeAnalysis_12 k;
	k.OpenWithList("inputFiles_179828_v5.txt");
	gROOT->ProcessLine("k->run(-1, 444, 294, 313)");
	gROOT->ProcessLine("delete k");
	UpgradeAnalysis_12 l;
	l.OpenWithList("inputFiles_179828_v5.txt");
	gROOT->ProcessLine("l->run(-1, 444, 314, 320)");
	gROOT->ProcessLine("delete l");
	UpgradeAnalysis_12 m;
	m.OpenWithList("inputFiles_179828_v5.txt");
	gROOT->ProcessLine("m->run(-1, 444, 323, 333)");
	gROOT->ProcessLine("delete m");
	UpgradeAnalysis_12 n;
	n.OpenWithList("inputFiles_179828_v5.txt");
	gROOT->ProcessLine("n->run(-1, 444, 337, 352)");
	gROOT->ProcessLine("delete n");
	UpgradeAnalysis_12 o;
	o.OpenWithList("inputFiles_179828_v5.txt");
	gROOT->ProcessLine("o->run(-1, 444, 353, 365)");
	gROOT->ProcessLine("delete o");
	UpgradeAnalysis_12 p;
	p.OpenWithList("inputFiles_179828_v5.txt");
	gROOT->ProcessLine("p->run(-1, 444, 367, 388)");
	gROOT->ProcessLine("delete p");
	UpgradeAnalysis_12 q;
	q.OpenWithList("inputFiles_179828_v5.txt");
	gROOT->ProcessLine("q->run(-1, 444, 389, 401)");
	gROOT->ProcessLine("delete q");
	
}
