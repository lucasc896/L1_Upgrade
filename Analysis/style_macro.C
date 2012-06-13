{

gROOT->SetStyle("Plain");

 // set the styles

 // Paper size

  gStyle->SetPaperSize(TStyle::kUSLetter);

  // Canvas

  gStyle->SetOptTitle(0);
  gStyle->SetCanvasColor     (0);
  gStyle->SetCanvasBorderSize(0);
  gStyle->SetCanvasBorderMode(0);
//gStyle->SetCanvasDefH      (600);
//gStyle->SetCanvasDefW      (600);
//gStyle->SetCanvasDefX      (10);
//gStyle->SetCanvasDefY      (10);

  // Pads

  gStyle->SetPadColor       (0);
  gStyle->SetPadBorderSize  (0);
  gStyle->SetPadBorderMode  (0);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadTopMargin   (0.08);
  gStyle->SetPadLeftMargin  (0.18);
  gStyle->SetPadRightMargin (0.05);
  gStyle->SetPadGridX       (0);
  gStyle->SetPadGridY       (0);
  gStyle->SetPadTickX       (1);
  gStyle->SetPadTickY       (1);

  // Frames

  gStyle->SetFrameFillStyle ( 0);
  gStyle->SetFrameFillColor ( 0);
  gStyle->SetFrameLineColor ( 1);
  gStyle->SetFrameLineStyle ( 0);
  gStyle->SetFrameLineWidth ( 2);
  gStyle->SetFrameBorderSize( 0);
  gStyle->SetFrameBorderMode( 0);


  // Histograms
  
  //gStyle->SetErrorX(0);
  //gStyle->SetHistFillColor(1);
  //gStyle->SetHistFillStyle(1);
  gStyle->SetHistLineColor(1);
  gStyle->SetLineColor(1);
  gStyle->SetLineWidth(1.5);
  gStyle->SetHistLineStyle(0);
  gStyle->SetHistLineWidth(1);
  gStyle->SetMarkerColor(1);
  gStyle->SetMarkerStyle(8);


  // Functions

  gStyle->SetFuncColor(2);
  gStyle->SetFuncStyle(0);
  gStyle->SetFuncWidth(2);

//   //Legends 

  gStyle->SetStatBorderSize(1);
  gStyle->SetStatFont      (42);
//gStyle->SetOptStat       (0);
gStyle->SetOptStat       (1111110);
  gStyle->SetOptFit        (1);
  gStyle->SetStatColor     (0);
//   gStyle->SetStatX         (0.93);
//   gStyle->SetStatY         (0.90);
//   gStyle->SetStatFontSize  (1.0);
//   gStyle->SetStatW         (0.2);
//   gStyle->SetStatH         (0.15);

//   // Labels, Ticks, and Titles

  gStyle->SetTickLength ( 0.015,"X");
  gStyle->SetTitleFont  (42,"xyz");
  gStyle->SetTitleSize(0.1);
  gStyle->SetTitleSize  ( 0.050,"X");
  gStyle->SetTitleOffset( 1.200,"X");
  gStyle->SetLabelOffset( 0.012,"X");
  gStyle->SetLabelSize  ( 0.04,"X");
  gStyle->SetStatFont   (42);
  gStyle->SetLabelFont  (42,"xyz");

  gStyle->SetTickLength ( 0.015,"Y");
  gStyle->SetTitleSize  ( 0.050,"Y");
  gStyle->SetTitleOffset( 1.200,"Y");
  gStyle->SetTitleBorderSize(0);
  gStyle->SetLabelOffset( 0.012,"Y");
  gStyle->SetLabelSize  ( 0.04,"Y");
  gStyle->SetTitleColor  (1);

  // Options


  gROOT->ForceStyle();

//TBrowser b;
}