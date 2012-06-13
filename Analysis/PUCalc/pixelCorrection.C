#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <algorithm>
#include <math.h>
#include <sstream>


double correctedLumi(double rawLumi, int nBX, int correctionMode){

   // rawLumi is assumed to be the lumi for one LS given in units of ub.

   // nBX is the number of colliding bunches (typically 1317 or 1331 for 2011)
   // nBX = 10 for Run 179828 (high PU run)

   // correction mode = 0 ==> no correction
   // correction mode = 2 ==> lumiCalc2.py  no good for high pile
   // correction mode = 3 ==> correction based on High PU scan

  

   double val = rawLumi;
   if(correctionMode == 0)return val;
   if(correctionMode != 2 && correctionMode != 3)return val;  

   double gamma = 1.141;
   double beta = 1.;
   if(nBX >= 1317){
      beta = 0.9748;
   } else if (nBX >= 1179) {
      beta = 0.9768;
   } else if (nBX >= 1041) {
      beta = 0.9788;
   } else if (nBX >= 873) {
      beta = 0.9788;
   } else if (nBX >= 700) {
      beta = 0.9836;
   } else if (nBX >= 597) {
      beta = 0.9851;
   } else if (nBX >= 423) {
      beta = 0.9881;
   } else if (nBX >= 321) {
      beta = 0.9902;
   } else if (nBX >= 213) {
      beta = 0.9916;
   }

   double fOrbit = 11246.;
   double secondsPerLS = pow(2, 18) / fOrbit;
   double Ldot = rawLumi/nBX/secondsPerLS;
   double alpha = 0.076;
   val = rawLumi*gamma*beta*(1. - alpha*Ldot); 
   if(correctionMode == 2)return val;

   // from vertex
   double alpha_1 = 0.0700;
   double alpha_2 = -0.00406;
   alpha_2 = -0.0045;
   val = rawLumi*gamma*beta/(1. + alpha_1*Ldot + alpha_2*Ldot*Ldot);
   return val;
}

int runCorrection(){
   
   std::vector<string> instLumis;

   ifstream ifs( "getLumi_out_179828_v5_unCorr.txt" );
   string temp;

   while( getline( ifs, temp ) )
      instLumis.push_back( temp );
   
   ofstream myfile;
   myfile.open ("getLumi_out_179828_v5_Corr.txt");
   
   for(unsigned int i=0; i < instLumis.size();i+=2){
      string lS;
      vector<string>::iterator it;
	  double numb, pixLumi;
	  
	  lS = instLumis.at(i);
	  
      // iterator to vector element:
      it = find (instLumis.begin(), instLumis.end(), lS);
      ++it;

      std::istringstream ( *it ) >> numb;
      
      pixLumi = correctedLumi(numb, 10, 3);
      myfile << lS << "\n" << pixLumi << "\n" << pixLumi*0.02714 << std::endl;
      
   }
   
   myfile.close();  
   
return 0;
}