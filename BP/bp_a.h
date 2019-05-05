#ifndef BP_A_H
#define BP_A_H

#include "dbx_a.h"
#include "ReadCard.h"
#include <iostream>
#include "analysis_core.h"
#include "myParticle.h"
#include "Node.h"
#include <list>

class BPdbxA : public dbxA {
  public: 
      BPdbxA(char *aname) : dbxA ( aname)
         {
         sprintf (cname,"%s",aname); // keep the current analysis name in the class variable
//       int r=dbxA::setDir(cname);  // make the relevant root directory
//       if (r)  std::cout <<"Root Directory Set Failure in:"<<cname<<std::endl;
         //grl_cut=false;
         }

      int initGRL();
      int readAnalysisParams();
      int plotVariables(int sel);
      int printEfficiencies();
      int bookAdditionalHistos();
      int makeAnalysis(AnalysisObjects ao ); 
      int saveHistos() {
        int r = dbxA::saveHistos();
        return r;
      }

   private:
      bool grl_cut;
      char cname[CHMAX];
      char algoname[CHMAX];


      // Sezen's handmade histograms
      //TH1D* h_eff = new TH1D("h_eff", "efficiencies", 20, 0, 20);
      TH1D* mWHh1;
      TH1D* mWHh2;
      TH1D* mTopHh1;
      TH1D* mTopHh2;
      TH1D* WHbRh1;
      TH1D* WHbRh2;
      TH1D* xWHbRh1;
      TH1D* xWHbRh2;
      TH1D* Zlm;

      Double_t minwidthw, minwidthz, minmassmm, mindrmm, maxetconem;
      Int_t chargeveto, minnjetsw;  
      int TRGe, TRGm;
 
//DBXcut relevant variables
        std::vector< std::vector<std::string> > myopelist;
        std::vector< std::vector<std::string> > terOpelistT;
        std::vector< std::vector<std::string> > terOpelistF;
        std::vector<   int                    > isTernary;
        std::vector< pair< int , int >        > forbidthese;
  
        map < string, vector<string> > obj_names;

        std::vector<int> ret_i;
        std::vector<int> ret_t;
        vector <int>  found_type_vecs, found_idx_vecs, found_idx_origs;
//-------------------------anna
        list<string> parts; //for def of particles as given by user
        map<string,Node*> NodeVars;//for variable defintion
        map<string,vector<myParticle*> > ListParts;//for particle definition
        map<int,Node*> NodeCuts;//cuts and histos
        map<string,Node*> ObjectCuts;//cuts for user defined objects

};
#endif
