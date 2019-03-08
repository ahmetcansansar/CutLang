#include "TParameter.h"
#include <TRandom.h>
#include "TText.h"

#include "bp_a.h"
#include "analysis_core.h"
#include "dbx_a.h"

//#define __VERBOSE3__
//#define _CLV_

#ifdef _CLV_
#define DEBUG(a) std::cout<<a
#else
#define DEBUG(a)
#endif

extern int yyparse(list<string> *parts,map<string,Node*>* NodeVars,map<string,vector<myParticle*> >* ListParts,map<int,Node*>* NodeCuts, map<string,Node*>* ObjectCuts, vector<double>* PtEtaInitializations , vector<double>* btagValues);

extern FILE* yyin;
extern int cutcount;


int BPdbxA::plotVariables(int sel) {
 return 0;  
}

//--------------------------
int BPdbxA:: readAnalysisParams() {
  int retval=0;

  dbxA::ChangeDir(cname);
  TString CardName=cname;
          CardName+="-card.ini";

  ifstream cardfile(CardName);
  if ( ! cardfile.good()) {
    cerr << "The cardfile " << CardName << " file has problems... " << endl;
    return -1;
  }

// ---------------------------DBX style defs, objects and cuts
    int kk=1;

    string tempLine;
    string tempS1, tempS2;
    string subdelimiter = " ";
    string hashdelimiter = "#";
    size_t found;
    bool foundInFile(false);
    TString DefList2file="\n";
    TString CutList2file="\n";
    TString ObjList2file="\n";
    std:vector<TString> effCL;


    bool algorithmnow=false;

    while ( ! cardfile.eof() ) {
       getline( cardfile, tempLine );
       if ( tempLine[0] == '#' ) continue; // skip comment lines
//---------obj
       found = tempLine.find("obj ");
       if (found!=std::string::npos) {
           ObjList2file+=tempLine;
           ObjList2file+="\n";
           continue;
       }

//---------algo
       found = tempLine.find("algo ");
       if (found!=std::string::npos) {
           algorithmnow=true;
           continue;
       }
//---------defs
       found = tempLine.find("def ");
       if (found!=std::string::npos) {
           DefList2file+=tempLine;
           DefList2file+="\n";
           continue;
       }
//---------cmds
       found=tempLine.find("cmd ");
       if (found!=std::string::npos) {
           if (algorithmnow) {
              CutList2file+=tempLine;
              CutList2file+="\n";
              size_t apos=tempLine.find(hashdelimiter);
              tempS1 = tempLine.substr(4, apos-4);
              tempS1.erase(remove_if(tempS1.begin(), tempS1.end(), ::isspace), tempS1.end());
              cout <<tempS1<<"\n";
              effCL.push_back(tempS1);
           } else {
              ObjList2file+=tempLine;
              ObjList2file+="\n";
           }
           continue;
       }
//---------histos
       found=tempLine.find("histo ");
       if (found!=std::string::npos) {
           CutList2file+=tempLine;
           CutList2file+="\n";
           size_t apos=tempLine.find(hashdelimiter);
           tempS1 = tempLine.substr(6, apos-6); // without the comments
           apos=tempS1.find_first_of('"');
           size_t bpos=tempS1.find_last_of('"');
           tempS1 = tempS1.substr(apos+1, bpos-apos-1); // without the comments
           tempS2 = "[Histo] ";
           tempS2 += tempS1;
           cout <<tempS2<<"\n";
           effCL.push_back(tempS2);
           continue;
       }

    } 

//----------put into output file as text
    TText cinfo(0,0,CutList2file.Data());
          cinfo.SetName("CLA2cuts");
          cinfo.Write();

    TText info(0,0,DefList2file.Data());
          info.SetName("CLA2defs");
          info.Write();

    TText oinfo(0,0,ObjList2file.Data());
          oinfo.SetName("CLA2Objs");
          oinfo.Write();




// ****************************************
// ---------------------------DBX style cuts
       eff->GetXaxis()->SetBinLabel(1,"all Events"); // this is hard coded.
       int kFillHistos=0;
    
       std::vector<double> PtEtaInitializations(11);
       PtEtaInitializations={15., 15., 15., 15., 15., 2.5, 2.5, 2.5, 2.5, 30, 1, 0};
       vector<double> btagValues=vector<double>(6);

       yyin=fopen(CardName,"r");
       if (yyin==NULL) { cout << "Cardfile "<<CardName<<" has problems, please check\n";}
       cutcount=0;
       retval=yyparse(&parts,&NodeVars,&ListParts,&NodeCuts, &ObjectCuts, &PtEtaInitializations, &btagValues);
       if (retval){
         cout << "\nSYNTAX error check the input file\n";
         exit (99); 
       }
       cout << "\nWe have "<<NodeCuts.size() << " CutLang Cuts and "<<ObjectCuts.size()  <<" CutLang objects cuts\n";

   minpte  = PtEtaInitializations[0];
   minptm  = PtEtaInitializations[1];
   minptj  = PtEtaInitializations[2];
   minptg  = PtEtaInitializations[3];
   maxetae = PtEtaInitializations[4];
   maxetam = PtEtaInitializations[5];
   maxetaj = PtEtaInitializations[6];
   maxetag = PtEtaInitializations[7];
   maxmet  = PtEtaInitializations[8];
   TRGe    = PtEtaInitializations[9];
   TRGm    = PtEtaInitializations[10];

    eff->GetXaxis()->SetBinLabel(1,"all Events"); // this is hard coded.

    DEBUG("CL CUTS: \n");
    std::map<int, Node*>::iterator iter = NodeCuts.begin();
    while(iter != NodeCuts.end())
    {
            DEBUG(" CUT "<<iter->first<<" ");
            DEBUG("--->"<<iter->second->getStr()<<"\n");

//           TString newLabels=iter->second->getStr();
           TString newLabels=effCL[ iter->first -1];
/*
            TString newLabels="CUT";
                    newLabels+=iter->first;
 */
           eff->GetXaxis()->SetBinLabel(iter->first+1,newLabels); // labels

            DEBUG(std::endl);
            iter++; 
    }

#ifdef _CLV__
     cout<<"\n Particle Lists: \n";

     for (map<string,vector<myParticle*> >::iterator it1 = ListParts.begin(); it1 != ListParts.end(); it1++)
         {
         cout << (*it1)->first << ": ";
         for (vector<myParticle*>::iterator lit = it1->second.begin(); lit  != it1->second.end(); lit++)
         cout << (*lit)->type << "_" << (*lit)->index << " ";
         cout << "\n";
         }

    cout<<"\n Particles defintions as given by user: \n";

    std::list<std::string>::iterator it = parts.begin();
    while(it != parts.end())
    {
            std::cout<<(*it)<<std::endl;
            it++;
    }

    cout<<"\n Variables results: \n";
    map<string,Node* >::iterator itv = NodeVars.begin();
    while(itv != NodeVars.end())
    {
            std::cout<<"**************************** "<<itv->first<<endl;
            itv->second->display();
            std::cout<<std::endl;
            itv++;
    }

#endif


// PUT ANALYSIS PARAMETERS INTO .ROOT //////////////
	
    TParameter<double> *minpte_tmp=new TParameter<double> ("minpte", minpte);
    TParameter<double> *maxetae_tmp=new TParameter<double> ("maxetae", maxetae);
    TParameter<double> *minptm_tmp=new TParameter<double> ("minptm", minptm);
    TParameter<double> *maxetam_tmp=new TParameter<double> ("maxetam", maxetam);
    TParameter<double> *minptj_tmp=new TParameter<double> ("minptj", minptj);
    TParameter<double> *maxetaj_tmp=new TParameter<double> ("maxetaj", maxetaj);
    TParameter<double> *TRGe_tmp=new TParameter<double> ("TRGe", TRGe);
    TParameter<double> *TRGm_tmp=new TParameter<double> ("TRGm", TRGm);

    minpte_tmp->Write("minpte");
    maxetae_tmp->Write("maxetae");
    minptm_tmp->Write("minptm");
    maxetam_tmp->Write("maxetam");
    minptj_tmp->Write("minptj");
    maxetaj_tmp->Write("maxetaj");
    TRGe_tmp->Write("TRGe");
    TRGm_tmp->Write("TRGm");

  return retval;
}

int BPdbxA:: printEfficiencies() {
  int retval=0;
  PrintEfficiencies(eff);
  return retval;
}

int BPdbxA:: initGRL() {
  int retval=0;
  grl_cut=true;
  return retval;
}

int BPdbxA:: bookAdditionalHistos() {
        int retval=0;
        dbxA::ChangeDir(cname);

#ifdef __VERBOSE3__
	// Sezen's handmade histograms
	mWHh1 = new TH1D("mWHh1", "Hadronic W best combi (GeV)", 50, 50, 150);
	mWHh2 = new TH1D("mHWh2", "Hadronic W best combi (GeV)", 50, 50, 150);
	mTopHh1 = new TH1D("mTopHh1", "Hadronic top combi (GeV)", 70, 0, 700);
	mTopHh2 = new TH1D("mTopHh2", "Hadronic top combi (GeV)", 70, 0, 700);
	WHbRh1 = new TH1D("WHbRh1", "Angular distance between W1 and bjet", 70, 0, 7);
	WHbRh2 = new TH1D("WHbRh2", "Angular distance between W2 and bjet", 70, 0, 7);
	xWHbRh1 = new TH1D("xWHbRh1", "Hadronic top combi (GeV) after angular cut", 70, 0, 700);
	xWHbRh2 = new TH1D("xWHbRh2", "Hadronic top combi (GeV) after angular cut", 70, 0, 700);
#endif

// ---------------------------DBX style defs from the main file

  return retval;
}

/////////////////////////
int BPdbxA::makeAnalysis(vector<dbxMuon> muons, vector<dbxElectron> electrons, vector <dbxPhoton> photons,
                         vector<dbxJet> jets, TVector2 met, evt_data anevt) {
  int retval=0;

  vector<dbxElectron>  goodElectrons;
  vector<dbxMuon>      goodMuons;
  vector<dbxJet>       goodJets;
  vector<dbxPhoton>    goodPhotons;

  DEBUG("-------------------------------------------------------------------- "<<cname<<"\n");
//----------------------selection of good gams-----------------
        for (UInt_t i=0; i<photons.size(); i++) {
               TLorentzVector gam4p = photons.at(i).lv();
               if (    (gam4p.Pt()  > minptg)
                    && (fabs(gam4p.Eta()) < maxetag)
                  )
                  goodPhotons.push_back(photons.at(i));
        }
//----------------------selection of good electrons-----------------
        for (UInt_t i=0; i<electrons.size(); i++) {
               if ( (electrons.at(i).lv().Pt()  > minpte)    // the electrons should have a minimum PT
                  &&(electrons.at(i).lv().Eta() < maxetae )  // and maximum eta.
                  )
                  goodElectrons.push_back( electrons.at(i) );
        }

//----------------------selection of good muons-----------------
        for (UInt_t i=0; i<muons.size(); i++) {
               TLorentzVector mu4p = muons.at(i).lv();
               if (    (mu4p.Pt()  > minptm)
                    && (fabs(mu4p.Eta()) < maxetam)
                  )
                  goodMuons.push_back(muons.at(i));
        }

//------------selection of good jets----------------------------------
        for (UInt_t i=0; i<jets.size(); i++) {
               TLorentzVector jet4p = jets.at(i).lv();
               if (   (fabs(jet4p.Pt())  > minptj ) // this corresponds to min PT cut
                    && (jet4p.E() >= 0)
                    && (fabs(jet4p.Eta())<= maxetaj) 
                   )
                   goodJets.push_back(jets.at(i) );
        }

///////
        double theLeptonWeight = 1;
        double theFourJetWeight = 1;
        unsigned int njets;
        double evt_weight = 1;

        if(TRGe==2 || TRGm== 2) evt_weight = anevt.weight_mc*anevt.weight_pileup*anevt.weight_jvt;//                                                                                                                                                                 
// --------- INITIAL  # events  ====> C0
        eff->Fill(1, 1);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    AnalysisObjects a0={goodMuons, goodElectrons, goodPhotons, goodJets, met, anevt};
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DEBUG("------------------------------------------------- Event ID:"<<anevt.event_no<<" \n");

//    std::cout<<"\n--------------Starting New Event: "<<anevt.event_no<<"  ";

// *************************************
/// CutLang execution starts-------here*
// *************************************

    unsigned int ternaryCount=0;
    std::map<int, Node*>::iterator iter = NodeCuts.begin();
    DEBUG("Start resetting cuts:"<< NodeCuts.size() <<"\n");
//----------------------reset 

    while(iter != NodeCuts.end())
    {   
        iter->second->Reset();
        iter++;
    }

    DEBUG("RESet ALL cuts\n");
    iter = NodeCuts.begin();

//----------------------execute
    while(iter != NodeCuts.end())
    {   
        a0={goodMuons, goodElectrons, goodPhotons, goodJets, met, anevt}; // we start from good ones.

        DEBUG("Selecting: "<<iter->first<<" |");
        double d=iter->second->evaluate(&a0); // execute the selection cut
        DEBUG(" Result : " << d << std::endl);
        if (d==0) return iter->first;         // quits the event.
        eff->Fill(iter->first+1, evt_weight); // filling starts from 1 which is already filled.
        iter++; //moves on to the next cut
    } // loop over all cutlang cuts
    DEBUG("   EOE\n     ");
return 1;
}
