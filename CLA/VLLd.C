#define VLLd_cxx
#include "VLLd.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <signal.h>

#include "dbx_electron.h"
#include "dbx_muon.h"
#include "dbx_jet.h"
#include "dbx_tau.h"
#include "dbx_a.h"
#include "DBXNtuple.h"
#include "analysis_core.h"
#include "AnalysisController.h"
#include <iostream>
#include <sstream>      // std::istringstream
#include <string>
#include "TTreeReader.h"
using namespace std;

extern void _fsig_handler (int) ;
extern bool fctrlc;
extern map<string, TTreeReader*> ttr_map;


//#define _CLV_
#ifdef _CLV_
#define DEBUG(a) std::cout<<a
#else
#define DEBUG(a)
#endif


void VLLd::GetPhysicsObjects( Long64_t j, AnalysisObjects *a0 )
{
       fChain->GetEntry(j);

       vector<dbxMuon>     muons;
       vector<dbxElectron> electrons;
       vector<dbxPhoton>   photons;
       vector<dbxJet>      jets;
       vector<dbxTau>      taus;
       vector<dbxJet>    ljets;
       vector<dbxTruth>   truth;
       vector<dbxTrack>    track;
       vector<dbxParticle> combos;
       vector<dbxParticle> constis;

       map<string, vector<dbxMuon>     > muos_map;
       map<string, vector<dbxElectron> > eles_map;
       map<string, vector<dbxTau>      > taus_map;
       map<string, vector<dbxPhoton>   > gams_map;
       map<string, vector<dbxJet>      > jets_map;
       map<string, vector<dbxJet>     >ljets_map;
       map<string, vector<dbxTruth>    >truth_map;
       map<string, vector<dbxTrack>    >track_map;
       map<string, vector<dbxParticle> >combo_map;
       map<string, vector<dbxParticle> >constits_map;
       map<string, TVector2            >  met_map;

       evt_data anevt;

//temporary variables
       TLorentzVector  alv, alv0, alv1, alv2, alv3, alv4, alv5;
       TVector2 met;
       dbxJet      *adbxj;
       dbxElectron *adbxe;
       dbxMuon     *adbxm;

DEBUG("Begin Filling\n");

//JETS
   unsigned int jet_n=jets_pt->size();
        for (unsigned int i=0; i<jet_n; i++) {
                alv.SetPtEtaPhiE( jets_pt->at(i)*0.001, jets_eta->at(i), jets_phi->at(i), jets_e->at(i)*0.001 ); // all in GeV
                adbxj= new dbxJet(alv);
                adbxj->setCharge(-99);
                adbxj->setParticleIndx(i);
//                adbxj->setJVtxf(jets_DFCommonJets_fJvt->at(i));
                adbxj->setFlavor(jets_btagFlag_DL1r_FixedCutBEff_85->at(i) );
                adbxj->set_isbtagged_77( jets_btagFlag_DL1r_FixedCutBEff_85->at(i) ); // 5 is btag
                jets.push_back(*adbxj);
                delete adbxj;
        }
DEBUG("Jets ok\n");

//LEPTONS
      for (int ii=1; ii<2; ii++){
        int nlep=0;
        //if( (int)onelep_type >0 ) nlep=1;
        if( (int)dilep_type >0  ) nlep=2;
        if( (int)trilep_type >0 ) nlep=3;
        if( (int)quadlep_type>0 ) nlep=4;
        //if( (int)fivelep_type>0 ) nlep=5;

        if (nlep==0) break;
        alv0.SetPtEtaPhiE( lep_Pt_0*0.001, lep_Eta_0, lep_Phi_0, lep_E_0*0.001 ); // all in GeV       
        if (abs(lep_ID_0)==13) { // muons
            adbxm= new dbxMuon(alv0);
            adbxm->setCharge(lep_ID_0 / 13 );
            adbxm->setPdgID( lep_ID_0 );
            adbxm->setParticleIndx(0);
            adbxm->setZ0(lep_Z0SinTheta_0 );
            adbxm->setIsLoose(lep_isLooseLH_0);
            adbxm->setIsMedium(lep_isMedium_0);
            muons.push_back(*adbxm);
            delete adbxm;
        }
        if (abs(lep_ID_0)==11) { // electrons
            adbxe= new dbxElectron(alv0);
            adbxe->setCharge(lep_ID_0 / 11 );
            adbxe->setPdgID( lep_ID_0 );
            adbxe->setParticleIndx(0);
            adbxe->setZ0(lep_Z0SinTheta_0 );
            adbxe->setIsLoose(lep_isLooseLH_0);
            adbxe->setIsTight(lep_isTightLH_0);
            electrons.push_back(*adbxe);
            delete adbxe;
        }
        if (nlep==1) break;
        alv1.SetPtEtaPhiE( lep_Pt_1*0.001, lep_Eta_1, lep_Phi_1, lep_E_1*0.001 ); // all in GeV       
//1
        if (abs(lep_ID_1)==13) { // muons
            adbxm= new dbxMuon(alv1);
            adbxm->setCharge(lep_ID_1 / 13 );
            adbxm->setPdgID( lep_ID_1 );
            adbxm->setParticleIndx(1);
            adbxm->setZ0(lep_Z0SinTheta_1 );
            adbxm->setIsLoose(lep_isLooseLH_1);
            adbxm->setIsMedium(lep_isMedium_1);
            muons.push_back(*adbxm);
            delete adbxm;
        }
        if (abs(lep_ID_1)==11) { // electrons
            adbxe= new dbxElectron(alv1);
            adbxe->setCharge(lep_ID_1 / 11 );
            adbxe->setPdgID( lep_ID_1 );
            adbxe->setParticleIndx(1);
            adbxe->setZ0(lep_Z0SinTheta_1 );
            adbxe->setIsLoose(lep_isLooseLH_1);
            adbxe->setIsTight(lep_isTightLH_1);
            electrons.push_back(*adbxe);
            delete adbxe;
        }
        if (nlep==2) break;
        alv2.SetPtEtaPhiE( lep_Pt_2*0.001, lep_Eta_2, lep_Phi_2, lep_E_2*0.001 ); // all in GeV       
//2
        if (abs(lep_ID_2)==13) { // muons
            adbxm= new dbxMuon(alv2);
            adbxm->setCharge(lep_ID_2 / 13 );
            adbxm->setPdgID( lep_ID_2 );
            adbxm->setParticleIndx(2);
            adbxm->setZ0(lep_Z0SinTheta_2 );
            adbxm->setIsLoose(lep_isLooseLH_2);
            adbxm->setIsMedium(lep_isMedium_2);
            muons.push_back(*adbxm);
            delete adbxm;
        }
        if (abs(lep_ID_2)==11) { // electrons
            adbxe= new dbxElectron(alv2);
            adbxe->setCharge(lep_ID_2 / 11 );
            adbxe->setPdgID( lep_ID_2 );
            adbxe->setParticleIndx(2);
            adbxe->setZ0(lep_Z0SinTheta_2 );
            adbxe->setIsLoose(lep_isLooseLH_2);
            adbxe->setIsTight(lep_isTightLH_2);
            electrons.push_back(*adbxe);
            delete adbxe;
        }
        if (nlep==3) break;
        alv3.SetPtEtaPhiE( lep_Pt_3*0.001, lep_Eta_3, lep_Phi_3, lep_E_3*0.001 ); // all in GeV       
//3
        if (abs(lep_ID_3)==13) { // muons
            adbxm= new dbxMuon(alv3);
            adbxm->setCharge(lep_ID_3 / 13 );
            adbxm->setPdgID( lep_ID_3 );
            adbxm->setParticleIndx(3);
            adbxm->setZ0(lep_Z0SinTheta_3 );
            adbxm->setIsLoose(lep_isLooseLH_3);
            adbxm->setIsMedium(lep_isMedium_3);
            muons.push_back(*adbxm);
            delete adbxm;
        }
        if (abs(lep_ID_3)==11) { // electrons
            adbxe= new dbxElectron(alv3);
            adbxe->setCharge(lep_ID_3 / 11 );
            adbxe->setPdgID( lep_ID_3 );
            adbxe->setParticleIndx(3);
            adbxe->setZ0(lep_Z0SinTheta_3 );
            adbxe->setIsLoose(lep_isLooseLH_3);
            adbxe->setIsTight(lep_isTightLH_3);
            electrons.push_back(*adbxe);
            delete adbxe;
        }
       }
//MET
        met.SetMagPhi( met_met*0.001,  met_phi);


//------------ auxiliary information -------
        anevt.run_no=runNumber;
	anevt.ChannelNo=mcChannelNumber;
	anevt.RunYear=RunYear;
        anevt.correction_weight=1.0;
        anevt.luminosity_weight=1.0;
        anevt.weight_xsec=1.0;
        anevt.user_evt_weight=1.0;
        anevt.lumiblk_no=1;
        anevt.top_hfor_type=0;
        anevt.event_no=eventNumber;
        anevt.TRG_e= 1;
        anevt.TRG_m= 1;
        anevt.TRG_j= 0;
        anevt.vxp_maxtrk_no= 9;
        anevt.badjet=0;
        anevt.weight_mc=weight_mc;
	anevt.m_HF_Classification=m_HF_Classification;
        anevt.weight_pileup=weight_pileup;
	anevt.weight_jvt=jvtSF_customOR;
        anevt.z_vtx_weight = 1.0;
        anevt.weight_bTagSF_77 = bTagSF_weight_DL1r_Continuous;
        anevt.weight_leptonSF = custTrigSF_TightElMediumMuID_FCLooseIso_SLTorDLT;
        anevt.vxpType=0;
        anevt.lar_Error=0;
        anevt.core_Flags=0;
        //anevt.maxEvents=nentries;



DEBUG("Filling finished"<<std::endl);
        muos_map.insert( pair <string,vector<dbxMuon>     > ("MUO",         muons) );
        eles_map.insert( pair <string,vector<dbxElectron> > ("ELE",     electrons) );
        taus_map.insert( pair <string,vector<dbxTau>      > ("TAU",          taus) );
        gams_map.insert( pair <string,vector<dbxPhoton>   > ("PHO",       photons) );
        jets_map.insert( pair <string,vector<dbxJet>      > ("JET",          jets) );
       ljets_map.insert( pair <string,vector<dbxJet>      > ("FJET",        ljets) );
       truth_map.insert( pair <string,vector<dbxTruth>    > ("Truth",       truth) );
       track_map.insert( pair <string,vector<dbxTrack>    > ("Track",       track) );
       combo_map.insert( pair <string,vector<dbxParticle> > ("Combo",      combos) );
         met_map.insert( pair <string,TVector2>             ("MET",           met) );
    if (constits_map.size() < 1) // we only add this if it was previously empty...
    constits_map.insert( pair <string,vector<dbxParticle> > ("Constits",  constis) );

        *a0={muos_map, eles_map, taus_map, gams_map, jets_map, ljets_map, truth_map,track_map, combo_map, constits_map, met_map, anevt};
}
 //--------------------------------------------------------LOOP
void VLLd::Loop( analy_struct aselect, char *extname)
{
// Signal HANDLER
  signal (SIGINT, _fsig_handler); // signal handler has issues with CINT
   TFile *afile= ((TChain *)fChain)->GetFile();

   if (fChain == 0) {
          cout <<"Error opening the data file"<<endl; return;
   }
   int verboseFreq(aselect.verbfreq);
   bool  doSystematics(aselect.dosystematics);
   map < string, syst_struct > systematics; // contains all
   map < string, string > syst_names; // contains all
   map < string, VLLd*> syst_objects;

   if (doSystematics) {
       string tempLine;
       cout << "Reading available systematics from ini file...\n";
       TString CardName="BP_1-card.ini";
       ifstream cardfile(CardName);
       if ( ! cardfile.good()) {
         cerr << "The cardfile " << CardName << " file has problems... " << endl;
       }
       int systindex=0;
       while ( ! cardfile.eof() ) {
          getline( cardfile, tempLine );
          if ( tempLine[0] == '#' ) continue; // skip comment lines starting with #
          if (tempLine.find_first_of("#") != std::string::npos ){ // skip anything after #
            tempLine.erase(tempLine.find_first_of("#"));
          }
          if (tempLine.size() < 3) continue; // skip the junk

          std::istringstream iss(tempLine);
          std::vector<std::string> resultstr((std::istream_iterator<std::string>(iss)), std::istream_iterator<std::string>());
          if (resultstr.size() < 1) continue;
          string firstword=resultstr[0];
          for(auto& c : firstword) { c = tolower(c); } // convert to lowercase

//---------do we have a systematic or NOT ?
         if (firstword == "systematic" ) {
              string ison=resultstr[1];
              for(auto& c : ison) { c = tolower(c); } // convert to lowercase
              if (ison == "on") {
                resultstr[2].erase(remove( resultstr[2].begin(), resultstr[2].end(), '\"' ),resultstr[2].end());
                resultstr[3].erase(remove( resultstr[3].begin(), resultstr[3].end(), '\"' ),resultstr[3].end());
                DEBUG("--> Syst  "<< resultstr[2] << " &&  " << resultstr[3]<<" "<<resultstr[4] <<"\n");
                if ( resultstr[4] != "ttree" ) { // maybe sshould contain weight
		 for (int ri=2; ri<4; ri++){
                  size_t findex = resultstr[ri].find_first_of("[");
                  syst_struct asyst;
                  if (findex == std::string::npos ){ // without any []
                      asyst.index=systindex;
                      asyst.vartype=resultstr[4];
                      asyst.varname=resultstr[ri];
                      asyst.varid=0;
                      asyst.systname=resultstr[ri];// I am double counting
                      systematics[ resultstr[ri] ] = asyst;  // add a simple variable
                      syst_names[ resultstr[ri] ] = resultstr[4];  // add a simple variable
                  } else { // here we have [] lets find the numbers
                      size_t lindex = resultstr[ri].find_first_of("]");
                      string numsection=resultstr[ri].substr(findex+1, lindex-findex-1);
                      string delimiters=",";
                      resultstr[ri].erase(findex);
                      stringstream jss(numsection);
                      string intermediate;
                      while(getline(jss, intermediate, ',')) {
                            cout <<" select subset @" << intermediate << " of "<< resultstr[ri] <<"\n";
                            string asysname=resultstr[ri]+"["+intermediate+"]";
                            //cout <<"==========>"<< asysname <<"\n";
                            asyst.vartype=resultstr[4];
                            asyst.varname=resultstr[ri];
                            asyst.systname=asysname;
                            asyst.varid=stoi(intermediate);
                            asyst.index=systindex;
                            systematics[asysname] = asyst ; // with []
                            syst_names[asysname] = resultstr[4] ; // with []
                      }
                  }
                  systindex++;
                  freaders.push_back(TTreeReaderArray<Float_t>( *(ttr_map["nominal"]), resultstr[ri].c_str() ) ); //push only the generic name
             //     cout <<"XXXXXXXXXXXXXXXXX:"<<ri<<" "<< resultstr[ri]<< " @:"<< ttr_map["nominal"] <<" finished\n";
                 }//2 3 counting
               } else { // tree
                    string xxx="XXX";
                    syst_names[resultstr[2]] = resultstr[4] ; // with []
                    syst_names[resultstr[3]] = resultstr[4] ; // with []
                    
             //       cout << "B r2:"<< resultstr[2]<<" r3:"<<resultstr[3]<<" r4:"<<resultstr[4]<<"\n";                   

                    syst_struct asyst;
                    asyst.vartype=resultstr[4];
                    asyst.varname=resultstr[2];
                    asyst.systname=resultstr[2];
                    asyst.varid=-1;
                    asyst.index=systindex;
                    syst_names[resultstr[2]] = resultstr[4] ; // with []
                    systematics[resultstr[2]] = asyst ; // with []
                    syst_objects[resultstr[2]] = new VLLd((char*)xxx.c_str() ,(TChain *)afile->Get(resultstr[2].c_str()) );
                    systindex++;
                    syst_struct bsyst;
                    bsyst.vartype=resultstr[4];
                    bsyst.varname=resultstr[3];
                    bsyst.systname=resultstr[3];
                    bsyst.varid=-1;
                    bsyst.index=systindex;
                    syst_names[resultstr[3]] = resultstr[4] ; // with []
                    syst_objects[resultstr[3]]=new VLLd((char*)xxx.c_str(),(TChain *)afile->Get(resultstr[3].c_str()) );
                    systematics[resultstr[3]] = bsyst ; // with []
                    systindex++;
                    cout << "A r2:"<< resultstr[2]<<" r3:"<<resultstr[3]<<" r4:"<<resultstr[4]<<"\n";                   

                    // TFile * _afile = TFile::Open(fileList[0].c_str());
                    // ttreader =  new TTreeReader(leafname.c_str(), _afile);
               }
              }// end of systematics on
             continue;
          }
          
       }// end of while reading the file
   }// end of do systematics

   //--------------start stop event ids etc
   map < string,   AnalysisObjects > analysis_objs_map;

   AnalysisController aCtrl(&aselect, syst_names);
   aCtrl.Initialize(extname);
   cout << "End of analysis initialization"<<endl;

   Long64_t nentries = fChain->GetEntriesFast();
   if (aselect.maxEvents>0 ) nentries=aselect.maxEvents;
   cout << "number of entries " << nentries << endl;
   Long64_t startevent = 0;
   if (aselect.startpt>0 ) startevent=aselect.startpt;
   cout << "starting entry " << startevent << endl;
   Long64_t lastevent = startevent + nentries;
   if (lastevent > fChain->GetEntriesFast() ) { lastevent=fChain->GetEntriesFast();
       cout << "Interval exceeds tree. Analysis is done on max available events starting from event : " << startevent << endl;
   }

// ******************************************************************
   for (Long64_t j=startevent; j<lastevent; ++j) { // event loop here

     //  if ( fctrlc ) { cout << "Processed " << j << " events\n"; break; }
       if (0 > LoadTree (j)) break;
       if ( j%verboseFreq == 0 ) cout << "Processing event " << j << endl;

       AnalysisObjects a0;
       GetPhysicsObjects(j, &a0);
       ttr_map["nominal"]->SetEntry(j);
       evt_data oldevt=a0.evt;
       double wvalue=0;
       


       for (map<string,syst_struct>::iterator it = systematics.begin(); it != systematics.end(); it++) {

 //        cout << "it nam:"<<it->first<<"\t"<<it->second.vartype<<"\t" ;
         int jsyst=it->second.index;
         int jid=it->second.varid;
 //        cout <<" jsyst:"<<jsyst <<" jid:"<< jid<<"\n";
         if (jid >=0){
           wvalue = freaders.at(jsyst)[jid]; //may not be zero?
 //          cout <<"Wv:"<<wvalue<<"\n";
           if (jsyst>0) a0.evt=oldevt;
         }

         if (it->second.vartype ==  "weight_jvt" ) {
                a0.evt.weight_jvt=wvalue; // cout <<  "jvt\n";
         	analysis_objs_map[it->first] = a0;
         } else if (it->second.vartype == "weight_pileup" ) {
                a0.evt.weight_pileup=wvalue; // cout <<  "pileup\n";
         	analysis_objs_map[it->first] = a0;
         } else if (it->second.vartype == "weight_leptonSF" ) {
                a0.evt.weight_leptonSF=wvalue; //  cout <<  "SF\n";
         	analysis_objs_map[it->first] = a0;
         } else if (it->second.vartype == "weight_BTagSF" ) {
                a0.evt.weight_bTagSF_77=wvalue; //  cout <<  "w BTAG:" << wvalue <<"\n";
         	analysis_objs_map[it->first] = a0;
         } else if (it->second.vartype == "ttree" ) {
		AnalysisObjects ad;
   //             cout<<"NAme D:"<<it->first<<"\n";
             	syst_objects[it->first]->GetPhysicsObjects(j, &ad);
        	analysis_objs_map[it->first] = ad;
  	      it++;
		AnalysisObjects au;
     //           cout<<"NAme U:"<<it->first<<"\n";
                syst_objects[it->first]->GetPhysicsObjects(j, &au);
        	analysis_objs_map[it->first] = au;

             for ( map<string, TTreeReader*>::iterator itr=ttr_map.begin(); itr!= ttr_map.end(); itr++){
                     (itr->second)->SetEntry(j);
             }


         } else {
	        cout << "problem with: "<<it->second.vartype<< ", no such systematics type.\n";
         }
       } // end of loop over syst subsyst name

       a0.evt=oldevt;
       aCtrl.RunTasks(a0, analysis_objs_map);

   }// end of event loop
   aCtrl.Finalize();

}
