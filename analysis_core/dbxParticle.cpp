#include "dbxParticle.h"

#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>

#ifndef __DBX_PARTICLE_C__
ClassImp(dbxParticle)

dbxParticle:: dbxParticle() : TObject() {
  p_charge=0; // not initialized
  p_pdgID= 0.;
  p_lvector.SetPtEtaPhiM(0, 0, 0, 0);
  p_istight=0;
  p_ismedium=0;
  p_isloose=0;
  p_attribute.clear();
}
dbxParticle:: dbxParticle (TLorentzVector lv){
  p_istight = 0;
  p_lvector=lv;
  p_istight=0;
  p_ismedium=0;
  p_isloose=0;
  p_attribute.clear();
}

dbxParticle:: ~dbxParticle() {}

dbxParticle:: dbxParticle (TLorentzVector lv, int q){
  p_pdgID = 0.;
  p_charge=q; // initalized
  p_lvector=lv;
  p_istight=0;
  p_ismedium=0;
  p_isloose=0;
  p_attribute.clear();
}
dbxParticle  dbxParticle::operator+ (dbxParticle& p)
{
     dbxParticle result (this->lv()+p.lv(), this->q()+p.q() );
//    dbxParticle result ;
//    result.setTlv(this->lv()+p.lv());
//    result.setCharge(this->q()+p.q());
      return result;
}

int dbxParticle::setPdgID(int pdgID)
{
	p_pdgID = pdgID;
	return 0;
}

double dbxParticle::deltaR(dbxParticle p1,dbxParticle p2){
    
  return p1.lv().DeltaR(p2.lv());
}


double dbxParticle::deltaPhi(dbxParticle p1,dbxParticle p2){
  
  return p1.lv().DeltaPhi(p2.lv());
}

int dbxParticle:: setCharge (int q){
 p_charge=q;
 return 0;
}

int dbxParticle:: setFlavor (double fla){
 p_flavor=fla;
 return 0;
}

int dbxParticle:: setEtCone (double iso){
 p_et_cone=iso;
 return 0;
}

int dbxParticle:: setPtCone (double iso){
 p_pt_cone=iso;
 return 0;
}

int dbxParticle:: setIsTight (int indx){
 p_istight=indx;
 return 0;
}
int dbxParticle:: setIsMedium (bool indx){
 p_ismedium=indx;
 return 0;
}
int dbxParticle:: setIsLoose (bool indx){
 p_isloose=indx;
 return 0;
}

int dbxParticle:: setParticleIndx (int indx){
 p_particleindx=indx;
 return 0;
}

int dbxParticle:: setTlv (TLorentzVector lv){
 p_lvector=lv;
 return 0;
}

int dbxParticle:: scaleLorentzVector ( double scale ){
 p_lvector*=scale;
 return 0;
}

int dbxParticle:: scalePt ( double scale ){
 p_lvector.SetPerp(p_lvector.Perp() * scale);
 return 0;
}

int dbxParticle:: scaleE ( double scale ){
 p_lvector.SetE( p_lvector.E() * scale);
 return 0;
}

int dbxParticle:: setZ0(double q){
  p_z0=q;
 return 0;
}

void dbxParticle:: dump (){
 std::cout << "Px="<<p_lvector.Px()<< "  Py="<<p_lvector.Py()<< "  Pz="<<p_lvector.Pz()<< "  E="<<p_lvector.E()<<std::endl;
}


void dbxParticle:: dump_b (){
 std::cout << "PT="<<p_lvector.Pt()<< "  Eta="<<p_lvector.Eta()<< "  Phi="<<p_lvector.Phi()<< "  M="<<p_lvector.M()<<std::endl;
}

void dbxParticle:: dumpLHCO (std::ofstream& fn){
 using namespace std;
//  #  typ      eta    phi      pt    jmas  ntrk  btag   had/em  sv0  dum2
 fn << fixed <<setprecision(1)<<p_lvector.Eta()<< setw(9)<<setprecision(1)<<p_lvector.Phi()<< setw(9)<<setprecision(1)<<p_lvector.Pt()<< setw(9)<<setprecision(1)<<p_lvector.M();
 fn << fixed<<setw(9)<<setprecision(1)<<p_charge<< setw(9)<<0.0<< setw(9)<<setprecision(1)<< 0.0 << setw(9)<<setprecision(2);
 if (fabs(p_charge)==1) fn << 0.0;
 else fn<<p_flavor;
 fn<<setw(9)<<0.0<<std::endl;
}


#define __DBX_PARTICLE_C__
#endif
