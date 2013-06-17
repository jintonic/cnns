#include "XMASS835kg.h"

#include <TH1D.h>

#include <MAD/NaturalXe.h>
#include <MAD/LiquidXenon.h>
using namespace MAD;

#include <UNIC/Units.h>
using namespace UNIC;

//______________________________________________________________________________
//

XMASS835kg::XMASS835kg(const char *name, const char *title) :
   LXeDetector(name, title), fHEff(0)
{
   fNatXe = new NaturalXe;
   fLXe = new LiquidXenon;
   fLXe->AddElement(fNatXe,1);

   SetTargetMaterial(fLXe);
   SetTargetMass(835*kg);
   SetLightYield(14.7*PE/keV);
   SetThreshold(0*PE);
}

//______________________________________________________________________________
//

XMASS835kg::~XMASS835kg()
{
   if (fNatXe) delete fNatXe;
   if (fLXe) delete fLXe;
   if (fHEff) delete fHEff;
}

//______________________________________________________________________________
//

TH1D* XMASS835kg::HEff()
{
   if (fHEff) return fHEff;

   // define bins
   Int_t nbinse=0;
   Double_t e=0, de, ebins[200];
   while (e<9.9999) {
      if (e<4.9999) de=0.1;
      else de=0.25;
      ebins[nbinse]=e;
      nbinse++;
      e+=de;
   }
   ebins[nbinse]=e;

   // create histogram
   fHEff = new TH1D("hEff","",nbinse,ebins);

   // data from simulation
   Double_t content[70] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.002915452, 0.002915452, 0.007575758, 0.01567398, 0.03643725,
      0.04020101, 0.06019656, 0.1010962, 0.1197982, 0.1404151, 0.194702,
      0.2350792, 0.2728365, 0.3432099, 0.3808323, 0.421519, 0.5163728,
      0.5214564, 0.5541069, 0.6514523, 0.6504425, 0.6865022, 0.734255,
      0.7682584, 0.7925697, 0.7879257, 0.832084, 0.8411215, 0.8576329,
      0.8799314, 0.8912685, 0.9057971, 0.9087719, 0.9212454, 0.9216418,
      0.9265537, 0.9217877, 0.9487705, 0.9663659, 0.9707965, 0.9835069,
      0.9817658, 0.9883495, 0.9886598, 0.9918946, 0.997897, 0.9952886,
      0.9963548, 0.9987245, 0.997426, 1, 1, 1, 0.9984301, 1, 1, 0.9982993, 1}; 

   Double_t errors[70] = {0.002919699, 0.002919699, 0.002919699, 0.002919699,
      0.002919699, 0.002919699, 0.002919699, 0.002919699, 0.002919699,
      0.002919699, 0.002919699, 0.002919699, 0.002919699, 0.002919699,
      0.0038022, 0.004995241, 0.007138964, 0.00724804, 0.008854557, 0.01164417,
      0.01300644, 0.01398288, 0.0175526, 0.01880541, 0.02043036, 0.02385664,
      0.02575139, 0.02754047, 0.03140324, 0.03212001, 0.03350731, 0.03857493,
      0.03979145, 0.04099254, 0.0442272, 0.04368039, 0.04689651, 0.0466983,
      0.0478072, 0.04911375, 0.05227528, 0.05326735, 0.05269708, 0.0559222,
      0.05516543, 0.05693545, 0.0574824, 0.05798005, 0.05743558, 0.06155325,
      0.03948218, 0.04114772, 0.0411509, 0.04321119, 0.04368004, 0.04502128,
      0.04474117, 0.04578668, 0.04836415, 0.04916159, 0.05045931, 0.05063668,
      0.05311191, 0.05471757, 0.05517373, 0.0559672, 0.05647825, 0.05832118,
      0.05824679, 0.06074567};

   for (Int_t ix=1; ix<=nbinse; ix++) {
      fHEff->SetBinContent(ix,content[ix-1]);
      fHEff->SetBinError(ix,errors[ix-1]);
   }

   fHEff->SetXTitle("nuclear recoil energy [keV]");
   fHEff->SetYTitle("acceptance");
   fHEff->SetStats(0);

   return fHEff;
}

//______________________________________________________________________________
//

Double_t XMASS835kg::Acceptance(Double_t Enr)
{
   return HEff()->Interpolate(Enr/keV);
}
