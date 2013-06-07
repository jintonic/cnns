#include "Experiment.h"

#include "Source.h"

#include <MAD/Element.h>
#include <MAD/Material.h>
using namespace MAD;

#include <TF1.h>
#include <TAxis.h>
#include <TMath.h>
using namespace TMath;

#include <UNIC/Units.h>
#include <UNIC/Constants.h>
using namespace UNIC;

//______________________________________________________________________________
//

Experiment::Experiment(Material *material, Source *source) : 
   TObject(), fSource(source), fMaterial(material), fMass(0), fThreshold(0)
{
   for (UShort_t i=0; i<7; i++) {
      fXSxNe[i]=0;
      fNevt[i]=0;
   }
}

//______________________________________________________________________________
//

Double_t Experiment::XSxNe(Double_t *x, Double_t *parameter)
{
   Double_t Ev = x[0]*MeV; // neutrino energy
   Double_t Er = parameter[0]*keV; // nuclear recoil energy
   UShort_t type = static_cast<UShort_t>(parameter[1]); // type of neutrino

   Element *element = (Element*) fMaterial->GetElement();
   Double_t dXS = element->CNNSdXS(Er, Ev);

   if (type<1 || type>6)
      return dXS*(fSource->Ne(1,Ev) + fSource->Ne(2,Ev) + 4*fSource->Ne(3,Ev));

   return dXS * fSource->Ne(type,Ev);
}

//______________________________________________________________________________
//

Double_t Experiment::FuncN(Double_t *x, Double_t *parameter)
{
   if (!fMaterial) {
      Warning("Nevt", "Please set targe material!");
      return 0;
   }
   if (fMaterial->Nelements()!=1) {
      Warning("Nevt", "Can only handle material with one element!");
      return 0;
   }
   if (!fSource) {
      Warning("Nevt", "Please set source of neutrino!");
      return 0;
   }

   Double_t Er = x[0]*keV; // nuclear recoil energy
   Double_t maxEv = parameter[0]*MeV; // max neutrino energy
   UShort_t type = static_cast<UShort_t>(parameter[1]); // type of neutrino

   Element *element = (Element*) fMaterial->GetElement();
   Double_t atomicMass  = element->A();
   Double_t nNuclei = fMass/atomicMass*Avogadro;
   Double_t distance = fSource->Distance()/hbarc;
   Double_t area = 4*pi*distance*distance;
   Double_t minEv = (Er + Sqrt(2*element->M()*Er))/2;

   //Printf("Er: %f keV, minEv: %f MeV, maxEv: %f MeV ==========", x[0], minEv, maxEv);
   TF1 *f = FXSxNe(type,Er,minEv,maxEv);
   return nNuclei/area*1e50*f->Integral(minEv/MeV, maxEv/MeV)*keV;
}

//______________________________________________________________________________
//

Double_t Experiment::Nevt(UShort_t type,
      Double_t nuclearRecoilEnergy, Double_t maxNeutrinoEnergy)
{
   Double_t x[1] = {nuclearRecoilEnergy/keV};
   Double_t p[2] = {maxNeutrinoEnergy/MeV, type};
   return FuncN(x,p);
}

//______________________________________________________________________________
//

TF1* Experiment::FNevt(UShort_t type,
      Double_t maxNuclearRecoilEnergy, Double_t maxNeutrinoEnergy)
{
   if (fNevt[type]) {
      Double_t min, max;
      fNevt[type]->GetRange(min,max);
      if (max!=maxNuclearRecoilEnergy/keV) {
         Info("FNevt","Reset range of recoil energy.");
         fNevt[type]->SetRange(0,maxNuclearRecoilEnergy/keV);
      }
      if (fNevt[type]->GetParameter(0)!=maxNeutrinoEnergy/MeV) {
         Info("FNevt","Reset maximal neutrino energy.");
         fNevt[type]->SetParameter(0,maxNeutrinoEnergy/MeV);
      }
      return fNevt[type];
   }

   fNevt[type] = new TF1(Form("fNevt_%s_%s_%f_%d",
            fMaterial->GetName(), fSource->GetName(), fMass, type),
         this, &Experiment::FuncN, 0, maxNuclearRecoilEnergy/keV,2);
   fNevt[type]->SetParameter(0,maxNeutrinoEnergy/MeV);
   fNevt[type]->SetParameter(1,type);

   if (type==0) {
      fNevt[type]->SetTitle(Form(
               "%s;nuclear recoil energy [keV];total events/(keV #times %.0f kg)",
               fSource->GetTitle(), fMass/kg));
      fNevt[type]->SetLineColor(kGray+2);
      fNevt[type]->SetLineWidth(2);
   } else {
      fNevt[type]->SetTitle(Form(
               "%s, type of neutrino: %d;nuclear recoil energy [keV];total events/(keV #times %.0f kg)",
               fSource->GetTitle(), type, fMass/kg));
      fNevt[type]->SetLineColor(type);
   }
   return fNevt[type];
}

//______________________________________________________________________________
//

TF1* Experiment::FXSxNe(UShort_t type, Double_t nuclearRecoilEnergy,
      Double_t minNeutrinoEnergy, Double_t maxNeutrinoEnergy)
{
   if (fXSxNe[type]) {
      fXSxNe[type]->SetRange(minNeutrinoEnergy/MeV,maxNeutrinoEnergy/MeV);
      fXSxNe[type]->SetParameter(0,nuclearRecoilEnergy/keV);
      return fXSxNe[type];
   }

   fXSxNe[type] = new TF1(Form("fXSxNe%s%s%f%d",
            fSource->GetName(), fMaterial->GetName(), fMass, type), this, 
         &Experiment::XSxNe, minNeutrinoEnergy/MeV, maxNeutrinoEnergy/MeV,2);
   fXSxNe[type]->SetParameter(0,nuclearRecoilEnergy/keV);
   fXSxNe[type]->SetParameter(1,type);

   if (type==0) {
      fXSxNe[type]->SetTitle(Form(
               "%s, target: %s, recoil energy: %.1f keV;neutrino energy [MeV];1/MeV^{4}", 
               fSource->GetName(), fMaterial->GetTitle(), 
               nuclearRecoilEnergy/keV));
      fXSxNe[type]->SetLineColor(kGray+2);
      fXSxNe[type]->SetLineWidth(2);
   } else {
      fXSxNe[type]->SetTitle(Form(
               "neutrino %d from %s, target: %s;neutrino energy [MeV];1/MeV^{4}", 
               type, fSource->GetName(), fMaterial->GetTitle()));
      fXSxNe[type]->SetLineColor(type);
   }
   return fXSxNe[type];
}

//______________________________________________________________________________
//

TH1D* Experiment::HNevt(UShort_t type,
      Double_t maxNuclearRecoilEnergy, Double_t maxNeutrinoEnergy)
{
   TH1D *h = (TH1D*) 
      FNevt(type,maxNuclearRecoilEnergy,maxNeutrinoEnergy)->GetHistogram();
   return h;
}

//______________________________________________________________________________
//

TH1D* Experiment::HXSxNe(UShort_t type, Double_t nuclearRecoilEnergy,
      Double_t minNeutrinoEnergy, Double_t maxNeutrinoEnergy)
{
   TH1D *h = (TH1D*) FXSxNe(type, nuclearRecoilEnergy, 
         minNeutrinoEnergy, maxNeutrinoEnergy)->GetHistogram();
   return h;
}

//______________________________________________________________________________
//

void Experiment::Clear(Option_t *option)
{
   for (UShort_t i=0; i<7; i++) {
      if (fNevt[i]) {
         delete fNevt[i];
         fNevt[i]=NULL;
      }
      if (fXSxNe[i]) {
         delete fXSxNe[i];
         fXSxNe[i]=NULL;
      }
   }
}

