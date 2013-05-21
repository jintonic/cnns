#include "Experiment.h"

#include "Source.h"

#include <Element.h>

#include <TAxis.h>
#include <TGeoMaterial.h>
#include <TMath.h>
using namespace TMath;

#include <CLHEP/Units/PhysicalConstants.h>
using namespace CLHEP;

//______________________________________________________________________________
//

Experiment::Experiment(TGeoMaterial *material, Source *source) : 
   TObject(), fSource(source), fMaterial(material), fMass(0), fThreshold(0)
{
   for (UShort_t i=0; i<7; i++) {
      fXSxNe[i]=0;
      fNevt[i]=0;
   }
}

//______________________________________________________________________________
//

Experiment::~Experiment() 
{
   for (UShort_t i=0; i<7; i++) {
      if (fNevt[i]) delete fNevt[i];
      if (fXSxNe[i]) delete fXSxNe[i];
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

   if (type<1 || type>6) {
      Double_t total = 0;
      for (UShort_t i=1; i<=6; i++) total += fSource->Ne(i,Ev);
      return dXS * total;
   }

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
   if (fMaterial->GetNelements()!=1) {
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
   Double_t atomicMass  = element->A()*gram/mole;
   Double_t nNuclei = fMass/atomicMass*Avogadro;
   Double_t distance = fSource->Distance()/hbarc;
   Double_t area = 4*pi*distance*distance;
   Double_t minEv = (Er + Sqrt(2*element->M()*Er))/2;

   TF1 *f = FXSxNe(type,Er,maxEv);
   return nNuclei/area*f->Integral(minEv/MeV, maxEv/MeV)*keV;
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
      fNevt[type]->SetRange(0,maxNuclearRecoilEnergy/keV);
      fNevt[type]->SetParameter(0,maxNeutrinoEnergy/MeV);
      fNevt[type]->SetParameter(1,type);
      return fNevt[type];
   }

   fNevt[type] = new TF1(Form("fNevt_%s_%s_%f_%d",
            fMaterial->GetName(), fSource->GetName(), fMass, type),
         this, &Experiment::FuncN, 0, maxNuclearRecoilEnergy/keV,2);
   fNevt[type]->SetParameter(0,maxNeutrinoEnergy/MeV);
   fNevt[type]->SetParameter(1,type);
   fNevt[type]->GetXaxis()->SetTitle("nuclear recoil energy [keV]");
   fNevt[type]->GetYaxis()->
      SetTitle(Form("events/(keV #times %.0f kg #times 20 s)", fMass/kg));

   if (type==0) {
      fNevt[type]->SetTitle(Form("%s", fSource->GetTitle()));
      fNevt[type]->SetLineColor(kGray+2);
      fNevt[type]->SetLineWidth(2);
   } else {
      fNevt[type]->SetTitle(
            Form("%s, type of neutrino: %d", fSource->GetTitle(), type));
      fNevt[type]->SetLineColor(type);
   }
   return fNevt[type];
}

//______________________________________________________________________________
//

TF1* Experiment::FXSxNe(UShort_t type,
      Double_t nuclearRecoilEnergy, Double_t maxNeutrinoEnergy)
{
   if (fXSxNe[type]) {
      fXSxNe[type]->SetRange(0,maxNeutrinoEnergy/MeV);
      fXSxNe[type]->SetParameter(0,nuclearRecoilEnergy/keV);
      fXSxNe[type]->SetParameter(1,type);
      return fXSxNe[type];
   }

   fXSxNe[type] = new TF1(Form("fXSxNe_%s_%s_%f_%d",
            fSource->GetName(), fMaterial->GetName(), fMass, type),
         this, &Experiment::XSxNe, 0, maxNeutrinoEnergy/MeV,2);
   fXSxNe[type]->SetParameter(0,nuclearRecoilEnergy/keV);
   fXSxNe[type]->SetParameter(1,type);
   fXSxNe[type]->GetXaxis()->SetTitle("neutrino energy [MeV]");
   fXSxNe[type]->GetYaxis()->SetTitle("1/MeV^{4}");

   if (type==0) {
      fXSxNe[type]->SetTitle(Form("neutrinos from %s, target: %s", 
               fSource->GetName(), fMaterial->GetTitle()));
      fXSxNe[type]->SetLineColor(kGray+2);
      fXSxNe[type]->SetLineWidth(2);
   } else {
      fXSxNe[type]->SetTitle(Form("neutrino %d from %s, target: %s", 
               type, fSource->GetName(), fMaterial->GetTitle()));
      fXSxNe[type]->SetLineColor(type);
   }
   return fXSxNe[type];
}

