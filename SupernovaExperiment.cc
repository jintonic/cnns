#include "SupernovaExperiment.h"

#include <MAD/Element.h>
#include <MAD/Material.h>
using namespace MAD;

#include <NEUS/SupernovaModel.h>
using namespace NEUS;

#include <UNIC/Units.h>
#include <UNIC/Constants.h>
using namespace UNIC;

#include <TF1.h>
#include <TH2D.h>
#include <TAxis.h>
#include <TMath.h>
using namespace TMath;

//______________________________________________________________________________
//

SupernovaExperiment::SupernovaExperiment(
      Material *material, SupernovaModel *model) : TNamed(),
   fMaterial(material), fModel(model), fMass(0), fDistance(0), fThreshold(0)
{
   for (UShort_t i=0; i<7; i++) {
      fXSxNe[i]=0;
      fXSxN2[i]=0;
      fNevt[i]=0;
      fHN2[i]=0;
      fHNt[i]=0;
   }
}

//______________________________________________________________________________
//

Double_t SupernovaExperiment::XSxNe(Double_t *x, Double_t *parameter)
{
   Double_t Ev = x[0]*MeV; // neutrino energy
   Double_t Er = parameter[0]*keV; // nuclear recoil energy
   UShort_t type = static_cast<UShort_t>(parameter[1]); // type of neutrino

   Element *element = (Element*) fMaterial->GetElement();
   Double_t dXS = element->CNNSdXS(Er, Ev);

   if (type<1 || type>6)
      return dXS*(fModel->Ne(1,Ev/MeV) + fModel->Ne(2,Ev/MeV) + 4*fModel->Ne(3,Ev/MeV))/MeV;

   return dXS * fModel->Ne(type,Ev/MeV)/MeV;
}

//______________________________________________________________________________
//

Double_t SupernovaExperiment::XSxN2(Double_t *x, Double_t *parameter)
{
   Double_t Ev = x[0]*MeV; // neutrino energy
   Double_t Er = parameter[0]*keV; // nuclear recoil energy
   UShort_t type = static_cast<UShort_t>(parameter[1]); // type of neutrino
   Double_t time = parameter[2]*second;

   Element *element = (Element*) fMaterial->GetElement();
   Double_t dXS = element->CNNSdXS(Er, Ev);

   if (type<1 || type>6)
      return dXS*(fModel->N2(1,time/second,Ev/MeV)/second/MeV
            + fModel->N2(2,time/second,Ev/MeV)/second/MeV 
            + 4*fModel->N2(3,time/second,Ev/MeV))/second/MeV;

   return dXS * fModel->N2(type,time/second,Ev/MeV)/second/MeV;
}

//______________________________________________________________________________
//

Double_t SupernovaExperiment::FuncN(Double_t *x, Double_t *parameter)
{
   if (!fMaterial) {
      Warning("Nevt", "Please set targe material!");
      return 0;
   }
   if (fMaterial->Nelements()!=1) {
      Warning("Nevt", "Can only handle material with one element!");
      return 0;
   }
   if (!fModel) {
      Warning("Nevt", "Please set source of neutrino!");
      return 0;
   }

   Double_t Er = x[0]*keV; // nuclear recoil energy
   Double_t maxEv = parameter[0]*MeV; // max neutrino energy
   UShort_t type = static_cast<UShort_t>(parameter[1]); // type of neutrino

   Element *element = (Element*) fMaterial->GetElement();
   Double_t atomicMass  = element->A();
   Double_t nNuclei = fMass/atomicMass*Avogadro;
   Double_t area = 4*pi*fDistance/hbarc*fDistance/hbarc;
   Double_t minEv = (Er + Sqrt(2*element->M()*Er))/2;

   TF1 *f = FXSxNe(type,Er,minEv,maxEv);
   return nNuclei/area*1e50*f->Integral(minEv/MeV, maxEv/MeV)*keV;
}

//______________________________________________________________________________
//

Double_t SupernovaExperiment::N2(UShort_t type, Double_t time, Double_t Enr)
{
   if (!fMaterial) {
      Warning("Nevt", "Please set targe material!");
      return 0;
   }
   if (fMaterial->Nelements()!=1) {
      Warning("Nevt", "Can only handle material with one element!");
      return 0;
   }
   if (!fModel) {
      Warning("Nevt", "Please set source of neutrino!");
      return 0;
   }

   Element *element = (Element*) fMaterial->GetElement();
   Double_t atomicMass  = element->A();
   Double_t nNuclei = fMass/atomicMass*Avogadro;
   Double_t area = 4*pi*fDistance*fDistance;

   Double_t minEv = (Enr + Sqrt(2*element->M()*Enr))/2;
   Double_t maxEv = 82.5*MeV; // max neutrino energy

   TF1 *f = FXSxN2(type,time,Enr,minEv,maxEv);
   return nNuclei/area*1e50*f->Integral(minEv/MeV, maxEv/MeV)*keV*second;
}

//______________________________________________________________________________
//

Double_t SupernovaExperiment::Nevt(UShort_t type,
      Double_t nEr, Double_t maxEv)
{
   Double_t x[1] = {nEr/keV};
   Double_t p[2] = {maxEv/MeV, type};
   return FuncN(x,p);
}

//______________________________________________________________________________
//

TF1* SupernovaExperiment::FNevt(UShort_t type, Double_t maxEr, Double_t maxEv)
{
   if (fNevt[type]) {
      Double_t min, max;
      fNevt[type]->GetRange(min,max);
      if (max!=maxEr/keV) {
         Info("FNevt","Reset range of recoil energy.");
         fNevt[type]->SetRange(0,maxEr/keV);
      }
      if (fNevt[type]->GetParameter(0)!=maxEv/MeV) {
         Info("FNevt","Reset maximal neutrino energy.");
         fNevt[type]->SetParameter(0,maxEv/MeV);
      }
      return fNevt[type];
   }

   fNevt[type] = new TF1(Form("fNevt_%s_%s_%f_%d",
            fMaterial->GetName(), fModel->GetName(), fMass, type),
         this, &SupernovaExperiment::FuncN, 0, maxEr/keV,2);
   fNevt[type]->SetParameter(0,maxEv/MeV);
   fNevt[type]->SetParameter(1,type);

   if (type==0) {
      fNevt[type]->SetTitle(Form(
               "%s;nuclear recoil energy [keV];total events/(keV #times %.0f kg)",
               fModel->GetTitle(), fMass/kg));
      fNevt[type]->SetLineColor(kGray+2);
      fNevt[type]->SetLineWidth(2);
   } else {
      fNevt[type]->SetTitle(Form(
               "%s, type of neutrino: %d;nuclear recoil energy [keV];total events/(keV #times %.0f kg)",
               fModel->GetTitle(), type, fMass/kg));
      fNevt[type]->SetLineColor(type);
   }
   return fNevt[type];
}

//______________________________________________________________________________
//

TF1* SupernovaExperiment::FXSxNe(UShort_t type, Double_t nEr,
      Double_t minEv, Double_t maxEv)
{
   if (fXSxNe[type]) {
      fXSxNe[type]->SetRange(minEv/MeV,maxEv/MeV);
      fXSxNe[type]->SetParameter(0,nEr/keV);
      return fXSxNe[type];
   }

   fXSxNe[type] = new TF1(Form("fXSxNe%s%s%f%d",
            fModel->GetName(), fMaterial->GetName(), fMass, type), this, 
         &SupernovaExperiment::XSxNe, minEv/MeV, maxEv/MeV,2);
   fXSxNe[type]->SetParameter(0,nEr/keV);
   fXSxNe[type]->SetParameter(1,type);

   if (type==0) {
      fXSxNe[type]->SetTitle(Form(
               "%s, target: %s, recoil energy: %.1f keV;neutrino energy [MeV];1/MeV^{4}", 
               fModel->GetName(), fMaterial->GetTitle(), 
               nEr/keV));
      fXSxNe[type]->SetLineColor(kGray+2);
      fXSxNe[type]->SetLineWidth(2);
   } else {
      fXSxNe[type]->SetTitle(Form(
               "neutrino %d from %s, target: %s;neutrino energy [MeV];1/MeV^{4}", 
               type, fModel->GetName(), fMaterial->GetTitle()));
      fXSxNe[type]->SetLineColor(type);
   }
   return fXSxNe[type];
}

//______________________________________________________________________________
//

TF1* SupernovaExperiment::FXSxN2(UShort_t type, Double_t time, 
      Double_t Enr, Double_t minEv, Double_t maxEv)
{
   if (fXSxN2[type]) {
      fXSxN2[type]->SetRange(minEv/MeV,maxEv/MeV);
      fXSxN2[type]->SetParameter(0,Enr/keV);
      fXSxN2[type]->SetParameter(2,time/second);
      return fXSxN2[type];
   }

   fXSxN2[type] = new TF1(Form("fXSxN2%s%s%f%d",
            fModel->GetName(), fMaterial->GetName(), fMass, type), this, 
         &SupernovaExperiment::XSxN2, minEv/MeV, maxEv/MeV,3);
   fXSxN2[type]->SetParameter(0,Enr/keV);
   fXSxN2[type]->SetParameter(1,type);
   fXSxN2[type]->SetParameter(2,time/second);

   if (type==0) {
      fXSxN2[type]->SetTitle(Form(
               "%s, target: %s, recoil energy: %.1f keV, time: %.1f second;neutrino energy [MeV];1/MeV^{4}", 
               fModel->GetName(), fMaterial->GetTitle(), Enr/keV, time/second));
      fXSxN2[type]->SetLineColor(kGray+2);
      fXSxN2[type]->SetLineWidth(2);
   } else {
      fXSxN2[type]->SetTitle(Form(
               "neutrino %d from %s, target: %s, recoil energy : %.1f, time: %.1f second;neutrino energy [MeV];1/MeV^{4}", 
               type, fModel->GetName(), fMaterial->GetTitle(),
               Enr/keV, time/second));
      fXSxN2[type]->SetLineColor(type);
   }
   return fXSxN2[type];
}

//______________________________________________________________________________
//

TH1D* SupernovaExperiment::HNevt(UShort_t type, Double_t maxEr, Double_t maxEv)
{
   TH1D *h = (TH1D*) FNevt(type,maxEr,maxEv)->GetHistogram();
   return h;
}

//______________________________________________________________________________
//

TH1D* SupernovaExperiment::HXSxNe(UShort_t type, Double_t nEr,
      Double_t minEv, Double_t maxEv)
{
   TH1D *h = (TH1D*) FXSxNe(type, nEr, minEv, maxEv)->GetHistogram();
   return h;
}

//______________________________________________________________________________
//

void SupernovaExperiment::Clear(Option_t *option)
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
      if (fXSxN2[i]) {
         delete fXSxN2[i];
         fXSxN2[i]=NULL;
      }
      if (fHN2[i]) {
         delete fHN2[i];
         fHN2[i]=NULL;
      }
      if (fHNt[i]) {
         delete fHNt[i];
         fHNt[i]=NULL;
      }
   }
}

//______________________________________________________________________________
//

TH2D* SupernovaExperiment::HN2(UShort_t type)
{
   if (fHN2[type]) return fHN2[type];

   Int_t nbinst=fModel->HN2()->GetXaxis()->GetNbins();
   const Double_t *tbins = fModel->HN2()->GetXaxis()->GetXbins()->GetArray();
   fHN2[type] = new TH2D(Form("hN2%d",type),"",nbinst,tbins,100,0,50);

   for (Int_t ix=1; ix<=fHN2[type]->GetNbinsX(); ix++) {
      for (Int_t iy=1; iy<=fHN2[type]->GetNbinsY(); iy++) {
         Double_t t = fHN2[type]->GetXaxis()->GetBinCenter(ix);
         Double_t e = fHN2[type]->GetYaxis()->GetBinCenter(iy);
         Double_t nevt = N2(type,t*second,e*keV);
         fHN2[type]->SetBinContent(ix, iy, nevt);
      }
   }
   fHN2[type]->SetStats(0);
   fHN2[type]->GetXaxis()->SetTitle("time [second]");
   fHN2[type]->GetYaxis()->SetTitle("nuclear recoil energy [keV]");
   fHN2[type]->SetTitle(Form("number of events / (%.0f kg)",fMass/kg));

   return fHN2[type];
}

//______________________________________________________________________________
//

TH1D* SupernovaExperiment::HNt(UShort_t type)
{
   if (fHNt[type]) return fHNt[type];

   Int_t nbinst=fModel->HN2()->GetXaxis()->GetNbins();
   const Double_t *tbins = fModel->HN2()->GetXaxis()->GetXbins()->GetArray();
   fHNt[type] = new TH1D(Form("hNt%d",type),"",nbinst,tbins);

   for (Int_t ix=1; ix<=fHNt[type]->GetNbinsX(); ix++) {
      Double_t nevt=0;
      for (Int_t iy=1; iy<=HN2(type)->GetNbinsY(); iy++) {
         Double_t dn = HN2(type)->GetBinContent(ix,iy);
         Double_t er = HN2(type)->GetYaxis()->GetBinWidth(iy);
         nevt+=dn*er;
      }
      fHNt[type]->SetBinContent(ix, nevt);
   }
   fHNt[type]->SetStats(0);
   fHNt[type]->GetXaxis()->SetTitle("time [second]");
   fHNt[type]->GetYaxis()->SetTitle(Form(
            "rate of events [Hz/(%.0f kg)]",fMass/kg));

   return fHNt[type];
}

