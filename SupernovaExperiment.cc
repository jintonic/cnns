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
      fFNevtE[i]=0;
      fHNevt2[i]=0;
      fHNevtT[i]=0;
   }
}

//______________________________________________________________________________
//

Double_t SupernovaExperiment::XSxNe(Double_t *x, Double_t *parameter)
{
   Double_t Ev = x[0]; // neutrino energy
   Double_t Er = parameter[0]; // nuclear recoil energy
   UShort_t type = static_cast<UShort_t>(parameter[1]); // type of neutrino

   Element *element = fMaterial->GetElement();
   Double_t dXS = element->CNNSdXS(Er*keV, Ev*MeV);

   if (type==0)
      return dXS*(fModel->Ne(1,Ev) + fModel->Ne(2,Ev) + 4*fModel->Ne(3,Ev))/MeV;

   return dXS * fModel->Ne(type,Ev)/MeV;
}

//______________________________________________________________________________
//

Double_t SupernovaExperiment::XSxN2(Double_t *x, Double_t *parameter)
{
   Double_t Ev = x[0]; // neutrino energy
   Double_t Er = parameter[0]; // nuclear recoil energy
   UShort_t type = static_cast<UShort_t>(parameter[1]); // type of neutrino
   Double_t time = parameter[2];

   Element *element = fMaterial->GetElement();
   Double_t dXS = element->CNNSdXS(Er*keV, Ev*MeV);

   if (type==0)
      return dXS*(fModel->N2(1,time,Ev) + fModel->N2(2,time,Ev)
            + 4*fModel->N2(3,time,Ev))/second/MeV;

   return dXS * fModel->N2(type,time,Ev)/second/MeV;
}

//______________________________________________________________________________
//

Double_t SupernovaExperiment::NevtE(Double_t *x, Double_t *parameter)
{
   if (!fMaterial) {
      Warning("NevtE", "Please set targe material!");
      return 0;
   }
   if (fMaterial->Nelements()!=1) {
      Warning("NevtE", "Can only handle material with one element!");
      return 0;
   }
   if (!fModel) {
      Warning("NevtE", "Please set supernova model!");
      return 0;
   }

   Double_t Er = x[0]*keV; // nuclear recoil energy
   UShort_t type = static_cast<UShort_t>(parameter[0]); // type of neutrino

   Double_t maxEv = fModel->EMax()*MeV; // max neutrino energy

   Element *element = fMaterial->GetElement();
   Double_t atomicMass  = element->A();
   Double_t nNuclei = fMass/atomicMass*Avogadro;
   Double_t area = 4*pi*fDistance/hbarc*fDistance/hbarc;
   Double_t minEv = (Er + Sqrt(2*element->M()*Er))/2;
   if (minEv<fModel->EMin()*MeV) {
      Warning("NevtE","Requested neutrino energy is too small.");
      Warning("NevtE","Reset it to the minimal energy provided by NEUS.");
      minEv=fModel->EMin()*MeV;
   }

   TF1 *f = FXSxNe(type,x[0]);
   return nNuclei/area*1e50*f->Integral(minEv/MeV, maxEv/MeV)*keV;
}

//______________________________________________________________________________
//

Double_t SupernovaExperiment::Nevt2(UShort_t type, Double_t time, Double_t Enr)
{
   if (!fMaterial) {
      Warning("Nevt2", "Please set targe material!");
      return 0;
   }
   if (fMaterial->Nelements()!=1) {
      Warning("Nevt2", "Can only handle material with one element!");
      return 0;
   }
   if (!fModel) {
      Warning("Nevt2", "Please set supernova model!");
      return 0;
   }

   Element *element = fMaterial->GetElement();
   Double_t atomicMass  = element->A();
   Double_t nNuclei = fMass/atomicMass*Avogadro;
   Double_t area = 4*pi*fDistance/hbarc*fDistance/hbarc;

   Double_t maxEv = fModel->EMax()*MeV; // max neutrino energy
   Double_t minEv = (Enr + Sqrt(2*element->M()*Enr))/2;
   if (minEv<fModel->EMin()*MeV) {
      Warning("Nevt2","Requested neutrino energy is too small.");
      Warning("Nevt2","Reset it to the minimal energy provided by NEUS.");
      minEv=fModel->EMin()*MeV;
   }

   TF1 *f = FXSxN2(type,time/second,Enr/keV);
   return nNuclei/area*1e50*f->Integral(minEv/MeV, maxEv/MeV)*keV*second;
}

//______________________________________________________________________________
//

TF1* SupernovaExperiment::FNevtE(UShort_t type, Double_t maxEr)
{
   if (fFNevtE[type]) {
      Double_t min, max;
      fFNevtE[type]->GetRange(min,max);
      if (max!=maxEr) {
         Info("FNevtE","Reset range of recoil energy.");
         fFNevtE[type]->SetRange(0,maxEr);
      }
      return fFNevtE[type];
   }

   fFNevtE[type] = new TF1(Form("fNevtE_%s_%s_%f_%d",
            fMaterial->GetName(), fModel->GetName(), fMass, type),
         this, &SupernovaExperiment::NevtE, 0, maxEr,1);
   fFNevtE[type]->SetParameter(0,type);

   if (type==0) {
      fFNevtE[type]->SetTitle(Form(
               "%s;nuclear recoil energy [keV];total events/(keV #times %.0f kg)",
               fModel->GetTitle(), fMass/kg));
      fFNevtE[type]->SetLineColor(kGray+2);
      fFNevtE[type]->SetLineWidth(2);
   } else {
      fFNevtE[type]->SetTitle(Form(
               "%s, type of neutrino: %d;nuclear recoil energy [keV];total events/(keV #times %.0f kg)",
               fModel->GetTitle(), type, fMass/kg));
      fFNevtE[type]->SetLineColor(type);
   }
   return fFNevtE[type];
}

//______________________________________________________________________________
//

TF1* SupernovaExperiment::FXSxNe(UShort_t type, Double_t nEr)
{
   if (fXSxNe[type]) {
      fXSxNe[type]->SetParameter(0,nEr);
      return fXSxNe[type];
   }

   fXSxNe[type] = new TF1(Form("fXSxNe%s%s%f%d",
            fModel->GetName(), fMaterial->GetName(), fMass, type), this, 
         &SupernovaExperiment::XSxNe, 0., fModel->EMax(),2);
   fXSxNe[type]->SetParameter(0,nEr);
   fXSxNe[type]->SetParameter(1,type);

   if (type==0) {
      fXSxNe[type]->SetTitle(Form(
               "%s, target: %s, recoil energy: %.1f keV;neutrino energy [MeV];1/MeV^{4}", 
               fModel->GetName(), fMaterial->GetTitle(), nEr));
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

TF1* SupernovaExperiment::FXSxN2(UShort_t type, Double_t time, Double_t Enr)
{
   if (fXSxN2[type]) {
      fXSxN2[type]->SetParameter(0,Enr);
      fXSxN2[type]->SetParameter(2,time);
      return fXSxN2[type];
   }

   fXSxN2[type] = new TF1(Form("fXSxN2%s%s%f%d",
            fModel->GetName(), fMaterial->GetName(), fMass, type), this, 
         &SupernovaExperiment::XSxN2, 0., fModel->EMax(),3);
   fXSxN2[type]->SetParameter(0,Enr);
   fXSxN2[type]->SetParameter(1,type);
   fXSxN2[type]->SetParameter(2,time);

   if (type==0) {
      fXSxN2[type]->SetTitle(Form(
               "%s, target: %s, recoil energy: %.1f keV, time: %.1f second;neutrino energy [MeV];1/MeV^{4}", 
               fModel->GetName(), fMaterial->GetTitle(), Enr, time));
      fXSxN2[type]->SetLineColor(kGray+2);
      fXSxN2[type]->SetLineWidth(2);
   } else {
      fXSxN2[type]->SetTitle(Form(
               "neutrino %d from %s, target: %s, recoil energy : %.1f, time: %.1f second;neutrino energy [MeV];1/MeV^{4}", 
               type, fModel->GetName(), fMaterial->GetTitle(),
               Enr, time));
      fXSxN2[type]->SetLineColor(type);
   }
   return fXSxN2[type];
}

//______________________________________________________________________________
//

TH1D* SupernovaExperiment::HNevtE(UShort_t type, Double_t maxEr)
{
   if (type>6) {
      Warning("HNevtE","Type of neutrinos must be in 0, 1, 2, 3, 4, 5, 6!");
      Warning("HNevtE","Return NULL pointer!");
      return 0;
   }
   FNevtE(type,maxEr/keV)->SetNpx(200);
   TH1D *h = (TH1D*) FNevtE(type,maxEr/keV)->GetHistogram();
   return h;
}

//______________________________________________________________________________
//

Double_t SupernovaExperiment::Nevt(Double_t maxEr)
{
   TH1D* h = HNevtE(0,maxEr);
   Double_t minEr = fMaterial->Enr(Threshold());
   Int_t startBin=0, endBin=0;
   for (Int_t i=1; i<=h->GetNbinsX(); i++) {
      startBin=i;
      if (h->GetBinLowEdge(i)*keV>=minEr) break;
   }
   for (Int_t i=1; i<=h->GetNbinsX(); i++) {
      if (h->GetBinLowEdge(i)*keV>=maxEr) break;
      endBin=i;
   }
   return h->Integral(startBin,endBin,"width");
}

//______________________________________________________________________________
//

TH1D* SupernovaExperiment::HXSxNe(UShort_t type, Double_t nEr)
{
   if (type>6) {
      Warning("HXSxNe","Type of neutrinos must be in 0, 1, 2, 3, 4, 5, 6!");
      Warning("HXSxNe","Return NULL pointer!");
      return 0;
   }
   TH1D *h = (TH1D*) FXSxNe(type, nEr/keV)->GetHistogram();
   return h;
}

//______________________________________________________________________________
//

void SupernovaExperiment::Clear(Option_t *option)
{
   for (UShort_t i=0; i<7; i++) {
      if (fFNevtE[i]) {
         delete fFNevtE[i];
         fFNevtE[i]=NULL;
      }
      if (fXSxNe[i]) {
         delete fXSxNe[i];
         fXSxNe[i]=NULL;
      }
      if (fXSxN2[i]) {
         delete fXSxN2[i];
         fXSxN2[i]=NULL;
      }
      if (fHNevt2[i]) {
         delete fHNevt2[i];
         fHNevt2[i]=NULL;
      }
      if (fHNevtT[i]) {
         delete fHNevtT[i];
         fHNevtT[i]=NULL;
      }
   }
}

//______________________________________________________________________________
//

TH2D* SupernovaExperiment::HNevt2(UShort_t type)
{
   if (type>6) {
      Warning("HNevt2","Type of neutrinos must be in 0, 1, 2, 3, 4, 5, 6!");
      Warning("HNevt2","Return NULL pointer!");
      return 0;
   }

   TString name = Form("hNevt2-%d-%f", type, fThreshold);
   if (fHNevt2[type]) {
      if (name.CompareTo(fHNevt2[type]->GetName())==0) return fHNevt2[type];
      else delete fHNevt2[type];
   }

   // define bins
   Int_t nbinst=fModel->HN2()->GetXaxis()->GetNbins();
   const Double_t *tbins = fModel->HN2()->GetXaxis()->GetXbins()->GetArray();
   Int_t nbinse=0;
   Double_t e=0, de, ebins[200];
   while (e<50.) {
      if (e<4.9999) de=0.1;
      else if (e<20.) de=1;
      else de=2;
      ebins[nbinse]=e;
      nbinse++;
      e+=de;
   }
   ebins[nbinse]=e;

   // create histogram
   Double_t minEr = fMaterial->Enr(Threshold());
   Info("HNevt2","Create HNevt2 with threshold %.3f keVnr",minEr/keV);
   fHNevt2[type] = new TH2D(name.Data(),"",nbinst,tbins,nbinse,ebins);

   // fill histogram
   for (Int_t ix=1; ix<=fHNevt2[type]->GetNbinsX(); ix++) {
      for (Int_t iy=1; iy<=fHNevt2[type]->GetNbinsY(); iy++) {
         Double_t t = fHNevt2[type]->GetXaxis()->GetBinCenter(ix);
         Double_t e = fHNevt2[type]->GetYaxis()->GetBinCenter(iy);
         if (e*keV<minEr) continue; // skip events below threshold
         Double_t nevt = Nevt2(type,t*second,e*keV);
         fHNevt2[type]->SetBinContent(ix, iy, nevt);
      }
   }
   fHNevt2[type]->SetStats(0);
   fHNevt2[type]->GetXaxis()->SetTitle("time [second]");
   fHNevt2[type]->GetYaxis()->SetTitle("nuclear recoil energy [keV]");
   fHNevt2[type]->SetTitle(Form("number of events / (%.0f kg)",fMass/kg));

   return fHNevt2[type];
}

//______________________________________________________________________________
//

TH1D* SupernovaExperiment::HNevtT(UShort_t type)
{
   TH2D *h = HNevt2(type);

   TString name = Form("hNevtT-%d-%f", type, fThreshold);
   if (fHNevtT[type]) {
      if (name.CompareTo(fHNevtT[type]->GetName())==0) return fHNevtT[type];
      else delete fHNevtT[type];
   }

   // create histogram
   Int_t nbinst=h->GetXaxis()->GetNbins();
   const Double_t *tbins = h->GetXaxis()->GetXbins()->GetArray();
   fHNevtT[type] = new TH1D(name.Data(),"",nbinst,tbins);

   // fill histogram
   for (Int_t ix=1; ix<=fHNevtT[type]->GetNbinsX(); ix++) {
      Double_t nevt=0;
      for (Int_t iy=1; iy<=h->GetNbinsY(); iy++) {
         Double_t dn = h->GetBinContent(ix,iy);
         Double_t er = h->GetYaxis()->GetBinWidth(iy);
         nevt+=dn*er;
      }
      fHNevtT[type]->SetBinContent(ix, nevt);
   }
   fHNevtT[type]->SetStats(0);
   fHNevtT[type]->SetXTitle("time [second]");
   fHNevtT[type]->SetYTitle(Form("rate of events [Hz/(%.0f kg)]",fMass/kg));

   return fHNevtT[type];
}

