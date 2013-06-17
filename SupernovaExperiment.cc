#include "SupernovaExperiment.h"
#include "Detector.h"

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
      Detector *detector, SupernovaModel *model) : TNamed(),
   fDetector(detector), fModel(model), fDistance(0)
{
   for (UShort_t i=0; i<SupernovaModel::fgNtype; i++) {
      fFXSxNe[i]=0;
      fFXSxN2[i]=0;
      fHNevt2[i]=0;
      fHNevtT[i]=0;
      fHNevtE[i]=0;
   }
}

//______________________________________________________________________________
//

Double_t SupernovaExperiment::XSxNe(Double_t *x, Double_t *parameter)
{
   Double_t Ev = x[0]; // neutrino energy
   Double_t Er = parameter[0]; // nuclear recoil energy
   UShort_t type = static_cast<UShort_t>(parameter[1]); // type of neutrino

   Element *element = fDetector->TargetMaterial()->GetElement();
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

   Element *element = fDetector->TargetMaterial()->GetElement();
   Double_t dXS = element->CNNSdXS(Er*keV, Ev*MeV);

   if (type==0)
      return dXS*(fModel->N2(1,time,Ev) + fModel->N2(2,time,Ev)
            + 4*fModel->N2(3,time,Ev))/second/MeV;

   return dXS * fModel->N2(type,time,Ev)/second/MeV;
}

//______________________________________________________________________________
//

Double_t SupernovaExperiment::NevtE(UShort_t type, Double_t Enr)
{
   if (!fDetector->TargetMaterial()) {
      Warning("NevtE", "Please set targe material!");
      return 0;
   }
   if (fDetector->TargetMaterial()->Nelements()!=1) {
      Warning("NevtE", "Can only handle material with one element!");
      return 0;
   }
   if (!fModel) {
      Warning("NevtE", "Please set supernova model!");
      return 0;
   }

   Element *element = fDetector->TargetMaterial()->GetElement();
   Double_t atomicMass  = element->A();
   Double_t nNuclei = fDetector->TargetMass()/atomicMass*Avogadro;
   Double_t area = 4*pi*fDistance/hbarc*fDistance/hbarc;

   Double_t maxEv = fModel->EMax()*MeV; // max neutrino energy
   Double_t minEv = (Enr + Sqrt(2*element->M()*Enr))/2;
   if (minEv<fModel->EMin()*MeV) {
      Warning("NevtE","Requested neutrino energy is too small.");
      Warning("NevtE","Reset it to the minimal energy provided by NEUS.");
      minEv=fModel->EMin()*MeV;
   }

   TF1 *f = FXSxNe(type,Enr/keV);
   return nNuclei/area*1e50*f->Integral(minEv/MeV, maxEv/MeV)*keV;
}

//______________________________________________________________________________
//

Double_t SupernovaExperiment::Nevt2(UShort_t type, Double_t time, Double_t Enr)
{
   if (!fDetector->TargetMaterial()) {
      Warning("Nevt2", "Please set targe material!");
      return 0;
   }
   if (fDetector->TargetMaterial()->Nelements()!=1) {
      Warning("Nevt2", "Can only handle material with one element!");
      return 0;
   }
   if (!fModel) {
      Warning("Nevt2", "Please set supernova model!");
      return 0;
   }

   Element *element = fDetector->TargetMaterial()->GetElement();
   Double_t atomicMass  = element->A();
   Double_t nNuclei = fDetector->TargetMass()/atomicMass*Avogadro;
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

TF1* SupernovaExperiment::FXSxNe(UShort_t type, Double_t Enr)
{
   if (fFXSxNe[type]) {
      fFXSxNe[type]->SetParameter(0,Enr);
      return fFXSxNe[type];
   }

   fFXSxNe[type] = new TF1(Form("fFXSxNe%s%s%f%d", fModel->GetName(),
            fDetector->TargetMaterial()->GetName(), fDetector->TargetMass(),
            type), this, &SupernovaExperiment::XSxNe, 0., fModel->EMax(),2);
   fFXSxNe[type]->SetParameter(0,Enr);
   fFXSxNe[type]->SetParameter(1,type);

   if (type==0) {
      fFXSxNe[type]->SetTitle(Form(
               "%s, target: %s, recoil energy: %.1f keV;neutrino energy [MeV];1/MeV^{4}", 
               fModel->GetName(), fDetector->TargetMaterial()->GetTitle(), Enr));
      fFXSxNe[type]->SetLineColor(kGray+2);
      fFXSxNe[type]->SetLineWidth(2);
   } else {
      fFXSxNe[type]->SetTitle(Form(
               "neutrino %d from %s, target: %s;neutrino energy [MeV];1/MeV^{4}", 
               type,fModel->GetName(),fDetector->TargetMaterial()->GetTitle()));
      fFXSxNe[type]->SetLineColor(type);
   }
   return fFXSxNe[type];
}

//______________________________________________________________________________
//

TF1* SupernovaExperiment::FXSxN2(UShort_t type, Double_t time, Double_t Enr)
{
   if (fFXSxN2[type]) {
      fFXSxN2[type]->SetParameter(0,Enr);
      fFXSxN2[type]->SetParameter(2,time);
      return fFXSxN2[type];
   }

   fFXSxN2[type] = new TF1(Form("fFXSxN2%s%s%f%d", fModel->GetName(),
            fDetector->TargetMaterial()->GetName(), fDetector->TargetMass(),
            type), this, &SupernovaExperiment::XSxN2, 0., fModel->EMax(),3);
   fFXSxN2[type]->SetParameter(0,Enr);
   fFXSxN2[type]->SetParameter(1,type);
   fFXSxN2[type]->SetParameter(2,time);

   if (type==0) {
      fFXSxN2[type]->SetTitle(Form(
               "%s, target: %s, recoil energy: %.1f keV, time: %.1f second;neutrino energy [MeV];1/MeV^{4}", 
               fModel->GetName(), 
               fDetector->TargetMaterial()->GetTitle(), Enr, time));
      fFXSxN2[type]->SetLineColor(kGray+2);
      fFXSxN2[type]->SetLineWidth(2);
   } else {
      fFXSxN2[type]->SetTitle(Form(
               "neutrino %d from %s, target: %s, recoil energy : %.1f, time: %.1f second;neutrino energy [MeV];1/MeV^{4}", 
               type, fModel->GetName(), 
               fDetector->TargetMaterial()->GetTitle(), Enr, time));
      fFXSxN2[type]->SetLineColor(type);
   }
   return fFXSxN2[type];
}

//______________________________________________________________________________
//

Double_t SupernovaExperiment::Nevt()
{
   TH1D* h = HNevtE(0);
   Double_t minEr = fDetector->Threshold();
   Int_t startBin=0, endBin=0;
   for (Int_t i=1; i<=h->GetNbinsX(); i++) {
      startBin=i;
      if (h->GetBinLowEdge(i)*keV>=minEr) break;
   }
   for (Int_t i=1; i<=h->GetNbinsX(); i++) {
      endBin=i;
   }
   return h->Integral(startBin,endBin,"width");
}

//______________________________________________________________________________
//

TH1D* SupernovaExperiment::HXSxNe(UShort_t type, Double_t Enr)
{
   if (type>6) {
      Warning("HXSxNe","Type of neutrinos must be in 0, 1, 2, 3, 4, 5, 6!");
      Warning("HXSxNe","Return NULL pointer!");
      return 0;
   }
   TH1D *h = (TH1D*) FXSxNe(type, Enr/keV)->GetHistogram();
   return h;
}

//______________________________________________________________________________
//

void SupernovaExperiment::Clear(Option_t *option)
{
   for (UShort_t i=0; i<SupernovaModel::fgNtype; i++) {
      if (fFXSxNe[i]) {
         delete fFXSxNe[i];
         fFXSxNe[i]=NULL;
      }
      if (fFXSxN2[i]) {
         delete fFXSxN2[i];
         fFXSxN2[i]=NULL;
      }
      if (fHNevt2[i]) {
         delete fHNevt2[i];
         fHNevt2[i]=NULL;
      }
      if (fHNevtT[i]) {
         delete fHNevtT[i];
         fHNevtT[i]=NULL;
      }
      if (fHNevtE[i]) {
         delete fHNevtE[i];
         fHNevtE[i]=NULL;
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

   TString name = Form("hNevt2-%d-%f", type, fDetector->Threshold());
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
      else if (e<9.9999) de=0.25;
      else if (e<20.) de=0.5;
      else de=1;
      ebins[nbinse]=e;
      nbinse++;
      e+=de;
   }
   ebins[nbinse]=e;

   // create histogram
   Double_t minEr = fDetector->Threshold();
   //Info("HNevt2","Create HNevt2-%d with threshold %.3f keVnr",type,minEr/keV);
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
   fHNevt2[type]->SetTitle(Form("number of events / (%.0f kg)",
            fDetector->TargetMass()/kg));

   return fHNevt2[type];
}

//______________________________________________________________________________
//

TH1D* SupernovaExperiment::HNevtT(UShort_t type, Bool_t detectableOnly)
{
   TH2D *h = HNevt2(type);

   TString name = Form("hNevtT-%d-%f-%d", 
         type, fDetector->Threshold(), detectableOnly);
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
         Double_t de = h->GetYaxis()->GetBinWidth(iy);
         if (detectableOnly) nevt+=dn*de*fDetector->Efficiency(
               h->GetYaxis()->GetBinCenter(iy)*keV);
         else nevt+=dn*de;
      }
      fHNevtT[type]->SetBinContent(ix, nevt);
   }
   fHNevtT[type]->SetStats(0);
   fHNevtT[type]->SetTitle(Form("%s",fModel->GetTitle()));
   fHNevtT[type]->SetXTitle("time [second]");
   fHNevtT[type]->SetYTitle(Form("rate of events [Hz/(%.0f kg)]",
            fDetector->TargetMass()/kg));
   fHNevtT[type]->GetYaxis()->SetTitleOffset(1.3);

   return fHNevtT[type];
}

//______________________________________________________________________________
//

TH1D* SupernovaExperiment::HNevtE(UShort_t type)
{
   if (type>6) {
      Warning("HNevtE","Type of neutrinos must be in 0, 1, 2, 3, 4, 5, 6!");
      Warning("HNevtE","Return NULL pointer!");
      return 0;
   }
   if (!fDetector->TargetMaterial()) {
      Warning("HNevtE", "Please set targe material!");
      return 0;
   }

   TString name = Form("hNevtE-%d-%f", type, fDetector->Threshold());
   if (fHNevtE[type]) {
      if (name.CompareTo(fHNevtE[type]->GetName())==0) return fHNevtE[type];
      else delete fHNevtE[type];
   }

   // define bins
   Int_t nbinse=0;
   Double_t e=0, de, ebins[200];
   while (e<50.) {
      if (e<4.9999) de=0.1;
      else if (e<9.9999) de=0.25;
      else if (e<20.) de=0.5;
      else de=1;
      ebins[nbinse]=e;
      nbinse++;
      e+=de;
   }
   ebins[nbinse]=e;

   // create histogram
   Double_t minEr = fDetector->Threshold();
   //Info("HNevtE","Create HNevtE-%d with threshold %.3f keVnr",type,minEr/keV);
   fHNevtE[type] = new TH1D(name.Data(),"",nbinse,ebins);

   // fill histogram
   for (Int_t ix=1; ix<=fHNevtE[type]->GetNbinsX(); ix++) {
      Double_t e = fHNevtE[type]->GetXaxis()->GetBinCenter(ix);
      if (e*keV<minEr) continue; // skip events below threshold
      Double_t nevt = NevtE(type,e*keV);
      fHNevtE[type]->SetBinContent(ix, nevt);
   }
   fHNevtE[type]->SetStats(0);
   fHNevtE[type]->SetTitle(Form("%s",fModel->GetTitle()));
   fHNevtE[type]->SetXTitle("nuclear recoil energy [keV]");
   fHNevtE[type]->SetYTitle(Form(
            "number of events / (keV #times %.0f kg)",
            fDetector->TargetMass()/kg));
   fHNevtE[type]->GetYaxis()->SetTitleOffset(1.3);
   if (type==0) fHNevtE[type]->SetLineColor(kGray+2);
   else fHNevtE[type]->SetLineColor(type);

   return fHNevtE[type];
}
