#include "Experiment.h"

#include "Supernova.h"

#include <NEUS/NakazatoModel.h>
#include <NEUS/LivermoreModel.h>
#include <MAD/NaturalXe.h>
#include <MAD/LiquidXenon.h>

using namespace UNIC;
using namespace MAD;
using namespace NEUS;

#include <TF1.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>

#include <iostream>
using namespace std;

int main ()
{
   // the weakest sn in Nakazato model
   NakazatoModel *model2001 = new NakazatoModel(20,0.02,100);
   model2001->SetDataLocation("../neus");
   model2001->LoadIntegratedData();
   model2001->LoadFullData();

   // the brightest sn in Nakazato model
   NakazatoModel *model3003 = new NakazatoModel(30,0.02,300);
   model3003->SetDataLocation("../neus");
   model3003->LoadIntegratedData();
   model3003->LoadFullData();

   // black hole in Nakazato model
   NakazatoModel *blackHole = new NakazatoModel(30,0.004);
   blackHole->SetDataLocation("../neus");
   blackHole->LoadIntegratedData();

   // Totani's Livermore model
   LivermoreModel *totani = new LivermoreModel;
   totani->SetDataLocation("../total");

   // Divari's approximation
   LivermoreModel *divari = new LivermoreModel;
   divari->SetDataLocation("../total");
   divari->UseDivariData();


   // set up experiment
   NaturalXe *natXe = new NaturalXe;
   LiquidXenon *LXe = new LiquidXenon;
   LXe->AddElement(natXe,1);

   Supernova *sn = new Supernova();
   //sn->SetDistance(196.22*pc); // Betelgeuse
   sn->SetDistance(10*kpc);

   Experiment *xmass = new Experiment(LXe, sn);
   xmass->SetTargetMass(835*kg);
   xmass->SetThreshold(0.3*keV);


   // draw results
   TCanvas *can = new TCanvas;
   can->Print("SNvInXMASS.ps[");

   sn->SetModel(divari);
   TH1D *h0 = xmass->HXSxNe(0);
   TH1D *h1 = xmass->HXSxNe(1);
   TH1D *h2 = xmass->HXSxNe(2);
   TH1D *h3 = xmass->HXSxNe(3);

   h0->Draw();
   h1->Draw("same");
   h2->Draw("same");
   h3->Draw("same");

   TLegend *leg = new TLegend(0.5,0.6,0.88,0.88);
   leg->SetFillColor(0);
   leg->SetBorderSize(0);
   leg->AddEntry(h0, "all flavors","l");
   leg->AddEntry(h1, "#nu_{e}","l");
   leg->AddEntry(h2, "#bar{#nu}_{e}","l");
   leg->AddEntry(h3, "#nu_{x}","l");
   leg->Draw();

   can->Print("SNvInXMASS.ps");

   //can->SetLogy();
   //TF1 *f = xmass->FNevt(0,50*keV);
   ////f->SetNpx(100);
   //f->Draw();
   ////f->GetHistogram()->Draw("histsame");
   //TF1 *fv1 = xmass->FNevt(1,50*keV);
   //fv1->Draw("same");
   //TF1 *fv2 = xmass->FNevt(2,50*keV);
   //fv2->Draw("same");
   //TF1 *fv3 = xmass->FNevt(3,50*keV);
   //fv3->Draw("same");

   //can->Print("SNvInXMASS.ps");

   //Printf("number of events: %f",f->GetHistogram()->Integral(7,100)*0.5);
   
   //can->SetLogy(0);
   //xmass1->FXSxNe(0)->Draw();
   //xmass1->FXSxNe(1)->Draw("same");
   //xmass1->FXSxNe(2)->Draw("same");
   //xmass1->FXSxNe(3)->Draw("same");
   //can->Print("SNvInXMASS.ps");

   //can->SetLogy();
   //TF1 *f = xmass1->FNevt(0,50*keV);
   //f->SetNpx(100);
   //f->Draw();
   //f->GetHistogram()->Draw("hist");
   //fv1 = xmass1->FNevt(1,50*keV);
   //fv1->Draw("same");
   //fv2 = xmass1->FNevt(2,50*keV);
   //fv2->Draw("same");
   //fv3 = xmass1->FNevt(3,50*keV);
   //fv3->Draw("same");

   //leg->Draw();

   //can->Print("SNvInXMASS.ps");

   //Printf("number of events: %f",f->GetHistogram()->Integral(7,100)*0.5);

   can->Print("SNvInXMASS.ps]");
}
