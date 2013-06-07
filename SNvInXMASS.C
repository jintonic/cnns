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

   // Divari approximation
   h0 = xmass->HNevt(0,50*keV);
   h1 = xmass->HNevt(1,50*keV);
   h2 = xmass->HNevt(2,50*keV);
   h3 = xmass->HNevt(3,50*keV);

   h0->GetYaxis()->SetRangeUser(0,10);
   h0->Draw();
   h2->Draw("same");
   h1->Draw("same");
   h3->Draw("same");

   leg->Draw();
   can->Print("SNvInXMASS.ps");

   Printf("number of events in Divari approximation: %f",
         h0->Integral(3*2+1,50*2)*0.5);

   // Totani's Livermore model
   sn->SetModel(totani);
   xmass->Clear(); // clear internal functions
   h0 = xmass->HNevt(0,50*keV);
   h1 = xmass->HNevt(1,50*keV);
   h2 = xmass->HNevt(2,50*keV);
   h3 = xmass->HNevt(3,50*keV);

   h0->GetYaxis()->SetRangeUser(0,10);
   h0->Draw();
   h2->Draw("same");
   h1->Draw("same");
   h3->Draw("same");

   leg->Draw();
   can->Print("SNvInXMASS.ps");

   Printf("number of events in Livermore model: %f",
         h0->Integral(3*2+1,50*2)*0.5);

   // weakest Nakazato Model
   sn->SetModel(model2001);
   xmass->Clear(); // clear internal functions
   h0 = xmass->HNevt(0,50*keV);
   h1 = xmass->HNevt(1,50*keV);
   h2 = xmass->HNevt(2,50*keV);
   h3 = xmass->HNevt(3,50*keV);

   h0->GetYaxis()->SetRangeUser(0,10);
   h0->Draw();
   h2->Draw("same");
   h1->Draw("same");
   h3->Draw("same");

   leg->Draw();
   can->Print("SNvInXMASS.ps");

   Printf("number of events in Nakazato model 2001: %f",
         h0->Integral(3*2+1,50*2)*0.5);

   // weakest Nakazato Model
   sn->SetModel(model2001);
   xmass->Clear(); // clear internal functions
   h0 = xmass->HNevt(0,50*keV);
   h1 = xmass->HNevt(1,50*keV);
   h2 = xmass->HNevt(2,50*keV);
   h3 = xmass->HNevt(3,50*keV);

   h0->GetYaxis()->SetRangeUser(0,10);
   h0->Draw();
   h2->Draw("same");
   h1->Draw("same");
   h3->Draw("same");

   leg->Draw();
   can->Print("SNvInXMASS.ps");

   Printf("number of events in Nakazato model 2001: %f",
         h0->Integral(3*2+1,50*2)*0.5);

   // brightest Nakazato Model
   sn->SetModel(model3003);
   xmass->Clear(); // clear internal functions
   h0 = xmass->HNevt(0,50*keV);
   h1 = xmass->HNevt(1,50*keV);
   h2 = xmass->HNevt(2,50*keV);
   h3 = xmass->HNevt(3,50*keV);

   h0->GetYaxis()->SetRangeUser(0,10);
   h0->Draw();
   h2->Draw("same");
   h1->Draw("same");
   h3->Draw("same");

   leg->Draw();
   can->Print("SNvInXMASS.ps");

   Printf("number of events in Nakazato model 3003: %f",
         h0->Integral(3*2+1,50*2)*0.5);

   // black hole in Nakazato Model
   sn->SetModel(blackHole);
   xmass->Clear(); // clear internal functions
   h0 = xmass->HNevt(0,50*keV);
   h1 = xmass->HNevt(1,50*keV);
   h2 = xmass->HNevt(2,50*keV);
   h3 = xmass->HNevt(3,50*keV);

   h0->GetYaxis()->SetRangeUser(0,10);
   h0->Draw();
   h2->Draw("same");
   h1->Draw("same");
   h3->Draw("same");

   leg->Draw();
   can->Print("SNvInXMASS.ps");

   Printf("number of events in black hole: %f",
         h0->Integral(3*2+1,50*2)*0.5);

   can->Print("SNvInXMASS.ps]");
}
