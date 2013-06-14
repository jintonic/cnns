#include "SupernovaExperiment.h"

#include <NEUS/NakazatoModel.h>
#include <NEUS/LivermoreModel.h>
using namespace NEUS;

#include <MAD/NaturalXe.h>
#include <MAD/LiquidXenon.h>
using namespace MAD;

#include <UNIC/Units.h>
using namespace UNIC;

#include <TF1.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>

int main ()
{
   // the weakest sn in Nakazato model
   NakazatoModel *model2001 = new NakazatoModel(20,0.02,100);
   model2001->LoadData("../neus");

   // the brightest sn in Nakazato model
   NakazatoModel *model3003 = new NakazatoModel(30,0.02,300);
   model3003->LoadData("../neus");

   // black hole in Nakazato model
   NakazatoModel *blackHole = new NakazatoModel(30,0.004);
   blackHole->LoadData("../neus");

   // Totani's Livermore model
   LivermoreModel *totani = new LivermoreModel;
   totani->LoadData("../total");

   // Divari's approximation
   LivermoreModel *divari = new LivermoreModel;
   divari->UseDivariData();


   // target material
   NaturalXe *natXe = new NaturalXe;
   LiquidXenon *LXe = new LiquidXenon;
   LXe->AddElement(natXe,1);


   // set up experiment
   SupernovaExperiment *xmass = new SupernovaExperiment(LXe);
   xmass->SetDistance(196.22*pc); // Betelgeuse
   xmass->SetDistance(10*kpc); // galaxy center
   xmass->SetTargetMass(835*kg);
   xmass->SetThreshold(0.3*keV);


   // draw results
   TCanvas *can = new TCanvas;
   can->Print("XMASS.ps[");

   // dXS()
   xmass->SetSupernovaModel(divari);
   TH1D *h0 = xmass->HXSxNe(0,5*keV);
   TH1D *h1 = xmass->HXSxNe(1,5*keV);
   TH1D *h2 = xmass->HXSxNe(2,5*keV);
   TH1D *h3 = xmass->HXSxNe(3,5*keV);

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

   can->Print("XMASS.ps");

   // Divari approximation
   Double_t maxEr = 50*keV;
   h0 = xmass->HNevt(0,maxEr);
   h1 = xmass->HNevt(1,maxEr);
   h2 = xmass->HNevt(2,maxEr);
   h3 = xmass->HNevt(3,maxEr);

   h0->GetYaxis()->SetRangeUser(0,10);
   h0->Draw();
   h2->Draw("same");
   h1->Draw("same");
   h3->Draw("same");

   leg->Draw();
   can->Print("XMASS.ps");

   Printf("number of events in Divari approximation: %.1f",
         xmass->TotalNevt(maxEr));

   // Totani's Livermore model
   xmass->SetSupernovaModel(totani);
   xmass->Clear(); // clear internal functions
   h0 = xmass->HNevt(0,maxEr);
   h1 = xmass->HNevt(1,maxEr);
   h2 = xmass->HNevt(2,maxEr);
   h3 = xmass->HNevt(3,maxEr);

   h0->GetYaxis()->SetRangeUser(0,10);
   h0->Draw();
   h2->Draw("same");
   h1->Draw("same");
   h3->Draw("same");

   leg->Draw();
   can->Print("XMASS.ps");

   Printf("number of events in Livermore model: %.1f",
         xmass->TotalNevt(maxEr));

   // weakest Nakazato Model
   xmass->SetSupernovaModel(model2001);
   xmass->Clear(); // clear internal functions
   h0 = xmass->HNevt(0,maxEr);
   h1 = xmass->HNevt(1,maxEr);
   h2 = xmass->HNevt(2,maxEr);
   h3 = xmass->HNevt(3,maxEr);

   h0->GetYaxis()->SetRangeUser(0,10);
   h0->Draw();
   h2->Draw("same");
   h1->Draw("same");
   h3->Draw("same");

   leg->Draw();
   can->Print("XMASS.ps");

   Printf("number of events in Nakazato model 2001: %.1f",
         xmass->TotalNevt(maxEr));

   // brightest Nakazato Model
   xmass->SetSupernovaModel(model3003);
   xmass->Clear(); // clear internal functions
   h0 = xmass->HNevt(0,maxEr);
   h1 = xmass->HNevt(1,maxEr);
   h2 = xmass->HNevt(2,maxEr);
   h3 = xmass->HNevt(3,maxEr);

   h0->GetYaxis()->SetRangeUser(0,10);
   h0->Draw();
   h2->Draw("same");
   h1->Draw("same");
   h3->Draw("same");

   leg->Draw();
   can->Print("XMASS.ps");

   Printf("number of events in Nakazato model 3003: %.1f",
         xmass->TotalNevt(maxEr));

   // black hole in Nakazato Model
   xmass->SetSupernovaModel(blackHole);
   xmass->Clear(); // clear internal functions
   h0 = xmass->HNevt(0,maxEr);
   h1 = xmass->HNevt(1,maxEr);
   h2 = xmass->HNevt(2,maxEr);
   h3 = xmass->HNevt(3,maxEr);

   h0->GetYaxis()->SetRangeUser(0,10);
   h0->Draw();
   h2->Draw("same");
   h1->Draw("same");
   h3->Draw("same");

   leg->Draw();
   can->Print("XMASS.ps");

   Printf("number of events in black hole: %.1f",
         xmass->TotalNevt(maxEr));

   // time dependent event rate
   xmass->SetSupernovaModel(totani);
   xmass->Clear();
   TH2D *hT = xmass->HN2(0);

   hT->GetXaxis()->SetRangeUser(1.2e-2,17.9012);
   //hT->GetXaxis()->SetRangeUser(1.2e-2,0.9012);
   //hT->GetYaxis()->SetRangeUser(0,0.5);
   can->SetLogx();
   can->SetLogz();
   hT->Draw("colz");
   can->Print("XMASS.ps");

   //can->SetLogy();
   xmass->HNt(0)->GetXaxis()->SetRangeUser(1.8e-2,17.9012);
   xmass->HNt(0)->Draw();
   can->Print("XMASS.ps");

   can->Print("XMASS.ps]");
}
