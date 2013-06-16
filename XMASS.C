#include "SupernovaExperiment.h"

#include <NEUS/NakazatoModel.h>
#include <NEUS/LivermoreModel.h>
using namespace NEUS;

#include <MAD/NaturalXe.h>
#include <MAD/LiquidXenon.h>
using namespace MAD;

#include <UNIC/Units.h>
using namespace UNIC;

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
   xmass->SetThreshold(0.*keV);


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
   h0 = xmass->HNevtE(0);
   h1 = xmass->HNevtE(1);
   h2 = xmass->HNevtE(2);
   h3 = xmass->HNevtE(3);

   h0->Draw();
   h2->Draw("same");
   h1->Draw("same");
   h3->Draw("same");

   leg->Draw();
   can->Print("XMASS.ps");

   // Totani's Livermore model
   xmass->SetSupernovaModel(totani);
   h0 = xmass->HNevtE(0);
   h1 = xmass->HNevtE(1);
   h2 = xmass->HNevtE(2);
   h3 = xmass->HNevtE(3);

   h0->Draw();
   h2->Draw("same");
   h1->Draw("same");
   h3->Draw("same");

   leg->Draw();
   can->Print("XMASS.ps");

   // weakest Nakazato Model
   xmass->SetSupernovaModel(model2001);
   h0 = xmass->HNevtE(0);
   h1 = xmass->HNevtE(1);
   h2 = xmass->HNevtE(2);
   h3 = xmass->HNevtE(3);

   h0->Draw();
   h2->Draw("same");
   h1->Draw("same");
   h3->Draw("same");

   leg->Draw();
   can->Print("XMASS.ps");

   // brightest Nakazato Model
   xmass->SetSupernovaModel(model3003);
   h0 = xmass->HNevtE(0);
   h1 = xmass->HNevtE(1);
   h2 = xmass->HNevtE(2);
   h3 = xmass->HNevtE(3);

   h0->Draw();
   h2->Draw("same");
   h1->Draw("same");
   h3->Draw("same");

   leg->Draw();
   can->Print("XMASS.ps");

   // black hole in Nakazato Model
   xmass->SetSupernovaModel(blackHole);
   h0 = xmass->HNevtE(0);
   h1 = xmass->HNevtE(1);
   h2 = xmass->HNevtE(2);
   h3 = xmass->HNevtE(3);

   h0->Draw();
   h2->Draw("same");
   h1->Draw("same");
   h3->Draw("same");

   leg->Draw();
   can->Print("XMASS.ps");

   // time dependent event rate
   xmass->SetSupernovaModel(totani);
   TH2D *hN = xmass->HNevt2(0);

   hN->GetXaxis()->SetRangeUser(1.2e-2,17.9012);
   can->SetLogx();
   can->SetLogz();
   hN->Draw("colz");
   can->Print("XMASS.ps");

   xmass->HNevtT(0)->GetXaxis()->SetRangeUser(1.8e-2,17.9012);
   TH1 *hc = xmass->HNevtT(0)->DrawCopy();

   xmass->SetThreshold(0.3*keV);
   TH1D *hT = xmass->HNevtT(0);
   hT->SetLineColor(kBlue);
   hT->Draw("same");

   leg->Clear();
   leg->SetHeader("Detector threshold:");
   leg->AddEntry(hc,"0 keVee","l");
   leg->AddEntry(hT,"0.3 keVee","l");
   leg->Draw();
   can->Print("XMASS.ps");

   can->Print("XMASS.ps]");

   xmass->SetSupernovaModel(divari);
   Printf("number of events in Divari approximation: %.1f", xmass->Nevt());
   xmass->SetSupernovaModel(totani);
   Printf("number of events in Livermore model: %.1f", xmass->Nevt());
   xmass->SetSupernovaModel(model2001);
   Printf("number of events in Nakazato model 2001: %.1f", xmass->Nevt());
   xmass->SetSupernovaModel(model3003);
   Printf("number of events in Nakazato model 3003: %.1f", xmass->Nevt());
   xmass->SetSupernovaModel(blackHole);
   Printf("number of events in black hole: %.1f", xmass->Nevt()); 
}
