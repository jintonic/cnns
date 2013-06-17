#include "SupernovaExperiment.h"
#include "XMASS835kg.h"

#include <NEUS/NakazatoModel.h>
#include <NEUS/LivermoreModel.h>
using namespace NEUS;

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

   // set up detector
   XMASS835kg *xmass = new XMASS835kg;

   // set up experiment
   SupernovaExperiment *xmass4sn = new SupernovaExperiment(xmass);
   xmass4sn->SetDistance(196.22*pc); // Betelgeuse
   xmass4sn->SetDistance(10*kpc); // galaxy center


   // draw results
   TCanvas *can = new TCanvas;
   can->Print("XMASS.ps[");

   // dXS()
   xmass4sn->SetSupernovaModel(divari);
   TH1D *h0 = xmass4sn->HXSxNe(0,5*keV);
   TH1D *h1 = xmass4sn->HXSxNe(1,5*keV);
   TH1D *h2 = xmass4sn->HXSxNe(2,5*keV);
   TH1D *h3 = xmass4sn->HXSxNe(3,5*keV);

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
   h0 = xmass4sn->HNevtE(0);
   h1 = xmass4sn->HNevtE(1);
   h2 = xmass4sn->HNevtE(2);
   h3 = xmass4sn->HNevtE(3);

   TH1 *hN0[5], *hN1[5];
   hN0[0] = (TH1*)h0->Clone("hN00");
   hN1[0] = (TH1*)h0->Clone("hN10");
   h0->Draw();
   h2->Draw("same");
   h1->Draw("same");
   h3->Draw("same");

   leg->Draw();
   can->Print("XMASS.ps");

   // Totani's Livermore model
   xmass4sn->SetSupernovaModel(totani);
   h0 = xmass4sn->HNevtE(0);
   h1 = xmass4sn->HNevtE(1);
   h2 = xmass4sn->HNevtE(2);
   h3 = xmass4sn->HNevtE(3);

   hN0[1] = (TH1*)h0->Clone("hN01");
   hN1[1] = (TH1*)h0->Clone("hN11");
   h0->Draw();
   h2->Draw("same");
   h1->Draw("same");
   h3->Draw("same");

   leg->Draw();
   can->Print("XMASS.ps");

   // weakest Nakazato Model
   xmass4sn->SetSupernovaModel(model2001);
   h0 = xmass4sn->HNevtE(0);
   h1 = xmass4sn->HNevtE(1);
   h2 = xmass4sn->HNevtE(2);
   h3 = xmass4sn->HNevtE(3);

   hN0[2] = (TH1*)h0->Clone("hN02");
   hN1[2] = (TH1*)h0->Clone("hN12");
   h0->Draw();
   h2->Draw("same");
   h1->Draw("same");
   h3->Draw("same");

   leg->Draw();
   can->Print("XMASS.ps");

   // brightest Nakazato Model
   xmass4sn->SetSupernovaModel(model3003);
   h0 = xmass4sn->HNevtE(0);
   h1 = xmass4sn->HNevtE(1);
   h2 = xmass4sn->HNevtE(2);
   h3 = xmass4sn->HNevtE(3);

   hN0[3] = (TH1*)h0->Clone("hN03");
   hN1[3] = (TH1*)h0->Clone("hN13");
   h0->Draw();
   h2->Draw("same");
   h1->Draw("same");
   h3->Draw("same");

   leg->Draw();
   can->Print("XMASS.ps");

   // black hole in Nakazato Model
   xmass4sn->SetSupernovaModel(blackHole);
   h0 = xmass4sn->HNevtE(0);
   h1 = xmass4sn->HNevtE(1);
   h2 = xmass4sn->HNevtE(2);
   h3 = xmass4sn->HNevtE(3);

   hN0[4] = (TH1*)h0->Clone("hN04");
   hN1[4] = (TH1*)h0->Clone("hN14");
   h0->Draw();
   h2->Draw("same");
   h1->Draw("same");
   h3->Draw("same");

   leg->Draw();
   can->Print("XMASS.ps");

   // fold in detection efficiency
   can->SetLogx(0);
   xmass->HEff()->Draw("E0");
   can->Print("XMASS.ps");

   leg->Clear();
   Double_t nevt[5] = {0}, content;
   for (Int_t j=0; j<5; j++) { // loop over models
      for (Int_t i=1; i<=hN0[j]->GetNbinsX(); i++) { // loop over Enr
         content = hN1[j]->GetBinContent(i) // nevt(Enr)/keV x eff(Enr)
            * xmass->Efficiency(hN1[j]->GetBinCenter(i)*keV);
         hN1[j]->SetBinContent(i, content);
         nevt[j]+=content*hN1[j]->GetBinWidth(i);
      }
      hN0[j]->Draw();
      hN1[j]->SetLineColor(kRed);
      hN1[j]->Draw("same");
      if (j==0) {
         leg->AddEntry(hN0[j],"All events","l");
         leg->AddEntry(hN1[j],"Detected events","l");
      }
      leg->Draw();
      can->Print("XMASS.ps");
   }

   Printf("number of events in Divari approximation: %.1f", nevt[0]);
   Printf("number of events in Livermore model: %.1f", nevt[1]);
   Printf("number of events in Nakazato model 2001: %.1f", nevt[2]);
   Printf("number of events in Nakazato model 3003: %.1f", nevt[3]);
   Printf("number of events in black hole: %.1f", nevt[4]); 

   // time dependent event rate
   xmass4sn->SetSupernovaModel(totani);
   TH2D *hN = xmass4sn->HNevt2(0);

   hN->GetXaxis()->SetRangeUser(1.2e-2,17.9012);
   can->SetLogx();
   can->SetLogz();
   hN->Draw("colz");
   can->Print("XMASS.ps");

   xmass4sn->HNevtT(0)->GetXaxis()->SetRangeUser(1.8e-2,17.9012);
   TH1 *hc = xmass4sn->HNevtT(0)->DrawCopy();

   TH1D *hT = xmass4sn->HNevtT(0);
   hT->SetLineColor(kBlue);
   hT->Draw("same");

   leg->Clear();
   leg->SetHeader("Detector threshold:");
   leg->AddEntry(hc,"0 keVee","l");
   leg->AddEntry(hT,"0.3 keVee","l");
   leg->Draw();
   can->Print("XMASS.ps");

   can->Print("XMASS.ps]");
}
