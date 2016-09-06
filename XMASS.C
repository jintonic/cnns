#include "SupernovaExperiment.h"
#include "XMASS835kg.h"
using namespace CNNS;

#include <NEUS/NakazatoModel.h>
#include <NEUS/LivermoreModel.h>
using namespace NEUS;

#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>

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

   // model close to Betelgeuse
   NakazatoModel *betelgeuse = new NakazatoModel(13,0.02,100);
   betelgeuse->LoadData("../neus");

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
   xmass4sn->Distance=10*kpc; // galaxy center
   xmass4sn->Distance=196.22*pc; // Betelgeuse


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
      hN0[j]->SetLineColor(kBlue);
      hN1[j]->SetLineColor(kRed);

      hN0[j]->Draw();
      hN1[j]->Draw("same");

      if (j==0) {
         leg->AddEntry(hN0[j],"All events above 1 keV","l");
         leg->AddEntry(hN1[j],"Observable events","l");
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

   xmass4sn->SetSupernovaModel(model2001);
   xmass4sn->HNevtT(0)->GetXaxis()->SetRangeUser(1.8e-2,17.9012);
   TH1 *hT0 = xmass4sn->HNevtT(0)->DrawCopy();
   hT0->SetLineColor(kBlue);
   TH1D *hT1 = xmass4sn->HNevtT(0,kTRUE);
   hT1->SetLineColor(kRed);
   hT1->Draw("same");

   leg->Draw();
   can->Print("XMASS.ps");

   xmass4sn->SetSupernovaModel(model3003);
   xmass4sn->HNevtT(0)->GetXaxis()->SetRangeUser(1.8e-2,17.9012);
   hT0 = xmass4sn->HNevtT(0)->DrawCopy();
   hT0->SetLineColor(kBlue);
   hT1 = xmass4sn->HNevtT(0,kTRUE);
   hT1->SetLineColor(kRed);
   hT1->Draw("same");

   leg->Draw();
   can->Print("XMASS.ps");

   xmass4sn->SetSupernovaModel(betelgeuse);
   xmass4sn->HNevtT(0)->GetXaxis()->SetRangeUser(1.8e-2,17.9012);
   hT0 = xmass4sn->HNevtT(0)->DrawCopy();
   hT0->SetLineColor(kBlue);
   hT1 = xmass4sn->HNevtT(0,kTRUE);
   hT1->SetLineColor(kRed);
   hT1->Draw("same");

   leg->Draw();
   can->Print("XMASS.ps");

   can->Print("XMASS.ps]");

   // generate plot for svx paper
   gROOT->SetStyle("Plain");
   gStyle->SetLegendBorderSize(0);
   gStyle->SetLegendFont(22);
   gStyle->SetTitleBorderSize(0);
   gStyle->SetLabelFont(22,"XYZ");
   gStyle->SetTitleFont(22,"H");
   gStyle->SetTitleFont(22,"XYZ");
   gStyle->SetLabelSize(0.05,"XYZ");
   gStyle->SetTitleSize(0.05,"XYZ");
   gStyle->SetTitleOffset(1.1,"Y");
   gStyle->SetTitleOffset(-0.5,"Z");
   gStyle->SetPadRightMargin(0.02);
   gStyle->SetPadLeftMargin(0.12);
   gStyle->SetPadTopMargin(0.02);
   gStyle->SetPadBottomMargin(0.11);
   gROOT->ForceStyle();
   TCanvas *c = new TCanvas;

   xmass4sn->SetSupernovaModel(totani);
   xmass4sn->Distance=196.22*pc; // Betelgeuse

   hT0 = xmass4sn->HNevtT(0)->DrawCopy();
   hT0->GetXaxis()->SetRangeUser(0,10.5);
   hT0->GetYaxis()->SetTitle("(number of events)/second/(832 kg)");
   hT0->SetTitle("");
   hT0->SetLineColor(kBlue);
   hT1 = xmass4sn->HNevtT(0,kTRUE);
   hT1->SetLineColor(kRed);
   hT1->Draw("same");

   TLegend *l = new TLegend(0.4,0.6,0.88,0.88);
   l->SetHeader("Livermore Model");
   l->AddEntry(hT0,"All events above 1 keV","l");
   l->AddEntry(hT1,"Observable events","l");
   l->Draw();
   c->Print("rate0.eps");
   c->Print("rate0.pdf");
}
