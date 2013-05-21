#include "Experiment.h"

#include "Supernova.h"

#include <NakazatoModel.h>
#include <NaturalXe.h>
#include <LiquidXenon.h>

#include <TCanvas.h>

using namespace CLHEP;

#include <iostream>
using namespace std;

int main ()
{
   NakazatoModel *model1311 = new NakazatoModel(13,0.004,100);
   model1311->LoadIntegratedData("../snvd/integdata");
   model1311->LoadFullData("../snvd/intpdata");

   Supernova *sn = new Supernova(model1311);
   //sn->SetDistance(196.22*pc); // Betelgeuse
   sn->SetDistance(10*1000*pc);

   NaturalXe *natXe = new NaturalXe;
   LiquidXenon *LXe = new LiquidXenon;
   LXe->SetElement(natXe);

   Experiment *xmass = new Experiment(LXe, sn);
   xmass->SetTargetMass(835*kg);
   xmass->SetThreshold(0*keV);

   TCanvas *can = new TCanvas;
   can->Print("SNvInXMASS.ps[");

   xmass->FXSxNe(0)->Draw();
   xmass->FXSxNe(1)->Draw("same");
   xmass->FXSxNe(2)->Draw("same");
   xmass->FXSxNe(3)->Draw("same");
   can->Print("SNvInXMASS.ps");

   can->SetLogy();
   TF1 *f = xmass->FNevt(0,50*keV);
   //f->SetNpx(100);
   f->Draw();
   f->GetHistogram()->Draw("histsame");
   xmass->FNevt(1,50*keV)->Draw("same");
   xmass->FNevt(2,50*keV)->Draw("same");
   xmass->FNevt(3,50*keV)->Draw("same");
   can->Print("SNvInXMASS.ps");

   Printf("number of events: %f",f->GetHistogram()->Integral(7,100)*0.5);
   
   can->Print("SNvInXMASS.ps]");
}
