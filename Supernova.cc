#include "Supernova.h"

#include <CLHEP/Units/SystemOfUnits.h>
using namespace CLHEP;

//______________________________________________________________________________
//

Double_t Supernova::N2(UShort_t type, Double_t energy, Double_t time)
{
   if (!fModel) return 0;
   TH2D *h = fModel->H2N(type);
   return h->Interpolate(energy/MeV, time/second);
}

//______________________________________________________________________________
//

Double_t Supernova::Ne(UShort_t type, Double_t energy)
{
   if (!fModel) return 0;
   TH1D *h = fModel->HNe(type);
   return h->Interpolate(energy/MeV);
}

//______________________________________________________________________________
//
