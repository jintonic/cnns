#include "Supernova.h"

#include <TH1D.h>

#include <NEUS/SupernovaModel.h>
using namespace NEUS;

#include <UNIC/Units.h>
using namespace UNIC;

//______________________________________________________________________________
//

Double_t Supernova::N2(UShort_t type, Double_t time, Double_t energy)
{
   if (!fModel) return 0;
   return fModel->N2(type, time/second, energy/MeV);
}

//______________________________________________________________________________
//

Double_t Supernova::Ne(UShort_t type, Double_t energy)
{
   if (!fModel) return 0;
   //return fModel->Ne(type, energy/MeV); // too slow for integration
   if (fModel->GetName()[0]=='D')
      return fModel->HNeFD(type)->Interpolate(energy/MeV);
   else
      return fModel->HNe(type, fModel->TMax())->Interpolate(energy/MeV);
}

//______________________________________________________________________________
//

Double_t Supernova::Nt(UShort_t type, Double_t time)
{
   if (!fModel) return 0;
   //return fModel->Ne(type, time/second); // too slow for integration
   return fModel->HNt(type)->Interpolate(time/second);
}

//______________________________________________________________________________
//

void Supernova::SetModel(SupernovaModel *model)
{
   fModel = model;
   fName  = model->GetName();
   fTitle = model->GetTitle();
}
