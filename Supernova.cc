#include <NEUS/SupernovaModel.h>
using namespace NEUS;

#include "Supernova.h"

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
   //Printf("type: %d energy: %f MeV, Ne: %e", type, energy/MeV,
         //fModel->Ne(type, energy/MeV));
   return fModel->Ne(type, energy/MeV);
}

//______________________________________________________________________________
//

Double_t Supernova::Nt(UShort_t type, Double_t time)
{
   if (!fModel) return 0;
   return fModel->Ne(type, time/second);
}

//______________________________________________________________________________
//

void Supernova::SetModel(SupernovaModel *model)
{
   fModel = model;
   fName  = model->GetName();
   fTitle = model->GetTitle();
}
