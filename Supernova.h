#ifndef SUPERNOVA_H
#define SUPERNOVA_H

#include "Source.h"

namespace NEUS { class SupernovaModel; }

class Supernova : public Source
{
   protected:
      NEUS::SupernovaModel *fModel;
      Bool_t fUseFermiDiracApproximation;

   public:
      Supernova(NEUS::SupernovaModel *model=0) 
         : Source(), fModel(0), fUseFermiDiracApproximation(kFALSE)
      { if (model) SetModel(model); }
      virtual ~Supernova() {};

      void SetModel(NEUS::SupernovaModel *model);
      NEUS::SupernovaModel* Model() { return fModel; }

      void UseFermiDiracApproximation() { fUseFermiDiracApproximation=kTRUE; }

      Double_t N2(UShort_t type, Double_t time, Double_t energy);
      Double_t Ne(UShort_t type, Double_t energy);
      Double_t Nt(UShort_t type, Double_t time);

      ClassDef(Supernova,1);
};

#endif
