#ifndef SUPERNOVA_H
#define SUPERNOVA_H

#include "Source.h"

class TArrayD;

namespace NEUS { class SupernovaModel; }

class Supernova : public Source
{
   protected:
      NEUS::SupernovaModel *fModel;

   public:
      Supernova(NEUS::SupernovaModel *model=0) 
         : Source(), fModel(0) { if (model) SetModel(model); }
      virtual ~Supernova() {};

      void SetModel(NEUS::SupernovaModel *model);
      NEUS::SupernovaModel* Model() { return fModel; }

      Double_t N2(UShort_t type, Double_t time, Double_t energy);
      Double_t Ne(UShort_t type, Double_t energy);

      Int_t NbinsT();
      const TArrayD* BinEdgesT();

      ClassDef(Supernova,1);
};

#endif
