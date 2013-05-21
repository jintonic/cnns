#ifndef SUPERNOVA_HH
#define SUPERNOVA_HH

#include <NakazatoModel.h>

#include "Source.h"

class Supernova : public Source
{
   protected:
      NakazatoModel *fModel;

   public:
      Supernova(NakazatoModel *model=0) 
         : Source(), fModel(0) { if (model) SetModel(model); }
      virtual ~Supernova() {};

      void SetModel(NakazatoModel *model);
      NakazatoModel* Model() { return fModel; }

      Double_t N2(UShort_t type, Double_t energy, Double_t time);
      Double_t Ne(UShort_t type, Double_t energy);

      ClassDef(Supernova,1);
};

#endif // SUPERNOVA_HH
