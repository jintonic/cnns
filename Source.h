#ifndef SOURCE_H
#define SOURCE_H

#include <TNamed.h>

class Source : public TNamed
{
   protected:
      Float_t fDistance; // distance between detector and source

   public:
      Source() : TNamed(), fDistance(0) {};
      virtual ~Source() {};

      virtual Double_t N2(UShort_t type, Double_t time, Double_t energy)
      { return 0; }
      virtual Double_t Ne(UShort_t type, Double_t energy) { return 0; }
      virtual Double_t Nt(UShort_t type, Double_t time) { return 0; }

      void SetDistance(Double_t distance) { fDistance=distance; }
      Double_t Distance() { return fDistance; }

      ClassDef(Source,1);
};

#endif
