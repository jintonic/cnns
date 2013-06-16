#ifndef SCINTILLATIONDETECTOR_H
#define SCINTILLATIONDETECTOR_H

#include "Detector.h"

class ScintillationDetector : public Detector
{
   protected:
      Double_t fLightYield;

   public:
      ScintillationDetector() : Detector(), fLightYield(0) {};
      ScintillationDetector(const char *name, const char *title) :
         Detector(name, title), fLightYield(0) {};

      virtual ~ScintillationDetector() {};

      void SetLightYield(Double_t yield) { fLightYield=yield; }
      Double_t LightYield() { return fLightYield; }

      ClassDef(ScintillationDetector,1);
};

#endif
