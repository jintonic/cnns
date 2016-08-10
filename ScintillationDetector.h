#ifndef CNNS_SCINTILLATIONDETECTOR_H
#define CNNS_SCINTILLATIONDETECTOR_H

#include "Detector.h"

namespace CNNS { class ScintillationDetector; }

class CNNS::ScintillationDetector : public CNNS::Detector
{
   public:
      Double_t LightYield;

   public:
      ScintillationDetector() : Detector(), LightYield(0) {};
      ScintillationDetector(const char *name, const char *title) :
         Detector(name, title), LightYield(0) {};

      virtual ~ScintillationDetector() {};

      ClassDef(ScintillationDetector,1);
};

#endif
