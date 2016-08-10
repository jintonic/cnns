#ifndef CNNS_LXEDETECTOR_H
#define CNNS_LXEDETECTOR_H

#include "ScintillationDetector.h"

namespace CNNS { class LXeDetector; }

class CNNS::LXeDetector : public CNNS::ScintillationDetector
{
   public:
      LXeDetector() : ScintillationDetector() {};
      LXeDetector(const char *name, const char *title) :
         ScintillationDetector(name, title) {};

      virtual ~LXeDetector() {};

      void SetThreshold(Double_t threshold);

      ClassDef(LXeDetector,1);
};

#endif
