#ifndef LXEDETECTOR_H
#define LXEDETECTOR_H

#include "ScintillationDetector.h"

class LXeDetector : public ScintillationDetector
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
