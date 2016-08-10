#ifndef CNNS_XMASS835KG_H
#define CNNS_XMASS835KG_H

#include "LXeDetector.h"

class TH1D;

namespace CNNS { class XMASS835kg; }

class CNNS::XMASS835kg : public CNNS::LXeDetector
{
   private:
      TH1D *fHEff;

   public:
      XMASS835kg(const char *name="XMASS835kg",
            const char *title="XMASS 835kg");
      ~XMASS835kg();

      TH1D* HEff();
      Double_t Efficiency(Double_t Enr);

      ClassDef(XMASS835kg,1);
};

#endif
