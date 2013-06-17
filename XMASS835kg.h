#ifndef XMASS835KG_H
#define XMASS835KG_H

#include "LXeDetector.h"

class TH1D;
namespace MAD {
   class NaturalXe;
   class LiquidXenon;
}

class XMASS835kg : public LXeDetector
{
   private:
      TH1D *fHTrgEff;
      MAD::NaturalXe *fNatXe;
      MAD::LiquidXenon *fLXe;

   public:
      XMASS835kg(const char *name="XMASS835kg",
            const char *title="XMASS 835kg");
      ~XMASS835kg();

      TH1D* HTrgEff();
      Double_t TriggerEfficiency(Double_t Enr);

      ClassDef(XMASS835kg,1);
};

#endif
