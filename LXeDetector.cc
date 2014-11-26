#include "LXeDetector.h"

#include <MAD/LiquidXenon.h>
using namespace MAD;

#include <UNIC/Units.h>
using namespace UNIC;

#include <TGraphErrors.h>

//______________________________________________________________________________
//

void LXeDetector::SetThreshold(Double_t threshold)
{
   if (!fMaterial) {
      Warning("SetThreshold","No target material available!");
      Warning("SetThreshold","Simply set fThreshold to threshold.");
      fThreshold = threshold;
      return;
   }
   TString material(fMaterial->GetName());
   if (material.CompareTo("LXe")!=0) {
      Warning("SetThreshold","Target material is not LXe!");
      Warning("SetThreshold","Simply set fThreshold to threshold.");
      fThreshold = threshold;
      return;
   }
   if (LightYield()==0.) {
      Warning("SetThreshold","Light yield equals to zero!");
      Warning("SetThreshold","Simply set fThreshold to threshold.");
      fThreshold = threshold;
      return;
   }
   LiquidXenon *LXe = (LiquidXenon*) fMaterial;
   fThreshold = LXe->EnrPE(LightYield())->Eval(threshold)*keV;
}

//______________________________________________________________________________
//
