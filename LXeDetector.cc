#include "LXeDetector.h"
using namespace CNNS;

#include <MAD/LiquidXenon.h>
using namespace MAD;

#include <TGraphErrors.h>

//______________________________________________________________________________
//

void LXeDetector::SetThreshold(Double_t threshold)
{
   if (!TargetMaterial) {
      Warning("SetThreshold","No target material available!");
      Warning("SetThreshold","Simply set EnergyThreshold to %f.", threshold);
      EnergyThreshold = threshold;
      return;
   }
   TString material(TargetMaterial->GetName());
   if (material.CompareTo("LXe")!=0) {
      Warning("SetThreshold","Target material is not LXe!");
      Warning("SetThreshold","Simply set EnergyThreshold to %f.", threshold);
      EnergyThreshold = threshold;
      return;
   }
   if (LightYield==0.) {
      Warning("SetThreshold","Light yield equals to zero!");
      Warning("SetThreshold","Simply set EnergyThreshold to %f.", threshold);
      EnergyThreshold = threshold;
      return;
   }
   LiquidXenon *LXe = (LiquidXenon*) TargetMaterial;
   EnergyThreshold = LXe->EnrPE(LightYield)->Eval(threshold)*keV;
}

//______________________________________________________________________________
//
