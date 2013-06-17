#ifndef DETECTOR_H
#define DETECTOR_H

#include <TNamed.h>

namespace MAD { class Material; }

class Detector : public TNamed
{
   protected:
      MAD::Material *fMaterial; // target material
      Double_t fMass;      // mass of target material
      Double_t fThreshold; // energy threshold of detector

   public:
      Detector() : TNamed(), fMaterial(0), fMass(0), fThreshold(0) {};
      Detector(const char *name, const char *title) : TNamed(name, title),
      fMaterial(0), fMass(0), fThreshold(0) {};
      virtual ~Detector() {};

      void SetTargetMaterial(MAD::Material *material) { fMaterial=material; }
      MAD::Material* TargetMaterial() { return fMaterial; }

      void SetTargetMass(Double_t mass) { fMass=mass; }
      Double_t TargetMass() { return fMass; }

      virtual void SetThreshold(Double_t threshold) { fThreshold=threshold; }
      Double_t Threshold() { return fThreshold; }

      virtual Double_t TriggerEfficiency(Double_t Enr) { return 1.; }

      ClassDef(Detector,1);
};

#endif
