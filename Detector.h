#ifndef CNNS_DETECTOR_H
#define CNNS_DETECTOR_H

#include <TNamed.h>

namespace MAD { class Material; }
namespace CNNS {
   class Detector;

   static const Double_t MeV=1.;
   static const Double_t keV=MeV/1000.;

   static const Double_t mm = 1.;
   static const Double_t cm = 10.*mm;
   static const Double_t meter = 1e3*mm;
   static const Double_t pc = 3.0856775807e+16*meter;
   static const Double_t kpc = 1e3*pc;
   static const Double_t ly = 9460730472580800000.*mm;

   static const Double_t ns = 1.;
   static const Double_t sec = 1.e+9 *ns;

   static const Double_t joule = 6.24150e+12*MeV;
   static const Double_t kg = joule*sec*sec/meter/meter;   

   static const Double_t number_of_photoelectrons = 1;
   static const Double_t PE = number_of_photoelectrons;

   static const Double_t Avogadro = 6.02214179e+23;
   static const Double_t pi = 3.14159265;
   static const Double_t hbarc = 197.32705e-12*MeV*mm;
}

class CNNS::Detector : public TNamed
{
   public:
      MAD::Material *TargetMaterial;
      Double_t TargetMass;
      Double_t EnergyThreshold;

   public:
      Detector() : TNamed(), TargetMaterial(0), TargetMass(0), EnergyThreshold(0) {};
      Detector(const char *name, const char *title) : TNamed(name, title),
      TargetMaterial(0), TargetMass(0), EnergyThreshold(0) {};
      virtual ~Detector() {};

      virtual Double_t Efficiency(Double_t Enr) { return 1.; }

      ClassDef(Detector,1);
};

#endif
