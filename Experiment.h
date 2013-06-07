#ifndef EXPERIMENT_H
#define EXPERIMENT_H

#include <TObject.h>
class TF1;
class TH1D;

#include <UNIC/Units.h>

class Source;

namespace MAD { class Material; }

class Experiment : public TObject
{
   protected:
      Source *fSource; // source of neutrino
      MAD::Material *fMaterial; // target material
      Double_t fMass;      // mass of target material
      Double_t fThreshold; // energy threshold of detector

      TF1 *fXSxNe[7]; // dXS(Ev) * Ne(Ev)
      TF1 *fNevt[7]; // Nevt in TF1 format

      Double_t XSxNe(Double_t *x, Double_t *parameter); // function of dXS * Ne
      Double_t FuncN(Double_t *x, Double_t *parameter); // function of Nevt

   public:
      Experiment(MAD::Material *material=0, Source *source=0);
      virtual ~Experiment() { Clear(); } 

      // in case of reactor, it gives Nevt/keVnr/second
      // in case of supernova, it gives Nevt/keVnr in 20 second
      Double_t Nevt(UShort_t type,
            Double_t nuclearRecoilEnergy, 
            Double_t maxNeutrinoEnergy=82.5*UNIC::MeV);

      void SetTargetMass(Double_t m) { fMass=m; }
      Double_t TargetMass() { return fMass; }

      void SetThreshold(Double_t threshold) { fThreshold=threshold; }
      Double_t Threshold() { return fThreshold; }

      void SetTargetMaterial(MAD::Material *material) { fMaterial=material; }
      MAD::Material* TargetMaterial() { return fMaterial; }

      void SetNeutrinoSource(Source *source) { fSource=source; }
      Source* NeutrinoSource() { return fSource; }

      TF1* FXSxNe(UShort_t type=1, 
            Double_t nuclearRecoilEnergy=5*UNIC::keV,
            Double_t minNeutrinoEnergy= 2.5*UNIC::MeV,
            Double_t maxNeutrinoEnergy=82.5*UNIC::MeV);

      TF1* FNevt(UShort_t type=1,
            Double_t maxNuclearRecoilEnergy=50*UNIC::keV,
            Double_t maxNeutrinoEnergy=82.5*UNIC::MeV);

      TH1D* HXSxNe(UShort_t type=1,
            Double_t nuclearRecoilEnergy=5*UNIC::keV,
            Double_t minNeutrinoEnergy= 2.5*UNIC::MeV,
            Double_t maxNeutrinoEnergy=82.5*UNIC::MeV);

      TH1D* HNevt(UShort_t type=1,
            Double_t maxNuclearRecoilEnergy=50*UNIC::keV,
            Double_t maxNeutrinoEnergy=82.5*UNIC::MeV);

      /**
       * Delete internal TF1 objects.
       */
      void Clear(Option_t *option="");

      ClassDef(Experiment,1);
};

#endif
