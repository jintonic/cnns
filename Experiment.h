#ifndef EXPERIMENT_HH
#define EXPERIMENT_HH

class TGeoMaterial;
class Source;
#include <TF1.h>

#include <CLHEP/Units/SystemOfUnits.h>

class Experiment : public TObject
{
   protected:
      Source *fSource; // source of neutrino
      TGeoMaterial *fMaterial; // target material
      Double_t fMass;      // mass of target material
      Double_t fThreshold; // energy threshold of detector

      TF1 *fXSxNe[7]; // dXS(Ev) * Ne(Ev)
      TF1 *fNevt[7]; // Nevt in TF1 format

      Double_t XSxNe(Double_t *x, Double_t *parameter); // function of dXS * Ne
      Double_t FuncN(Double_t *x, Double_t *parameter); // function of Nevt

   public:
      Experiment(TGeoMaterial *material=0, Source *source=0);
      virtual ~Experiment(); 

      // in case of reactor, it gives Nevt/keVnr/second
      // in case of supernova, it gives Nevt/keVnr in 20 second
      Double_t Nevt(UShort_t type,
            Double_t nuclearRecoilEnergy, 
            Double_t maxNeutrinoEnergy=100*CLHEP::MeV);

      void SetTargetMass(Double_t m) { fMass=m; }
      Double_t TargetMass() { return fMass; }

      void SetThreshold(Double_t threshold) { fThreshold=threshold; }
      Double_t Threshold() { return fThreshold; }

      void SetTargetMaterial(TGeoMaterial *material) { fMaterial=material; }
      TGeoMaterial* TargetMaterial() { return fMaterial; }

      void SetSource(Source *source) { fSource=source; }
      Source* NeutrinoSource() { return fSource; }

      TF1* FXSxNe(UShort_t type=1, 
            Double_t nuclearRecoilEnergy=5*CLHEP::keV,
            Double_t maxNeutrinoEnergy=100*CLHEP::MeV);

      TF1* FNevt(UShort_t type=1,
            Double_t maxNuclearRecoilEnergy=25*CLHEP::keV,
            Double_t maxNeutrinoEnergy=100*CLHEP::MeV);

      ClassDef(Experiment,1);
};

#endif //EXPERIMENT_HH
