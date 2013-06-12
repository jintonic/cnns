#ifndef SUPERNOVAEXPERIMENT_H
#define SUPERNOVAEXPERIMENT_H

#include <TNamed.h>
class TF1;
class TH1D;
class TH2D;

#include <UNIC/Units.h>

namespace MAD { class Material; }
namespace NEUS { class SupernovaModel; }

class SupernovaExperiment : public TNamed
{
   protected:
      MAD::Material *fMaterial; // target material
      NEUS::SupernovaModel *fModel; // supernova model

      Double_t fMass;      // mass of target material
      Double_t fDistance; // distance between detector and Supernova
      Double_t fThreshold; // energy threshold of detector

      TF1 *fXSxNe[7]; // dXS(Ev) * Ne(Ev)
      TF1 *fXSxN2[7]; // dXS(Ev) * N2(time,Ev)
      TF1 *fNevt[7]; // Nevt in TF1 format
      TH2D *fHN2[7]; // Nevt(t, Enr)
      TH1D *fHNt[7]; // Nevt(t)

      Double_t XSxNe(Double_t *x, Double_t *parameter); // function of dXS * Ne
      Double_t XSxN2(Double_t *x, Double_t *parameter); // function of dXS * N2
      Double_t FuncN(Double_t *x, Double_t *parameter); // function of Nevt

   public:
      SupernovaExperiment(MAD::Material *material=0, NEUS::SupernovaModel *model=0);
      virtual ~SupernovaExperiment() { Clear(); } 

      // in case of reactor, it gives Nevt/keVnr/second
      // in case of supernova, it gives Nevt/keVnr in 20 second
      Double_t Nevt(UShort_t type,
            Double_t nEr, 
            Double_t maxEv=82.5*UNIC::MeV);

      void SetTargetMass(Double_t m) { fMass=m; }
      Double_t TargetMass() { return fMass; }

      void SetThreshold(Double_t threshold) { fThreshold=threshold; }
      Double_t Threshold() { return fThreshold; }

      void SetTargetMaterial(MAD::Material *material) { fMaterial=material; }
      MAD::Material* TargetMaterial() { return fMaterial; }

      void SetModel(NEUS::SupernovaModel *model) { fModel = model; }
      NEUS::SupernovaModel* Model() { return fModel; }

      void SetDistance(Double_t distance) { fDistance=distance; }
      Double_t Distance() { return fDistance; }

      TF1* FXSxNe(UShort_t type=1, 
            Double_t nEr=5*UNIC::keV,
            Double_t minEv= 2.5*UNIC::MeV,
            Double_t maxEv=82.5*UNIC::MeV);

      TF1* FNevt(UShort_t type=1,
            Double_t maxEr=50*UNIC::keV,
            Double_t maxEv=82.5*UNIC::MeV);

      TH1D* HXSxNe(UShort_t type=1,
            Double_t nEr=5*UNIC::keV,
            Double_t minEv= 2.5*UNIC::MeV,
            Double_t maxEv=82.5*UNIC::MeV);

      TH1D* HNevt(UShort_t type=1,
            Double_t maxEr=50*UNIC::keV,
            Double_t maxEv=82.5*UNIC::MeV);

      Double_t N2(UShort_t type, Double_t time, Double_t Enr);
      TH2D* HN2(UShort_t type); // Nevt(t, Enr)
      TH1D* HNt(UShort_t type); // Nevt(t)
      TF1* FXSxN2(UShort_t type=1, Double_t time=0, Double_t Enr=5*UNIC::keV,
            Double_t minEv= 2.5*UNIC::MeV, Double_t maxEv=82.5*UNIC::MeV);


      /**
       * Delete internal objects.
       */
      void Clear(Option_t *option="");

      ClassDef(SupernovaExperiment,1);
};

#endif
