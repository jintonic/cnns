#ifndef SUPERNOVAEXPERIMENT_H
#define SUPERNOVAEXPERIMENT_H

#include <TNamed.h>
class TF1;
class TH1D;
class TH2D;

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

      TF1 *fFXSxNe[7]; // dXS(Ev) * Ne(Ev)
      TF1 *fFXSxN2[7]; // dXS(Ev) * N2(time,Ev)
      TF1 *fFNevtE[7]; // Nevt(Enr) in TF1 format
      TH2D *fHNevt2[7]; // Nevt(t, Enr)
      TH1D *fHNevtT[7]; // Nevt(t)

      Double_t XSxNe(Double_t *x, Double_t *parameter); // function of dXS * Ne
      Double_t XSxN2(Double_t *x, Double_t *parameter); // function of dXS * N2
      Double_t NevtE(Double_t *x, Double_t *parameter); // function of NevtE

   public:
      SupernovaExperiment(MAD::Material *material=0, NEUS::SupernovaModel *model=0);
      virtual ~SupernovaExperiment() { Clear(); } 

      void SetTargetMass(Double_t m) { fMass=m; }
      Double_t TargetMass() { return fMass; }

      void SetThreshold(Double_t threshold) { fThreshold=threshold; }
      Double_t Threshold() { return fThreshold; }

      void SetTargetMaterial(MAD::Material *material) { fMaterial=material; }
      MAD::Material* TargetMaterial() { return fMaterial; }

      void SetSupernovaModel(NEUS::SupernovaModel *model) { fModel = model; }
      NEUS::SupernovaModel* Model() { return fModel; }

      void SetDistance(Double_t distance) { fDistance=distance; }
      Double_t Distance() { return fDistance; }

      TF1* FXSxNe(UShort_t type, Double_t Enr);

      TF1* FNevtE(UShort_t type, Double_t maxEnr);

      TH1D* HXSxNe(UShort_t type, Double_t Enr);

      TH1D* HNevtE(UShort_t type, Double_t maxEnr);

      Double_t Nevt(Double_t maxEnr); // total number of events

      TF1* FXSxN2(UShort_t type, Double_t time, Double_t Enr);
      Double_t Nevt2(UShort_t type, Double_t time, Double_t Enr);
      TH2D* HNevt2(UShort_t type); // Nevt(t, Enr)
      TH1D* HNevtT(UShort_t type); // Nevt(t)

      /**
       * Delete internal objects.
       */
      void Clear(Option_t *option="");

      ClassDef(SupernovaExperiment,1);
};

#endif
