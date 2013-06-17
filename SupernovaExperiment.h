#ifndef SUPERNOVAEXPERIMENT_H
#define SUPERNOVAEXPERIMENT_H

#include <TNamed.h>
class TF1;
class TH1D;
class TH2D;

class Detector;

namespace NEUS { class SupernovaModel; }

class SupernovaExperiment : public TNamed
{
   protected:
      Detector* fDetector;
      NEUS::SupernovaModel *fModel; // supernova model
      Double_t fDistance; // distance between detector and Supernova

      TF1 *fFXSxN2[7]; // dXS(Ev) * N2(time,Ev)
      TF1 *fFXSxNe[7]; // dXS(Ev) * Ne(Ev)
      TH2D *fHNevt2[7]; // Nevt(t, Enr)
      TH1D *fHNevtT[7]; // Nevt(t)
      TH1D *fHNevtE[7]; // Nevt(Enr)

      Double_t XSxNe(Double_t *x, Double_t *parameter); // function of dXS * Ne
      Double_t XSxN2(Double_t *x, Double_t *parameter); // function of dXS * N2

   public:
      SupernovaExperiment(Detector *detector=0, NEUS::SupernovaModel *model=0);
      virtual ~SupernovaExperiment() { Clear(); } 

      void SeteDetector(Detector *detector) { Clear(); fDetector = detector; }

      void SetSupernovaModel(NEUS::SupernovaModel *model)
      { Clear(); fModel = model; }
      NEUS::SupernovaModel* Model() { return fModel; }

      void SetDistance(Double_t distance) { fDistance=distance; }
      Double_t Distance() { return fDistance; }

      TF1* FXSxNe(UShort_t type, Double_t Enr);
      TH1D* HXSxNe(UShort_t type, Double_t Enr);
      Double_t NevtE(UShort_t type, Double_t Enr);
      TH1D* HNevtE(UShort_t type); // Nevt(Enr)

      Double_t Nevt(); // total number of events

      TF1* FXSxN2(UShort_t type, Double_t time, Double_t Enr);
      Double_t Nevt2(UShort_t type, Double_t time, Double_t Enr);
      TH2D* HNevt2(UShort_t type); // Nevt(t, Enr)
      TH1D* HNevtT(UShort_t type, Bool_t detectableOnly=kFALSE); // Nevt(t)

      /**
       * Delete internal objects.
       */
      void Clear(Option_t *option="");

      ClassDef(SupernovaExperiment,1);
};

#endif
