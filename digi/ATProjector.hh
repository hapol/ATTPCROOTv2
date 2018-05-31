// -----------------------------------------------------------------------------
// -----                  ATProjector header file                     -----
// -----                  Created 27/03/18  by H.Alvarez                   -----
// -----                Last modification 23/05/18                         -----
// -----------------------------------------------------------------------------

/**  ATProjector.hh 
 * Projects the electrons derived from the energy loss in the gas into
 * a virtual pad plane of variable size pads, including time information
 **/
#ifndef ATPROJECTOR_H
#define ATPROJECTOR_H

#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"

#include "FairTask.h"
#include "FairMCPoint.h"

#include "TClonesArray.h"
#include "AtTpcPoint.h"
#include "ATDigiPar.hh"
#include "TH2Poly.h"
#include "AtTpcMap.h"
#include "ATRawEvent.hh"
#include "ATGas.hh"

/**
 * GTPC point projector task
 * @author Héctor Alvarez Pol
 *
 * For each event, get the AtTpcPoints and determine the projection on the pad plane
 *   Input:  Branch AtTpcPoints = TClonesArray("AtTpcPoint")
 *   Output: Branch fRawEventArray = TClonesArray("fRawEventArray")
 */

class ATProjector : public FairTask {
 public:

  /** Default constructor **/
  ATProjector();

  /** Destructor **/
  ~ATProjector();

  /** Virtual method Exec **/
  void Exec(Option_t*);

  /** Set parameters -- To be removed when parameter containers are ready **/
  void SetDriftParameters(Double_t ion, Double_t driftv,
			  Double_t tDiff, Double_t lDiff,
			  Double_t fanoFactor);

 protected:

  /** Virtual method Init **/
  InitStatus Init();

  /** Virtual method Finish **/
  void Finish();

  /** Virtual method SetParContainers **/
  void SetParContainers();

 private:

  TClonesArray* fAtTpcPoints;          //!< Incoming points (NEW _2)
  TClonesArray* fRawEventArray;        //!< Raw Event array(only one)
  ATRawEvent* fRawEvent;               //!< Raw Event Object

  ATGas*     fGas;                     //!< Gas parameter container.
  ATDigiPar* fPar;                     //!< Base parameter container.

  Double_t fGain;                      //!< Gain.

  Double_t fEIonization;   //!< Effective ionization energy of gas [eV]
  Double_t fDriftVelocity; //!< Drift velocity in gas [cm/ns]
  Double_t fTransDiff;     //!< Transversal diffusion coefficient [cm^2/ns]
  Double_t fLongDiff;      //!< Longitudinal diffusion coefficient [cm^2/ns]
  Double_t fFanoFactor;    //!< Fano factor to calculate electron number fluctuations

  TH2Poly *fPadPlane;                  //!< pad plane
  AtTpcMap *fMap;                      //!<ATTPC map

  ClassDef(ATProjector, 1)
};

#endif // ATPROJECTOR_H
