/********************************************************************************
 *    Copyright (C) 2014 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    *
 *                                                                              *
 *              This software is distributed under the terms of the             *
 *         GNU Lesser General Public Licence version 3 (LGPL) version 3,        *
 *                  copied verbatim in the file "LICENSE"                       *
 ********************************************************************************/
#ifndef ATTPC_H
#define ATTPC_H

#include "FairDetector.h"
#include "TVector3.h"
#include "TLorentzVector.h"

class AtTpcPoint;
class FairVolume;
class TClonesArray;

class AtTpc: public FairDetector {

 public:

  /**      Name :  Detector Name
   *       Active: kTRUE for active detectors (ProcessHits() will be called)
   *               kFALSE for inactive detectors
   */
  AtTpc(const char* Name, Bool_t Active);

  /**      default constructor    */
  AtTpc();

  /**       destructor     */
  virtual ~AtTpc();

  /**      Initialization of the detector is done here    */
  virtual void   Initialize();

  /**       this method is called for each step during simulation
   *       (see FairMCApplication::Stepping())
   */
  virtual Bool_t ProcessHits( FairVolume* v=0);

  /**       Registers the produced collections in FAIRRootManager.     */
  virtual void   Register();

  /** Gets the produced collections */
  virtual TClonesArray* GetCollection(Int_t iColl) const ;

  /**      has to be called after each event to reset the containers      */
  virtual void   Reset();

  /** Virtual method Print
   **
   ** Screen output of hit collection.
   **/
  virtual void Print(Option_t* option = "") const;

  /** Virtual method EndOfEvent
   **
   ** If verbosity level is set, print hit collection at the
   ** end of the event and resets it afterwards.
   **/
  virtual void EndOfEvent();

  /**      Create the detector geometry        */
  void ConstructGeometry();

  Bool_t CheckIfSensitive(std::string name);

  std::pair<Int_t,Int_t> DecodePdG(Int_t PdG_Code);

  Int_t GetTrackStatus(bool NewTrack,
		       bool TrackDisappeared,
		       bool TrackStop,
		       bool TrackAlive,
		       bool TrackEntering,
		       bool TrackExiting,
		       bool TrackInside,
		       bool TrackOut);

 private:

    /** container for data points */
    TClonesArray*  fAtTpcPointCollection;  //!

    AtTpc(const AtTpc&);
    AtTpc& operator=(const AtTpc&);

    ClassDef(AtTpc,1)
};

#endif //NEWDETECTOR_H
