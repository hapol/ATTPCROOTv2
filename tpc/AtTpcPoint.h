/********************************************************************************
 *    Copyright (C) 2014 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    *
 *                                                                              *
 *              This software is distributed under the terms of the             *
 *         GNU Lesser General Public Licence version 3 (LGPL) version 3,        *
 *                  copied verbatim in the file "LICENSE"                       *
 ********************************************************************************/
#ifndef NEWDETECTORPOINT_H
#define NEWDETECTORPOINT_H 1


#include "FairMCPoint.h"
#include "TObject.h"
#include "TVector3.h"


class AtTpcPoint : public FairMCPoint {
  
 public:
  
  /** Default constructor **/
  AtTpcPoint();
  
  
  /** Constructor with arguments
   *@param trackID  Index of MCTrack
   *@param detID    Detector ID
   *@param pos      Ccoordinates at entrance to active volume [cm]
   *@param mom      Momentum of track at entrance [GeV]
   *@param tof      Time since event start [ns]
   *@param length   Track length since creation [cm]
   *@param eLoss    Energy deposit [GeV]
   *@param eventId            MC event id
   *@param trackStatus        Status of the track
   *@param detCopyID          Detector Copy ID 
   *@param volName            Volume name for a given volume identifier id
   **/
  AtTpcPoint(Int_t trackID, Int_t detID, TVector3 pos,
	     TVector3 mom, Double_t tof, Double_t length, 
	     Double_t eLoss, const UInt_t EventId,
	     Int_t trackStatus, Int_t detCopyID, TString volName);
  
  /** Destructor **/
  virtual ~AtTpcPoint();
  
  /** Accessors **/
  Int_t    GetDetCopyID() const { return fDetCopyID; } // added by Marc
  TString  GetVolName()   const { return fVolName; }
  Int_t GetTrackStatus() const { return fTrackStatus; }

  
  /** Modifiers **/
  void SetDetCopyID(Int_t id)       { fDetCopyID = id; }; // added by Marc
  void SetVolName(TString name)  { fVolName = name; }    
  void SetTrackStatus(Int_t code) { fTrackStatus = code; }          
    
  /** Output to screen **/
  virtual void Print(const Option_t* opt) const;
  
 private:
  /** Copy constructor **/
  AtTpcPoint(const AtTpcPoint& point);
  AtTpcPoint operator=(const AtTpcPoint& point);
  
 protected:
  
  Int_t fDetCopyID;
  TString fVolName;
  Int_t fTrackStatus;           ///< Status of the track   
  
  ClassDef(AtTpcPoint,1)
    
};

#endif
