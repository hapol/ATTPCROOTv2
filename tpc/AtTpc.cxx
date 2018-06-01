/********************************************************************************
 *    Copyright (C) 2014 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    *
 *                                                                              *
 *              This software is distributed under the terms of the             *
 *         GNU Lesser General Public Licence version 3 (LGPL) version 3,        *
 *                  copied verbatim in the file "LICENSE"                       *
 ********************************************************************************/
#include "AtTpc.h"

#include "AtTpcPoint.h"
#include "AtTpcGeo.h"
#include "AtTpcGeoPar.h"
#include "ATVertexPropagator.h"

#include "FairVolume.h"
#include "FairGeoVolume.h"
#include "FairGeoNode.h"
#include "FairRootManager.h"
#include "FairGeoLoader.h"
#include "FairGeoInterface.h"
#include "FairRun.h"
#include "FairRuntimeDb.h"
#include "AtDetectorList.h"
#include "AtStack.h"

#include "TClonesArray.h"
#include "TVirtualMC.h"
#include "TGeoManager.h"
#include "TGeoBBox.h"
#include "TGeoCompositeShape.h"
#include "TGeoTube.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TParticle.h"
#include "TRandom.h"
#include "TRandom3.h"

#include <iostream>
using std::cout;
using std::endl;

AtTpc::AtTpc()
  : FairDetector("AtTpc", kTRUE, kAtTpc),
    fAtTpcPointCollection(new TClonesArray("AtTpcPoint")) {
}

AtTpc::AtTpc(const char* name, Bool_t active)
  : FairDetector(name, active, kAtTpc),
    fAtTpcPointCollection(new TClonesArray("AtTpcPoint")) {
}

AtTpc::~AtTpc() {
  if (fAtTpcPointCollection) {
    fAtTpcPointCollection->Delete();
    delete fAtTpcPointCollection;
  }
}

void AtTpc::Initialize() {
  FairDetector::Initialize();
  FairRuntimeDb* rtdb= FairRun::Instance()->GetRuntimeDb();
  AtTpcGeoPar* par=(AtTpcGeoPar*)(rtdb->getContainer("AtTpcGeoPar"));
}

Bool_t  AtTpc::ProcessHits(FairVolume* vol) {
  
  TLorentzVector pos;
  TLorentzVector mom;
  gMC->TrackPosition(pos);
  gMC->TrackMomentum(mom);


  Int_t theTrackStatus = GetTrackStatus(gMC->IsNewTrack(),
					gMC->IsTrackDisappeared(),
					gMC->IsTrackStop(),
					gMC->IsTrackAlive(),
					gMC->IsTrackEntering(),
					gMC->IsTrackExiting(),
					gMC->IsTrackInside(),
					gMC->IsTrackOut());
  
  
  Int_t parentTrackID = gMC->GetStack()->GetCurrentParentTrackNumber();
  TString particleName = gMC->GetStack()->GetCurrentTrack()->GetName();
  //_______________only care about primary particle and decayed particle
  if(gMC->CurrentVolName()=="drift_volume" && 
     (parentTrackID==-1 || (parentTrackID==0 && particleName!="e-"))) {
    Int_t size = fAtTpcPointCollection->GetEntriesFast();
    new((*fAtTpcPointCollection)[size]) 
      AtTpcPoint(gMC->GetStack()->GetCurrentTrackNumber(), //trackID
		   vol->getModId(), //check if getModId or CurrentVolOffID(1,modID)
		   pos.Vect(), //pos from gMC->TrackPosition(pos);
		   mom.Vect(), //mom from gMC->TrackMomentum(pos);
		   gMC->TrackTime(), //time in s
		   gMC->TrackLength(), //length
		   gMC->Edep(), //eloss
		   gMC->CurrentEvent(), //EventID
		   theTrackStatus, //trackStatus from GetTrackStatus(...)
		   vol->getCopyNo(), //detCopyID
		   vol->GetName() //volName
		   );
  }
  
  // Increment number of Points for this track
  AtStack* stack = (AtStack*) gMC->GetStack();
  stack->AddPoint(kAtTpc);
  
  return kTRUE;
}

void AtTpc::EndOfEvent(){ 
  fAtTpcPointCollection->Clear(); 
}



void AtTpc::Register() {

  /** This will create a branch in the output tree called
      AtTpcPoint, setting the last parameter to kFALSE means:
      this collection will not be written to the file, it will exist
      only during the simulation.
  */
  
  FairRootManager::Instance()->Register("AtTpcPoint", "AtTpc",
                                        fAtTpcPointCollection, kTRUE);
}


TClonesArray* AtTpc::GetCollection(Int_t iColl) const {
  if (iColl == 0) { return fAtTpcPointCollection; }
  else { return NULL; }
}

void AtTpc::Reset() {
  fAtTpcPointCollection->Clear();
}

void AtTpc::Print(Option_t* option) const {
  Int_t nHits = fAtTpcPointCollection->GetEntriesFast();
  LOG(INFO) << "ATTPC: " << nHits << " points registered in this event" << FairLogger::endl;
}


void AtTpc::ConstructGeometry() {
  TString fileName=GetGeometryFileName();
  if (fileName.EndsWith(".geo")) {
    LOG(INFO)<<"Constructing ATTPC geometry from ASCII file "<<fileName<<FairLogger::endl;
    //ConstructASCIIGeometry();
  } else if (fileName.EndsWith(".root")) {
    LOG(INFO)<<"Constructing ATTPC geometry from ROOT file "<<fileName<<FairLogger::endl;
    ConstructRootGeometry();
  } else {
    std::cout << "Geometry format not supported." << std::endl;
  }
}

Bool_t AtTpc::CheckIfSensitive(std::string name) {
  
  TString tsname = name;
  if (tsname.Contains("drift_volume") || tsname.Contains("window")) {
    LOG(INFO)<<" ATTPC geometry: Sensitive volume found: "<<tsname<<FairLogger::endl;
    return kTRUE;
  }
  return kFALSE;
}

std::pair<Int_t,Int_t> AtTpc::DecodePdG(Int_t PdG_Code) {
  Int_t A = PdG_Code/10%1000;
  Int_t Z = PdG_Code/10000%1000;
  
  std::pair<Int_t,Int_t> nucleus;
  
  if(PdG_Code==2212){
    nucleus.first  = 1;
    nucleus.second = 1;
  }else if(PdG_Code==2112)
    {
      nucleus.first  = 1;
      nucleus.second = 0;
    }else{
    nucleus.first  = A;
    nucleus.second = Z;
  }
  return nucleus;
}


Int_t AtTpc::GetTrackStatus(bool NewTrack,
			      bool TrackDisappeared,
			      bool TrackStop,
			      bool TrackAlive,
			      bool TrackEntering,
			      bool TrackExiting,
			      bool TrackInside,
			      bool TrackOut) {
  int trackstatus =
    1*TrackOut+10*TrackInside+
    100*TrackExiting+1000*TrackEntering+
    10000*TrackAlive + 100000*TrackStop +
    1000000*TrackDisappeared + 10000000*NewTrack;
  
  return trackstatus;
  
}

ClassImp(AtTpc)
