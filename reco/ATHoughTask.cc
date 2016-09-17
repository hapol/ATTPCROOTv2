#include "ATHoughSpace.hh"
#include "ATHoughSpaceLine.hh"
#include "ATHoughSpaceCircle.hh"


// FAIRROOT classes
#include "FairRootManager.h"
#include "FairRun.h"
#include "FairRuntimeDb.h"

#include <iostream>

#include "ATHoughTask.hh"

ClassImp(ATHoughTask);

ATHoughTask::ATHoughTask():fAtPadCoord(boost::extents[10240][3][2])
{
  fLogger = FairLogger::GetLogger();
  fPar = NULL;

  fIsPersistence = kFALSE;
  fIsLinear = kFALSE;
  fIsCircular = kFALSE;
  fIsPhiReco = kFALSE;
  fRadThreshold=0.0;
  fHoughThreshold = 0.0;
  fHoughDistance = 5.0;

  fEvent = NULL;
  fProtoevent = NULL;

  fInternalID = 0;

  fIsEnableMap = kFALSE;



}

ATHoughTask::~ATHoughTask()
{
}


void   ATHoughTask::SetPersistence(Bool_t value)           { fIsPersistence = value; }
void   ATHoughTask::SetThreshold(Double_t threshold)       { fThreshold = threshold; }
void   ATHoughTask::SetRadiusThreshold(Float_t value)      { fRadThreshold = value; }
void   ATHoughTask::SetLinearHough()                       { fIsLinear = kTRUE;fIsCircular = kFALSE;}
void   ATHoughTask::SetCircularHough()                     { fIsCircular = kTRUE;fIsLinear = kFALSE;}
void   ATHoughTask::SetPhiReco()                           { fIsPhiReco = kTRUE;}
void   ATHoughTask::SetHoughThreshold(Double_t value)      { fHoughThreshold = value;}
void   ATHoughTask::SetHoughDistance(Double_t value)       { fHoughDistance = value;}
void   ATHoughTask::SetEnableMap()                         { fIsEnableMap = kTRUE;}
void   ATHoughTask::SetMap(Char_t const *map)              { fMap = map; }

InitStatus
ATHoughTask::Init()
{


  if(fIsLinear) fHoughArray = new TClonesArray("ATHoughSpaceLine");
  else if(fIsCircular) fHoughArray = new TClonesArray("ATHoughSpaceCircle");
  else{

    fLogger -> Error(MESSAGE_ORIGIN, "-I- ATHoughTask : Hough Space Calculation NOT Set. Please choose a Hough Space Topology");
    return kERROR;

  }

  FairRootManager *ioMan = FairRootManager::Instance();
  if (ioMan == 0) {
    fLogger -> Error(MESSAGE_ORIGIN, "Cannot find RootManager!");
    return kERROR;
  }

  fEventHArray = (TClonesArray *) ioMan -> GetObject("ATEventH");
  if (fEventHArray == 0) {
    fLogger -> Error(MESSAGE_ORIGIN, "Cannot find ATEvent array!");
    return kERROR;
  }

  if(fIsPhiReco){ //Find the Array of ProtoEvents
      fProtoEventHArray = (TClonesArray *) ioMan -> GetObject("ATProtoEvent");
      if (fProtoEventHArray == 0) {
        fLogger -> Error(MESSAGE_ORIGIN, "Cannot find ATProtoEvent array! If SetPhiReco method is enabled, Phi Reconstruction is needed");
        return kERROR;
      }
  }

  // Pointer to the Pad Plane map for digitization during the MC
  if(fIsEnableMap){
  fAtMapPtr = new AtTpcMap();
  fAtMapPtr->GenerateATTPC();
  fPadPlane = fAtMapPtr->GetATTPCPlane();
  Bool_t MapIn = fAtMapPtr->ParseXMLMap(fMap);
  fLogger -> Info(MESSAGE_ORIGIN, "ATTPC Map enabled");
  if(!MapIn) std::cerr<<" -E- ATHoughTask - : Map was enabled but not found ! "<<std::endl;
  fAtPadCoord = fAtMapPtr->GetPadCoordArr();
  }



  ioMan -> Register("ATHough", "ATTPC", fHoughArray, fIsPersistence);




  return kSUCCESS;
}

void
ATHoughTask::SetParContainers()
{



  FairRun *run = FairRun::Instance();
  if (!run)
    fLogger -> Fatal(MESSAGE_ORIGIN, "No analysis run!");

  FairRuntimeDb *db = run -> GetRuntimeDb();
  if (!db)
    fLogger -> Fatal(MESSAGE_ORIGIN, "No runtime database!");

  fPar = (ATDigiPar *) db -> getContainer("ATDigiPar");
  if (!fPar)
    fLogger -> Fatal(MESSAGE_ORIGIN, "ATDigiPar not found!!");
}

void
ATHoughTask::Exec(Option_t *opt)
{
   fHoughArray -> Delete();



    if (fEventHArray -> GetEntriesFast() == 0)
     return;

     if(fIsPhiReco){
         if (fProtoEventHArray -> GetEntriesFast() == 0)
          return;
        }


     //ATEvent *Event = (ATEvent *) fEventHArray -> At(0);
     fEvent  = (ATEvent *) fEventHArray -> At(0);
     if(fIsPhiReco) fProtoevent = (ATProtoEvent *) fProtoEventHArray -> At(0);
     fInternalID++;
     if(fInternalID%100==0) std::cout << "  -I- ATHoughTask -  Event Number :  " << fEvent -> GetEventID()<<" Internal ID : "<<fInternalID<< std::endl;

    if(fIsLinear){ // TODO: Solve this dirty way with a dynamic cast and make global pointers

            ATHoughSpaceLine *HoughSpace = (ATHoughSpaceLine *) new ((*fHoughArray)[0]) ATHoughSpaceLine();
            HoughSpace->SetRadiusThreshold(fRadThreshold);
            HoughSpace->SetHoughDistance(fHoughDistance);
            if(fIsPhiReco) HoughSpace ->CalcHoughSpace(fProtoevent,kTRUE,kTRUE,kTRUE,kTRUE);
            else HoughSpace->CalcHoughSpace(fEvent);


    }
    else if(fIsCircular){
            ATHoughSpaceCircle *HoughSpace = (ATHoughSpaceCircle *) new ((*fHoughArray)[0]) ATHoughSpaceCircle();
            HoughSpace ->SetThreshold(fHoughThreshold);
            //if(fIsEnableMap) HoughSpace ->CalcHoughSpace(fEvent,fPadPlane);
            if(fIsEnableMap) HoughSpace ->CalcHoughSpace(fEvent,fPadPlane,fAtPadCoord);
            else HoughSpace ->CalcHoughSpace(fEvent,kTRUE,kTRUE,kTRUE);
    }



  //(ATHoughSpaceLine *) new ((*fHoughArray)[0]) ATHoughSpaceLine();
  //event -> SetEventID(event -> GetEventID());
   /* event -> SetEventID(rawEvent -> GetEventID());

  if (!(rawEvent -> IsGood()))
    event -> SetIsGood(kFALSE);
  else {
    fPSA -> Analyze(rawEvent, event);
    event -> SetIsGood(kTRUE);
  }*/
}

/*void
ATHoughTask::FinishEvent()
{

  if (fEventHArray -> GetEntriesFast() == 0) return;

  fEvent  = (ATEvent *) fEventHArray -> At(0);
  if (fEvent == NULL)
  {
    fLogger -> Info(MESSAGE_ORIGIN, "End of file. Terminating FairRun.");
    FairRootManager::Instance() -> SetFinishRun();
  }
}*/
