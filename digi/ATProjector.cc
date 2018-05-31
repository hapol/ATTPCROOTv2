#include "ATProjector.hh"

#include "TClonesArray.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"
#include "ATPad.hh"

// STL class headers
#include <cmath>
#include <iostream>
#include <iomanip>

#include "TRandom.h"
#include "TMath.h"
#include "TF1.h"

ATProjector::ATProjector() : FairTask("ATProjector") {
  fEIonization = 15.e-9;  // [GeV]
  fDriftVelocity = 0.0005; // [cm/ns]
  fTransDiff = 0.00000001; // [cm^2/ns]
  fLongDiff  = 0.0000001;  // [cm^2/ns]
  fFanoFactor = 2; 
}

ATProjector::~ATProjector(){
  fAtTpcPoints->Clear();
  fRawEventArray->Clear();
}

void ATProjector::SetParContainers() {
  FairRunAna* run = FairRunAna::Instance();
  if (!run){
    fLogger -> Error(MESSAGE_ORIGIN, "ATProjector::SetParContainers: No analysis run");
    return;
  }
  FairRuntimeDb* rtdb = run->GetRuntimeDb();
  if (!rtdb) {
    fLogger -> Error(MESSAGE_ORIGIN, "ATProjector::SetParContainers: No runtime database");
    return;
  }

  fPar = (ATDigiPar*) rtdb->getContainer("ATDigiPar");

}

InitStatus ATProjector::Init() {
  FairRootManager* ioman = FairRootManager::Instance();
  if (!ioman) {
    fLogger -> Error(MESSAGE_ORIGIN, "ATProjector::Init: No FairRootManager");
    return kFATAL;
  }
  // Input: TClonesArray of AtTpcPoints
  if ((TClonesArray*)ioman->GetObject("AtTpcPoint") == nullptr) {
    fLogger -> Error(MESSAGE_ORIGIN, "ATProjector::Init No AtTpcPoint!");
    return kFATAL;
  }
  fAtTpcPoints = (TClonesArray*)ioman->GetObject("AtTpcPoint");
  
  //Output:
  fRawEventArray  = new TClonesArray("ATRawEvent", 100);        //!< Raw Event array(only one)
  ioman -> Register("ATRawEvent", "cbmsim", fRawEventArray, kTRUE);
  

  ////PARS
  fEIonization = fPar->GetEIonize();  // [GeV]
  fDriftVelocity = fPar->GetDriftVelocity()/1000; // [cm/us]
  fTransDiff = fPar->GetCoefDiffusionTrans()*sqrt(10.); // [cm^2/ns] ???
  fLongDiff  = fPar->GetCoefDiffusionLong()*sqrt(10.);  // [cm^2/ns] ???
  fFanoFactor = 2; 

  fGain = fPar->GetGain();
  std::cout<<"Gain: "<<fGain<<std::endl;

  // ***************Create ATTPC Pad Plane***************************
  TString scriptfile = "Lookup20150611.xml";
  TString dir = getenv("VMCWORKDIR");
  TString scriptdir = dir + "/scripts/"+ scriptfile;
  
  fMap = new AtTpcMap();
  fMap->GenerateATTPC();
  Bool_t MapIn = fMap->ParseXMLMap(scriptdir);
  fPadPlane = fMap->GetATTPCPlane();
  
  fRawEvent = NULL;

  return kSUCCESS;
}

void ATProjector::SetDriftParameters(Double_t ion, Double_t driftv, 
				     Double_t tDiff, Double_t lDiff, 
				     Double_t fanoFactor) {
  fEIonization = ion;       // [eV]    15/1000000000; [GeV]
  fDriftVelocity = driftv;  // [mm/ns]
  fTransDiff = tDiff;       // [mm^(-1/2)] ~0.010
  fLongDiff = lDiff;        // [mm^(-1/2)] ~0.025
  fFanoFactor = fanoFactor; //
}

struct vPad{
  Double_t RawADC[512];
  Int_t padnumb;
};

void ATProjector::Exec(Option_t*) {
  fLogger->Debug(MESSAGE_ORIGIN,"Exec of ATProjector");
  Int_t nPoints = fAtTpcPoints->GetEntries();
  if(nPoints<2){
    fLogger->Warning(MESSAGE_ORIGIN, "Not enough hits for digitization! (<2)");
    return;
  }  
  
  fRawEventArray -> Delete();
  fRawEvent = NULL;
  fRawEvent = (ATRawEvent*)fRawEventArray->ConstructedAt(0);

  AtTpcPoint* aPoint; 
  //ATProjPoint* aProjPoint; 
  Int_t presentTrackID = -10; //control of the point trackID 
  Double_t xPre, yPre, zPre; Double_t xPost, yPost, zPost; Double_t projX, projY, projTime;
  Double_t energyDep = 0.; Double_t timeBeforeDrift = 0.;
  Bool_t readyToProject = kFALSE; Bool_t virtualPadFound = kFALSE;
  Int_t electrons = 0; Int_t flucElectrons = 0; Int_t generatedElectrons = 0;
  Double_t stepX, stepY, stepZ;
  Double_t driftDistance, driftTime;
  Double_t sigmaLongAtPadPlane;  Double_t sigmaTransvAtPadPlane; 
  Int_t evtID;
  
  Double_t tau = 1; //shaping time (us)
  Double_t samplingtime = 60;
  Double_t samplingrate = 0.080; //us
  Double_t timeBucket[512] = {0};
  Int_t counter = 0;
  Double_t output = 0;
  std::vector<vPad> padarray;
  Double_t output_sum = 0;
  Int_t tbcounter = 0;
  Double_t c = 10;
  Int_t  vsize;
  Int_t cc = 0;
  Double_t pBin, g, xElectron, yElectron, eTime, clusterNum, eventID;
  Int_t padNumber;
  TVector3 coord;
  std::vector<Float_t> PadCenterCoord;
  //Polya distribution of gain
  TF1 *gain =  new TF1("gain", "4*(x/[0])*pow(2.718, -2*(x/[0]))", 80, 120);
  gain->SetParameter(0, fGain);  
  // ***************Create Time Buckets*******************************
  for(Double_t d = 40; d<81; d+=samplingrate) {
    timeBucket[counter] = d;
    counter++;
  }
  
  for(Int_t i=0; i<nPoints; i++) {
    aPoint = (AtTpcPoint*) fAtTpcPoints->At(i);
    evtID=aPoint->GetEventID();
    if(aPoint->GetTrackStatus() == 11000 || 
       aPoint->GetTrackStatus() == 10010010  || 
       aPoint->GetTrackStatus() == 10010000) { 
      //entering the gas volume or new track inside the gas (is 10010010 or 10010000??)
      presentTrackID = aPoint->GetTrackID();
      xPre = aPoint->GetX(); yPre = aPoint->GetY(); zPre = aPoint->GetZ();
      readyToProject = kTRUE;
      continue; //no energy deposited in this point, just taking in entrance coordinates
    }
    else { //any other case
      if(presentTrackID != aPoint->GetTrackID()) {
	fLogger->Error(MESSAGE_ORIGIN, "ATProjector::Exec: Problem 2 in point logic");
	break;
      }
      if(readyToProject!=kTRUE) {
	fLogger->Error(MESSAGE_ORIGIN, "ATProjector::Exec: Problem 3 in point logic");
	break;
      }
      if (aPoint->GetTrackStatus() == 10100 || 
	  aPoint->GetTrackStatus() == 1000000 ) { //exiting the gas volume or dissappearing
	readyToProject = kFALSE;
      }	
      
      xPost = aPoint->GetX(); yPost = aPoint->GetY(); zPost = aPoint->GetZ();   
      energyDep = aPoint->GetEnergyLoss();
      timeBeforeDrift = aPoint->GetTime(); 
    }
    
    electrons = energyDep/fEIonization;
    //electron number fluctuates as the square root of the 
    //Fano factor times the number of electrons
    flucElectrons = pow(fFanoFactor*electrons, 0.5); 
    generatedElectrons = gRandom->Gaus(electrons, flucElectrons); //generated electrons
    
    //step in each direction for an homogeneous electron creation position along the track
    stepX = (xPost-xPre) / generatedElectrons;
    stepY = (yPost-yPre) / generatedElectrons;
    stepZ = (zPost-zPre) / generatedElectrons;
    
    //taken a mean driftDistance for the calculation of the sigmaLong and sigmaTrans
    //improve (make the calculation individual for electron) if needed, but probably slower
    Double_t zApprox = (zPost+zPre)/2;
    driftDistance = abs(1000-zApprox*10); //in mm!!
    sigmaLongAtPadPlane = sqrt(driftDistance*2*fLongDiff/fDriftVelocity);
    sigmaTransvAtPadPlane = sqrt(driftDistance*2*fTransDiff/fDriftVelocity);
    
    for(Int_t ele=1;ele<=generatedElectrons;ele++){
      driftTime = (1000-(zPre+stepZ*ele))/fDriftVelocity;
      projX =  10*gRandom->Gaus(xPre+stepX*ele,sigmaTransvAtPadPlane); //in mm!!
      projY =  10*gRandom->Gaus(yPre+stepY*ele,sigmaTransvAtPadPlane); //in mm!!
      projTime = gRandom->Gaus(driftTime+timeBeforeDrift,sigmaLongAtPadPlane/fDriftVelocity);
      
      xElectron                     = projX; //mm
      yElectron                     = projY; //mm
      eTime                         = projTime; //us
      counter                       = 0;
      pBin                          = fPadPlane->Fill(xElectron,yElectron,c);
      padNumber                     = pBin-1;
      Double_t pointmem[1000][3]    = {0};
      Double_t digital[512]         = {0};
      
      //*******Create new element in padarray if there's a new pad******//
      TString check = kTRUE;
      vsize  = padarray.size();
      for(Int_t r = 0; r<vsize; r++){
	if(padNumber == padarray[r].padnumb) check = kFALSE;
      }
      if(check == kTRUE){
	padarray.push_back(vPad());
	padarray[vsize].padnumb = padNumber;
      }
      
      // *********Pulse Generation for each electron************
      for(Double_t j = eTime; j<eTime+10; j+=samplingrate/5){
	output = pow(2.718,-3*((j-eTime)/tau))*sin((j-eTime)/tau)*pow((j-eTime)/tau,3);
	pointmem[counter][0]  = j;
	pointmem[counter][1]  = output;
	counter++;
	
	// Once a point is assigned a height in time, it assigns time to a time bucket
	for(Int_t k = 0; k<512; k++){//go through all time buckets
	  if(j>=timeBucket[k] && j<timeBucket[k+1]){//if point on pulse is in this time bucket, assign it that time bucket
	    pointmem[counter][2]  = k;
	    break;
	  }//end if for time buckets
	}//end assigning pulse to time buckets
      }//end plotting pulse function

      //*********Once pulse is generated, it adds points to the running average********
      Int_t pTimebucket  = pointmem[0][2];
      Int_t A            = 0;
      Int_t nOPoints     = 0;
      Double_t acum[512] = {0};
      while (A<1000){
	if(pTimebucket == pointmem[A][2]){
	  acum[pTimebucket]+= pointmem[A][1];
	  nOPoints++;
	  A++;
	}
	else{
	  digital[pTimebucket] = acum[pTimebucket]/nOPoints++;
	  pTimebucket          = pointmem[A][2];
	  nOPoints             = 0;
	}
      }
      
      // ********Adds pulse to output array**************
      vsize  = padarray.size();
      for(Int_t y = 0; y<vsize; y++){
	if(padarray[y].padnumb == padNumber){
	  for(Int_t del = 0; del<512; del++){//go through every time bucket
	    if(digital[del] != 0)  padarray[y].RawADC[del] += digital[del];
	  }
	  break;
	}
      }
    }// end through all electrons
  }//end through all Points    
  
  // ***************Set Pad and add to event**************
  vsize = padarray.size();
  Int_t thepad;
  for(Int_t q = 0; q<vsize; q++){
    g = gain->GetRandom();
    ATPad *pad = new ATPad();
    thepad = padarray[q].padnumb;
    if(thepad<10240 && thepad>0){
      pad->SetPad(thepad);
      PadCenterCoord = fMap->CalcPadCenter(thepad);
      pad->SetValidPad(kTRUE);
      pad->SetPadXCoord(PadCenterCoord[0]);
      pad->SetPadYCoord(PadCenterCoord[1]);
      pad->SetPedestalSubtracted(kTRUE);
      for(Int_t p = 0; p<512; p++){
	pad->SetADC(p, padarray[q].RawADC[p]*g);
      }
      fRawEvent->SetPad(pad);
      fRawEvent->SetEventID(aPoint->GetEventID());
    }
  }
  padarray.clear();
  //return;
}


void ATProjector::Finish() {
  
}

ClassImp(ATProjector)
