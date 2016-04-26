#include "FairRootManager.h"

#include "ATEventDrawTaskProto.hh"


#include "TEveManager.h"
#include "TEveGeoShape.h"
#include "TEveTrans.h"
#include "TGeoSphere.h"
#include "TEveTrans.h"
#include "TPaletteAxis.h"
#include "TStyle.h"
#include "TRandom.h"
#include "TColor.h"

#include "AtTpcMap.h"
#include "AtTpcProtoMap.h"
#include "TH2Poly.h"
#include "TF1.h"

#ifndef __CINT__ // Boost
#include <boost/multi_array.hpp>
#endif //__CINT__

#include <iostream>

#define cRED "\033[1;31m"
#define cYELLOW "\033[1;33m"
#define cNORMAL "\033[0m"
#define cGREEN "\033[1;32m"


using namespace std;

ClassImp(ATEventDrawTaskProto);


ATEventDrawTaskProto::ATEventDrawTaskProto():
fHitColor(kPink),
fHitSize(1),
fHitStyle(kFullDotMedium),
fHitSet(0),
fhitBoxSet(0)
{
  Char_t padhistname[256];
  fMultiHit=10;

  for(Int_t i=0;i<300;i++){ // TODO: Full-scale must be accomodated
      sprintf(padhistname,"pad_%d",i);
        fPadAll[i] = new TH1I(padhistname,padhistname,512,0,511);
  }

  Char_t phihistname[256];

          for(Int_t i=0;i<5;i++){
           sprintf(phihistname,"PhiDistr_%d",i);
           fPhiDistr[i] = new TH1D(phihistname,phihistname,180.0,-180.0,180.0);
           if(i==0) fPhiDistr[i]->SetLineColor(kRed);
           else if(i==1) fPhiDistr[i]->SetLineColor(kBlue);
           else if(i==2) fPhiDistr[i]->SetLineColor(kGreen);
           else if(i==3) fPhiDistr[i]->SetLineColor(kCyan);
                 else if(i==4) fPhiDistr[i]->SetLineColor(kMagenta);
                 fPhiDistr[i]->SetLineWidth(2);
                 fPhiDistr[i]->GetYaxis()->SetRangeUser(0., 20.);
          }

          fIsCircularHough=kFALSE;
          fIsLinearHough=kTRUE;
          fIsRawData=kFALSE;
          fHoughLinearFit =new TF1("HoughLinearFit"," (  (-TMath::Cos([0])/TMath::Sin([0]))*x ) + [1]/TMath::Sin([0])",0,500);

          f3DHitStyle=1;


}

ATEventDrawTaskProto::~ATEventDrawTaskProto()
{

    delete fHoughLinearFit;

}

InitStatus
ATEventDrawTaskProto::Init()
{


  std::cout<<" =====  ATEventDrawTaskProto::Init ====="<<std::endl;

  gROOT->Reset();
  FairRootManager* ioMan = FairRootManager::Instance();
  fEventManager = ATEventManagerProto::Instance();
  fDetmap  =  new AtTpcProtoMap();
  fDetmap -> SetProtoMap(fMap.Data());
  fDetmap -> SetGeoFile("proto_geo_hires.root");
  fDetmap -> SetName("fMap");
  gROOT->GetListOfSpecials()->Add(fDetmap);

  fHitArray = (TClonesArray*) ioMan->GetObject("ATEventH"); // TODO: Why this confusing name? It should be fEventArray
  if(fHitArray) LOG(INFO)<<"Hit Array Found."<<FairLogger::endl;

  fRawEventArray = (TClonesArray*) ioMan->GetObject("ATRawEvent");
  if(fRawEventArray){
       LOG(INFO)<<"Raw Event Array  Found."<<FairLogger::endl;
       fIsRawData=kTRUE;
  }

  fHoughSpaceArray =  (TClonesArray*) ioMan->GetObject("ATHough");
  if(fHoughSpaceArray) LOG(INFO)<<"Hough Array Found."<<FairLogger::endl;


  fProtoEventArray =  (TClonesArray*) ioMan->GetObject("ATProtoEvent");
  if(fProtoEventArray) LOG(INFO)<<"Prototype Event Array Found."<<FairLogger::endl;

  //Drawing histograms

  gStyle -> SetPalette(55);
  fCvsPadWave = fEventManager->GetCvsPadWave();
  fCvsPadWave->SetName("fCvsPadWave");
  gROOT->GetListOfSpecials()->Add(fCvsPadWave);
  DrawPadWave();
  fCvsPadPlane = fEventManager->GetCvsPadPlane();// There is a problem if the pad plane is drawn first
  fCvsPadPlane -> ToggleEventStatus();
  fCvsPadPlane->AddExec("ex","ATEventDrawTaskProto::SelectPad(\"fRawEvent\")");
  DrawPadPlane();

}

void
ATEventDrawTaskProto::Exec(Option_t* option)
{
    Reset();
  //ResetPadAll();
  //ResetPhiDistr();

    if(fHitArray) DrawHitPoints();

    gEve -> Redraw3D(kFALSE);

    UpdateCvsPadWave();
    UpdateCvsPadPlane();

}

void
ATEventDrawTaskProto::Reset()
{

  if(fHitSet) {
    fHitSet->Reset();
    gEve->RemoveElement(fHitSet, fEventManager);

  }

   if(fhitBoxSet) {
    fhitBoxSet->Reset();
    gEve->RemoveElement(fhitBoxSet, fEventManager);

  }

  if(fPadPlane!=NULL)
    fPadPlane->Reset(0);


}

void
ATEventDrawTaskProto::DrawHitPoints()
{

  Float_t *MeshArray;
  //fMesh->Reset(0);
  //f3DHist->Reset(0);
  //TRandom r(0);


  std::ofstream dumpEvent;
  dumpEvent.open ("event.dat");

  std::vector<Double_t> fPosXMin;
  std::vector<Double_t> fPosYMin;
  std::vector<Double_t> fPosZMin;


  //fQEventHist_H->Reset(0);
  ATEvent* event = (ATEvent*) fHitArray->At(0); // TODO: Why this confusing name? It should be fEventArray
  Double_t Qevent=event->GetEventCharge();
  Double_t RhoVariance=event->GetRhoVariance();
  MeshArray = event->GetMesh();
  Int_t eventID=event->GetEventID();
  TString TSevt =" Event ID : ";
  TString TSpad =" Pad ID : ";
  dumpEvent<<TSevt<<eventID<<std::endl;

  //if(fEventManager->GetEraseQEvent()){
	//fQEventHist->Reset();
  //      fRhoVariance->Reset();
  //}

  //fQEventHist->Fill(Qevent);
  //fQEventHist_H->Fill(Qevent);
  //fRhoVariance->Fill(RhoVariance);

    for(Int_t i=0;i<512;i++){

		//fMesh->SetBinContent(i,MeshArray[i]);

	}


  if(fIsRawData){
  fRawevent = (ATRawEvent*) fRawEventArray->At(0);
  fRawevent->SetName("fRawEvent");
  gROOT->GetListOfSpecials()->Add(fRawevent);
  }


  Int_t nHits = event->GetNumHits();
  fHitSet = new TEvePointSet("Hit",nHits, TEvePointSelectorConsumer::kTVT_XYZ);
  fHitSet->SetOwnIds(kTRUE);
  fHitSet->SetMarkerColor(fHitColor);
  fHitSet->SetMarkerSize(fHitSize);
  fHitSet->SetMarkerStyle(fHitStyle);
  std::cout<<cYELLOW<<" Number of hits : "<<nHits<<cNORMAL<<std::endl;


  //////////////////////////////////////////////

  fhitBoxSet = new TEveBoxSet("hitBox");
  fhitBoxSet->Reset(TEveBoxSet::kBT_AABox, kTRUE, 64);

  for(Int_t iHit=0; iHit<nHits; iHit++)
  {

    ATHit hit = event->GetHitArray()->at(iHit);
    Int_t PadNumHit = hit.GetHitPadNum();
    Int_t PadMultHit = event->GetHitPadMult(PadNumHit);
    Double_t BaseCorr = hit.GetBaseCorr();
    Int_t Atbin = -1;



    //if(hit.GetCharge()<fThreshold) continue;
    //if(PadMultHit>fMultiHit) continue;
    TVector3 position = hit.GetPosition();
    TVector3 positioncorr = hit.GetPositionCorr();


    fHitSet->SetMarkerColor(fHitColor);
    fHitSet->SetNextPoint(position.X()/10.,position.Y()/10.,position.Z()/10.); // Convert into cm
    fHitSet->SetPointId(new TNamed(Form("Hit %d",iHit),""));
    Atbin = fPadPlane->Fill(position.X(), position.Y(), hit.GetCharge());

    Bool_t fValidPad;


    if(fIsRawData){
    ATPad *RawPad = fRawevent->GetPad(PadNumHit,fValidPad);
    Double_t *adc = RawPad->GetADC();
        for(Int_t i=0;i<512;i++){

            //f3DThreshold = fEventManager->Get3DThreshold();
            //if(adc[i]>f3DThreshold)
            //f3DHist->Fill(position.X()/10.,position.Y()/10.,i,adc[i]);

          }
    }

    //if(fSaveTextData)
      //dumpEvent<<position.X()<<" "<<position.Y()<<" "<<position.Z()<<" "<<hit.GetTimeStamp()<<" "<<hit.GetCharge()<<std::endl;

  }
    //////////////////////// Colored Box Drawing ////////////////

    fPadPlane -> Draw("zcol");
    gPad ->Update();
    fPadPlanePal
    = (TPaletteAxis *) fPadPlane->GetListOfFunctions()->FindObject("palette");





    for(Int_t iHit=0; iHit<nHits; iHit++)
    {

    ATHit hit = event->GetHitArray()->at(iHit);
    TVector3 position = hit.GetPosition();
    TVector3 positioncorr = hit.GetPositionCorr();

        if(f3DHitStyle==0){

          Float_t HitBoxYDim = hit.GetCharge()*0.001;
          Float_t HitBoxZDim = 0.05;
          Float_t HitBoxXDim = 0.05;


          fhitBoxSet->AddBox(position.X()/10. - HitBoxXDim/2.0, position.Y()/10., position.Z()/10. - HitBoxZDim/2.0,
          HitBoxXDim,HitBoxYDim,HitBoxZDim); //This coordinates are x,y,z in our system


        }else if(f3DHitStyle==1){

         Float_t HitBoxYDim = hit.GetCharge()*0.0002;
         Float_t HitBoxZDim = hit.GetCharge()*0.0002;
         Float_t HitBoxXDim = hit.GetCharge()*0.0002;



         fhitBoxSet->AddBox(position.X()/10. - HitBoxXDim/2.0, position.Y()/10. - HitBoxYDim/2.0, position.Z()/10. - HitBoxZDim/2.0,
                  HitBoxXDim,HitBoxYDim,HitBoxZDim); //This coordinates are x,y,z in our system


      }

    Float_t xrgb=255,yrgb=0,zrgb=0;
    if(fPadPlanePal){

        Int_t cHit = fPadPlanePal->GetValueColor(hit.GetCharge());
        TColor *hitBoxColor = gROOT->GetColor(cHit);
        hitBoxColor->GetRGB(xrgb,yrgb,zrgb);

    }

         fhitBoxSet->DigitColor(xrgb*255,yrgb*255,zrgb*255, 0);

    }

     /////////////////////// End of colored box drawing ////////////////////////////

    fhitBoxSet->RefitPlex();
    TEveTrans& tHitBoxPos = fhitBoxSet->RefMainTrans();
    tHitBoxPos.SetPos(0.0, 0.0, 0.0);

   //for(Int_t i=0;i<hitSphereArray.size();i++) gEve->AddElement(hitSphereArray[i]);


 if(fIsRawData){
    Int_t nPads = fRawevent->GetNumPads();
    std::cout<<"Num of pads : "<<nPads<<std::endl;

        for(Int_t iPad = 0;iPad<nPads;iPad++){


            ATPad *fPad = fRawevent->GetPad(iPad);
            //std::cout<<"Pad num : "<<iPad<<" Is Valid? : "<<fPad->GetValidPad()<<" Pad num in pad object :"<<fPad->GetPadNum()<<std::endl;
            Int_t *rawadc = fPad->GetRawADC();
            Double_t *adc = fPad->GetADC();
	   // dumpEvent<<TSpad<<fPad->GetPadNum()<<std::endl;

            for(Int_t j=0;j<512;j++){ // TODO: This is limited to 256 pads only. Increment the size of the array and put another option for ATTPC

                if (fPad->GetValidPad() && iPad<256){


                    fPadAll[iPad]->SetBinContent(j,adc[j]);
		    //if(fSaveTextData) dumpEvent<<adc[j]<<"     "<<j<<"     "<<fPad->GetPadNum()<<std::endl;

                }



            }


        }

  }


    gEve -> AddElement(fHitSet);
    gEve -> AddElement(fhitBoxSet);

}


void
ATEventDrawTaskProto::DrawPadWave()
{


        fPadWave = new TH1I("fPadWave","fPadWave",512,0,511);
        gROOT->GetListOfSpecials()->Add(fPadWave);
        fCvsPadWave->cd();
        fPadWave -> Draw();
}

void
ATEventDrawTaskProto::DrawPadPlane()
{

  /*if(fPadPlane)
  {
    fPadPlane->Reset(0);
    return;
  }*/

    fPadPlane = fDetmap->GetATTPCPlane("ATTPC_Proto");
    fCvsPadPlane -> cd();
    fPadPlane -> Draw("zcol");




}

void
ATEventDrawTaskProto::UpdateCvsPadWave()
{
    fCvsPadWave -> Modified();
    fCvsPadWave -> Update();


}

void
ATEventDrawTaskProto::UpdateCvsPadPlane()
{
  fCvsPadPlane -> Modified();
  fCvsPadPlane -> Update();

}

void
ATEventDrawTaskProto::SelectPad(const char *rawevt)
{
    int event = gPad->GetEvent();
    if (event != 11) return; //may be comment this line
    TObject *select = gPad->GetSelected();
    if (!select) return;
    if (select->InheritsFrom(TH2Poly::Class())) {
        TH2Poly *h = (TH2Poly*)select;
        gPad->GetCanvas()->FeedbackMode(kTRUE);
         ATRawEvent* tRawEvent = NULL;
         tRawEvent = (ATRawEvent*)gROOT->GetListOfSpecials()->FindObject(rawevt);
         if(tRawEvent == NULL){
		std::cout<<" = ATEventDrawTaskProto::SelectPad NULL pointer for the ATRawEvent! Please select an event first "<<std::endl;
		return;
	}

        int pyold = gPad->GetUniqueID();
        int px = gPad->GetEventX();
        int py = gPad->GetEventY();
        float uxmin = gPad->GetUxmin();
        float uxmax = gPad->GetUxmax();
        int pxmin = gPad->XtoAbsPixel(uxmin);
        int pxmax = gPad->XtoAbsPixel(uxmax);
        if(pyold) gVirtualX->DrawLine(pxmin,pyold,pxmax,pyold);
        gVirtualX->DrawLine(pxmin,py,pxmax,py);
        gPad->SetUniqueID(py);
        Float_t upx = gPad->AbsPixeltoX(px);
        Float_t upy = gPad->AbsPixeltoY(py);
        Double_t x = gPad->PadtoX(upx);
        Double_t y = gPad->PadtoY(upy);
        Int_t bin = h->FindBin(x,y);
        const char *bin_name = h->GetBinName(bin);
        //std::cout<<" X : "<<x<<"  Y: "<<y<<std::endl;
        //std::cout<<bin_name<<std::endl;
        std::cout<<" =========================="<<std::endl;
        std::cout<<" Bin number selected : "<<bin<<" Bin name :"<<bin_name<<std::endl;
        Bool_t IsValid = kFALSE;

        AtTpcMap *tmap = NULL;
        tmap = (AtTpcMap*)gROOT->GetListOfSpecials()->FindObject("fMap");
        //new AtTpcProtoMap();
        //TString map = "/Users/yassidayyad/fair_install/ATTPCROOT_v2_06042015/scripts/proto.map";
        //tmap->SetProtoMap(map.Data());
        Int_t tPadNum =tmap->BinToPad(bin);
        std::cout<<" Bin : "<<bin<<" to Pad : "<<tPadNum<<std::endl;
        ATPad *tPad = tRawEvent->GetPad(tPadNum,IsValid);
        std::cout<<" Event ID (Select Pad) : "<<tRawEvent->GetEventID()<<std::endl;
        std::cout<<" Raw Event Pad Num "<<tPad->GetPadNum()<<" Is Valid? : "<<IsValid<<std::endl;
        std::cout<<std::endl;
        //TH1D* tPadWaveSub = NULL;
        //tPadWaveSub = new TH1D("tPadWaveSub","tPadWaveSub",512.0,0.0,511.0);
        //tPadWaveSub->SetLineColor(kRed);
        TH1I* tPadWave = NULL;
        tPadWave = (TH1I*)gROOT->GetListOfSpecials()->FindObject("fPadWave");
        Int_t *rawadc = tPad->GetRawADC();
        Double_t *adc = tPad->GetADC();
        if(tPadWave == NULL){
            std::cout<<" = ATEventDrawTask::SelectPad NULL pointer for the TH1I! Please select an event first "<<std::endl;
            return;
	     }
         tPadWave->Reset();
         //tPadWaveSub->Reset();
        for(Int_t i=0;i<512;i++){

			      // tPadWave->SetBinContent(i,rawadc[i]);
                   tPadWave->SetBinContent(i,adc[i]);
         		   //tPadWaveSub->SetBinContent(i,adc[i]);

		    }



        TCanvas *tCvsPadWave = NULL;
        tCvsPadWave = (TCanvas*)gROOT->GetListOfSpecials()->FindObject("fCvsPadWave");
        if(tCvsPadWave == NULL){
            std::cout<<" = ATEventDrawTask::SelectPad NULL pointer for the TCanvas! Please select an event first "<<std::endl;
            return;
        }
        tCvsPadWave->cd();
        tPadWave->Draw();
        //tPadWaveSub->Draw("SAME");
        tCvsPadWave->Update();
    }


}

void
ATEventDrawTaskProto::SetHitAttributes(Color_t color, Size_t size, Style_t style)
{
  fHitColor = color;
  fHitSize = size;
  fHitStyle = style;
}

void
ATEventDrawTaskProto::Set3DHitStyleBar() {f3DHitStyle=0;}

void
ATEventDrawTaskProto::Set3DHitStyleBox() {f3DHitStyle=1;}
