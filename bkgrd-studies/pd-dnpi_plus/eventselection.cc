/***********************************************
* Software developement for WASA-at-COSY
* (c) 2005-2020 The WASA-at-COSY Collaboration
* Aleksander K.                 2016-09
* This software is distributed under the terms
  of the GNU General Public Licence v3.0
*
* Modified 2020-05
***********************************************/

//Selection cuts for the analysis of the simulation results of the background pd -> dnpi+ reaction

#include "Wasa.hh"
#include "CDataManager.hh"
#include "CHistoManager.hh"
#include "CParameterManager.hh"
#include "CLog.hh"
#include "CConst.hh"
#include "EmsEvent.hh"

#include "WHitBank.hh"
#include "WHitScint.hh"
#include "WVertex.hh"
#include "WVertexBank.hh"
#include "WCluster.hh"
#include "WClusterBank.hh"
#include "WClusterChamb.hh"
#include <WClusterFinder.hh>
#include "WTrack.hh"
#include "WTrackBank.hh"
#include "CDTracksSimple.hh"
#include "FDFTHTracks.hh"

#include "TString.h"
#include "TMath.h"
#include <TVector3.h>
#include <TLorentzVector.h>
#include <Riostream.h>
#include <iostream>
#include <fstream>

#include "eventselection.hh"

ClassImp(eventselection);

eventselection::eventselection() {}

eventselection::eventselection(const char * name):CAnalysisModule(name) {

    fDetectorTable = dynamic_cast<CCardWDET*>(gParameterManager->GetParameterObject("CCardWDET","default"));
    //FD table
    kFTH1_old = fDetectorTable->GetDet(CConst::kFTH)->Get1stPlane();
    kFRH1_old = fDetectorTable->GetDet(CConst::kFRH)->Get1stPlane();
    //kFVH1_old = fDetectorTable->GetDet(CConst::kFVH)->Get1stPlane();
    kFWC1_old = fDetectorTable->GetDet(CConst::kFWC)->Get1stPlane();
    //printf("Plane Numbers FTH,FRH,FVH,FWC: %i,%i,%i,%i \n",kFTH1_old,kFRH1_old,kFVH1_old,kFWC1_old);
    printf("Plane Numbers FTH,FRH,FVH,FWC: %i,%i,%i \n",kFTH1_old,kFRH1_old,kFWC1_old);

    //CD table
    //Get PlaneNumbers of first and last Planes of PS and SE
    kPSfirst_old = fDetectorTable->GetDet(CConst::kPSB)->Get1stPlane();     //141
    kPSlast_old = fDetectorTable->GetDet(CConst::kPSF)->Get1stPlane();      //143
    kSEfirst_old = fDetectorTable->GetDet(CConst::kSEB)->Get1stPlane();     //151
    kSElast_old = fDetectorTable->GetDet(CConst::kSEF)->Get1stPlane() + fDetectorTable->GetDet(CConst::kSEF)->GetWasaPlanes()-1;  //174
    gScreen<<kPSfirst_old<<"\t"<<kPSlast_old<<"\t"<<kSEfirst_old<<"\t"<<kSElast_old<<CLog::endl;

    fCDTrackFinder = dynamic_cast<CDTracksSimple*>(gDataManager->GetAnalysisModule("CDTracksSimple","default"));    //CD
    if(fCDTrackFinder!=0) fCDTrackBank = fCDTrackFinder->GetTrackBank();

    fFDTrackFinder = dynamic_cast<FDFTHTracks*>(gDataManager->GetAnalysisModule("FDFTHTracks","default"));          //FD
    if(fFDTrackFinder!=0) fFDTrackBank = fFDTrackFinder->GetTrackBank();

    WTrackFinder *MCTrf = dynamic_cast<WTrackFinder*>(gDataManager->GetAnalysisModule("MCTrackFinder","default"));
    fMCTrackBank  = MCTrf->GetTrackBank();
    fMCVertexBank = MCTrf->GetVertexBank();

    fEventHeader = dynamic_cast<REventWmcHeader*>(gDataManager->GetDataObject("REventWmcHeader","EventHeader"));    //WMC Event header
    fHeader = dynamic_cast<REventHeader*>(gDataManager->GetDataObject("REventHeader","Header"));                    //DATA Event Header

    SetupSpectra(name);

}

eventselection::~eventselection() {}

void eventselection::ProcessEvent() {    //01//

    if (fProcessed) return;
    fProcessed = kTRUE;

    //Int_t RunNumber = SorterOption::GetIntValue("RunNumber");

/////////////////////////////////////////ANALYSIS START/////////////////////////////////////////
/*
    //get event weight. For real data and Pluto, weight should be always one,
    //but for MC data, histograms have to be filled with according event weight
    Double_t ww=1;
    if (gWasa->IsAnalysisMode(Wasa::kMCRaw)||
            gWasa->IsAnalysisMode(Wasa::kMCReco)||
            gWasa->IsAnalysisMode(Wasa::kMC))
        ww=fEventHeader->GetWeight();
    cout<<"ww= "<<ww<<endl;
    //Int_t RunNumber = SorterOption::GetIntValue("RunNumber");
*/

    Double_t beamMom = 0.;
    Double_t Q = 0.;

    ////PARTICLE MASSES////
    const Double_t m_target = 1.875613; //deuteron target mass  [GeV]
    const Double_t m_beam = 0.938272;   //proton beam mass      [GeV]
    const Double_t m_3He = 2.808950;    //3He mass              [GeV]
    const Double_t m_n = 0.939565;      //neutron mass          [GeV]
    const Double_t m_p = 0.938272;      //proton mass           [GeV]
    const Double_t m_d = 1.875613;      //deuteron mass         [GeV]
    const Double_t m_pi = 0.13957;      //charged pion mass     [GeV]
    const Double_t m_pi0 = 0.13497;     //neutral pion mass     [GeV]
    const Double_t m_eta = 0.547853;    //eta mass              [GeV]

    Double_t s_thr = m_3He + m_eta;     //invariant mass on threshold [GeV]

///////////////////////////////GENERATED (TRUE) EVENTS FROM PLUTO///////////////////////////////
//////////////////////////////////////////KINEMATICS////////////////////////////////////////////

    if (gWasa->IsAnalysisMode(Wasa::kMCRaw)||
            gWasa->IsAnalysisMode(Wasa::kMCReco)||
            gWasa->IsAnalysisMode(Wasa::kMC)) {     //A01//

        TVector3 vec_d_MC;
        TVector3 vec_n_MC;
        TVector3 vec_pi_MC;

        TLorentzVector P_d_MC;
        TLorentzVector P_n_MC;
        TLorentzVector P_pi_MC;

        Double_t Ekin_d_lab_MC;
        Double_t E_d_lab_MC;
        Double_t p_d_lab_MC;
        Double_t Theta_d_lab_MC;
        Double_t Phi_d_lab_MC;

        Double_t Ekin_n_lab_MC;
        Double_t E_n_lab_MC;
        Double_t p_n_lab_MC;
        Double_t Theta_n_lab_MC;
        Double_t Phi_n_lab_MC;

        Double_t Ekin_pi_lab_MC;
        Double_t E_pi_lab_MC;
        Double_t p_pi_lab_MC;
        Double_t Theta_pi_lab_MC;
        Double_t Phi_pi_lab_MC;

        ////

        Int_t PType;
        WParticle *part = 0;

        WVertexIter vIt(fMCVertexBank);
        Int_t NrVertex = 0;

        while (WVertex* vert = dynamic_cast<WVertex*>(vIt.Next())) {    //A02//

            NrVertex++;

            for (Int_t particleindex = 0; particleindex < vert->NumberOfParticles(); particleindex++) {     //A03//

                part = vert->GetParticle(particleindex);
                PType = part->GetType();

                //cout<<"NrVertex: "<<NrVertex<<endl;
                //cout<<"particleindex: "<<particleindex<<endl;
                //cout<<"PType: "<<PType<<endl;

                ////LAB////

                //deutron
                if ((NrVertex == 2) && (PType == 45)) {

                    Ekin_d_lab_MC = part->GetEkin();
                    Theta_d_lab_MC = part->GetTheta();
                    Phi_d_lab_MC = part->GetPhi();

                    //cout<<"Ekin_d_lab_MC = "<<Ekin_d_lab_MC<<endl;
                    //cout<<"Theta_d_lab_MC = "<<Theta_d_lab_MC*TMath::RadToDeg()<<endl;
                    //cout<<"Phi_d_lab_MC = "<<Phi_d_lab_MC*TMath::RadToDeg()<<endl;
                    //cout<<"m_d = "<<part->GetMass();<<endl;

                    p_d_lab_MC = TMath::Sqrt(Ekin_d_lab_MC*(Ekin_d_lab_MC + 2*m_d));                //E^2 = m^2 + p^2, E = Ekin + m => p^2 = Ekin^2 + 2*m*Ekin
                    //E_d_lab_MC = TMath::Sqrt(p_d_lab_MC*p_d_lab_MC + m_d*m_d);                    //total energy
                    E_d_lab_MC = TMath::Sqrt(TMath::Power(p_d_lab_MC, 2) + TMath::Power(m_d, 2));   //total energy

                    vec_d_MC.SetMagThetaPhi(p_d_lab_MC,Theta_d_lab_MC,Phi_d_lab_MC);
                    P_d_MC.SetVectM(vec_d_MC,m_d);

                    //histograms
                    hEkin_d_lab_MC->Fill(Ekin_d_lab_MC);
                    hp_d_lab_MC->Fill(p_d_lab_MC);
                    hE_d_lab_MC->Fill(E_d_lab_MC);
                    hTheta_d_lab_MC->Fill(Theta_d_lab_MC*TMath::RadToDeg());
                    hPhi_d_lab_MC->Fill(Phi_d_lab_MC*TMath::RadToDeg());
                    hEkin_vs_Theta_d_lab_MC->Fill(Ekin_d_lab_MC,Theta_d_lab_MC*TMath::RadToDeg());

                }

                //neutron
                if ((NrVertex == 2) && (PType == 13)) {

                    Ekin_n_lab_MC = part->GetEkin();
                    Theta_n_lab_MC = part->GetTheta();
                    Phi_n_lab_MC = part->GetPhi();

                    //cout<<"Ekin_n_lab_MC = "<<Ekin_n_lab_MC<<endl;
                    //cout<<"Theta_n_lab_MC = "<<Theta_n_lab_MC*TMath::RadToDeg()<<endl;
                    //cout<<"Phi_n_lab_MC = "<<Phi_n_lab_MC*TMath::RadToDeg()<<endl;
                    //cout<<"m_n = "<<part->GetMass();<<endl;

                    p_n_lab_MC = TMath::Sqrt(Ekin_n_lab_MC*(Ekin_n_lab_MC + 2*m_n));                //E^2 = m^2 + p^2, E = Ekin + m => p^2 = Ekin^2 + 2*m*Ekin
                    //E_n_lab_MC = TMath::Sqrt(p_n_lab_MC*p_n_lab_MC + m_n*m_n);                    //total energy
                    E_n_lab_MC = TMath::Sqrt(TMath::Power(p_n_lab_MC, 2) + TMath::Power(m_n, 2));   //total energy

                    vec_n_MC.SetMagThetaPhi(p_n_lab_MC,Theta_n_lab_MC,Phi_n_lab_MC);
                    P_n_MC.SetVectM(vec_n_MC,m_n);

                    //histograms
                    hEkin_n_lab_MC->Fill(Ekin_n_lab_MC);
                    hp_n_lab_MC->Fill(p_n_lab_MC);
                    hE_n_lab_MC->Fill(E_n_lab_MC);
                    hTheta_n_lab_MC->Fill(Theta_n_lab_MC*TMath::RadToDeg());
                    hPhi_n_lab_MC->Fill(Phi_n_lab_MC*TMath::RadToDeg());
                    hEkin_vs_Theta_n_lab_MC->Fill(Ekin_n_lab_MC,Theta_n_lab_MC*TMath::RadToDeg());

                }

                //cgarged pion
                if ((NrVertex == 2) && (PType == 8)) {

                    Ekin_pi_lab_MC = part->GetEkin();
                    Theta_pi_lab_MC = part->GetTheta();
                    Phi_pi_lab_MC = part->GetPhi();

                    //cout<<"Ekin_pi_lab_MC = "<<Ekin_pi_lab_MC<<endl;
                    //cout<<"Theta_pi_lab_MC = "<<Theta_pi_lab_MC*TMath::RadToDeg()<<endl;
                    //cout<<"Phi_pi_lab_MC = "<<Phi_pi_lab_MC*TMath::RadToDeg()<<endl;
                    //cout<<"m_pi = "<<part->GetMass();<<endl;

                    p_pi_lab_MC = TMath::Sqrt(Ekin_pi_lab_MC*(Ekin_pi_lab_MC + 2*m_pi));                //E^2 = m^2 + p^2, E = Ekin + m => p^2 = Ekin^2 + 2*m*Ekin
                    //E_pi_lab_MC = TMath::Sqrt(p_pi_lab_MC*p_pi_lab_MC + m_pi*m_pi);                    //total energy
                    E_pi_lab_MC = TMath::Sqrt(TMath::Power(p_pi_lab_MC, 2) + TMath::Power(m_pi, 2));   //total energy

                    vec_pi_MC.SetMagThetaPhi(p_pi_lab_MC,Theta_pi_lab_MC,Phi_pi_lab_MC);
                    P_pi_MC.SetVectM(vec_pi_MC,m_pi);

                    //histograms
                    hEkin_pi_lab_MC->Fill(Ekin_pi_lab_MC);
                    hp_pi_lab_MC->Fill(p_pi_lab_MC);
                    hE_pi_lab_MC->Fill(E_pi_lab_MC);
                    hTheta_pi_lab_MC->Fill(Theta_pi_lab_MC*TMath::RadToDeg());
                    hPhi_pi_lab_MC->Fill(Phi_pi_lab_MC*TMath::RadToDeg());
                    hEkin_vs_Theta_pi_lab_MC->Fill(Ekin_pi_lab_MC,Theta_pi_lab_MC*TMath::RadToDeg());

                }

            }   //A03//

        }   //A02//


        //BEAM KINETIC ENERGY//
        //FOUR-VECTORS//

        TVector3 vec_beam_MC = vec_d_MC + vec_n_MC + vec_pi_MC;

        beamMom = vec_beam_MC.Mag();

        hp_beam_MC->Fill(beamMom);

        TLorentzVector P_b_MC;      //4-vector of the beam
        P_b_MC.SetVectM(vec_beam_MC,m_beam);

        TVector3 vec_target_MC;
        vec_target_MC.SetMagThetaPhi(0.,0.,0.);

        TLorentzVector P_t_MC;      //4-vector of the target
        P_t_MC.SetVectM(vec_target_MC,m_target);

        TLorentzVector P_tot_MC = P_b_MC + P_t_MC;  //total 4-vector

        //Double_t s = TMath::Sqrt(m_beam*m_beam + m_target*m_target + 2*m_target*TMath::Sqrt(m_beam*m_beam + beamMom*beamMom));  //corresponds to Q from Q=-70MeV to Q=30MeV
        Double_t sqMass = TMath::Power(m_beam, 2) + TMath::Power(m_target, 2);
        Double_t beamEnergy = TMath::Sqrt(TMath::Power(beamMom, 2) + TMath::Power(m_beam, 2));
        Double_t s = TMath::Sqrt(sqMass + 2*m_target*beamEnergy);   //corresponds to Q from Q=-70MeV to Q=30MeV

        //cout<<"beamMom = "<<beamMom<<endl;
        //cout<<"s = "<<s<<endl;
        //cout<<"s_thr = "<<s_thr<<endl;

        Q = 1000*(s - s_thr);   //excess energy

        //cout<<"Q= "<<Q<<endl;

        hGenerated_Q->Fill(Q);

////////////////////////////////////////////CM FRAME////////////////////////////////////////////

        TVector3 vec_d_cm_MC;
        TVector3 vec_n_cm_MC;
        TVector3 vec_pi_cm_MC;

        Double_t p_d_cm_MC;
        Double_t E_d_cm_MC;
        Double_t Ekin_d_cm_MC;
        Double_t Theta_d_cm_MC;
        Double_t Phi_d_cm_MC;

        Double_t p_n_cm_MC;
        Double_t E_n_cm_MC;
        Double_t Ekin_n_cm_MC;
        Double_t Theta_n_cm_MC;
        Double_t Phi_n_cm_MC;

        Double_t p_pi_cm_MC;
        Double_t E_pi_cm_MC;
        Double_t Ekin_pi_cm_MC;
        Double_t Theta_pi_cm_MC;
        Double_t Phi_pi_cm_MC;

        ////boost to CM////
        TVector3 b_MC;
        b_MC = P_tot_MC.BoostVector();  //boost to LAB

        ////beam and target////
        P_b_MC.Boost(-b_MC);    //to CM
        P_t_MC.Boost(-b_MC);    //to CM

        ////deuteron////
        P_d_MC.Boost(-b_MC);    //to CM

        vec_d_cm_MC = P_d_MC.Vect();

        p_d_cm_MC = vec_d_cm_MC.Mag();
        //E_d_cm_MC = TMath::Sqrt(p_d_cm_MC*p_d_cm_MC + m_d*m_d);
        E_d_cm_MC = P_d_MC.E();
        Ekin_d_cm_MC = E_d_cm_MC - m_d;
        Theta_d_cm_MC = P_d_MC.Theta();
        Phi_d_cm_MC = P_d_MC.Phi();

        //histograms
        hp_d_cm_MC->Fill(p_d_cm_MC);
        hE_d_cm_MC->Fill(E_d_cm_MC);
        hEkin_d_cm_MC->Fill(Ekin_d_cm_MC);
        hTheta_d_cm_MC->Fill(Theta_d_cm_MC*TMath::RadToDeg());
        hPhi_d_cm_MC->Fill(Phi_d_cm_MC*TMath::RadToDeg());
        hEkin_vs_Theta_d_cm_MC->Fill(Ekin_d_cm_MC,Theta_d_cm_MC*TMath::RadToDeg());

        ////neutron////
        P_n_MC.Boost(-b_MC);    //to CM

        vec_n_cm_MC = P_n_MC.Vect();
        p_n_cm_MC = vec_n_cm_MC.Mag();
        //E_n_cm_MC = TMath::Sqrt(p_n_cm_MC*p_n_cm_MC + m_n*m_n);
        E_n_cm_MC = P_n_MC.E();
        Ekin_n_cm_MC = E_n_cm_MC - m_n;
        Theta_n_cm_MC = P_n_MC.Theta();
        Phi_n_cm_MC = P_n_MC.Phi();

        //histograms
        hp_n_cm_MC->Fill(p_n_cm_MC);
        hE_n_cm_MC->Fill(E_n_cm_MC);
        hEkin_n_cm_MC->Fill(Ekin_n_cm_MC);
        hTheta_n_cm_MC->Fill(Theta_n_cm_MC*TMath::RadToDeg());
        hPhi_n_cm_MC->Fill(Phi_n_cm_MC*TMath::RadToDeg());
        hEkin_vs_Theta_n_cm_MC->Fill(Ekin_n_cm_MC,Theta_n_cm_MC*TMath::RadToDeg());

        ////pion////
        P_pi_MC.Boost(-b_MC);    //to CM

        vec_pi_cm_MC = P_pi_MC.Vect();
        p_pi_cm_MC = vec_pi_cm_MC.Mag();
        //E_pi_cm_MC = TMath::Sqrt(p_pi_cm_MC*p_pi_cm_MC + m_pi*m_pi);
        E_pi_cm_MC = P_pi_MC.E();
        Ekin_pi_cm_MC = E_pi_cm_MC - m_pi;
        Theta_pi_cm_MC = P_pi_MC.Theta();
        Phi_pi_cm_MC = P_pi_MC.Phi();

        //histograms
        hp_pi_cm_MC->Fill(p_pi_cm_MC);
        hE_pi_cm_MC->Fill(E_pi_cm_MC);
        hEkin_pi_cm_MC->Fill(Ekin_pi_cm_MC);
        hTheta_pi_cm_MC->Fill(Theta_pi_cm_MC*TMath::RadToDeg());
        hPhi_pi_cm_MC->Fill(Phi_pi_cm_MC*TMath::RadToDeg());
        hEkin_vs_Theta_pi_cm_MC->Fill(Ekin_pi_cm_MC,Theta_pi_cm_MC*TMath::RadToDeg());

    }   //A01//

//////////////////////////////////////RECONSTRUCTED EVENTS//////////////////////////////////////

    ///////LEVEL 0///////

    hStatistics[0]->Fill(0);

    hp_beam[0][0]->Fill(beamMom);

    hQ[0][0]->Fill(Q);

////////////////////////////////////////////////////////////////////////////////////////////////

    //Forward Detector//
    Int_t NumTrackFD = fFDTrackBank->GetEntries();
    Int_t NumNeutTrackFD = fFDTrackBank->GetEntries(1);     //Neutral Tracks in FD have type 1
    Int_t NumCharTrackFD = fFDTrackBank->GetEntries(2);     //Charged Tracks in FD have type 2

    hTracksFD[0][0]->Fill(NumTrackFD);
    hNeutralTracksFD[0][0]->Fill(NumNeutTrackFD);
    hChargedTracksFD[0][0]->Fill(NumCharTrackFD);

    WTrackIter TrackIterFD(fFDTrackBank);   //define an Iterator for th FD TrackBank
    TrackIterFD.SetType(2);                 //prepare Iterator to look at charged tracks only

    while (WTrack *FDTrack = dynamic_cast<WTrack*> (TrackIterFD.Next())) {

        if(FDTrack->Theta() != 0.125) {

            hTime_FDC[0][0]->Fill(FDTrack->Time());
            hTheta_FDC[0][0]->Fill(FDTrack->Theta()*TMath::RadToDeg());
            hPhi_FDC[0][0]->Fill(FDTrack->Phi()*TMath::RadToDeg());

            //energy deposited in different FD layers
            hEdepFWC1vsFRH1[0][0]->Fill(FDTrack->Edep(kFWC1),FDTrack->Edep(kFRH1));
            hEdepFWC2vsFRH1[0][0]->Fill(FDTrack->Edep(kFWC2),FDTrack->Edep(kFRH1));
            hEdepFTH1vsFRH1[0][0]->Fill(FDTrack->Edep(kFTH1),FDTrack->Edep(kFRH1));
            hEdepFRH1vsFRH2[0][0]->Fill(FDTrack->Edep(kFRH1),FDTrack->Edep(kFRH2));
            hEdepFRH2vsFRH3[0][0]->Fill(FDTrack->Edep(kFRH2),FDTrack->Edep(kFRH3));
            hEdepFWC1vsFRH1FRH2FRH3[0][0]->Fill(FDTrack->Edep(kFWC1),(FDTrack->Edep(kFRH1)+FDTrack->Edep(kFRH2)+FDTrack->Edep(kFRH3)));

        }

    }

    //Central Detector//
    Int_t NumTrackCD = fCDTrackBank->GetEntries();
    Int_t NumNeutTrackCD = fCDTrackBank->GetEntries(11);    //Neutral Tracks in CD have Type 11
    Int_t NumCharTrackCD = fCDTrackBank->GetEntries(12);    //Charged Tracks in CD have Type 12

    hTracksCD[0][0]->Fill(NumTrackCD);
    hNeutralTracksCD[0][0]->Fill(NumNeutTrackCD);
    hChargedTracksCD[0][0]->Fill(NumCharTrackCD);

    WTrackIter TrackIterCD(fCDTrackBank);   //define an Iterator for th CD TrackBank
    TrackIterCD.SetType(12);                //prepare Iterator to look at charged tracks only

    while (WTrack *CDTrack = dynamic_cast<WTrack*> (TrackIterCD.Next())) {

        hTime_CDC[0][0]->Fill(CDTrack->Time());
        hTheta_CDC[0][0]->Fill(CDTrack->Theta()*TMath::RadToDeg());
        hPhi_CDC[0][0]->Fill(CDTrack->Phi()*TMath::RadToDeg());
        hMom_CDC[0][0]->Fill(CDTrack->Momentum());

        hEdepPSBvsSEC[0][0]->Fill((CDTrack->Edep(151,174)),(CDTrack->GetSpecELossPS()));
        hEdepPSBvsSigMom[0][0]->Fill((CDTrack->Momentum()*CDTrack->Charge()),(CDTrack->GetSpecELossPS()));
        hEdepSECvsSigMom[0][0]->Fill((CDTrack->Momentum()*CDTrack->Charge()),(CDTrack->Edep(151,174)));

    }

    TrackIterCD.Reset();        //reset Iterator and
    TrackIterCD.SetType(11);    //prepare for looking at neutral tracks

    while (WTrack *CDTrack = dynamic_cast<WTrack*> (TrackIterCD.Next())) {

        hTime_CDN[0][0]->Fill(CDTrack->Time());
        hTheta_CDN[0][0]->Fill(CDTrack->Theta()*TMath::RadToDeg());
        hPhi_CDN[0][0]->Fill(CDTrack->Phi()*TMath::RadToDeg());
        hMom_CDN[0][0]->Fill(CDTrack->Momentum());
    }

////////////////////////////////////////////////////////////////////////////////////////////////

    //////LEVEL 1//////
    //exactly 1 charged particle in FD & exactly 1 charged particle in CD//
    if((NumCharTrackCD == 1) && (NumCharTrackFD == 1))  {    //B01//

        //hStatistics[0]->Fill(1);

        hNeutralTracksCD[1][0]->Fill(NumNeutTrackCD);
        hChargedTracksCD[1][0]->Fill(NumCharTrackCD);
        hChargedTracksFD[1][0]->Fill(NumCharTrackFD);

        //Forward Detector//
        Double_t EdepFWC1, EdepFWC2;
        Double_t EdepFTH1;
        Double_t EdepFRH1, EdepFRH2, EdepFRH3;
        Double_t TimeFD;
        Double_t ThetaFD_lab;
        Double_t PhiFD_lab;

        WTrackIter FDTrackIter(fFDTrackBank);   //define an iterator for FD TrackBank
        FDTrackIter.SetType(2);                 //prepare Iterator to look at charged tracks only

        while (WTrack *TrackFD = dynamic_cast<WTrack*> (FDTrackIter.Next())) {

            TimeFD = TrackFD->Time();
            ThetaFD_lab = TrackFD->Theta();
            PhiFD_lab = TrackFD->Phi();

            EdepFWC1 = TrackFD->Edep(kFWC1);
            EdepFWC2 = TrackFD->Edep(kFWC2);
            EdepFTH1 = TrackFD->Edep(kFTH1);
            EdepFRH1 = TrackFD->Edep(kFRH1);
            EdepFRH2 = TrackFD->Edep(kFRH2);
            EdepFRH3 = TrackFD->Edep(kFRH3);

        }

        //Central Detector//
        Double_t EdepSEC;
        Double_t EdepPSB;
        Double_t SgnMom;

        Double_t TimeCD;
        Double_t ThetaCD_lab;
        Double_t PhiCD_lab;
        Double_t MomCD_lab;

        WTrackIter CDTrackIter(fCDTrackBank);   //define an Iterator for CD TrackBank
        CDTrackIter.SetType(12);                //prepare Iterator to look at charged tracks only

        while (WTrack *TrackCD = dynamic_cast<WTrack*> (CDTrackIter.Next())) {

            TimeCD = TrackCD->Time();
            ThetaCD_lab = TrackCD->Theta();
            PhiCD_lab = TrackCD->Phi();
            MomCD_lab = TrackCD->Momentum();

            EdepSEC = TrackCD->Edep(151,174);
            EdepPSB = TrackCD->GetSpecELossPS();
            SgnMom = TrackCD->Momentum()*(TrackCD->Charge());

        }

        Double_t Delta_Phi = (PhiFD_lab - PhiCD_lab)*TMath::RadToDeg();
        Double_t Delta_Phi_abs = (TMath::Abs(PhiFD_lab - PhiCD_lab))*TMath::RadToDeg();
        Double_t Delta_Phi_sym = (fmod((2*TMath::Pi() + (PhiFD_lab - PhiCD_lab)), (2*TMath::Pi())))*TMath::RadToDeg();

//////////////////////////////////////////LEVELS 1-2////////////////////////////////////////////

        Bool_t lev[7][1];

        //main level
        lev[1][0] = (ThetaFD_lab != 0.125);

        //cut on scattering angle in FD (3,18) [deg]
        lev[2][0] = (lev[1][0] && (ThetaFD_lab >= 0.052) && (ThetaFD_lab <= 0.314));

        //positively charged particles in CD
        lev[3][0] = (lev[2][0] &&  (SgnMom >= 0.));

        //cut on scattering angle in CD > 50 [deg]
        lev[4][0] = (lev[3][0] && (ThetaCD_lab > 0.873));

	//cut on scattering angle in CD > 60 [deg]
        lev[5][0] = (lev[3][0] && (ThetaCD_lab > 1.047));

        //loop
        for (Int_t l = 1; l < 6; l++) {     //C01//

            if (lev[l][0]) {            //C02//

                hStatistics[0]->Fill(l);

                hQ[l][0]->Fill(Q);

                hTime_FDC[l][0]->Fill(TimeFD);
                hTheta_FDC[l][0]->Fill(ThetaFD_lab*TMath::RadToDeg());
                hPhi_FDC[l][0]->Fill(PhiFD_lab*TMath::RadToDeg());

                hEdepFWC1vsFRH1[l][0]->Fill(EdepFWC1,EdepFRH1);
                hEdepFWC2vsFRH1[l][0]->Fill(EdepFWC2,EdepFRH1);
                hEdepFTH1vsFRH1[l][0]->Fill(EdepFTH1,EdepFRH1);
                hEdepFRH1vsFRH2[l][0]->Fill(EdepFRH1,EdepFRH2);
                hEdepFRH2vsFRH3[l][0]->Fill(EdepFRH2,EdepFRH3);
                hEdepFWC1vsFRH1FRH2FRH3[l][0]->Fill(EdepFWC1,(EdepFRH1 + EdepFRH2 + EdepFRH3));

                hTime_CDC[l][0]->Fill(TimeCD);
                hMom_CDC[l][0]->Fill(MomCD_lab);
                hTheta_CDC[l][0]->Fill(ThetaCD_lab*TMath::RadToDeg());
                hPhi_CDC[l][0]->Fill(PhiCD_lab*TMath::RadToDeg());

                hEdepPSBvsSEC[l][0]->Fill(EdepSEC,EdepPSB);
                hEdepPSBvsSigMom[l][0]->Fill(SgnMom,EdepPSB);
                hEdepSECvsSigMom[l][0]->Fill(SgnMom,EdepSEC);

                hTheta_FDvsTheta_CD[l][0]->Fill(ThetaFD_lab*TMath::RadToDeg(),ThetaCD_lab*TMath::RadToDeg());
                hDeltaTime[l][0]->Fill(TimeFD-TimeCD);

                hDelta_Phi[l][0]->Fill(Delta_Phi);
                hDelta_Phi_abs[l][0]->Fill(Delta_Phi_abs);
                hDelta_Phi_sym[l][0]->Fill(Delta_Phi_sym);

            }   //C02//

        }       //C01//

    }   //B01//

    return;

}   //01//

////////////////////////////////////////////////////////////////////////////////////////////////

void eventselection::SetupSpectra(const char * lpath) {   //02//

    TString h_mc = Form("WMC");

    TString g[7][1];
    g[0][0] = Form("DATA_lev0");
    g[1][0] = Form("DATA_lev1");
    g[2][0] = Form("DATA_lev2");
    g[3][0] = Form("DATA_lev3");
    g[4][0] = Form("DATA_lev4");
    g[5][0] = Form("DATA_lev5");

    TString h_st = Form("Statistics");

////////////////////////////////////TRUE EVENTS (MONTE CARLO)///////////////////////////////////

    hp_beam_MC = new TH1F("hp_beam_MC","",200,1.41,1.65);
    hp_beam_MC->GetXaxis()->SetTitle("p_{beam} [GeV/c]");
    hp_beam_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hp_beam_MC,h_mc);

    hGenerated_Q = new TH1F("hGenerated_Q","",40,-70.,30.);
    hGenerated_Q->GetXaxis()->SetTitle("Q [MeV]");
    hGenerated_Q->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hGenerated_Q,h_mc);

    //deuteron//
    hEkin_d_lab_MC = new TH1F("hEkin_d_lab_MC","",500,0.,1.);
    hEkin_d_lab_MC->GetXaxis()->SetTitle("E^{lab}_{kin} [GeV]");
    hEkin_d_lab_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hEkin_d_lab_MC,h_mc);

    hEkin_d_cm_MC = new TH1F("hEkin_d_cm_MC","",500,0.,2.);
    hEkin_d_cm_MC->GetXaxis()->SetTitle("E^{cm}_{kin} [GeV]");
    hEkin_d_cm_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hEkin_d_cm_MC,h_mc);

    hE_d_lab_MC = new TH1F("hE_d_lab_MC","",500,1.8,3.0);
    hE_d_lab_MC->GetXaxis()->SetTitle("E^{lab} [GeV]");
    hE_d_lab_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hE_d_lab_MC,h_mc);

    hE_d_cm_MC = new TH1F("hE_d_cm_MC","",500,0.,2.1);
    hE_d_cm_MC->GetXaxis()->SetTitle("E^{cm} [GeV]");
    hE_d_cm_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hE_d_cm_MC,h_mc);

    hp_d_lab_MC = new TH1F("hp_d_lab_MC","",500,0.,2.3);
    hp_d_lab_MC->GetXaxis()->SetTitle("p^{lab} [GeV/c]");
    hp_d_lab_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hp_d_lab_MC,h_mc);

    hp_d_cm_MC = new TH1F("hp_d_cm_MC","",500,0.,1.);
    hp_d_cm_MC->GetXaxis()->SetTitle("p^{cm} [GeV/c]");
    hp_d_cm_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hp_d_cm_MC,h_mc);

    hTheta_d_lab_MC = new TH1F("hTheta_d_lab_MC","",360,0.,180.);
    hTheta_d_lab_MC->GetXaxis()->SetTitle("#theta^{lab} [deg]");
    hTheta_d_lab_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hTheta_d_lab_MC,h_mc);

    hTheta_d_cm_MC = new TH1F("hTheta_d_cm_MC","",360,0.,180.);
    hTheta_d_cm_MC->GetXaxis()->SetTitle("#theta^{cm} [deg]");
    hTheta_d_cm_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hTheta_d_cm_MC,h_mc);

    hPhi_d_lab_MC = new TH1F("hPhi_d_lab_MC","",360,0.,360.);
    hPhi_d_lab_MC->GetXaxis()->SetTitle("#phi^{lab} [deg]");
    hPhi_d_lab_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hPhi_d_lab_MC,h_mc);

    hPhi_d_cm_MC = new TH1F("hPhi_d_cm_MC","",360,-180.,180.);
    hPhi_d_cm_MC->GetXaxis()->SetTitle("#phi^{cm} [deg]");
    hPhi_d_cm_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hPhi_d_cm_MC,h_mc);

    hEkin_vs_Theta_d_lab_MC = new TH2F("hEkin_vs_Theta_d_lab_MC","",500,0.,1.0,500,0.,180.);
    hEkin_vs_Theta_d_lab_MC->GetXaxis()->SetTitle("E^{lab}_{kin} [GeV]");
    hEkin_vs_Theta_d_lab_MC->GetYaxis()->SetTitle("#theta^{lab} [deg]");
    gHistoManager->Add(hEkin_vs_Theta_d_lab_MC,h_mc);

    hEkin_vs_Theta_d_cm_MC = new TH2F("hEkin_vs_Theta_d_cm_MC","",500,0.,2.,500,0.,180.);
    hEkin_vs_Theta_d_cm_MC->GetXaxis()->SetTitle("E^{cm}_{kin} [GeV]");
    hEkin_vs_Theta_d_cm_MC->GetYaxis()->SetTitle("#theta^{cm} [deg]");
    gHistoManager->Add(hEkin_vs_Theta_d_cm_MC,h_mc);

    //neutron//
    hEkin_n_lab_MC = new TH1F("hEkin_n_lab_MC","",500,0.,1.);
    hEkin_n_lab_MC->GetXaxis()->SetTitle("E^{lab}_{kin} [GeV]");
    hEkin_n_lab_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hEkin_n_lab_MC,h_mc);

    hEkin_n_cm_MC = new TH1F("hEkin_n_cm_MC","",500,0.,0.3);
    hEkin_n_cm_MC->GetXaxis()->SetTitle("E^{cm}_{kin} [GeV]");
    hEkin_n_cm_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hEkin_n_cm_MC,h_mc);

    hE_n_lab_MC = new TH1F("hE_n_lab_MC","",500,0.8,2.);
    hE_n_lab_MC->GetXaxis()->SetTitle("E^{lab} [GeV]");
    hE_n_lab_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hE_n_lab_MC,h_mc);

    hE_n_cm_MC = new TH1F("hE_n_cm_MC","",500,0.8,1.3);
    hE_n_cm_MC->GetXaxis()->SetTitle("E^{cm} [GeV]");
    hE_n_cm_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hE_n_cm_MC,h_mc);

    hp_n_lab_MC = new TH1F("hp_n_lab_MC","",500,0.,1.5);
    hp_n_lab_MC->GetXaxis()->SetTitle("p^{lab} [GeV/c]");
    hp_n_lab_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hp_n_lab_MC,h_mc);

    hp_n_cm_MC = new TH1F("hp_n_cm_MC","",500,0.,0.8);
    hp_n_cm_MC->GetXaxis()->SetTitle("p^{cm} [GeV/c]");
    hp_n_cm_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hp_n_cm_MC,h_mc);

    hTheta_n_lab_MC = new TH1F("hTheta_n_lab_MC","",360,0.,180.);
    hTheta_n_lab_MC->GetXaxis()->SetTitle("#theta^{lab} [deg]");
    hTheta_n_lab_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hTheta_n_lab_MC,h_mc);

    hTheta_n_cm_MC = new TH1F("hTheta_n_cm_MC","",360,0.,180.);
    hTheta_n_cm_MC->GetXaxis()->SetTitle("#theta^{cm} [deg]");
    hTheta_n_cm_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hTheta_n_cm_MC,h_mc);

    hPhi_n_lab_MC = new TH1F("hPhi_n_lab_MC","",360,0.,360.);
    hPhi_n_lab_MC->GetXaxis()->SetTitle("#phi^{lab} [deg]");
    hPhi_n_lab_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hPhi_n_lab_MC,h_mc);

    hPhi_n_cm_MC = new TH1F("hPhi_n_cm_MC","",360,-180.,180.);
    hPhi_n_cm_MC->GetXaxis()->SetTitle("#phi^{cm} [deg]");
    hPhi_n_cm_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hPhi_n_cm_MC,h_mc);

    hEkin_vs_Theta_n_lab_MC = new TH2F("hEkin_vs_Theta_n_lab_MC","",500,0.,1.,500,0.,180.);
    hEkin_vs_Theta_n_lab_MC->GetXaxis()->SetTitle("E^{lab}_{kin} [GeV]");
    hEkin_vs_Theta_n_lab_MC->GetYaxis()->SetTitle("#theta^{lab} [deg]");
    gHistoManager->Add(hEkin_vs_Theta_n_lab_MC,h_mc);

    hEkin_vs_Theta_n_cm_MC = new TH2F("hEkin_vs_Theta_n_cm_MC","",500,0.,0.3,500,0.,180.);
    hEkin_vs_Theta_n_cm_MC->GetXaxis()->SetTitle("E^{cm}_{kin} [GeV]");
    hEkin_vs_Theta_n_cm_MC->GetYaxis()->SetTitle("#theta^{cm} [deg]");
    gHistoManager->Add(hEkin_vs_Theta_n_cm_MC,h_mc);

    //pion//
    hEkin_pi_lab_MC = new TH1F("hEkin_pi_lab_MC","",500,0.,0.75);
    hEkin_pi_lab_MC->GetXaxis()->SetTitle("E^{lab}_{kin} [GeV]");
    hEkin_pi_lab_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hEkin_pi_lab_MC,h_mc);

    hEkin_pi_cm_MC = new TH1F("hEkin_pi_cm_MC","",500,0.,0.5);
    hEkin_pi_cm_MC->GetXaxis()->SetTitle("E^{cm}_{kin} [GeV]");
    hEkin_pi_cm_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hEkin_pi_cm_MC,h_mc);

    hE_pi_lab_MC = new TH1F("hE_pi_lab_MC","",500,0.,1.);
    hE_pi_lab_MC->GetXaxis()->SetTitle("E^{lab} [GeV]");
    hE_pi_lab_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hE_pi_lab_MC,h_mc);

    hE_pi_cm_MC = new TH1F("hE_pi_cm_MC","",500,0.,0.7);
    hE_pi_cm_MC->GetXaxis()->SetTitle("E^{cm} [GeV]");
    hE_pi_cm_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hE_pi_cm_MC,h_mc);

    hp_pi_lab_MC = new TH1F("hp_pi_lab_MC","",500,0.,1.);
    hp_pi_lab_MC->GetXaxis()->SetTitle("p^{lab} [GeV/c]");
    hp_pi_lab_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hp_pi_lab_MC,h_mc);

    hp_pi_cm_MC = new TH1F("hp_pi_cm_MC","",500,0.,0.6);
    hp_pi_cm_MC->GetXaxis()->SetTitle("p^{cm} [GeV/c]");
    hp_pi_cm_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hp_pi_cm_MC,h_mc);

    hTheta_pi_lab_MC = new TH1F("hTheta_pi_lab_MC","",360,0.,180.);
    hTheta_pi_lab_MC->GetXaxis()->SetTitle("#theta^{lab} [deg]");
    hTheta_pi_lab_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hTheta_pi_lab_MC,h_mc);

    hTheta_pi_cm_MC = new TH1F("hTheta_pi_cm_MC","",360,0.,180.);
    hTheta_pi_cm_MC->GetXaxis()->SetTitle("#theta^{cm} [deg]");
    hTheta_pi_cm_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hTheta_pi_cm_MC,h_mc);

    hPhi_pi_lab_MC = new TH1F("hPhi_pi_lab_MC","",360,0.,360.);
    hPhi_pi_lab_MC->GetXaxis()->SetTitle("#phi^{lab} [deg]");
    hPhi_pi_lab_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hPhi_pi_lab_MC,h_mc);

    hPhi_pi_cm_MC = new TH1F("hPhi_pi_cm_MC","",360,-180.,180.);
    hPhi_pi_cm_MC->GetXaxis()->SetTitle("#phi^{cm} [deg]");
    hPhi_pi_cm_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hPhi_pi_cm_MC,h_mc);

    hEkin_vs_Theta_pi_lab_MC = new TH2F("hEkin_vs_Theta_pi_lab_MC","",500,0.,0.75,500,0.,180.);
    hEkin_vs_Theta_pi_lab_MC->GetXaxis()->SetTitle("E^{lab}_{kin} [GeV]");
    hEkin_vs_Theta_pi_lab_MC->GetYaxis()->SetTitle("#theta^{lab} [deg]");
    gHistoManager->Add(hEkin_vs_Theta_pi_lab_MC,h_mc);

    hEkin_vs_Theta_pi_cm_MC = new TH2F("hEkin_vs_Theta_pi_cm_MC","",500,0.,0.5,500,0.,180.);
    hEkin_vs_Theta_pi_cm_MC->GetXaxis()->SetTitle("E^{cm}_{kin} [GeV]");
    hEkin_vs_Theta_pi_cm_MC->GetYaxis()->SetTitle("#theta^{cm} [deg]");
    gHistoManager->Add(hEkin_vs_Theta_pi_cm_MC,h_mc);

//////////////////////////////////////////RECONSTRUCTED/////////////////////////////////////////

    ////DATA: LEVEL 0////
    hp_beam[0][0] = new TH1F("hp_beam_lev0","",200,1.41,1.65);
    hp_beam[0][0]->GetXaxis()->SetTitle("p_{beam} [GeV/c]");
    hp_beam[0][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hp_beam[0][0],g[0][0]);

    hQ[0][0] = new TH1F("hQ_lev0","",40,-70.,30.);
    hQ[0][0]->GetXaxis()->SetTitle("Q [MeV]");
    hQ[0][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hQ[0][0],g[0][0]);

    //FD
    hTracksFD[0][0] = new TH1F("hTracksFD_lev0","",11,-0.5,10.5);
    hTracksFD[0][0]->GetXaxis()->SetTitle("Tracks in FD");
    hTracksFD[0][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hTracksFD[0][0],g[0][0]);

    hNeutralTracksFD[0][0] = new TH1F("hNeutralTracksFD_lev0","",11,-0.5,10.5);
    hNeutralTracksFD[0][0]->GetXaxis()->SetTitle("Neutral Tracks in FD");
    hNeutralTracksFD[0][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hNeutralTracksFD[0][0],g[0][0]);

    hChargedTracksFD[0][0] = new TH1F("hChargedTracksFD_lev0","",11,-0.5,10.5);
    hChargedTracksFD[0][0]->GetXaxis()->SetTitle("Charged Tracks in FD");
    hChargedTracksFD[0][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hChargedTracksFD[0][0],g[0][0]);

    hEdepFWC1vsFRH1[0][0] = new TH2F("hEdepFWC1vsFRH1_lev0","",1000,0.,0.01,1000,0.,0.5);
    hEdepFWC1vsFRH1[0][0]->GetXaxis()->SetTitle("Edep(FWC1) [GeV]");
    hEdepFWC1vsFRH1[0][0]->GetYaxis()->SetTitle("Edep(FRH1) [GeV]");
    gHistoManager->Add(hEdepFWC1vsFRH1[0][0],g[0][0]);

    hEdepFWC2vsFRH1[0][0] = new TH2F("hEdepFWC2vsFRH1_lev0","",1000,0.,0.01,1000,0.,0.5);
    hEdepFWC2vsFRH1[0][0]->GetXaxis()->SetTitle("Edep(FWC2) [GeV]");
    hEdepFWC2vsFRH1[0][0]->GetYaxis()->SetTitle("Edep(FRH1) [GeV]");
    gHistoManager->Add(hEdepFWC2vsFRH1[0][0],g[0][0]);

    hEdepFTH1vsFRH1[0][0] = new TH2F("hEdepFTH1vsFRH1_lev0","",1000,0.,0.02,1000,0.,0.5);
    hEdepFTH1vsFRH1[0][0]->GetXaxis()->SetTitle("Edep(FTH1) [GeV]");
    hEdepFTH1vsFRH1[0][0]->GetYaxis()->SetTitle("Edep(FRH1) [GeV]");
    gHistoManager->Add(hEdepFTH1vsFRH1[0][0],g[0][0]);

    hEdepFRH1vsFRH2[0][0] = new TH2F("hEdepFRH1vsFRH2_lev0","",1000,0.,0.5,1000,0.,0.5);
    hEdepFRH1vsFRH2[0][0]->GetXaxis()->SetTitle("Edep(FRH1) [GeV]");
    hEdepFRH1vsFRH2[0][0]->GetYaxis()->SetTitle("Edep(FRH2) [GeV]");
    gHistoManager->Add(hEdepFRH1vsFRH2[0][0],g[0][0]);

    hEdepFRH2vsFRH3[0][0] = new TH2F("hEdepFRH2vsFRH3_lev0","",1000,0.,0.5,1000,0.,0.5);
    hEdepFRH2vsFRH3[0][0]->GetXaxis()->SetTitle("Edep(FRH2) [GeV]");
    hEdepFRH2vsFRH3[0][0]->GetYaxis()->SetTitle("Edep(FRH3) [GeV]");
    gHistoManager->Add(hEdepFRH2vsFRH3[0][0],g[0][0]);

    hEdepFWC1vsFRH1FRH2FRH3[0][0] = new TH2F("hEdepFWC1vsFRH1FRH2FRH3_lev0","",1000,0.,0.01,1000,0.,1.);
    hEdepFWC1vsFRH1FRH2FRH3[0][0]->GetXaxis()->SetTitle("Edep(FWC1) [GeV]");
    hEdepFWC1vsFRH1FRH2FRH3[0][0]->GetYaxis()->SetTitle("Edep(FRH1+FRH2+FRH3) [GeV]");
    gHistoManager->Add(hEdepFWC1vsFRH1FRH2FRH3[0][0],g[0][0]);

    hTheta_FDC[0][0] = new TH1F("hThetaFDC_lev0","",250,0.,25.);
    hTheta_FDC[0][0]->GetXaxis()->SetTitle("#theta^{lab}_{FD} [deg]");
    hTheta_FDC[0][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hTheta_FDC[0][0],g[0][0]);

    hPhi_FDC[0][0] = new TH1F("hPhiFDC_lev0","",360,-180.,180.);
    hPhi_FDC[0][0]->GetXaxis()->SetTitle("#phi^{lab}_{FD} [deg]");
    hPhi_FDC[0][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hPhi_FDC[0][0],g[0][0]);

    hTime_FDC[0][0] = new TH1F("hTimeFDC_lev0","",1000,0.,2500.);
    hTime_FDC[0][0]->GetXaxis()->SetTitle("Time_{FD} [ns]");
    hTime_FDC[0][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hTime_FDC[0][0],g[0][0]);

    //CD
    hTracksCD[0][0] = new TH1F("hTracksCD_lev0","",11,-0.5,10.5);
    hTracksCD[0][0]->GetXaxis()->SetTitle("Tracks in CD");
    hTracksCD[0][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hTracksCD[0][0],g[0][0]);

    hNeutralTracksCD[0][0] = new TH1F("hNeutralTracksCD_lev0","",11,-0.5,10.5);
    hNeutralTracksCD[0][0]->GetXaxis()->SetTitle("Neutral Tracks in CD");
    hNeutralTracksCD[0][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hNeutralTracksCD[0][0],g[0][0]);

    hChargedTracksCD[0][0] = new TH1F("hChargedTracksCD_lev0","",11,-0.5,10.5);
    hChargedTracksCD[0][0]->GetXaxis()->SetTitle("Charged Tracks in CD");
    hChargedTracksCD[0][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hChargedTracksCD[0][0],g[0][0]);

    //charged in CD
    hEdepPSBvsSEC[0][0] = new TH2F("hEdepPSBvsSEC_lev0","",1000,0.,0.75,1000,0.,0.025);
    hEdepPSBvsSEC[0][0]->GetXaxis()->SetTitle("E_{dep(SEC)} [GeV]");
    hEdepPSBvsSEC[0][0]->GetYaxis()->SetTitle("E_{dep(PSB)} [GeV]");
    gHistoManager->Add(hEdepPSBvsSEC[0][0],g[0][0]);

    hEdepPSBvsSigMom[0][0] = new TH2F("hEdepPSBvsSigMom_lev0","",1000,-2.5,2.5,1000,0.,0.025);
    hEdepPSBvsSigMom[0][0]->GetXaxis()->SetTitle("Momentum [GeV/c]");
    hEdepPSBvsSigMom[0][0]->GetYaxis()->SetTitle("E_{dep(PSB)} [GeV]");
    gHistoManager->Add(hEdepPSBvsSigMom[0][0],g[0][0]);

    hEdepSECvsSigMom[0][0] = new TH2F("hEdepSECvsSigMom_lev0","",1000,-2.5,2.5,1000,0.,0.75);
    hEdepSECvsSigMom[0][0]->GetXaxis()->SetTitle("Momentum [GeV/c]");
    hEdepSECvsSigMom[0][0]->GetYaxis()->SetTitle("E_{dep(PSB)} [GeV]");
    gHistoManager->Add(hEdepSECvsSigMom[0][0],g[0][0]);

    hMom_CDC[0][0] = new TH1F("hMomentumCDC_lev0","",750,0.,2.5);
    hMom_CDC[0][0]->GetXaxis()->SetTitle("p_{CD} [GeV/c]");
    hMom_CDC[0][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hMom_CDC[0][0],g[0][0]);

    hTheta_CDC[0][0] = new TH1F("hThetaCDC_lev0","",360,0.,180.);
    hTheta_CDC[0][0]->GetXaxis()->SetTitle("#theta_{CD} [deg]");
    hTheta_CDC[0][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hTheta_CDC[0][0],g[0][0]);

    hPhi_CDC[0][0] = new TH1F("hPhiCDC_lev0","",360,-180.,180.);
    hPhi_CDC[0][0]->GetXaxis()->SetTitle("#phi_{CD} [deg]");
    hPhi_CDC[0][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hPhi_CDC[0][0],g[0][0]);

    hTime_CDC[0][0] = new TH1F("hTimeCDC_lev0","",1000,0.,2500.);
    hTime_CDC[0][0]->GetXaxis()->SetTitle("Time_{CD} [ns]");
    hTime_CDC[0][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hTime_CDC[0][0],g[0][0]);

    //neutral in CD
    hMom_CDN[0][0] = new TH1F("hMomentumCDN_lev0","",500,0.,0.8);
    hMom_CDN[0][0]->GetXaxis()->SetTitle("p_{CD} [GeV/c]");
    hMom_CDN[0][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hMom_CDN[0][0],g[0][0]);

    hTheta_CDN[0][0]= new TH1F("hThetaCDN_lev0","",360,0.,180.);
    hTheta_CDN[0][0]->GetXaxis()->SetTitle("#theta_{CD} [deg]");
    hTheta_CDN[0][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hTheta_CDN[0][0],g[0][0]);

    hPhi_CDN[0][0] = new TH1F("hPhiCDN_lev0","",360,-180.,180.);
    hPhi_CDN[0][0]->GetXaxis()->SetTitle("#phi_{CD} [deg]");
    hPhi_CDN[0][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hPhi_CDN[0][0],g[0][0]);

    hTime_CDN[0][0] = new TH1F("hTimeCDN_lev0","",1000,0.,2500.);
    hTime_CDN[0][0]->GetXaxis()->SetTitle("Time_{CD} [ns]");
    hTime_CDN[0][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hTime_CDN[0][0],g[0][0]);

////DATA: LEVELS 1-4////

    hNeutralTracksCD[1][0] = new TH1F("hNeutralTracksCD_lev1","",11,-0.5,10.5);
    hNeutralTracksCD[1][0]->GetXaxis()->SetTitle("Neutral Tracks in CD");
    hNeutralTracksCD[1][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hNeutralTracksCD[1][0],g[1][0]);

    hChargedTracksCD[1][0] = new TH1F("hChargedTracksCD_lev1","",11,-0.5,10.5);
    hChargedTracksCD[1][0]->GetXaxis()->SetTitle("Charged Tracks in CD");
    hChargedTracksCD[1][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hChargedTracksCD[1][0],g[1][0]);

    hChargedTracksFD[1][0] = new TH1F("hChargedTracksFD_lev1","",11,-0.5,10.5);
    hChargedTracksFD[1][0]->GetXaxis()->SetTitle("Charged Tracks in FD");
    hChargedTracksFD[1][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hChargedTracksFD[1][0],g[1][0]);

    for (Int_t l = 1; l < 6; l++) {
        TString hname = Form("hQ_lev%d",l);
        hQ[l][0] = new TH1F(hname,"",40,-70.,30.);
        hQ[l][0]->GetXaxis()->SetTitle("Q [MeV]");
        hQ[l][0]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hQ[l][0],g[l][0]);
    }

    //FD//
    for (Int_t l = 1; l < 6; l++) {
        TString hname = Form("hEdepFWC1vsFRH1_lev%d",l);
        hEdepFWC1vsFRH1[l][0] =new TH2F(hname,"",1000,0.,0.01,1000,0.,0.5);
        hEdepFWC1vsFRH1[l][0]->GetXaxis()->SetTitle("Edep(FWC1) [GeV]");
        hEdepFWC1vsFRH1[l][0]->GetYaxis()->SetTitle("Edep(FRH1) [GeV]");
        gHistoManager->Add(hEdepFWC1vsFRH1[l][0],g[l][0]);
    }

    for (Int_t l = 1; l < 6; l++) {
        TString hname = Form("hEdepFWC2vsFRH1_lev%d",l);
        hEdepFWC2vsFRH1[l][0] =new TH2F(hname,"",1000,0.,0.01,1000,0.,0.5);
        hEdepFWC2vsFRH1[l][0]->GetXaxis()->SetTitle("Edep(FWC2) [GeV]");
        hEdepFWC2vsFRH1[l][0]->GetYaxis()->SetTitle("Edep(FRH1) [GeV]");
        gHistoManager->Add(hEdepFWC2vsFRH1[l][0],g[l][0]);
    }

    for (Int_t l = 1; l < 6; l++) {
        TString hname = Form("hEdepFTH1vsFRH1_lev%d",l);
        hEdepFTH1vsFRH1[l][0] =new TH2F(hname,"",1000,0.,0.02,1000,0.,0.5);
        hEdepFTH1vsFRH1[l][0]->GetXaxis()->SetTitle("Edep(FTH1) [GeV]");
        hEdepFTH1vsFRH1[l][0]->GetYaxis()->SetTitle("Edep(FRH1) [GeV]");
        gHistoManager->Add(hEdepFTH1vsFRH1[l][0],g[l][0]);
    }

    for (Int_t l = 1; l < 6; l++) {
        TString hname = Form("hEdepFRH1vsFRH2_lev%d",l);
        hEdepFRH1vsFRH2[l][0] =new TH2F(hname,"",1000,0.,0.5,1000,0.,0.5);
        hEdepFRH1vsFRH2[l][0]->GetXaxis()->SetTitle("Edep(FRH1) [GeV]");
        hEdepFRH1vsFRH2[l][0]->GetYaxis()->SetTitle("Edep(FRH2) [GeV]");
        gHistoManager->Add(hEdepFRH1vsFRH2[l][0],g[l][0]);
    }

    for (Int_t l = 1; l < 6; l++) {
        TString hname = Form("hEdepFRH2vsFRH3_lev%d",l);
        hEdepFRH2vsFRH3[l][0] =new TH2F(hname,"",1000,0.,0.5,1000,0.,0.5);
        hEdepFRH2vsFRH3[l][0]->GetXaxis()->SetTitle("Edep(FRH2) [GeV]");
        hEdepFRH2vsFRH3[l][0]->GetYaxis()->SetTitle("Edep(FRH3) [GeV]");
        gHistoManager->Add(hEdepFRH2vsFRH3[l][0],g[l][0]);
    }

    for (Int_t l = 1; l < 6; l++) {
        TString hname = Form("hEdepFWC1vsFRH1FRH2FRH3_lev%d",l);
        hEdepFWC1vsFRH1FRH2FRH3[l][0] =new TH2F(hname,"",1000,0.,0.01,1000,0.,1.);
        hEdepFWC1vsFRH1FRH2FRH3[l][0]->GetXaxis()->SetTitle("Edep(FWC1) [GeV]");
        hEdepFWC1vsFRH1FRH2FRH3[l][0]->GetYaxis()->SetTitle("Edep(FRH1+FRH2+FRH3) [GeV]");
        gHistoManager->Add(hEdepFWC1vsFRH1FRH2FRH3[l][0],g[l][0]);
    }

    for (Int_t l = 1; l < 6; l++) {
        TString hname = Form("hThetaFD_lev%d",l);
        hTheta_FDC[l][0] = new TH1F(hname,"",250,0.,25.);
        hTheta_FDC[l][0]->GetXaxis()->SetTitle("#theta^{lab}_{FD} [deg]");
        hTheta_FDC[l][0]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hTheta_FDC[l][0],g[l][0]);
    }

    for (Int_t l = 1; l < 6; l++) {
        TString hname = Form("hPhiFD_lev%d",l);
        hPhi_FDC[l][0] = new TH1F(hname,"",360,-180.,180.);
        hPhi_FDC[l][0]->GetXaxis()->SetTitle("#phi^{lab}_{FD} [deg]");
        hPhi_FDC[l][0]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hPhi_FDC[l][0],g[l][0]);
    }

    for (Int_t l = 1; l < 6; l++) {
        TString hname = Form("hTimeFD_lev%d",l);
        hTime_FDC[l][0] = new TH1F(hname,"",1000,0.,2500.);
        hTime_FDC[l][0]->GetXaxis()->SetTitle("Time_{FD} [ns]");
        hTime_FDC[l][0]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hTime_FDC[l][0],g[l][0]);
    }

    //charged in CD//
    for (Int_t l = 1; l < 6; l++) {
        TString hname = Form("hEdepPSBvsSEC_lev%d",l);
        hEdepPSBvsSEC[l][0] = new TH2F(hname,"",1000,0.,0.75,1000,0.,0.025);
        hEdepPSBvsSEC[l][0]->GetXaxis()->SetTitle("E_{dep(SEC)} [GeV]");
        hEdepPSBvsSEC[l][0]->GetYaxis()->SetTitle("E_{dep(PSB)} [GeV]");
        gHistoManager->Add(hEdepPSBvsSEC[l][0],g[l][0]);
    }

    for (Int_t l = 1; l < 6; l++) {
        TString hname = Form("hEdepPSBvsSigMom_lev%d",l);
        hEdepPSBvsSigMom[l][0] = new TH2F(hname,"",1000,-2.5,2.5,1000,0.,0.025);
        hEdepPSBvsSigMom[l][0]->GetXaxis()->SetTitle("Momentum [GeV/c]");
        hEdepPSBvsSigMom[l][0]->GetYaxis()->SetTitle("E_{dep(SEC)} [GeV]");
        gHistoManager->Add(hEdepPSBvsSigMom[l][0],g[l][0]);
    }

    for (Int_t l = 1; l < 6; l++) {
        TString hname = Form("hEdepSECvsSigMom_lev%d",l);
        hEdepSECvsSigMom[l][0] = new TH2F(hname,"",1000,-2.5,2.5,1000,0.,0.75);
        hEdepSECvsSigMom[l][0]->GetXaxis()->SetTitle("Momentum [GeV/c]");
        hEdepSECvsSigMom[l][0]->GetYaxis()->SetTitle("E_{dep(PSB)} [GeV]");
        gHistoManager->Add(hEdepSECvsSigMom[l][0],g[l][0]);
    }

    for (Int_t l = 1; l < 6; l++) {
        TString hname = Form("hMomentumCD_lev%d",l);
        hMom_CDC[l][0] = new TH1F(hname,"",750,0.,2.5);
        hMom_CDC[l][0]->GetXaxis()->SetTitle("p^{lab}_{p} [GeV/c]");
        hMom_CDC[l][0]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hMom_CDC[l][0],g[l][0]);
    }

    for (Int_t l = 1; l < 6; l++) {
        TString hname = Form("hThetaCD_lev%d",l);
        hTheta_CDC[l][0] = new TH1F(hname,"",360,0.,180.);
        hTheta_CDC[l][0]->GetXaxis()->SetTitle("#theta^{lab}_{CD} [deg]");
        hTheta_CDC[l][0]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hTheta_CDC[l][0],g[l][0]);
    }

    for (Int_t l = 1; l < 6; l++) {
        TString hname = Form("hPhiCD_lev%d",l);
        hPhi_CDC[l][0] = new TH1F(hname,"",360,-180.,180.);
        hPhi_CDC[l][0]->GetXaxis()->SetTitle("#phi^{lab}_{CD} [deg]");
        hPhi_CDC[l][0]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hPhi_CDC[l][0],g[l][0]);
    }

    for (Int_t l = 1; l < 6; l++) {
        TString hname = Form("hTimeCD_lev%d",l);
        hTime_CDC[l][0] = new TH1F(hname,"",1000,0.,2500.);
        hTime_CDC[l][0]->GetXaxis()->SetTitle("Time_{CD} [ns]");
        hTime_CDC[l][0]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hTime_CDC[l][0],g[l][0]);
    }

    for (Int_t l = 1; l < 6; l++) {
        TString hname = Form("hThetaFDvsThetaCD_lev%d",l);
        hTheta_FDvsTheta_CD[l][0] = new TH2F(hname,"",540,0.,20.,540,0.,180.);
        hTheta_FDvsTheta_CD[l][0]->GetXaxis()->SetTitle("#theta_{FD} [deg]");
        hTheta_FDvsTheta_CD[l][0]->GetYaxis()->SetTitle("#theta_{CD} [deg]");
        gHistoManager->Add(hTheta_FDvsTheta_CD[l][0],g[l][0]);
    }

    for (Int_t l = 1; l < 6; l++) {
        TString hname = Form("hTimeDifference_lev%d",l);
        hDeltaTime[l][0] = new TH1F(hname,"",1000,-50.,50.);
        hDeltaTime[l][0]->GetXaxis()->SetTitle("Time difference [ns]");
        hDeltaTime[l][0]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hDeltaTime[l][0],g[l][0]);
    }

    for (Int_t l = 1; l < 6; l++) {
        TString hname = Form("hDeltaPhi_lev%d",l);
        hDelta_Phi[l][0] = new TH1F(hname,"",360,0.,360.);
        hDelta_Phi[l][0]->GetXaxis()->SetTitle("#Delta#phi [deg]");
        hDelta_Phi[l][0]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hDelta_Phi[l][0],g[l][0]);
    }

    for (Int_t l = 1; l < 6; l++) {
        TString hname = Form("hDeltaPhi_abs_lev%d",l);
        hDelta_Phi_abs[l][0] = new TH1F(hname,"",360,0.,360.);
        hDelta_Phi_abs[l][0]->GetXaxis()->SetTitle("|#Delta#phi| [deg]");
        hDelta_Phi_abs[l][0]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hDelta_Phi_abs[l][0],g[l][0]);
    }

    for (Int_t l = 1; l < 6; l++) {
        TString hname = Form("hDeltaPhi_sym_lev%d",l);
        hDelta_Phi_sym[l][0] = new TH1F(hname,"",360,0.,360.);
        hDelta_Phi_sym[l][0]->GetXaxis()->SetTitle("(2#pi+#Delta#phi)mod2#pi [deg]");
        hDelta_Phi_sym[l][0]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hDelta_Phi_sym[l][0],g[l][0]);
    }


    ////statistics [h_st]/////
    hStatistics[0] = new TH1F("Statistics","",11,-0.5,10.5);
    hStatistics[0]->GetXaxis()->SetTitle("cut");
    hStatistics[0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hStatistics[0],h_st);

    return;

}   //02//

void eventselection::Clear(Option_t *option){
    fProcessed = kFALSE;
    return;
}

void eventselection::Print(Option_t *option){
    return;
}

void eventselection::UserCommand(CCommand * command){
    return;
}
