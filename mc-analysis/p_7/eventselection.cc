/***********************************************
* Software developement for WASA-at-COSY
* (c) 2005-2020 The WASA-at-COSY Collaboration
* Aleksander K.                 2016-09
* This software is distributed under the terms
  of the GNU General Public Licence v3.0
*
* Modified 2020-05
***********************************************/

//Selection cuts for the analysis of the simulation results of the quasi-free pd -> ppn_spectator reaction

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

////////////////////////////file with beam momentum and scattering angle////////////////////////////

    beam_file.open("/data7/users/khreptak/SIMULATION_OUTPUT/PLUTO_OUTPUT/ppn_qf/ProtonVariables/ProtonVariables-PARIS-7.dat");
    //beam_file.open("/data7/users/khreptak/SIMULATION_OUTPUT/PLUTO_OUTPUT/ppn_qf/ProtonVariables/ProtonVariables-CDBONN-7.dat");
    Double_t p_pp_beam_MC = 0.;
    Double_t Theta_scatt_MC = 0.;
    for (Int_t i = 1; i < 1000001; i++) {
        beam_file>>p_pp_beam_MC>>Theta_scatt_MC;
        p_pp_beam_lab[i] = p_pp_beam_MC;
        Theta_scatt_cm[i] = Theta_scatt_MC;
    }

//////////////file with differential cross section from SAID for Theta_scatt=(0,180)deg//////////////

    for (Int_t th = 0; th < 181; th++) {
        TString name = Form("../SAID_data/Theta_%d.txt", th);
        said_file = fopen(name, "r");
        gCrossSection[th] = new TGraph();
        Float_t p_pbeam = 0.;
        Float_t XS = 0.;
        Int_t i = -1;
        while (!feof(said_file)) {
            i++;
            fscanf(said_file, "%f %f\n", &p_pbeam, &XS);
            p_pbeam = p_pbeam/1000.;
            gCrossSection[th]->SetPoint(i, p_pbeam, XS);
        }
    }
/*
/////////////////////////////////time in cycle to beam momentum/////////////////////////////////

    BMfile = fopen("Time.2.PBeam.calibration.dat", "r");    //file with beam momentum
    gBeamMomentum = new TGraph();
    Int_t time_cycle;
    Float_t beam_mom;
    Int_t i = -1;
    while (!feof(BMfile)) {
        i++;
        fscanf(BMfile, "%d %f\n", &time_cycle, &beam_mom);
        beam_mom = beam_mom/1000.;
        gBeamMomentum->SetPoint(i, time_cycle, beam_mom);
    }
*/
////////////////////////////////////////////////////////////////////////////////////////////////

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

Int_t counter = 0;

void eventselection::ProcessEvent() {    //01//

    if (fProcessed) return;
    fProcessed = kTRUE;

    counter++;
    //cout<<"counter = "<<counter<<endl;
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

    //beam momentum offset for DATA analysis
    //const Double_t pbeam_offset = 0.004;      //(from O. Rundel)
    //const Double_t pbeam_offset = 0.;

    Double_t E_calib = 1.;      //energy correction factor for WMC
    //Double_t E_calib = 1.526;   //energy correction factor for DATA

    ////PARTICLE MASSES////
    const Double_t m_target = 1.875613; //deuteron target mass  [GeV]
    const Double_t m_beam = 0.938272;   //proton beam mass      [GeV]
    const Double_t m_3He = 2.808950;    //3He mass              [GeV]
    const Double_t m_n = 0.939565;      //neutron mass          [GeV]
    const Double_t m_p = 0.938272;      //proton mass           [GeV]
    const Double_t m_d = 1.875613;      //deuteron mass         [GeV]
    const Double_t m_pi0 = 0.13497;     //pion mass             [GeV]
    const Double_t m_eta = 0.547853;    //eta mass              [GeV]

    Double_t s_thr = m_3He + m_eta;     //invariant mass on threshold [GeV]

    //cross section
    Double_t XS_bilin = 0.; //for bilinear interpolation
    Double_t XS_lin = 0.;   //for linear interpolation
    Double_t XS_prx = 0.;   //for proximate points searching
    Double_t XS_sq = 0.;    //square of XS (bilinear) for statistical uncertainties calculations

///////////////////////////////GENERATED (TRUE) EVENTS FROM PLUTO///////////////////////////////
//////////////////////////////////////////KINEMATICS////////////////////////////////////////////

    if (gWasa->IsAnalysisMode(Wasa::kMCRaw)||
            gWasa->IsAnalysisMode(Wasa::kMCReco)||
            gWasa->IsAnalysisMode(Wasa::kMC)) {     //A01//

        TVector3 vec_n_MC;
        TVector3 vec_p_MC[2];

        TLorentzVector P_n_MC;
        TLorentzVector P_p_MC[2];

        Double_t Ekin_n_lab_MC;
        Double_t E_n_lab_MC;
        Double_t p_n_lab_MC;
        Double_t Theta_n_lab_MC;
        Double_t Phi_n_lab_MC;

        Double_t Ekin_p_lab_MC[2];
        Double_t E_p_lab_MC[2];
        Double_t p_p_lab_MC[2];
        Double_t Theta_p_lab_MC[2];
        Double_t Phi_p_lab_MC[2];

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

                //neutron
                if ((NrVertex == 1) && (PType == 13) && (particleindex == 0)) {

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

                //proton 1//
                if ((NrVertex == 1) && (PType == 14)&& (particleindex == 1)) {

                    Ekin_p_lab_MC[0] = part->GetEkin();
                    Theta_p_lab_MC[0] = part->GetTheta();
                    Phi_p_lab_MC[0] = part->GetPhi();

                    //cout<<"Ekin_p1_lab_MC = "<<Ekin_p_lab_MC[0]<<endl;
                    //cout<<"Theta_p1_lab_MC = "<<Theta_p_lab_MC[0]*TMath::RadToDeg()<<endl;
                    //cout<<"Phi_p1_lab_MC = "<<Phi_p_lab_MC[0]*TMath::RadToDeg()<<endl;
                    //cout<<"m_p1 = "<<part->GetMass();<<endl;

                    p_p_lab_MC[0] = TMath::Sqrt(Ekin_p_lab_MC[0]*(Ekin_p_lab_MC[0] + 2*m_p));           //E^2 = m^2 + p^2, E = Ekin + m => p^2 = Ekin^2 + 2*m*Ekin
                    //E_p_lab_MC[0] = TMath::Sqrt(p_p_lab_MC[0]*p_p_lab_MC[0] + m_p*m_p);               //total energy
                    E_p_lab_MC[0] = TMath::Sqrt(TMath::Power(p_p_lab_MC[0], 2) + TMath::Power(m_p, 2)); //total energy

                    vec_p_MC[0].SetMagThetaPhi(p_p_lab_MC[0],Theta_p_lab_MC[0],Phi_p_lab_MC[0]);
                    P_p_MC[0].SetVectM(vec_p_MC[0],m_p);

                    //histograms
                    hEkin_p_lab_MC[0]->Fill(Ekin_p_lab_MC[0]);
                    hp_p_lab_MC[0]->Fill(p_p_lab_MC[0]);
                    hE_p_lab_MC[0]->Fill(E_p_lab_MC[0]);
                    hTheta_p_lab_MC[0]->Fill(Theta_p_lab_MC[0]*TMath::RadToDeg());
                    hPhi_p_lab_MC[0]->Fill(Phi_p_lab_MC[0]*TMath::RadToDeg());
                    hEkin_vs_Theta_p_lab_MC[0]->Fill(Ekin_p_lab_MC[0],Theta_p_lab_MC[0]*TMath::RadToDeg());

                }

                //proton 2//
                if ((NrVertex == 1) && (PType == 14)&& (particleindex == 2)) {

                    Ekin_p_lab_MC[1] = part->GetEkin();
                    Theta_p_lab_MC[1] = part->GetTheta();
                    Phi_p_lab_MC[1] = part->GetPhi();

                    //cout<<"Ekin_p2_lab_MC = "<<Ekin_p_lab_MC[1]<<endl;
                    //cout<<"Theta_p2_lab_MC = "<<Theta_p_lab_MC[1]*TMath::RadToDeg()<<endl;
                    //cout<<"Phi_p2_lab_MC = "<<Phi_p_lab_MC[1]*TMath::RadToDeg()<<endl;
                    //cout<<"m_p2 = "<<part->GetMass();<<endl;

                    p_p_lab_MC[1] = TMath::Sqrt(Ekin_p_lab_MC[1]*(Ekin_p_lab_MC[1] + 2*m_p));           //E^2 = m^2 + p^2, E = Ekin + m => p^2 = Ekin^2 + 2*m*Ekin
                    //E_p_lab_MC[1] = TMath::Sqrt(p_p_lab_MC[1]*p_p_lab_MC[1] + m_p*m_p);               //total energy
                    E_p_lab_MC[1] = TMath::Sqrt(TMath::Power(p_p_lab_MC[1], 2) + TMath::Power(m_p, 2)); //total energy

                    vec_p_MC[1].SetMagThetaPhi(p_p_lab_MC[1],Theta_p_lab_MC[1],Phi_p_lab_MC[1]);
                    P_p_MC[1].SetVectM(vec_p_MC[1],m_p);

                    //histograms
                    hEkin_p_lab_MC[1]->Fill(Ekin_p_lab_MC[1]);
                    hp_p_lab_MC[1]->Fill(p_p_lab_MC[1]);
                    hE_p_lab_MC[1]->Fill(E_p_lab_MC[1]);
                    hTheta_p_lab_MC[1]->Fill(Theta_p_lab_MC[1]*TMath::RadToDeg());
                    hPhi_p_lab_MC[1]->Fill(Phi_p_lab_MC[1]*TMath::RadToDeg());
                    hEkin_vs_Theta_p_lab_MC[1]->Fill(Ekin_p_lab_MC[1],Theta_p_lab_MC[1]*TMath::RadToDeg());

                }

            }   //A03//

        }   //A02//

        ////Angles////
        Double_t OpeningAngle_p1_p2_lab_MC = (vec_p_MC[0].Angle(vec_p_MC[1]))*TMath::RadToDeg();    //angle between protons in CD & FD [deg]

        Double_t Delta_Phi_lab_MC = (Phi_p_lab_MC[1] - Phi_p_lab_MC[0])*TMath::RadToDeg();          //coplanarity definition

        hOpeningAngle_p1_p2_lab_MC->Fill(OpeningAngle_p1_p2_lab_MC);
        hDelta_Phi_lab_MC->Fill(Delta_Phi_lab_MC);

        //theta for protons
        hTheta_p1_vs_Theta_p2_lab_MC->Fill(Theta_p_lab_MC[0]*TMath::RadToDeg(),Theta_p_lab_MC[1]*TMath::RadToDeg());

        if ((Theta_p_lab_MC[0] >= 0.052) && (Theta_p_lab_MC[0] <= 0.314) && (Theta_p_lab_MC[1] >= 0.349) && (Theta_p_lab_MC[1] <= 2.950)) {
            hTheta_FDvsTheta_CD_lab_MC->Fill(Theta_p_lab_MC[0]*TMath::RadToDeg(), Theta_p_lab_MC[1]*TMath::RadToDeg());
        }

        if ((Theta_p_lab_MC[1] >= 0.052) && (Theta_p_lab_MC[1] <= 0.314) && (Theta_p_lab_MC[0] >= 0.349) && (Theta_p_lab_MC[0] <= 2.950)) {
            hTheta_FDvsTheta_CD_lab_MC->Fill(Theta_p_lab_MC[1]*TMath::RadToDeg(), Theta_p_lab_MC[0]*TMath::RadToDeg());
        }

        //BEAM KINETIC ENERGY//
        //FOUR-VECTORS//

        TVector3 vec_beam_MC = vec_n_MC + vec_p_MC[0] + vec_p_MC[1];

        beamMom = vec_beam_MC.Mag();

        hp_beam_MC->Fill(beamMom);

        TLorentzVector P_b_MC;      //4-vector of the beam
        P_b_MC.SetVectM(vec_beam_MC,m_beam);

        TVector3 vec_target_MC;
        vec_target_MC.SetMagThetaPhi(0.,0.,0.);

        TLorentzVector P_t_MC;      //4-vector of the target
        P_t_MC.SetVectM(vec_target_MC,m_target);

        TLorentzVector P_tot_MC = P_b_MC + P_t_MC;  //total 4-vector

        //proton-proton total energy
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

        TVector3 vec_n_cm_MC;
        TVector3 vec_p_cm_MC[2];

        Double_t p_n_cm_MC;
        Double_t E_n_cm_MC;
        Double_t Ekin_n_cm_MC;
        Double_t Theta_n_cm_MC;
        Double_t Phi_n_cm_MC;

        Double_t p_p_cm_MC[2];
        Double_t E_p_cm_MC[2];
        Double_t Ekin_p_cm_MC[2];
        Double_t Theta_p_cm_MC[2];
        Double_t Phi_p_cm_MC[2];

        ////boost to CM////
        TVector3 b_MC;
        b_MC = P_tot_MC.BoostVector();  //boost to LAB

        ////beam and target////
        P_b_MC.Boost(-b_MC);    //to CM
        P_t_MC.Boost(-b_MC);    //to CM

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

        ////proton 1////
        P_p_MC[0].Boost(-b_MC);    //to CM

        vec_p_cm_MC[0] = P_p_MC[0].Vect();
        p_p_cm_MC[0] = vec_p_cm_MC[0].Mag();
        //E_p_cm_MC[0] = TMath::Sqrt(p_p_cm_MC[0]*p_p_cm_MC[0] + m_p*m_p);
        E_p_cm_MC[0] = P_p_MC[0].E();
        Ekin_p_cm_MC[0] = E_p_cm_MC[0] - m_p;
        Theta_p_cm_MC[0] = P_p_MC[0].Theta();
        Phi_p_cm_MC[0] = P_p_MC[0].Phi();

        //histograms
        hp_p_cm_MC[0]->Fill(p_p_cm_MC[0]);
        hE_p_cm_MC[0]->Fill(E_p_cm_MC[0]);
        hEkin_p_cm_MC[0]->Fill(Ekin_p_cm_MC[0]);
        hTheta_p_cm_MC[0]->Fill(Theta_p_cm_MC[0]*TMath::RadToDeg());
        hPhi_p_cm_MC[0]->Fill(Phi_p_cm_MC[0]*TMath::RadToDeg());
        hEkin_vs_Theta_p_cm_MC[0]->Fill(Ekin_p_cm_MC[0],Theta_p_cm_MC[0]*TMath::RadToDeg());

        ////proton 2////
        P_p_MC[1].Boost(-b_MC);    //to CM

        vec_p_cm_MC[1] = P_p_MC[1].Vect();
        p_p_cm_MC[1] = vec_p_cm_MC[1].Mag();
        //E_p_cm_MC[1] = TMath::Sqrt(p_p_cm_MC[1]*p_p_cm_MC[1] + m_p*m_p);
        E_p_cm_MC[1] = P_p_MC[1].E();
        Ekin_p_cm_MC[1] = E_p_cm_MC[1] - m_p;
        Theta_p_cm_MC[1] = P_p_MC[1].Theta();
        Phi_p_cm_MC[1] = P_p_MC[1].Phi();

        //histograms
        hp_p_cm_MC[1]->Fill(p_p_cm_MC[1]);
        hE_p_cm_MC[1]->Fill(E_p_cm_MC[1]);
        hEkin_p_cm_MC[1]->Fill(Ekin_p_cm_MC[1]);
        hTheta_p_cm_MC[1]->Fill(Theta_p_cm_MC[1]*TMath::RadToDeg());
        hPhi_p_cm_MC[1]->Fill(Phi_p_cm_MC[1]*TMath::RadToDeg());
        hEkin_vs_Theta_p_cm_MC[1]->Fill(Ekin_p_cm_MC[1],Theta_p_cm_MC[1]*TMath::RadToDeg());

        ////
        Double_t OpeningAngle_p1_p2_cm_MC = (vec_p_cm_MC[0].Angle(vec_p_cm_MC[1]))*TMath::RadToDeg();   //angle between protons in CM [deg]

        Double_t Delta_Phi_cm_MC = (Phi_p_cm_MC[1] - Phi_p_cm_MC[0])*TMath::RadToDeg();  //coplanarity

        hOpeningAngle_p1_p2_cm_MC->Fill(OpeningAngle_p1_p2_cm_MC);
        hDelta_Phi_cm_MC->Fill(Delta_Phi_cm_MC);

        hTheta_p1_vs_Theta_p2_cm_MC->Fill(Theta_p_cm_MC[0]*TMath::RadToDeg(),Theta_p_cm_MC[1]*TMath::RadToDeg());

////////////////////////////////////////////////////////////////////////////////////////////////

        //cross section calculation for proton beam momentum in pp LAB frame (p_pp_beam)
        //and CM scattering angle (Theta_scatt) based on SAID data base

        Double_t p_pp_beam = 0.;  //effective beam momentum
        Double_t Theta_scatt = 0.;

        //beam momentum of first proton
        //(in the lab system of two protons=>assumption that first proton is moving and second is at rest)
        p_pp_beam = p_pp_beam_lab[counter];     //[GeV/c]

        //scattering angle in the CM
        Theta_scatt = Theta_scatt_cm[counter]*TMath::RadToDeg();    //[deg]

        //cout<<"p_pp_beam = "<<p_pp_beam<<endl;
        //cout<<"Theta_scatt = "<<Theta_scatt<<endl;

        hp_pp_beam->Fill(p_pp_beam);
        hTheta_scatt_cm->Fill(Theta_scatt);

////////////////////////////////////////////////////////////////////////////////////////////////

        //interpolation for above p_pp_beam and Theta_scatt
        Double_t XS_th1 = 0.;
        Double_t XS_th2 = 0.;

        Double_t XS_p1th1 = 0.;
        Double_t XS_p1th2 = 0.;
        Double_t XS_p2th1 = 0.;
        Double_t XS_p2th2 = 0.;

        //in TGraph
        Int_t Theta = 0;
        Int_t Theta_1 = 0;
        Int_t Theta_2 = 0;

        Double_t pbeam = 0.;
        Double_t pbeam_1 = 0.;
        Double_t pbeam_2 = 0.;

        Double_t t = 0.;
        Double_t u = 0.;

        for (Int_t Theta_s = 0; Theta_s < 181; Theta_s++) {

            Theta_1 = Theta_s;
            Theta_2 = Theta_s + 1;

            if (!((Theta_scatt >= Theta_1) && (Theta_scatt < Theta_2))) continue;

            XS_th1 = gCrossSection[Theta_1]->Eval(p_pp_beam);
            XS_th2 = gCrossSection[Theta_2]->Eval(p_pp_beam);

            if ((TMath::Abs(Theta_scatt - Theta_1)) <= (TMath::Abs(Theta_scatt - Theta_2))) {
                Theta = Theta_1;
            } else {
                Theta = Theta_2;
            }

            //cout<<"Theta_scatt = "<<Theta_scatt<<endl;
            //cout<<"Theta_1 = "<<Theta_1<<endl;
            //cout<<"Theta_2 = "<<Theta_2<<endl;
            //cout<<"Theta = "<<Theta<<endl;

            for (Int_t pbeam_s = 700; pbeam_s < 3050; pbeam_s += 50) {
                pbeam_1 = pbeam_s/1000.;
                pbeam_2 = (pbeam_s + 50.)/1000.;

                if (!((p_pp_beam >= pbeam_1) && (p_pp_beam < pbeam_2))) continue;

                if ((TMath::Abs(p_pp_beam - pbeam_1)) <= (TMath::Abs(p_pp_beam - pbeam_2))) {
                    pbeam = pbeam_1;
                } else {
                    pbeam = pbeam_2;
                }

                //cout<<"p_pp_beam = "<<p_pp_beam<<endl;
                //cout<<"pbeam_1 = "<<pbeam_1<<endl;
                //cout<<"pbeam_2 = "<<pbeam_2<<endl;
                //cout<<"pbeam = "<<pbeam<<endl;

                XS_p1th1 = gCrossSection[Theta_1]->Eval(pbeam_1);
                XS_p2th1 = gCrossSection[Theta_1]->Eval(pbeam_2);
                XS_p1th2 = gCrossSection[Theta_2]->Eval(pbeam_1);
                XS_p2th2 = gCrossSection[Theta_2]->Eval(pbeam_2);

                //cout<<"XS_p1th1 = "<<XS_p1th1<<endl;
                //cout<<"XS_p2th1 = "<<XS_p2th1<<endl;
                //cout<<"XS_p1th2 = "<<XS_p1th2<<endl;
                //cout<<"XS_p2th2 = "<<XS_p2th2<<endl;

                //bilinear interpolation
                t = p_pp_beam - pbeam_1;
                u = Theta_scatt - Theta_1;

                XS_bilin = (1-t)*(1-u)*XS_p1th1 + t*(1-u)*XS_p2th1 + t*u*XS_p2th2+(1-t)*u*XS_p1th2;

                //proximate points
                XS_prx = gCrossSection[Theta]->Eval(pbeam);

                //linear extrapolation
                if (XS_th2 > XS_th1) {
                    XS_lin = ((XS_th2 - XS_th1)/(Theta_2 - Theta_1))*(Theta_scatt - Theta_1) + XS_th1;
                } else {
                    XS_lin = ((XS_th1 - XS_th2)/(Theta_2 - Theta_1))*(Theta_2 - Theta_scatt) + XS_th2;
                }
            }
        }

        XS_sq = TMath::Power(XS_bilin, 2);

        //cout<<"XS_bilin = "<<XS_bilin<<endl;
        //cout<<"XS_lin = "<<XS_linr<<endl;
        //cout<<"XS_prx = "<<XS_prx<<endl;
        //cout<<"XS_sq = "<<XS_sq<<endl;

        hCrossSection_bilin->Fill(XS_bilin);
        hCrossSection_prx->Fill(XS_prx);
        hCrossSection_lin->Fill(XS_lin);

        hGenerated_Q_wght_bilin->Fill(Q, XS_bilin);
        hGenerated_Q_wght_prx->Fill(Q, XS_prx);
        hGenerated_Q_wght_lin->Fill(Q, XS_lin);
        hGenerated_Q_wght_sq->Fill(Q, XS_sq);

    }   //A01//

//////////////////////////////////////RECONSTRUCTED EVENTS//////////////////////////////////////

    ///////LEVEL 0: TRIGGER///////

    //if (!(fHeader->TriggerNumSet(21)))  return; //h[0]: level 0 (trigger #21)

    hStatistics[0]->Fill(0);
/*
    //BEAM MOMENTUM//
    Double_t t_incycle = 1000*fHeader->GetTimeInCycle();        //time in cycle (t_axis - t_start)
    beamMom = (gBeamMomentum->Eval(t_incycle)) + pbeam_offset;  //beam momentum [GeV/c]

    Double_t sqMass = TMath::Power(m_beam, 2) + TMath::Power(m_target, 2);
    Double_t beamEnergy = TMath::Sqrt(TMath::Power(beamMom, 2) + TMath::Power(m_beam, 2));
    Double_t s = TMath::Sqrt(sqMass + 2*m_target*beamEnergy);  //corresponds to Q from Q=-70MeV to Q=30MeV
    //Double_t s = TMath::Sqrt(m_beam*m_beam + m_target*m_target +2*m_target*TMath::Sqrt(m_beam*m_beam + beamMom*beamMom)); //corresponds to Q from Q=-70MeV to Q=30MeV

    Q = 1000*(s - s_thr);   //excess energy [MeV]
*/
    hp_beam[0][0]->Fill(beamMom);

    hQ[0][0]->Fill(Q);
    hQ_wght_bilin[0][0]->Fill(Q, XS_bilin);
    hQ_wght_prx[0][0]->Fill(Q, XS_prx);
    hQ_wght_lin[0][0]->Fill(Q, XS_lin);
    hQ_wght_sq[0][0]->Fill(Q, XS_sq);

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

        hEdepPSBvsSEC[0][0]->Fill((CDTrack->Edep(151,174)),(CDTrack->GetSpecELossPS()*E_calib));
        hEdepPSBvsSigMom[0][0]->Fill((CDTrack->Momentum()*CDTrack->Charge()),(CDTrack->GetSpecELossPS()*E_calib));
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
            EdepPSB = TrackCD->GetSpecELossPS()*E_calib;
            SgnMom = TrackCD->Momentum()*(TrackCD->Charge());

        }

        //coplanarity definition
        //correction factor
        Double_t Delta_Phi_offset = 0.0392699;  //for MC
        //Double_t Delta_Phi_offset = 0.0436332;  //for DATA

        Double_t Delta_Phi = (PhiFD_lab - PhiCD_lab + Delta_Phi_offset)*TMath::RadToDeg();
        Double_t Delta_Phi_abs = (TMath::Abs(PhiFD_lab - PhiCD_lab))*TMath::RadToDeg();
        Double_t Delta_Phi_sym = (fmod((2*TMath::Pi() + (PhiFD_lab - PhiCD_lab + Delta_Phi_offset)), (2*TMath::Pi())))*TMath::RadToDeg();

/////////////////////////////////////////LEVELS 1-5(6)//////////////////////////////////////////

        Bool_t lev[7][10];

        //main level
        lev[1][0] = (ThetaFD_lab != 0.125);

        //cut on scattering angle in FD (3,18) [deg]
        lev[2][0] = (lev[1][0] && (ThetaFD_lab >= 0.052) && (ThetaFD_lab <= 0.314));

        //positively charged particles in CD
        lev[3][0] = (lev[2][0] &&  (SgnMom >= 0.));

        //graphical cut on spectrum of energy loss in the PSB vs. energy deposited in the SEC
        lev[4][0] = (lev[3][0] && (EdepPSB >= -0.00618182*EdepSEC + 0.003400));

        //cut on scattering angle in CD (40,100) [deg]
        lev[5][0] = (lev[4][0] && (ThetaCD_lab >= 0.698) && (ThetaCD_lab <= 1.745));

        //without condion #3
        lev[6][0] = (lev[2][0] && (EdepPSB >= -0.00618182*EdepSEC + 0.003400) && (ThetaCD_lab >= 0.698) && (ThetaCD_lab <= 1.745));

        //loop
        for (Int_t l = 1; l < 7; l++) {     //C01//

            if (lev[l][0]) {            //C02//

                if (l < 6) { hStatistics[0]->Fill(l); }

                hQ[l][0]->Fill(Q);

                hQ_wght_bilin[l][0]->Fill(Q, XS_bilin);
                hQ_wght_prx[l][0]->Fill(Q, XS_prx);
                hQ_wght_lin[l][0]->Fill(Q, XS_lin);
                hQ_wght_sq[l][0]->Fill(Q, XS_sq);

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

                for(Int_t i = 1; i < 41; i++) {
                    Double_t Q_min = -70.0 + (i-1)*2.5;
                    Double_t Q_max = -67.5 + (i-1)*2.5;

                    if((Q >= Q_min) && (Q < Q_max)) {
                        hDelta_Phi_sym_Q[l][0][i]->Fill(Delta_Phi_sym);
                        hTheta_CDC_Q[l][0][i]->Fill(ThetaCD_lab*TMath::RadToDeg());

                        for(Int_t k = 1; k < 16; k++) {
                            Double_t ThetaFDC_min_i = 3.0 + (k-1);
                            Double_t ThetaFDC_max_i = 4.0 + (k-1);

                            if((ThetaFD_lab*TMath::RadToDeg()>=ThetaFDC_min_i) && (ThetaFD_lab*TMath::RadToDeg()<ThetaFDC_max_i)) {
                                hTheta_CDC_Q_Th[l][0][i][k]->Fill(ThetaCD_lab*TMath::RadToDeg());
                            }
                        }
                    }
                }

                for(Int_t k = 1; k < 16; k++) {
                    Double_t ThetaFDC_min = 3.0 + (k-1);
                    Double_t ThetaFDC_max = 4.0 + (k-1);

                    if((ThetaFD_lab*TMath::RadToDeg()>=ThetaFDC_min) && (ThetaFD_lab*TMath::RadToDeg()<ThetaFDC_max)) {
                        hTheta_CDC_Th[l][0][k]->Fill(ThetaCD_lab*TMath::RadToDeg());
                    }
                }

            }   //C02//

        }       //C01//

//////////////////////////////////////////SYSTEMATICS///////////////////////////////////////////

        lev[5][1] = (lev[3][0] && (EdepPSB >= -0.00618182*EdepSEC + 0.002980) && (ThetaCD_lab >= 0.698) && (ThetaCD_lab <= 1.745)); //graphical cut (+10%)
        lev[5][2] = (lev[3][0] && (EdepPSB >= -0.00618182*EdepSEC + 0.003922) && (ThetaCD_lab >= 0.698) && (ThetaCD_lab <= 1.745)); //graphical cut (-10%)

        lev[5][3] = (lev[3][0] && (EdepPSB >= -0.00618182*EdepSEC + 0.003400) && (ThetaCD_lab >= 0.646) && (ThetaCD_lab <= 1.798)); //CD scattering angle cut (37,103) [deg]
        lev[5][4] = (lev[3][0] && (EdepPSB >= -0.00618182*EdepSEC + 0.003400) && (ThetaCD_lab >= 0.751) && (ThetaCD_lab <= 1.693)); //CD scattering angle cut (43,97) [deg]

        lev[5][5] = (lev[3][0] && (EdepPSB >= -0.00618182*EdepSEC + 0.003400) && (ThetaCD_lab >= 0.594) && (ThetaCD_lab <= 1.745)); //CD scattering angle cut (34,100) [deg]
        lev[5][6] = (lev[3][0] && (EdepPSB >= -0.00618182*EdepSEC + 0.003400) && (ThetaCD_lab >= 0.802) && (ThetaCD_lab <= 1.745)); //CD scattering angle cut (46,100) [deg]
        lev[5][7] = (lev[3][0] && (EdepPSB >= -0.00618182*EdepSEC + 0.003400) && (ThetaCD_lab >= 0.698) && (ThetaCD_lab <= 1.850)); //CD scattering angle cut (40,106) [deg]
        lev[5][8] = (lev[3][0] && (EdepPSB >= -0.00618182*EdepSEC + 0.003400) && (ThetaCD_lab >= 0.698) && (ThetaCD_lab <= 1.641)); //CD scattering angle cut (40,94) [deg]

        //loop
        for (Int_t c = 1; c < 9; c++) { //D01//

            if (lev[5][c]) {            //D02//

                hQ[5][c]->Fill(Q);

                hQ_wght_bilin[5][c]->Fill(Q, XS_bilin);
                hQ_wght_prx[5][c]->Fill(Q, XS_prx);
                hQ_wght_lin[5][c]->Fill(Q, XS_lin);
                hQ_wght_sq[5][c]->Fill(Q, XS_sq);

                hTime_FDC[5][c]->Fill(TimeFD);
                hTheta_FDC[5][c]->Fill(ThetaFD_lab*TMath::RadToDeg());
                hPhi_FDC[5][c]->Fill(PhiFD_lab*TMath::RadToDeg());

                hEdepFWC1vsFRH1[5][c]->Fill(EdepFWC1,EdepFRH1);
                hEdepFWC2vsFRH1[5][c]->Fill(EdepFWC2,EdepFRH1);
                hEdepFTH1vsFRH1[5][c]->Fill(EdepFTH1,EdepFRH1);
                hEdepFRH1vsFRH2[5][c]->Fill(EdepFRH1,EdepFRH2);
                hEdepFRH2vsFRH3[5][c]->Fill(EdepFRH2,EdepFRH3);
                hEdepFWC1vsFRH1FRH2FRH3[5][c]->Fill(EdepFWC1,(EdepFRH1 + EdepFRH2 + EdepFRH3));

                hTime_CDC[5][c]->Fill(TimeCD);
                hMom_CDC[5][c]->Fill(MomCD_lab);
                hTheta_CDC[5][c]->Fill(ThetaCD_lab*TMath::RadToDeg());
                hPhi_CDC[5][c]->Fill(PhiCD_lab*TMath::RadToDeg());

                hEdepPSBvsSEC[5][c]->Fill(EdepSEC,EdepPSB);
                hEdepPSBvsSigMom[5][c]->Fill(SgnMom,EdepPSB);
                hEdepSECvsSigMom[5][c]->Fill(SgnMom,EdepSEC);

                hTheta_FDvsTheta_CD[5][c]->Fill(ThetaFD_lab*TMath::RadToDeg(),ThetaCD_lab*TMath::RadToDeg());
                hDeltaTime[5][c]->Fill(TimeFD-TimeCD);

                hDelta_Phi[5][c]->Fill(Delta_Phi);
                hDelta_Phi_abs[5][c]->Fill(Delta_Phi_abs);
                hDelta_Phi_sym[5][c]->Fill(Delta_Phi_sym);

                for(Int_t i = 1; i < 41; i++) {
                    Double_t Q_min = -70.0 + (i-1)*2.5;
                    Double_t Q_max = -67.5 + (i-1)*2.5;

                    if((Q >= Q_min) && (Q < Q_max)) {
                        hDelta_Phi_sym_Q[5][c][i]->Fill(Delta_Phi_sym);
                        hTheta_CDC_Q[5][c][i]->Fill(ThetaCD_lab*TMath::RadToDeg());

                        for(Int_t k = 1; k < 16; k++) {
                            Double_t ThetaFDC_min_i = 3.0 + (k-1);
                            Double_t ThetaFDC_max_i = 4.0 + (k-1);

                            if((ThetaFD_lab*TMath::RadToDeg()>=ThetaFDC_min_i) && (ThetaFD_lab*TMath::RadToDeg()<ThetaFDC_max_i)) {
                                hTheta_CDC_Q_Th[5][c][i][k]->Fill(ThetaCD_lab*TMath::RadToDeg());
                            }
                        }
                    }
                }

                for(Int_t k = 1; k < 16; k++) {
                    Double_t ThetaFDC_min = 3.0 + (k-1);
                    Double_t ThetaFDC_max = 4.0 + (k-1);

                    if((ThetaFD_lab*TMath::RadToDeg()>=ThetaFDC_min) && (ThetaFD_lab*TMath::RadToDeg()<ThetaFDC_max)) {
                        hTheta_CDC_Th[5][c][k]->Fill(ThetaCD_lab*TMath::RadToDeg());
                    }
                }

            }   //D02//

        }       //D01//

    }           //B01//

    return;

}   //01//

////////////////////////////////////////////////////////////////////////////////////////////////

void eventselection::SetupSpectra(const char * lpath) {   //02//

    TString h_mc = Form("WMC");

    TString g[7][10];
    g[0][0] = Form("DATA_lev0");
    g[1][0] = Form("DATA_lev1");
    g[2][0] = Form("DATA_lev2");
    g[3][0] = Form("DATA_lev3");
    g[4][0] = Form("DATA_lev4");
    g[5][0] = Form("DATA_lev5");
    g[6][0] = Form("DATA_lev6");

    g[5][1] = Form("DATA_syst_lev5_E1");
    g[5][2] = Form("DATA_syst_lev5_E2");
    g[5][3] = Form("DATA_syst_lev5_T1");
    g[5][4] = Form("DATA_syst_lev5_T2");
    g[5][5] = Form("DATA_syst_lev5_T3");
    g[5][6] = Form("DATA_syst_lev5_T4");
    g[5][7] = Form("DATA_syst_lev5_T5");
    g[5][8] = Form("DATA_syst_lev5_T6");

    TString h_st = Form("Statistics");

////////////////////////////////////TRUE EVENTS (MONTE CARLO)///////////////////////////////////

    hp_pp_beam = new TH1F("hp_pp_beam","",1000,0.,3.);
    hp_pp_beam->GetXaxis()->SetTitle("p^{pp}_{beam} [GeV/c]");
    hp_pp_beam->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hp_pp_beam,h_mc);

    hTheta_scatt_cm = new TH1F("hTheta_scatt_cm","",360,0.,180.);
    hTheta_scatt_cm->GetXaxis()->SetTitle("#theta^{cm}_{scatt} [deg]");
    hTheta_scatt_cm->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hTheta_scatt_cm,h_mc);

    hp_beam_MC = new TH1F("hp_beam_MC","",200,1.41,1.65);
    hp_beam_MC->GetXaxis()->SetTitle("p_{beam} [GeV/c]");
    hp_beam_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hp_beam_MC,h_mc);

    hGenerated_Q = new TH1F("hGenerated_Q","",40,-70.,30.);
    hGenerated_Q->GetXaxis()->SetTitle("Q [MeV]");
    hGenerated_Q->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hGenerated_Q,h_mc);

    hGenerated_Q_wght_bilin =new TH1F("hGenerated_Q_wght_bilin","",40,-70.,30.);
    hGenerated_Q_wght_bilin->GetXaxis()->SetTitle("Q [MeV]");
    hGenerated_Q_wght_bilin->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hGenerated_Q_wght_bilin,h_mc);

    hGenerated_Q_wght_prx =new TH1F("hGenerated_Q_wght_prx","",40,-70.,30.);
    hGenerated_Q_wght_prx->GetXaxis()->SetTitle("Q [MeV]");
    hGenerated_Q_wght_prx->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hGenerated_Q_wght_prx,h_mc);

    hGenerated_Q_wght_lin =new TH1F("hGenerated_Q_wght_lin","",40,-70.,30.);
    hGenerated_Q_wght_lin->GetXaxis()->SetTitle("Q [MeV]");
    hGenerated_Q_wght_lin->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hGenerated_Q_wght_lin,h_mc);

    hGenerated_Q_wght_sq =new TH1F("hGenerated_Q_wght_sq","",40,-70.,30.);
    hGenerated_Q_wght_sq->GetXaxis()->SetTitle("Q [MeV]");
    hGenerated_Q_wght_sq->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hGenerated_Q_wght_sq,h_mc);

    //neutron//
    hEkin_n_lab_MC = new TH1F("hEkin_n_lab_MC","",500,0.,0.1);
    hEkin_n_lab_MC->GetXaxis()->SetTitle("E^{lab}_{kin} [GeV]");
    hEkin_n_lab_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hEkin_n_lab_MC,h_mc);

    hEkin_n_cm_MC = new TH1F("hEkin_n_cm_MC","",500,0.,0.4);
    hEkin_n_cm_MC->GetXaxis()->SetTitle("E^{cm}_{kin} [GeV]");
    hEkin_n_cm_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hEkin_n_cm_MC,h_mc);

    hE_n_lab_MC = new TH1F("hE_n_lab_MC","",500,0.9,1.04);
    hE_n_lab_MC->GetXaxis()->SetTitle("E^{lab} [GeV]");
    hE_n_lab_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hE_n_lab_MC,h_mc);

    hE_n_cm_MC = new TH1F("hE_n_cm_MC","",500,0.9,1.3);
    hE_n_cm_MC->GetXaxis()->SetTitle("E^{cm} [GeV]");
    hE_n_cm_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hE_n_cm_MC,h_mc);

    hp_n_lab_MC = new TH1F("hp_n_lab_MC","",500,0.,0.4);
    hp_n_lab_MC->GetXaxis()->SetTitle("p^{lab} [GeV/c]");
    hp_n_lab_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hp_n_lab_MC,h_mc);

    hp_n_cm_MC = new TH1F("hp_n_cm_MC","",500,0.,1.);
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

    hEkin_vs_Theta_n_lab_MC = new TH2F("hEkin_vs_Theta_n_lab_MC","",500,0.,0.08,500,0.,180.);
    hEkin_vs_Theta_n_lab_MC->GetXaxis()->SetTitle("E^{lab}_{kin} [GeV]");
    hEkin_vs_Theta_n_lab_MC->GetYaxis()->SetTitle("#theta^{lab} [deg]");
    gHistoManager->Add(hEkin_vs_Theta_n_lab_MC,h_mc);

    hEkin_vs_Theta_n_cm_MC = new TH2F("hEkin_vs_Theta_n_cm_MC","",500,0.,0.4,500,0.,180.);
    hEkin_vs_Theta_n_cm_MC->GetXaxis()->SetTitle("E^{cm}_{kin} [GeV]");
    hEkin_vs_Theta_n_cm_MC->GetYaxis()->SetTitle("#theta^{cm} [deg]");
    gHistoManager->Add(hEkin_vs_Theta_n_cm_MC,h_mc);

    //protons//
    for (Int_t j = 0; j < 2; j++) {
        TString hname = Form("hEkin_p%d_lab_MC",j+1);
        hEkin_p_lab_MC[j] = new TH1F(hname,"",500,0.,1.);
        hEkin_p_lab_MC[j]->GetXaxis()->SetTitle("E^{lab}_{kin} [GeV]");
        hEkin_p_lab_MC[j]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hEkin_p_lab_MC[j],h_mc);
    }

    for (Int_t j = 0; j < 2; j++) {
        TString hname = Form("hEkin_p%d_cm_MC",j+1);
        hEkin_p_cm_MC[j] = new TH1F(hname,"",500,0.,0.5);
        hEkin_p_cm_MC[j]->GetXaxis()->SetTitle("E^{cm}_{kin} [GeV]");
        hEkin_p_cm_MC[j]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hEkin_p_cm_MC[j],h_mc);
    }

    for (Int_t j = 0; j < 2; j++) {
        TString hname = Form("hE_p%d_lab_MC",j+1);
        hE_p_lab_MC[j] = new TH1F(hname,"",500,0.9,2.);
        hE_p_lab_MC[j]->GetXaxis()->SetTitle("E^{lab} [GeV]");
        hE_p_lab_MC[j]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hE_p_lab_MC[j],h_mc);
    }

    for (Int_t j = 0; j < 2; j++) {
        TString hname = Form("hE_p%d_cm_MC",j+1);
        hE_p_cm_MC[j] = new TH1F(hname,"",500,0.9,1.4);
        hE_p_cm_MC[j]->GetXaxis()->SetTitle("E^{cm} [GeV]");
        hE_p_cm_MC[j]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hE_p_cm_MC[j],h_mc);
    }

    for (Int_t j = 0; j < 2; j++) {
        TString hname = Form("hp_p%d_lab_MC",j+1);
        hp_p_lab_MC[j] = new TH1F(hname,"",500,0.,1.8);
        hp_p_lab_MC[j]->GetXaxis()->SetTitle("p^{lab} [GeV/c]");
        hp_p_lab_MC[j]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hp_p_lab_MC[j],h_mc);
    }

    for (Int_t j = 0; j < 2; j++) {
        TString hname = Form("hp_p%d_cm_MC",j+1);
        hp_p_cm_MC[j] = new TH1F(hname,"",500,0.,1.2);
        hp_p_cm_MC[j]->GetXaxis()->SetTitle("p^{cm} [GeV/c]");
        hp_p_cm_MC[j]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hp_p_cm_MC[j],h_mc);
    }

    for (Int_t j = 0; j < 2; j++) {
        TString hname = Form("hTheta_p%d_lab_MC",j+1);
        hTheta_p_lab_MC[j] = new TH1F(hname,"",360,0.,180.);
        hTheta_p_lab_MC[j]->GetXaxis()->SetTitle("#theta^{lab} [deg]");
        hTheta_p_lab_MC[j]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hTheta_p_lab_MC[j],h_mc);
    }

    for (Int_t j = 0; j < 2; j++) {
        TString hname = Form("hTheta_p%d_cm_MC",j+1);
        hTheta_p_cm_MC[j] = new TH1F(hname,"",360,0.,180.);
        hTheta_p_cm_MC[j]->GetXaxis()->SetTitle("#theta^{cm} [deg]");
        hTheta_p_cm_MC[j]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hTheta_p_cm_MC[j],h_mc);
    }

    for (Int_t j = 0; j < 2; j++) {
        TString hname = Form("hPhi_p%d_lab_MC",j+1);
        hPhi_p_lab_MC[j] = new TH1F(hname,"",360,0.,360.);
        hPhi_p_lab_MC[j]->GetXaxis()->SetTitle("#phi^{lab} [deg]");
        hPhi_p_lab_MC[j]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hPhi_p_lab_MC[j],h_mc);
    }

    for (Int_t j = 0; j < 2; j++) {
        TString hname = Form("hPhi_p%d_cm_MC",j+1);
        hPhi_p_cm_MC[j] = new TH1F(hname,"",360,-180.,180.);
        hPhi_p_cm_MC[j]->GetXaxis()->SetTitle("#phi^{cm} [deg]");
        hPhi_p_cm_MC[j]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hPhi_p_cm_MC[j],h_mc);
    }

    for (Int_t j = 0; j < 2; j++) {
        TString hname = Form("hEkin_vs_Theta_p%d_lab_MC",j+1);
        hEkin_vs_Theta_p_lab_MC[j] = new TH2F(hname,"",500,0.,1.,500,0.,180.);
        hEkin_vs_Theta_p_lab_MC[j]->GetXaxis()->SetTitle("E^{lab}_{kin} [GeV]");
        hEkin_vs_Theta_p_lab_MC[j]->GetYaxis()->SetTitle("#theta^{lab} [deg]");
        gHistoManager->Add(hEkin_vs_Theta_p_lab_MC[j],h_mc);
    }

    for (Int_t j = 0; j < 2; j++) {
        TString hname = Form("hEkin_vs_Theta_p%d_cm_MC",j+1);
        hEkin_vs_Theta_p_cm_MC[j] = new TH2F(hname,"",500,0.,0.5,500,0.,180.);
        hEkin_vs_Theta_p_cm_MC[j]->GetXaxis()->SetTitle("E^{cm}_{kin} [GeV]");
        hEkin_vs_Theta_p_cm_MC[j]->GetYaxis()->SetTitle("#theta^{cm} [deg]");
        gHistoManager->Add(hEkin_vs_Theta_p_cm_MC[j],h_mc);
    }

    ////
    hOpeningAngle_p1_p2_lab_MC=new TH1F("hOpeningAngle_p1_p2_lab_MC","",360,0.,180.);
    hOpeningAngle_p1_p2_lab_MC->GetXaxis()->SetTitle("#theta^{lab}_{p_{1},p_{2}} [deg]");
    hOpeningAngle_p1_p2_lab_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hOpeningAngle_p1_p2_lab_MC,h_mc);

    hOpeningAngle_p1_p2_cm_MC=new TH1F("hOpeningAngle_p1_p2_cm_MC","",360,0.,180.);
    hOpeningAngle_p1_p2_cm_MC->GetXaxis()->SetTitle("#theta^{cm}_{p_{1},p_{2}} [deg]");
    hOpeningAngle_p1_p2_cm_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hOpeningAngle_p1_p2_cm_MC,h_mc);

    hTheta_FDvsTheta_CD_lab_MC = new TH2F("hTheta_FDvsTheta_CD_lab_MC","",540,0.,20.,540,0.,180.);
    hTheta_FDvsTheta_CD_lab_MC->GetXaxis()->SetTitle("#theta^{lab}_{FD} [deg]");
    hTheta_FDvsTheta_CD_lab_MC->GetYaxis()->SetTitle("#theta^{lab}_{CD} [deg]");
    gHistoManager->Add(hTheta_FDvsTheta_CD_lab_MC,h_mc);

    hTheta_p1_vs_Theta_p2_lab_MC = new TH2F("hTheta_p1_vs_Theta_p2_lab_MC","",540,0.,180.,540,0.,180.);
    hTheta_p1_vs_Theta_p2_lab_MC->GetXaxis()->SetTitle("#theta^{lab}_{p_{1}} [deg]");
    hTheta_p1_vs_Theta_p2_lab_MC->GetYaxis()->SetTitle("#theta^{lab}_{p_{2}} [deg]");
    gHistoManager->Add(hTheta_p1_vs_Theta_p2_lab_MC,h_mc);

    hTheta_p1_vs_Theta_p2_cm_MC = new TH2F("hTheta_p1_vs_Theta_p2_cm_MC","",540,0.,180.,540,0.,180.);
    hTheta_p1_vs_Theta_p2_cm_MC->GetXaxis()->SetTitle("#theta^{cm}_{p_{1}} [deg]");
    hTheta_p1_vs_Theta_p2_cm_MC->GetYaxis()->SetTitle("#theta^{cm}_{p_{2}} [deg]");
    gHistoManager->Add(hTheta_p1_vs_Theta_p2_cm_MC,h_mc);

    ////
    hDelta_Phi_lab_MC = new TH1F("hDelta_Phi_lab_MC","",360,0.,360.);
    hDelta_Phi_lab_MC->GetXaxis()->SetTitle("#Delta#phi^{lab} [deg]");
    hDelta_Phi_lab_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hDelta_Phi_lab_MC,h_mc);

    hDelta_Phi_cm_MC = new TH1F("hDelta_Phi_cm_MC","",360,0.,360.);
    hDelta_Phi_cm_MC->GetXaxis()->SetTitle("#Delta#phi^{cm} [deg]");
    hDelta_Phi_cm_MC->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hDelta_Phi_cm_MC,h_mc);

    ////
    hCrossSection_bilin = new TH1F("hCrossSection_bilin","Bilinear",500,0.,30.);
    hCrossSection_bilin->GetXaxis()->SetTitle("#sigma_{pp} [nb]");
    hCrossSection_bilin->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hCrossSection_bilin,h_mc);

    hCrossSection_prx = new TH1F("hCrossSection_prx","Approximation",500,0.,30.);
    hCrossSection_prx->GetXaxis()->SetTitle("#sigma_{pp} [nb]");
    hCrossSection_prx->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hCrossSection_prx,h_mc);

    hCrossSection_lin = new TH1F("hCrossSection_lin","Linear",500,0.,30.);
    hCrossSection_lin->GetXaxis()->SetTitle("#sigma_{pp} [nb]");
    hCrossSection_lin->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hCrossSection_lin,h_mc);

//////////////////////////////////////////RECONSTRUCTED/////////////////////////////////////////

    ////DATA level 0: trigger #21 (g[0][0])////
    hp_beam[0][0] = new TH1F("hp_beam_lev0","",200,1.41,1.65);
    hp_beam[0][0]->GetXaxis()->SetTitle("p_{beam} [GeV/c]");
    hp_beam[0][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hp_beam[0][0],g[0][0]);

    hQ[0][0] = new TH1F("hQ_lev0","",40,-70.,30.);
    hQ[0][0]->GetXaxis()->SetTitle("Q [MeV]");
    hQ[0][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hQ[0][0],g[0][0]);

    hQ_wght_bilin[0][0] = new TH1F("hQ_wght_bilin_lev0","",40,-70.,30.);
    hQ_wght_bilin[0][0]->GetXaxis()->SetTitle("Q [MeV]");
    hQ_wght_bilin[0][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hQ_wght_bilin[0][0],g[0][0]);

    hQ_wght_prx[0][0] = new TH1F("hQ_wght_prx_lev0","",40,-70.,30.);
    hQ_wght_prx[0][0]->GetXaxis()->SetTitle("Q [MeV]");
    hQ_wght_prx[0][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hQ_wght_prx[0][0],g[0][0]);

    hQ_wght_lin[0][0] = new TH1F("hQ_wght_lin_lev0","",40,-70.,30.);
    hQ_wght_lin[0][0]->GetXaxis()->SetTitle("Q [MeV]");
    hQ_wght_lin[0][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hQ_wght_lin[0][0],g[0][0]);

    hQ_wght_sq[0][0] = new TH1F("hQ_wght_sq_lev0","",40,-70.,30.);
    hQ_wght_sq[0][0]->GetXaxis()->SetTitle("Q [MeV]");
    hQ_wght_sq[0][0]->GetYaxis()->SetTitle("counts");
    gHistoManager->Add(hQ_wght_sq[0][0],g[0][0]);

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

////DATA: LEVELS 1-5(6)////

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

    for (Int_t l = 1; l < 7; l++) {
        TString hname = Form("hQ_lev%d",l);
        hQ[l][0] = new TH1F(hname,"",40,-70.,30.);
        hQ[l][0]->GetXaxis()->SetTitle("Q [MeV]");
        hQ[l][0]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hQ[l][0],g[l][0]);
    }

    for (Int_t l = 1; l < 7; l++) {
        TString hname = Form("hQ_wght_bilin_lev%d",l);
        hQ_wght_bilin[l][0] = new TH1F(hname,"",40,-70.,30.);
        hQ_wght_bilin[l][0]->GetXaxis()->SetTitle("Q [MeV]");
        hQ_wght_bilin[l][0]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hQ_wght_bilin[l][0],g[l][0]);
    }

    for (Int_t l = 1; l < 7; l++) {
        TString hname = Form("hQ_wght_prx_lev%d",l);
        hQ_wght_prx[l][0] = new TH1F(hname,"",40,-70.,30.);
        hQ_wght_prx[l][0]->GetXaxis()->SetTitle("Q [MeV]");
        hQ_wght_prx[l][0]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hQ_wght_prx[l][0],g[l][0]);
    }

    for (Int_t l = 1; l < 7; l++) {
        TString hname = Form("hQ_wght_lin_lev%d",l);
        hQ_wght_lin[l][0] = new TH1F(hname,"",40,-70.,30.);
        hQ_wght_lin[l][0]->GetXaxis()->SetTitle("Q [MeV]");
        hQ_wght_lin[l][0]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hQ_wght_lin[l][0],g[l][0]);
    }

    for (Int_t l = 1; l < 7; l++) {
        TString hname = Form("hQ_wght_sq_lev%d",l);
        hQ_wght_sq[l][0] = new TH1F(hname,"",40,-70.,30.);
        hQ_wght_sq[l][0]->GetXaxis()->SetTitle("Q [MeV]");
        hQ_wght_sq[l][0]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hQ_wght_sq[l][0],g[l][0]);
    }

    //FD//
    for (Int_t l = 1; l < 7; l++) {
        TString hname = Form("hEdepFWC1vsFRH1_lev%d",l);
        hEdepFWC1vsFRH1[l][0] =new TH2F(hname,"",1000,0.,0.01,1000,0.,0.5);
        hEdepFWC1vsFRH1[l][0]->GetXaxis()->SetTitle("Edep(FWC1) [GeV]");
        hEdepFWC1vsFRH1[l][0]->GetYaxis()->SetTitle("Edep(FRH1) [GeV]");
        gHistoManager->Add(hEdepFWC1vsFRH1[l][0],g[l][0]);
    }

    for (Int_t l = 1; l < 7; l++) {
        TString hname = Form("hEdepFWC2vsFRH1_lev%d",l);
        hEdepFWC2vsFRH1[l][0] =new TH2F(hname,"",1000,0.,0.01,1000,0.,0.5);
        hEdepFWC2vsFRH1[l][0]->GetXaxis()->SetTitle("Edep(FWC2) [GeV]");
        hEdepFWC2vsFRH1[l][0]->GetYaxis()->SetTitle("Edep(FRH1) [GeV]");
        gHistoManager->Add(hEdepFWC2vsFRH1[l][0],g[l][0]);
    }

    for (Int_t l = 1; l < 7; l++) {
        TString hname = Form("hEdepFTH1vsFRH1_lev%d",l);
        hEdepFTH1vsFRH1[l][0] =new TH2F(hname,"",1000,0.,0.02,1000,0.,0.5);
        hEdepFTH1vsFRH1[l][0]->GetXaxis()->SetTitle("Edep(FTH1) [GeV]");
        hEdepFTH1vsFRH1[l][0]->GetYaxis()->SetTitle("Edep(FRH1) [GeV]");
        gHistoManager->Add(hEdepFTH1vsFRH1[l][0],g[l][0]);
    }

    for (Int_t l = 1; l < 7; l++) {
        TString hname = Form("hEdepFRH1vsFRH2_lev%d",l);
        hEdepFRH1vsFRH2[l][0] =new TH2F(hname,"",1000,0.,0.5,1000,0.,0.5);
        hEdepFRH1vsFRH2[l][0]->GetXaxis()->SetTitle("Edep(FRH1) [GeV]");
        hEdepFRH1vsFRH2[l][0]->GetYaxis()->SetTitle("Edep(FRH2) [GeV]");
        gHistoManager->Add(hEdepFRH1vsFRH2[l][0],g[l][0]);
    }

    for (Int_t l = 1; l < 7; l++) {
        TString hname = Form("hEdepFRH2vsFRH3_lev%d",l);
        hEdepFRH2vsFRH3[l][0] =new TH2F(hname,"",1000,0.,0.5,1000,0.,0.5);
        hEdepFRH2vsFRH3[l][0]->GetXaxis()->SetTitle("Edep(FRH2) [GeV]");
        hEdepFRH2vsFRH3[l][0]->GetYaxis()->SetTitle("Edep(FRH3) [GeV]");
        gHistoManager->Add(hEdepFRH2vsFRH3[l][0],g[l][0]);
    }

    for (Int_t l = 1; l < 7; l++) {
        TString hname = Form("hEdepFWC1vsFRH1FRH2FRH3_lev%d",l);
        hEdepFWC1vsFRH1FRH2FRH3[l][0] =new TH2F(hname,"",1000,0.,0.01,1000,0.,1.);
        hEdepFWC1vsFRH1FRH2FRH3[l][0]->GetXaxis()->SetTitle("Edep(FWC1) [GeV]");
        hEdepFWC1vsFRH1FRH2FRH3[l][0]->GetYaxis()->SetTitle("Edep(FRH1+FRH2+FRH3) [GeV]");
        gHistoManager->Add(hEdepFWC1vsFRH1FRH2FRH3[l][0],g[l][0]);
    }

    for (Int_t l = 1; l < 7; l++) {
        TString hname = Form("hThetaFD_lev%d",l);
        hTheta_FDC[l][0] = new TH1F(hname,"",250,0.,25.);
        hTheta_FDC[l][0]->GetXaxis()->SetTitle("#theta^{lab}_{FD} [deg]");
        hTheta_FDC[l][0]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hTheta_FDC[l][0],g[l][0]);
    }

    for (Int_t l = 1; l < 7; l++) {
        TString hname = Form("hPhiFD_lev%d",l);
        hPhi_FDC[l][0] = new TH1F(hname,"",360,-180.,180.);
        hPhi_FDC[l][0]->GetXaxis()->SetTitle("#phi^{lab}_{FD} [deg]");
        hPhi_FDC[l][0]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hPhi_FDC[l][0],g[l][0]);
    }

    for (Int_t l = 1; l < 7; l++) {
        TString hname = Form("hTimeFD_lev%d",l);
        hTime_FDC[l][0] = new TH1F(hname,"",1000,0.,2500.);
        hTime_FDC[l][0]->GetXaxis()->SetTitle("Time_{FD} [ns]");
        hTime_FDC[l][0]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hTime_FDC[l][0],g[l][0]);
    }

    //charged in CD//
    for (Int_t l = 1; l < 7; l++) {
        TString hname = Form("hEdepPSBvsSEC_lev%d",l);
        hEdepPSBvsSEC[l][0] = new TH2F(hname,"",1000,0.,0.75,1000,0.,0.025);
        hEdepPSBvsSEC[l][0]->GetXaxis()->SetTitle("E_{dep(SEC)} [GeV]");
        hEdepPSBvsSEC[l][0]->GetYaxis()->SetTitle("E_{dep(PSB)} [GeV]");
        gHistoManager->Add(hEdepPSBvsSEC[l][0],g[l][0]);
    }

    for (Int_t l = 1; l < 7; l++) {
        TString hname = Form("hEdepPSBvsSigMom_lev%d",l);
        hEdepPSBvsSigMom[l][0] = new TH2F(hname,"",1000,-2.5,2.5,1000,0.,0.025);
        hEdepPSBvsSigMom[l][0]->GetXaxis()->SetTitle("Momentum [GeV/c]");
        hEdepPSBvsSigMom[l][0]->GetYaxis()->SetTitle("E_{dep(SEC)} [GeV]");
        gHistoManager->Add(hEdepPSBvsSigMom[l][0],g[l][0]);
    }

    for (Int_t l = 1; l < 7; l++) {
        TString hname = Form("hEdepSECvsSigMom_lev%d",l);
        hEdepSECvsSigMom[l][0] = new TH2F(hname,"",1000,-2.5,2.5,1000,0.,0.75);
        hEdepSECvsSigMom[l][0]->GetXaxis()->SetTitle("Momentum [GeV/c]");
        hEdepSECvsSigMom[l][0]->GetYaxis()->SetTitle("E_{dep(PSB)} [GeV]");
        gHistoManager->Add(hEdepSECvsSigMom[l][0],g[l][0]);
    }

    for (Int_t l = 1; l < 7; l++) {
        TString hname = Form("hMomentumCD_lev%d",l);
        hMom_CDC[l][0] = new TH1F(hname,"",750,0.,2.5);
        hMom_CDC[l][0]->GetXaxis()->SetTitle("p^{lab}_{p} [GeV/c]");
        hMom_CDC[l][0]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hMom_CDC[l][0],g[l][0]);
    }

    for (Int_t l = 1; l < 7; l++) {
        TString hname = Form("hThetaCD_lev%d",l);
        hTheta_CDC[l][0] = new TH1F(hname,"",360,0.,180.);
        hTheta_CDC[l][0]->GetXaxis()->SetTitle("#theta^{lab}_{CD} [deg]");
        hTheta_CDC[l][0]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hTheta_CDC[l][0],g[l][0]);
    }

    for (Int_t l = 1; l < 7; l++) {
        TString hname = Form("hPhiCD_lev%d",l);
        hPhi_CDC[l][0] = new TH1F(hname,"",360,-180.,180.);
        hPhi_CDC[l][0]->GetXaxis()->SetTitle("#phi^{lab}_{CD} [deg]");
        hPhi_CDC[l][0]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hPhi_CDC[l][0],g[l][0]);
    }

    for (Int_t l = 1; l < 7; l++) {
        TString hname = Form("hTimeCD_lev%d",l);
        hTime_CDC[l][0] = new TH1F(hname,"",1000,0.,2500.);
        hTime_CDC[l][0]->GetXaxis()->SetTitle("Time_{CD} [ns]");
        hTime_CDC[l][0]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hTime_CDC[l][0],g[l][0]);
    }

    for (Int_t l = 1; l < 7; l++) {
        TString hname = Form("hThetaFDvsThetaCD_lev%d",l);
        hTheta_FDvsTheta_CD[l][0] = new TH2F(hname,"",540,0.,20.,540,0.,180.);
        hTheta_FDvsTheta_CD[l][0]->GetXaxis()->SetTitle("#theta_{FD} [deg]");
        hTheta_FDvsTheta_CD[l][0]->GetYaxis()->SetTitle("#theta_{CD} [deg]");
        gHistoManager->Add(hTheta_FDvsTheta_CD[l][0],g[l][0]);
    }

    for (Int_t l = 1; l < 7; l++) {
        TString hname = Form("hTimeDifference_lev%d",l);
        hDeltaTime[l][0] = new TH1F(hname,"",1000,-50.,50.);
        hDeltaTime[l][0]->GetXaxis()->SetTitle("Time difference [ns]");
        hDeltaTime[l][0]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hDeltaTime[l][0],g[l][0]);
    }

    for (Int_t l = 1; l < 7; l++) {
        TString hname = Form("hDeltaPhi_lev%d",l);
        hDelta_Phi[l][0] = new TH1F(hname,"",360,0.,360.);
        hDelta_Phi[l][0]->GetXaxis()->SetTitle("#Delta#phi [deg]");
        hDelta_Phi[l][0]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hDelta_Phi[l][0],g[l][0]);
    }

    for (Int_t l = 1; l < 7; l++) {
        TString hname = Form("hDeltaPhi_abs_lev%d",l);
        hDelta_Phi_abs[l][0] = new TH1F(hname,"",360,0.,360.);
        hDelta_Phi_abs[l][0]->GetXaxis()->SetTitle("|#Delta#phi| [deg]");
        hDelta_Phi_abs[l][0]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hDelta_Phi_abs[l][0],g[l][0]);
    }

    for (Int_t l = 1; l < 7; l++) {
        TString hname = Form("hDeltaPhi_sym_lev%d",l);
        hDelta_Phi_sym[l][0] = new TH1F(hname,"",360,0.,360.);
        hDelta_Phi_sym[l][0]->GetXaxis()->SetTitle("(2#pi+#Delta#phi)mod2#pi [deg]");
        hDelta_Phi_sym[l][0]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hDelta_Phi_sym[l][0],g[l][0]);
    }

    for (Int_t l = 1; l < 7; l++) {
        TString hpath = g[l][0] + Form("/DeltaPhi_lev%d",l);
        for (Int_t i = 1; i < 41; i++) {
            TString name = Form("hDeltaPhi_lev%d_Q%d",l,i);
            TString descr = Form("Q#in(%G,%G)GeV",-70.+(i-1)*2.5,-67.5+(i-1)*2.5);
            hDelta_Phi_sym_Q[l][0][i] = new TH1D(name,descr,360,0.,360.);
            hDelta_Phi_sym_Q[l][0][i]->GetXaxis()->SetTitle("(2#pi+#Delta#phi)mod2#pi [deg]");
            hDelta_Phi_sym_Q[l][0][i]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hDelta_Phi_sym_Q[l][0][i],hpath);
        }
    }

    for (Int_t l = 1; l < 7; l++) {
        TString hpath = g[l][0] + Form("/ThetaCD_lev%d",l);
        for (Int_t j = 1; j < 16; j++) {
            TString name = Form("hThetaCD_lev%d_Th%d",l,j);
            TString descr = Form("#theta_{FD}#in(%G,%G)deg",3.+(j-1),4.+(j-1));
            hTheta_CDC_Th[l][0][j] = new TH1D(name,descr,360,0.,180.);
            hTheta_CDC_Th[l][0][j]->GetXaxis()->SetTitle("#theta_{CD} [deg]");
            hTheta_CDC_Th[l][0][j]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hTheta_CDC_Th[l][0][j],hpath);
        }
    }

    for (Int_t l = 1; l < 7; l++) {
        TString hpath = g[l][0] + Form("/ThetaCD_lev%d/ThetaCD_Q",l);
        for (Int_t i = 1; i < 41; i++) {
            TString name = Form("hThetaCD_lev%d_Q%d",l,i);
            TString descr = Form("Q#in(%G,%G)GeV",-70.+(i-1)*2.5,-67.5+(i-1)*2.5);
            hTheta_CDC_Q[l][0][i] = new TH1D(name,descr,360,0.,180.);
            hTheta_CDC_Q[l][0][i]->GetXaxis()->SetTitle("#theta_{CD} [deg]");
            hTheta_CDC_Q[l][0][i]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hTheta_CDC_Q[l][0][i],hpath);
        }
    }

    for (Int_t l = 1; l < 7; l++) {
        for (Int_t i = 1; i < 41; i++) {
            TString hpath = g[l][0] + Form("/ThetaCD_lev%d/ThetaCD_Q%d",l,i);
            for (Int_t j = 1; j < 16; j++) {
                TString name = Form("hThetaCD_lev%d_Q%d_Th%d",l,i,j);
                TString descr = Form("Q#in(%G,%G)GeV, #theta_{FD}#in(%G,%G)deg",-70.+(i-1)*2.5,-67.5+(i-1)*2.5,3.+(j-1),4.+(j-1));
                hTheta_CDC_Q_Th[l][0][i][j] = new TH1D(name,descr,360,0.,180.);
                hTheta_CDC_Q_Th[l][0][i][j]->GetXaxis()->SetTitle("#theta_{CD} [deg]");
                hTheta_CDC_Q_Th[l][0][i][j]->GetYaxis()->SetTitle("counts");
                gHistoManager->Add(hTheta_CDC_Q_Th[l][0][i][j],hpath);
            }
        }
    }

////SYSTEMATICS////

    for (Int_t c = 1; c < 9; c++) {
        TString hname = Form("hQ_syst_lev5_cut%d",c);
        hQ[5][c] = new TH1F(hname,"",40,-70.,30.);
        hQ[5][c]->GetXaxis()->SetTitle("Q [MeV]");
        hQ[5][c]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hQ[5][c],g[5][c]);
    }

    for (Int_t c = 1; c < 9; c++) {
        TString hname = Form("hQ_wght_bilin_syst_lev5_cut%d",c);
        hQ_wght_bilin[5][c] = new TH1F(hname,"",40,-70.,30.);
        hQ_wght_bilin[5][c]->GetXaxis()->SetTitle("Q [MeV]");
        hQ_wght_bilin[5][c]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hQ_wght_bilin[5][c],g[5][c]);
    }

    for (Int_t c = 1; c < 9; c++) {
        TString hname = Form("hQ_wght_prx_syst_lev5_cut%d",c);
        hQ_wght_prx[5][c] = new TH1F(hname,"",40,-70.,30.);
        hQ_wght_prx[5][c]->GetXaxis()->SetTitle("Q [MeV]");
        hQ_wght_prx[5][c]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hQ_wght_prx[5][c],g[5][c]);
    }

    for (Int_t c = 1; c < 9; c++) {
        TString hname = Form("hQ_wght_lin_syst_lev5_cut%d",c);
        hQ_wght_lin[5][c] = new TH1F(hname,"",40,-70.,30.);
        hQ_wght_lin[5][c]->GetXaxis()->SetTitle("Q [MeV]");
        hQ_wght_lin[5][c]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hQ_wght_lin[5][c],g[5][c]);
    }

    for (Int_t c = 1; c < 9; c++) {
        TString hname = Form("hQ_wght_sq_syst_lev5_cut%d",c);
        hQ_wght_sq[5][c] = new TH1F(hname,"",40,-70.,30.);
        hQ_wght_sq[5][c]->GetXaxis()->SetTitle("Q [MeV]");
        hQ_wght_sq[5][c]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hQ_wght_sq[5][c],g[5][c]);
    }

    //FD//
    for (Int_t c = 1; c < 9; c++) {
        TString hname = Form("hEdepFWC1vsFRH1_syst_lev5_cut%d",c);
        hEdepFWC1vsFRH1[5][c] =new TH2F(hname,"",1000,0.,0.01,1000,0.,0.5);
        hEdepFWC1vsFRH1[5][c]->GetXaxis()->SetTitle("Edep(FWC1) [GeV]");
        hEdepFWC1vsFRH1[5][c]->GetYaxis()->SetTitle("Edep(FRH1) [GeV]");
        gHistoManager->Add(hEdepFWC1vsFRH1[5][c],g[5][c]);
    }

    for (Int_t c = 1; c < 9; c++) {
        TString hname = Form("hEdepFWC2vsFRH1_syst_lev5_cut%d",c);
        hEdepFWC2vsFRH1[5][c] =new TH2F(hname,"",1000,0.,0.01,1000,0.,0.5);
        hEdepFWC2vsFRH1[5][c]->GetXaxis()->SetTitle("Edep(FWC2) [GeV]");
        hEdepFWC2vsFRH1[5][c]->GetYaxis()->SetTitle("Edep(FRH1) [GeV]");
        gHistoManager->Add(hEdepFWC2vsFRH1[5][c],g[5][c]);
    }

    for (Int_t c = 1; c < 9; c++) {
        TString hname = Form("hEdepFTH1vsFRH1_syst_lev5_cut%d",c);
        hEdepFTH1vsFRH1[5][c] =new TH2F(hname,"",1000,0.,0.02,1000,0.,0.5);
        hEdepFTH1vsFRH1[5][c]->GetXaxis()->SetTitle("Edep(FTH1) [GeV]");
        hEdepFTH1vsFRH1[5][c]->GetYaxis()->SetTitle("Edep(FRH1) [GeV]");
        gHistoManager->Add(hEdepFTH1vsFRH1[5][c],g[5][c]);
    }

    for (Int_t c = 1; c < 9; c++) {
        TString hname = Form("hEdepFRH1vsFRH2_syst_lev5_cut%d",c);
        hEdepFRH1vsFRH2[5][c] =new TH2F(hname,"",1000,0.,0.5,1000,0.,0.5);
        hEdepFRH1vsFRH2[5][c]->GetXaxis()->SetTitle("Edep(FRH1) [GeV]");
        hEdepFRH1vsFRH2[5][c]->GetYaxis()->SetTitle("Edep(FRH2) [GeV]");
        gHistoManager->Add(hEdepFRH1vsFRH2[5][c],g[5][c]);
    }

    for (Int_t c = 1; c < 9; c++) {
        TString hname = Form("hEdepFRH2vsFRH3_syst_lev5_cut%d",c);
        hEdepFRH2vsFRH3[5][c] =new TH2F(hname,"",1000,0.,0.5,1000,0.,0.5);
        hEdepFRH2vsFRH3[5][c]->GetXaxis()->SetTitle("Edep(FRH2) [GeV]");
        hEdepFRH2vsFRH3[5][c]->GetYaxis()->SetTitle("Edep(FRH3) [GeV]");
        gHistoManager->Add(hEdepFRH2vsFRH3[5][c],g[5][c]);
    }

    for (Int_t c = 1; c < 9; c++) {
        TString hname = Form("hEdepFWC1vsFRH1FRH2FRH3_syst_lev5_cut%d",c);
        hEdepFWC1vsFRH1FRH2FRH3[5][c] =new TH2F(hname,"",1000,0.,0.01,1000,0.,1.);
        hEdepFWC1vsFRH1FRH2FRH3[5][c]->GetXaxis()->SetTitle("Edep(FWC1) [GeV]");
        hEdepFWC1vsFRH1FRH2FRH3[5][c]->GetYaxis()->SetTitle("Edep(FRH1+FRH2+FRH3) [GeV]");
        gHistoManager->Add(hEdepFWC1vsFRH1FRH2FRH3[5][c],g[5][c]);
    }

    for (Int_t c = 1; c < 9; c++) {
        TString hname = Form("hThetaFD_syst_lev5_cut%d",c);
        hTheta_FDC[5][c] = new TH1F(hname,"",250,0.,25.);
        hTheta_FDC[5][c]->GetXaxis()->SetTitle("#theta^{lab}_{FD} [deg]");
        hTheta_FDC[5][c]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hTheta_FDC[5][c],g[5][c]);
    }

    for (Int_t c = 1; c < 9; c++) {
        TString hname = Form("hPhiFD_syst_lev5_cut%d",c);
        hPhi_FDC[5][c] = new TH1F(hname,"",360,-180.,180.);
        hPhi_FDC[5][c]->GetXaxis()->SetTitle("#phi^{lab}_{FD} [deg]");
        hPhi_FDC[5][c]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hPhi_FDC[5][c],g[5][c]);
    }

    for (Int_t c = 1; c < 9; c++) {
        TString hname = Form("hTimeFD_syst_lev5_cut%d",c);
        hTime_FDC[5][c] = new TH1F(hname,"",1000,0.,2500.);
        hTime_FDC[5][c]->GetXaxis()->SetTitle("Time_{FD} [ns]");
        hTime_FDC[5][c]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hTime_FDC[5][c],g[5][c]);
    }

    //charged in CD//
    for (Int_t c = 1; c < 9; c++) {
        TString hname = Form("hEdepPSBvsSEC_syst_lev5_cut%d",c);
        hEdepPSBvsSEC[5][c] = new TH2F(hname,"",1000,0.,0.75,1000,0.,0.025);
        hEdepPSBvsSEC[5][c]->GetXaxis()->SetTitle("E_{dep(SEC)} [GeV]");
        hEdepPSBvsSEC[5][c]->GetYaxis()->SetTitle("E_{dep(PSB)} [GeV]");
        gHistoManager->Add(hEdepPSBvsSEC[5][c],g[5][c]);
    }

    for (Int_t c = 1; c < 9; c++) {
        TString hname = Form("hEdepPSBvsSigMom_syst_lev5_cut%d",c);
        hEdepPSBvsSigMom[5][c] = new TH2F(hname,"",1000,-2.5,2.5,1000,0.,0.025);
        hEdepPSBvsSigMom[5][c]->GetXaxis()->SetTitle("Momentum [GeV/c]");
        hEdepPSBvsSigMom[5][c]->GetYaxis()->SetTitle("E_{dep(SEC)} [GeV]");
        gHistoManager->Add(hEdepPSBvsSigMom[5][c],g[5][c]);
    }

    for (Int_t c = 1; c < 9; c++) {
        TString hname = Form("hEdepSECvsSigMom_syst_lev5_cut%d",c);
        hEdepSECvsSigMom[5][c] = new TH2F(hname,"",1000,-2.5,2.5,1000,0.,0.75);
        hEdepSECvsSigMom[5][c]->GetXaxis()->SetTitle("Momentum [GeV/c]");
        hEdepSECvsSigMom[5][c]->GetYaxis()->SetTitle("E_{dep(PSB)} [GeV]");
        gHistoManager->Add(hEdepSECvsSigMom[5][c],g[5][c]);
    }

    for (Int_t c = 1; c < 9; c++) {
        TString hname = Form("hMomentumCD_syst_lev5_cut%d",c);
        hMom_CDC[5][c] = new TH1F(hname,"",750,0.,2.5);
        hMom_CDC[5][c]->GetXaxis()->SetTitle("p^{lab}_{p} [GeV/c]");
        hMom_CDC[5][c]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hMom_CDC[5][c],g[5][c]);
    }

    for (Int_t c = 1; c < 9; c++) {
        TString hname = Form("hThetaCD_syst_lev5_cut%d",c);
        hTheta_CDC[5][c] = new TH1F(hname,"",360,0.,180.);
        hTheta_CDC[5][c]->GetXaxis()->SetTitle("#theta^{lab}_{CD} [deg]");
        hTheta_CDC[5][c]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hTheta_CDC[5][c],g[5][c]);
    }

    for (Int_t c = 1; c < 9; c++) {
        TString hname = Form("hPhiCD_syst_lev5_cut%d",c);
        hPhi_CDC[5][c] = new TH1F(hname,"",360,-180.,180.);
        hPhi_CDC[5][c]->GetXaxis()->SetTitle("#phi^{lab}_{CD} [deg]");
        hPhi_CDC[5][c]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hPhi_CDC[5][c],g[5][c]);
    }

    for (Int_t c = 1; c < 9; c++) {
        TString hname = Form("hTimeCD_syst_lev5_cut%d",c);
        hTime_CDC[5][c] = new TH1F(hname,"",1000,0.,2500.);
        hTime_CDC[5][c]->GetXaxis()->SetTitle("Time_{CD} [ns]");
        hTime_CDC[5][c]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hTime_CDC[5][c],g[5][c]);
    }

    for (Int_t c = 1; c < 9; c++) {
        TString hname = Form("hThetaFDvsThetaCD_syst_lev5_cut%d",c);
        hTheta_FDvsTheta_CD[5][c] = new TH2F(hname,"",540,0.,20.,540,0.,180.);
        hTheta_FDvsTheta_CD[5][c]->GetXaxis()->SetTitle("#theta_{FD} [deg]");
        hTheta_FDvsTheta_CD[5][c]->GetYaxis()->SetTitle("#theta_{CD} [deg]");
        gHistoManager->Add(hTheta_FDvsTheta_CD[5][c],g[5][c]);
    }

    for (Int_t c = 1; c < 9; c++) {
        TString hname = Form("hTimeDifference_syst_lev5_cut%d",c);
        hDeltaTime[5][c] = new TH1F(hname,"",1000,-50.,50.);
        hDeltaTime[5][c]->GetXaxis()->SetTitle("Time difference [ns]");
        hDeltaTime[5][c]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hDeltaTime[5][c],g[5][c]);
    }

    for (Int_t c = 1; c < 9; c++) {
        TString hname = Form("hDeltaPhi_syst_lev5_cut%d",c);
        hDelta_Phi[5][c] = new TH1F(hname,"",360,0.,360.);
        hDelta_Phi[5][c]->GetXaxis()->SetTitle("#Delta#phi [deg]");
        hDelta_Phi[5][c]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hDelta_Phi[5][c],g[5][c]);
    }

    for (Int_t c = 1; c < 9; c++) {
        TString hname = Form("hDeltaPhi_abs_syst_lev5_cut%d",c);
        hDelta_Phi_abs[5][c] = new TH1F(hname,"",360,0.,360.);
        hDelta_Phi_abs[5][c]->GetXaxis()->SetTitle("|#Delta#phi| [deg]");
        hDelta_Phi_abs[5][c]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hDelta_Phi_abs[5][c],g[5][c]);
    }

    for (Int_t c = 1; c < 9; c++) {
        TString hname = Form("hDeltaPhi_sym_syst_lev5_cut%d",c);
        hDelta_Phi_sym[5][c] = new TH1F(hname,"",360,0.,360.);
        hDelta_Phi_sym[5][c]->GetXaxis()->SetTitle("(2#pi+#Delta#phi)mod2#pi [deg]");
        hDelta_Phi_sym[5][c]->GetYaxis()->SetTitle("counts");
        gHistoManager->Add(hDelta_Phi_sym[5][c],g[5][c]);
    }

    for (Int_t c = 1; c < 9; c++) {
        TString hpath = g[5][c] + Form("/DeltaPhi_syst_lev5_cut%d",c);
        for (Int_t i = 1; i < 41; i++) {
            TString name = Form("hDeltaPhi_syst_lev5_cut%d_Q%d",c,i);
            TString descr = Form("Q#in(%G,%G)GeV",-70.+(i-1)*2.5,-67.5+(i-1)*2.5);
            hDelta_Phi_sym_Q[5][c][i] = new TH1D(name,descr,360,0.,360.);
            hDelta_Phi_sym_Q[5][c][i]->GetXaxis()->SetTitle("(2#pi+#Delta#phi)mod2#pi [deg]");
            hDelta_Phi_sym_Q[5][c][i]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hDelta_Phi_sym_Q[5][c][i],hpath);
        }
    }

    for (Int_t c = 1; c < 9; c++) {
        TString hpath = g[5][c] + Form("/ThetaCD_syst_lev5_cut%d",c);
        for (Int_t j = 1; j < 16; j++) {
            TString name = Form("hThetaCD_syst_lev5_cut%d_Th%d",c,j);
            TString descr = Form("#theta_{FD}#in(%G,%G)deg",3.+(j-1),4.+(j-1));
            hTheta_CDC_Th[5][c][j] = new TH1D(name,descr,360,0.,180.);
            hTheta_CDC_Th[5][c][j]->GetXaxis()->SetTitle("#theta_{CD} [deg]");
            hTheta_CDC_Th[5][c][j]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hTheta_CDC_Th[5][c][j],hpath);
        }
    }

    for (Int_t c = 1; c < 9; c++) {
        TString hpath = g[5][c] + Form("/ThetaCD_syst_lev5_cut%d/ThetaCD_Q",c);
        for (Int_t i = 1; i < 41; i++) {
            TString name = Form("hThetaCD_syst_lev5_cut%d_Q%d",c,i);
            TString descr = Form("Q#in(%G,%G)GeV",-70.+(i-1)*2.5,-67.5+(i-1)*2.5);
            hTheta_CDC_Q[5][c][i] = new TH1D(name,descr,360,0.,180.);
            hTheta_CDC_Q[5][c][i]->GetXaxis()->SetTitle("#theta_{CD} [deg]");
            hTheta_CDC_Q[5][c][i]->GetYaxis()->SetTitle("counts");
            gHistoManager->Add(hTheta_CDC_Q[5][c][i],hpath);
        }
    }

    for (Int_t c = 1; c < 9; c++) {
        for (Int_t i = 1; i < 41; i++) {
            TString hpath = g[5][c] + Form("/ThetaCD_syst_lev5_cut%d/ThetaCD_Q%d",c,i);
            for (Int_t j = 1; j < 16; j++) {
                TString name = Form("hThetaCD_syst_lev5_cut%d_Q%d_Th%d",c,i,j);
                TString descr = Form("Q#in(%G,%G)GeV, #theta_{FD}#in(%G,%G)deg",-70.+(i-1)*2.5,-67.5+(i-1)*2.5,3.+(j-1),4.+(j-1));
                hTheta_CDC_Q_Th[5][c][i][j] = new TH1D(name,descr,360,0.,180.);
                hTheta_CDC_Q_Th[5][c][i][j]->GetXaxis()->SetTitle("#theta_{CD} [deg]");
                hTheta_CDC_Q_Th[5][c][i][j]->GetYaxis()->SetTitle("counts");
                gHistoManager->Add(hTheta_CDC_Q_Th[5][c][i][j],hpath);
            }
        }
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
