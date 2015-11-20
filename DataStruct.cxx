#include <cmath>
#include "DataStruct.h"

ClassImp(Event)
ClassImp(Dimuon)
ClassImp(Spill)
ClassImp(Track)

#define SPILLID_MIN_57 303215
#define SPILLID_MAX_57 370100
#define SPILLID_MIN_59 370110
#define SPILLID_MAX_59 388580
#define SPILLID_MIN_61 388611
#define SPILLID_MAX_61 391100
#define SPILLID_MIN_62 394308
#define SPILLID_MAX_62 482573
#define SPILLID_MIN_67 484946
#define SPILLID_MAX_67 676224
#define SPILLID_MIN_70 676498
#define SPILLID_MAX_70 696455

#define PEDESTAL 36.791

Event::Event() : runID(-1), spillID(-1), eventID(-1), status(-1), MATRIX1(-1), weight(0.)
{
    for(int i = 0; i < 33; ++i) intensity[i] = 0.;
}

bool Event::goodEvent()
{
    return MATRIX1 > 0 && status == 0;
}

float Event::weightedIntensity(float unit)
{
    double weight[] = {0.000814246430361413, 0.0028662467149288, 0.00597015326639906, 0.0121262946061061,
                       0.0300863195179747, 0.0777262437180552, 0.159446650644417, 0.259932709364831,
                       0.36718876894966, 0.488159093692654, 0.678969311099113, 0.847788074599439, 0.956475273764143,
                       1.,
                       0.989173954042814, 0.897678016090413, 0.767828869998712, 0.647167321489559, 0.533894756174369,
                       0.448848741080746, 0.356435437171761, 0.263693103645649, 0.177964720504253,
                       0.108504562083177, 0.0540099990325891, 0.019218568399343, 0.00308302089003216};

    double sum = 0.;
    double wsum = 0.;
    for(int i = -13; i <= 13; ++i)
    {
        sum += (weight[i+13]*(intensity[i+16] - PEDESTAL));
        wsum += weight[i+13];
    }

    return unit*sum/wsum;
}

int Event::branch_mapping(TTree *newtree, const char *prefix)
{
    newtree->Branch(Form("%s%s", prefix, "runID"), &runID, Form("%s%s", prefix, "runID/I"));
    newtree->Branch(Form("%s%s", prefix, "spillID"), &spillID, Form("%s%s", prefix, "spillID/I"));
    newtree->Branch(Form("%s%s", prefix, "eventID"), &eventID, Form("%s%s", prefix, "eventID/I"));
    newtree->Branch(Form("%s%s", prefix, "status"), &status, Form("%s%s", prefix, "status/I"));
    newtree->Branch(Form("%s%s", prefix, "MATRIX1"), &MATRIX1, Form("%s%s", prefix, "MATRIX1/I"));
    newtree->Branch(Form("%s%s", prefix, "weight"), &weight, Form("%s%s", prefix, "weight/F"));
    newtree->Branch(Form("%s%s", prefix, "intensity"), intensity, Form("%s%s", prefix, "intensity[33]/F"));
    return 0;
}

int Event::SetBranchAddress(TTree *newtree, const char *prefix)
{
    newtree->SetBranchAddress(Form("%s%s", prefix, "runID"), &runID);
    newtree->SetBranchAddress(Form("%s%s", prefix, "spillID"), &spillID);
    newtree->SetBranchAddress(Form("%s%s", prefix, "eventID"), &eventID);
    newtree->SetBranchAddress(Form("%s%s", prefix, "status"), &status);
    newtree->SetBranchAddress(Form("%s%s", prefix, "MATRIX1"), &MATRIX1);
    newtree->SetBranchAddress(Form("%s%s", prefix, "weight"), &weight);
    newtree->SetBranchAddress(Form("%s%s", prefix, "intensity"), intensity);
    return 0;
}

bool Dimuon::goodDimuon(int polarity)
{
    if(fabs(dx) > 2. || fabs(dy) > 2.) return false;
    if(dz < -300. || dz > 200.) return false;
    if(fabs(dpx) > 3. || fabs(dpy) > 3.) return false;
    if(dpz < 30. || dpz > 120.) return false;
    if(x1 < 0. || x1 > 1.) return false;
    if(x2 < 0. || x2 > 1.) return false;
    if(xF < -1. || xF > 1.) return false;
    if(fabs(trackSeparation) > 250.) return false;
    if(chisq_dimuon > 15.) return false;
    if(polarity*px1 < 0. || polarity*px2 > 0.) return false;

    return true;
}

bool Dimuon::targetDimuon()
{
    if(dz > -60. || dz < -300.) return false;
    return true;
}

bool Dimuon::dumpDimuon()
{
    if(dz < 0. || dz > 150.) return false;
    return true;
}

int Dimuon::branch_mapping(TTree *newtree, const char *prefix)
{
    newtree->Branch(Form("%s%s", prefix, "dimuonID"), &dimuonID, Form("%s%s", prefix, "dimuonID/I"));
    newtree->Branch(Form("%s%s", prefix, "posTrackID"), &posTrackID, Form("%s%s", prefix, "posTrackID/I"));
    newtree->Branch(Form("%s%s", prefix, "negTrackID"), &negTrackID, Form("%s%s", prefix, "negTrackID/I"));
    newtree->Branch(Form("%s%s", prefix, "chisq_dimuon"), &chisq_dimuon, Form("%s%s", prefix, "chisq_dimuon/F"));
    newtree->Branch(Form("%s%s", prefix, "trackSeparation"), &trackSeparation, Form("%s%s", prefix, "trackSeparation/F"));
    newtree->Branch(Form("%s%s", prefix, "dx"), &dx, Form("%s%s", prefix, "dx/F"));
    newtree->Branch(Form("%s%s", prefix, "dy"), &dy, Form("%s%s", prefix, "dy/F"));
    newtree->Branch(Form("%s%s", prefix, "dz"), &dz, Form("%s%s", prefix, "dz/F"));
    newtree->Branch(Form("%s%s", prefix, "dpx"), &dpx, Form("%s%s", prefix, "dpx/F"));
    newtree->Branch(Form("%s%s", prefix, "dpy"), &dpy, Form("%s%s", prefix, "dpy/F"));
    newtree->Branch(Form("%s%s", prefix, "dpz"), &dpz, Form("%s%s", prefix, "dpz/F"));
    newtree->Branch(Form("%s%s", prefix, "px1"), &px1, Form("%s%s", prefix, "px1/F"));
    newtree->Branch(Form("%s%s", prefix, "py1"), &py1, Form("%s%s", prefix, "py1/F"));
    newtree->Branch(Form("%s%s", prefix, "pz1"), &pz1, Form("%s%s", prefix, "pz1/F"));
    newtree->Branch(Form("%s%s", prefix, "px2"), &px2, Form("%s%s", prefix, "px2/F"));
    newtree->Branch(Form("%s%s", prefix, "py2"), &py2, Form("%s%s", prefix, "py2/F"));
    newtree->Branch(Form("%s%s", prefix, "pz2"), &pz2, Form("%s%s", prefix, "pz2/F"));
    newtree->Branch(Form("%s%s", prefix, "mass"), &mass, Form("%s%s", prefix, "mass/F"));
    newtree->Branch(Form("%s%s", prefix, "xF"), &xF, Form("%s%s", prefix, "xF/F"));
    newtree->Branch(Form("%s%s", prefix, "x1"), &x1, Form("%s%s", prefix, "x1/F"));
    newtree->Branch(Form("%s%s", prefix, "x2"), &x2, Form("%s%s", prefix, "x2/F"));
    newtree->Branch(Form("%s%s", prefix, "pT"), &pT, Form("%s%s", prefix, "pT/F"));
    newtree->Branch(Form("%s%s", prefix, "costh"), &costh, Form("%s%s", prefix, "costh/F"));
    newtree->Branch(Form("%s%s", prefix, "phi"), &phi, Form("%s%s", prefix, "phi/F"));
    return 0;
}

int Dimuon::SetBranchAddress(TTree *newtree, const char *prefix)
{
    newtree->SetBranchAddress(Form("%s%s", prefix, "dimuonID"), &dimuonID);
    newtree->SetBranchAddress(Form("%s%s", prefix, "posTrackID"), &posTrackID);
    newtree->SetBranchAddress(Form("%s%s", prefix, "negTrackID"), &negTrackID);
    newtree->SetBranchAddress(Form("%s%s", prefix, "chisq_dimuon"), &chisq_dimuon);
    newtree->SetBranchAddress(Form("%s%s", prefix, "trackSeparation"), &trackSeparation);
    newtree->SetBranchAddress(Form("%s%s", prefix, "dx"), &dx);
    newtree->SetBranchAddress(Form("%s%s", prefix, "dy"), &dy);
    newtree->SetBranchAddress(Form("%s%s", prefix, "dz"), &dz);
    newtree->SetBranchAddress(Form("%s%s", prefix, "dpx"), &dpx);
    newtree->SetBranchAddress(Form("%s%s", prefix, "dpy"), &dpy);
    newtree->SetBranchAddress(Form("%s%s", prefix, "dpz"), &dpz);
    newtree->SetBranchAddress(Form("%s%s", prefix, "px1"), &px1);
    newtree->SetBranchAddress(Form("%s%s", prefix, "py1"), &py1);
    newtree->SetBranchAddress(Form("%s%s", prefix, "pz1"), &pz1);
    newtree->SetBranchAddress(Form("%s%s", prefix, "px2"), &px2);
    newtree->SetBranchAddress(Form("%s%s", prefix, "py2"), &py2);
    newtree->SetBranchAddress(Form("%s%s", prefix, "pz2"), &pz2);
    newtree->SetBranchAddress(Form("%s%s", prefix, "mass"), &mass);
    newtree->SetBranchAddress(Form("%s%s", prefix, "xF"), &xF);
    newtree->SetBranchAddress(Form("%s%s", prefix, "x1"), &x1);
    newtree->SetBranchAddress(Form("%s%s", prefix, "x2"), &x2);
    newtree->SetBranchAddress(Form("%s%s", prefix, "pT"), &pT);
    newtree->SetBranchAddress(Form("%s%s", prefix, "costh"), &costh);
    newtree->SetBranchAddress(Form("%s%s", prefix, "phi"), &phi);
    return 0;
}

Spill::Spill() : spillID(-1), quality(-1), targetPos(-1), TARGPOS_CONTROL(-1), nEvents(0), nTracks(0), nDimuons(0), skipflag(false)
{}

bool Spill::goodSpill()
{
    return skipflag || quality == 0; //(goodTargetPos() && goodTSGo() && goodScaler() && goodBeam() && goodBeamDAQ());
}

bool Spill::goodTargetPos()
{
    if(targetPos != TARGPOS_CONTROL) return false;
    if(targetPos < 1 || targetPos > 7) return false;

    return true;
}

bool Spill::goodTSGo()
{
    int trigSet = triggerSet();
    if(trigSet < 0)
    {
        return false;
    }
    else if(trigSet <= 61)
    {
        if(TSGo < 1E3 || TSGo > 8E3) return false;
    }
    else if(trigSet <= 70)
    {
        if(TSGo < 1E2 || TSGo > 6E3) return false;
    }
    else
    {
        return false;
    }

    return true;
}

bool Spill::goodScaler()
{
    int trigSet = triggerSet();
    if(trigSet < 0)
    {
        return false;
    }
    else if(trigSet <= 61)
    {
        if(acceptedMatrix1 < 1E3 || acceptedMatrix1 > 8E3) return false;
        if(afterInhMatrix1 < 1E3 || afterInhMatrix1 > 3E4) return false;
        if(acceptedMatrix1/afterInhMatrix1 < 0.2 || acceptedMatrix1/afterInhMatrix1 > 0.9) return false;
    }
    else if(trigSet <= 70)
    {
        if(acceptedMatrix1 < 1E2 || acceptedMatrix1 > 6E3) return false;
        if(afterInhMatrix1 < 1E2 || afterInhMatrix1 > 1E4) return false;
        if(acceptedMatrix1/afterInhMatrix1 < 0.2 || acceptedMatrix1/afterInhMatrix1 > 1.05) return false;
    }
    else
    {
        return false;
    }

    return true;
}

bool Spill::goodBeam()
{
    int trigSet = triggerSet();
    if(trigSet < 0)
    {
        return false;
    }
    else if(trigSet <= 61)
    {
        //if(NM3ION < 2E12 || NM3ION > 1E13) return false;
        if(G2SEM < 2E12 || G2SEM > 1E13) return false;
        if(dutyFactor < 15. || dutyFactor > 60.) return false;
    }
    else if(trigSet <= 70)
    {
        //if(NM3ION < 2E12 || NM3ION > 1E13) return false;
        if(G2SEM < 2E12 || G2SEM > 1E13) return false;
        if(dutyFactor < 10. || dutyFactor > 60.) return false;
    }
    else
    {
        return false;
    }

    return true;
}

bool Spill::goodBeamDAQ()
{
    int trigSet = triggerSet();
    if(trigSet < 0)
    {
        return false;
    }
    else if(trigSet <= 61)
    {
        if(QIESum < 4E10 || QIESum > 1E12) return false;
        if(inhibitSum < 4E9 || inhibitSum > 1E11) return false;
        if(busySum < 4E9 || busySum > 1E11) return false;
    }
    else if(trigSet <= 70)
    {
        if(QIESum < 4E10 || QIESum > 1E12) return false;
        if(inhibitSum < 4E9 || inhibitSum > 2E11) return false;
        if(busySum < 4E9 || busySum > 1E11) return false;
    }
    else
    {
        return false;
    }

    return true;
}

int Spill::triggerSet()
{
    if(spillID >= 303215 && spillID <= 310954) return -1;           //bad QIE range
    if(spillID >= 371870 && spillID <= 376533) return -1;           //timing shift in #59
    if(spillID >= 378366 && spillID <= 379333) return -1;           //timing shift in #59
    if(spillID >= 416207 && spillID <= 424180) return -1;           //manual target movement in #62
    if(spillID >= SPILLID_MIN_62 && spillID <= 409540) return -1;   //unstable trigger timing in #62
    if(spillID >= SPILLID_MIN_57 && spillID <= SPILLID_MAX_57) return 57;
    if(spillID >= SPILLID_MIN_59 && spillID <= SPILLID_MAX_59) return 59;
    if(spillID >= SPILLID_MIN_61 && spillID <= SPILLID_MAX_61) return 61;
    if(spillID >= SPILLID_MIN_62 && spillID <= SPILLID_MAX_62) return 62;
    if(spillID >= SPILLID_MIN_67 && spillID <= SPILLID_MAX_67) return 67;
    if(spillID >= SPILLID_MIN_70 && spillID <= SPILLID_MAX_70) return 70;
    return -1;
}

float Spill::QIEUnit()
{
    return G2SEM/(QIESum - 588*360000*PEDESTAL);
}

float Spill::liveG2SEM()
{
    return (QIESum - inhibitSum - busySum)*QIEUnit();
}

void Spill::print()
{
    using namespace std;
    cout << " trigge set:        " << triggerSet() << endl;
    cout << " targetPos:         " << targetPos << "  " << TARGPOS_CONTROL << endl;
    cout << " TSGo:              " << TSGo << endl;
    cout << " acceptedMatrix1:   " << acceptedMatrix1 << endl;
    cout << " afterInhMatrix1:   " << afterInhMatrix1 << endl;
    cout << " NM3ION:            " << NM3ION << endl;
    cout << " G2SEM:             " << G2SEM << endl;
    cout << " QIESum:            " << QIESum << endl;
    cout << " inhibitSum:        " << inhibitSum << endl;
    cout << " busySum:           " << busySum << endl;
    cout << " dutyFactor:        " << dutyFactor << endl;
    cout << " liveG2SEM:         " << liveG2SEM() << endl;
    cout << " liveProton:        " << liveProton << endl;
    cout << " QIEUnit:           " << QIEUnit() << endl;
}

int Spill::branch_mapping(TTree *newtree, const char *prefix)
{
    newtree->Branch(Form("%s%s", prefix, "TSGo"), &TSGo, Form("%s%s", prefix, "TSGo/F"));
    newtree->Branch(Form("%s%s", prefix, "acceptedMatrix1"), &acceptedMatrix1, Form("%s%s", prefix, "acceptedMatrix1/F"));
    newtree->Branch(Form("%s%s", prefix, "afterInhMatrix1"), &afterInhMatrix1, Form("%s%s", prefix, "afterInhMatrix1/F"));
    newtree->Branch(Form("%s%s", prefix, "NM3ION"), &NM3ION, Form("%s%s", prefix, "NM3ION/F"));
    newtree->Branch(Form("%s%s", prefix, "G2SEM"), &G2SEM, Form("%s%s", prefix, "G2SEM/F"));
    newtree->Branch(Form("%s%s", prefix, "QIESum"), &QIESum, Form("%s%s", prefix, "QIESum/F"));
    newtree->Branch(Form("%s%s", prefix, "inhibitSum"), &inhibitSum, Form("%s%s", prefix, "inhibitSum/F"));
    newtree->Branch(Form("%s%s", prefix, "busySum"), &busySum, Form("%s%s", prefix, "busySum/F"));
    newtree->Branch(Form("%s%s", prefix, "dutyFactor"), &dutyFactor, Form("%s%s", prefix, "dutyFactor/F"));
    newtree->Branch(Form("%s%s", prefix, "liveProton"), &liveProton, Form("%s%s", prefix, "liveProton/F"));
    newtree->Branch(Form("%s%s", prefix, "spillID"), &spillID, Form("%s%s", prefix, "spillID/I"));
    newtree->Branch(Form("%s%s", prefix, "quality"), &quality, Form("%s%s", prefix, "quality/I"));
    newtree->Branch(Form("%s%s", prefix, "targetPos"), &targetPos, Form("%s%s", prefix, "targetPos/I"));
    newtree->Branch(Form("%s%s", prefix, "TARGPOS_CONTROL"), &TARGPOS_CONTROL, Form("%s%s", prefix, "TARGPOS_CONTROL/I"));
    newtree->Branch(Form("%s%s", prefix, "nEvents"), &nEvents, Form("%s%s", prefix, "nEvents/I"));
    newtree->Branch(Form("%s%s", prefix, "nTracks"), &nTracks, Form("%s%s", prefix, "nTracks/I"));
    newtree->Branch(Form("%s%s", prefix, "nDimuons"), &nDimuons, Form("%s%s", prefix, "nDimuons/I"));
    newtree->Branch(Form("%s%s", prefix, "KMAG"), &KMAG, Form("%s%s", prefix, "KMAG/F"));
    newtree->Branch(Form("%s%s", prefix, "MATRIX3Prescale"), &MATRIX3Prescale, Form("%s%s", prefix, "MATRIX3Prescale/I"));
    newtree->Branch(Form("%s%s", prefix, "skipflag"), &skipflag, Form("%s%s", prefix, "skipflag/O"));
    return 0;
}

int Spill::SetBranchAddress(TTree *newtree, const char *prefix)
{
    newtree->SetBranchAddress(Form("%s%s", prefix, "TSGo"), &TSGo);
    newtree->SetBranchAddress(Form("%s%s", prefix, "acceptedMatrix1"), &acceptedMatrix1);
    newtree->SetBranchAddress(Form("%s%s", prefix, "afterInhMatrix1"), &afterInhMatrix1);
    newtree->SetBranchAddress(Form("%s%s", prefix, "NM3ION"), &NM3ION);
    newtree->SetBranchAddress(Form("%s%s", prefix, "G2SEM"), &G2SEM);
    newtree->SetBranchAddress(Form("%s%s", prefix, "QIESum"), &QIESum);
    newtree->SetBranchAddress(Form("%s%s", prefix, "inhibitSum"), &inhibitSum);
    newtree->SetBranchAddress(Form("%s%s", prefix, "busySum"), &busySum);
    newtree->SetBranchAddress(Form("%s%s", prefix, "dutyFactor"), &dutyFactor);
    newtree->SetBranchAddress(Form("%s%s", prefix, "liveProton"), &liveProton);
    newtree->SetBranchAddress(Form("%s%s", prefix, "spillID"), &spillID);
    newtree->SetBranchAddress(Form("%s%s", prefix, "quality"), &quality);
    newtree->SetBranchAddress(Form("%s%s", prefix, "targetPos"), &targetPos);
    newtree->SetBranchAddress(Form("%s%s", prefix, "TARGPOS_CONTROL"), &TARGPOS_CONTROL);
    newtree->SetBranchAddress(Form("%s%s", prefix, "nEvents"), &nEvents);
    newtree->SetBranchAddress(Form("%s%s", prefix, "nTracks"), &nTracks);
    newtree->SetBranchAddress(Form("%s%s", prefix, "nDimuons"), &nDimuons);
    newtree->SetBranchAddress(Form("%s%s", prefix, "KMAG"), &KMAG);
    newtree->SetBranchAddress(Form("%s%s", prefix, "MATRIX3Prescale"), &MATRIX3Prescale);
    newtree->SetBranchAddress(Form("%s%s", prefix, "skipflag"), &skipflag);
    return 0;
}

bool Track::goodTrack()
{
    if(nHits <= 14) return false;
    if(chisq/(nHits - 5) > 5.) return false;
    if(z0 < -400. || z0 > 200.) return false;
    if(roadID == 0) return false;
    if(nHits < 18 && pz1 < 18.) return false;

    return true;
}

bool Track::targetTrack()
{
    if(z0 <= -300. || z0 >= 0.) return false;
    if(chisq_dump - chisq_target < 10.) return false;

    return true;
}

bool Track::dumpTrack()
{
    if(z0 <= 0. && z0 >= 150.) return false;
    if(chisq_target - chisq_dump < 10.) return false;

    return true;
}

int Track::branch_mapping(TTree *newtree, const char *prefix)
{
    newtree->Branch(Form("%s%s", prefix, "trackID"), &trackID, Form("%s%s", prefix, "trackID/I"));
    newtree->Branch(Form("%s%s", prefix, "roadID"), &roadID, Form("%s%s", prefix, "roadID/I"));
    newtree->Branch(Form("%s%s", prefix, "charge"), &charge, Form("%s%s", prefix, "charge/I"));
    newtree->Branch(Form("%s%s", prefix, "nHits"), &nHits, Form("%s%s", prefix, "nHits/I"));
    newtree->Branch(Form("%s%s", prefix, "nHitsSt1"), &nHitsSt1, Form("%s%s", prefix, "nHitsSt1/I"));
    newtree->Branch(Form("%s%s", prefix, "nHitsSt2"), &nHitsSt2, Form("%s%s", prefix, "nHitsSt2/I"));
    newtree->Branch(Form("%s%s", prefix, "nHitsSt3"), &nHitsSt3, Form("%s%s", prefix, "nHitsSt3/I"));
    newtree->Branch(Form("%s%s", prefix, "nHitsSt4H"), &nHitsSt4H, Form("%s%s", prefix, "nHitsSt4H/I"));
    newtree->Branch(Form("%s%s", prefix, "nHitsSt4V"), &nHitsSt4V, Form("%s%s", prefix, "nHitsSt4V/I"));
    newtree->Branch(Form("%s%s", prefix, "chisq"), &chisq, Form("%s%s", prefix, "chisq/F"));
    newtree->Branch(Form("%s%s", prefix, "chisq_dump"), &chisq_dump, Form("%s%s", prefix, "chisq_dump/F"));
    newtree->Branch(Form("%s%s", prefix, "chisq_target"), &chisq_target, Form("%s%s", prefix, "chisq_target/F"));
    newtree->Branch(Form("%s%s", prefix, "chisq_upstream"), &chisq_upstream, Form("%s%s", prefix, "chisq_upstream/F"));
    newtree->Branch(Form("%s%s", prefix, "x1"), &x1, Form("%s%s", prefix, "x1/F"));
    newtree->Branch(Form("%s%s", prefix, "y1"), &y1, Form("%s%s", prefix, "y1/F"));
    newtree->Branch(Form("%s%s", prefix, "z1"), &z1, Form("%s%s", prefix, "z1/F"));
    newtree->Branch(Form("%s%s", prefix, "x3"), &x3, Form("%s%s", prefix, "x3/F"));
    newtree->Branch(Form("%s%s", prefix, "y3"), &y3, Form("%s%s", prefix, "y3/F"));
    newtree->Branch(Form("%s%s", prefix, "z3"), &z3, Form("%s%s", prefix, "z3/F"));
    newtree->Branch(Form("%s%s", prefix, "x0"), &x0, Form("%s%s", prefix, "x0/F"));
    newtree->Branch(Form("%s%s", prefix, "y0"), &y0, Form("%s%s", prefix, "y0/F"));
    newtree->Branch(Form("%s%s", prefix, "z0"), &z0, Form("%s%s", prefix, "z0/F"));
    newtree->Branch(Form("%s%s", prefix, "xT"), &xT, Form("%s%s", prefix, "xT/F"));
    newtree->Branch(Form("%s%s", prefix, "yT"), &yT, Form("%s%s", prefix, "yT/F"));
    newtree->Branch(Form("%s%s", prefix, "zT"), &zT, Form("%s%s", prefix, "zT/F"));
    newtree->Branch(Form("%s%s", prefix, "xD"), &xD, Form("%s%s", prefix, "xD/F"));
    newtree->Branch(Form("%s%s", prefix, "yD"), &yD, Form("%s%s", prefix, "yD/F"));
    newtree->Branch(Form("%s%s", prefix, "zD"), &zD, Form("%s%s", prefix, "zD/F"));
    newtree->Branch(Form("%s%s", prefix, "px0"), &px0, Form("%s%s", prefix, "px0/F"));
    newtree->Branch(Form("%s%s", prefix, "py0"), &py0, Form("%s%s", prefix, "py0/F"));
    newtree->Branch(Form("%s%s", prefix, "pz0"), &pz0, Form("%s%s", prefix, "pz0/F"));
    newtree->Branch(Form("%s%s", prefix, "px1"), &px1, Form("%s%s", prefix, "px1/F"));
    newtree->Branch(Form("%s%s", prefix, "py1"), &py1, Form("%s%s", prefix, "py1/F"));
    newtree->Branch(Form("%s%s", prefix, "pz1"), &pz1, Form("%s%s", prefix, "pz1/F"));
    newtree->Branch(Form("%s%s", prefix, "px3"), &px3, Form("%s%s", prefix, "px3/F"));
    newtree->Branch(Form("%s%s", prefix, "py3"), &py3, Form("%s%s", prefix, "py3/F"));
    newtree->Branch(Form("%s%s", prefix, "pz3"), &pz3, Form("%s%s", prefix, "pz3/F"));
    newtree->Branch(Form("%s%s", prefix, "pxT"), &pxT, Form("%s%s", prefix, "pxT/F"));
    newtree->Branch(Form("%s%s", prefix, "pyT"), &pyT, Form("%s%s", prefix, "pyT/F"));
    newtree->Branch(Form("%s%s", prefix, "pzT"), &pzT, Form("%s%s", prefix, "pzT/F"));
    newtree->Branch(Form("%s%s", prefix, "pxD"), &pxD, Form("%s%s", prefix, "pxD/F"));
    newtree->Branch(Form("%s%s", prefix, "pyD"), &pyD, Form("%s%s", prefix, "pyD/F"));
    newtree->Branch(Form("%s%s", prefix, "pzD"), &pzD, Form("%s%s", prefix, "pzD/F"));
    newtree->Branch(Form("%s%s", prefix, "pxv"), &pxv, Form("%s%s", prefix, "pxv/F"));
    newtree->Branch(Form("%s%s", prefix, "pyv"), &pyv, Form("%s%s", prefix, "pyv/F"));
    newtree->Branch(Form("%s%s", prefix, "pzv"), &pzv, Form("%s%s", prefix, "pzv/F"));
    newtree->Branch(Form("%s%s", prefix, "tx_PT"), &tx_PT, Form("%s%s", prefix, "tx_PT/F"));
    newtree->Branch(Form("%s%s", prefix, "ty_PT"), &ty_PT, Form("%s%s", prefix, "ty_PT/F"));
    newtree->Branch(Form("%s%s", prefix, "thbend"), &thbend, Form("%s%s", prefix, "thbend/F"));
    return 0;
}

int Track::SetBranchAddress(TTree *newtree, const char *prefix)
{
    newtree->SetBranchAddress(Form("%s%s", prefix, "trackID"), &trackID);
    newtree->SetBranchAddress(Form("%s%s", prefix, "roadID"), &roadID);
    newtree->SetBranchAddress(Form("%s%s", prefix, "charge"), &charge);
    newtree->SetBranchAddress(Form("%s%s", prefix, "nHits"), &nHits);
    newtree->SetBranchAddress(Form("%s%s", prefix, "nHitsSt1"), &nHitsSt1);
    newtree->SetBranchAddress(Form("%s%s", prefix, "nHitsSt2"), &nHitsSt2);
    newtree->SetBranchAddress(Form("%s%s", prefix, "nHitsSt3"), &nHitsSt3);
    newtree->SetBranchAddress(Form("%s%s", prefix, "nHitsSt4H"), &nHitsSt4H);
    newtree->SetBranchAddress(Form("%s%s", prefix, "nHitsSt4V"), &nHitsSt4V);
    newtree->SetBranchAddress(Form("%s%s", prefix, "chisq"), &chisq);
    newtree->SetBranchAddress(Form("%s%s", prefix, "chisq_dump"), &chisq_dump);
    newtree->SetBranchAddress(Form("%s%s", prefix, "chisq_target"), &chisq_target);
    newtree->SetBranchAddress(Form("%s%s", prefix, "chisq_upstream"), &chisq_upstream);
    newtree->SetBranchAddress(Form("%s%s", prefix, "x1"), &x1);
    newtree->SetBranchAddress(Form("%s%s", prefix, "y1"), &y1);
    newtree->SetBranchAddress(Form("%s%s", prefix, "z1"), &z1);
    newtree->SetBranchAddress(Form("%s%s", prefix, "x3"), &x3);
    newtree->SetBranchAddress(Form("%s%s", prefix, "y3"), &y3);
    newtree->SetBranchAddress(Form("%s%s", prefix, "z3"), &z3);
    newtree->SetBranchAddress(Form("%s%s", prefix, "x0"), &x0);
    newtree->SetBranchAddress(Form("%s%s", prefix, "y0"), &y0);
    newtree->SetBranchAddress(Form("%s%s", prefix, "z0"), &z0);
    newtree->SetBranchAddress(Form("%s%s", prefix, "xT"), &xT);
    newtree->SetBranchAddress(Form("%s%s", prefix, "yT"), &yT);
    newtree->SetBranchAddress(Form("%s%s", prefix, "zT"), &zT);
    newtree->SetBranchAddress(Form("%s%s", prefix, "xD"), &xD);
    newtree->SetBranchAddress(Form("%s%s", prefix, "yD"), &yD);
    newtree->SetBranchAddress(Form("%s%s", prefix, "zD"), &zD);
    newtree->SetBranchAddress(Form("%s%s", prefix, "px0"), &px0);
    newtree->SetBranchAddress(Form("%s%s", prefix, "py0"), &py0);
    newtree->SetBranchAddress(Form("%s%s", prefix, "pz0"), &pz0);
    newtree->SetBranchAddress(Form("%s%s", prefix, "px1"), &px1);
    newtree->SetBranchAddress(Form("%s%s", prefix, "py1"), &py1);
    newtree->SetBranchAddress(Form("%s%s", prefix, "pz1"), &pz1);
    newtree->SetBranchAddress(Form("%s%s", prefix, "px3"), &px3);
    newtree->SetBranchAddress(Form("%s%s", prefix, "py3"), &py3);
    newtree->SetBranchAddress(Form("%s%s", prefix, "pz3"), &pz3);
    newtree->SetBranchAddress(Form("%s%s", prefix, "pxT"), &pxT);
    newtree->SetBranchAddress(Form("%s%s", prefix, "pyT"), &pyT);
    newtree->SetBranchAddress(Form("%s%s", prefix, "pzT"), &pzT);
    newtree->SetBranchAddress(Form("%s%s", prefix, "pxD"), &pxD);
    newtree->SetBranchAddress(Form("%s%s", prefix, "pyD"), &pyD);
    newtree->SetBranchAddress(Form("%s%s", prefix, "pzD"), &pzD);
    newtree->SetBranchAddress(Form("%s%s", prefix, "pxv"), &pxv);
    newtree->SetBranchAddress(Form("%s%s", prefix, "pyv"), &pyv);
    newtree->SetBranchAddress(Form("%s%s", prefix, "pzv"), &pzv);
    newtree->SetBranchAddress(Form("%s%s", prefix, "tx_PT"), &tx_PT);
    newtree->SetBranchAddress(Form("%s%s", prefix, "ty_PT"), &ty_PT);
    newtree->SetBranchAddress(Form("%s%s", prefix, "thbend"), &thbend);
    return 0;
}