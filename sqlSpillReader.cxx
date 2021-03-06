#include <iostream>
#include <cmath>
#include <stdlib.h>

#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>

#include "DataStruct.h"

using namespace std;

int main(int argc, char* argv[])
{
    // define the output file structure
    Spill* p_spill = new Spill; Spill& spill = *p_spill;
    spill.skipflag = false;

    TFile* saveFile = new TFile(argv[2], "recreate");
    TTree* saveTree = new TTree("save", "save");

    saveTree->Branch("spill", &p_spill, 256000, 99);

    //Connect to server
    TSQLServer* server = TSQLServer::Connect(Form("mysql://%s:%d", argv[3], atoi(argv[4])), "seaguest", "qqbar2mu+mu-");
    server->Exec(Form("USE %s", argv[1]));
    cout << "Reading schema " << argv[1] << " and save to " << argv[2] << endl;

    char query[2000];
    //sprintf(query, "SELECT spillID FROM Spill WHERE runID in (SELECT run FROM summary.production WHERE ktracked=1) ORDER BY spillID");
    sprintf(query, "SELECT runID,spillID FROM Spill ORDER BY spillID");

    TSQLResult* res = server->Query(query);
    int nSpillsRow = res->GetRowCount();

    int nGoodSpill = 0;
    int nBadSpill_record = 0;
    int nBadSpill_quality = 0;
    for(int i = 0; i < nSpillsRow; ++i)
    {
        if(i % 100 == 0)
        {
            cout << "Converting spill " << i << "/" << nSpillsRow << endl;
        }

        //basic spillID info
        TSQLRow* row = res->Next();
        int runID = atoi(row->GetField(0));
        spill.spillID = atoi(row->GetField(1));
        delete row;

        //run configuration
        sprintf(query, "SELECT value FROM Run WHERE runID=%d AND name IN ('KMAG-Avg','MATRIX3Prescale') ORDER BY name", runID);
        TSQLResult* res_spill = server->Query(query);
        if(res_spill->GetRowCount() != 2)
        {
            ++nBadSpill_record;
            spill.log("lacks Run table info");

            delete res_spill;
            continue;
        }

        TSQLRow* row_spill = res_spill->Next();  spill.KMAG = atof(row_spill->GetField(0)); delete row_spill;
        row_spill = res_spill->Next(); spill.MATRIX3Prescale = atoi(row_spill->GetField(0)); delete row_spill;
        delete res_spill;

        //target position
        sprintf(query, "SELECT a.targetPos,b.value,a.dataQuality,a.liveProton FROM Spill AS a,Target AS b WHERE a.spillID"
            "=%d AND b.spillID=%d AND b.name='TARGPOS_CONTROL' AND a.liveProton IS NOT NULL", spill.spillID, spill.spillID);
        res_spill = server->Query(query);
        if(res_spill->GetRowCount() != 1)
        {
            ++nBadSpill_record;
            spill.log("lacks target position info");

            delete res_spill;
            continue;
        }

        row_spill = res_spill->Next();
        spill.targetPos       = atoi(row_spill->GetField(0));
        spill.TARGPOS_CONTROL = atoi(row_spill->GetField(1));
        spill.quality         = atoi(row_spill->GetField(2));
        spill.liveProton      = atof(row_spill->GetField(3));
        delete row_spill;
        delete res_spill;

        //Reconstruction info
        sprintf(query, "SELECT (SELECT COUNT(*) FROM Event WHERE spillID=%d),"
                              "(SELECT COUNT(*) FROM kDimuon WHERE spillID=%d),"
                              "(SELECT COUNT(*) FROM kTrack WHERE spillID=%d)", spill.spillID, spill.spillID, spill.spillID);
        res_spill = server->Query(query);
        if(res_spill->GetRowCount() != 1)
        {
            ++nBadSpill_record;
            spill.log("lacks reconstructed tables");

            delete res_spill;
            continue;
        }

        row_spill = res_spill->Next();
        spill.nEvents = atoi(row_spill->GetField(0));
        spill.nDimuons = atoi(row_spill->GetField(1));
        spill.nTracks = atoi(row_spill->GetField(2));
        delete row_spill;
        delete res_spill;

        //Beam/BeamDAQ
        sprintf(query, "SELECT a.value,b.NM3ION,b.QIESum,b.inhibit_block_sum,b.trigger_sum_no_inhibit,"
            "b.dutyFactor53MHz FROM Beam AS a,BeamDAQ AS b WHERE a.spillID=%d AND b.spillID=%d AND "
            "a.name='S:G2SEM'", spill.spillID, spill.spillID);
        res_spill = server->Query(query);
        if(res_spill->GetRowCount() != 1)
        {
            ++nBadSpill_record;
            spill.log("lacks Beam/BeamDAQ info");

            delete res_spill;
            continue;
        }

        row_spill = res_spill->Next();
        spill.G2SEM      = atof(row_spill->GetField(0));
        spill.NM3ION     = atof(row_spill->GetField(1));
        spill.QIESum     = atof(row_spill->GetField(2));
        spill.inhibitSum = atof(row_spill->GetField(3));
        spill.busySum    = atof(row_spill->GetField(4));
        spill.dutyFactor = atof(row_spill->GetField(5));

        delete row_spill;
        delete res_spill;

        //Scalar table
        sprintf(query, "SELECT value FROM Scaler WHERE spillType='EOS' AND spillID=%d AND scalerName "
            "in ('TSGo','AcceptedMatrix1','AfterInhMatrix1') ORDER BY scalerName", spill.spillID);
        res_spill = server->Query(query);
        if(res_spill->GetRowCount() != 3)
        {
            ++nBadSpill_record;
            spill.log("lacks scaler info");

            delete res_spill;
            continue;
        }

        row_spill = res_spill->Next(); spill.acceptedMatrix1 = atof(row_spill->GetField(0)); delete row_spill;
        row_spill = res_spill->Next(); spill.afterInhMatrix1 = atof(row_spill->GetField(0)); delete row_spill;
        row_spill = res_spill->Next(); spill.TSGo            = atof(row_spill->GetField(0)); delete row_spill;
        delete res_spill;

        if(spill.goodSpill())
        {
            saveTree->Fill();
            ++nGoodSpill;

            if(nGoodSpill % 100 == 0) saveTree->AutoSave("self");
        }
        else
        {
            ++nBadSpill_quality;
            spill.log("spill bad");
            spill.print();
        }
    }
    delete res;

    cout << "sqlReader finished successfully." << endl;
    cout << nGoodSpill << " good spills, " << nBadSpill_record << " spills have insufficient info in database, " << nBadSpill_quality << " rejected because of bad quality." << endl;

    saveFile->cd();
    saveTree->Write();
    saveFile->Close();

    return EXIT_SUCCESS;
}
