void flat_tree_convert()
{
    gSystem->Load("libAnaUtil.so");

    // define the output file structure
    Dimuon* p_dimuon = new Dimuon; Dimuon& dimuon = *p_dimuon;
    Spill* p_spill = new Spill; Spill& spill = *p_spill;
    Event* p_event = new Event; Event& event = *p_event;
    Track* p_posTrack = new Track; Track& posTrack = *p_posTrack;
    Track* p_negTrack = new Track; Track& negTrack = *p_negTrack;

    TFile* dataFile = new TFile("/sp8data11/shyao/rootfile/data_67.root", "READ");
    TTree* dataTree = (TTree*)dataFile->Get("save");

    dataTree->SetBranchAddress("dimuon", &p_dimuon);
    dataTree->SetBranchAddress("event", &p_event);
    dataTree->SetBranchAddress("spill", &p_spill);
    dataTree->SetBranchAddress("posTrack", &p_posTrack);
    dataTree->SetBranchAddress("negTrack", &p_negTrack);

    TFile *newfile = new TFile("/sp8data11/shyao/rootfile/clone2.root","recreate");
    TTree *newtree = new TTree("save2","pick up branch");
 
    dimuon.branch_mapping(newtree,"d_");
    event.branch_mapping(newtree,"e_");
    spill.branch_mapping(newtree,"s_");
    posTrack.branch_mapping(newtree,"p_");
    negTrack.branch_mapping(newtree,"n_");
	
    //loop over all the events
    for(int i = 0; i < dataTree->GetEntries(); ++i)
    {
        dataTree->GetEntry(i);

        //real analysis here
        newtree->Fill();
    }
    newtree->Print();
    newtree->Write();
    newfile->Close();

}
