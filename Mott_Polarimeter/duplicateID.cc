void duplicateID(){

  TFile* mottFile = new TFile("MottSim.root");
  TTree* mottTree = (TTree*) mottFile->Get("Mott");

  Int_t nentries = GetEntriesFast();
  
  Int_t event, event1, event2;
  
  mottTree->SetBranchStatus("*",0);
  mottTree->SetBranchStatus("Event_ID",1);
  mottTree->SetBranchAddress("Event_ID", &event);
  
  for(Int_t i=0; i<nentries; i++) {
    mottTree->GetEntry(i);
    event1 = event;
    mottTree->GetEntry(i+1);
    event2 = event;
    if (event1==event2) cout << event1 << endl;
  }

  return;
}
