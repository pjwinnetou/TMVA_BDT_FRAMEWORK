#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include "TROOT.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include <ctime>
long _ts;

bool BDTClassifier_Function(bool IDvar = false, bool MoreVar = true, bool IDonly= true ){
  
  
  //**********************************************
  //Load Library
  //**********************************************
  TMVA::Tools::Instance();
  //**********************************************

  
  
  //**********************************************
  //THESE ARE JUST FOR TIME STEP AND LOG. YOU CAN REMOVE THIS
  //**********************************************
  std::time_t tstamp = std::time(nullptr);
  _ts = (long) tstamp;
  std::cout <<"time stamp---> " <<  tstamp << std::endl;
  std::system(Form("cat BDTClassifier.C >> ./.past_source/_BDTClassifier_%ld.old",(long) tstamp));
  ofstream log;
  log.open("BDT_description.log", std::ios_base::out|std::ios_base::app);
  log << tstamp;
  char logbuf[2000];
  std::cout << "Write down description for this run :";
  std::cin.getline(logbuf,2000);
  log << ", " << logbuf << std::endl;
  //**********************************************

  TString mainDIR= gSystem->ExpandPathName(gSystem->pwd());
  TString BDTDir = mainDIR + ("/BDTResult");
  void* dirp = gSystem->OpenDirectory(BDTDir.Data());
  if(dirp) gSystem->FreeDirectory(dirp);
  else gSystem->mkdir(BDTDir.Data(),kTRUE);
  
  
  
  //**********************************************
  //INPUT & OUTPUT Call : input files for signal and background
  //**********************************************
  TFile* inputDATA = new TFile("INPUT_BACKGROUND.root","read");
  TFile* inputMC   = new TFile("INPUT_SIGNAL.root","read");
  TFile* output    = new TFile(Form("%s/BDTresult_%ld_IDv%d_MoreVar%d.root",BDTDir.Data(),(long) tstamp, (int) IDvar, (int) MoreVar),"recreate");
  //**********************************************
	
  
  //**********************************************
  //PRESELECTION
  //**********************************************
  string CutId = "CUT_FOR_PARTCLES"; // could be tracks for instance
  string rejectNAN = "&&!TMath::IsNaN(ctau)&&!TMath::IsNaN(ctau3D)&&!TMath::IsNaN(cosAlpha)&&!TMath::IsNaN(cosAlpha3D)"; // to reject NAN values of variables -- put here your training variables 
  TCut cut1 = Form("%s%s",CutId.c_str(),rejectNAN.c_str());
  TCut cut2 = Form("%s%s",CutId.c_str(),rejectNAN.c_str());
  //**********************************************

  
  //**********************************************
  //Assign signal and background tree + training / testing tree
  //**********************************************
  TMVA::DataLoader *loader = new TMVA::DataLoader("dataset");
  TTree* SigTree =(TTree*) inputMC->Get("tree"); // name of the tree in the input signal file
  TTree* BkgTreeTest =(TTree*) inputDATA->Get("tree"); // name of the tree in the input background file <-- suggest to make them consistent
  TTree* BkgTreeTrain = BkgTreeTest->CopyTree("IN_CASE_YOU_NEED_ADDITIONAL_CUTS_PUT_HERE");
  std::cout << "Number of Events in Trees (Sig, BkgTest, BkgTrain) : ( " << SigTree->GetEntries(cut1) << ", "<< BkgTreeTest->GetEntries(cut2) << ", " << BkgTreeTrain->GetEntries(cut2) << " )" << std::endl;
  //**********************************************
  
  
  //**********************************************
  //Factory Call
  //**********************************************
  TMVA::Factory *factory = new TMVA::Factory("TMVA_BDT_Classifier", output, "!V:Silent:Color:DrawProgressBar:Transformations=I;D;P;G;D:AnalysisType=Classification");
  loader->AddVariable("QQMassErr", "Dimu Mass error", "F");
  loader->AddVariable("QQVtxProb", "Vtx prob", "F");
  loader->AddVariable("QQdca", "QQdca", "F");
  loader->AddVariable("cosAlpha", "cos alpha for trajectory angle", "F");
  loader->AddVariable("cosAlpha3D", "cos alpha for trajectory angle 3D", "F");
  //**********************************************

  
  //**********************************************
  //Spectator Call, Will NOT Use For Training 
  //**********************************************
  loader->AddSpectator("cBin","Centrality times 2 bin", "F");
  loader->AddSpectator("mass","Dimuon mass", "F");
  loader->AddSpectator("pt","Dimuon pt","F");
  loader->AddSpectator("y","Dimuon y",  "F");
  //**********************************************


  //**********************************************
  //WEIGHTS -- Careful in the application
  //**********************************************
  
  
  //**********************************************
  // In case you want to apply constant weight for all events (barely the case)
  //**********************************************
  Double_t signalWeight     = 1.0;
  Double_t backgroundWeight = 1.0;
  loader->AddSignalTree     (SigTree, signalWeight);
  loader->AddBackgroundTree (BkgTreeTrain, backgroundWeight);

  //Assign weights from the branch inside the tree : you NEED to have the weight branch already in your input tree 
  //Then apply the weight to signal or (and) background
  loader->SetSignalWeightExpression("weight");
  //loader->SetBackgroundWeightExpression("weight");
  //**********************************************


  //**********************************************
  //Preselection Cut -> Conventional Kinematics -- You can choose the number of signal and background events to be trained here 
  //**********************************************
  loader->PrepareTrainingAndTestTree( cut1, cut2, "nTrain_Signal=200000:nTrain_Background=200000:SplitMode=Random:Random:NormMode=NumEvents:!V");
  //**********************************************


  //**********************************************
  //Book Training BDT Method -- MOST IMPORTANT PART!!!!!!
  //NTrees -- # of trees to be used in the BDT training -- typical suggestion is 200-500 
  //**********************************************
  factory->BookMethod( loader, TMVA::Types::kBDT, TString::Format("BDT_train_%ld", (long) tstamp ),
  	"!H:!V:NTrees=200:MaxDepth=4:MinNodeSize=5%:BoostType=AdaBoost:AdaBoostBeta=0.6:UseBaggedBoost:SeparationType=GiniIndex:PruneMethod=CostComplexity:PruneStrength=1:PruningValFraction=0.3:UseRandomisedTrees=True:UseNvars=2:BaggedSampleFraction=0.4:nCuts=20000:CreateMVAPdfs");
  //**********************************************

  
  //**********************************************
  //Train Test Evaluate
  //**********************************************
  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();
  auto c1 = factory->GetROCCurve(loader);
  c1->Draw();
  c1->SaveAs(Form("%s/ROC_%ld_ID%d_MV%d.pdf",BDTDir.Data(), (long) tstamp,(int) IDvar,(int) MoreVar));
  output->Close();
  //**********************************************

  ofstream out;
  out.open("timestamp.tmp");
  out << tstamp ;
  out.close();
  return true;
}

//Main Function
void BDTClassifier(bool _IDvar = true, bool _MoreVar = false, bool _IDonly= false ){
  bool res = BDTClassifier_Function(_IDvar, _MoreVar, _IDonly);
  if(!res)std::system(Form("rm ./.past_source/_BDTClassifier_%ld.old",(long) _ts));
//  BDTClassifier_Function(true, false, false);
//  BDTClassifier_Function(false, false, false);
}
