#!/usr/bin/env python3

import ROOT
from ROOT import TChain, RDataFrame, TCanvas
ROOT.gROOT.SetBatch(True)

ROOT.TH1.SetDefaultSumw2()

#ROOT::VecOps::RVec<float>(0.10566, 0.10566)

# Defining custom function to use in PyROOT
# Todo: add to arrays (see below)

#Ideally this should return the indices that are the specific value.
ElementSel_code = '''
ROOT::VecOps::RVec<int> ElementSel(ROOT::VecOps::RVec<float> v, float value)
{
ROOT::VecOps::RVec<int> r;
const auto size = v.size();
r.reserve(size);
for(int i=0; i<size; i++) {
      if(v[i] == value) {
         r.emplace_back(i);
      }
   }
   return r;
}
'''
ROOT.gInterpreter.Declare(ElementSel_code)

LeptonMass_code = '''
ROOT::VecOps::RVec<float> LeptonMass(ROOT::VecOps::RVec<float> e_pt, ROOT::VecOps::RVec<float> m_pt,
 ROOT::VecOps::RVec<float> l_pt)
{
  ROOT::VecOps::RVec<float> l_m(l_pt.size());
  for (int i=0; i < l_pt.size(); ++i) {
    for (int j=0; j < std::max(e_pt.size(), m_pt.size()); ++j) {
      if(j < e_pt.size() && l_pt[i] == e_pt[j]) {
        l_m[i] = 0.000511;
        break;
      } else if (j < m_pt.size() && l_pt[i] == m_pt[j]) {
        l_m[i] = 0.10566;
        break;
      } else {
        l_m[i] = -1;
      }
    }
  }
  return l_m;
}
'''
ROOT.gInterpreter.Declare(LeptonMass_code)


ComputeInvariantMass_code = '''
double ComputeInvariantMass(ROOT::VecOps::RVec<float> pt, ROOT::VecOps::RVec<float> eta, ROOT::VecOps::RVec<float> phi,
 ROOT::VecOps::RVec<float> m)
{
  if(2 > pt.size()) return std::numeric_limits<double>::min();
  TLorentzVector p[2];
  for(int i : {0, 1}) p[i].SetPtEtaPhiM(pt[i], eta[i], phi[i], m[i]);
  return (p[0]+p[1]).M();
}
'''

ROOT.gInterpreter.Declare(ComputeInvariantMass_code)


ComputeInvariantMassLeptons_code = '''
double ComputeInvariantMassLeptons(
  ROOT::VecOps::RVec<float> l_pt, ROOT::VecOps::RVec<float> l_eta, ROOT::VecOps::RVec<float> l_phi,
  ROOT::VecOps::RVec<float> l_m
) {
  if(2 <= l_pt.size()) {
    return ComputeInvariantMass(l_pt, l_eta, l_phi, l_m);
  } else {
    return std::numeric_limits<double>::min();
  }
}
'''
ROOT.gInterpreter.Declare(ComputeInvariantMassLeptons_code)


ComputeDPhijj_code = '''
double ComputeDPhijj(
  ROOT::VecOps::RVec<float> jet_phi
) {
  if(2 <= jet_phi.size()) {
    Double_t dPhi = fabs(jet_phi[0] - jet_phi[1]);
    while(true) {
      if(dPhi < TMath::Pi()) return dPhi;
      dPhi = fabs(dPhi - 2.0 * TMath::Pi());
    }
  } else {
    return std::numeric_limits<double>::lowest();
  }
}
'''
ROOT.gInterpreter.Declare(ComputeDPhijj_code)

ComputeDRjl_code = '''                                                                                                                                                                                                                          
double ComputeDRjl(
  ROOT::VecOps::RVec<float> jet_phi, ROOT::VecOps::RVec<float> jet_eta, ROOT::VecOps::RVec<float> lep_phi, ROOT::VecOps::RVec<float> lep_eta, int x, int y
) {
  if(jet_phi.size() >=2 && lep_phi.size() >=2) {
    Double_t dPhilj = fabs(jet_phi[x] - lep_phi[y]);
    Double_t dEtalj = fabs(jet_eta[x] - lep_eta[y]);
    while(true) {
      if(dPhilj < TMath::Pi()) return sqrt(dPhilj * dPhilj + dEtalj * dEtalj);
      dPhilj = fabs(dPhilj - 2.0 * TMath::Pi());
    }
  }
  else {
    return std::numeric_limits<double>::lowest();
  }
}
'''
ROOT.gInterpreter.Declare(ComputeDRjl_code)


ComputeDRjj_code = '''
double ComputeDRjj(
  ROOT::VecOps::RVec<float> jet_phi, ROOT::VecOps::RVec<float> jet_eta
) {
  if(jet_phi.size() >= 2) {
    Double_t dPhijj = fabs(jet_phi[0] - jet_phi[1]);
    Double_t dEtajj = fabs(jet_eta[0] - jet_eta[1]); 
    while(true) {
      if(dPhijj < TMath::Pi()) return sqrt(dPhijj * dPhijj + dEtajj * dEtajj);
      dPhijj = fabs(dPhijj - 2.0 * TMath::Pi());
    }
  }
  else {
    return std::numeric_limits<double>::lowest();
  }
}
'''
ROOT.gInterpreter.Declare(ComputeDRjj_code)

#This Code should return the reconcstructed W or Z mass aswell as any other parameters if needed 
#int type, 0 = return W 4 vec, 1 = return Z 4 vec
ComputeWZMass_code = '''
TLorentzVector ComputeWZMass(int type,int nleptons, int nelectrons , int nmuons, ROOT::VecOps::RVec<float> L_pt, 
ROOT::VecOps::RVec<float> L_eta, ROOT::VecOps::RVec<float> L_phi, ROOT::VecOps::RVec<float> L_m, float Met_met, 
float Met_eta, float Met_phi
){
  int result_type = pow(2,type);
  if(nleptons == 3){
    // Same flavor leptons
    if(abs(nelectrons - nmuons) == 3){
      TLorentzVector v1;
      return v1;
    }
    //Different flavor leptons
    else if (abs(nelectrons - nmuons) == 1){
      if(nelectrons > nmuons){
        switch(result_type){
          case 1: {
            TLorentzVector p;
            TLorentzVector m;
            auto muon_idx = ElementSel(L_m,0.10566);
            p.SetPtEtaPhiM(L_pt[muon_idx[0]],L_eta[muon_idx[0]],L_phi[muon_idx[0]],L_m[muon_idx[0]]);
            m.SetPtEtaPhiM(Met_met, Met_eta, Met_phi, 0.0);
            return(p + m);
            break; 
          }
          case 2:{
            auto electron_idx = ElementSel(L_m, 0.000511);
            TLorentzVector p[2];
            for(int i : {0, 1}) p[i].SetPtEtaPhiM(L_pt[electron_idx[i]],L_eta[electron_idx[i]],L_phi[electron_idx[i]],L_m[electron_idx[i]]);
            return (p[0]+p[1]);
            break;
          }
        }
      }
      else if (nmuons > nelectrons){
        switch(result_type){
          case 1:{
            TLorentzVector p;
            TLorentzVector m;
            auto electron_idx = ElementSel(L_m,0.000511);
            p.SetPtEtaPhiM(L_pt[electron_idx[0]],L_eta[electron_idx[0]],L_phi[electron_idx[0]],L_m[electron_idx[0]]);
            m.SetPtEtaPhiM(Met_met, Met_eta, Met_phi, 0.0);
            return(p + m);
            break;
          }
          case 2:{
            auto muon_idx = ElementSel(L_m, 0.10566);
            TLorentzVector p[2];
            for(int i : {0, 1}) p[i].SetPtEtaPhiM(L_pt[muon_idx[i]],L_eta[muon_idx[i]],L_phi[muon_idx[i]],L_m[muon_idx[i]]);
            return (p[0]+p[1]);
            break;
          }
        }
      }
    }
  }
  else {
    TLorentzVector v1;
    return v1;
  }
}
'''
ROOT.gInterpreter.Declare(ComputeWZMass_code)

samples = [
#  ("ssWWjjEW", "Samples/PROC_100TeV_ssWWjj_EW/*root", "Delphes"),
## ("ssWWjjEWLL", "Samples/PROC_100TeV_ssWWjj_EW_LL/*root", "Delphes"),
## ("ssWWjjEWLT", "Samples/PROC_100TeV_ssWWjj_EW_LT/*root", "Delphes"),
## ("ssWWjjEWTT", "Samples/PROC_100TeV_ssWWjj_EW_TT/*root", "Delphes"),
##  ("ssWWjjEWLLWWcmf", "Samples/PROC_100TeV_ssWWjj_EW_WWcmf_LL/*root", "Delphes"),
##  ("ssWWjjEWLTWWcmf", "Samples/PROC_100TeV_ssWWjj_EW_WWcmf_LT/*root", "Delphes"),
##  ("ssWWjjEWTTWWcmf", "Samples/PROC_100TeV_ssWWjj_EW_WWcmf_TT/*root", "Delphes"),
##  ("ssWWjjQCD", "Samples/PROC_100TeV_ssWWjj_QCD/*root", "Delphes"),
##  ("tZq3l", "Samples/PROC_100TeV_tZq_3l/*root", "Delphes"),
##  ("WlljjEW4F", "Samples/PROC_100TeV_Wlljj_EW_4F/*root", "Delphes"),
 ## ("WlljjQCD2jnomatch", "Samples/PROC_100TeV_Wlljj_QCD_2j_nomatch/*root", "Delphes"),
#("ssWWSM","Samples/ssWW_SM/events.root","Delphes"),
("WZQuad","Samples/sigWZ/WZquad/T2/events_WZquad.root","Delphes"),

]
# Units: pb
cross_section = {
  "ssWWjjEW": 0.7778,
  "ssWWjjEWLL": 0.04538,
  "ssWWjjEWLT": 0.26998,
  "ssWWjjEWTT": 0.46618,
  "ssWWjjEWLLWWcmf": 0.08144,
  "ssWWjjEWLTWWcmf": 0.2781,
  "ssWWjjEWTTWWcmf": 0.5067,
  "ssWWjjQCD": 1.232,
  "tZq3l": 1.115,
  "WlljjEW4F": 0.4533,
  "WlljjQCD2jnomatch": 6.9,
  "ssWWSM":1,
  "WZQuad":1,
}

# Assuming 30 ab-1
lumi = 1# pb-1

# Adding these definitions to provoke copying the values (avoids the non-split branch warning), ultimately more efficient
# This should be used to improve performance if needed, and a new TTree should be created
aliases = [

  ("el_pt", "Electron.PT"),
  ("el_eta", "Electron.Eta"),
  ("el_phi", "Electron.Phi"),
  ("mu_pt", "Muon.PT"),
  ("mu_eta", "Muon.Eta"),
  ("mu_phi", "Muon.Phi"),
  ("jet_pt", "Jet.PT"),
  ("jet_eta", "Jet.Eta"),
  ("jet_phi", "Jet.Phi"),
  ("jet_mass", "Jet.Mass"),

]

definitions = [
  ("argsMuons", "ROOT::VecOps::Reverse(ROOT::VecOps::Argsort(Muon.PT));"),
  ("argsElectron", "ROOT::VecOps::Reverse(ROOT::VecOps::Argsort(Electron.PT));"),
  ("leptons", "ROOT::VecOps::Concatenate(Electron.PT, Muon.PT);"),
  ("argLeptons", "ROOT::VecOps::Reverse(ROOT::VecOps::Argsort(leptons));"),
  ("argJets", "ROOT::VecOps::Reverse(ROOT::VecOps::Argsort(Jet.PT));"),


  ("muonPt", "ROOT::VecOps::Take(Muon.PT, argsMuons);"),
  ("electronPt", "ROOT::VecOps::Take(Electron.PT, argsElectron);"),
  ("leptonPt", "ROOT::VecOps::Take(leptons, argLeptons);"),
  ("jetPt", " ROOT::VecOps::Take(Jet.PT, argJets);"),

  ("met","MissingET.MET[0]"),

  ("leptonMass", "LeptonMass(electronPt, muonPt, leptonPt);"),
  ("leptonEta", "ROOT::VecOps::Take(ROOT::VecOps::Concatenate(Electron.Eta, Muon.Eta), argLeptons);"),
  ("leptonPhi", "ROOT::VecOps::Take(ROOT::VecOps::Concatenate(Electron.Phi, Muon.Phi), argLeptons);"),
  ("leptonCharge", "ROOT::VecOps::Take(ROOT::VecOps::Concatenate(Electron.Charge, Muon.Charge), argLeptons);"),
  ("jetPhi", "ROOT::VecOps::Take(Jet.Phi, argJets);"),
  ("jetEta", "ROOT::VecOps::Take(Jet.Eta, argJets);"),
  ("jetMass", "ROOT::VecOps::Take(Jet.Mass, argJets);"),

  ("goodElectrons", "electronPt > 15"),
  ("goodMuons", "muonPt > 15"),
  ("goodLeptons", "leptonPt > 15"),
  ("nLeptons", "Sum(goodElectrons) + Sum(goodMuons)"),
  ("goodJets", "jetPt >= 50"),
  ("nJets", "Sum(goodJets)"),

  #("gd_el_pt", "el_pt[goodElectrons]"),
  #("gd_el_eta", "el_pt[goodElectrons]"),
  #("gd_el_phi", "el_phi[goodElectrons]"),
  #("gd_mu_pt", "mu_pt[goodMuons]"),
  #("gd_mu_eta", "mu_pt[goodMuons]"),
  #("gd_mu_phi", "mu_phi[goodMuons]"),
  #("gd_jet_pt","jet_pt[goodJets]"),
  #("gd_jet_eta","jet_eta[goodJets]"),
  #("gd_jet_phi","jet_phi[goodJets]"),
  #("gd_jet_mass","jet_mass[goodJets]"),
  #("mll", "ComputeInvariantMassLeptons( gd_el_pt, gd_el_eta, gd_el_phi, gd_mu_pt, gd_mu_eta, gd_mu_phi )"),
  #("mjj", "ComputeInvariantMass( gd_jet_pt, gd_jet_eta, gd_jet_phi, gd_jet_mass )"),

  ("mll", "ComputeInvariantMassLeptons(leptonPt[goodLeptons], leptonEta[goodLeptons],"
          " leptonPhi[goodLeptons], leptonMass[goodLeptons])"),

  ("mjj", "ComputeInvariantMass(jetPt[goodJets], jetEta[goodJets], jetPhi[goodJets], jetMass[goodJets])"),
  ("dPhijj", "ComputeDPhijj(jetPhi[goodJets])"),
  ("dEtajj", "jetPhi[goodJets].size() >= 2 ? std::fabs(jetEta[goodJets][0] - jetEta[goodJets][1]) "
             ": std::numeric_limits<float>::lowest()"),
  ("pTj1", "jetPt[goodJets].size() >= 2 ? jetPt[goodJets][0] : std::numeric_limits<double>::lowest();"),
  ("pTj2", "jetPt[goodJets].size() >= 2 ? jetPt[goodJets][1] : std::numeric_limits<double>::lowest();"),
  ("pTj3", "jetPt[goodJets].size() >= 3 ? jetPt[goodJets][2] : std::numeric_limits<double>::lowest();"),
  ("pTl1", "leptonPt[goodLeptons].size() >= 2 ? leptonPt[goodLeptons][0] : "
            "std::numeric_limits<double>::lowest();"),
  ("pTl2", "leptonPt[goodLeptons].size() >= 2 ? leptonPt[goodLeptons][1] : "
            "std::numeric_limits<double>::lowest();"),
  ("etal1", "leptonEta[goodLeptons].size() >= 2 ? leptonEta[goodLeptons][0] : "
              "std::numeric_limits<double>::lowest();"),
  ("etal2", "leptonEta[goodLeptons].size() >= 2 ? leptonEta[goodLeptons][1] : "
              "std::numeric_limits<double>::lowest();"),
  ("ptrel", "(pTl1*pTl2)/(pTj1*pTj2)"),
  ("dRj1l1", "ComputeDRjl(jetPhi[goodJets], jetEta[goodJets], leptonPhi[goodLeptons], leptonEta[goodLeptons], 0, 0)"),
  ("dRj2l2", "ComputeDRjl(jetPhi[goodJets], jetEta[goodJets], leptonPhi[goodLeptons], leptonEta[goodLeptons], 1, 1)"),
  ("dRjj", "ComputeDRjj(jetPhi[goodJets], jetEta[goodJets])"),
  ("Wmass", "ComputeWZMass(0, nLeptons,Sum(goodElectrons), Sum(goodMuons),leptonPt[goodLeptons], leptonEta[goodLeptons],"
          " leptonPhi[goodLeptons], leptonMass[goodLeptons], MissingET.MET[0], MissingET.Eta[0],MissingET.Phi[0])"),
  ("Zmass", "ComputeWZMass(1, nLeptons,Sum(goodElectrons), Sum(goodMuons),leptonPt[goodLeptons], leptonEta[goodLeptons],"
          " leptonPhi[goodLeptons], leptonMass[goodLeptons], MissingET.MET[0], MissingET.Eta[0],MissingET.Phi[0])")

]
#("dPhij1l1", "jetPhi[goodJets].size() >= 2 ? leptonPhi[goodLeptons].size() >= 2 ? jetPhi[goodJets][0] - leptonPhi[goodLeptons][0]"),
#("dEtaj1l1", "jetEta[goodJets].size() >= 2 ? leptonEta[goodLeptons].size() >= 2 ? jetEta[goodJets][0] - leptonEta[goodLeptons][0]")

#I'm pretty sure filters are applied sequentially so make sure that you don't have two mutually exclusive filters -WS
filters = [
  ("nLeptons == 3", "LeptonCut"),
  #("Wmass == 1", "CodeCheck"), #From me checking that the Caveman code works
  #("leptonCharge[goodLeptons][0] * leptonCharge[goodLeptons][1] > 0", "SameSignCut"),
  ("nJets >= 2", "JetCut"),
  ("mll >= 60", "MllCut"),
  ("mjj >= 500", "MjjCut"),
  ("dRjj > 2.5", "dRjjCut"),
  ("dRj1l1 > 0.4 && dRj2l2 > 0.4", "dRllCut"),
  ("Sum(MissingET.MET) >= 50.0", "METCut"),
]

# Regions: these are defined based on a certain cut level
regions = [
  ("METCut", ("Sum(goodElectrons)==2", "ee")),
  ("METCut", ("Sum(goodMuons)==2", "mm")),
  ("METCut", ("Sum(goodElectrons)==1 && Sum(goodMuons)==1", "emme")),
]

# These are the definitions of the histograms
# ((name, title, nbins, minx, miny), variable)
# Todo: expand to multiple dimensions

histograms = [
#  (("njets", "nJets", 10, 0, 10), "nJets"),
#  (("nleptons", "nLeptons", 10, 0, 10), "nLeptons"),
  (("mll", "Mll", 50, 0, 5000.), "mll"),
  (("mjj", "Mjj", 25, 0, 10000.), "mjj"),
  (("dPhijj", "DPhijj", 16, 0, 3.14159), "dPhijj"),
  (("dEtajj", "DEtajj", 60, 0, 15.), "dEtajj"),
  (("pTj1", "pTj1", 50, 0, 500.), "pTj1"),
  (("pTj2", "pTj2", 50, 0, 500.), "pTj2"),
  (("pTl1", "pTl1", 50, 0, 500.), "pTl1"),
  (("pTl2", "pTl2", 50, 0, 500.), "pTl2"),
  (("etal1", "etal1", 60, -6.0, 6.0), "etal1"),
  (("etal2", "etal2", 60, -6.0, 6.0), "etal2"),
  (("met", "met", 100, 0.0, 1000.0), "met"),
  (("dPhijj_vs_mll", "dPhijj_vs_mll", 16, 0, 3.14159, 20, 0, 1000), "dPhijj", "mll"),
  (("ptrel", "PtRel", 20, 0, 1.), "ptrel"),
  (("dRj1l1", "DRl1j1", 40, 0, 10.), "dRj1l1"),
  (("dRj2l2", "DRl2j2", 40, 0, 10.), "dRj2l2"),
  (("dRjj", "DRjj", 40, 0, 10.), "dRjj"),
]


def ana():
  output_root_file = "output.root"

  import os
  # It's unlikely that you have '/DelphesAnalysis/' unless you're in the container
  # If you do, then it's still OK as likely that's what you want
  inContainer = os.path.isdir("/DelphesAnalysis/")
  if inContainer: os.chdir("/DelphesAnalysis/")

  oF = ROOT.TFile.Open(output_root_file, "RECREATE")

  for s in samples:
    sampleName = s[0]
    sampleGlob = s[1]

    print(f"Processing sample {sampleName} at {s[1]}")
    t = TChain(s[2])
    t.Add(sampleGlob)
    print(f"TChain created with {t.GetEntries()} entries")
    if t.GetEntries() == 0: continue

    # WARNING: assumes input data is NOT filtered
    tot_events = t.GetEntries()
  
    rd = RDataFrame(t)
    # rdNodes contains pointers to the various nodes (currently used only for filters)
    rdNodes = { "nofilter": rd }
  
  
    # Consider whether to move the function elsewhere, or make a class (better)
    def make_histos(rd, regName, chName = "all"):
      for h in histograms:
        print(h)
        # Todo: handle different sizes
        if len(h) == 2: histo = rd.Histo1D(*h)
        elif len(h) == 3: histo = rd.Histo2D(*h)
        histo.Scale( cross_section[sampleName] / tot_events  * lumi )
        varName = h[0][1]
        hName = "_".join([varName, regName, sampleName, chName])
        histo.SetName(hName)
        nBins = histo.GetNbinsX()
        entries = histo.GetEntries()
        # Handling overflow bin ... to lazy for underflow ; is it needed ?
        lastBinContent = histo.GetBinContent(nBins) + histo.GetBinContent(nBins+1)
        lastBinError = (histo.GetBinError(nBins)**2 + histo.GetBinError(nBins+1)**2)**0.5
        histo.SetBinContent(nBins, lastBinContent)
        histo.SetBinContent(nBins+1, 0.)
        histo.SetBinError(nBins, lastBinError)
        histo.SetBinError(nBins+1, 0.)
        histo.SetEntries(entries)
        c = TCanvas(hName)
        histo.Draw()
        #c.Print(".png")
        oF.cd()
        histo.Write()
  
    if False:
      for a in aliases:
        print(a)
        rd = rd.Define(*a)
  
    for d in definitions:
      print(d)
      rd = rd.Define(*d)
  
    # cols = ROOT.vector('string')(["nLeptons", "nElectrons", "nMuons"])
    # rd.Display(cols).Print()
  
    make_histos(rd, regName = "nofilter")
    for f in filters:
      print(f)
      rd = rd.Filter(*f)
      rdNodes[f[1]] = rd
      make_histos(rd, regName = f[1])

    for r in regions:
      print(r)
      rd = rdNodes[r[0]]
      rd = rd.Filter(*r[1])
      rdNodes[r[1][1]] = rd
      make_histos(rd, regName = r[0], chName = r[1][1])

    rd.Report().Print()

  oF.ls()
  oF.Close()

if __name__ == '__main__':
  ana()

