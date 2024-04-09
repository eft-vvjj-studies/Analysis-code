#!/usr/bin/env python3

import ROOT
import numpy as np

CollectHistos = False

ROOT.TH1.AddDirectory(False)

regions = [ "MjjCut" ]
channels = [ "ee", "mm", "all" ]
channels = [ "all" ]
samples = [ "ssWWjjEWLL", "ssWWjjEWLT", "ssWWjjEWTT", "ssWWjjQCD" ]
signalSample = "ssWWjjEWLL"

meas = ROOT.RooStats.HistFactory.Measurement("ssWWjjSnowmassMeas", "ssWWjj Snowmass measurement")
meas.SetOutputFilePrefix("output/out")
meas.SetPOI("mu")

meas.AddConstantParam("Lumi")
meas.SetLumi( 1.0 )
meas.SetLumiRelErr( 0.02 )

meas.SetExportOnly(False)

inFile = ROOT.TFile.Open("output.root")

varName = "Mll"

import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("Logger")

inFileRebinnedName = "output_rebinnned.root"
inFileRebinned = ROOT.TFile.Open(inFileRebinnedName, "recreate")

for k in inFile.GetListOfKeys():
  oName = k.GetName()
  if oName.find(varName+"_MjjCut") == 0:
    histo = k.ReadObj()
    histo = histo.Rebin(4, oName, np.array([0., 100., 200., 500., 1000.]))
    inFileRebinned.cd()
    histo.Write()

inFileRebinned.Close()
inFile.Close()

inFile = ROOT.TFile.Open("output.root")

def getHistName(varName, regName, sampleName, chName):
  return "_".join([varName, regName, sampleName, chName])

def getHist(varName, regName, sampleName, chName):
  hName = "_".join([varName, regName, sampleName, chName])
  inFile.cd()
  inFile.ls()
  logger.info(f"Getting histogram {hName}")
  histo = inFile.Get(hName)
  histo.SetDirectory(0)
  histo = histo.Rebin(4, hName, np.array([0., 100., 200., 500., 1000.]))
  if not histo: raise ValueError(f"Cannot find histogram {hName}")
  return histo

def makeSample(regName, sampleName, chName):
  sName = "_".join([varName, regName, sampleName, chName])
  logger.info(f"Building sample {sName}")
  if not CollectHistos:
    sample = ROOT.RooStats.HistFactory.Sample(sName)
    nomHist = getHist(varName, regName, sampleName, chName)
    sample.SetHisto(nomHist)
  else: sample = ROOT.RooStats.HistFactory.Sample( sName, getHistName(varName, regName, sampleName, chName), inFileRebinnedName )
  if sampleName == signalSample:
    sample.AddNormFactor( "mu", 1, 0, 3 )
  else:
    pass
    sample.ActivateStatError()


  return sample

def makeChannel(regName, chName):
  cName = f"{regName}_{chName}"
  logger.info(f"Building channel {cName}")
  chan  = ROOT.RooStats.HistFactory.Channel(cName)

  if not CollectHistos: chan.SetData( getHist(varName, regName, "ssWWjjEW", chName) )
  else: chan.SetData( getHistName(varName, regName, "ssWWjjEW", chName), inFileRebinnedName )
  chan.SetStatErrorConfig( 0.05, "Poisson" )

  for sampleName in samples:
    chan.AddSample( makeSample(regName, sampleName, chName) )
  return chan

# Looping over each region and channel to create fit channels

for regName in regions:
  for chName in channels:
    meas.AddChannel( makeChannel(regName, chName) )

if CollectHistos: meas.CollectHistograms()
meas.PrintTree()

ws = ROOT.RooStats.HistFactory.MakeModelAndMeasurementFast( meas )

ws.Print()

oF = ROOT.TFile("output/out_combined_meas_model.root", "update")
ws.Write(ws.GetName(), ROOT.TObject.kOverwrite)
oF.Close()

#import sys
#sys.exit(0)

mu = ws.var("mu")
pdf = ws.pdf("simPdf")
sbModel = ws.obj("ModelConfig")
data = ws.data("obsData")
asimov = ws.data("asimovData")

ROOT.Math.MinimizerOptions.SetDefaultMinimizer("Minuit2")
ROOT.Math.MinimizerOptions.SetDefaultStrategy(0)
ROOT.Math.MinimizerOptions.SetDefaultPrintLevel(1)

params = ROOT.RooArgSet( sbModel.GetNuisanceParameters(), sbModel.GetParametersOfInterest() )

theData = asimov

nll = pdf.createNLL( theData, ROOT.RooFit.Constrain(params), ROOT.RooFit.GlobalObservables(sbModel.GetGlobalObservables()), ROOT.RooFit.Offset(1) )


def minimize(fcn, save = False, retry_mode = 3):
  printLevel = ROOT.Math.MinimizerOptions.DefaultPrintLevel()
  msgLevel = ROOT.RooMsgService.instance().globalKillBelow()
  if printLevel < 0:
      ROOT.RooMsgService.instance().globalKillBelow(ROOT.RooFit.FATAL)
  
  strategy = ROOT.Math.MinimizerOptions.DefaultStrategy()
  save_def_strategy = strategy

  minimizer = ROOT.RooMinimizer(fcn)
  minimizer.optimizeConst(2)
  minimizer.setStrategy(strategy)
  minimizer.setPrintLevel(printLevel)
  minimizer.setMinimizerType(ROOT.Math.MinimizerOptions.DefaultMinimizerType())

  status = minimizer.minimize( ROOT.Math.MinimizerOptions.DefaultMinimizerType(), ROOT.Math.MinimizerOptions.DefaultMinimizerAlgo() )

  # Possibly re-trying if the fit didn't work
  if retry_mode == 0:
    if status != 0 and status != 1 and strategy < 2:
      strategy += 1
      logger.warning( f"Fit failed with status {status}. Retrying with strategy {strategy}" )
      minimizer.setStrategy(strategy)
      status = minimizer.minimize( ROOT.Math.MinimizerOptions.DefaultMinimizerType(), ROOT.Math.MinimizerOptions.DefaultMinimizerAlgo() )

    if status != 0 and status != 1 and strategy < 2:
      strategy += 1
      logger.warning( f"Fit failed with status {status}. Retrying with strategy {strategy}" )
      minimizer.setStrategy(strategy)
      status = minimizer.minimize( ROOT.Math.MinimizerOptions.DefaultMinimizerType(), ROOT.Math.MinimizerOptions.DefaultMinimizerAlgo() )

  else:
    for i in range(retry_mode):
      if status == 0 or status == 1: break
      logger.warning( f"Fit failed with status {status}. Retrying with strategy {strategy}." )
      minimizer.setStrategy(strategy)
      status = minimizer.minimize( ROOT.Math.MinimizerOptions.DefaultMinimizerType(), ROOT.Math.MinimizerOptions.DefaultMinimizerAlgo() )

  if printLevel < 0:
    ROOT.RooMsgService.insurance().setGlobalKillerBelow(msgLevel)
  ROOT.Math.MinimizerOptions.SetDefaultStrategy(save_def_strategy)

  if save:
    fitRes = minimizer.save( f"fitresult_{fcn.GetName()}", f"fitresult_{fcn.GetName()}" )

  return fitRes


res = minimize(nll, True, 0)

res.SaveAs("result.root")
