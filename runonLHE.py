#!/usr/bin/env python

if __name__ == "__main__":
  #check the inputs before importing mela, which takes a while
  import argparse
  parser = argparse.ArgumentParser()
  parser.add_argument("lhefile", nargs="+")
  group = parser.add_mutually_exclusive_group(required=True)
  group.add_argument("--VBF", action="store_true")
  group.add_argument("--ZH", action="store_true")
  group.add_argument("--WH", action="store_true")
  group.add_argument("--ggH", action="store_true")
  parser.add_argument("--firstevent", type=int, required=True)
  parser.add_argument("--angles", action="store_true")

  args = parser.parse_args()


import array
import os
import re
import ROOT

from lhefile import LHEFile_JHUGenVBFVH
from ZZMatrixElement.MELA.mela import TVar

def runonLHE(lhefile, productionmode, firsteventnumber, doangles=False):
  # constant definitions, YM
  # VBF
  a3 = 0.297979018705
  a2 = 0.27196538
  L1 = -2158.21307286
  L1Zg = -4091.0
  if productionmode == "ZH":
    a3 = 0
    a2 = 0
    L1 = 0
    L1Zg = -642.9534
  if productionmode == "WH":
    a3 = 0
    a2 = 0
    L1 = 0
    L1Zg = 0.0


  rootfile = lhefile.replace(".lhe", ".root")
  inputg1 = inputg2 = inputg4 = inputg1prime2 = 0
  with open(lhefile) as f:
    for line in f:
      if "<init>" in line: break
      line = line.strip()
      if productionmode == "ggH":
        match = re.match("ghg2= *([0-9E+-.]*) *[0-9E+-.]*i", line)
        if match: inputg2 = float(match.group(1))
        match = re.match("ghg4= *([0-9E+-.]*) *[0-9E+-.]*i", line)
        if match: inputg4 = float(match.group(1))
      else:
        match = re.match("ghz1= *([0-9E+-.]*) *[0-9E+-.]*i", line)
        if match: inputg1 = float(match.group(1))
        match = re.match("ghz2= *([0-9E+-.]*) *[0-9E+-.]*i", line)
        if match: inputg2 = float(match.group(1))
        match = re.match("ghz4= *([0-9E+-.]*) *[0-9E+-.]*i", line)
        if match: inputg4 = float(match.group(1))
        match = re.match("ghz1_prime2= *([0-9E+-.]*) *[0-9E+-.]*i", line)
        if match: inputg1prime2 = float(match.group(1))

  print "Here are the couplings (please check):"
  print "For "+lhefile+":"
  print "        g1 =", inputg1
  print "        g2 =", inputg2
  print "        g4 =", inputg4
  print "  g1prime2 =", inputg1prime2


  if os.path.exists(rootfile): raise IOError(rootfile + " already exists")
  rootf = ROOT.TFile(rootfile, "CREATE")
  t = ROOT.TTree("tree", "tree")

  if productionmode <> "ggH":
    branches_float = ["wt_a1", "wt_a2", "wt_a3", "wt_L1", "wt_L1Zg", "wt_a1a2", "wt_a1a3", "wt_a1L1", "wt_a1L1Zg", \
                        "wt_a3int", "wt_a2int", "wt_L1int", "wt_L1Zgint"]
  else:
    branches_float = ["wt_a2", "wt_a3", "wt_a2a3"]
  branches_int = ["eventnumber"]

  if doangles:
    branches_float += ["costheta1", "costheta2", "Phi", "costhetastar", "Phi1"]
    if productionmode in ("VBF", "VH"):
      branches_float += ["q1", "q2"]

  branches = {name: array.array("f", [0]) for name in branches_float}
  branches.update({name: array.array("i", [0]) for name in branches_int})

  for branchname in branches_float:
    t.Branch(branchname, branches[branchname], branchname+"/F")
  for branchname in branches_int:
    t.Branch(branchname, branches[branchname], branchname+"/I")

  with LHEFile_JHUGenVBFVH(lhefile) as lhef:
    for i, event in enumerate(lhef):
      branches["eventnumber"][0] = firsteventnumber+i
      if i>0 and i % 1000 == 0: print "Processed", i, "events"
      if productionmode == "VBF":
        process = TVar.JJVBF
      elif productionmode == "ggH":
        process = TVar.JJQCD
      elif productionmode == "ZH":
        if all(11 <= abs(particle.first) <= 16 for particle in event.associated):
          process = TVar.Lep_ZH
        elif all(1 <= abs(particle.first) <= 5 for particle in event.associated):
          process = TVar.Had_ZH
        else:
          assert False
      elif productionmode == "WH":
        if all(11 <= abs(particle.first) <= 16 for particle in event.associated):
          process = TVar.Lep_WH
        elif all(1 <= abs(particle.first) <= 5 for particle in event.associated):
          process = TVar.Had_WH
        else:
          assert False

      event.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, process)
      event.ghz1 = inputg1
      event.ghz2 = event.ghg2 = inputg2
      event.ghz4 = event.ghg4 = inputg4
      event.ghz1_prime2 = inputg1prime2
      p_base = event.computeProdP(False)

      if productionmode == "ggH":
        event.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, process)
        event.ghg2 = 1
        p_a2 = event.computeProdP(False)

        event.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, process)
        event.ghg4 = 1
        p_a3 = event.computeProdP(False)

        event.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, process)
        event.ghg2 = event.ghg4 = 1
        p_a2a3 = event.computeProdP(False) - p_a2 - p_a3
      else:
        event.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, process)
        event.ghz1 = 1
        p_a1 = event.computeProdP(False)

        event.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, process)
        event.ghz2 = 1
        p_a2 = event.computeProdP(False)

        event.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, process)
        event.ghz4 = 1
        p_a3 = event.computeProdP(False)

        event.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, process)
        event.ghz1_prime2 = 1
        p_L1 = event.computeProdP(False)

        event.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, process)
        event.ghzgs1_prime2 = 1
        p_L1Zg = event.computeProdP(False)

        event.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, process)
        event.ghz1 = event.ghz2 = 1
        p_a1a2 = event.computeProdP(False) - p_a1 - p_a2

        event.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, process)
        event.ghz1 = event.ghz4 = 1
        p_a1a3 = event.computeProdP(False) - p_a1 - p_a3

        event.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, process)
        event.ghz1 = 1
        event.ghz1_prime2 = 1e4  #otherwise you get big + small - big and lose precision
        p_a1L1 = (event.computeProdP(False) - p_a1 - p_L1*1e8) / 1e4

        event.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, process)
        event.ghz1 = 1
        event.ghzgs1_prime2 = 1e4  #otherwise you get big + small - big and lose precision
        p_a1L1Zg = (event.computeProdP(False) - p_a1 - p_L1Zg*1e8) / 1e4

      if productionmode != "ggH":
        branches["wt_a1"][0] = p_a1 / p_base

      branches["wt_a2"][0] = p_a2 / p_base
      branches["wt_a3"][0] = p_a3 / p_base

      if productionmode == "ggH":
        branches["wt_a2a3"] = p_a2a3 / p_base
      else:
        branches["wt_L1"][0] = p_L1 / p_base
        branches["wt_L1Zg"][0] = p_L1Zg / p_base
        branches["wt_a1a2"][0] = p_a1a2 / p_base
        branches["wt_a1a3"][0] = p_a1a3 / p_base
        branches["wt_a1L1"][0] = p_a1L1 / p_base
        branches["wt_a1L1Zg"][0] = p_a1L1Zg / p_base
        branches["wt_a3int"][0] = (p_a1 + a3*p_a1a3 + a3*a3*p_a3) / p_base
        branches["wt_a2int"][0] = (p_a1 + a2*p_a1a2 + a2*a2*p_a2) / p_base
        branches["wt_L1int"][0] = (p_a1 + L1*p_a1L1 + L1*L1*p_L1) / p_base
        branches["wt_L1Zgint"][0] = (p_a1 + L1Zg*p_a1L1Zg + L1Zg*L1Zg*p_L1Zg) / p_base

      if doangles:
        if productionmode in ("VBF", "HJJ"):
          (Q2V1, Q2V2, branches["costheta1"][0], branches["costheta2"][0],
           branches["Phi"][0], branches["costhetastar"][0], branches["Phi1"][0]) = event.computeVBFAngles()
          branches["q1"][0] = Q2V1**.5
          branches["q2"][0] = Q2V2**.5
        elif productionmode in ("ZH", "WH"):
          (branches["costheta1"][0], branches["costheta2"][0],
           branches["Phi"][0], branches["costhetastar"][0], branches["Phi1"][0]) = event.computeVHAngles()

      t.Fill()

    print "Processed", i+1, "events"
    rootf.Write()
    rootf.Close()

if __name__ == "__main__":
  for filename in args.lhefile:
    if args.VBF: productionmode = "VBF"
    if args.ZH: productionmode = "ZH"
    if args.WH: productionmode = "WH"
    if args.ggH: productionmode = "ggH"
    runonLHE(filename, productionmode, args.firstevent, args.angles)

