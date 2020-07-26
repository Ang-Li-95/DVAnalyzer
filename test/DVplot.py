import ROOT as R
import math as M
import argparse
import subprocess


parser = argparse.ArgumentParser()
parser.add_argument('--inputDir',dest='inputDir',default="")
parser.add_argument('--inputeosDir',dest='inputeosDir',default="")
parser.add_argument('--pattern',dest='pattern')
parser.add_argument('--output',dest="output")

args = parser.parse_args()

if (args.inputDir != ""):
  print (args.inputDir)
  chain=R.TChain("DVAnalyzer/tree_DV")
  files = []
  if args.inputDir[-1] != '/':
    args.inputDir += '/'
  print ('>> Creating list of files from : \n'+args.inputDir)
  command = '/bin/find '+args.inputDir+' -type f | grep root | grep -v failed | grep '+args.pattern
  str_files = subprocess.check_output(command,shell=True).splitlines()
  files.extend(['file:'+ifile for ifile in str_files])
  for file in files:
    print ">>Adding "+file
    chain.Add(file)

elif (args.inputeosDir != ""):
  print (args.inputeosDir)
  chain=R.TChain("DVAnalyzer/tree_DV")
  files = []
  if args.inputeosDir[-1] != '/':
    args.inputeosDir += '/'
  command = 'xrdfs root://cmseos.fnal.gov ls '+args.inputeosDir
  paths = subprocess.check_output(command,shell=True).splitlines()
  for path in paths:
    print ('>> Creating list of files from : \n'+path)
    command = 'xrdfs root://cmseos.fnal.gov ls -u '+path+"| grep '\.root' | grep "+args.pattern
    str_files = subprocess.check_output(command,shell=True).splitlines()
    files.extend(str_files)
  for file in files:
    print ">>Adding "+file
    chain.Add(file)

else:
  print('!!! no inputDir!!!')

histos = {}

histos["vtx_sigma_dBV"]=R.TH1F("vtx_sigma_dBV","vtx_sigma_dBV",100,0,0.05)
histos["vtx_dBV"]=R.TH1F("vtx_dBV","vtx_dBV",80,0,0.4)
histos["vtx_tkSize"]=R.TH1F("vtx_tkSize","vtx_tkSize",40,0,40)
histos["vtx_xy"]=R.TH2F("vtx_xy","vtx_xy",80,-4,4,80,-4,4)
histos["nvtx_per_event"]=R.TH1F("nvtx_per_event","nvtx_per_event",8,0,8);

#loop over all events
for ievt,evt in enumerate(chain):
  if(ievt%100000==0): print ('analyzing event {0}'.format(ievt))
  nVtx = 0
  for iv in range(0,len(evt.vtx_dBV)):
    sigma_dBV = evt.vtx_sigma_dBV[iv]
    dBV = evt.vtx_dBV[iv]
    tkSize = evt.vtx_track_size[iv]
    x = evt.vtx_x[iv]
    y = evt.vtx_y[iv]
    if(tkSize<3):
      continue
      #print("tkSize=0: {}".format(ievt))
      #print(evt.evt)
    if(x*x+y*y<2.09 and dBV>0.01 and sigma_dBV<0.0025):
      histos["vtx_tkSize"].Fill(tkSize)
    if(tkSize>=5 and dBV>0.01 and sigma_dBV<0.0025):
      histos["vtx_xy"].Fill(x, y)
    if(tkSize>=5 and x*x+y*y<2.09 and sigma_dBV<0.0025):
      histos["vtx_dBV"].Fill(dBV)
    if(tkSize>=5 and x*x+y*y<2.09 and dBV>0.01):
      histos["vtx_sigma_dBV"].Fill(sigma_dBV)
    if(tkSize>=5 and x*x+y*y<2.09 and dBV>0.01 and sigma_dBV<0.0025):
      nVtx += 1
  histos["nvtx_per_event"].Fill(nVtx)
    #histos["vtx_sigma_dBV"].Fill(evt.vtx_sigma_dBV[iv])
    #histos["vtx_dBV"].Fill(evt.vtx_dBV[iv])
    #histos["vtx_tkSize"].Fill(evt.vtx_track_size[iv])
    #histos["vtx_xy"].Fill(evt.vtx_x[iv], evt.vtx_y[iv])



fOut=R.TFile(args.output,"RECREATE")
for hn, histo in histos.iteritems():
  histo.Write()
fOut.Close()
print "Saved histos in "+args.output
