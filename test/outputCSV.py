import ROOT as R
import math as M
import argparse
import subprocess
import csv


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
    print (">>Adding "+file)
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
    print (">>Adding "+file)
    chain.Add(file)

else:
  print('!!! no inputDir!!!')

var = ['vtx_track_size','vtx_dBV','vtx_sigma_dBV','vtx_x','vtx_y','vtx_z','evt'] #title of csv
#nr = 0
with open(args.output,'w') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(var)
    #loop over all events
    for ievt,evt in enumerate(chain):
      if(ievt%100000==0): print ('analyzing event {0}'.format(ievt))
      for iv in range(0,len(evt.vtx_dBV)):
        sigma_dBV = evt.vtx_sigma_dBV[iv]
        dBV = evt.vtx_dBV[iv]
        tkSize = evt.vtx_track_size[iv]
        x = evt.vtx_x[iv]
        y = evt.vtx_y[iv]
        z = evt.vtx_z[iv]
        n_evt = evt.evt
        writer.writerow([tkSize, dBV, sigma_dBV, x,y,z, n_evt])
        #nr += 1
        #if(nr>3770 and nr<3780):
        #  print([tkSize, dBV, sigma_dBV, x,y,z])

print ("Saved data in "+args.output)
