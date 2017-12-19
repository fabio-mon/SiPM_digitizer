#!/usr/bin/env python
import os, re
import commands
import math, time
import sys
import argparse
import subprocess

parser = argparse.ArgumentParser(description='This script creates the script for the digitizing of simulated light')

parser.add_argument("-l", "--label",          required=True,     type=str,  help="job label")
parser.add_argument("-e", "--exe",            required=True,     type=str,  help="executable")
parser.add_argument("-o", "--outputFolder",   required=True,     type=str,  help="folder where to store output files")
#parser.add_argument("-f", "--outputFileName", required=True,     type=str,  help="base name of output files [outputFileName]_i.root")
parser.add_argument("-c", "--configFile",     required=True,     type=str,  help="config file to be run")
#parser.add_argument("-g", "--gpsFile",        required=True,     type=str,  help="gps.mac file to be run")
parser.add_argument("-i", "--inputFileName",   required=True,     type=str,  help="temp name of input root files [inputFileName]")
parser.add_argument("-n", "--nJobs",           required=True,     type=int,  help="number of jobs")
parser.add_argument("-t", "--thread",          required=True,     type=int,  help="number of thread per job")
parser.add_argument("-q", "--queue",          default="1nd",     type=str,  help="hercules queue to use")
parser.add_argument("-s", "--submit",                                       help="submit jobs", action='store_true')
parser.add_argument("-v", "--verbose",                                      help="increase output verbosity", action='store_true')


args = parser.parse_args()


print 
print 'START'
print 

currDir = os.getcwd()

print

try:
   subprocess.check_output(['mkdir','jobs'])
except subprocess.CalledProcessError as e:
   print e.output
try:
   subprocess.check_output(['mkdir','jobs/'+args.label])
except subprocess.CalledProcessError as e:
   print e.output
try:
   subprocess.check_output(['mkdir',args.outputFolder+"/"+args.label+"/"])
except subprocess.CalledProcessError as e:
   print e.output


##### loop for creating and sending jobs #####
for x in range(1, args.nJobs+1):

   ##### creates directory and file list for job #######
   jobDir = currDir+'/jobs/'+args.label+'/job_'+str(x)
   os.system('mkdir '+jobDir)
   os.chdir(jobDir)

   ##### copy executable to the jobDir ######
   os.system('cp '+args.exe+' '+jobDir+"/executable.exe")

   ##### creates Geant4 config file #######
   with open(args.configFile) as fi:
      contents = fi.read()
      filename = args.inputFileName.replace('_seed%', str("_seed"+str(x)))
      filename = filename.replace('_t%.', str(".") )
      replaced_contents = contents.replace('FILENAME', str(filename))
   with open(jobDir+"/config.cfg", "w") as fo:
      fo.write(replaced_contents)

   ##### creates job #######
   with open('job_'+str(x)+'.sh', 'w') as fout:
      fout.write("#!/bin/sh\n")
      fout.write("echo\n")
      fout.write("echo 'START---------------'\n")
      fout.write("echo 'current dir: ' ${PWD}\n")
      fout.write("source /afs/cern.ch/sw/lcg/external/gcc/4.9/x86_64-slc6-gcc49-opt/setup.sh\n")
      fout.write("source /afs/cern.ch/sw/lcg/app/releases/ROOT/6.06.06/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh\n")
      inputreplaced_contents = args.inputFileName.replace('_seed%', str("_seed"+str(x)))
      inputreplaced_contents = inputreplaced_contents.replace('_t%', str("_t*"))
      outputreplaced_contents = args.inputFileName.replace('_seed%', str("_seed"+str(x)))
      outputreplaced_contents = outputreplaced_contents.replace('_t%.', str(".") )
      fout.write(str("hadd " + outputreplaced_contents + " " + inputreplaced_contents + "\n"))
      fout.write("mkdir "+str(args.outputFolder+"/"+args.label)+"\n")
      fout.write("cd "+str(jobDir)+"\n")
      fout.write("echo 'current dir: ' ${PWD}\n")
      fout.write("./executable.exe config.cfg\n")#+args.outputFolder+"/"+args.label+"/"+args.outputFileName+"_"+str(x)+"\n")
      fout.write("mv *.root " +str(args.outputFolder+"/"+args.label)+"\n")
      fout.write("echo 'STOP---------------'\n")
      fout.write("echo\n")
      fout.write("echo\n")
   os.system("chmod 755 job_"+str(x)+".sh")
   os.system("chmod +x executable.exe")

   ###### sends bjobs ######
   if args.submit:
      os.system("bsub -q "+args.queue+" job_"+str(x)+".sh")
      print "job nr. " + str(x) + " submitted"
   
   os.chdir("../..")
   
print
print "your jobs:"
os.system("bjobs")
print
print 'END'
print



