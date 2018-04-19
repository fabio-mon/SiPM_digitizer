# SiPM_digitizer

Produce sampled output of SiPM given the photon arrival time and some other parameters

#instruction to compile on lxplus (ROOT 6 is required)
   source /afs/cern.ch/sw/lcg/external/gcc/4.9/x86_64-slc6-gcc49-opt/setup.sh

   source /afs/cern.ch/sw/lcg/app/releases/ROOT/6.06.06/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh

   c++ ConfigFile.cc utilities.cc -o digitizer digitizer.cpp `root-config --glibs --cflags`	


