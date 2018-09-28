#!/bin/bash

g++ -std=c++14 validation_histograms.cxx -L$(root-config --libdir) -lCore -lVMC -lHist -lGpad -lGraf -I$(root-config --incdir) -I/home/bvolkel/alice/software/dev/concurrentEngines/sw_test/MCStepLogger/include/ -L/home/bvolkel/alice/software/dev/concurrentEngines/sw_test/MCStepLogger/lib -lMCStepLogger -o validationHistograms
