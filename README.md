# Long exercise on ttbar at CMS DAS in Bari

This repository contains code for the long exercise on ttbar at the CMS data analysis school in Bari in 2015. The code is organised into three independent modules described in further details below.

To get started, initialise a read-only copy of the repository, checkout the recommended tag, and set up the environment:
```
git clone https://github.com/andrey-popov/CMSDASTTbar.git
cd CMSDASTTbar/
git checkout 1.2
source env.sh
```
The environment should be set up each time you log in.


## Reader

The module provides a set of C++ classes to read source ROOT files that contain basic properties of events, such as momenta of leptons and jets. The user is provided with an access to these properties. The code allows to evaluate systematical variations in jet energy corrections (JEC) and b-tagging.

The classes are utilised in an example program that produces a set of histograms of the MtW observable and stores them in a new ROOT file. The program can be compiled and executed with the following commands:
```
cd Reader/
make
./produceExampleHist
```
The source trees are pretty large, and the execution takes several minutes.


## Plotter

The module provides a C++ class to build a stacked distribution of simulated processes and compare it against data. It starts from a ROOT file that containst a set of histograms, with a structure similar to the file produced by the `Reader/produceExampleHist` program. An example program to plot a distribution can be compiled and executed with the following commands:
```
cd Plotter/
ln -s ../Reader/MtW.root  # The output file produced by the produceExampleHist program
make
./produceExamplePlot
```


## Fit

Provides an example script to fit data and extract experimental cross section of the ttbar production. Histograms created by the Reader are used as the input.
