#!/bin/bash
g++ main.cpp  -lLHAPDF -O2  `gsl-config --cflags` `gsl-config --libs` ` /opt/homebrew/Cellar/lhapdf/6.5.4/bin/lhapdf-config --cflags` ` /opt/homebrew/Cellar/lhapdf/6.5.4/bin/lhapdf-config --cflags --ldflags`
