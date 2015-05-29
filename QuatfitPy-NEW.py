#! /usr/bin/env python2.7

# The MIT License (MIT)
#
# Copyright (c) 2015 Thomas J. L. Mustard, Paul Ha-Yeon Cheong
#
#   PHYC Group
#   Oregon State University
#   College of Science
#   Department of Chemistry
#   153 Gilbert Hall
#   Corvallis OR, 97331
#   E-mail:  mustardt@onid.orst.edu
#   Ph.  (541)-737-2081
#   http://phyc.chem.oregonstate.edu/
#
#   Permission is hereby granted, free of charge, to any person obtaining a copy
#   of this software and associated documentation files (the "Software"), to deal
#   in the Software without restriction, including without limitation the rights
#   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#   copies of the Software, and to permit persons to whom the Software is
#   furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in all
#   copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#   SOFTWARE.

import sys
import getopt
import shutil
import os
import math
import quatfit

### --- Arguments --- ###
program = "QuatfitPy-CLASSICAL.py"
reffile = ''
fitfile = ''
ofile = ''
pairsfile = ''
statfile = ''
debug = 0

# Read command line args
try:
  myopts, args = getopt.getopt(sys.argv[1:], "r:f:p:o:s:hd")
except getopt.GetoptError:
  print program + " -r <ref> -f <fit> -p <pairs> -o <out> -s <stat> -h -d"
  sys.exit(2)
###############################
# o == option
# a == argument passed to the o
###############################
for o, a in myopts:
  if o == '-r':
    reffile = a
  elif o == '-f':
    fitfile = a
  elif o == '-o':
    ofile = a
  elif o == '-p':
    pairsfile = a
  elif o == '-s':
    statfile = a
  elif o == '-d':
    debug += 1
  elif o == '-h':
    print program + " -r <ref> -f <fit> -p <pairs> -o <out> -s <stat> -h -d"
    sys.exit(0)
  else:
    print("Usage: %s -r <ref> -f <fit> -p <pairs> -o <out> -s <stat> -h -d" % sys.argv[0])
    sys.exit(0)

reffilelol = quatfit.parseXYZ(reffile)
fitfilelol = quatfit.parseXYZ(fitfile)

pairs, weights = quatfit.parsePairs(pairsfile)

fitfilelol, rms = quatfit.quatfitGetMolecule(reffilelol, fitfilelol, pairs, weights)

print "Weighted root mean square = " + str(rms)

quatfit.outputXYZ(ofile, fitfilelol, 1)

######################################################################
### END OF SCRIPT
######################################################################