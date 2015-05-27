#! /usr/bin/env python2.7

#   The MIT License (MIT)
#
#   Copyright (c) 2015 Thomas J. L. Mustard, Paul Ha-Yeon Cheong
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
refFile = ''
fitFile = ''
ofile = ''
pairsFile = ''
statFile = ''
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

#Create the xyz coord in Molecule class
reffilelol = quatfit.parseXYZ(reffile)
fitfilelol = quatfit.parseXYZ(fitfile)

#print reffile
#quatfit.printMolecule(reffilelol)
#print fitfile
#quatfit.printMolecule(fitfilelol)

#Parse the pairs file
pairs, weights = quatfit.parsePairs(pairsfile)
#print pairs
#print weights

#Create the reference coords based on the pairs
ref_xyz, fit_xyz = quatfit.createRefGeom(reffilelol, fitfilelol, pairs)

#Center the reference coords around 0,0,0
refcenter, ref_xyz = quatfit.center(ref_xyz, weights, 1, [float(0),float(0),float(0)])
fitcenter, fit_xyz = quatfit.center(fit_xyz, weights, 1, [float(0),float(0),float(0)])

#print "ref_xyz"
#print refcenter
#quatfit.printMolecule(ref_xyz)
#print "fit_xyz"
#print fitcenter
#quatfit.printMolecule(fit_xyz)

#fit the specified atom coords of the fit to reference
quaternion, rotmat, maxsweeps = quatfit.qtrfit(fit_xyz, ref_xyz, weights, 30)

#print quaternion
#print rotmat
#print maxsweeps

#subtract coordinates of the center of fitted atoms of the fitted molecule
#from all atom coordinates of the fitted molecule (note that weight is
#a dummy parameter)
fitcenter, fitfilelol = quatfit.center(fitfilelol, weights, 2, fitcenter)

# rotate the fitted molecule by the rotation matrix u
fitfilelol = quatfit.rotmol(fitfilelol, rotmat)

# same with set of fitted atoms of the fitted molecule
fit_xyz = quatfit.rotmol(fit_xyz, rotmat)

### Haven't yet fully translated this section
# if modes given in fitted molecule, rotate the modes too
#if n_fields_f > 4:
#  rotmol(nat_f, modes_f, modes_f, u)
#  # calculate dot product of reference and fitted molecule modes
#  if n_fields_r > 4:
#    dotm = 0.0
#    #for (i = 1; i <= npairs; i++) {
#    for i in range(1,len(pairs)):
#      #for (j = 1; j <= 3; j++) {
#      for j in range(1,3):
#        dotm += modes_r[j][atoms_r[i]]*modes_f[j][atoms_f[i]]

# translate atoms of the fitted molecule to the center
# of fitted atoms of the reference molecule
refcenter, fitfilelol = quatfit.center(fitfilelol, weights, 3, refcenter)

# same with set of fitted atoms of the fitted molecule
refcenter, fit_xyz = quatfit.center(fit_xyz, weights, 3, refcenter)

# translate fitted atoms of reference molecule to their orig. location
refcenter, ref_xyz = quatfit.center(ref_xyz, weights, 3, refcenter)


# write modified XYZ file for fitted molecule */
#quatfit.printMolecule(fitfilelol)
quatfit.outputXYZ("OUT-" + fitfile, fitfilelol, 1)

# find distances between fitted and reference atoms and print them in out file
print "Distances and weighted distances between fitted atoms"
print "Ref.At. Fit.At.  Distance  Dist*sqrt(weight)  weight"

rms = 0.0
wnorm = 0.0
for i in range(len(pairs)):
  d = 0.0
  for j in range(3):
    if j == 0:
      s = ref_xyz[i].x - fit_xyz[i].x
    elif j == 1:
      s = ref_xyz[i].y - fit_xyz[i].y
    elif j == 2:
      s = ref_xyz[i].z - fit_xyz[i].z
    d += s*s
  rms += d
  print "  " + str("{0: <4}".format(ref_xyz[i].e)) + "    " + str("{0: <5}".format(fit_xyz[i].e)) + "  "\
        + str("{:.6f}".format(d)) + "      " + str("{:.6f}".format(weights[i]*d)) + "      "+ str("{:.6f}".format(weights[i]))

rms = math.sqrt(rms/len(pairs))
print  "\nWeighted root mean square = " + str("{:.6f}".format(rms))

print  "\nCenter of reference molecule fitted atoms"
print  "Xc = " + str("{:.6f}".format(refcenter[0])) + " Yc = " + str("{:.6f}".format(refcenter[1])) + " Zc = " + str("{:.6f}".format(refcenter[2]))

print "\nCenter of fitted molecule fitted atoms"
print "Xc = " + str("{:.6f}".format(fitcenter[0])) + " Yc = " + str("{:.6f}".format(fitcenter[1])) + " Zc = " + str("{:.6f}".format(fitcenter[2]))

print "\nLeft rotation matrix"
for i in range(3):
  print " " + str("{:10.6f}".format(rotmat[i][0])) + "  " + str("{:10.6f}".format(rotmat[i][1])) + "  " + str("{:10.6f}".format(rotmat[i][2]))


#/* {
# double rmsd = 0.0;
# double x,y,z;
# for(int i = 1; i <= npairs; i++)
# {
#   x = ref_xyz[1][i] - fit_xyz[1][i];
#   y = ref_xyz[2][i] - fit_xyz[2][i];
#   z = ref_xyz[3][i] - fit_xyz[3][i];
#
#   rmsd += x*x;
#   rmsd += y*y;
#   rmsd += z*z;
#   //rmsd += (pow(ref_xyz[1][i] - fit_xyz[1][i], 2) + pow(ref_xyz[2][i] - fit_xyz[2][i], 2) + pow(ref_xyz[3][i] - fit_xyz[3][i], 2));
# }
# //rmsd /= static_cast<double>(npairs);
# rms = sqrt(rmsd/npairs);
# }*/
#if((n_fields_f > 4) && (n_fields_r > 4)) {
#  fprintf(statfile,
#  "\nDot product of normal modes on fitted atom pairs =%11.6f\n", dotm);
#  }

######################################################################
### END OF SCRIPT
######################################################################