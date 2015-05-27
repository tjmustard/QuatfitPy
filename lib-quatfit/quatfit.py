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

##############################################################################
# The program to superimpose atoms of two molecules by quaternion method
#
# David J. Heisterberg
# The Ohio Supercomputer Center
# 1224 Kinnear Rd.
# Columbus, OH  43212-1163
# (614)292-6036
# djh@ccl.net    djh@ohstpy.bitnet    ohstpy::djh
#
# Translated to C from fitest.f program and interfaced with Xmol program
# by Jan Labanowski,  jkl@ccl.net   jkl@ohstpy.bitnet   ohstpy::jkl
#
# Translated to python from quatfit.c
# by Thomas J. L. Mustard, mustardt@onid.orst.edu
#
# Copyright: Ohio Supercomputer Center, David J. Heisterberg, 1990.
# The program can be copied and distributed freely, provided that
# this copyright in not removed. You may acknowledge the use of the
# program in published material as:
# David J. Heisterberg, 1990, unpublished results.
#
#
# #include <iostream>
# #include <fstream>
# #include <cstdlib>
# #include <vector>
# #include <string>
# #include <cassert>
# #include <cmath>
# #include <algorithm>
# #include <sstream>
# #using namespace std
#
# #define MAXPOINTS     400
# #define MAXLINELEN    250
# #define OPTIONS       "r:f:p:o:s:"
#
#  options
#  -r refmol    reference molecule xmol file. If this option is not given, the
#               information is read from standard input.
#
#  -f fitmol    fitted molecule input xmol file If this option is not given, the
#               information is read from standard input.
#
#  -p pairs     file with the list of fitted atom pairs and weights. If this
#               option is not specified, pairs and weights are taken from
#               stdin. If file name "none" is used (i.e. -p none), atoms of the
#               fitted molecule are fitted to atoms of the reference
#               molecule with the same number, and all weights are assumed 1.0.
#               If molecules do not have the same number of atoms, the
#               smaller number of atoms is fitted.
#
#  -o outmol    output file for transformed fitted molecule. If this option
#               is not given, the information is written to standard output.
#
#  -s statfile  file with fit statistics (distances between fitted atoms, etc).
#               If this option is not given, the information is written to
#               standard output.
#
#  If any files are read from stdin, the order is: refmol, fitmol, pairs.
#  If any files are written to stdout, the order is: outmol, statfile.
#
#The file formats are:
# The refmol, fitmol and outmol files are in the XYZ format used with xmol
# program from Minnesota Supercomputer Institute. The format is:
#   1st line: number of atoms
#   2nd line: title
#   3rd and next lines have the format depending on the kind of information:
#       T  X  Y  Z                (total of 4 columns)
#       T  X  Y  Z  C             (total of 5 columns)
#       T  X  Y  Z  Mx My Mz      (total of 7 columns)
#       T  X  Y  Z  C Mx My Mz    (total of 8 columns)
#     where T is atom type (usually elemnt symbol), X, Y, Z are cartesian
#     coordinates of an atom in Angstroms, C is atomic charge, and Mx, My, Mz
#     are normal modes.
#
# The pairs file format is:
#   1st line: number of pairs
#   2nd and next lines:
#       Ar   Af    W
#     where Ar is the atom number of the reference molecule, Af is the atom
#     number of fitted molecule, w is the statistical weight. Weights W are
#     related to the square of expected deviation "sigma" between the reference
#     and fitted molecule atoms and allow to make fit of some atom pairs more
#     tight. W is proportional to 1/sigma^2. The larger the weight, the more
#     tight will be the resulting fit for the given pair.
#
# The statfile lists results of the fit with explanation.
#
# there was a typo in the formula for RMS. In the part calculating the
# RMS, there was a line:
#    wnorm += s
# while it should be:
#    wnorm += s*s
# Martin Lema mlema[-at-]unq.edu.ar was kind to find it. THANKS!!! jkl.
# It did not affect fitting too much. But the RMS value was not right when
# weights were different than 1.
#
####################################################################################

import os
from sys import *
import math
import sys
from decimal import *


###==================================================================================================================###
### --- Molecule class --- ###
class Molecule(object):
  def __init__(self):
    #Hold the string for the element
    self.e = ''
    #Hold the int for the element
    self.en = int(0)
    #Hold the float for the x coordinate
    self.x = float(0)
    #Hold the float for the y coordinate
    self.y = float(0)
    #Hold the float for the z coordinate
    self.z = float(0)
    #Hold an int/float for the charge of this atom
    self.charge = float(0)
    #Hold the float for the x coordinate
    self.mx = float(0)
    #Hold the float for the y coordinate
    self.my = float(0)
    #Hold the float for the z coordinate
    self.mz = float(0)



### --- Parse the input XYZ file and convert it into a Molecule Class list of lists --- ###
def parseXYZ(ifile):
  ### --- Open parent file, read it in and close it --- ###
  f = open(ifile, "r")
  ifileList = f.readlines()
  f.close()
  #Find out the final length of the ifilelol
  ifileLength = len(ifileList)
  #Fill the ifilelol with 0's so that there is a place to put the atomic data
  ifilelol = [0] * (ifileLength - 2)
  #### --- Input/parse the input file into a list of lists called ifilelol --- ###
  for i in range(2,ifileLength):
    line = ifileList[i].rstrip().split()
    ifilelol[i-2] = Molecule()
    ifilelol[i-2].e = line[0]
    ifilelol[i-2].x = float(line[1])
    ifilelol[i-2].y = float(line[2])
    ifilelol[i-2].z = float(line[3])
    if len(line) >= 5:
      ifilelol[i-2].charge = float(line[4])
    if len(line) == 8:
      ifilelol[i-2].mx = float(line[5])
      ifilelol[i-2].my = float(line[6])
      ifilelol[i-2].mz = float(line[7])
  #If the XYZ file is longer or shorter than the list of data should be, warn the user.
  if int(ifileList[0]) != len(ifilelol):
    print "Your file is the wrong length!"
    print "Make sure you have the right number of atoms."
    exit(0)
  return ifilelol

### --- Output the structural data in an XYZ format
def outputXYZ(ofile, ofilelol, printlevel):
  #Open the file for writing
  f = open(ofile, "w+")
  f.write(str(len(ofilelol)) + "\n")
  f.write(ofile + "\n")
  #iterate through the whole ifilelol and write each line in XYZ format to the file
  for i in range(len(ofilelol)):
    #if the instance is in Atom form print the 'element X-Coord Y-Coord Z-Coord'
    if isinstance(ofilelol[i], Molecule):
      if printlevel == 1:
        line = ofilelol[i].e + "  " + str("{:.6f}".format(ofilelol[i].x)) + "  " + str("{:.6f}".format(ofilelol[i].y)) + "  " + str("{:.6f}".format(ofilelol[i].z))
      elif printlevel == 2:
        line = ofilelol[i].e + "  " + str("{:.6f}".format(ofilelol[i].x)) + "  " + str("{:.6f}".format(ofilelol[i].y)) + "  " + str("{:.6f}".format(ofilelol[i].z)) + "  " + str("{:.6f}".format(ofilelol[i].charge))
      elif printlevel == 3:
        line = ofilelol[i].e + "  " + str("{:.6f}".format(ofilelol[i].x)) + "  " + str("{:.6f}".format(ofilelol[i].y)) + "  " + str("{:.6f}".format(ofilelol[i].z)) + "  " + str("{:.6f}".format(ofilelol[i].charge)) + "  " + str("{:.6f}".format(ofilelol[i].mx)) + "  " + str("{:.6f}".format(ofilelol[i].my)) + "  " + str("{:.6f}".format(ofilelol[i].mz))
    f.write(line + "\n")
  f.close
  return

### --- Parse the pairs file and return both the pairs and the weights --- ###
def parsePairs(pairsfile):
  ### --- Open parent file, read it in and close it --- ###
  f = open(pairsfile, "r")
  pairsfilelist = f.readlines()
  f.close()
  #Find out the final length of the ifilelol
  pairsfilelength = len(pairsfile)
  if int(pairsfilelist[0]) != len(pairsfilelist) - 1:
    print "The pairs file is the incorrect length.\nIt should be " \
          + str(int(pairsfilelist[0])) + " pairs long and it is " + str(len(pairsfilelist) - 1) + " pairs long."
  else:
    pairs = []
    weights=[]
    for i in range(len(pairsfilelist)):
      if i == 0:
        pairnum = int(pairsfilelist[i])
      elif i >= 1:
        line = pairsfilelist[i].strip().split()
        pairs.append([int(line[0]),int(line[1])])
        weights.append(int(line[2]))
  return pairs, weights

### --- Create the reference XYZ files for future alignment --- ###
def createRefGeom(reflol, fitlol, pairs):
  #Create ref_xyz and fit_xyz
  ref_xyz = [0] * len(pairs)
  fit_xyz = [0] * len(pairs)
  #Copy the pair atoms from the original molecule into the reference xyz file
  for i in range(len(pairs)):
    #Copy the ref pairs
    ref_xyz[i] = Molecule()
    ref_xyz[i].e = reflol[pairs[i][0]-1].e
    ref_xyz[i].x = reflol[pairs[i][0]-1].x
    ref_xyz[i].y = reflol[pairs[i][0]-1].y
    ref_xyz[i].z = reflol[pairs[i][0]-1].z
    ref_xyz[i].charge = reflol[pairs[i][0]-1].charge
    ref_xyz[i].mx = reflol[pairs[i][0]-1].mx
    ref_xyz[i].my = reflol[pairs[i][0]-1].my
    ref_xyz[i].mz = reflol[pairs[i][0]-1].mz
    #Copy the fit pairs
    fit_xyz[i] = Molecule()
    fit_xyz[i].e = fitlol[pairs[i][1]-1].e
    fit_xyz[i].x = fitlol[pairs[i][1]-1].x
    fit_xyz[i].y = fitlol[pairs[i][1]-1].y
    fit_xyz[i].z = fitlol[pairs[i][1]-1].z
    fit_xyz[i].charge = fitlol[pairs[i][1]-1].charge
    fit_xyz[i].mx = fitlol[pairs[i][1]-1].mx
    fit_xyz[i].my = fitlol[pairs[i][1]-1].my
    fit_xyz[i].mz = fitlol[pairs[i][1]-1].mz
  return ref_xyz, fit_xyz

### --- Print the molecule class object --- ###
def printMolecule(ofilelol):
  for i in range(len(ofilelol)):
    print ofilelol[i].e + "  " + str("{:.6f}".format(ofilelol[i].x)) + "  " + str("{:.6f}".format(ofilelol[i].y)) + "  " + str("{:.6f}".format(ofilelol[i].z)) + "  " + str("{:.6f}".format(ofilelol[i].charge)) + "  " + str("{:.6f}".format(ofilelol[i].mx)) + "  " + str("{:.6f}".format(ofilelol[i].my)) + "  " + str("{:.6f}".format(ofilelol[i].mz))
  return

### --- center the coordinates, or translate them to some xyz --- ###
def center(filelol, weights, centerswitch, centerxyz):
  '''====================================================================
  CENTER
   center or translate a molecule.
   atomnum (n) - number of atoms
   filelol (x) - on input  - original xyz coordinates of a molecule
       on output - moved xyz coordinates (see io for modes).

   weights (w) - if centerswitch=1, weights of atoms
       if centerswitch=2 or 3, unused

   centerswitch (io) - 1 weighted geometric center of the molecule will be at (0,0,0)
        2 molecule will be moved by a vector -center (i.e., components of a vector center
          will be subtracted from atom coordinates).
        3 molecule will be moved by a vector +center (i.e., components of a vector center
          will be added atom coordinates).

   centerxyz (o) - if centerswitch=1, output, center of original coordinates
       if centerswitch=2, input, vector center will be subtracted from atomic coordinates
       if centerswitch=3, input, vector center will be added to atomic coordinates

  ====================================================================='''
  wnorm = float(0.00)
  modif = float(0.00)
  #int i
  if (centerswitch == 2):
    modif = -1.0
  elif (centerswitch == 3):
    modif = 1.0
  else:
    modif = -1.0
    centerxyz[0] = float(0.0)
    centerxyz[1] = float(0.0)
    centerxyz[2] = float(0.0)
    wnorm = 0.0
    for i in range(len(filelol)):
      centerxyz[0] += filelol[i].x * math.sqrt(weights[i])
      centerxyz[1] += filelol[i].y * math.sqrt(weights[i])
      centerxyz[2] += filelol[i].z * math.sqrt(weights[i])
      wnorm += math.sqrt(weights[i])
    centerxyz[0] = centerxyz[0] / wnorm
    centerxyz[1] = centerxyz[1] / wnorm
    centerxyz[2] = centerxyz[2] / wnorm
  for i in range(len(filelol)):
    filelol[i].x = filelol[i].x + modif*centerxyz[0]
    filelol[i].y = filelol[i].y + modif*centerxyz[1]
    filelol[i].z = filelol[i].z + modif*centerxyz[2]
  return centerxyz, filelol


#void rotmol (int n, double x[4][MAXPOINTS], double y[4][MAXPOINTS], double u[4][4])
#{
def rotmol(filelol, rotmat):
  '''
  ROTMOL
  rotate a molecule
  n - number of atoms
  filelol (x) - input coordinates
  filelol (y) - rotated coordinates y = u * x
  rotmat (u) - left rotation matrix
  '''
  yx = float(0.0)
  yy = float(0.0)
  yz = float(0.0)

  for i in range(len(filelol)):
    yx = rotmat[0][0] * filelol[i].x + rotmat[1][0] * filelol[i].y + rotmat[2][0] * filelol[i].z
    yy = rotmat[0][1] * filelol[i].x + rotmat[1][1] * filelol[i].y + rotmat[2][1] * filelol[i].z
    yz = rotmat[0][2] * filelol[i].x + rotmat[1][2] * filelol[i].y + rotmat[2][2] * filelol[i].z
    filelol[i].x = yx #x
    filelol[i].y = yy #y
    filelol[i].z = yz #z
  return filelol


#void jacobi (double a[4][4], double d[4], double v[4][4], int nrot)
#{
def jacobi(matrix, maxsweeps):
  '''
  JACOBI
  Jacobi diagonalizer with sorted output. It is only good for 4x4 matrices.
  (was too lazy to do pointers...)
  matrix (a) - input: matrix to diagonalize
  eigenvect (v) - output: eigenvectors
  eigenval (d) - output: eigenvalues
  maxsweeps (nrot) - input: maximum number of sweeps
  '''

  eigenvect = [[float(0.0) for x in range(4)] for x in range(4)]
  eigenval = [float(0.0) for x in range(4)]
  onorm = float(0.0)
  dnorm = float(0.0)
  b = float(0.0)
  dma = float(0.0)
  q = float(0.0)
  t = float(0.0)
  c = float(0.0)
  s = float(0.0)
  atemp = float(0.0)
  vtemp = float(0.0)
  dtemp = float(0.0)

  for j in range(4):
    #for i in range(4):
    #  eigenvect[i][j] = 0.0
    eigenvect[j][j] = 1.0
    eigenval[j] = matrix[j][j]

  for m in range(maxsweeps):
    dnorm = 0.0
    onorm = 0.0
    for j in range(4):
      dnorm = dnorm + math.fabs(eigenval[j])
      for i in range(j):
        onorm = onorm + math.fabs(matrix[i][j])
    if (onorm/dnorm) <= 1.0e-12: break #goto Exit_now;
    for j in range(1,4):
      for i in range(j):
        b = matrix[i][j]
        if math.fabs(b) > 0.0:
          dma = eigenval[j] - eigenval[i]
          if (math.fabs(dma) + math.fabs(b)) <=  math.fabs(dma):
            t = b / dma
          else:
            q = 0.5 * dma / b
            t = 1.0/(math.fabs(q) + math.sqrt(1.0+q*q))
            if q < 0.0:
              t = -t
          c = 1.0/math.sqrt(t * t + 1.0)
          s = t * c
          matrix[i][j] = 0.0
          for k in range(i):
            atemp = c * matrix[k][i] - s * matrix[k][j]
            matrix[k][j] = s * matrix[k][i] + c * matrix[k][j]
            matrix[k][i] = atemp
          for k in range(i+1,j):
            atemp = c * matrix[i][k] - s * matrix[k][j]
            matrix[k][j] = s * matrix[i][k] + c * matrix[k][j]
            matrix[i][k] = atemp
          for k in range(j+1,4):
            atemp = c * matrix[i][k] - s * matrix[j][k]
            matrix[j][k] = s * matrix[i][k] + c * matrix[j][k]
            matrix[i][k] = atemp
          for k in range(4):
            vtemp = c * eigenvect[k][i] - s * eigenvect[k][j]
            eigenvect[k][j] = s * eigenvect[k][i] + c * eigenvect[k][j]
            eigenvect[k][i] = vtemp
          dtemp = c*c*eigenval[i] + s*s*eigenval[j] - 2.0*c*s*b
          eigenval[j] = s*s*eigenval[i] + c*c*eigenval[j] +  2.0*c*s*b
          eigenval[i] = dtemp

  ###Exit_now:

  maxsweeps = m

  for j in range(3):
    k = j
    dtemp = eigenval[k]
    for i in range(j+1,4):
      if eigenval[i] < dtemp:
        k = i
        dtemp = eigenval[k]

    if k > j:
      eigenval[k] = eigenval[j]
      eigenval[j] = dtemp
      for i in range(4):
        dtemp = eigenvect[i][k]
        eigenvect[i][k] = eigenvect[i][j]
        eigenvect[i][j] = dtemp

  return eigenvect, eigenval, maxsweeps



#void q2mat (double q[4], double u[4][4])
#{
def q2mat(quaternion):
  '''
  Q2MAT
  Generate a left rotation matrix from a normalized quaternion

  INPUT
    quaternion (q)      - normalized quaternion

  OUTPUT
    rotmat (u)      - the rotation matrix
  '''
  rotmat = [[float(0.0) for x in range(3)] for x in range(3)]
  rotmat[0][0] = quaternion[0]*quaternion[0] + quaternion[1]*quaternion[1] - quaternion[2]*quaternion[2] - quaternion[3]*quaternion[3]
  rotmat[1][0] = 2.0 * (quaternion[1] * quaternion[2] - quaternion[0] * quaternion[3])
  rotmat[2][0] = 2.0 * (quaternion[1] * quaternion[3] + quaternion[0] * quaternion[2])
  rotmat[0][1] = 2.0 * (quaternion[2] * quaternion[1] + quaternion[0] * quaternion[3])
  rotmat[1][1] = quaternion[0]*quaternion[0] - quaternion[1]*quaternion[1] + quaternion[2]*quaternion[2] - quaternion[3]*quaternion[3]
  rotmat[2][1] = 2.0 * (quaternion[2] * quaternion[3] - quaternion[0] * quaternion[1])
  rotmat[0][2] = 2.0 *(quaternion[3] * quaternion[1] - quaternion[0] * quaternion[2])
  rotmat[1][2] = 2.0 * (quaternion[3] * quaternion[2] + quaternion[0] * quaternion[1])
  rotmat[2][2] = quaternion[0]*quaternion[0] - quaternion[1]*quaternion[1] - quaternion[2]*quaternion[2] + quaternion[3]*quaternion[3]
  return rotmat


#void qtrfit (int n, double x[4][MAXPOINTS], double y[4][MAXPOINTS], double w[MAXPOINTS], double q[4], double u[4][4], int nr)
#{
def qtrfit(fit_xyz, ref_xyz, weights, maxsweeps):
  '''
   QTRFIT
   Find the quaternion, q,[and left rotation matrix, u] that minimizes

     |qTXq - Y| ^ 2  [|uX - Y| ^ 2]

   This is equivalent to maximizing Re (qTXTqY).

   This is equivalent to finding the largest eigenvalue and corresponding
   eigenvector of the matrix

   [A2   AUx  AUy  AUz ]
   [AUx  Ux2  UxUy UzUx]
   [AUy  UxUy Uy2  UyUz]
   [AUz  UzUx UyUz Uz2 ]

   where

     A2   = Xx Yx + Xy Yy + Xz Yz
     Ux2  = Xx Yx - Xy Yy - Xz Yz
     Uy2  = Xy Yy - Xz Yz - Xx Yx
     Uz2  = Xz Yz - Xx Yx - Xy Yy
     AUx  = Xz Yy - Xy Yz
     AUy  = Xx Yz - Xz Yx
     AUz  = Xy Yx - Xx Yy
     UxUy = Xx Yy + Xy Yx
     UyUz = Xy Yz + Xz Yy
     UzUx = Xz Yx + Xx Yz

   The left rotation matrix, u, is obtained from q by

     u = qT1q

   INPUT
     n      - number of points
     fit_xyz (x)      - fitted molecule coordinates
     ref_xyz (y)      - reference molecule coordinates
     weights (w)      - weights

   OUTPUT
     quaternion (q)      - the best-fit quaternion
     rotmat (u)      - the best-fit left rotation matrix
     maxsweeps (nr)     - max number of jacobi sweeps

  '''


  #Create variables/lists/matrixes
  matrix = [[float(0.0) for x in range(4)] for x in range(4)] #double c[4][4]
  quaternion = [float(0.0) for x in range(4)]

  # generate the upper triangle of the quadratic form matrix
  xxyx = float(0.0)
  xxyy = float(0.0)
  xxyz = float(0.0)
  xyyx = float(0.0)
  xyyy = float(0.0)
  xyyz = float(0.0)
  xzyx = float(0.0)
  xzyy = float(0.0)
  xzyz = float(0.0)

  for i in range(len(fit_xyz)):
    xxyx = xxyx + fit_xyz[i].x * ref_xyz[i].x * weights[i]
    xxyy = xxyy + fit_xyz[i].x * ref_xyz[i].y * weights[i]
    xxyz = xxyz + fit_xyz[i].x * ref_xyz[i].z * weights[i]
    xyyx = xyyx + fit_xyz[i].y * ref_xyz[i].x * weights[i]
    xyyy = xyyy + fit_xyz[i].y * ref_xyz[i].y * weights[i]
    xyyz = xyyz + fit_xyz[i].y * ref_xyz[i].z * weights[i]
    xzyx = xzyx + fit_xyz[i].z * ref_xyz[i].x * weights[i]
    xzyy = xzyy + fit_xyz[i].z * ref_xyz[i].y * weights[i]
    xzyz = xzyz + fit_xyz[i].z * ref_xyz[i].z * weights[i]

  matrix[0][0] = xxyx + xyyy + xzyz
  matrix[0][1] = xzyy - xyyz
  matrix[1][1] = xxyx - xyyy - xzyz
  matrix[0][2] = xxyz - xzyx
  matrix[1][2] = xxyy + xyyx
  matrix[2][2] = xyyy - xzyz - xxyx
  matrix[0][3] = xyyx - xxyy
  matrix[1][3] = xzyx + xxyz
  matrix[2][3] = xyyz + xzyy
  matrix[3][3] = xzyz - xxyx - xyyy

  # diagonalize c
  eigenvect, eigenval, maxsweeps = jacobi(matrix, maxsweeps)

  # extract the desired quaternion
  quaternion[0] = eigenvect[0][3]
  quaternion[1] = eigenvect[1][3]
  quaternion[2] = eigenvect[2][3]
  quaternion[3] = eigenvect[3][3]

  # generate the rotation matrix
  rotmat = q2mat(quaternion)

  return quaternion, rotmat, maxsweeps


'''=======================================================
 FITEST
 rigid fit test driver
 reads in data, fits, and writes out
'''

#int main(int argc, char *argv[])
#{
# int n_fields_r;                 no of fields in xmol file for ref. molec
# int n_fields_f;                 no of fields in xmol file for fit. molec
# int nat_r;                      number of all atoms in reference molecule
# int nat_f;                      number of all atoms in fitted molecule
# char title_r[MAXLINELEN];       title line in reference mol. xmol file
# char title_f[MAXLINELEN];       title line in fitted mol. xmol file
# char symb_r[MAXPOINTS][10];     atom type symbols for reference molecule
# char symb_f[MAXPOINTS][10];     atom type symbols for fitted molecule
# char line[MAXLINELEN];          scratch space for line
# double xyz_r[4][MAXPOINTS];     coordinates for reference molecule
# double xyz_f[4][MAXPOINTS];     coordinates for fitted molecule
# double charge_r[MAXPOINTS];     charges for reference molecule
# double charge_f[MAXPOINTS];     charges for fitted molecule
# double modes_r[4][MAXPOINTS];   normal modes for reference
# double modes_f[4][MAXPOINTS];   normal modes for fitted
# int npairs;                     no of fitted atom pairs
# int atoms_r[MAXPOINTS];         atoms of ref. molecule to be superimposed
# int atoms_f[MAXPOINTS];         atoms of fit. molecule to be superimposed
# double ref_xyz[4][MAXPOINTS];   ref. molecule atom coordinates to fit
# double fit_xyz[4][MAXPOINTS];   fit. molecule atom coordinates to fit
# double weight[MAXPOINTS];       fitted atom pair weights
# double ref_center[4];           center of ref. molecule fitted atoms
# double fit_center[4];           center of ref. molecule fitted atoms
# double quatnorm[4];                    quaternion
# double leftrotmatrix[4][4];                 left rotation matrix for coordinates
# int i, j, ch;                   aux variables
# double s, d, wd, rms, wnorm, dotm;    aux variables
# FILE *refmol;                   file variable for reference mol xmol file
# FILE *fitmol;                   file variable for reference mol xmol file
# FILE *pairinp;                  file with fitted atom pairs and weights
# FILE *outfile;                  rotated and translated fit. mol. xmol file
# FILE *statfile;                 file with fit goodness values
# int max_sweeps;                 max number of iterations in jacobi
# int opt;                        option letter
# int read_pairs;                 1 - read pairs, 0 - make pairs, weights
# static char usage[]=
#  "Usage : quatfit [-r ref] [-f fit] [-p pairs] [-o out] [-s stat]\n"
# extern char *optarg;            option argument from getopt
#
# # set defaults
# max_sweeps = 30
# read_pairs = 1
# refmol = stdin
# fitmol = stdin
# pairinp = stdin
# outfile = stdout
# statfile = stdout
#
# while ((opt = getopt(argc, argv, OPTIONS)) != EOF) {
#   switch (opt) {
#     case 'r':
#       if((refmol = fopen(optarg, "r")) == NULL)  {
#         fprintf(stderr,"Error: Could not find ref-mol-file: %s\n", optarg)
#         return(1)
#         }
#       break
#
#     case 'f':
#       if((fitmol = fopen(optarg, "r")) == NULL)  {
#         fprintf(stderr,"Error: Could not find fit-mol-file: %s\n", optarg)
#         return(2)
#         }
#       break
#
#     case 'p':
#       if(strcmp(optarg,"none") == 0) {  if argument is "none"
#         read_pairs = 0
#         break
#         }
#       if((pairinp = fopen(optarg, "r")) == NULL)  {
#         fprintf(stderr,
#                 "Error: Could not find pairs-wieight-file: %s\n", optarg)
#         return(3)
#         }
#       break
#
#     case 's':
#       if((statfile = fopen(optarg, "w")) == NULL) {
#         fprintf(stderr, "Error: Could not open out-file %s\n", optarg)
#         return(4)
#         }
#       break
#
#     case 'o':
#       if((outfile = fopen(optarg, "w")) == NULL)  {
#         fprintf(stderr,
#                 "Error: Could not open fitted-mol-out-file: %s\n", optarg)
#         return(5)
#         }
#       break
#
#     case '?':
#       fprintf(stderr,"Error: %s\n", usage)
#       return(6)
#     }   end switch
#   }  end while
#
#
# # Now read in the ref molecule
# while (fgets(line, MAXLINELEN, refmol) != NULL) {  scan until n_at found
#   i = 0
#   while ((ch = line[i++]) != '\0') {
#     if(!isspace(ch)) {
#       i = -1
#       break
#       }
#     }
#   if(i == -1) {
#     break
#     }
#   }
#
# if(i != -1) {  if white space only to the end
#  fprintf(stderr, "Error: Error in line 1 of refmol file\n")
#  return(7)
#  }
#
# if(sscanf(line, "%d\n", &nat_r) != 1) {   n  atoms
#   fprintf(stderr, "Error: Error in line 1 of refmol file.\n")
#   return(8)
#   }
#
# if(nat_r >= MAXPOINTS) {
#   fprintf(stderr,
#   "Error: Molecule too big. Recompile program with larger MAXPOINTS\n")
#   return(8)
#   }
#
# if(fgets(title_r,  MAXLINELEN, refmol) == NULL) {   title line
#   fprintf(stderr, "Error: Error in line 2 of refmol file.\n")
#   return(9)
#   }
#
# # read coordinates of ref molecule
# n_fields_r = 0
# for (i = 1; i <=  nat_r; i++) {
#   if(scan_line(refmol, &n_fields_r, &symb_r[i][0], &xyz_r[1][i], &xyz_r[2][i],
#                &xyz_r[3][i], &charge_r[i], &modes_r[1][i], &modes_r[2][i],
#                &modes_r[3][i]) != 0) {
#     fprintf(stderr,"Error: Error in line %d of refmol file.\n", i+2)
#     return(10)
#     }
#   }
#
# # Now read in the fitted molecule
# while (fgets(line, MAXLINELEN, fitmol) != NULL) {  scan until n_at found
#   i = 0
#   while ((ch = line[i++]) != '\0') {
#     if(!isspace(ch)) {
#       i = -1
#       break
#       }
#     }
#   if(i == -1) {
#     break
#     }
#   }
#
# if(i != -1) {  if white space only to the end
#  fprintf(stderr, "Error: Error in line 1 of fitmol file\n")
#  return(11)
#  }
#
# if(sscanf(line, "%d\n", &nat_f) != 1) {
#   fprintf(stderr, "Error: Error in line 1 of fitmol file.\n")
#   return(12)
#   }
#
# if(nat_f >= MAXPOINTS) {
#   fprintf(stderr,
#   "Error: Molecule too big. Recompile program with larger MAXPOINTS\n")
#   return(13)
#   }
#
# if(fgets(title_f,  MAXLINELEN, fitmol) == NULL) {
#   fprintf(stderr, "Error: Error in line 2 of fitmol file\n")
#   return(13)
#   }
#
# # read coordinates of fitted molecule
# n_fields_f = 0
# for (i = 1; i <=  nat_f; i++) {
#   if(scan_line(fitmol, &n_fields_f, &symb_f[i][0], &xyz_f[1][i], &xyz_f[2][i],
#                &xyz_f[3][i], &charge_f[i], &modes_f[1][i], &modes_f[2][i],
#                &modes_f[3][i]) != 0) {
#     fprintf(stderr,"Error: Error in line %d of fitmol file.\n", i+2)
#     return(14)
#     }
#   }
#
# # Read or set pairs and weights
# if(read_pairs == 1) {
#   while (fgets(line, MAXLINELEN, pairinp) != NULL) {  skip white space
#     i = 0
#     while ((ch = line[i++]) != '\0') {
#       if(!isspace(ch)) {
#         i = -1
#         break
#         }
#       }
#     if(i == -1) {
#       break
#       }
#     }
#
#   if(i != -1) {  if white space only to the end
#     fprintf(stderr, "Error: Error in line 1 of pairs and weights file\n")
#     return(15)
#     }
#
#   if(sscanf(line, "%d", &npairs) != 1) {
#     fprintf(stderr,
#             "Error: Error in line 1 of file with pairs and weights\n")
#     return(16)
#     }
#
#   if(npairs < 2) {
#     fprintf(stderr,
#             "Error: Cannot fit a single atom. Need at least 2\n")
#     return(16)
#     }
#
#   for (i = 1; i <= npairs; i++) {
#     if(fscanf(pairinp,"%d%d%le",&atoms_r[i], &atoms_f[i], &weight[i]) != 3) {
#       fprintf(stderr,
#               "Error: Error in line %d of pairs and weights file.\n", i+2)
#       return(17)
#       }
#     if((atoms_r[i] < 1) || (atoms_r[i] > nat_r) || (atoms_f[i] < 1) ||
#        (atoms_f[i] > nat_f) ) {
#       fprintf(stderr,
#       "Error: Error in line %d of pairs and weights (number out of range)\n",
#       i+2)
#       return(18)
#       }
#     }
#   }
# else {  # initialize pairs to consectutive numbers
#   npairs = nat_r
#   if(nat_r > nat_f) {
#     npairs = nat_f
#     }
#
#   for(i = 1; i <= npairs; i++) {
#     atoms_r[i] = i
#     atoms_f[i] = i
#     weight[i] = 1.0
#     }
#   }
#
# # extract fitted atoms to tables
# for (i = 1; i <= npairs; i++) {
#   for (j = 1; j <= 3; j++) {
#     ref_xyz[j][i] = xyz_r[j][atoms_r[i]]
#     fit_xyz[j][i] = xyz_f[j][atoms_f[i]]
#     }
#   }
#
# # ===  Atom coordinates are fit in both modes ===
# # center ref molecule fitted atoms around (0,0,0)
# center (npairs, ref_xyz, weight, 1, ref_center)
#
# # center fitted molecule fitted atoms around (0,0,0)
# center (npairs, fit_xyz, weight, 1, fit_center)
#
# # fit specified atoms of fit_molecule to those of ref_molecule
# qtrfit(npairs, fit_xyz, ref_xyz, weight, quatnorm, u, max_sweeps)
#
# ''' subtract coordinates of the center of fitted atoms of the fitted molecule
#    from all atom coordinates of the fitted molecule (note that weight is
#    a dummy parameter) '''
# center(nat_f, xyz_f, weight, 2, fit_center)
#
# # rotate the fitted molecule by the rotation matrix u
# rotmol(nat_f, xyz_f, xyz_f, u)
# # same with set of fitted atoms of the fitted molecule
# rotmol(npairs, fit_xyz, fit_xyz, u)
#
# # if modes given in fitted molecule, rotate the modes too
# if(n_fields_f > 4) {
#   rotmol(nat_f, modes_f, modes_f, u)
#
#   # calculate dot product of reference and fitted molecule modes
#   if(n_fields_r > 4) {
#     dotm = 0.0
#     for (i = 1; i <= npairs; i++) {
#       for (j = 1; j <= 3; j++) {
#         dotm += modes_r[j][atoms_r[i]]*modes_f[j][atoms_f[i]]
#       }
#     }
#   }
# }
#
# ''' translate atoms of the fitted molecule to the center
#      of fitted atoms of the reference molecule '''
# center(nat_f, xyz_f, weight, 3, ref_center)
# # same with set of fitted atoms of the fitted molecule
# center(npairs, fit_xyz, weight, 3, ref_center)
# # translate fitted atoms of reference molecule to their orig. location
# center(npairs, ref_xyz, weight, 3, ref_center)
#
# # write modified XYZ file for fitted molecule
# fprintf(outfile,"%d\n%s", nat_f, title_f)
# for (i = 1; i <= nat_f; i++) {
#   write_line(outfile, n_fields_f, symb_f[i], xyz_f[1][i], xyz_f[2][i],
#              xyz_f[3][i], charge_f[i], modes_f[1][i], modes_f[2][i],
#              modes_f[3][i])
#   }
#
# fflush(outfile)
# ''' find distances between fitted and reference atoms and print them in
#    out file '''
#
# fprintf(statfile,
#          "\nDistances and weighted distances between fitted atoms\n")
#fprintf(statfile,"Ref.At. Fit.At.  Distance  Dist*math.sqrt(weight)  weight\n")
# rms = 0.0
# wnorm = 0.0
# for (i = 1; i <= npairs; i++) {
#   d = 0.0
#   for (j = 1; j <= 3; j++) {
#     s = ref_xyz[j][i] - fit_xyz[j][i]
#     d += s*s
#     }
#   rms += d
#   //cout << "D is " << d << " s is " << s << " wd is " << wd << " wnorm is " << wnorm <<  " numPairs " << npairs << endl
#   fprintf(statfile, "  %3d    %3d  %11.6f   %11.6f    %11.6f\n", atoms_r[i], atoms_f[i], d, wd, weight[i])
#   }
#
#
# rms = math.sqrt(rms/npairs);//wnorm)
# fprintf(statfile, "\n\nWeighted root mean square=%10.6f\n\n", rms)
# fprintf(statfile, "\n\nCenter of reference molecule fitted atoms\n")
# fprintf(statfile, "Xc = %11.6f Yc = %11.6f Zc = %11.6f\n",
#         ref_center[1], ref_center[2], ref_center[3])
#
# ''' {
#  double rmsd = 0.0
#  double x,y,z
#  for(int i = 1; i <= npairs; i++)
#  {
#    x = ref_xyz[1][i] - fit_xyz[1][i]
#    y = ref_xyz[2][i] - fit_xyz[2][i]
#    z = ref_xyz[3][i] - fit_xyz[3][i]
#
#    rmsd += x*x
#    rmsd += y*y
#    rmsd += z*z
#    //rmsd += (pow(ref_xyz[1][i] - fit_xyz[1][i], 2) + pow(ref_xyz[2][i] - fit_xyz[2][i], 2) + pow(ref_xyz[3][i] - fit_xyz[3][i], 2))
#  }
#  //rmsd /= static_cast<double>(npairs)
#  rms = math.sqrt(rmsd/npairs)
#  }'''
#
# fprintf(statfile, "\n\nWeighted root mean square=%10.6f\n\n", rms)
# fprintf(statfile, "\n\nCenter of reference molecule fitted atoms\n")
# fprintf(statfile, "Xc = %11.6f Yc = %11.6f Zc = %11.6f\n",
#         ref_center[1], ref_center[2], ref_center[3])
#
# fprintf(statfile, "\n\nCenter of fitted molecule fitted atoms\n")
# fprintf(statfile, "Xc = %11.6f Yc = %11.6f Zc = %11.6f\n",
#         fit_center[1], fit_center[2], fit_center[3])
#
# fprintf(statfile,"\n\nLeft rotation matrix\n")
# for (i = 1; i <= 3; i++) {
#   fprintf(statfile, " %11.6f  %11.6f  %11.6f\n",
#           u[1][i], u[2][i], u[3][i])
#   }
#
# if((n_fields_f > 4) && (n_fields_r > 4)) {
#   fprintf(statfile,
#   "\nDot product of normal modes on fitted atom pairs =%11.6f\n", dotm)
#   }
#
# fclose(outfile)
# fclose(statfile)
# fclose(pairinp)
# fclose(fitmol)
# fclose(refmol)
#
# return(0)
#
# }
#