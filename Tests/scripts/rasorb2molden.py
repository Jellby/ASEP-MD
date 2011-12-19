#!/usr/bin/python

# Use with 2 arguments:
#  1 (read):       Input RasOrb file
#  2 (read/write): Input molden file

import sys
import os
import re
import math
import fileinput

#=============================
# Get input files

try:
  rasorb_input = sys.argv[1]
except IndexError:
  sys.exit("Missing input rasorb file")

try:
  molden_input = sys.argv[2]
except IndexError:
  sys.exit("Missing input molden file")

dir = os.path.dirname(sys.argv[0])

#=============================
# Read the orbitals from RasOrb

orbital = []
symmetry = []
energy = []
occupancy = []

file_rasorb = open(rasorb_input, "r")
for line in file_rasorb:
  if (re.search("#INFO", line, re.I)): break
file_rasorb.next()
file_rasorb.next()
nbas = int(file_rasorb.next())
norb = int(file_rasorb.next())

for line in file_rasorb:
  if (re.search("#ORB", line, re.I)): break
for oi in range(norb):
  orbital.append([])
  sym = re.search("\* ORBITAL\s+(\d+)\s+(\d+)", file_rasorb.next(), re.I).group(1)
  symmetry.append(int(sym))
  num = (int(nbas)-1)/4+1
  text = ""
  for i in range(num):
    text += file_rasorb.next().rstrip()
  for i in range(nbas):
    orbital[oi].append(float(text[18*i:18*(i+1)]))

for line in file_rasorb:
  if (re.search("#OCC", line, re.I)): break
file_rasorb.next()
num = (int(norb)-1)/4+1
text = ""
for i in range(num):
  text += file_rasorb.next().rstrip()
for i in range(norb):
  occupancy.append(float(text[18*i:18*(i+1)]))

for line in file_rasorb:
  if (re.search("#ONE", line, re.I)): break
try:
  file_rasorb.next()
  num = (int(norb)-1)/4+1
  text = ""
  for i in range(num):
    text += file_rasorb.next().rstrip()
  for i in range(norb):
    energy.append(float(text[18*i:18*(i+1)]))
except:
  for i in range(norb):
    energy.append(0.0)

file_rasorb.close()

#=============================
# Read the basis set from molden

basis = []
numd = 6
numf = 10
numg = 15

file_molden = open(molden_input, "r")
# Cartesian or spherical
for line in file_molden:
  if (re.search("\[5D\]", line, re.I)): numd = 5; numf = 7
  if (re.search("\[5D7F\]", line, re.I)): numd = 5; numf = 7
  if (re.search("\[5D10F\]", line, re.I)): numd = 5
  if (re.search("\[7F\]", line, re.I)): numf = 7
  if (re.search("\[9G\]", line, re.I)): numg = 9
file_molden.seek(0)
# Locate the beginning of the basis set
for line in file_molden:
  if (re.search("\[GTO\]", line, re.I)): break
# Get the basis set for each atom as a string like "sssppd"
for line in file_molden:
  if (re.search("\[", line)): break
  num = re.search("^\s*(\d+)\s*$", line)
  if num:
    num = int(num.group(1))-1
    basis.append("")
    while True:
      shell = re.search("^\s*([spdfg])\s*(\d+)\s*$", file_molden.next())
      if (not shell): break
      basis[num] += shell.group(1)
      for i in range(int(shell.group(2))):
        file_molden.next()
file_molden.close()

numbas = 0
num = 0
order = []
for shell in basis:
  funcs = re.search("(s*)(p*)(d*)(f*)(g*)", shell)
  # Count total basis functions
  numbas += len(funcs.group(1)) + 3*len(funcs.group(2)) + numd*len(funcs.group(3)) + numf*len(funcs.group(4)) + numg*len(funcs.group(5))

if (numbas == 0):
  sys.exit("No basis set found in molden file")
if (numbas != nbas):
  sys.exit("Number of basis functions does not match (%i in RasOrb, %i in molden)\n\
\'Contaminants\' are not supported" % (nbas,numbas) )

for shell in basis:
  funcs = re.search("(s*)(p*)(d*)(f*)(g*)", shell)
  numsame = len(funcs.group(1))
  # Reorder s shells (easy)
  for i in range(numsame):
    order.append(num+i)
  num += numsame
  # Reorder p shells
  numsame = len(funcs.group(2))
  for i in range(numsame):
    order.append(num+i+0*numsame)
    order.append(num+i+1*numsame)
    order.append(num+i+2*numsame)
  num += 3*numsame
  # Reorder d shells
  numsame = len(funcs.group(3))
  if (numd == 5):
    for i in range(numsame):
      order.append(num+i+2*numsame)
      order.append(num+i+3*numsame)
      order.append(num+i+1*numsame)
      order.append(num+i+4*numsame)
      order.append(num+i+0*numsame)
    num += numd*numsame
  else:
    for i in range(numsame):
      order.append(num+i+0*numsame)
      order.append(num+i+3*numsame)
      order.append(num+i+5*numsame)
      order.append(num+i+1*numsame)
      order.append(num+i+2*numsame)
      order.append(num+i+4*numsame)
      # Refer to normalized functions
      for j in range(norb):
        orbital[j][num+i+0*numsame] *= math.sqrt(3) # x^2 -> sqrt(1*3)
        orbital[j][num+i+3*numsame] *= math.sqrt(3)
        orbital[j][num+i+5*numsame] *= math.sqrt(3)
    num += numd*numsame
  # Reorder f shells
  numsame = len(funcs.group(4))
  if (numf == 7):
    for i in range(numsame):
      order.append(num+i+3*numsame)
      order.append(num+i+4*numsame)
      order.append(num+i+2*numsame)
      order.append(num+i+5*numsame)
      order.append(num+i+1*numsame)
      order.append(num+i+6*numsame)
      order.append(num+i+0*numsame)
    num += numf*numsame
  else:
    for i in range(numsame):
      order.append(num+i+0*numsame)
      order.append(num+i+6*numsame)
      order.append(num+i+9*numsame)
      order.append(num+i+3*numsame)
      order.append(num+i+1*numsame)
      order.append(num+i+2*numsame)
      order.append(num+i+5*numsame)
      order.append(num+i+8*numsame)
      order.append(num+i+7*numsame)
      order.append(num+i+4*numsame)
      # Refer to normalized functions
      for j in range(norb):
        orbital[j][num+i+0*numsame] *= math.sqrt(15) # x^3 -> sqrt(1*3*5)
        orbital[j][num+i+6*numsame] *= math.sqrt(15)
        orbital[j][num+i+9*numsame] *= math.sqrt(15)
        orbital[j][num+i+3*numsame] *= math.sqrt(3)  # x^2 -> sqrt(1*3)
        orbital[j][num+i+1*numsame] *= math.sqrt(3)
        orbital[j][num+i+2*numsame] *= math.sqrt(3)
        orbital[j][num+i+5*numsame] *= math.sqrt(3)
        orbital[j][num+i+8*numsame] *= math.sqrt(3)
        orbital[j][num+i+7*numsame] *= math.sqrt(3)
    num += numf*numsame
  # Reorder g shells
  numsame = len(funcs.group(5))
  if (numg == 9):
    for i in range(numsame):
      order.append(num+i+4*numsame)
      order.append(num+i+5*numsame)
      order.append(num+i+3*numsame)
      order.append(num+i+6*numsame)
      order.append(num+i+2*numsame)
      order.append(num+i+7*numsame)
      order.append(num+i+1*numsame)
      order.append(num+i+8*numsame)
      order.append(num+i+0*numsame)
    num += numg*numsame
  else:
    for i in range(numsame):
      order.append(num+i+ 0*numsame)
      order.append(num+i+10*numsame)
      order.append(num+i+14*numsame)
      order.append(num+i+ 1*numsame)
      order.append(num+i+ 2*numsame)
      order.append(num+i+ 6*numsame)
      order.append(num+i+11*numsame)
      order.append(num+i+ 9*numsame)
      order.append(num+i+13*numsame)
      order.append(num+i+ 3*numsame)
      order.append(num+i+ 5*numsame)
      order.append(num+i+12*numsame)
      order.append(num+i+ 4*numsame)
      order.append(num+i+ 7*numsame)
      order.append(num+i+ 8*numsame)
      # Refer to normalized functions
      for j in range(norb):
        orbital[j][num+i+ 0*numsame] *= math.sqrt(105) # x^4 -> sqrt(1*3*5*7)
        orbital[j][num+i+10*numsame] *= math.sqrt(105)
        orbital[j][num+i+14*numsame] *= math.sqrt(105)
        orbital[j][num+i+ 1*numsame] *= math.sqrt(15)  # x^3 -> sqrt(1*3*5)
        orbital[j][num+i+ 2*numsame] *= math.sqrt(15)
        orbital[j][num+i+ 6*numsame] *= math.sqrt(15)
        orbital[j][num+i+11*numsame] *= math.sqrt(15)
        orbital[j][num+i+ 9*numsame] *= math.sqrt(15)
        orbital[j][num+i+13*numsame] *= math.sqrt(15)
        orbital[j][num+i+ 3*numsame] *= 3              # x^2*y^2 -> sqrt(1*3)*sqrt(1*3)
        orbital[j][num+i+ 5*numsame] *= 3
        orbital[j][num+i+12*numsame] *= 3
        orbital[j][num+i+ 4*numsame] *= math.sqrt(3)   # x^2 -> sqrt(1*3)
        orbital[j][num+i+ 7*numsame] *= math.sqrt(3)
        orbital[j][num+i+ 8*numsame] *= math.sqrt(3)
    num += numg*numsame

#=============================
# Write the orbitals in molden format

file_molden = fileinput.input(molden_input, inplace=1)
for line in file_molden:
  if (re.search("\[MO\]", line, re.I)): break
  print line.rstrip()

print " [MO]"
for i in range(norb):
  print "Sym= %3ia" % (i+1)
  print "Ene= %10.4f" % energy[i]
  print "Spin= Alpha"
  print "Occup= %10.5f" % occupancy[i]
  for j in range(nbas):
    print "%4i %18.8f" % (j+1, orbital[i][order[j]])

file_molden.close()
