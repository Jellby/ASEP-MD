#!/usr/bin/python

# Use 1 or 2 arguments:
#  1 (read):  Gaussian output file
#  2 (read):  Gaussian fchk file
#  3 (write): generic output file (default: append ".gen" to Gaussian output)

import sys
import re

#=============================
# Get input files

try:
  gau_output = sys.argv[1]
except IndexError:
  sys.exit("Missing input file")

try:
  gau_fchk = sys.argv[2]
except IndexError:
  sys.exit("Missing input file")

try:
  gen_output = sys.argv[3]
except IndexError:
  gen_output = gau_output + ".gen"

#=============================
# Read the data from the output file(s)

Q_natoms = 0
Q_charge = 0
Q_multiplicity = 0
Q_energy = 0
Q_energy_lower = 0
Q_selfenergy = 0
Q_dipole = [0, 0, 0]
Q_mulliken = []
Q_esp = []
Q_gradient = []
Q_gradient_lower = []
Q_hessian = []
Q_potfile = ""
Q_potential = []

# Read the data available in the fchk file, in atomic units
file_gau = open(gau_fchk, "r")
for line in file_gau:

  if re.match("Number of atoms", line):
    Q_natoms = int(line.rstrip().split()[4])

  elif re.match("Charge", line):
    Q_charge = int(line.rstrip().split()[2])

  elif re.match("Multiplicity", line):
    Q_multiplicity = int(line.rstrip().split()[2])

  elif re.match("Total Energy", line):
    Q_energy = float(line.rstrip().split()[3])

  elif re.match("Dipole Moment", line):
    Q_dipole = map(float, file_gau.next().split())

  elif re.match("Cartesian Gradient", line):
    num = (int(line.split()[4])-1)/5+1
    line = ""
    for i in range(num):
      line += file_gau.next()
    Q_gradient = map(float, line.split())

  elif re.match("Cartesian Force Constants", line):
    num = (int(line.split()[5])-1)/5+1
    line = ""
    for i in range(num):
      line += file_gau.next()
    Q_hessian = map(float, line.split())

file_gau.close()

# Read the data in the Gaussian output file
file_gau = open(gau_output, "r")
for line in file_gau:

  # Self energy of the charges
  if re.search("Self energy of the charges", line):
    Q_selfenergy = float(line.rstrip().split()[6])
    Q_energy -= Q_selfenergy

  # Energy difference and gradients for a conical intersection
  elif re.search("Energy difference", line):
    Q_energy_lower = Q_energy + float(line.split()[2])
  elif re.match("\s*Gradient of iOther State", line):
    del Q_gradient_lower[:]
    for i in range(Q_natoms):
      Q_gradient_lower.extend(map(float, file_gau.next().split()))
  elif re.match("\s*Gradient of iVec State", line):
    del Q_gradient[:]
    for i in range(Q_natoms):
      Q_gradient.extend(map(float, file_gau.next().split()))

  # Mulliken charges
  elif (re.match("\s*Total atomic charges", line) or re.match("\s*Mulliken atomic charges", line)):
    file_gau.next()
    for i in range(Q_natoms):
      Q_mulliken.append(float(file_gau.next().split()[2]))

  # ESP charges
  elif re.search("Charges from ESP fit,", line):
    file_gau.next()
    file_gau.next()
    del Q_esp[:]
    for i in range(Q_natoms):
      Q_esp.append(float(file_gau.next().split()[2]))

  # ESP charges calculated externally (added by a script)
  elif re.search("ESP Charges \(Molden\)", line):
    del Q_esp[:]
    for i in range(Q_natoms):
      Q_esp.append(float(file_gau.next().split()[1]))

  # Fortran unit where the electrostatic potential is written
  elif re.search("Compute potential derivative range", line):
    Q_potfile = int(line.split()[11])

file_gau.close()

# Read the potential
if (Q_potfile):
  potfile = open("fort.%i" % Q_potfile, "r")
  for line in potfile:
    Q_potential.append(float(line.split()[3]))
  potfile.close()
  
#=============================
# Write the generic output file, in atomic units

file_gen = open(gen_output, "w")

print >> file_gen, "Number of atoms\n%4d\n" % Q_natoms
print >> file_gen, "Charge\n%4d\n" % Q_charge
print >> file_gen, "Multiplicity\n%4d\n" % Q_multiplicity
print >> file_gen, "Energy\n%20.12E\n" % Q_energy
if (Q_energy_lower):
  print >> file_gen, "Energy (lower state)\n%20.12E\n" % Q_energy_lower
print >> file_gen, "Dipole moment\n%20.12f %20.12f %20.12f\n" % tuple(Q_dipole)
if (Q_mulliken):
  print >> file_gen, "Mulliken charges"
  for i in range(0, len(Q_mulliken), 1):
    nums = tuple(Q_mulliken[i:i+1])
    print >> file_gen, (len(nums)*"%20.12f ")[:-1] % nums
  print >> file_gen, ""
if (Q_esp):
  print >> file_gen, "ESP charges"
  for i in range(0, len(Q_esp), 1):
    nums = tuple(Q_esp[i:i+1])
    print >> file_gen, (len(nums)*"%20.12f ")[:-1] % nums
  print >> file_gen, ""
if (Q_gradient):
  print >> file_gen, "Cartesian gradient"
  for i in range(0, len(Q_gradient), 5):
    nums = tuple(Q_gradient[i:i+5])
    print >> file_gen, (len(nums)*"%20.12E ")[:-1] % nums
  print >> file_gen, ""
if (Q_gradient_lower):
  print >> file_gen, "Cartesian gradient (lower state)"
  for i in range(0, len(Q_gradient_lower), 5):
    nums = tuple(Q_gradient_lower[i:i+5])
    print >> file_gen, (len(nums)*"%20.12E ")[:-1] % nums
  print >> file_gen, ""
if (Q_hessian):
  print >> file_gen, "Cartesian Hessian"
  for i in range(0, len(Q_hessian), 5):
    nums = tuple(Q_hessian[i:i+5])
    print >> file_gen, (len(nums)*"%20.12E ")[:-1] % nums
  print >> file_gen, ""
if (Q_potential):
  print >> file_gen, "Electrostatic potential"
  for i in range(0, len(Q_potential), 5):
    nums = tuple(Q_potential[i:i+5])
    print >> file_gen, (len(nums)*"%20.12E ")[:-1] % nums
  print >> file_gen, ""

file_gen.close()
