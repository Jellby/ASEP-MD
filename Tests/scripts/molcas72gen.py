#!/usr/bin/python

# Use 1 or 2 arguments:
#  1 (read):  Molcas output file
#  2 (write): generic output file (default: append ".gen" to Molcas output)

import sys
import re

#=============================
# Get input files

try:
  mol_output = sys.argv[1]
except IndexError:
  sys.exit("Missing input file")

try:
  gen_output = sys.argv[2]
except IndexError:
  gen_output = mol_output + ".gen"

#=============================
# Read the data from the output file(s)

# Change to appropriate value if searching a CI
lower_root = 2

debye = 0.3934303073
Q_natoms = 0
Q_charge = 0
Q_multiplicity = 0
Q_energy = 0
Q_energy_lower = 0
Q_dipole = [0, 0, 0]
Q_mulliken = []
Q_esp = []
Q_gradient = []
Q_gradient_lower = []
Q_hessian = []
Q_potfile = ""
Q_potential = []
Q_root = 0

file_mol = open(mol_output, "r")
rasscf_root = 0
caspt2_root = 0
module = ""
for line in file_mol:
  match = re.search("Start Module: (\w+)", line)
  if match: module = match.group(1)

  # Read the data available from each module's output
  # (later modules overwrite earlier ones)

  if (module == "gateway"):
    if (re.search("Cartesian Coordinates", line)):
      file_mol.next()
      file_mol.next()
      file_mol.next()
      Q_natoms = -1
      while (not re.match("\s*$", line)):
        Q_natoms += 1
        line = file_mol.next()

  elif (module == "scf"):
    if (re.match("\s*Molecular charge", line)):
      Q_charge = float(line.rstrip().split()[2])
    elif (re.match("\s*Total spin, S ", line)):
      Q_multiplicity = int(float(line.rstrip().split()[3])*2+1.5)
    elif (re.match("\s*Total (SCF|KS-DFT) energy", line)):
      Q_energy = float(line.rstrip().split()[3])
    elif (re.match("\s*Dipole Moment \(Debye\):", line)):
      file_mol.next()
      Q_dipole = map(lambda x: float(x)*debye, file_mol.next().rstrip().split()[1:6:2])
    elif (re.match("\s*Mulliken charges per cent(er|re) and basis function type", line)):
      del Q_mulliken[:]
      file_mol.next()
      file_mol.next()
      file_mol.next()
      while (not re.match("\s*Total electronic charge", line)):
        while (not re.match("\s*$", line)):
          line = file_mol.next()
        Q_mulliken.extend(map(float,file_mol.next().split()[1:]))
        file_mol.next()
        line = file_mol.next()
    elif (re.match("\s*Electric Potential:", line)):
      if (re.search("centre no.\s*(\d+)\s*\(", line).group(1) == "1"): del Q_potential[:]
      file_mol.next()
      file_mol.next()
      file_mol.next()
      file_mol.next()
      Q_potential.append(float(file_mol.next().split()[1]))
    elif (re.match("\s*ESPF analysis", line)):
      del Q_esp[:]
      while (not re.match("\s*Total ESPF QM/MM interaction energy", line)):
        line = file_mol.next()
        if (re.match("\s*Charge on", line)):
          Q_esp.append(float(line.split()[4]))
     
  elif (module == "rasscf"):
    if (re.match("\s*Total molecular charge", line)):
      Q_charge = float(line.rstrip().split()[3])
    elif (re.match("\s*Spin quantum number", line)):
      Q_multiplicity = int(float(line.rstrip().split()[3])*2+1.5)
    elif (re.match("\s*Root passed to geometry opt.", line)):
      Q_root = int(line.rstrip().split()[5])
    elif (re.match("\s*Final state energy\(ies\):", line)):
      file_mol.next()
      file_mol.next()
      while (not re.match("\s*$", line)):
        line = file_mol.next()
        match = re.search("root number\s+(\d+)", line)
        if (match and (int(match.group(1)) == Q_root)):
          Q_energy = float(line.rstrip().split()[5])
        if (match and (int(match.group(1)) == lower_root)):
          Q_energy_lower = float(line.rstrip().split()[5])
    # Read values only for the interesting root
    elif (re.match("\s*Mulliken population Analysis for root number:", line)):
      rasscf_root = int(line.rstrip().split()[6])
    if (rasscf_root != Q_root): continue
    if (re.match("\s*Dipole Moment \(Debye\):", line)):
      file_mol.next()
      Q_dipole = map(lambda x: float(x)*debye, file_mol.next().rstrip().split()[1:6:2])
    elif (re.match("\s*Mulliken charges per cent(er|re) and basis function type", line)):
      del Q_mulliken[:]
      file_mol.next()
      file_mol.next()
      file_mol.next()
      while (not re.match("\s*Total electronic charge", line)):
        while (not re.match("\s*$", line)):
          line = file_mol.next()
        Q_mulliken.extend(map(float,file_mol.next().split()[1:]))
        file_mol.next()
        line = file_mol.next()
    elif (re.match("\s*Electric Potential:", line)):
      if (re.search("centre no.\s*(\d+)\s*\(", line).group(1) == "1"): del Q_potential[:]
      file_mol.next()
      file_mol.next()
      file_mol.next()
      file_mol.next()
      Q_potential.append(float(file_mol.next().split()[1]))
    elif (re.match("\s*ESPF analysis", line)):
      del Q_esp[:]
      while (not re.match("\s*Total ESPF QM/MM interaction energy", line)):
        line = file_mol.next()
        if (re.match("\s*Charge on", line)):
          Q_esp.append(float(line.split()[4]))

  elif (module == "caspt2"):
    # Read values only for the interesting root
    if (re.match("\s*Single-state initialization phase begins for state", line)):
      caspt2_root = int(line.rstrip().split()[6])
#IFG revisar
    if ((caspt2_root == lower_root) and (re.match("\s*Total energy:", line))):
      Q_energy_lower = float(line.rstrip().split()[2])
    if (caspt2_root != Q_root): continue
    if (re.match("\s*Total energy:", line)):
      Q_energy = float(line.rstrip().split()[2])
    elif (re.match("\s*Dipole Moment \(Debye\):", line)):
      file_mol.next()
      Q_dipole = map(lambda x: float(x)*debye, file_mol.next().rstrip().split()[1:6:2])
    elif (re.match("\s*Mulliken charges per cent(er|re) and basis function type", line)):
      del Q_mulliken[:]
      file_mol.next()
      file_mol.next()
      file_mol.next()
      while (not re.match("\s*Total electronic charge", line)):
        while (not re.match("\s*$", line)):
          line = file_mol.next()
        Q_mulliken.extend(map(float,file_mol.next().split()[1:]))
        file_mol.next()
        line = file_mol.next()
    elif (re.match("\s*Electric Potential:", line)):
      if (re.search("centre no.\s*(\d+)\s*\(", line).group(1) == "1"): del Q_potential[:]
      file_mol.next()
      file_mol.next()
      file_mol.next()
      file_mol.next()
      Q_potential.append(float(file_mol.next().split()[1]))

  elif (module == "alaska"):
    if (re.match("\s*\*\s*Molecular gradients", line)):
      del Q_gradient[:]
      file_mol.next()
      file_mol.next()
      file_mol.next()
      file_mol.next()
      file_mol.next()
      for i in range(3*Q_natoms):
        Q_gradient.append(float(file_mol.next().split()[2]))

  # ESP charges calculated externally (added by a script)
  if re.search("ESP Charges \(Molden\)", line):
    del Q_esp[:]
    for i in range(Q_natoms):
      Q_esp.append(float(file_mol.next().split()[1]))

  # Hessian, copied from UnSym file
  elif re.match("\*BEGIN HESSIAN", line):
    del Q_hessian[:]
    num = (int(file_mol.next().split()[3])-1)/4+1
    for i in range(3*Q_natoms):
      file_mol.next()
      line = ""
      for j in range(num):
        line += file_mol.next()
      Q_hessian.extend(map(float, line.split()[0:i+1]))

file_mol.close()

# # Read the data in the Gaussian output file
# 
#   # Energy difference and gradients for a conical intersection
#   elif re.search("Energy difference", line):
#     Q_energy_lower = Q_energy + float(line.split()[2])
#   elif re.match("\s*Gradient of iOther State", line):
#     del Q_gradient_lower[:]
#     for i in range(Q_natoms):
#       Q_gradient_lower.extend(map(float, file_gau.next().split()))
#   elif re.match("\s*Gradient of iVec State", line):
#     del Q_gradient[:]
#     for i in range(Q_natoms):
#       Q_gradient.extend(map(float, file_gau.next().split()))
  
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
