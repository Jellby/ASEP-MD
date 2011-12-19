#!/usr/bin/python

# Generate a generic system input file from a Moldy save-file

# Use with 1 or 2 arguments:
#  1 (read): moldy save-file
#  2:        index of the solute molecule

import sys
import re

#=============================
# Get input arguments

try:
  save_input = sys.argv[1]
except IndexError:
  sys.exit("Missing input file")

try:
  sol = int(sys.argv[2])
except IndexError:
  sol = 1
sol -= 1

amu = 1.6605402e-27
angstrom = 1.0e-10
electron = 1.60217733e-19
hartree = 4.35974394e-18
bohr = 5.29177208e-11

atnum = {}
atnum["h "] = 1
atnum["he"] = 2
atnum["li"] = 3
atnum["be"] = 4
atnum["b "] = 5
atnum["c "] = 6
atnum["n "] = 7
atnum["o "] = 8
atnum["f "] = 9
atnum["ne"] = 10
atnum["na"] = 11
atnum["mg"] = 12
atnum["al"] = 13
atnum["si"] = 14
atnum["p "] = 15
atnum["s "] = 16
atnum["cl"] = 17
atnum["ar"] = 18

#=============================
# Read the save-file

mass_unit = 1.6605402e-27
length_unit = 1.0e-10
time_unit = 1.0e-13
charge_unit = 1.60217733e-19

# First read the system specification units
file_save = open(save_input, "r")
for line in file_save:
  if (re.search("mass-unit",line)):
    mass_unit = float(line.split("=")[1].lstrip().rstrip())

  if (re.search("length-unit",line)):
    length_unit = float(line.split("=")[1].lstrip().rstrip())

  if (re.search("time-unit",line)):
    time_unit = float(line.split("=")[1].lstrip().rstrip())

  if (re.search("charge-unit",line)):
    charge_unit = float(line.split("=")[1].lstrip().rstrip())

  # Stop at the first "end"
  if (re.match("\s*end",line)):
    break

# Calculate the energy units
energy_unit = mass_unit*length_unit*length_unit/time_unit/time_unit

mol_name = []
mol_number = []
mol_geom = []

line = file_save.next()
i = 0
maxid = 0
# Read molecules until the next "end" is reached
while (not(re.match("\s*end",line))):
  # Read name and number of molecules
  mol_name.append(line.split()[0].lstrip().rstrip())
  mol_number.append(int(line.split()[1].lstrip().rstrip()))
  mol_geom.append([])
  line = file_save.next()
  at = {}
  # Read the atoms in the correct units (coordinates always in angstrom)
  while (re.match("\s*\d+",line.split()[0])):
    tmp = dict(zip(("id","x","y","z","m","q","name"),line.split()))
    tmp["id"] = int(tmp["id"])
    maxid = max(maxid, tmp["id"])
    tmp["x"] = float(tmp["x"])*angstrom
    tmp["y"] = float(tmp["y"])*angstrom
    tmp["z"] = float(tmp["z"])*angstrom
    tmp["m"] = float(tmp["m"])*mass_unit
    tmp["q"] = float(tmp["q"])*charge_unit
    # Try to guess the atomic number from the name
    at = tmp["name"].lower().ljust(2)[0:2]
    if (not(re.search("[a-z]",at[1]))): at = at[0] + " "
    tmp["atnum"] = atnum[at]
    mol_geom[i].append(tmp)
    line = file_save.next()
  i += 1

# Abort if the solute molecule does not exist
if (sol > len(mol_geom)-1):
  sys.exit("There are only %i molecule types" % len(mol_geom))

# Read the Van der Waals parameters
line = file_save.next()
potential = line.split()[0].lstrip().rstrip()
if ((potential.lower() != "lennard-jones") and \
    (potential.lower() != "generic")):
  sys.exit("Potential type not supported: \"%s\"" % potential)
# First dimension the parameter matrix
params = [None]*maxid
for i in range(len(params)):
  params[i] = [None]*maxid
# Read the parameters in the correct units
line = file_save.next()
while (not(re.match("\s*end",line))):
  i = int(line.split()[0])-1
  j = int(line.split()[1])-1
  params[i][j] = map(float,line.split()[2:])
  if (potential.lower() == "lennard-jones"):
    params[i][j][0] = params[i][j][0]*energy_unit
    params[i][j][1] = params[i][j][1]*length_unit
  elif (potential.lower() == "generic"):
    params[i][j][0] = params[i][j][0]*energy_unit
    params[i][j][1] = params[i][j][1]/length_unit
    params[i][j][2] = params[i][j][2]*energy_unit*length_unit**12
    params[i][j][3] = params[i][j][3]*energy_unit*length_unit**4
    params[i][j][4] = params[i][j][4]*energy_unit*length_unit**6
    params[i][j][5] = params[i][j][5]*energy_unit*length_unit**8
  params[j][i] = params[i][j]
  line = file_save.next()

file_save.close()

#=============================
# Print the generic system file

print "# PLEASE CHECK EVERYTHING IS CORRECT"
print

# Print the solute data (in angstrom, amu, and e)
print "Solute"
print "  %s" % mol_name[sol]
print "  %i" % mol_number[sol]
print "  %i" % len(mol_geom[sol])
for at in mol_geom[sol]:
  print "  %2i %16s %3i %5.2f %10.6f %10.6f %10.6f %10.6f" % \
  (at["id"], at["name"].ljust(16), at["atnum"], at["m"]/amu, at["x"]/angstrom, at["y"]/angstrom, at["z"]/angstrom, at["q"]/electron)
print

# Calculate the number of solvent atoms
num = 0
for i in range(len(mol_geom)):
  num += mol_number[i]*len(mol_geom[i])
num -= len(mol_geom[sol])
# Print the solvent data (in e)
print "Solvent"
print "  %i" % num
for i in range(len(mol_geom)):
  if (i == sol): continue
  for j in range(mol_number[i]):
    for at in mol_geom[i]:
      print "  %2i %16s %10.6f" % (at["id"], at["name"].ljust(16), at["q"]/electron)
print

# Print only the Van der Waals parameters between solute and solvent
sol_ids = []
nosol_ids = []
for i in range(len(mol_geom)):
  for j in range(mol_number[i]):
    for k in range(len(mol_geom[i])):
      at = mol_geom[i][k]
      if ((i == sol) and (j == 0)):
        sol_ids.append(at["id"]-1)
      else:
        nosol_ids.append(at["id"]-1)
num = 0
# First count the number of pairs
for i in range(len(params)):
  for j in range(len(params[i])):
    if ((i in sol_ids) and (j in nosol_ids) and params[i][j]): num += 1
# Then print the values (in atomic units)
print "Non-Bonded"
print "  %s" % potential
print "  %i" % num
if (potential.lower() == "lennard-jones"):
  for i in range(len(params)):
    for j in range(len(params[i])):
      if ((i in sol_ids) and (j in nosol_ids) and params[i][j]):
        print "  %2i %2i %14.10f %14.10f" % (i+1, j+1, params[i][j][0]/hartree, params[i][j][1]/bohr)
elif (potential.lower() == "generic"):
  for i in range(len(params)):
    for j in range(len(params[i])):
      if ((i in sol_ids) and (j in nosol_ids) and params[i][j]):
        print "  %2i %2i %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f" % \
              (i+1, j+1, params[i][j][0]/hartree, params[i][j][1]*bohr, \
              params[i][j][2]/hartree/bohr**12, params[i][j][3]/hartree/bohr**4, \
              params[i][j][4]/hartree/bohr**6, params[i][j][5]/hartree/bohr**8)

