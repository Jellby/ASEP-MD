#!/usr/bin/python

# Use 1 or 2 arguments:
#  1 (read):  generic input file
#  2 (write): Gaussian input file (default: append ".gau" to generic input)

import sys

#=============================
# CHANGE THIS DATA AS NEEDED

header = """\
#P HF/6-31G*
GFInput iop(6/7=3) Pop=CHELPG
FCheck=All NoSymm"""

jobname = "Job Name"

chargemult = "0 1"

#=============================
# Get input files

try:
  gen_input = sys.argv[1]
except IndexError:
  sys.exit("Missing input file")

try:
  gau_input = sys.argv[2]
except IndexError:
  gau_input = gen_input + ".gau"

#=============================
# Read the data from the generic input file

angstrom = 1.88972613289
Q_potfile1 = 63  #default fortran unit for potential input
Q_potfile2 = 64  #default fortran unit for potential output
Q_derivative = 0
Q_natoms = 0
Q_molecule = []
Q_chargesfile = ""
Q_charges = 0
Q_charges_points = []
Q_potential = 0
Q_potential_points = []

file_gen = open(gen_input, "r")
for line in file_gen:

  # Read the order of the energy derivative needed
  if (line.rstrip() == "Derivative"):
    Q_derivative = int(file_gen.next())

  # Read the geometry
  #  atom name, atomic symbol, atomic number, coordinates in Angstrom
  elif (line.rstrip() == "Geometry"):
    Q_natoms = int(file_gen.next())
    for i in range(Q_natoms):
      tmp = dict(zip(("name","type","number","x","y","z"),file_gen.next().split()))
      tmp["number"] = int(tmp["number"])
      tmp["x"] = float(tmp["x"])
      tmp["y"] = float(tmp["y"])
      tmp["z"] = float(tmp["z"])
      Q_molecule.append(tmp)

  # Read the filename containing external charges, in Angstrom
  elif (line.rstrip() == "External charge file"):
    Q_chargesfile = file_gen.next().rstrip()

  # ... or read the external charges positions and values, in Angstrom
  elif (line.rstrip() == "External charges"):
    Q_charges = int(file_gen.next())
    for i in range(Q_charges):
      tmp = dict(zip(("x","y","z","q"),file_gen.next().split()))
      tmp["x"] = float(tmp["x"])
      tmp["y"] = float(tmp["y"])
      tmp["z"] = float(tmp["z"])
      tmp["q"] = float(tmp["q"])
      Q_charges_points.append(tmp)

  # Read the points where the electrostatic potential is to be calculated, in Angstrom
  elif (line.rstrip() == "Potential points"):
    Q_potential = int(file_gen.next())
    for i in range(Q_potential):
      tmp = dict(zip(("x","y","z"),file_gen.next().split()))
      tmp["x"] = float(tmp["x"])
      tmp["y"] = float(tmp["y"])
      tmp["z"] = float(tmp["z"])
      Q_potential_points.append(tmp)

file_gen.close()

#=============================
# Write the Gaussian input file

file_gau = open(gau_input, "w")

print >> file_gau, header
if (Q_derivative == 1): print >> file_gau, "Force"
if (Q_derivative == 2): print >> file_gau, "Freq"
if (Q_potential): print >> file_gau, "Prop=(Grid,Potential)"
if (Q_charges or Q_chargesfile): print >> file_gau, "Charge=Angstrom"
print >> file_gau, ""
print >> file_gau, jobname
print >> file_gau, ""
print >> file_gau, chargemult

# Molecular geometry, in Angstrom
for atom in Q_molecule:
  print >> file_gau, "%16s %20.12f %20.12f %20.12f" % (atom["type"], atom["x"], atom["y"], atom["z"])
print >> file_gau, ""

# External charges, in Angstrom
if (Q_charges):
  for charge in Q_charges_points:
    print >> file_gau, "%20.12f %20.12f %20.12f %20.12f" % (charge["x"], charge["y"], charge["z"], charge["q"])
  print >> file_gau, ""
elif (Q_chargesfile):
  print >> file_gau, "@%s/N\n" % Q_chargesfile

# Potential points
if (Q_potential):
  # Number of points and fortran units
  print >> file_gau, "%i 3 %i %i\n" % (Q_potential, Q_potfile1, Q_potfile2)
  # Coordinates in the fort.XX file, in Angstrom
  pot_gau = open("fort.%i" % Q_potfile1, "w")
  for point in Q_potential_points:
    print >> pot_gau, "%20.12f %20.12f %20.12f" % (point["x"], point["y"], point["z"])
  pot_gau.close()

file_gau.close()

