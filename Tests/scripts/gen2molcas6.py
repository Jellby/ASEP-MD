#!/usr/bin/python

# Use 1 or 2 arguments:
#  1 (read):  generic input file
#  2 (write): Molcas input file (default: append ".mol" to generic input)

import sys

#=============================
# CHANGE THIS DATA AS NEEDED

Q_root = 1

seward_header = """\
*** Usa el directorio de bases (6-31G* con funciones esfericas)
!ln -sf $BASISDIR BASLIB
&SEWARD
  Title
    N-fenilpirrol
  >>> coord inline basis sph-31G*
  NoSymm"""

main_block = """\
!ln -sf $SaveDir/$Project.JobIph JOBOLD
&RASSCF
  JobIph
  NActElectrons = 12 0 0
  Inactive = 32
  Ras2 = 11
  CIRoot = 5 5 1
  RlxRoot = %(root)i
  OutOrbitals = Natural
    5
End of input
!rm -f JOBOLD

!cp $Project.rasscf.%(root)i.molden $Project.esp.molden
!cp $Project.RasOrb%(root_b)s $Project.esp.RasOrb

&RASSI
  NrOfJobIphs = 1  5
    1 2 3 4 5
End of input
""" % {"root": Q_root, "root_b": Q_root if (Q_root > 1) else "" }

alaska = """\
&MCLR
  SALa = %(root)i
  Iterations = 200
End of input

&ALASKA
End of input
""" % {"root": Q_root}

mckinley = """\
&MCKINLEY
End of input
!cat $Project.UnSym
"""

potential = """\
!cp $Project.JobIph $Project.JobOld

&RASSCF
  JobIph
  CIRestart
  CIOnly
  NActElectrons = 12 0 0
  Inactive = 32
  Ras2 = 11
  CIRoot = 1 %(root)i
    %(root)i
End of input

!cp $Project.JobOld $Project.JobIph
""" % {"root": Q_root}

#=============================
# Get input files

try:
  gen_input = sys.argv[1]
except IndexError:
  sys.exit("Missing input file")

try:
  mol_input = sys.argv[2]
except IndexError:
  mol_input = gen_input + ".mol"

#=============================
# Read the data from the generic input file

angstrom = 1.88972613289
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
# Write the Molcas input file

file_mol = open(mol_input, "w")

# Write the SEWARD input
print >> file_mol, seward_header

# Write the molecular geometry, in Angstrom
print >> file_mol, "    %i\n" % Q_natoms
for atom in Q_molecule:
  print >> file_mol, "%16s %20.12f %20.12f %20.12f" % (atom["type"], atom["x"], atom["y"], atom["z"])
print >> file_mol, "  End of input"

# Write the external charges, in Bohr, or the charge file
if (Q_charges):
  print >> file_mol, "  XField"
  print >> file_mol, "    %i 1" % Q_charges
  for charge in Q_charges_points:
    print >> file_mol, "%20.12f %20.12f %20.12f %20.12f %3.1f %3.1f %3.1f" % \
                       (charge["x"]*angstrom, charge["y"]*angstrom, charge["z"]*angstrom, charge["q"], 0, 0, 0)
elif (Q_chargesfile):
  print >> file_mol, "  XField"
  print >> file_mol, "@$MOLCAS_SUBMIT_PWD/%s" % Q_chargesfile

# Write the potential points, in Bohr
# (Skip, to avoid printing the potential for every root)
# (use molden instead)
#if (Q_potential):
#  print >> file_mol, "  EPot"
#  print >> file_mol, "    %i" % Q_potential
#  for point in Q_potential_points:
#    print >> file_mol, "%20.12f %20.12f %20.12f" % (point["x"]*angstrom, point["y"]*angstrom, point["z"]*angstrom)

# End of SEWARD input
print >> file_mol, "End of input\n"

print >> file_mol, main_block
if (Q_derivative >= 1): print >> file_mol, alaska
if (Q_derivative == 2): print >> file_mol, mckinley

# Work-around to avoid useless printing of the potential for every root...
#print >> file_mol, seward_header
#print >> file_mol, "    %i\n" % Q_natoms
#for atom in Q_molecule:
#  print >> file_mol, "%16s %20.12f %20.12f %20.12f" % (atom["type"], atom["x"], atom["y"], atom["z"])
#print >> file_mol, "  End of input"
#if (Q_charges):
#  print >> file_mol, "  XField"
#  print >> file_mol, "    %i 1" % Q_charges
#  for charge in Q_charges_points:
#    print >> file_mol, "%20.12f %20.12f %20.12f %20.12f %3.1f %3.1f %3.1f" % \
#                       (charge["x"]*angstrom, charge["y"]*angstrom, charge["z"]*angstrom, charge["q"], 0, 0, 0)
#elif (Q_chargesfile):
#  print >> file_mol, "  XField"
#  print >> file_mol, "@$MOLCAS_SUBMIT_PWD/%s" % Q_chargesfile
#if (Q_potential):
#  print >> file_mol, "  EPot"
#  print >> file_mol, "    %i" % Q_potential
#  for point in Q_potential_points:
#    print >> file_mol, "%20.12f %20.12f %20.12f" % (point["x"]*angstrom, point["y"]*angstrom, point["z"]*angstrom)
#print >> file_mol, "End of input\n"
#print >> file_mol, potential

file_mol.close()

