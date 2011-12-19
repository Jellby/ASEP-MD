#!/usr/bin/python

# Use with 2 arguments:
#  1 (read): Input molden file
#  2 (read): Input generic file

import sys
import os
import subprocess

#=============================
# Get input files

try:
  molden_input = sys.argv[1]
except IndexError:
  sys.exit("Missing input molden file")

try:
  gen_input = sys.argv[2]
except IndexError:
  sys.exit("Missing input generic file")

dir = os.path.dirname(sys.argv[0])
molden_input = os.path.relpath(molden_input)
gen_input = os.path.relpath(gen_input)

#=============================
# Read and write the points to calculate the potential

angstrom = 1.88972613289
Q_potential = 0
Q_potential_points = []

file_gen = open(gen_input, "r")
for line in file_gen:
  if (line.rstrip() == "Potential points"):
    Q_potential = int(file_gen.next())
    for i in range(Q_potential):
      tmp = dict(zip(("x","y","z"),file_gen.next().split()))
      tmp["x"] = float(tmp["x"])
      tmp["y"] = float(tmp["y"])
      tmp["z"] = float(tmp["z"])
      Q_potential_points.append(tmp)
file_gen.close()

file_points = open("points.dat", "w")
for point in Q_potential_points:
  print >> file_points, "%20.12f %20.12f %20.12f" % (point["x"]*angstrom, point["y"]*angstrom, point["z"]*angstrom)
file_points.close()

#=============================
# Calculate the electrostatic potential with a patched molden

file_pot = open("pot.in", "w")

print >> file_pot, """\
potential
espch
file=%s
""" % molden_input

file_pot.close()

Q_potential = []

pot_output = subprocess.Popen(["%s/molden-pot" % dir, "pot.in"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()[0]

pot_output = open("potential.dat", "r")
for line in pot_output:
  Q_potential.append(float(line.split()[3]))
print "Electrostatic potential"
for i in range(0, len(Q_potential), 5):
  nums = tuple(Q_potential[i:i+5])
  print (len(nums)*"%20.12E ")[:-1] % nums
print ""

pot_output.close()

os.unlink("points.dat")
os.unlink("potential.dat")
os.unlink("pot.in")
os.unlink("esp.xyz")
os.unlink("plot")
