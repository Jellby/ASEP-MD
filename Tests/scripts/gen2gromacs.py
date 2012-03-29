#!/usr/bin/python

# Modify the solute geometry and charges in Gromacs .gro and .top files

# Use with 5 arguments:
#  1 (read): generic system file
#  2 (read): .top file
#  3 (read): .gro file
#  4 (write): modified .top file
#  5 (write): modified .gro file

import sys
import re

#=============================
# Get input arguments

try:
  system_input = sys.argv[1]
except IndexError:
  sys.exit("Missing input file")

try:
  top_input = sys.argv[2]
except IndexError:
  sys.exit("Missing input file")

try:
  gro_input = sys.argv[3]
except IndexError:
  sys.exit("Missing input file")

try:
  top_output = sys.argv[4]
except IndexError:
  sys.exit("Missing output file")

try:
  gro_output = sys.argv[5]
except IndexError:
  sys.exit("Missing output file")

#=============================
# Read the system file

# Skip the file until the solute is found is found
file_system = open(system_input, "r")
for line in file_system:
  if (re.match("Solute",line)):
    break

# Skip name and number of molecules
file_system.next()
file_system.next()

# Read coordinates and charges
mol = []
num = int(file_system.next())
for i in range(num):
  tmp = dict(zip(("x","y","z","q"),file_system.next().split()[4:8]))
  tmp["x"] = float(tmp["x"])
  tmp["y"] = float(tmp["y"])
  tmp["z"] = float(tmp["z"])
  tmp["q"] = float(tmp["q"])
  mol.append(tmp)

# Calculate center of solute molecule
sol_center = { "x": 0.0, "y": 0.0, "z": 0.0 }
for i in range(num):
  sol_center["x"] += mol[i]["x"]
  sol_center["y"] += mol[i]["y"]
  sol_center["z"] += mol[i]["z"]
sol_center["x"] = sol_center["x"]/num
sol_center["y"] = sol_center["y"]/num
sol_center["z"] = sol_center["z"]/num

file_system.close()

#=============================
# Read the topology file
# and write the modified charges

file_top = open(top_input, "r")
file_top_out = open(top_output, "w")

# Skip to the definition of the first molecule's atoms
for line in file_top:
  file_top_out.write(line)
  if (re.match("\[\s*atoms\s*\]",line)):
    break

# Replace the 7th word (the charge) with the new charge
for i in range(num):
  line = file_top.next()
  # Skip comment lines
  while (re.match("\s*;", line)):
    file_top_out.write(line)
    line = file_top.next()
  m = re.match("(\s*\S*\s+\S*\s+\S*\s+\S*\s+\S*\s+\S*\s+)(\S*)(.*)",line)
  file_top_out.write("%s%10.6f%s\n" % (m.group(1), mol[i]["q"], m.group(3)))

# Copy the rest of the file unchanged
for line in file_top:
  file_top_out.write(line)

file_top.close()
file_top_out.close()

#=============================
# Read the coordinates file
# and write the modified coordinates

coord_prec = "11.6"
veloc_prec = "11.6"
format_str = "%%5d%%5s%%5s%%5d%%%sf%%%sf%%%sf%%%sf%%%sf%%%sf\n" % (coord_prec, coord_prec, coord_prec, veloc_prec, veloc_prec, veloc_prec)

file_gro = open(gro_input, "r")
file_gro_out = open(gro_output, "w")

# First read the solute coordinates and calculate its center
file_gro.next()
file_gro.next()
gro_center = { "x": 0.0, "y": 0.0, "z": 0.0 }
for i in range(num):
  tmp = dict(zip(("x","y","z"),file_gro.next()[20:].split()))
  gro_center["x"] += float(tmp["x"])*10
  gro_center["y"] += float(tmp["y"])*10
  gro_center["z"] += float(tmp["z"])*10
gro_center["x"] = gro_center["x"]/num
gro_center["y"] = gro_center["y"]/num
gro_center["z"] = gro_center["z"]/num

# Modify the input coordinates to match centers
# (assuming orientation is the same)
for i in range(num):
  mol[i]["x"] += gro_center["x"]-sol_center["x"]
  mol[i]["y"] += gro_center["y"]-sol_center["y"]
  mol[i]["z"] += gro_center["z"]-sol_center["z"]

# Back to the top of the file
file_gro.seek(0)

# Copy title and total number of atoms
file_gro_out.write(file_gro.next())
numtot = int(file_gro.next())
file_gro_out.write("%5d\n" % numtot)

# Read the atom coordinates and velocities
for i in range(numtot):
  line = file_gro.next()
  tmp = dict(zip(("x","y","z","vx","vy","vz"),line[20:].split()))
  tmp["resnum"] = int(line[0:5])
  tmp["resname"] = line[5:10]
  tmp["atname"] = line[10:15]
  tmp["atnum"] = int(line[15:20])
  # For the solute, write the new coordinates, in nm
  if (i < num):
    tmp["x"] = 0.1*mol[i]["x"]
    tmp["y"] = 0.1*mol[i]["y"]
    tmp["z"] = 0.1*mol[i]["z"]
  else:
    tmp["x"] = float(tmp["x"])
    tmp["y"] = float(tmp["y"])
    tmp["z"] = float(tmp["z"])
  # Write the velocities if present
  if "vx" in tmp:
    tmp["vx"] = float(tmp["vx"])
    tmp["vy"] = float(tmp["vy"])
    tmp["vz"] = float(tmp["vz"])
  else:
    tmp["vx"] = 0.0
    tmp["vy"] = 0.0
    tmp["vz"] = 0.0
  file_gro_out.write(format_str % \
        (tmp["resnum"], tmp["resname"], tmp["atname"], tmp["atnum"], tmp["x"], tmp["y"], tmp["z"], tmp["vx"], tmp["vy"], tmp["vz"]))

# Copy the cell tensor
file_gro_out.write(file_gro.next())

file_gro.close()
file_gro_out.close()
