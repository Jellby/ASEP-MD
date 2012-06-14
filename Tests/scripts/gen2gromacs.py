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
import math
import copy

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
# Function to superpose molecules
# see:   Acta Chrystallogr. Sec. A 61 (2005), 478
#        J. Comput. Chem. 31 (2010), 1561

def superpose ( mol1, mol2 ):

  center1 = { "x": 0.0, "y": 0.0, "z": 0.0 }
  for i in range(len(mol1)):
    center1["x"] += mol1[i]["x"]
    center1["y"] += mol1[i]["y"]
    center1["z"] += mol1[i]["z"]
  center1["x"] = center1["x"]/len(mol1)
  center1["y"] = center1["y"]/len(mol1)
  center1["z"] = center1["z"]/len(mol1)
  for i in range(len(mol1)):
    mol1[i]["x"] -= center1["x"]
    mol1[i]["y"] -= center1["y"]
    mol1[i]["z"] -= center1["z"]

  G1 = 0
  for i in range(len(mol1)):
    G1 += mol1[i]["x"]**2+mol1[i]["y"]**2+mol1[i]["z"]**2

  center2 = { "x": 0.0, "y": 0.0, "z": 0.0 }
  for i in range(len(mol2)):
    center2["x"] += mol2[i]["x"]
    center2["y"] += mol2[i]["y"]
    center2["z"] += mol2[i]["z"]
  center2["x"] = center2["x"]/len(mol2)
  center2["y"] = center2["y"]/len(mol2)
  center2["z"] = center2["z"]/len(mol2)
  for i in range(len(mol2)):
    mol2[i]["x"] -= center2["x"]
    mol2[i]["y"] -= center2["y"]
    mol2[i]["z"] -= center2["z"]

  G2 = 0
  for i in range(len(mol2)):
    G2 += mol2[i]["x"]**2+mol2[i]["y"]**2+mol2[i]["z"]**2

  M = {}
  for i in ["x", "y", "z"]:
    for j in ["x", "y", "z"]:
      M[i+j] = 0
      for k in range(len(mol1)):
        M[i+j] += mol1[k][i] * mol2[k][j]

  K = []
  K.append( [ M["xx"]+M["yy"]+M["zz"], M["yz"]-M["zy"], M["zx"]-M["xz"], M["xy"]-M["yx"] ] )
  K.append( [ M["yz"]-M["zy"], M["xx"]-M["yy"]-M["zz"], M["xy"]+M["yx"], M["xz"]+M["zx"] ] )
  K.append( [ M["zx"]-M["xz"], M["xy"]+M["yx"], M["yy"]-M["xx"]-M["zz"], M["yz"]+M["zy"] ] )
  K.append( [ M["xy"]-M["yx"], M["xz"]+M["yz"], M["yz"]+M["zy"], M["zz"]-M["xx"]-M["yy"] ] )

  coef = []
  D = (M["xy"]**2+M["xz"]**2-M["yx"]**2-M["zx"]**2)**2
  E = (-M["xx"]**2+M["yy"]**2+M["zz"]**2+M["yz"]**2+M["zy"]**2-2*(M["yy"]*M["zz"]-M["yz"]*M["zy"]))*\
      (-M["xx"]**2+M["yy"]**2+M["zz"]**2+M["yz"]**2+M["zy"]**2+2*(M["yy"]*M["zz"]-M["yz"]*M["zy"]))
  F = (-(M["xz"]+M["zx"])*(M["yz"]-M["zy"])+(M["xy"]-M["yx"])*(M["xx"]-M["yy"]-M["zz"]))*\
      (-(M["xz"]-M["zx"])*(M["yz"]+M["zy"])+(M["xy"]-M["yx"])*(M["xx"]-M["yy"]+M["zz"]))
  G = (-(M["xz"]+M["zx"])*(M["yz"]+M["zy"])-(M["xy"]+M["yx"])*(M["xx"]+M["yy"]-M["zz"]))*\
      (-(M["xz"]-M["zx"])*(M["yz"]-M["zy"])-(M["xy"]+M["yx"])*(M["xx"]+M["yy"]+M["zz"]))
  H = ( (M["xy"]+M["yx"])*(M["yz"]+M["zy"])+(M["xz"]+M["zx"])*(M["xx"]-M["yy"]+M["zz"]))*\
      (-(M["xy"]-M["yx"])*(M["yz"]-M["zy"])+(M["xz"]+M["zx"])*(M["xx"]+M["yy"]+M["zz"]))
  I = ( (M["xy"]+M["yx"])*(M["yz"]-M["zy"])+(M["xz"]-M["zx"])*(M["xx"]-M["yy"]-M["zz"]))*\
      (-(M["xy"]-M["yx"])*(M["yz"]+M["zy"])+(M["xz"]-M["zx"])*(M["xx"]+M["yy"]-M["zz"]))
  coef.append( D+E+F+G+H+I )
  coef.append( -8.0*( M["xx"]*M["yy"]*M["zz"]+M["xy"]*M["yz"]*M["zx"]+M["xz"]*M["yx"]*M["zy"]
                     -M["xx"]*M["yz"]*M["zy"]-M["xy"]*M["yx"]*M["zz"]-M["xz"]*M["yy"]*M["zx"] ) )
  coef.append( -2.0*( M["xx"]**2+M["xy"]**2+M["xz"]**2+M["yx"]**2+M["yy"]**2+M["yz"]**2+M["zx"]**2+M["zy"]**2+M["zz"]**2 ) )
  coef.append( 0.0 )
  coef.append( 1.0 )
  
  root_old = 0.0
  root = 0.5*(G1+G2)
  while (math.fabs(root-root_old) > 1.0e-6):
    root_old = root
    P = root**4+coef[2]*root**2+coef[1]*root+coef[0]
    dP = 4*root**3+2*coef[2]*root+coef[1]
    root -= P/dP

  for i in range(len(K)):
    K[i][i] -= root

  for i in range(len(K)):
    vect = []
    for j in range(len(K)):
      adj = copy.deepcopy(K)
      del adj[i]
      for k in range(len(adj)):
        del adj[k][j]
      det = adj[0][0]*adj[1][1]*adj[2][2]+adj[0][1]*adj[1][2]*adj[2][0]+adj[0][2]*adj[1][0]*adj[2][1] \
           -adj[0][0]*adj[1][2]*adj[2][1]-adj[0][1]*adj[1][0]*adj[2][2]-adj[0][2]*adj[1][1]*adj[2][0]
      det *= (-1)**(i+j)
      vect.append(det)
    norm = math.sqrt(vect[0]**2+vect[1]**2+vect[2]**2+vect[3]**2)
    if (norm > 1.0e-6):
      vect[0] = -vect[0]/norm
      vect[1] = vect[1]/norm
      vect[2] = vect[2]/norm
      vect[3] = vect[3]/norm
      break

  M["xx"] =vect[0]**2+vect[1]**2-vect[2]**2-vect[3]**2
  M["yy"] =vect[0]**2-vect[1]**2+vect[2]**2-vect[3]**2
  M["zz"] =vect[0]**2-vect[1]**2-vect[2]**2+vect[3]**2
  M["xy"] =2.0*(vect[1]*vect[2]-vect[0]*vect[3])
  M["yx"] =2.0*(vect[1]*vect[2]+vect[0]*vect[3])
  M["yz"] =2.0*(vect[2]*vect[3]-vect[0]*vect[1])
  M["zy"] =2.0*(vect[2]*vect[3]+vect[0]*vect[1])
  M["zx"] =2.0*(vect[1]*vect[3]-vect[0]*vect[2])
  M["xz"] =2.0*(vect[1]*vect[3]+vect[0]*vect[2])

  old = copy.deepcopy(mol2)
  for i in range(len(mol2)):
    mol2[i]["x"] = M["xx"]*old[i]["x"]+M["xy"]*old[i]["y"]+M["xz"]*old[i]["z"]+center1["x"]
    mol2[i]["y"] = M["yx"]*old[i]["x"]+M["yy"]*old[i]["y"]+M["yz"]*old[i]["z"]+center1["y"]
    mol2[i]["z"] = M["zx"]*old[i]["x"]+M["zy"]*old[i]["y"]+M["zz"]*old[i]["z"]+center1["z"]

  return

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

# First read the solute coordinates
file_gro.next()
file_gro.next()
mol_gro = []
for i in range(num):
  line = file_gro.next()
  dots = [match.start() for match in re.finditer("\.", line[20:])]
  width = dots[1]-dots[0]
  tmp = dict(zip(("x","y","z"), [line[i:i+width] for i in range(20, len(line), width)]))
  tmp["x"] = float(tmp["x"])*10
  tmp["y"] = float(tmp["y"])*10
  tmp["z"] = float(tmp["z"])*10
  mol_gro.append(tmp)

# Modify the input coordinates to fit the original orientation
superpose ( mol_gro, mol )

# Back to the top of the file
file_gro.seek(0)

# Copy title and total number of atoms
file_gro_out.write(file_gro.next())
numtot = int(file_gro.next())
file_gro_out.write("%5d\n" % numtot)

# Read the atom coordinates and velocities
for i in range(numtot):
  line = file_gro.next()
  dots = [match.start() for match in re.finditer("\.", line[20:])]
  width = dots[1]-dots[0]
  tmp = dict(zip(("x","y","z","vx","vy","vz"), [line[i:i+width] for i in range(20, len(line), width)]))
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
