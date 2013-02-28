#!/usr/bin/python

# Generate a generic system input file from a Gromacs dump file

# Use with 1 or 2 arguments:
#  1 (read): gromacs dump (from gmxdump -s *.tpr)
#  2:        index of the solute molecule, or indexes (ini:fin) of the solute atoms (1-based)

import sys
import re

#=============================
# Get input arguments

try:
  dump_input = sys.argv[1]
except IndexError:
  sys.exit("Missing input file")

try:
  sol = sys.argv[2]
  match = re.match("(\d+):(\d+)",sol)
  # Set the solute molecule, or the initial and final solute atoms
  if (re.match("(\d+):(\d+)",sol)):
    sol_ini = int(match.group(1))-1
    sol_fin = int(match.group(2))-1
    sol = 1
  else:
    sol_ini = None
    sol_fin = None
  sol = int(sol)-1
except IndexError:
  sol_ini = None
  sol_fin = None
  sol = 0

amu = 1.6605402e-27
angstrom = 1.0e-10
electron = 1.60217733e-19
hartree = 4.35974394e-18
bohr = 5.29177208e-11

#=============================
# Read the dump file

# Units in Gromacs
mass_unit = 1.66053878316273e-27
length_unit = 1.0e-9
energy_unit = 1.66053878316273e-21
charge_unit = 1.602176487e-19

mol_name = []
mol_number = []
mol_geom = []

# Skip the file until "topology" is found
file_dump = open(dump_input, "r")
for line in file_dump:
  if (re.match("topology:",line)):
    break

# If solute given as a whole molecule, search the numbers of molecules
if (sol_ini == None):
  for line in file_dump:
    if (re.search("#molecules",line)):
      mol_number.append(int(line.split()[2].lstrip().rstrip()))
    if (re.search("moltype \(",line)): break
# Else, just assume everything is a single molecule
else:
  mol_number.append(1)

# Read the atom data for each molecule
maxid = 0
for line in file_dump:
  if (re.match("^\s*name=",line)):
    mol_name.append(re.search("\"(.*)\"",line).group(1))
  if (re.match("\s*atoms:",line)):
    line = file_dump.next()
    num = int(re.search("atom \((\d+)\)",line).group(1))
    mol = []
    for i in range(num):
      line = file_dump.next()
      tmp = {}
      tmp["nt"] = int(re.search("atom\[\s*(\d+)\]",line).group(1))
      tmp["id"] = int(re.search("type=([^,}]*)",line).group(1))+1
      tmp["m"] = float(re.search("m=([^,}]*)",line).group(1))*mass_unit
      tmp["q"] = float(re.search("q=([^,}]*)",line).group(1))*charge_unit
      tmp["atnum"] = int(re.search("atomnumber=([^,}]*)",line).group(1))
      mol.append(tmp)
    line = file_dump.next()
    num = int(re.search("atom \((\d+)\)",line).group(1))
    names = []
    # the atom names are stored apart
    for i in range(num):
      line = file_dump.next()
      names.append(re.search("\"(.*)\"",line).group(1))
    for at in mol:
      at["name"] = names[at["nt"]]
      maxid = max(maxid,at["id"])
    mol_geom.append(mol)
  if (re.match("x \(",line)): break

# Abort if solute indexes were given, but the dump is not for the whole system (with -sys)
if (sol_ini != None):
  if ((len(mol_geom) > 1) or (mol_number[0] > 1)):
    sys.exit("Give a single solute number, or generate the dump with -sys")

# Abort if solute indexes were not given, but the dump is for the whole system (with -nosys)
if (sol_ini == None):
  if (len(mol_number) < 1):
    sys.exit("Give two solute indexes, or generate the dump with -nosys")

# Abort if the solute molecule does not exist
if (sol > sum(mol_number)-1):
  sys.exit("There are only %i molecules in the system" % sum(mol_number))

# Read the coordinates of the solute molecule
k=-1
for i in range(len(mol_geom)):
  for j in range(mol_number[i]):
    k += 1
    for at in mol_geom[i]:
      line = file_dump.next()
      if (k == sol):
        pos = re.search("{([^,}]*),([^,}]*),([^,}]*)}",line)
        at["x"] = float(pos.group(1))*length_unit
        at["y"] = float(pos.group(2))*length_unit
        at["z"] = float(pos.group(3))*length_unit

file_dump.close()

# Dimension the Van der Waals parameter matrix
params = [None]*maxid
for i in range(len(params)):
  params[i] = [None]*maxid

# Read the parameters
file_dump = open(dump_input, "r")
for line in file_dump:
  if (re.match("\s*ntypes=",line)): break
potential = None
for i in range(maxid):
  for j in range(maxid):
    line = file_dump.next()
    params[i][j] = []
    if (re.search("LJ_SR",line)):
      if ((potential != None) and (potential.lower() != "lennard-jones")):
        sys.exit("Mixed parameter types not supported")
      potential = "Lennard-Jones"
      pars = re.search("c6=([^,}]*), c12=([^,}]*)",line)
      a = float(pars.group(1))
      b = float(pars.group(2))
      if ((a == 0.0) or (b == 0.0)):
        params[i][j] = None
      else:
        params[i][j].append(a*a/b*energy_unit)
        params[i][j].append((b/a)**(1.0/6)*length_unit)
    elif (re.search("BHAM",line)):
      if ((potential != None) and (potential.lower() != "generic")):
        sys.exit("Mixed parameter types not supported")
      potential = "Generic"
      pars = re.search("a=([^,}]*), b=([^,}]*), c=([^,}]*)",line)
      a = float(pars.group(1))
      b = float(pars.group(2))
      c = float(pars.group(3))
      if ((a == 0.0) and (c == 0.0)):
        params[i][j] = None
      else:
        params[i][j].append(a*energy_unit)
        params[i][j].append(b/length_unit)
        params[i][j].append(0.0)
        params[i][j].append(0.0)
        params[i][j].append(c*energy_unit*length_unit**6)
        params[i][j].append(0.0)
    else:
      sys.exit("Parameter type not supported")

file_dump.close()

# Read exclusions and 1-4 parameters if solute indexes are given
exclusions = {}
pairs = {}
if (sol_ini != None):

  file_dump = open(dump_input, "r")
  for line in file_dump:
    match = re.match("\s*excls\[(\d+)\]\[.*\]={(.*?)(}?)$",line)
    if (match):
      atom = int(match.group(1))
      # Read only solute exclusions
      if ((atom < sol_ini) or (atom > sol_fin)): continue
      excllist = match.group(2)
      while (match.group(3) == ""):
        match = re.match("(\s*)(.*?)(}?)$",file_dump.next())
        excllist += match.group(2)
      excllist = map(int,excllist.split(", "))
      # Convert the indexes to 1-based, without solute
      # and remove solute-solute exclusions
      for i in reversed(range(len(excllist))):
        if ((excllist[i] < sol_ini) or (excllist[i] > sol_fin)):
          if (excllist[i] > sol_fin): excllist[i] -= sol_fin-sol_ini
        else:
          del excllist[i]
      if (excllist): exclusions[atom-sol_ini+1] = excllist
  file_dump.close()

  file_dump = open(dump_input, "r")
  for line in file_dump:
    if (re.match("\s*LJ-14",line)): break
  num = int(file_dump.next().split()[1])/3
  file_dump.next()
  for i in range(num):
    pair = map(int, file_dump.next().split()[3:5])
    check = map(lambda x: (x < sol_ini) or (x > sol_fin), pair)
    # Remove solute-solute and solvent-solvent pairs
    # and convert the indexes to 1-based, without solute
    if (check[0]):
      if (not check[1]):
        a = pair[1]-sol_ini+1
        b = pair[0]
    elif (check[1]):
      a = pair[0]-sol_ini+1
      b = pair[1]
    if (check[0] != check[1]):
      if (b > sol_fin): b -= sol_fin-sol_ini
      if (not a in pairs): pairs[a] = []
      pairs[a].append(b)

  file_dump.close()

#=============================
# Print the generic system file

print "# PLEASE CHECK EVERYTHING IS CORRECT"
print

if (exclusions):
  print "# solute-solvent exclusions (not converted):"
  for i in sorted(exclusions.keys()):
    if i in exclusions: print "#", i, exclusions[i]
  print

if (pairs):
  print "# solute-solvent pairs (not converted):"
  for i in sorted(pairs.keys()):
    if i in pairs: print "#", i, pairs[i]
  print

# Print the solute data (in angstrom, amu, and e)
print "Solute"
k = -1
for i in range(len(mol_geom)):
  k += mol_number[i]
  if (k < sol): continue
  if (sol_ini == None): sol_ini = 0
  if (sol_fin == None): sol_fin = len(mol_geom[i])-1
  print "  %s" % re.sub(" ","_",mol_name[i])
  print "  %i" % 1
  print "  %i" % (sol_fin-sol_ini+1)
  for j in range(sol_ini,sol_fin+1):
    at = mol_geom[i][j]
    print "  %2i %16s %3i %7.4f %10.6f %10.6f %10.6f %10.6f" % \
    (at["id"], at["name"].ljust(16), at["atnum"], at["m"]/amu, at["x"]/angstrom, at["y"]/angstrom, at["z"]/angstrom, at["q"]/electron)
  break
print

# Calculate the number of solvent atoms
num = 0
for i in range(len(mol_geom)):
  num += mol_number[i]*len(mol_geom[i])
num -= sol_fin-sol_ini+1
# Print the solvent data (in e)
print "Solvent"
print "  %i" % num
l = -1
for i in range(len(mol_geom)):
  for j in range(mol_number[i]):
    l += 1
    for k in range(len(mol_geom[i])):
      if ((l == sol) and (k >= sol_ini) and (k <= sol_fin)): continue
      at = mol_geom[i][k]
      print "  %2i %16s %10.6f" % (at["id"], at["name"].ljust(16), at["q"]/electron)
print

# Print only the Van der Waals parameters between solute and solvent
sol_ids = []
nosol_ids = []
l = -1
for i in range(len(mol_geom)):
  for j in range(mol_number[i]):
    l += 1
    for k in range(len(mol_geom[i])):
      at = mol_geom[i][k]
      if ((l == sol) and (k >= sol_ini) and (k <= sol_fin)):
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
print "  %s" % potential.lower()
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

