#!/usr/bin/python

# Use with 2 or 3 arguments:
#  1 (read/write): moldy input file (ctr)
#  2:              bare extension of the current run
#  3:              extension of the previous save file

import sys
import re
import fileinput
import os.path

#=============================
# Get input arguments

try:
  ctr_input = sys.argv[1]
except IndexError:
  sys.exit("Missing input file")

try:
  extension = sys.argv[2]
except IndexError:
  sys.exit("Missing extension")

try:
  oldext = sys.argv[3]
except IndexError:
  num = int(re.search("(.*?)(\d+)$",extension).group(2))
  oldext = re.search("(.*?)(\d+)$",extension).group(1) + str(num-1)

#=============================
# Read the ctr file

lattice = 0

file_ctr = open(ctr_input, "r")
for line in file_ctr:
  if (re.search("sys-spec-file",line)):
    sys_name = line.split("=")[1].lstrip().rstrip()

  if (re.search("save-file",line)):
    save_name = line.split("=")[1].lstrip().rstrip()

  if (re.search("lattice-start",line)):
    lattice = int(line.split("=")[1])

file_ctr.close()

#=============================
# exit if the previous save file does not exist

save_name = save_name.rstrip(extension) + oldext
if (not os.path.exists(save_name)):
  sys.exit()

#=============================
# Add "lattice-start = 1" to the ctr file if needed

if (not lattice):
  file_ctr = fileinput.input(ctr_input, inplace=1)
  for line in file_ctr:
    if (line.lstrip().rstrip() == "end"):
      print " lattice-start = 1"

    print line.rstrip()

  file_ctr.close()

#=============================
# Append the start configuration to the sys-spec-file

file_sys = fileinput.input(sys_name, inplace=1)
file_save = open(save_name, "r")
end = 0
for line in file_sys:
  if (line.lstrip().rstrip() == "end"): end += 1
  print line.rstrip()
  if (end == 2): break

end = 0
for line in file_save:
  if (end >= 3): print line.rstrip()
  if (line.lstrip().rstrip() == "end"): end += 1

file_sys.close()
file_save.close()
