#!/usr/bin/python

# Use with 1 argument:
#  1 (read): Input molden file

import sys
import os
import fileinput
import subprocess

#=============================
# Get input file

try:
  molden_input = sys.argv[1]
except IndexError:
  sys.exit("Missing input molden file")

dir = os.path.dirname(sys.argv[0])
molden_input = os.path.relpath(molden_input)

#=============================
# Calculate ESP charges with molden

file_esp = open("esp.in", "w")

print >> file_esp, """\
esp molcas
espch numsurf=4 connsc=1.4 connincr=0.2 ptden=3.0
file=%s
""" % molden_input

file_esp.close()

esp_output = subprocess.Popen(["molden", "esp.in"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()[0]
esp_output = open("esp.xyz", "r")
esp_output.next()
esp_output.next()
print "ESP charges"
for line in esp_output:
  print "%20.12f" % float(line.split()[4])
print ""

esp_output.close()

os.unlink("esp.in")
os.unlink("esp.xyz")
os.unlink("plot")
