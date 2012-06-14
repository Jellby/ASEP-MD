#!/bin/bash

# $1: input name
# $2: output name
# $3: bare extension

DIR=$(dirname $0)

[ -n "$PBS_O_WORKDIR" ] || PBS_O_WORKDIR=$PWD
cd $PBS_O_WORKDIR

export Project=Test
export HomeDir=$PBS_O_WORKDIR
export WorkDir=$HomeDir/Work
export SaveDir=$HomeDir/Save

export MOLCASMEM=8000
export MOLCAS=/home/todos/Programas/molcas-6.4
export BASISDIR=/home/todos/BasesMolcas

$DIR/gen2molcas6.py $1 $1.mol

rm -f $WorkDir/$Project.RunFile
molcas $1.mol >& $2.mol
rm -f temp.coord

cp $WorkDir/$Project.JobIph $SaveDir
cp $WorkDir/$Project.JobIph $SaveDir/$Project.JobIph$3
cp $WorkDir/$Project.rasscf.molden $SaveDir/$Project.rasscf.molden$3
cp $WorkDir/$Project.esp.molden $SaveDir/$Project.esp.molden$3
cp $WorkDir/$Project.esp.RasOrb $SaveDir/$Project.esp.RasOrb$3

$DIR/molcas62gen.py $2.mol $2

$DIR/rasorb2molden.py $SaveDir/$Project.esp.RasOrb$3 $SaveDir/$Project.esp.molden$3
$DIR/pot-molcas.py $SaveDir/$Project.esp.molden$3 $1 >> $2
$DIR/esp-molcas.py $SaveDir/$Project.esp.molden$3 >> $2
