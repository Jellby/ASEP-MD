#!/bin/bash

# $1: input name
# $2: output name
# $3: bare extension

export Project=test
export WorkDir=$SCRDIR/$Project
rm -rf $WorkDir

molcas $1 >& $2
rm -f temp.coord

# calculate ESP charges with molden
if [ -f esp.in ] ; then
  if [ -f $WorkDir/$Project.rasscf.molden ] ; then
    sed -i "s?file=.*\$?file=$WorkDir/$Project.rasscf.molden?" esp.in
   elif [ -f $WorkDir/$Project.scf.molden ] ; then
    sed -i "s?file=.*\$?file=$WorkDir/$Project.scf.molden?" esp.in
  fi
  molden esp.in > esp.out 2> /dev/null
  echo "*** ESP Charges (Molden) ***"     >> $2
  awk '{if(NR>2)print $1 " " $5}' esp.xyz >> $2
  echo "****************************"     >> $2
  rm -f esp.out esp.xyz
fi
