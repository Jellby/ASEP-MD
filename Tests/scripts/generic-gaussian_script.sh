#!/bin/bash

# $1: input name
# $2: output name
# $3: bare extension

if [ -z "$CPUS" ] ; then CPUS=1 ; fi

DIR=$(dirname $0)

$DIR/gen2gaussian.py $1 $1.gau

gaussian98 $1.gau $2.gau
mv -f Test.FChk fchk$3

$DIR/gaussian2gen.py $2.gau fchk$3 $2
