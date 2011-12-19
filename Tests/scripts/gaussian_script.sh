#!/bin/bash

# $1: input name
# $2: output name
# $3: bare extension

gaussian98 $1 $2
mv -f Test.FChk fchk$3
