#!/bin/bash

# $1: input name
# $2: output name
# $3: bare extension

if [ $3 == "skip" ]; then
  exit
fi

DIR=$(dirname $0)

$DIR/moldy-conf.py $1 $3

moldy $1 $2
