#!/bin/bash

# $1: input name
# $2: output name (unused)
# $3: bare extension
# $4: DCD trajectory file

if [ -z "$CPUS" ] ; then CPUS=1 ; fi

DIR=$(dirname $0)

export FILENAME=PP
export OUTPUT=$FILENAME$3
export DCD=$4

if [ $3 != "skip" ]; then
  $DIR/gen2gromacs.py $1 $FILENAME.top conf.gro $OUTPUT.top $OUTPUT.init.gro &> /dev/null

  grompp -f $FILENAME.mdp -c $OUTPUT.init.gro -p $OUTPUT.top -o $OUTPUT.tpr &> /dev/null && \
  mdrun -nt $CPUS -deffnm $OUTPUT &> /dev/null
fi

echo -e 'SLT\nSystem' | trjconv -f $OUTPUT.cpt -s $OUTPUT.tpr -o conf.gro -pbc mol -center -boxcenter zero -ur compact -ndec 6 &> /dev/null

rm -f $DCD
echo -e 'SLT\nSystem' | trjconv -f $OUTPUT.trr -s $OUTPUT.tpr -o tmp.trr -pbc mol -center -boxcenter zero -ur compact -b 250.5 &> /dev/null
/soft/catdcd/catdcd -otype dcd -o $DCD -trr tmp.trr &> /dev/null

rm -f mdout.mdp tmp.trr
