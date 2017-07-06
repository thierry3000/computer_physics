#!/bin/bash
if [ $# -eq 0 ]; then
echo "Chose exercise: 0-2"
exit -1;
fi

## further information, see compile_infos.txt

gcc -I/usr/include/cairo -lm -lcairo -o mdbasic_$1 -O2 mdbasic_$1.c && echo "programm \"mdbasic_$1\" created"
