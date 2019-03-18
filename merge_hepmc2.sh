#!/bin/bash

filelist=("$@")

awk '/^$/,/HepMC::IO_GenEvent-START_EVENT_LISTING/' ${filelist[0]}
for file in ${filelist[@]}; do
    awk '/HepMC::IO_GenEvent-START_EVENT_LISTING/{flag=1;next}/HepMC::IO_GenEvent-END_EVENT_LISTING/{flag=0}flag' $file
done
echo "HepMC::IO_GenEvent-END_EVENT_LISTING"
