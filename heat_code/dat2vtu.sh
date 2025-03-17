#!/bin/bash

for i in {0..50}
do
    ../../bin/oomph-convert.py RESLT/soln$i.dat
done

pwd