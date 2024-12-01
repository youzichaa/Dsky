#!/bin/bash

numbers="100 500 1000 2000 4000"
for dt in {3..7}
do
    for tt in {0..4}
    do
        ./skyline2 n=1000 d=$dt t=$tt
        sleep 1
    done
done
for nt in $numbers;
do
    for tt in {0..4}
    do
        ./skyline2 n=$nt d=4 t=$tt
        sleep 1
    done
done

