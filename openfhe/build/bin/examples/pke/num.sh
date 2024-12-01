#!/bin/bash

numbers="100 500 1000 2000 4000"
for dv in {81..100}
do
    for dt in {3..8}
    do
        for nt in $numbers;
        do
            for it in {0..2}
            do
                ./skylineNum n=$nt d=$dt t=$it v=$dv
                sleep 1
            done
        done
    done
done

