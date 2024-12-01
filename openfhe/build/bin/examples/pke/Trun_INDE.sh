#!/bin/bash

numbers="100 500 1000 2000 4000"

for dt in {3..7}
do
    ./skylineTime-2T n=1000 d=$dt t=2
    sleep 1
done
for nt in $numbers;
do
    ./skylineTime-2T n=$nt d=4 t=2
    sleep 1
done


