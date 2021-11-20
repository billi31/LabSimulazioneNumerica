#!/bin/bash
rm output/*r100unif*
make all
./bin/main --type radius --state 100 --mode unif
