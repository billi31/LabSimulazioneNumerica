#!/bin/bash
rm output/*r210unif*
make all
./bin/main --type radius --state 210 --mode unif
