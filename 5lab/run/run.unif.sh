#!/bin/bash
rm output/*.r*unif*
make all
./bin/main --type radius --mode unif
