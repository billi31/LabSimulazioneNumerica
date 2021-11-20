#!/bin/bash
rm output/*r100gauss*
make all
./bin/main --type radius --state 100 --mode gauss
