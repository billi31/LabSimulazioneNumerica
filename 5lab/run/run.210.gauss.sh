#!/bin/bash
rm output/*r210gauss*
make all
./bin/main --type radius --state 210 --mode gauss
