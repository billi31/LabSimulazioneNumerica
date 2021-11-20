#!/bin/bash
rm output/*.r*gauss*
make all
./bin/main --type radius --mode gauss
