#!/bin/bash
rm output/*r100*
make all
./bin/main --type radius --state 100
