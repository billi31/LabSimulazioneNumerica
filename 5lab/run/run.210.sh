#!/bin/bash
rm output/*r210*
make all
./bin/main --type radius --state 210
