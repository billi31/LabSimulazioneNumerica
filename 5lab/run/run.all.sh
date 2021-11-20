#!/bin/bash
rm output/*.r*
make all
./bin/main --type radius
