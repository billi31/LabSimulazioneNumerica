rm -rf output/output*
rm -rf output/grid*
make all
./bin/main --grid --fixed
