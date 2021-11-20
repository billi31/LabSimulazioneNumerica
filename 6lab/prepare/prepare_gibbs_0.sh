rm -rf gibbs/output0/out*/output.*
rm -rf gibbs/output0/output.*
rm -rf seed.out
rm -rf config.0
cp config.final config.0
cp gibbs/output0/input.dat input.dat
make all
