cp config.fcc config.0
cp liquid/input.liquid input.dat
rm liquid/old.*
rm liquid/output.*.0
make all
./bin/main start
mv config.final liquid
mv old.final liquid
mv frames/*.xyz liquid/frames
mv *epot* liquid
mv *ekin* liquid
mv *etot* liquid
mv *temp* liquid
mv *pres* liquid
