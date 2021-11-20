cp config.fcc config.0
cp gas/input.gas input.dat
rm gas/old.*
rm gas/output.*.0
make all
./bin/main start
mv config.final gas
mv old.final gas
mv frames/*.xyz gas/frames
mv *epot* gas
mv *ekin* gas
mv *etot* gas
mv *temp* gas
mv *pres* gas
