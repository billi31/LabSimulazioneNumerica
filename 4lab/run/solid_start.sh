cp config.fcc config.0
cp solid/input.solid input.dat
rm solid/old.*
rm solid/output.*.0
make all
./bin/main start
mv config.final solid
mv old.final solid
mv frames/*.xyz solid/frames
mv *epot* solid
mv *ekin* solid
mv *etot* solid
mv *temp* solid
mv *pres* solid
