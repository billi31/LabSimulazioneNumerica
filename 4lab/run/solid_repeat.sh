cp solid/config.final config.0
cp solid/old.final old.0
cp solid/input.solid input.dat
mv solid/insta.* ./
#rm solid/output.*.0
make all
./bin/main repeat
mv config.final solid
mv old.final solid
mv frames/*.xyz solid/frames
mv *epot* solid
mv *ekin* solid
mv *etot* solid
mv *temp* solid
mv *pres* solid
