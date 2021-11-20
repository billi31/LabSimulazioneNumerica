cp gas/config.final config.0
cp gas/old.final old.0
cp gas/input.gas input.dat
mv gas/insta.* ./
#rm gas/output.*.0
make all
./bin/main repeat
mv config.final gas
mv old.final gas
mv frames/*.xyz gas/frames
mv *epot* gas
mv *ekin* gas
mv *etot* gas
mv *temp* gas
mv *pres* gas
