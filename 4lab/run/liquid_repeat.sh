cp liquid/config.final config.0
cp liquid/old.final old.0
cp liquid/input.liquid input.dat
mv liquid/insta.* ./
#rm liquid/output.*.0
make all
./bin/main repeat
mv config.final liquid
mv old.final liquid
mv frames/*.xyz liquid/frames
mv *epot* liquid
mv *ekin* liquid
mv *etot* liquid
mv *temp* liquid
mv *pres* liquid
