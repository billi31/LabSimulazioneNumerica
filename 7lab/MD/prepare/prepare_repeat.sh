rm -rf src/seed.out
cp config.final config.0
cp old.final old.0
rm -rf *.final
make all
