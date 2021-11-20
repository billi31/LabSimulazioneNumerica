sh prepare/clean.sh
sh prepare/prepare_solid.sh
sh prepare/prepare_start.sh
./bin/main --state solid --mode start
cp config.final ./output/solid/config.final
cp old.final ./output/solid/old.final
sh prepare/prepare_liquid.sh
sh prepare/prepare_repeat.sh
./bin/main --state liquid --mode repeat
cp config.final ./output/liquid/config.final
cp old.final ./output/liquid/old.final
sh prepare/prepare_gas.sh
sh prepare/prepare_repeat.sh
./bin/main --state gas --mode repeat
cp config.final ./output/liquid/config.final
cp old.final ./output/liquid/old.final
