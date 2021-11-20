sh prepare/prepare_metro_0.sh
./bin/main --mode start
sh prepare/prepare_gibbs_0.sh
./bin/main --mode repeat
sh prepare/prepare_metro_02.sh
./bin/main --mode repeat
sh prepare/prepare_gibbs_02.sh
./bin/main --mode repeat
