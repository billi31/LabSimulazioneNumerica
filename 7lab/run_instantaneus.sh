sh prepare/prepare_solid.sh
sh prepare/prepare_start.sh
rm -rf output/solid/instantaneus.*
./bin/main --state solid --insta
sh prepare/prepare_liquid.sh
sh prepare/prepare_repeat.sh
rm -rf output/liquid/instantaneus.*
./bin/main --state liquid --insta
sh prepare/prepare_gas.sh
sh prepare/prepare_repeat.sh
rm -rf output/gas/instantaneus.*
./bin/main --state gas --insta
