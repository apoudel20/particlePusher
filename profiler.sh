
sed -i "s/GRID_SIZE .*/GRID_SIZE $1/g" particlePusher.cu

sed -i "s/log_profiles\/log-[0-9]*-/log_profiles\/log-$1-/g" Makefile

make reset profile
