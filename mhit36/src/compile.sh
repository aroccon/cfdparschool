cd src
module purge
module load nvhpc/24.5
module load cuda/12.2
rm *.mod
rm mhit36
make clean
make
make

rm -rf output
mkdir output
rm *.o

#run the code
#./mhit36


