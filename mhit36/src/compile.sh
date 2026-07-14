cd src
module purge
module load nvhpc/25.11
module load cuda/13.0
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


