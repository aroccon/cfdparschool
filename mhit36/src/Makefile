# Compiler
FC = pgf90

# Compiler flags
FFLAGS = -fast -gpu=managed  -cuda -cudalib=cufft -acc -Mfree -Mpreprocess -mcmodel=medium
#FFLAGS = -fast -gpu=managed -acc -cudalib -Minfo=accel

# Executable name
EXEC = mhit36

# Source files
SRCS = module.f90  main.f90 poissonfast.f90 readwrite.f90 particles.f90

# Libraries
#LIBS    = -L/usr/local/cuda/lib64 -lcufft
LIBS    =  -lcufft


# Object files
OBJS = $(SRCS:.f90=.o)

# Rule to compile .f90 files to .o files
%.o: %.f90
	$(FC) $(FFLAGS) $(LIBS)  -c $<

# Rule to link object files into executable
$(EXEC): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^

# Clean rule
clean:
	rm -f $(OBJS) $(EXEC)
