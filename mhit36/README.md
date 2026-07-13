
~~~text
███    ███ ██   ██ ██ ████████ ██████   ██████  
████  ████ ██   ██ ██    ██         ██ ██       
██ ████ ██ ███████ ██    ██     █████  ███████  
██  ██  ██ ██   ██ ██    ██         ██ ██    ██ 
██      ██ ██   ██ ██    ██    ██████   ██████         
~~~


#### GPU-based Finite difference code for DNS of Multiphase Homogenous isotropic turbulence (phase-field method)

Main developer: A. Roccon 

Future developments:
* Conservative Allen-Cahn equation for phase-field modeling (CDI)

Performance (only NS)
* 64  x  64 x  64 | RTX5000 |   5 ms/timestep
* 128 x 128 x 128 | RTX5000 |  30 ms/timestep
* 256 x 256 x 256 | RTX5000 | 240 ms/timestep

#### Systems supported:
* Unix + nvfortran 

### To compile the code
* Use compile.sh in the src/folder

#### To run a simulation (on Leonardo)
* go to src folder and submit ./go.sh

#### Parallelization strategy
* The code is serial and exploit a single GPU.

#### Output and restart files.
* Files containing the Eulerian fields (u\_\*\*\*,v\_\*\*\*\*,w\_\*\*\*\*  etc.) and the particle positions (p\_\*\*\*) are stored in set_run/output

