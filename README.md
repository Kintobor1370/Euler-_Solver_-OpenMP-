# Euler-_Solver_-OpenMP-
Solver for integral-differential equations for the first and second spatial moments using the Euler method.\
These programs carry out experiments 14.4 and 14.6 (LawDieckmann200a.pdf: pages 13, 17)
\
\
_Tested on Ubuntu 22.04 and SteamOS 3.4.6._
\
\
HOW TO EXECUTE THE CODE:

1. Before compiling, navigate to the main body of Solver.cpp and uncomment the experiment of your choice.
2. Launch the command prompt from the folder containing the programs.
3. Input the following line:
```
g++ Solver.cpp -fopenmp -O3
```
```
./a.out
```
This will start the experiment of your choice.
