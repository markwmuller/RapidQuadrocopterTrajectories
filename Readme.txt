Quadrocopter trajectory generator 

Mark W. Mueller (mwm@ethz.ch)

This contains source code implementing a method for rapidly generating 
trajectories for quadrocopters. The two folders give implementations for both
C++ code, and Python. These are equivalent, but independent of one another.

The Python code allows to quickly plot some trajectories too.

The C++ code has comments which can be compiled with Doxygen.

The algorithm is described in the following paper: 
M.W. Mueller, M. Hehn, and R. D’Andrea, “A computationally efficient motion primitive for quadrocopter trajectory generation,” IEEE Transactions on Robotics (to appear), 2015.
The paper may be downloaded from www.mwm.im/research/publications

Quickstart
====
If you'd just quickly like to see something running, run the file
	Python/demo.py
This will make a plot, visualising a trajectory. You can then change the 
trajectory parameters, and see how this affects things.

