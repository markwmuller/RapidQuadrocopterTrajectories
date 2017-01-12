Quadrocopter trajectory generator 
================================

Mark W. Mueller (mwm@berkeley.edu)

This contains source code implementing a method for rapidly generating trajectories for quadrocopters. The two folders give implementations for both C++ code, and Python. These are equivalent, but independent of one another.

The Python code allows to quickly plot some trajectories too.

The C++ code has comments which can be compiled with Doxygen.

The algorithm is described in the following paper: 
M.W. Mueller, M. Hehn, and R. D'Andrea, "A computationally efficient motion primitive for quadrocopter trajectory generation," IEEE Transactions on Robotics Volume 31, no.8, pages 1294-1310, 2015.

The paper may be downloaded [here][paperLink].

[paperLink]: http://muellerlab.berkeley.edu/publications/

Quickstart
----------
If you'd just quickly like to see something running, run the file
	`Python/demo.py`
This will make a plot, visualising a trajectory. You can then change the 
trajectory parameters, and see how this affects things.

Licensing
---------
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.  
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
