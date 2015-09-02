"""

SYNOPSIS

    A simple demo for Rapid trajectory generation for quadrocopters

DESCRIPTION
    
    Generates a single trajectory, and runs input and position feasibility
    tests. Then some plots are generated to visualise the results.

AUTHOR
    
    Mark W. Mueller <mwm@mwm.im>

LICENSE

    Copyright 2014 by Mark W. Mueller <mwm@mwm.im>

    This code is free software: you can redistribute
    it and/or modify it under the terms of the GNU General Public
    License as published by the Free Software Foundation, either
    version 3 of the License, or (at your option) any later version.

    This code is distributed in the hope that it will
    be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the code.  If not, see <http://www.gnu.org/licenses/>.
    
VERSION 

    0.0

"""

from __future__ import print_function, division
import quadrocoptertrajectory as quadtraj

# Define the trajectory starting state:
pos0 = [0, 0, 2] #position
vel0 = [0, 0, 0] #velocity
acc0 = [0, 0, 0] #acceleration

# Define the goal state:
posf = [1, 0, 1]  # position
velf = [0, 0, 1]  # velocity
accf = [0, 9.81, 0]  # acceleration

# Define the duration:
Tf = 1

# Define the input limits:
fmin = 5  #[m/s**2]
fmax = 25 #[m/s**2]
wmax = 20 #[rad/s]
minTimeSec = 0.02 #[s]

# Define how gravity lies:
gravity = [0,0,-9.81]
 
traj = quadtraj.RapidTrajectory(pos0, vel0, acc0, gravity)
traj.set_goal_position(posf)
traj.set_goal_velocity(velf)
traj.set_goal_acceleration(accf)

# Note: if you'd like to leave some states free, there are two options to 
# encode this. As exmample, we will be leaving the velocity in `x` (axis 0)
# free:
#
# Option 1: 
# traj.set_goal_velocity_in_axis(1,velf_y);
# traj.set_goal_velocity_in_axis(2,velf_z);
# 
# Option 2:
# traj.set_goal_velocity([None, velf_y, velf_z])
 
# Run the algorithm, and generate the trajectory.
traj.generate(Tf)

# Test input feasibility
inputsFeasible = traj.check_input_feasibility(fmin, fmax, wmax, minTimeSec)

# Test whether we fly into the floor
floorPoint  = [0,0,0]  # a point on the floor
floorNormal = [0,0,1]  # we want to be in this direction of the point (upwards)
positionFeasible = traj.check_position_feasibility(floorPoint, floorNormal)
 
for i in range(3):
    print("Axis #" , i)
    print("\talpha = " ,traj.get_param_alpha(i), "\tbeta = "  ,traj.get_param_beta(i), "\tgamma = " ,traj.get_param_gamma(i))
print("Total cost = " , traj.get_cost())
print("Input feasibility result: ",    quadtraj.InputFeasibilityResult.to_string(inputsFeasible),   "(", inputsFeasible, ")")
print("Position feasibility result: ", quadtraj.StateFeasibilityResult.to_string(positionFeasible), "(", positionFeasible, ")")

###########################################
# Plot the trajectories, and their inputs #
###########################################

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

numPlotPoints = 100
time = np.linspace(0, Tf, numPlotPoints)
position = np.zeros([numPlotPoints, 3])
velocity = np.zeros([numPlotPoints, 3])
acceleration = np.zeros([numPlotPoints, 3])
thrust = np.zeros([numPlotPoints, 1])
ratesMagn = np.zeros([numPlotPoints,1])

for i in range(numPlotPoints):
    t = time[i]
    position[i, :] = traj.get_position(t)
    velocity[i, :] = traj.get_velocity(t)
    acceleration[i, :] = traj.get_acceleration(t)
    thrust[i] = traj.get_thrust(t)
    ratesMagn[i] = np.linalg.norm(traj.get_body_rates(t))

figStates, axes = plt.subplots(3,1,sharex=True)
gs = gridspec.GridSpec(6, 2)
axPos = plt.subplot(gs[0:2, 0])
axVel = plt.subplot(gs[2:4, 0])
axAcc = plt.subplot(gs[4:6, 0])

for ax,yvals in zip([axPos, axVel, axAcc], [position,velocity,acceleration]):
    cols = ['r','g','b']
    labs = ['x','y','z']
    for i in range(3):
        ax.plot(time,yvals[:,i],cols[i],label=labs[i])

axPos.set_ylabel('Pos [m]')
axVel.set_ylabel('Vel [m/s]')
axAcc.set_ylabel('Acc [m/s^2]')
axAcc.set_xlabel('Time [s]')
axPos.legend()
axPos.set_title('States')

infeasibleAreaColour = [1,0.5,0.5]
axThrust = plt.subplot(gs[0:3, 1])
axOmega  = plt.subplot(gs[3:6, 1])
axThrust.plot(time,thrust,'k', label='command')
axThrust.plot([0,Tf],[fmin,fmin],'r--', label='fmin')
axThrust.fill_between([0,Tf],[fmin,fmin],-1000,facecolor=infeasibleAreaColour, color=infeasibleAreaColour)
axThrust.fill_between([0,Tf],[fmax,fmax], 1000,facecolor=infeasibleAreaColour, color=infeasibleAreaColour)
axThrust.plot([0,Tf],[fmax,fmax],'r-.', label='fmax')

axThrust.set_ylabel('Thrust [m/s^2]')
axThrust.legend()

axOmega.plot(time, ratesMagn,'k',label='command magnitude')
axOmega.plot([0,Tf],[wmax,wmax],'r--', label='wmax')
axOmega.fill_between([0,Tf],[wmax,wmax], 1000,facecolor=infeasibleAreaColour, color=infeasibleAreaColour)
axOmega.set_xlabel('Time [s]')
axOmega.set_ylabel('Body rates [rad/s]')
axOmega.legend()

axThrust.set_title('Inputs')

#make the limits pretty:
axThrust.set_ylim([min(fmin-1,min(thrust)), max(fmax+1,max(thrust))])
axOmega.set_ylim([0, max(wmax+1,max(ratesMagn))])

plt.show()
