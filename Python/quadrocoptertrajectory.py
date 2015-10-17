"""

SYNOPSIS

    Rapid trajectory generation for quadrocopters

DESCRIPTION
    
    An implementation of the algorithm described in the paper "Rapid 
    quadrocopter trajectory generation". Generates trajectories from an 
    arbitrary initial state, to a final state described by (any combination 
    of) position, velocity, and acceleration.
    
    Please refer to the paper for more information.

EXAMPLES

    Please see attached `demo.py` scripts for an example on how to use the 
    trajectory generator. `demo.py` will generate a trajectory with given 
    constraints, and return whether it passes feasibility tests. Then, some 
    plots are generated to visualise the resulting trajectory.
    
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

    0.1

"""

import numpy as np


class SingleAxisTrajectory:
    """A trajectory along one axis.
    
    This is used to construct the optimal trajectory in one axis, planning
    in the jerk to achieve position, velocity, and/or acceleration final
    conditions. The trajectory is initialised with a position, velocity and
    acceleration. 
    
    The trajectory is optimal with respect to the integral of jerk squared.
    
    Do not use this in isolation, this useful through the "RapidTrajectory"
    class, which wraps three of these and allows to test input/state 
    feasibility.

    """

    def __init__(self, pos0, vel0, acc0):
        """Initialise the trajectory with starting state."""
        self._p0 = pos0
        self._v0 = vel0
        self._a0 = acc0
        self._pf = 0
        self._vf = 0
        self._af = 0
        self.reset()

    def set_goal_position(self, posf):
        """Define the goal position for a trajectory."""
        self._posGoalDefined = True
        self._pf = posf

    def set_goal_velocity(self, velf):
        """Define the goal velocity for a trajectory."""
        self._velGoalDefined = True
        self._vf = velf

    def set_goal_acceleration(self, accf):
        """Define the goal acceleration for a trajectory."""
        self._accGoalDefined = True
        self._af = accf

    def generate(self, Tf):
        """ Generate a trajectory of duration Tf.

        Generate a trajectory, using the previously defined goal end states 
        (such as position, velocity, and/or acceleration).

        """
        #define starting position:
        delta_a = self._af - self._a0
        delta_v = self._vf - self._v0 - self._a0*Tf
        delta_p = self._pf - self._p0 - self._v0*Tf - 0.5*self._a0*Tf*Tf

        #powers of the end time:
        T2 = Tf*Tf
        T3 = T2*Tf
        T4 = T3*Tf
        T5 = T4*Tf

        #solve the trajectories, depending on what's constrained:
        if self._posGoalDefined and self._velGoalDefined and self._accGoalDefined:
            self._a = ( 60*T2*delta_a - 360*Tf*delta_v + 720* 1*delta_p)/T5
            self._b = (-24*T3*delta_a + 168*T2*delta_v - 360*Tf*delta_p)/T5
            self._g = (  3*T4*delta_a -  24*T3*delta_v +  60*T2*delta_p)/T5
        elif self._posGoalDefined and self._velGoalDefined:
            self._a = (-120*Tf*delta_v + 320*   delta_p)/T5
            self._b = (  72*T2*delta_v - 200*Tf*delta_p)/T5
            self._g = ( -12*T3*delta_v +  40*T2*delta_p)/T5
        elif self._posGoalDefined and self._accGoalDefined:
            self._a = (-15*T2*delta_a + 90*   delta_p)/(2*T5)
            self._b = ( 15*T3*delta_a - 90*Tf*delta_p)/(2*T5)
            self._g = (- 3*T4*delta_a + 30*T2*delta_p)/(2*T5)
        elif self._velGoalDefined and self._accGoalDefined:
            self._a = 0
            self._b = ( 6*Tf*delta_a - 12*   delta_v)/T3
            self._g = (-2*T2*delta_a +  6*Tf*delta_v)/T3
        elif self._posGoalDefined:
            self._a =  20*delta_p/T5
            self._b = -20*delta_p/T4
            self._g =  10*delta_p/T3
        elif self._velGoalDefined:
            self._a = 0
            self._b =-3*delta_v/T3
            self._g = 3*delta_v/T2
        elif self._accGoalDefined:
            self._a = 0
            self._b = 0
            self._g = delta_a/Tf
        else:
            #Nothing to do!
            self._a = self._b = self._g = 0

        #Calculate the cost:
        self._cost =  (self._g**2) + self._b*self._g*Tf + (self._b**2)*T2/3.0 + self._a*self._g*T2/3.0 + self._a*self._b*T3/4.0 + (self._a**2)*T4/20.0
                
    def reset(self):
        """Reset the trajectory parameters."""
        self._cost = float("inf")
        self._accGoalDefined = self._velGoalDefined = self._posGoalDefined = False
        self._accPeakTimes = [None,None]
        pass
    
    def get_jerk(self, t):
        """Return the scalar jerk at time t."""
        return self._g  + self._b*t  + (1.0/2.0)*self._a*t*t
    
    def get_acceleration(self, t):
        """Return the scalar acceleration at time t."""
        return self._a0 + self._g*t  + (1.0/2.0)*self._b*t*t  + (1.0/6.0)*self._a*t*t*t

    def get_velocity(self, t):
        """Return the scalar velocity at time t."""
        return self._v0 + self._a0*t + (1.0/2.0)*self._g*t*t  + (1.0/6.0)*self._b*t*t*t + (1.0/24.0)*self._a*t*t*t*t

    def get_position(self, t):
        """Return the scalar position at time t."""
        return self._p0 + self._v0*t + (1.0/2.0)*self._a0*t*t + (1.0/6.0)*self._g*t*t*t + (1.0/24.0)*self._b*t*t*t*t + (1.0/120.0)*self._a*t*t*t*t*t

    def get_min_max_acc(self, t1, t2):
        """Return the extrema of the acceleration trajectory between t1 and t2."""
        if self._accPeakTimes[0] is None:
            #uninitialised: calculate the roots of the polynomial
            if self._a:
                #solve a quadratic
                det = self._b*self._b - 2*self._g*self._a
                if det<0:
                    #no real roots
                    self._accPeakTimes[0] = 0
                    self._accPeakTimes[1] = 0
                else:
                    self._accPeakTimes[0] = (-self._b + np.sqrt(det))/self._a
                    self._accPeakTimes[1] = (-self._b - np.sqrt(det))/self._a
            else:
                #_g + _b*t == 0:
                if self._b:
                    self._accPeakTimes[0] = -self._g/self._b
                    self._accPeakTimes[1] = 0
                else:
                    self._accPeakTimes[0] = 0
                    self._accPeakTimes[1] = 0

        #Evaluate the acceleration at the boundaries of the period:
        aMinOut = min(self.get_acceleration(t1), self.get_acceleration(t2))
        aMaxOut = max(self.get_acceleration(t1), self.get_acceleration(t2))

        #Evaluate at the maximum/minimum times:
        for i in [0,1]:
            if self._accPeakTimes[i] <= t1: continue
            if self._accPeakTimes[i] >= t2: continue
            
            aMinOut = min(aMinOut, self.get_acceleration(self._accPeakTimes[i]))
            aMaxOut = max(aMaxOut, self.get_acceleration(self._accPeakTimes[i]))
        return (aMinOut, aMaxOut)
 
    def get_max_jerk_squared(self,t1, t2):
        """Return the extrema of the jerk squared trajectory between t1 and t2."""
        jMaxSqr = max(self.get_jerk(t1)**2,self.get_jerk(t2)**2)
        
        if self._a:
            tMax = -self._b/self._a
            if(tMax>t1 and tMax<t2):
                jMaxSqr = max(pow(self.get_jerk(tMax),2),jMaxSqr)

        return jMaxSqr


    def get_param_alpha(self):
        """Return the parameter alpha which defines the trajectory."""
        return self._a

    def get_param_beta (self):
        """Return the parameter beta which defines the trajectory."""
        return self._b

    def get_param_gamma(self):
        """Return the parameter gamma which defines the trajectory."""
        return self._g

    def get_initial_acceleration(self):
        """Return the start acceleration of the trajectory."""
        return self._a0

    def get_initial_velocity(self):
        """Return the start velocity of the trajectory."""
        return self._v0

    def get_initial_position(self):
        """Return the start position of the trajectory."""
        return self._p0

    def get_cost(self):
        """Return the total cost of the trajectory."""
        return self._cost


#enums for feasibility results:
class InputFeasibilityResult:
    """An enumeration of the possible outcomes for the input feasiblity test.

    If the test does not return ``feasible``, it returns the outcome of the 
    first segment that fails. The different outcomes are:
        0: Feasible -- trajectory is feasible with respect to inputs
        1: Indeterminable -- a section's feasibility could not be determined
        2: InfeasibleThrustHigh -- a section failed due to max thrust constraint
        3: InfeasibleThrustLow -- a section failed due to min thrust constraint

    """
    Feasible, Indeterminable, InfeasibleThrustHigh, InfeasibleThrustLow = range(4)
    
    @classmethod
    def to_string(cls,ifr):
        """Return the name of the result."""
        if   ifr==InputFeasibilityResult.Feasible:
            return "Feasible"
        elif ifr==InputFeasibilityResult.Indeterminable:
            return "Indeterminable"
        elif  ifr==InputFeasibilityResult.InfeasibleThrustHigh:
            return "InfeasibleThrustHigh"
        elif ifr==InputFeasibilityResult.InfeasibleThrustLow:
            return "InfeasibleThrustLow"
        return "Unknown"


class StateFeasibilityResult:
    """An enumeration of the possible outcomes for the state feasiblity test.

    The result is either feasible (0), or infeasible (1).
    """
    Feasible, Infeasible = range(2)
    
    @classmethod
    def to_string(cls,ifr):
        """Return the name of the result."""
        if   ifr==StateFeasibilityResult.Feasible:
            return "Feasible"
        elif ifr==StateFeasibilityResult.Infeasible:
            return "Infeasible"
        return "Unknown"
            

class RapidTrajectory:
    """Rapid quadrocopter trajectory generator.

    A quadrocopter state interception trajectory. The trajectory starts at a
    state defined by the vehicle's position, velocity, and acceleration. The
    acceleration can be calculated directly from the quadrocopter's attitude
    and thrust value. The trajectory duration is fixed, and given by the user.

    The trajectory goal state can include any combination of components from
    the quadrocopter's position, velocity, and acceleration. The acceleration
    allows to encode the direction of the quadrocopter's thrust at the end time.

    The trajectories are generated without consideration for any constraints,
    and are optimal with respect to the integral of the jerk squared (which is
    equivalent to an upper bound on a product of the inputs).

    The trajectories can then be tested with respect to input constraints
    (thrust/body rates) with an efficient, recursive algorithm. Whether linear
    combinations of states along the trajectory remain within some bounds can
    also be tested efficiently.

		For more information, please see the publication 'A computationally 
		efficient motion primitive for quadrocopter trajectory generation', 
		avaialable here: http://www.mwm.im/research/publications/

    NOTE: in the publication, axes are 1-indexed, while here they are
    zero-indexed.

    """

    def __init__(self, pos0, vel0, acc0, gravity):
        """Initialise the trajectory.

        Initialise the trajectory with the initial quadrocopter state, and the
        orientation of gravity for this problem.

        The orientation of gravity is required for the feasibility tests.

        Args:
          pos0 (array(3)): Initial position
          vel0 (array(3)): Initial velocity
          acc0 (array(3)): Initial acceleration
          gravity (array(3)): The acceleration due to gravity, in the frame of 
              the trajectories (e.g. [0,0,-9.81] for an East-North-Up frame).

        """

        self._axis = [SingleAxisTrajectory(pos0[i],vel0[i],acc0[i]) for i in range(3)]
        self._grav = gravity
        self._tf = None
        self.reset()

    def set_goal_position(self, pos):
        """ Define the goal end position.

        Define the end position for all three axes. To leave components free, 
        list the end state as ``None`` in the argument, or use the function
        `set_goal_position_in_axis`.

        """
        for i in range(3):
            if pos[i] is None:
                continue
            self.set_goal_position_in_axis(i,pos[i])

    def set_goal_velocity(self, vel):
        """ Define the goal end velocity.

        Define the end velocity for all three axes. To leave components free, 
        list the end state as ``None`` in the argument, or use the function
        `set_goal_velocity_in_axis`.

        """
        for i in range(3):
            if vel[i] is None:
                continue
            self.set_goal_velocity_in_axis(i,vel[i])

    def set_goal_acceleration(self, acc):
        """ Define the goal end acceleration.

        Define the end acceleration for all three axes. To leave components
        free, list the end state as ``None`` in the argument, or use the
        function `set_goal_acceleration_in_axis`.

        """
        for i in range(3):
            if acc[i] is None:
                continue
            self.set_goal_acceleration_in_axis(i,acc[i])

    def set_goal_position_in_axis(self, axNum, pos):
        """ Define the goal end position in axis `axNum`."""
        self._axis[axNum].set_goal_position(pos)

    def set_goal_velocity_in_axis(self, axNum, vel):
        """ Define the goal end velocity in axis `axNum`."""
        self._axis[axNum].set_goal_velocity(vel)

    def set_goal_acceleration_in_axis(self, axNum, acc):
        """ Define the goal end acceleration in axis `axNum`."""
        self._axis[axNum].set_goal_acceleration(acc)

    def reset(self):
        """ Reset the trajectory generator.

        Removes all goal states, and resets the cost. Use this if you want to
        try multiple trajectories from one initial state.

        """
        for i in range(3):
            self._axis[i].reset()

    def generate(self, timeToGo):
        """ Calculate a trajectory of duration `timeToGo`.

        Calculates a trajectory of duration `timeToGo`, with the problem data
        defined so far. If something (e.g. goal position) has not been defined,
        it is assumed to be left free. 

        """
        self._tf = timeToGo
        for i in range(3):
            self._axis[i].generate(self._tf)

    def check_input_feasibility(self, fminAllowed, fmaxAllowed, wmaxAllowed,
                                minTimeSection):
        """ Run recursive input feasibility test on trajectory.

        Attempts to prove/disprove the feasibility of the trajectory with
        respect to input constraints. The result is one of three outcomes:
        (i):   the trajectory is definitely input feasible
        (ii):  the trajectory is definitely input infeasible
        (iii): input feasibility could not be determined

        If the feasibility is indeterminable, this should be treated as
        infeasible.

        Args:
            fminAllowed (float): minimum thrust allowed. [m/s**2]
            fmaxAllowed (float): maximum thrust allowed. [m/s**2]
            wmaxAllowed (float): maximum body rates allowed. [rad/s]
            minTimeSection (float): minimum time interval to be tested during
                the recursion. [s]

        Returns:
            An enumeration, of type InputFeasibilityResult.

        """
        
        return self._check_input_feasibility_section(fminAllowed, fmaxAllowed,
                                    wmaxAllowed, minTimeSection, 0, self._tf)

    def _check_input_feasibility_section(self, fminAllowed, fmaxAllowed, 
                             wmaxAllowed, minTimeSection, t1, t2):
        """Recursive test used by `check_input_feasibility`.

        Returns:
            An enumeration, of type InputFeasibilityResult.

        """

        if (t2-t1)<minTimeSection:
            return InputFeasibilityResult.Indeterminable

        #test the acceleration at the two limits:
        if max(self.get_thrust(t1), self.get_thrust(t2)) > fmaxAllowed:
            return InputFeasibilityResult.InfeasibleThrustHigh
        if min(self.get_thrust(t1), self.get_thrust(t2)) < fminAllowed:
            return InputFeasibilityResult.InfeasibleThrustLow

        fminSqr = 0
        fmaxSqr = 0
        jmaxSqr = 0

        #Test the limits of the box we're putting around the trajectory:
        for i in range(3):
            amin, amax = self._axis[i].get_min_max_acc(t1, t2)

            #distance from zero thrust point in this axis
            v1 = amin - self._grav[i] #left
            v2 = amax - self._grav[i] #right

            #definitely infeasible:
            if (max(v1**2, v2**2) > fmaxAllowed**2):
                return InputFeasibilityResult.InfeasibleThrustHigh

            if(v1*v2 < 0):
                #sign of acceleration changes, so we've gone through zero
                fminSqr += 0
            else:
                fminSqr += min(np.fabs(v1), np.fabs(v2))**2

            fmaxSqr += max(np.fabs(v1), np.fabs(v2))**2

            jmaxSqr += self._axis[i].get_max_jerk_squared(t1, t2)

        fmin = np.sqrt(fminSqr)
        fmax = np.sqrt(fmaxSqr)
        if fminSqr > 1e-6:
            wBound = np.sqrt(jmaxSqr / fminSqr)#the 1e-6 is a divide-by-zero protection
        else:
            wBound = float("inf")

        #definitely infeasible:
        if fmax < fminAllowed:
            return InputFeasibilityResult.InfeasibleThrustLow
        if fmin > fmaxAllowed:
            return InputFeasibilityResult.InfeasibleThrustHigh

        #possibly infeasible:
        if (fmin < fminAllowed) or (fmax > fmaxAllowed) or (wBound > wmaxAllowed):
            #indeterminate: must check more closely:
            tHalf = (t1 + t2) / 2.0
            r1 = self._check_input_feasibility_section(fminAllowed, fmaxAllowed, wmaxAllowed, minTimeSection, t1, tHalf)
            
            if r1 == InputFeasibilityResult.Feasible:
                #check the other half
                return self._check_input_feasibility_section(fminAllowed, fmaxAllowed, wmaxAllowed, minTimeSection, tHalf, t2)
            else:
                #infeasible, or indeterminable
                return r1

        #definitely feasible:
        return InputFeasibilityResult.Feasible

    def check_position_feasibility(self, boundaryPoint, boundaryNormal):
        """Test whether the position trajectory is allowable w.r.t. a plane.

        Test whether the position trajectory remains on the allowable side
        of a given plane. The plane is defined by giving a point on the plane,
        and the normal vector to the plane.
        
        The result is of the class StateFeasibilityResult, either Feasible or
        Infeasible.

        Args:
            boundaryPoint (array(3)): a point lying on the plane defining the 
                boundary.
            boundaryNormal (array(3)): a vector defining the normal of the 
                boundary. All points lying in the direction of the normal from
                the boundary are taken as feasible.

        Returns:
            An enumeration, of type StateFeasibilityResult.

        """

        boundaryNormal = np.array(boundaryNormal)
        boundaryPoint  = np.array(boundaryPoint)
        
        #make sure it's a unit vector:
        boundaryNormal = boundaryNormal/np.linalg.norm(boundaryNormal)

        #first, we will build the polynomial describing the velocity of the a
        #quadrocopter in the direction of the normal. Then we will solve for 
        #the zeros of this, which give us the times when the position is at a
        #critical point. Then we evaluate the position at these points, and at
        #the trajectory beginning and end, to see whether we are feasible. 
        
        coeffs = np.zeros(5)

        for i in range(3):
            coeffs[0] += boundaryNormal[i]*self._axis[i].get_param_alpha()/24.0            # t**4
            coeffs[1] += boundaryNormal[i]*self._axis[i].get_param_beta() /6.0              # t**3
            coeffs[2] += boundaryNormal[i]*self._axis[i].get_param_gamma()/2.0             # t**2
            coeffs[3] += boundaryNormal[i]*self._axis[i].get_initial_acceleration()/6.0    # t
            coeffs[4] += boundaryNormal[i]*self._axis[i].get_initial_velocity()            # 1
        
        #calculate the roots
        tRoots = np.roots(coeffs)
        
        #test these times, and the initial & end times:
        for t in np.append(tRoots,[0,self._tf]):
            distToPoint = np.dot(self.get_position(t) - boundaryPoint, boundaryNormal)
            if distToPoint <= 0:
                return StateFeasibilityResult.Infeasible
        
        #all points tested feasible:
        return StateFeasibilityResult.Feasible

    def get_jerk(self, t):
        """ Return the trajectory's 3D jerk value at time `t`."""
        return np.array([self._axis[i].get_jerk(t) for i in range(3)])

    def get_acceleration(self, t):
        """ Return the trajectory's 3D acceleration value at time `t`."""
        return np.array([self._axis[i].get_acceleration(t) for i in range(3)])

    def get_velocity(self, t):
        """ Return the trajectory's 3D velocity value at time `t`."""
        return np.array([self._axis[i].get_velocity(t) for i in range(3)])

    def get_position(self, t):
        ''' Return the trajectory's 3D position value at time `t`.'''
        return np.array([self._axis[i].get_position(t) for i in range(3)])

    def get_normal_vector(self, t):
        """ Return the vehicle's normal vector at time `t`.

        The vehicle's normal vector is that vector along which the thrust
        points, e_3. The required body rates to fly a trajectory can be 
        calculated by finding that angular velocity which rotates this 
        normal vector from one direction to another. Note that the result
        will be expressed in the planning frame, so that a rotation is
        necessary to the body frame.

        Args:
            t (float): time argument.

        Returns:
            np.array() containing a unit vector.

        """
        v = (self.get_acceleration(t) - self._grav)
        return v/np.linalg.norm(v)

    def get_thrust(self, t):
        """ Return the thrust input at time `t`.

        Returns the thrust required at time `t` along the trajectory, in units
        of acceleration. 

        Args:
            t (float): time argument.

        Returns:
            np.array() containing a unit vector.

        """
        return np.linalg.norm(self.get_acceleration(t) - self._grav)

    def get_body_rates(self, t, dt=1e-3):
        """ Return the body rates input at time `t`, in inertial frame.

        Returns the body rates required at time `t` along the trajectory, in 
        units of [rad/s]. This is done by discretizing the normal direction
        trajectory, with discretization `dt`.
        
        **To get (p,q,r) rates, rotate these with the vehicle's attitude.**

        Args:
            t (float): time argument.
            dt (float, optional): discretization time, default is 1ms

        Returns:
            np.array() containing the rates, in the inertial frame.

        """
        n0 = self.get_normal_vector(t)
        n1 = self.get_normal_vector(t + dt)

        crossProd = np.cross(n0,n1) #direction of omega, in inertial axes

        if np.linalg.norm(crossProd) > 1e-6:
            return  np.arccos(np.dot(n0,n1))/dt*(crossProd/np.linalg.norm(crossProd))
        else:
            return np.array([0,0,0])


    def get_cost(self):
        """ Return the total trajectory cost.

        Returns the total trajectory cost. Trajectories with higher cost will 
        tend to have more aggressive inputs (thrust and body rates), so that 
        this is a cheap way to compare two trajectories.

        """
        return self._axis[0].get_cost() + self._axis[1].get_cost() + self._axis[2].get_cost() 

    def get_param_alpha(self, axNum):
        """Return the three parameters alpha which defines the trajectory."""
        return self._axis[axNum].get_param_alpha()

    def get_param_beta(self, axNum):
        """Return the three parameters beta which defines the trajectory."""
        return self._axis[axNum].get_param_beta()

    def get_param_gamma(self, axNum):
        """Return the three parameters gamma which defines the trajectory."""
        return self._axis[axNum].get_param_gamma()

