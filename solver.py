# Mechanics - Kinematics Visualisation

import numpy as np
import sympy as sp
from math import acos
from db import Database

"""
    Given x, y and z Cartesian coordinates (in meters) of a particle P
    as a function of time (in seconds),

    Compute and visualise the change over time of various kinematic
    properties such as:
        - velocity vector V
        - unit tangent vector ut
        - direction angles (theta_x, theta_y, theta_z)
        - acceleration A
        - unit binormal vector ub
        - angles of phi_x and phi_y that A makes with x, y and z
        - radius of curvature of the path (rho)
        - coordinates of rho for a point in the trajectory
"""


class Solver:

    def __init__(self, t, R, t0, tf, dt, table, memory):
        """
        initialise kinematics solver
        param t: symbolic variable t defined using sympi
        param R: symbolic equation of change in position over time
        R = rxi + ryj + rzk
        param t0: initial time
        param tf: final time
        param dt: change in time
        param table: database name and table name
        param memory: (bool) indicating to use RAM or file for db
        """

        if not self.validate_time(t0, tf, dt):
            print('Invalid time input...')
            raise ValueError('Time must be positive and numeric')

        self.t = t  # symbolic variable t
        columns = self.get_columns()
        print('variables to write to db:', columns)
        self.db = Database(table, memory, columns)

        print('\n', '=' * 50, '\n\nInitialising solver...')

        self.times = np.arange(t0, tf, dt, dtype='float')
        self.n_steps = (tf - t0) / dt

        print('\nt0 = {}\ntf = {}\ndt = {}\nsteps = {}\n'
              .format(t0, tf, dt, self.n_steps))

        self.R = R
        self.V = self.vector_derivative(R, t, 'position', 'velocity')
        self.A = self.vector_derivative(self.V, t, 'velocity', 'acceleration')
        self.at = self.vector_magnitude(self.V).diff(t)

        self.system_equations_info()
        print('=' * 50)

    @staticmethod
    def validate_time(t0, tf, dt):
        """
        validate time inputs should be numeric and > 0
        param t0: initial time (s)
        param tf: final time (s)
        param dt: change in time (dt)
        """
        try:
            t0 = float(t0)
            tf = float(tf)
            dt = float(dt)
        except Exception:
            return False
        if t0 >= 0 and tf > 0 and dt > 0:
            return True
        return False

    def system_equations_info(self):
        """
        print system information to console
        """
        print('\nPosition, R = {}\n'.format(self.R))
        print('     rx =', self.R[0])
        print('     ry =', self.R[1])
        print('     rz =', self.R[2])
        print('\nVelocity, V = {}\n'.format(self.V))
        print('     vx =', self.V[0])
        print('     vy =', self.V[1])
        print('     vz =', self.V[2], '\n')
        print('\nAccekeration, A = {}\n'.format(self.A))
        print('     ax =', self.A[0])
        print('     ay =', self.A[1])
        print('     az =', self.A[2], '\n')

    def evaluate_state(self, ti):
        """
        evaluate positon, velocity and acceleration at time t
        param ti: current time step
        """
        ti_r = self.d2r(ti)
        r_mag, R = self.eval_v(self.R, ti_r)
        v_mag, V = self.eval_v(self.V, ti_r)
        a_mag, A = self.eval_v(self.A, ti_r)
        # v_theta = [self.r2d(tht) for tht in self.direction_angles(V, v_mag)]
        # a_theta = [self.r2d(tht) for tht in self.direction_angles(A, a_mag)]
        ut = self.unit_vector(V, v_mag)
        ub = self.unit_binormal(V, A)
        un = self.unit_normal(ub, ut)
        at = self.at.subs(self.t, ti_r)
        an = np.dot(A, un)
        rho = v_mag**2 / an
        Rc = R + (rho * un)
        self.insert_data(ti, R[0], R[1], R[2], r_mag,
                         V[0], V[1], V[2], v_mag,
                         A[0], A[1], A[2], a_mag, an, at,
                         rho, Rc[0], Rc[1], Rc[2])

    @staticmethod
    def get_columns():
        """
        table names for results database
        return: string of variable names
        """
        return """
            \'t\', \'rx\', \'ry\', \'rz\', \'r_mag\',
            \'vx\', \'vy\', \'vz\', \'v_mag\',
            \'ax\', \'ay\', \'az\', \'a_mag\', \'an\', \'at\',
            \'rho\', \'rcx\', \'rcy\', \'rcz\'
            """

    def insert_data(self, *args):
        """
        insert results into database
        param *args: parameters to save
        """
        args = [float(arg) for arg in args]
        self.db.insert(args)

    def propogate(self):
        """
        propogate system through time steps
        """
        print('\npropagating state...\n')
        for ti in self.times:
            self.evaluate_state(ti)
        print('...finished')

    def vector_derivative(self, v, wrt, diff=None, result=None):
        """
        differentiate vector components wrt a symbolic variable
        param v: vector to differentiate
        param wrt: symbolic variable
        param diff: (str) quantity being differentiated
        param result: (str) name of resultant quantity
        return V: velcoity vector and components in x, y, z
        """
        if diff is not None and result is not None:
            print('d({})/d{} = {}'
                  .format(diff, wrt, result))
        return [c.diff(wrt) for c in v]

    def eval_v(self, v, ti):
        """
        evaluate vector components and magnitude @ti
        param v: symbolic vector expression to evaluate @ti
        param ti: time step
        return mag, v_xyz: magnitude of vector and components evaluated @ti
        """
        v_xyz = [c.subs(self.t, ti).evalf() for c in v]
        mag = self.vector_magnitude(v_xyz)
        return mag, (v_xyz)

    def direction_angles(self, v, mag=None):
        """
        compute direction angles a vector makes with + x,y,z axes
        param v: vector with x, y, z components
        param mag: magnitude of vector
        """
        if mag is None:
            mag = self.vector_magnitude(v)
        angles = [acos(c / mag) for c in v]
        return angles

    def unit_vector(self, v, v_mag=None):
        """
        determine unit vector in the direction of v
        param v: vector with components x,y,z
        param v_mag: magnitude of vector
        """
        if v_mag is None:
            v_mag = self.vector_magnitude(v)
        return [c / v_mag for c in v]

    def unit_binormal(self, v=None, a=None, ut=None, un=None):
        """
        orthogonal unit vectors in the direction of velocity (ut)
        and in the direction towards the centre of curvature (un)
        for the osculating plane. The unit binormal vector is the
        vector perpendicular to the osculating plane.
        return ub: unit binormal vector
        """
        if ut is not None and un is not None:
            ub = np.cross(ut, un)
        else:
            b = np.cross(v, a)
            ub = self.unit_vector(b)
        return ub

    def unit_normal(self, ub, ut, Rcp=None, rho=None):
        """
        unit normal is the unit vector in the direction of the
        centre of curvature from position P
        param ub: unit binormal
        param ut: unit tangent
        param Rcp: position vector of the centre of curvature, C
        relative to P
        param rho: radius of curvature (m)
        return un: unit normal
        """
        if Rcp is not None:
            un = self.unit_vector(Rcp, rho)
        else:
            un = np.cross(ub, ut)
        return un

    @staticmethod
    def get_components(v):
        """
        get components from vector, v
        return [vx, vy, vz]: vector in x, y, z components
        """
        return v[0], v[1], v[2]

    @staticmethod
    def vector_magnitude(v):
        """
        compute magnitude of a vector
        param v: vector with components of Cartesian form
        return: magnitude of vector
        """
        return (v[0]**2 + v[1]**2 + v[2]**2)**(1/2)

    @staticmethod
    def d2r(d):
        """
        convert from degrees to radians using sympi pi module
        param d: degrees
        return: radians
        """
        return d * (sp.pi.evalf() / 180)

    @staticmethod
    def r2d(rad):
        """
        convert from radians to degrees
        param rad: radians
        return: degrees
        """
        return rad * (180 / sp.pi.evalf())

    def close(self):
        """
        cleanup code
        """
        self.db.close()
