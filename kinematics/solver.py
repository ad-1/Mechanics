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

    def __init__(self, db_dir, db_name, table, drop_on_init, t, R, t0, tf, dt):
        """
        initialise kinematics solver
        param db_dir: simulation results directory
        param db_name: simulation database name
        param table: table name
        param drop_on_int: remove existing database if exists
        param t: symbolic variable t - removes req to import sympy
        param R: symbolic equation for change in position over time
        vector in the form R = [rxi, ryj, rzk]
        param t0: initial time
        param tf: final time
        param dt: change in time
        """

        print('\n', '=' * 50, '\n\nInitialising solver...')
        columns = self.get_columns()
        print('\nvariables to write to db:', columns)
        self.db = Database(db_dir, db_name, table, drop_on_init, columns)
        self.t = t

        if not self.validate_time(t0, tf, dt):
            print('Invalid time input...')
            raise ValueError('Time must be positive and numeric')

        self.times = np.arange(t0, tf, dt, dtype='float')
        self.n_steps = (tf - t0) / dt

        print('\nt0 = {}\ntf = {}\ndt = {}\nsteps = {}\n'
              .format(t0, tf, dt, self.n_steps))

        self.R = R
        self.V = self.vector_derivative(R, t, 'position', 'velocity')
        self.A = self.vector_derivative(self.V, t, 'velocity', 'acceleration')
        self.at = self.vector_magnitude(self.V).diff(t)

        self.system_equations_info()
        print('\n', '=' * 50)

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
        print equation information to console
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
        print('     az =', self.A[2])

    def evaluate_state(self, ti):
        """
        evaluate position, velocity, acceleration etc. at time ti
        param ti: current time step
        """
        # NOTE: no. of parms in insert_data must match no. of columns
        ti_r = self.d2r(ti)
        r_mag, R = self.eval_v(self.R, ti_r)
        v_mag, V = self.eval_v(self.V, ti_r)
        a_mag, A = self.eval_v(self.A, ti_r)
        v_theta = [self.r2d(tht) for tht in self.direction_angles(V, v_mag)]
        a_psi = [self.r2d(psi) for psi in self.direction_angles(A, a_mag)]
        ut = self.unit_vector(v_vmag=(V, v_mag))
        ub = self.unit_vector(v_a=(V, A))
        un = self.unit_vector(u1_u2=(ub, ut))
        at = self.at.subs(self.t, ti_r)
        an = np.dot(A, un)
        rho = v_mag**2 / an
        Rc = R + (rho * un)
        rc_mag = self.vector_magnitude(Rc)
        self.insert_data(ti, R[0], R[1], R[2], r_mag,
                         V[0], V[1], V[2], v_mag,
                         A[0], A[1], A[2], a_mag, an, at,
                         rho, Rc[0], Rc[1], Rc[2], rc_mag)

    @staticmethod
    def get_columns():
        """
        table names for results database
        return: string of variable names separated by ,
        """
        # NOTE: these are the column names used in the db
        return """
            \'t\', \'rx\', \'ry\', \'rz\', \'r_mag\',
            \'vx\', \'vy\', \'vz\', \'v_mag\',
            \'ax\', \'ay\', \'az\', \'a_mag\', \'an\', \'at\',
            \'rho\', \'rcx\', \'rcy\', \'rcz\', \'rc_mag\'
            """

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
            print('d({})/d{} = {}'.format(diff, wrt, result))
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
        return [acos(c / mag) for c in v]

    def unit_vector(self, v_vmag=None, v_a=None, u1_u2=None):
        """
        Calculate a unit vector using one of three input parameters.
        Two ways to compute unit vector:

        1. using vector and vector magnitude
        e.g. unit tangent vector is given by v/v_mag
        param v: vector with components x,y,z
        param v_mag: magnitude of vector

        2. using vectors v and a.
        e.g. velocity and acceleration vectors can be used to
        determine the unit binormal. ub = cross(v,a)/mag(cross(v,a))

        3. using orthogonal unit vectors.
        e.g. unit tangent and unit normal vectors can be used to
        determine the unit binormal. ub = corss(ut, un)

        2. and 3. used generally used to calculate one of two vectors.

        a. unit binormal
        Orthogonal unit vectors in the direction of velocity (ut)
        and in the direction towards the centre of curvature (un)
        form the osculating plane. The unit binormal vector is the
        vector perpendicular to the osculating plane.
        param v_a: tuple of vectors velocity vector and acceleration vector
        param ua_ub: tuple of unit tangent and unit normal vectors
        return ub: unit binormal vector

        b. unit normal
        The unit normal is the unit vector in the direction of the
        centre of curvature from the particle at position P.
        param Rcp_rho: tuple position vector of C relative to P and rho
        param ub_ut: tuple of unit binormal and unit tangent vectors
        return un: unit normal vector
        """

        if v_vmag is not None:
            v, mag = v_vmag[0], v_vmag[1]
            if mag is None:
                mag = self.vector_magnitude(v)
            return [c / mag for c in v]

        if v_a is not None:
            v, a = v_a[0], v_a[1]
            b = np.cross(v, a)
            return self.unit_vector(v_vmag=(b, None))

        if u1_u2 is not None:
            u1, u2 = u1_u2[0], u1_u2[1]
            return np.cross(u1, u2)

    @staticmethod
    def vector_magnitude(v):
        """
        compute magnitude of a vector
        param v: vector with components of Cartesian form
        return: magnitude of vector
        """
        # NOTE: np.linalg.norm(v) computes Euclidean norm
        mag = 0
        for c in v:
            mag += c**2
        return mag**(1/2)

    @staticmethod
    def get_components(v):
        """
        get components from vector, v
        return [vx, vy, vz]: vector in x, y, z components
        """
        # TODO: Extend for more dimensions
        return v[0], v[1], v[2]

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
        print('\npropagating state...', end=' ')
        for ti in self.times:
            self.evaluate_state(ti)
        print('finished')

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
        cleanup code. close connection to database
        """
        self.db.close()
