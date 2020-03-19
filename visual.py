# Kinematics Visualisation

import sqlite3
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import mpl_toolkits.mplot3d.axes3d as p3

"""
Visualise the change over time of various
kinematic properties such as.
"""


class Visual:

    def __init__(self, db, table):
        """
        Used to visualise the simulation results from the
        kinematic propagation of a particle moving through
        space.
        Requires a connection to the results database and
        runs through the simulation.
        """
        self.ani = None
        self.df = self.query_db(db, table)
        print('Starting visualisation...\n')
        self.fig = plt.figure()
        self.ax = p3.Axes3D(self.fig)
        self.plot_trajectory()
        self.vel_vectors, = plt.plot([], [], 'r-', label='velocity')
        self.pos_vectors, = plt.plot([], [], 'g-', label='position')
        self.annotate_plot('r')  # using position limits for plot
        self.animate()

    @staticmethod
    def query_db(db, table):
        """
        query the results db and read into pandas dataframe
        """
        conn = sqlite3.connect(db)
        return pd.read_sql_query("SELECT * FROM {}".format(table), conn)

    def annotate_plot(self, lim_param):
        """
        Setting the axes properties
        param lim_param: limiting parameter
        """
        self.ax.set_xlim3d(self.get_limits('%sx' % lim_param))
        self.ax.set_ylim3d(self.get_limits('%sy' % lim_param))
        self.ax.set_zlim3d(self.get_limits('%sz' % lim_param))
        self.ax.set_xlabel('X')
        self.ax.set_ylabel('Y')
        self.ax.set_zlabel('Z')
        self.ax.set_title('Kinematics Visualisation')
        self.ax.legend()

    def plot_trajectory(self):
        """
        plot trajectory of particle in 3d space
        """
        plt.plot(self.df['rx'],
                 self.df['ry'],
                 self.df['rz'],
                 'b-',
                 label='trajectory')

    def plot_vectors(self, i):
        """
        plot velocity vector at each time step for the particle
        param i: time step
        """
        print(self.df['r_mag'][i], self.df['v_mag'][i])
        rx = self.df['rx'][i]
        ry = self.df['ry'][i]
        rz = self.df['rz'][i]
        vx = rx + self.df['vx'][i]
        vy = ry + self.df['vy'][i]
        vz = rz + self.df['vz'][i]
        self.vel_vectors.set_data((rx, vx), (ry, vy))
        self.vel_vectors.set_3d_properties((rz, vz))
        self.pos_vectors.set_data((0.0, rx), (0.0, ry))
        self.pos_vectors.set_3d_properties((0.0, rz))

    def get_limits(self, param):
        """
        get upper and lower limits for parameter
        param param: array of data
        """
        upper_lim = max(self.df[param])
        lower_lim = min(self.df[param])
        return (lower_lim, upper_lim)

    def animate(self):
        """
        animate drawing velocity vector as particle
        moves along trajectory
        """
        self.ani = FuncAnimation(self.fig,
                                 self.plot_vectors,
                                 frames=len(self.df['t']),
                                 init_func=None,
                                 blit=False,
                                 interval=100)
        plt.show()
