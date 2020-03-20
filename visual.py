# Kinematics Visualisation

import sqlite3
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import mpl_toolkits.mplot3d.axes3d as p3
from arrow_3d import Arrow3D
from matplotlib import animation

# Visualise the change over time of kinematic properties


class Visual:

    def __init__(self, table, save=False):
        """
        Used to visualise the simulation results from the
        kinematic propagation of a particle moving through
        space.
        Requires a connection to the results database and
        runs through the simulation.
        """
        print('Visualising...\n')
        plt.rcParams['animation.ffmpeg_path'] = '/usr/local/bin/ffmpeg'
        self.table = table
        self.df = self.query_db()
        self.ani = None
        self.fig = plt.figure(figsize=(11, 8), facecolor='black')
        self.ax = p3.Axes3D(self.fig)
        self.coc, = self.create_line_plot('Centre of Curvature', 'mo')
        self.pos_text = self.new_text(0, 0, 0, ' P', 'g')
        self.coc_text = self.new_text(0, 0, 0, ' C', 'm')
        self.vel_text = self.new_text(0, 0, 0, ' V', 'r')
        self.plot_trajectory()
        self.prev_artists = []
        self.config_plot()
        self.animate()
        if save:
            self.save_animation()
        print('\n... finished')

    def new_text(self, x, y, z, txt, color):
        """
        create new text artist for the plot
        param x: x coordinate of text
        param y: y coordinate of text
        param z: z coordinate of text
        param txt: text string
        param color: text color
        """
        return self.ax.text(x, y, z, txt, size=10, color=color)

    def new_quiver_plot(self, var, color):
        """
        define a new quiver plot for vector
        param var: label text for plot
        param color: color of plot
        """
        # NOTE: method previously used for plotting
        return self.ax.quiver([], [], [], [], [], [],
                              length=1.0,
                              colors=color,
                              normalize=False,
                              label=var)

    @staticmethod
    def create_line_plot(var, style):
        """
        define an empty plot for animating
        param var: label for plot
        param style: line style e.g. 'mo' or 'r-'
        """
        return plt.plot([], [], [], style, label=var)

    @staticmethod
    def vector_artist(x0, x1, y0, y1, z0, z1, color):
        """
        method to create a new arrow in 3d for vectors
        return Arrow3D: new arrow
        """
        return Arrow3D([x0, x1], [y0, y1], [z0, z1],
                       mutation_scale=10,
                       lw=1, arrowstyle='-|>', color=color)

    def query_db(self):
        """
        query the results db and read into pandas dataframe
        """
        conn = sqlite3.connect('./results/%s.db' % self.table)
        return pd.read_sql_query("SELECT * FROM {}".format(self.table), conn)

    def config_plot(self):
        """
        Setting the axes properties such as title, limits, labels
        """
        lim_params = ['r', 'v']
        x_limits = self.get_limits(lim_params, 'x')
        y_limits = self.get_limits(lim_params, 'y')
        z_limits = self.get_limits(lim_params, 'z')
        self.ax.set_xlim3d(x_limits)
        self.ax.set_ylim3d(y_limits)
        self.ax.set_zlim3d(z_limits)
        self.ax.plot([0, 0], [0, 0], [0, 0], 'ko', label='Origin')
        self.ax.plot([x_limits[0], x_limits[1]], [0, 0], [0, 0], 'k-', lw=1)
        self.ax.plot([0, 0], [y_limits[0], y_limits[1]], [0, 0], 'k-', lw=1)
        self.ax.plot([0, 0], [0, 0], [z_limits[0], z_limits[1]], 'k-', lw=1)
        self.new_text(x_limits[1], 0, 0, 'x', 'k')
        self.new_text(0, y_limits[1], 0, 'y', 'k')
        self.new_text(0, 0, z_limits[1], 'z', 'k')
        self.ax.set_xlabel('x')
        self.ax.set_ylabel('y')
        self.ax.set_zlabel('z')
        self.ax.set_title('Kinematics Visualisation')
        self.ax.legend()

    def get_limits(self, params, axis):
        """
        get upper and lower limits for parameter
        param axis: get limits for axis, i.e x or y
        param params: list of varaible names
        """
        lower_lim, upper_lim = 0, 0
        for p in params:
            m = max(self.df['%s%s' % (p, axis)])
            if m > upper_lim:
                upper_lim = m
            m = min(self.df['%s%s' % (p, axis)])
            if m < lower_lim:
                lower_lim = m
        return (lower_lim, upper_lim)

    def plot_trajectory(self):
        """
        plot trajectory of particle in 3d space
        """
        plt.plot(self.df['rx'], self.df['ry'], self.df['rz'], 'b-',
                 label='trajectory')

    def plot_vectors(self, i):
        """
        plot animation function for vector at each time step.
        visualise change in position, velocity, radius of curvature
        param i: time step
        """
        # position vector
        rx, ry, rz = self.get_vector(['rx', 'ry', 'rz'], i)

        # velocity vector
        # v_mag = self.df['v_mag'][i]
        v = self.get_vector(['vx', 'vy', 'vz'], i)
        vx, vy, vz = rx + v[0], ry + v[1], rz + v[2]

        # radius of curvature - (normalising rc vector for plot if /rho)
        # rho = self.df['rho'][i]
        rcx, rcy, rcz = self.get_vector(['rcx', 'rcy', 'rcz'], i)
        self.coc.set_data([[0, rcx], [0, rcy]])
        self.coc.set_3d_properties([0, rcz])

        # update label text
        self.pos_text.set_position((rx, ry))
        self.pos_text.set_3d_properties(rz, 'x')
        self.coc_text.set_position((rcx, rcy))
        self.coc_text.set_3d_properties(rcz, 'x')
        self.vel_text.set_position((vx, vy))
        self.vel_text.set_3d_properties(vz, 'x')

        # update vector artists - 3D Arrows for R, V and Rc
        artists = [self.vector_artist(0, rx, 0, ry, 0, rz, 'g'),
                   self.vector_artist(rx, vx, ry, vy, rz, vz, 'r'),
                   self.vector_artist(rcx, rx, rcy, ry, rcz, rz, 'm')]
        [self.ax.add_artist(a) for a in artists]
        if self.prev_artists:
            [artist.remove() for artist in self.prev_artists]
        self.prev_artists = artists

    def get_vector(self, params, i):
        """
        get dataframe values at index for list of parameters
        param params: list of values to retrieve
        param i: index of data
        """
        return [self.df[param][i] for param in params]

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

    def save_animation(self):
        """
        save simulation animation
        """
        print('Creating animation movie...')
        FFwriter = animation.FFMpegWriter(fps=30)
        self.ani.save('./animations/%s.mp4' % self.table, writer=FFwriter)
        print('... movie saved as %s.mp4' % self.table)
