# Kinematics Visualisation

import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.animation import FuncAnimation
from matplotlib import animation
from arrow_3d import Arrow3D
from db import Database

# Visualise the change over time of kinematic properties


class Visual:

    def __init__(self, db_dir, db_name, table, save_dir, save):
        """
        Used to visualise the simulation results from the
        kinematic propagation of a particle moving through
        space.
        Requires a connection to the results database and
        runs through the simulation.
        param db_dir: simulation results database directory
        param db_name: database name -- filename
        param table: database table name == database filename
        param save_dir: directory to save mp4 animation
        param save: bool indicating to save animation mp4
        """
        print('Visualising...', end=' ')
        self.table = table
        self.db = Database(db_dir, db_name, table)
        self.df = self.db.query()

        plt.rcParams['animation.ffmpeg_path'] = '/usr/local/bin/ffmpeg'
        self.anim_filename = '%s%s.mp4' % (save_dir, self.table)
        self.fig = plt.figure(figsize=(15, 8))
        self.fig.subplots_adjust(left=0.05,
                            bottom=None,
                            right=0.95,
                            top=None,
                            wspace=None,
                            hspace=0.28)
        gs = gridspec.GridSpec(2, 2)

        # axis 1 - 3d visualisation
        self.ax1 = self.fig.add_subplot(gs[:, 0], projection='3d')
        self.traj, = self.ax1.plot([], [], [], 'bo', markersize=1)
        self.coc, = self.ax1.plot([], [], [], 'mo')

        # list of artists which need to be removed and redrawn each frame
        self.vector_lines = []

        # text animations for points moving through 3d space
        self.pos_text = self.text_artist_3d(' P', 'g')
        self.coc_text = self.text_artist_3d(' C', 'm')
        self.vel_text = self.text_artist_3d(' V', 'r')

        # axis 2 - position and velocity vs time
        self.ax2 = self.fig.add_subplot(gs[0, 1])
        self.ax3 = self.ax2.twinx()
        self.rvt, = self.ax2.plot([], [], 'g-')
        self.vvt, = self.ax3.plot([], [], 'r-')

        # axis 3
        self.ax4 = self.fig.add_subplot(gs[1, 1])
        self.ax5 = self.ax4.twinx()
        self.ant, = self.ax4.plot([], [], 'c-')
        self.att, = self.ax5.plot([], [], 'y-')

        self.n_frames = len(self.df['t'])
        self.anim_counter = 0
        self.n_loops = 3 * self.n_frames

        self.anim = self.animate()
        plt.show()
        if save:
            self.save_animation()
        print('finished')

    def text_artist_3d(self, txt, color, x=0, y=0, z=0):
        """
        create new text artist for the plot
        param txt: text string
        param color: text color
        param x: x coordinate of text
        param y: y coordinate of text
        param z: z coordinate of text
        """
        return self.ax1.text(x, y, z, txt, size=12, color=color)

    @staticmethod
    def vector_arrow_3d(x0, x1, y0, y1, z0, z1, color):
        """
        method to create a new arrow in 3d for vectors
        return Arrow3D: new arrow
        """
        return Arrow3D([x0, x1], [y0, y1], [z0, z1],
                       mutation_scale=10, lw=1,
                       arrowstyle='-|>', color=color)

    def query_db(self):
        """
        query the results db and read into pandas dataframe
        """
        if os.path.isfile(self.db):
            self.db.query()
        else:
            print('Simulation results database does not exist...')

    def config_plots(self):
        """
        Setting the axes properties such as title, limits, labels
        """
        self.set_axes_titles()
        self.set_axes_limits()
        self.ax1.set_position([0, 0, 0.5, 1])
        self.ax1.set_aspect('auto')
        self.ax2.grid()
        self.ax4.grid()

    def set_axes_titles(self):
        """
        put titles on plots
        """
        self.ax1.set_title('Trajectory Visualisation')
        self.ax2.set_title('Position and Velocity Magnitudes vs Time')
        self.ax4.set_title('Normal and Tangential Acceleration vs Time')

    def set_axes_limits(self):
        """
        set the axis limits for each plot, label axes
        """
        # NOTE: Using ax1 limits to draw custom xyz axis
        lim_params = ['r', 'v']
        x_lims = self.get_limits(lim_params, 'x')
        y_lims = self.get_limits(lim_params, 'y')
        z_lims = self.get_limits(lim_params, 'z')
        self.ax1.set_xlim3d(x_lims)
        self.ax1.set_ylim3d(y_lims)
        self.ax1.set_zlim3d(z_lims)
        self.draw_xyz_axis(x_lims, y_lims, z_lims)

        t_lims = self.get_limits(['t'], '')

        self.ax2.set_xlim(t_lims)
        self.ax2.set_ylim(self.get_limits(['r_mag'], ''))
        self.ax2.set_xlabel('t (s)')
        self.ax2.set_ylabel('r_mag (m)', color='g')

        self.ax3.set_ylim(self.get_limits(['v_mag'], ''))
        self.ax3.set_ylabel('v_mag (m/s)', color='r')

        self.ax4.set_xlim(t_lims)
        self.ax4.set_ylim(self.get_limits(['an'], ''))
        self.ax4.set_xlabel('t (s)')
        self.ax4.set_ylabel('an (m/s^2)', color='c')

        self.ax5.set_ylim(self.get_limits(['at'], ''))
        self.ax5.set_ylabel('at (m/s^2)', color='y')

    def draw_xyz_axis(self, x_lims, y_lims, z_lims):
        """
        draw xyz axis on ax1 3d plot
        param x_lims: upper and lower x limits
        param y_lims: upper and lower y limits
        param z_lims: upper and lower z limits
        """
        self.ax1.plot([0, 0], [0, 0], [0, 0], 'ko', label='Origin')
        self.ax1.plot(x_lims, [0, 0], [0, 0], 'k-', lw=1)
        self.ax1.plot([0, 0], y_lims, [0, 0], 'k-', lw=1)
        self.ax1.plot([0, 0], [0, 0], z_lims, 'k-', lw=1)
        self.text_artist_3d('x', 'k', x_lims[1], 0, 0)
        self.text_artist_3d('y', 'k', 0, y_lims[1], 0)
        self.text_artist_3d('z', 'k', 0, 0, z_lims[1])
        self.ax1.set_xlabel('x')
        self.ax1.set_ylabel('y')
        self.ax1.set_zlabel('z')

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

    def visualise_data(self, i):
        """
        plot animation function for vector at each time step.
        visualise change in position, velocity, radius of curvature
        param i: time step
        """

        # if self.anim_counter >= self.n_loops:
        #    return

        # ax1

        # position vector
        rx, ry, rz = self.get_vector(['rx', 'ry', 'rz'], i)
        self.traj.set_data(self.df['rx'][:i], self.df['ry'][:i])
        self.traj.set_3d_properties(self.df['rz'][:i])

        # velocity vector
        v = self.get_vector(['vx', 'vy', 'vz'], i)
        vx, vy, vz = rx + v[0], ry + v[1], rz + v[2]

        # radius of curvature - (normalising rc vector)
        rc_mag = self.df['rc_mag'][i]
        rcx, rcy, rcz = self.get_vector(['rcx', 'rcy', 'rcz'], i) / rc_mag
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
        vectors = [self.vector_arrow_3d(0, rx, 0, ry, 0, rz, 'g'),
                   self.vector_arrow_3d(rx, vx, ry, vy, rz, vz, 'r'),
                   self.vector_arrow_3d(rcx, rx, rcy, ry, rcz, rz, 'm')]
        [self.ax1.add_artist(v) for v in vectors]
        if self.vector_lines:
            [vector.remove() for vector in self.vector_lines]
        self.vector_lines = vectors

        # ax2

        t = self.df['t'][:i]

        # magnitude of position vs time
        self.rvt.set_data(t, self.df['r_mag'][:i])
        # magnitude of velocity vs time
        self.vvt.set_data(t, self.df['v_mag'][:i])

        # ax3

        self.ant.set_data(t, self.df['an'][:i])
        self.att.set_data(t, self.df['at'][:i])

        self.anim_counter += 1

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
        return: animation
        """
        return FuncAnimation(self.fig,
                             self.visualise_data,
                             frames=self.n_frames,
                             init_func=self.config_plots(),
                             blit=False,
                             repeat=True,
                             interval=2)

    def save_animation(self):
        """
        save simulation animation
        """
        print('Creating animation movie...')
        FFwriter = animation.FFMpegWriter(fps=30)
        self.anim.save(self.anim_filename, writer=FFwriter)
        print('... movie saved as %s' % self.anim_filename)
