import os
import sympy as sp
from solver import Solver
from visual import Visual

""" Kinematics Visualisation """

# Define symbolic variable t
t = sp.symbols('t')

# Start, end and time step for simulation
t0, tf, dt = 0, 360, 1

# [x, y, z] positions of particle in space in terms of time
R = [sp.sin(3*t), sp.cos(t), sp.cos(2*t)]  # Sample orbit
# R = [sp.cos(t), sp.sin(t), t/5]  # Spiral: z=t, Circle: z=0*t
print('\nR =', R, '\n')


#######################################################


print('*' * 50, '\n\nRun configuration...')

db_dir, db_name, table = './results/', 'Orbit', 'orbit'
anim_dir = './animations/'
solve = True
drop = True
visualise = True
save_anim = True

print('\nDatabase: %s.db\tTable: %s\n' % (db_name, table))

if solve:
    print('generating new simulation results...')
    print('* if db already exists results will be overwritten *')
else:
    print('using existing database...')

if solve and drop:
    print('simulation database will NOT be saved...')
else:
    print('simulation database will be saved...')

if visualise:
    print('running simulation visualisation...')
elif solve and not visualise:
    print('no visualisation, only solving...')
else:
    print('Doing nothing...')

if visualise and save_anim:
    print('will save animation mp4 file...')
else:
    print('will NOT save animation mp4 file...')

print('\n', '*' * 50, '\n')


#######################################################


def make_dir(_dir):
    """
    make directory if doesn't exist
    param _dir: directory string
    """
    if not os.path.exists(_dir):
        os.mkdir(_dir)


# Program driver
if __name__ == '__main__':

    make_dir(db_dir)
    make_dir(anim_dir)

    if solve:
        sv = Solver(db_dir, db_name, table, drop, t, R, t0, tf, dt)
        sv.propogate()
        sv.close()

    if visualise:
        vz = Visual(db_dir, db_name, table, anim_dir, save_anim)

    if solve and drop:
        sv.db.drop()
