import sympy as sp
from solver import Solver
from visual import Visual

""" Kinematics Visualisation """

# Define symbolic variable t
t = sp.symbols('t')

# Start, end and time step for simulation
t0, tf, dt = 0, 360 * 6, 2

# [x, y, z] positions of particle in space in terms of time
R1 = [sp.sin(3*t), sp.cos(t), sp.sin(2*t)]  # Sample orbit
R = [sp.cos(t), sp.sin(t), t]  # Spiral: z=t, Circle: z=0*t
print('\nR =', R, '\n')


#######################################################


print('*' * 50, '\n\nRun configuration...')

# NOTE: table name == db name
db_dir, table, memory = './results/', 'spiral_2', False
solve = False
drop = False
visualise = True
save = False

print('\nDatabase: %s.db\tTable: %s\n' % (table, table))

if solve:
    print('generating new simulation results...')
    print('* if db already exists results will be overwritten *')
else:
    print('using existing database...')

if drop:
    print('simulation results will NOT be saved...')
else:
    print('simulation results will be saved...')

if visualise:
    print('running simulation visualisation...')
elif solve and not visualise:
    print(solve, visualise, 'no visualisation, only solving...')
else:
    print('Doing nothing...')

if save:
    print('will save animation mp4 file...')
else:
    print('will NOT save animation mp4 file...')

print('\n', '*' * 50, '\n')


#######################################################


# Program driver
if __name__ == '__main__':

    if solve:
        sv = Solver(db_dir, table, memory, t, R, t0, tf, dt)
        sv.propogate()
        sv.close()

    if visualise:
        vz = Visual(table, db_dir)

    if drop:
        sv.db.drop()
