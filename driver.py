import sympy as sp
from solver import Solver
from visual import Visual

""" Kinematics Visualisation """

# Program driver
if __name__ == '__main__':

    t = sp.symbols('t')
    tf, dt = 360, 0.5

    # xyz positions of particle in space in terms of time
    x, y, z = sp.sin(3*t), sp.cos(t), sp.sin(2*t)
    # x, y, z = 8*t**2 + 7*t + 6, 5*t**3 + 4, 0.3*t**4 + 2*t**2 + 1
    R = [x, y, z]

    # Database name
    db = 'knmtx.db'
    table = db[:-3]
    print('\nDatabase:', db, '\nTable:', table, '\n')

    # Run configuration
    solve = False
    drop = False
    visualise = True

    if solve:
        sv = Solver(t, R, 0, tf, dt, db, table, memory=False)
        sv.propogate()
        sv.close()

    if visualise:
        vz = Visual(db, table)

    if solve and drop:
        sv.db.drop()
