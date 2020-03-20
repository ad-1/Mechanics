import sympy as sp
from solver import Solver
from visual import Visual

""" Kinematics Visualisation """

# Program driver
if __name__ == '__main__':

    t = sp.symbols('t')
    tf, dt = 360 * 6, 2

    # xyz positions of particle in space in terms of time
    # x, y, z = sp.sin(3*t), sp.cos(t), sp.sin(2*t)  # Random orbit
    # x, y, z = sp.cos(t), sp.sin(t), t  # Spiral
    # x, y, z = sp.cos(t), sp.sin(t), 0*t  # Circle
    x, y, z = sp.cos(t), sp.sin(t), 10*t  # Convergent Spiral
    R = [x, y, z]

    # Database name == table
    table = 'convergent_spiral'
    print('\nDatabase: %s.db\tTable:%s\n', table, table)

    # Run configuration
    solve = True
    drop = False
    visualise = True
    save = False

    if solve:
        sv = Solver(t, R, 0, tf, dt, table, memory=False)
        sv.propogate()
        sv.close()

    if visualise:
        vz = Visual(table, save)

    if solve and drop:
        sv.db.drop()
