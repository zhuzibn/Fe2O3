1.main:
    the main function, run it to start simulation.
2.distrib:
    the hexagonal atom distribution of alp-Fe2O3
3.systemgeneration:
    initialize the magnetic moment directions
4.integrate_llg；
    calculate effective fields, call solver and save data
5.exchangej_:
    calculate exchange effective field
6.dmi:
    calculate dmi effective field
7.dipole_:
    calculate dipole effective field, usually we don't active dipole cause it costs too much time
8.constantfile:
    const parameters
9.atomgpu(/pc4/rk4):
    differential equation solver
10.One_data_curve；
    use to visualize the result, it won't be called automatically by main 
11.Relax_(~):
    Relaxation results with different angles(with respect to +x axis) of layer1's magnetic moment
