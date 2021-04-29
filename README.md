# fenics_cahn-hilliard
Extension of the demo provided in fenics-dolfin to solve the cahn-hilliard equations in 3d using FEM for the purpose of shape generation

Requirement to run this code is a recent version of fenics-dolfin.
On debian-like systems, this can be installed as: apt install python3-dolfin

You can compile the routine by:

cmake . && make -j2

After this, the program can be run as, for open top boundary condition:
./cahn-hilliard-3d

and for zero c at the top boundary condition (Dirichlet):

./cahn-hilliard-3d-dirichlet

Any changes to the parameters can be done in open.cpp and dirichlet.cpp
