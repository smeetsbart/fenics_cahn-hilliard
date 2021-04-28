//Bart Smeets: Demo for CahnHilliard solving for simulating initial conditions for bacterial biofilms
//Parameters are still a bit all over, but you have to change the following:

//c0: Sets the base level of the initial concentration (aka ~ cell density)
//cr: Initial randomness
//dx,dy,dz: Resolution of the mesh. Size of the thing in x will be 1 and y, z will be adjusted accordingly
//lambda: This is the same as the 'gamma' parameter on the Cahn-Hilliard Wiki


#include <dolfin.h>
#include "CahnHilliard2D.h"
#include "CahnHilliard3D.h"

using namespace dolfin;

// Initial conditions
class InitialConditions : public Expression
{
public:

  InitialConditions(MPI_Comm comm) : Expression(2)
  { dolfin::seed(2 + dolfin::MPI::rank(comm)); }

  //To do: Maybe generate 'values' based on some kind of 'seeding' algorithm (we can make it 'x' dependent!)
  void eval(Array<double>& values, const Array<double>& x) const
  {
    double c0 = 0.3;//Sets the base concentration
    double cr = 0.05;//Sets the initial random perturbation

    double randsymm = 2*(1-dolfin::rand());//Symmetric random number between -1 and 1
    values[0]= c0 + cr*randsymm;
    values[1]= 0.0;
  }

};

class DirichletBoundary : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
//     return x[2] < DOLFIN_EPS or x[2] > 0.25 - DOLFIN_EPS;
    return x[2] < DOLFIN_EPS;
  }
};

// User defined nonlinear problem
class CahnHilliardEquation : public NonlinearProblem
{
  public:

    // Constructor
    CahnHilliardEquation(std::shared_ptr<const Form> F,
                         std::shared_ptr<const Form> J) : _F(F), _J(J) {}

    // User defined residual vector
    void F(GenericVector& b, const GenericVector& x)
    { assemble(b, *_F); }

    // User defined assemble of Jacobian
    void J(GenericMatrix& A, const GenericVector& x)
    { assemble(A, *_J); }

  private:

    // Forms
    std::shared_ptr<const Form> _F;
    std::shared_ptr<const Form> _J;
};


int main(int argc, char* argv[])
{
  dolfin::init(argc, argv);

  //Resolution of the mesh in x,y,z:
  unsigned dx = 24;  //Will be length 1 (ONE)
  unsigned dy = 24;  //Will be length dy/dx
  unsigned dz = 6;   //Will be length dz/dx

  // Mesh
  auto mesh = std::make_shared<Mesh>(
    UnitCubeMesh::create({{dx, dy, dz}}, CellType::Type::tetrahedron));

  //Get a reference to the vertex coordinates of the mesh
  std::vector<double> &x = mesh->coordinates();
  unsigned jv = 0;//Index of the vertices
  //Since this loops over the 'flattened' array, we use modulo to select for x, y and z.
  for ( auto i = x.begin(); i != x.end(); ++i ) {
      double dn = (jv%3==0) ? dx : ( (jv%3==1) ? dy : dz  );
      (*i) *= dn / static_cast<double>(dx);//Rescale the unit cube to desired size based on resolution
      jv++;
  }

  // Create function space and forms, depending on spatial dimension
  // of the mesh
  std::shared_ptr<FunctionSpace> V;
  std::shared_ptr<Form> F, J;
  if (mesh->geometry().dim() == 2)
  {
    V = std::make_shared<CahnHilliard2D::FunctionSpace>(mesh);
    F = std::make_shared<CahnHilliard2D::ResidualForm>(V);
    J = std::make_shared<CahnHilliard2D::JacobianForm>(V, V);
  }
  else if(mesh->geometry().dim() == 3)
  {
    V = std::make_shared<CahnHilliard3D::FunctionSpace>(mesh);
    F = std::make_shared<CahnHilliard3D::ResidualForm>(V);
    J = std::make_shared<CahnHilliard3D::JacobianForm>(V, V);
  }
  else
    error("This demo only supports two or three spatial dimensions.");


  // Create solution Functions (at t_n and t_{n+1})
  auto u0 = std::make_shared<Function>(V);
  auto u = std::make_shared<Function>(V);

  auto boundary = std::make_shared<DirichletBoundary>();
  DirichletBC bc(V, u0, boundary);

  // Set solution to intitial condition
  InitialConditions u_initial(mesh->mpi_comm());
  *u0 = u_initial;
  *u = u_initial;

  // Time stepping and model parameters
  auto dt = std::make_shared<Constant>(5.0e-6);//Simulation timestep
  auto theta = std::make_shared<Constant>(0.5);//Determines integration type. Set to 0.5
  auto lambda = std::make_shared<Constant>(20e-3);//This is Gamma from the wiki of the Cahn-Hilliard equations!

  // Collect coefficient into groups
  std::map<std::string, std::shared_ptr<const GenericFunction>> coefficients
    = {{"u", u}, {"lmbda", lambda}, {"dt", dt}, {"theta", theta}};

  // Add extra coefficient for residual
  std::map<std::string, std::shared_ptr<const GenericFunction>> coefficients_F = coefficients;
  coefficients_F.insert({"u0", u0});

  // Attach coefficients to form
  J->set_coefficients(coefficients);
  F->set_coefficients(coefficients_F);

  double t = 0.0;
  double T = 200*(*dt);

  // Create user-defined nonlinear problem
  CahnHilliardEquation cahn_hilliard(F, J);

  // Create nonlinear solver and set parameters
  NewtonSolver newton_solver;
  newton_solver.parameters["linear_solver"] = "lu";
  newton_solver.parameters["convergence_criterion"] = "incremental";
  newton_solver.parameters["maximum_iterations"] = 10;
  newton_solver.parameters["relative_tolerance"] = 5e-6;
  newton_solver.parameters["absolute_tolerance"] = 1e-13;

  // Save initial condition to file
  File file("cahn_hilliard.pvd", "compressed");
  file << (*u)[0];

  // Solve
  while (t < T)
  {
    // Update for next time step
    t += *dt;
    *u0->vector() = *u->vector();

    // Solve
    newton_solver.solve(cahn_hilliard, *u->vector() );

    // Save function to file
    file << std::pair<const Function*, double>(&((*u)[0]), t);
  }

  return 0;
}
