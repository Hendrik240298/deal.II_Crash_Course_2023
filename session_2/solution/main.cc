/*
 * Authors: Hendrik Fischer, 2023, Leibniz University Hannover, Germany
 *          Thomas Wick, 2014, ICES, UT Austin, USA; and 2023, Leibniz University Hannover, Germany
 *
 * Features of this program:
 * -------------------------
 *    Computes the nonstationary heat equation in 2d and 3d
 *    A) Constant and/or nonconstant diffusion coefficient
 *    B) Constant and/or nonconstant right hand side f
 *    C) Homogeneous or nonhomogeneous Dirichlet conditions
 *    D) Projection of initial data
 *    E) Declaration of the 2d and 3d classes
 *    F) Possible extension: implement different timestepping schemes
 *       with a suitable choice of theta
 *
 * Nov 2021
 * Running under deal.II version dealii-9.1.1
 *
 * This program is originally based on the deal.II steps 4 and 5
 * in version 8.1.0
 * for implementing the Poisson problem
 * and we refer for Author and License information to those tutorial steps.
 * The time step loop implementation is a straightforward
 * extension; here it is based on the ideas presented in
 * T. Wick; ANS, Vol. 1, No. 1, pp. 1-19.
 *
 * License of this program: GNU Lesser General
 * Public License as published by the Free Software Foundation
 *
 *
 */

////////////////////////////////////////////////
// For explanations of the includes, we refer
// to the deal.II documentation
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

using namespace dealii;

//////////////////////////////////////////////////////
// Main class
template <int dim>
class Step_Heat {
 public:
  Step_Heat();
  void run();

 private:
  void set_runtime_parameters();
  void make_grid();
  void setup_system();
  void assemble_system();
  void solve();
  void output_results(const unsigned int cycle) const;

  void compute_functionals();

  Triangulation<dim> triangulation;
  FE_Q<dim> fe;
  DoFHandler<dim> dof_handler;

  SparsityPattern sparsity_pattern;
  SparseMatrix<double> system_matrix;

  Vector<double> solution, old_timestep_solution;
  Vector<double> system_rhs;

  double timestep, theta, time;
  std::string time_stepping_scheme;
  unsigned int max_no_timesteps, timestep_number;
};

////////////////////////////////////////////////////////
// 1. Initial values
// 2. Boundary values
// 3. Right hand side f
// 4. Function for coefficients

////////////////////////////////////////////////////////
// Initial values

template <int dim>
class InitialValues : public Function<dim> {
 public:
  InitialValues() : Function<dim>(1) {}

  virtual double value(const Point<dim> &p,
                       const unsigned int component = 0) const;

  virtual void vector_value(const Point<dim> &p, Vector<double> &value) const;
};

template <int dim>
double InitialValues<dim>::value(const Point<dim> &p,
                                 const unsigned int component) const {
  return (std::sin(p[0]) * std::sin(p[1]));
}

template <int dim>
void InitialValues<dim>::vector_value(const Point<dim> &p,
                                      Vector<double> &values) const {
  for (unsigned int comp = 0; comp < this->n_components; ++comp)
    values(comp) = InitialValues<dim>::value(p, comp);
}

///////////////////////////////////////////////////////////
// Right hand side
// Info: not fully implemented, but rather directly in the assembly
template <int dim>
class BoundaryValues : public Function<dim> {
 public:
  BoundaryValues() : Function<dim>() {}

  virtual double value(const Point<dim> &p,
                       const unsigned int component = 0) const;
};

template <int dim>
double BoundaryValues<dim>::value(const Point<dim> &p,
                                  const unsigned int /*component*/) const {
  return p.square();
}

///////////////////////////////////////////////////////////
// Right hand side
// Info: not fully implemented, but rather directly in the assembly
template <int dim>
class RightHandSide : public Function<dim> {
 public:
  RightHandSide() : Function<dim>() {}

  virtual double value(const Point<dim> &p,
                       const unsigned int component = 0) const;
};

template <int dim>
double RightHandSide<dim>::value(const Point<dim> &p,
                                 const unsigned int /*component*/) const {
  double return_value = 0;
  for (unsigned int i = 0; i < dim; ++i) return_value += 4 * std::pow(p(i), 4);

  return return_value;
}

///////////////////////////////////////////////////////////
// Non-constant coefficient values for Laplacian
template <int dim>
class Coefficient : public Function<dim> {
 public:
  Coefficient() : Function<dim>() {}

  virtual double value(const Point<dim> &p,
                       const unsigned int component = 0) const;

  virtual void value_list(const std::vector<Point<dim> > &points,
                          std::vector<double> &values,
                          const unsigned int component = 0) const;
};

template <int dim>
double Coefficient<dim>::value(const Point<dim> &p,
                               const unsigned int /*component*/) const {
  if (p[1] >= 0.0)
    return 20;
  else
    return 1;
}

template <int dim>
void Coefficient<dim>::value_list(const std::vector<Point<dim> > &points,
                                  std::vector<double> &values,
                                  const unsigned int component) const {
  Assert(values.size() == points.size(),
         ExcDimensionMismatch(values.size(), points.size()));

  Assert(component == 0, ExcIndexRange(component, 0, 1));

  const unsigned int n_points = points.size();

  for (unsigned int i = 0; i < n_points; ++i) {
    if (points[i][1] >= 0.0)
      values[i] = 20;
    else
      values[i] = 1;
  }
}

///////////////////////////////////////////////////////////
// Constructor
template <int dim>
Step_Heat<dim>::Step_Heat() : fe(1), dof_handler(triangulation) {}

///////////////////////////////////////////////////////////
// Definiting several model parameters
template <int dim>
void Step_Heat<dim>::set_runtime_parameters() {
  // Timestepping schemes
  // BE, FE, CN
  time_stepping_scheme = "CN";

  // BE = 1.0,
  // FE = 0.0,
  // CN = 0.5,
  theta = 0.5;

  // Timestep size:
  timestep = 0.125;

  // Maximum number of timesteps:
  max_no_timesteps = 200;

  // A variable to count the number of time steps
  timestep_number = 0;

  // Counts total time
  time = 0;
}

///////////////////////////////////////////////////////////
// Initialize the mesh
template <int dim>
void Step_Heat<dim>::make_grid() {
  // Mesh generation
  GridGenerator::hyper_cube(triangulation, 0, M_PI);

  triangulation.refine_global(6);

  std::cout << "   Number of active cells: " << triangulation.n_active_cells()
            << std::endl
            << "   Total number of cells: " << triangulation.n_cells()
            << std::endl;
}

///////////////////////////////////////////////////////////
// Setting up sparsity pattern for system matrix, rhs, solution vector
template <int dim>
void Step_Heat<dim>::setup_system() {
  dof_handler.distribute_dofs(fe);

  std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs()
            << std::endl;

  // CompressedSparsityPattern c_sparsity(dof_handler.n_dofs());
  DynamicSparsityPattern c_sparsity(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, c_sparsity);
  sparsity_pattern.copy_from(c_sparsity);

  system_matrix.reinit(sparsity_pattern);

  solution.reinit(dof_handler.n_dofs());
  old_timestep_solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());
}

///////////////////////////////////////////////////////////
// This is the heart in which the
// discretized equations are locally (per element) assembled.
// Therein, local integrals are evaluated with Gauss quadrature.
template <int dim>
void Step_Heat<dim>::assemble_system() {
  QGauss<dim> quadrature_formula(2);

  FEValues<dim> fe_values(fe, quadrature_formula,
                          update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);

  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_q_points = quadrature_formula.size();

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double> cell_rhs(dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  // For previous time step solutions
  std::vector<double> old_timestep_solution_values(n_q_points);

  std::vector<std::vector<Tensor<1, dim> > > old_timestep_solution_grads(
      n_q_points, std::vector<Tensor<1, dim> >(1));

  // Right hand side
  const RightHandSide<dim> right_hand_side;

  // Class for nonconstant coefficients
  const Coefficient<dim> coefficient;
  std::vector<double> coefficient_values(n_q_points);

  double density = 1.0;

  // Next, we again have to loop over all elements and assemble local
  // contributions.  Note, that an element is a quadrilateral in two space
  // dimensions, but a hexahedron in 3D.
  typename DoFHandler<dim>::active_cell_iterator cell =
                                                     dof_handler.begin_active(),
                                                 endc = dof_handler.end();

  for (; cell != endc; ++cell) {
    fe_values.reinit(cell);
    cell_matrix = 0;
    cell_rhs = 0;

    // Old_timestep_solution values are interpolated on the present cell
    fe_values.get_function_values(old_timestep_solution,
                                  old_timestep_solution_values);
    fe_values.get_function_gradients(old_timestep_solution,
                                     old_timestep_solution_grads);

    // Update of nonconstant coefficients
    coefficient.value_list(fe_values.get_quadrature_points(),
                           coefficient_values);

    // Loop over all quadrature points (evaluating integrals numerically)
    for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) {
      // Get FE-solution values for previous time step at present quadrature
      // point
      const double old_timestep_u = old_timestep_solution_values[q_point];

      Tensor<1, dim> old_timestep_grad_u;
      old_timestep_grad_u[0] = old_timestep_solution_grads[q_point][0][0];
      old_timestep_grad_u[1] = old_timestep_solution_grads[q_point][0][1];

      for (unsigned int i = 0; i < dofs_per_cell; ++i) {
        for (unsigned int j = 0; j < dofs_per_cell; ++j) {
          // Mass term from the time derivative
          // Change indices such that (j,i) because test
          // function determines the row in system matrix A
          // See lecture notes in Chapter 6
          cell_matrix(j, i) += density * (fe_values.shape_value(i, q_point) *
                                          fe_values.shape_value(j, q_point) *
                                          fe_values.JxW(q_point));

          cell_matrix(j, i) +=
              timestep * theta *
              (  // A) Constant diffusion a=1:
                  1.0 *
                  // Nonconstant diffusion coefficient
                  // coefficient_values[q_point] *
                  fe_values.shape_grad(i, q_point) *
                  fe_values.shape_grad(j, q_point) * fe_values.JxW(q_point));
        }


        // rhs of PDE
        cell_rhs(i) +=
            (fe_values.shape_value(i, q_point) *
            // B) Constant right hand side f=0:
            0.0
            // Non constant right hand side:
            // right_hand_side.value (fe_values.quadrature_point (q_point))
             ) *
            fe_values.JxW(q_point);

        // contribution of old time step solution to rhs
        cell_rhs(i) += 
            (
              old_timestep_u * fe_values.shape_value(i, q_point)
              - timestep * (1.0 - theta) * 1.0 *
            (old_timestep_grad_u * fe_values.shape_grad(i, q_point))
            ) 
            * fe_values.JxW(q_point);
      }
    }

    // Putting local element contributions into global system matrix A
    // and global right hand side vector b
    cell->get_dof_indices(local_dof_indices);
    for (unsigned int i = 0; i < dofs_per_cell; ++i) {
      for (unsigned int j = 0; j < dofs_per_cell; ++j)
        system_matrix.add(local_dof_indices[i], local_dof_indices[j],
                          cell_matrix(i, j));

      system_rhs(local_dof_indices[i]) += cell_rhs(i);
    }
  }

  // Apply Dirichlet boundary conditions
  std::map<types::global_dof_index, double> boundary_values;
  VectorTools::interpolate_boundary_values(
      dof_handler, 0,
      // C) Homogeneous Dirichlet conditions:
      ZeroFunction<dim>(),
      // Non-homogeneous Dirichlet conditions:
      // BoundaryValues<dim>(),
      boundary_values);
  MatrixTools::apply_boundary_values(boundary_values, system_matrix, solution,
                                     system_rhs);
}

////////////////////////////////////////////////////////
// Linear solution using the conjugate gradient
// method without preconditioning.
// See deal.II tutorials for extensions with preconditioners
// See also Thomas Wick's lecture notes for extensions towards geometric
// multigrid
template <int dim>
void Step_Heat<dim>::solve() {
  SolverControl solver_control(1000, 1e-12);
  SolverCG<> solver(solver_control);
  solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());

  std::cout << "   " << solver_control.last_step()
            << " CG iterations needed to obtain convergence." << std::endl;
}

/////////////////////////////////////////////////////////
// Goal functional evaluations
// -> correctness of solution
// -> possible interest in terms of applications
template <int dim>
void Step_Heat<dim>::compute_functionals() {
  // Compute global L2 norm
  Vector<float> difference_per_cell(triangulation.n_active_cells());
  VectorTools::integrate_difference(dof_handler, solution, ZeroFunction<dim>(),
                                    difference_per_cell, QGauss<dim>(4),
                                    VectorTools::L2_norm);
  const double L2_error = VectorTools::compute_global_error(
      triangulation, difference_per_cell, VectorTools::L2_norm);

  // Compute point value
  Point<dim> p = Point<dim>(M_PI / 2.0, M_PI / 2.0);
  unsigned int component = 0;
  Vector<double> tmp_vector(1);
  VectorTools::point_value(dof_handler, solution, p, tmp_vector);

  // Output
  std::cout << "L2 norm:     " << time << "   " << std::scientific << L2_error
            << std::endl;
  std::cout << "Point value: " << time << "   " << p << "   "
            << tmp_vector(component) << std::endl;
  std::cout << std::endl;
}

/////////////////////////////////////////////////////////
// Graphical output: here *.vtk
template <int dim>
void Step_Heat<dim>::output_results(const unsigned int cycle) const {
  DataOut<dim> data_out;

  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "solution");

  data_out.build_patches();

  std::string filename_basis = "results/solution_";
  std::ostringstream filename;
  std::cout << "------------------" << std::endl;
  std::cout << "Write solution" << std::endl;
  std::cout << "------------------" << std::endl;
  std::cout << std::endl;

  filename << filename_basis << "2d_" << Utilities::int_to_string(cycle, 2)
            << ".vtk";

  std::ofstream output(filename.str().c_str());
  data_out.write_vtk(output);
}

/////////////////////////////////////////////////////////
// Final method to run the program and designing the
// time step loop
template <int dim>
void Step_Heat<dim>::run() {
  std::cout << "Solving problem in " << dim << " space dimensions."
            << std::endl;

  set_runtime_parameters();
  make_grid();
  setup_system();

  // D) Project (possibly nonhomogeneous) initial values
  {
    AffineConstraints<double> constraints;
    constraints.close();

    VectorTools::project(dof_handler, constraints, QGauss<dim>(3),
                         InitialValues<dim>(), solution);

    output_results(0);
  }

  // Timestep loop
  do {
    std::cout << "Timestep " << timestep_number << " (" << time_stepping_scheme
              << ")"
              << ": " << time << " (" << timestep << ")"
              << "\n=============================="
              << "=====================================" << std::endl;

    std::cout << std::endl;

    // Solve for next time step
    old_timestep_solution = solution;
    assemble_system();
    solve();

    // Write solutions
    output_results(timestep_number + 1);

    // Compute goal functionals
    compute_functionals();

    time += timestep;
    ++timestep_number;
  } while (timestep_number <= max_no_timesteps);
}

/////////////////////////////////////////////////////////
// Final method to start the program
int main() {
  deallog.depth_console(0);

  // initialize instance of step-heat class and call run method
  Step_Heat<2> laplace_problem_2d;
  laplace_problem_2d.run();

  return 0;
}
