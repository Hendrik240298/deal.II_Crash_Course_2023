/* ---------------------------------------------------------------------
 *
 * Copyright (C) 1999 - 2023 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------
 
 *
 * Authors: Wolfgang Bangerth, 1999,
 *          Guido Kanschat, 2011
 */
 
 
 
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
 
#include <deal.II/fe/fe_q.h>
 
#include <deal.II/dofs/dof_tools.h>
 
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
 
#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
 
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
 
#include <deal.II/numerics/data_out.h>
#include <fstream>
#include <iostream>
 
using namespace dealii;
 
 
 
class Poisson
{
public:
  Poisson();
 
  void run();
 
 
private:
  void make_grid();
  void setup_system();
  void assemble_system();
  void solve();
  void output_results() const;
 
  Triangulation<2> triangulation;
  FE_Q<2>          fe;
  DoFHandler<2>    dof_handler;
 
  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;
 
  Vector<double> solution;
  Vector<double> system_rhs;
};
 
 
Poisson::Poisson()
  : fe(1)
  , dof_handler(triangulation)
{}
 
 
 
void Poisson::make_grid()
{
  /* MISSING CODE: Use GridGenerator to generate quadratic mesh*/
  
  triangulation.refine_global(5);
 
  std::cout << "Number of active cells: " << triangulation.n_active_cells()
            << std::endl;
}
 
 
 
 
void Poisson::setup_system()
{
  dof_handler.distribute_dofs(fe);
  std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs()
            << std::endl;
 
  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp);
  sparsity_pattern.copy_from(dsp);
 
  system_matrix.reinit(sparsity_pattern);
 
  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());
}
 
 
 
void Poisson::assemble_system()
{
  QGauss<2> quadrature_formula(fe.degree + 1);

  FEValues<2> fe_values(fe,
                        quadrature_formula,
                        /* MISSING CODE: Which information is needed for the assembly? */);
 
  const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
 
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);
 
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      fe_values.reinit(cell);
 
      /* MISSING CODE:
      Reset the local cell's contributions (cell_matrix and cell_rhs)
      to global matrix (cell_matrix) and global rhs
      */
 
      for (const unsigned int q_index : fe_values.quadrature_point_indices())
        {
          // system matrix
          for (const unsigned int i : fe_values.dof_indices())
            for (const unsigned int j : fe_values.dof_indices())
                /* MISSING CODE:
                Assemble: grad phi_i(x_q) * grad phi_j(x_q) * dx
                */
 
          // system rhs
          for (const unsigned int i : fe_values.dof_indices())
              /* MISSING CODE:
              Assemble: phi_i(x_q) * f(x_q) * dx
              */            
        }

      // transfer the local elements to the global matrix.
      cell->get_dof_indices(local_dof_indices);
      for (const unsigned int i : fe_values.dof_indices())
        for (const unsigned int j : fe_values.dof_indices())
          system_matrix.add(/* MISSING CODE: ADD local to global*/);
 
      for (const unsigned int i : fe_values.dof_indices())
        /* MISSING CODE: ADD local to global*/
    }
 
 

  /* BONUS 1:
  Try out non-zero BC instead of homogeneous Dirichlet BC.
  */
  std::map<types::global_dof_index, double> boundary_values;
  VectorTools::interpolate_boundary_values(dof_handler,
                                           0,
                                           /* Missing Code: Set zero Dirichlet BC*/,
                                           boundary_values);
  MatrixTools::apply_boundary_values(boundary_values,
                                     system_matrix,
                                     solution,
                                     system_rhs);
}
 
 
 
void Poisson::solve()
{
  SolverControl            solver_control(1000, 1e-6 * system_rhs.l2_norm());
  SolverCG<Vector<double>> solver(solver_control);
  solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());
}
 
 
 
void Poisson::output_results() const
{
  DataOut<2> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "solution");
  data_out.build_patches();
 
  std::ofstream output("solution.vtk");
  data_out.write_vtk(output);
}
 
 
 
void Poisson::run()
{
  /* MISSING CODE:
  The run method is used to call the other methods of the class in order to:
  - create the grid
  - setup the system
  - assemble the system
  - solve the system
  - postprocessing and plotting
  */

 /* Bonus 2:
 Try to solve the problem for different mesh sizes by 
 successively refining the mesh in a loop
 */
}
 
 
 
int main()
{
  deallog.depth_console(2);
 
  Poisson laplace_problem;
  laplace_problem.run();
 
  return 0;
}
