//                                Solution of distributed control problem
//
// Compile with: make distributed_control
//
// Sample runs:
//    distributed_control -r 3
//    distributed_control -m ../../data/star.mesh -r 3
//
// Description:  This examples solves the following PDE-constrained
//               optimization problem:
//
//         min J(f) = 1/2 \| u - w \|^2_{L^2} + \alpha/2 \| f \|^2_{L^2} 
//
//         subject to   - \Delta u = f    in \Omega
//                               u = 0    on \partial\Omega
//         and            a <= f(x) <= b
//
//                where w = / 1   if x^2 + y^2 <= 0.5
//                          \ 0   otherwise
//
//

#include "mfem.hpp"
#include <memory>
#include <iostream>
#include <fstream>

using namespace std;
using namespace mfem;

/** The Lagrangian for this problem is
 *    
 *    L(u,f,p) = 1/2 (u - w,u-w) + \alpha/2 (f,f)
 *             - (\nabla u, \nabla p) + (f,p)
 * 
 *      u, p \in H^1_0(\Omega)
 *      f \in L^2(\Omega)
 * 
 *  Note that
 * 
 *    \partial_p L = 0        (1)
 *  
 *  delivers the state equation
 *    
 *    (\nabla u, \nabla v) = (f,v)  for all v in H^1_0(\Omega)
 * 
 *  and
 *  
 *    \partial_u L = 0        (2)
 * 
 *  delivers the adjoint equation
 * 
 *    (\nabla p, \nabla v) = (u-w,v)  for all v in H^1_0(\Omega)
 *    
 *  and at the solutions u and p(u) of (1) and (2), respectively,
 * 
 *  D_f J = D_f L = \partial_u L \partial_f u + \partial_p L \partial_f p
 *                + \partial_f L
 *                = \partial_f L
 *                = (\alpha f + p, \cdot)
 * 
 * We update the control f_k with projected gradient descent via
 * 
 *  f_{k+1} = P ( f_k - \gamma R_{L^2}^{-1} D_f J )
 * 
 * where P is the projection operator enforcing a <= u(x) <= b, \gamma is
 * a specified step length and R_{L^2} is the L2-Riesz operator. In other
 * words, we have that
 * 
 * f_{k+1} = max { a, min { b, f_k - \gamma (\alpha f_k + p) } }
 * 
 */

double indicator_function(const Vector & x)
{
   double r = sqrt(x[0]*x[0] + x[1]*x[1]);
   if (r <= 0.5)
   {
      return 1.0;
   }
   else
   {
      return 0.0;
   }
}

double compute_energy(const GridFunction & u, FunctionCoefficient & w_coeff, const GridFunction & f, double alpha)
{
   ConstantCoefficient zero(0.0);
   double energy = f.ComputeL2Error(zero);
   energy *= energy*alpha;
   double diff = u.ComputeL2Error(w_coeff);
   energy += diff * diff;
   return energy/2.0;
}

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   const char *mesh_file = "../../data/inline-quad.mesh";
   int ref_levels = 2;
   int order = 2;
   bool visualization = true;
   double alpha = 1e-4;
   double step_length = 1e0;
   double epsilon = 1.0;
   int max_it = 1e3;
   double tol = 1e-4;
   bool momentum = false;
   double momentum_param = 0.9;

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&ref_levels, "-r", "--refine",
                  "Number of times to refine the mesh uniformly.");
   args.AddOption(&order, "-o", "--order",
                  "Order (degree) of the finite elements.");
   args.AddOption(&step_length, "-sl", "--step-length",
                  "Step length for gradient descent.");
   args.AddOption(&max_it, "-mi", "--max-it",
                  "Maximum number of gradient descent iterations.");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.AddOption(&momentum, "-mom", "--momentum", "-no-mom",
                  "--no-momentum",
                  "Enable gradient descent with momentum.");
   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(cout);
      return 1;
   }
   args.PrintOptions(cout);

   // 2. Read the mesh from the given mesh file. We can handle triangular,
   //    quadrilateral, tetrahedral and hexahedral meshes with the same code.
   Mesh mesh(mesh_file, 1, 1);
   int dim = mesh.Dimension();

   // 3. Define the target function w.
   FunctionCoefficient w_coeff(indicator_function);
   ConstantCoefficient negative_one(-1.0);
   ProductCoefficient negative_w_coeff(w_coeff, negative_one);

   // 4. Refine the mesh to increase the resolution. In this example we do
   //    'ref_levels' of uniform refinement, where 'ref_levels' is a
   //    command-line parameter.
   for (int lev = 0; lev < ref_levels; lev++)
   {
      mesh.UniformRefinement();
   }

   // 5. Define the vector finite element spaces representing the state variable u,
   //    adjoint variable p, and the control variable f.
   H1_FECollection state_fec(order, dim);
   L2_FECollection control_fec(order-1, dim);
   FiniteElementSpace state_fes(&mesh, &state_fec);
   FiniteElementSpace control_fes(&mesh, &control_fec);

   int state_size = state_fes.GetTrueVSize();
   int control_size = control_fes.GetTrueVSize();
   cout << "Number of state unknowns: " << state_size << endl;
   cout << "Number of control unknowns: " << control_size << endl;

   // 6. All boundary attributes will be used for essential (Dirichlet) BC.
   MFEM_VERIFY(mesh.bdr_attributes.Size() > 0,
               "Boundary attributes required in the mesh.");
   Array<int> ess_bdr(mesh.bdr_attributes.Max());
   ess_bdr = 1;
   Array<int> ess_tdof_list;
   state_fes.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

   // 7. Set the initial guess for f and the boundary conditions for u and p.
   GridFunction u(&state_fes);
   GridFunction p(&state_fes);
   GridFunction f(&control_fes);
   u = 0.0;
   p = 0.0;
   f = 0.0;

   // 8. Set up the bilinear form a(.,.) for the state and adjoint equation.
   BilinearForm a(&state_fes);
   ConstantCoefficient diffusion_coeff(epsilon);
   ConstantCoefficient zero(0.0);
   a.AddDomainIntegrator(new DiffusionIntegrator(diffusion_coeff));
   a.Assemble();
   OperatorPtr A;
   // Vector B, C, X;

   // 9. Define the gradient function
   GridFunction grad(&control_fes);
   grad = 0.0;

   // 10. Connect to GLVis. Prepare for VisIt output.
   char vishost[] = "localhost";
   int  visport   = 19916;
   socketstream sout_u,sout_p,sout_f;
   if (visualization)
   {
      sout_u.open(vishost, visport);
      sout_p.open(vishost, visport);
      sout_f.open(vishost, visport);
      sout_u.precision(8);
      sout_p.precision(8);
      sout_f.precision(8);
   }

   // State equation system
   GridFunctionCoefficient f_coeff;
   GridFunctionCoefficient u_coeff;
   a.FormSystemMatrix(ess_tdof_list, A);
   GSSmoother M((SparseMatrix&)(*A));
   CGSolver pcg;
   double rtol = 1e-12;
   double atol = 0.;
   int max_cgit = 200;
   int print_level = 0;

   pcg.SetMaxIter(max_cgit);
   pcg.SetRelTol(rtol);
   pcg.SetAbsTol(atol);
   pcg.SetPrintLevel(print_level);
   pcg.SetOperator(*A);
   pcg.SetPreconditioner(M);

   // 11. Perform projected gradient descent
   for (int k = 1; k < max_it; k++)
   {
      // A. Form state equation
      LinearForm b(&state_fes);
      f_coeff.SetGridFunction(&f);
      b.AddDomainIntegrator(new DomainLFIntegrator(f_coeff));
      b.Assemble();
      a.EliminateVDofsInRHS(ess_tdof_list,u,b);
      
      // B. Solve state equation
      pcg.Mult(b,u);


      // D. Send the solution by socket to a GLVis server.
      // if (visualization && (k % int(max_it/5) == 0) )
      if (visualization)
      {
         sout_u << "solution\n" << mesh << u
                << "window_title 'State u'" << flush;
      }

      // E. Form adjoint equation
      LinearForm c(&state_fes);
      u_coeff.SetGridFunction(&u);
      c.AddDomainIntegrator(new DomainLFIntegrator(u_coeff));
      c.AddDomainIntegrator(new DomainLFIntegrator(negative_w_coeff));
      c.Assemble();
      a.EliminateVDofsInRHS(ess_tdof_list,p,c);
      
      // F. Solve adjoint equation
      pcg.Mult(c,p);

      if (visualization)
      {
         sout_p << "solution\n" << mesh << p 
                << "window_title 'Adjoint p'" << flush;
      }

      // H. Constuct gradient function (i.e., \alpha f + p)
      GridFunction p_L2(&control_fes);
      // p_L2.ProjectGridFunction(p);

      GridFunctionCoefficient pcoeff(&p);
      
      LinearForm d(&control_fes);
      d.AddDomainIntegrator(new DomainLFIntegrator(pcoeff));
      d.Assemble();
      BilinearForm L2proj(&control_fes);
      InverseIntegrator * m = new InverseIntegrator(new MassIntegrator());
      L2proj.AddDomainIntegrator(m);
      L2proj.Assemble();
      Array<int> empty_list;
      OperatorPtr invM;
      L2proj.FormSystemMatrix(empty_list,invM);   
      invM->Mult(d,p_L2);




      if (momentum)
      {
         grad *= momentum_param;
      }
      else
      {
         grad = 0.0;
      }
      grad += p_L2;
      // grad.ProjectGridFunction(p);
      grad /= alpha;
      grad += f;
      grad *= alpha;

      // I. Compute norm of gradient.
      double norm = grad.ComputeL2Error(zero);
      double energy = compute_energy(u, w_coeff, f, alpha);

      // J. Update control.
      grad *= step_length;
      f -= grad;


      // "Project" to [a,b]
      double cmin = -1., cmax = 1.;
      for (int i = 0; i < f.Size(); i++)
      {
         if (f[i] > cmax) 
         {
            f[i] = cmax;
         }
         else if (f[i] < cmin)
         {
            f[i] = cmin;
         }
         else
         { // do nothing
         }
      }

      if (visualization)
      {
         sout_f << "solution\n" << mesh << f 
                << "window_title 'Control f'" << flush;
      }

      // K. Exit if norm of grad is small enough.
      mfem::out << "norm of gradient = " << norm << endl;
      mfem::out << "energy = " << energy << endl;
      if (norm < tol)
      {
         break;
      }

   }

   return 0;
}