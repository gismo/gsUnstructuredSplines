/** @file gsBiharmonicAssembler2.h

    @brief Provides assembler for a (planar) Biharmonic equation.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/

#pragma once

#include <gsAssembler/gsExprHelper.h>

namespace gismo
{

template <class T>
  class gsBiharmonicAssembler2
  {
  public:

    /// Default empty constructor
    gsBiharmonicAssembler2() { };

    gsBiharmonicAssembler2(gsMultiPatch<T> & mp,
      gsMultiBasis<T> & mb,
      bool & second = false)
    : m_mp(mp), m_mb(mb)
    {
      gsFunctionExpr<>source("1*pi*pi*pi*pi*(4*cos(1*pi*x)*cos(1*pi*y) - cos(1*pi*x) - cos(1*pi*y))",2);
      m_f.swap(source);
      gsInfo << "Source function " << m_f << "\n";

      gsFunctionExpr<> solution("(cos(1*pi*x) - 1) * (cos(1*pi*y) - 1)",2);
      m_ms.swap(solution);
      gsInfo << "Exact function  " << m_ms << "\n";

      gsBoundaryConditions<> bc;	
      for (gsMultiPatch<>::const_biterator bit = mp.bBegin(); bit != mp.bEnd(); ++bit) 
      {
        // Laplace
        gsFunctionExpr<> laplace ("-1*pi*pi*(2*cos(1*pi*x)*cos(1*pi*y) - cos(1*pi*x) - cos(1*pi*y))",2);

        // Neumann
        gsFunctionExpr<> sol1der( "-1*pi*(cos(1*pi*y) - 1)*sin(1*pi*x)",
                                  "-1*pi*(cos(1*pi*x) - 1)*sin(1*pi*y)", 2);

        bc.addCondition(*bit, condition_type::dirichlet, m_ms);
        if (second)
          bc.addCondition(*bit, condition_type::laplace, laplace);
        else
          bc.addCondition(*bit, condition_type::neumann, sol1der);  
      }
      bc.setGeoMap(mp);
      m_bc.swap(bc);
      gsInfo << "Boundary conditions:\n" << m_bc << "\n";
   };

   gsBiharmonicAssembler2(gsMultiPatch<T> & mp,
    gsMultiBasis<T> & mb,
    const gsBoundaryConditions<T> & bc,
    const gsFunctionExpr<T> & f,
    const gsFunctionExpr<T> & ms)
   : m_mp(mp), m_mb(mb), m_bc(bc), m_f(f), m_ms(ms)
   {

   };


 public:

  void initialize()
  {
      // Expression Assembler (1x1 Block Matrix)
    m_A = gsExprAssembler<real_t>(1,1);

      // Elements used for numerical integration
    m_A.setIntegrationElements(m_mb);
    m_ev = gsExprEvaluator<real_t>(m_A);
  };

  void assemble(gsMappedBasis<2,T> & bb2)
  {
      // Set the geometry map
    auto G = m_A.getMap(m_mp);

      // Set the source term
      auto ff = m_A.getCoeff(m_f, G); // Laplace example

      // Set the discretization space
      auto u = m_A.getSpace(bb2);
      
      // Solution vector and solution variable
      gsMatrix<real_t> solVector;
      auto u_sol = m_A.getSolution(u, solVector);

      // Setup the mapper
      gsDofMapper map;
      setMapperForBiharmonic(m_bc, bb2, map);

      // Setup the system
      u.setupMapper(map);
      gsDirichletNeumannValuesL2Projection(m_mp, m_mb, m_bc, bb2, u);

      // Initialize the system
      m_A.initSystem();

      // Compute the system matrix and right-hand side
      m_A.assemble(ilapl(u, G) * ilapl(u, G).tr() * meas(G),u * ff * meas(G));
      
      // Enforce Laplace conditions to right-hand side
      auto g_L = m_A.getBdrFunction(G); // Set the laplace bdy value
      //auto g_L = A.getCoeff(laplace, G);
      m_A.assembleBdr(m_bc.get("Laplace"), (igrad(u, G) * nv(G)) * g_L.tr() );
    };

    void solve()
    {
     gsSparseSolver<real_t>::SimplicialLDLT solver;
     solver.compute( m_A.matrix() );
     m_solVector = solver.solve(m_A.rhs());
   };

   void computeError(gsMappedBasis<2,T> & bb2,
    T & l2err,
    T & h1err,
    T & h2err)
   {	
	// Set the geometry map
     auto G = m_A.getMap(m_mp);

     auto u = m_A.getSpace(bb2);
     auto u_sol = m_A.getSolution(u, m_solVector);

	// Recover manufactured solution
     auto u_ex = m_ev.getVariable(m_ms, G);

     l2err = math::sqrt( m_ev.integral( (u_ex - u_sol).sqNorm() * meas(G) ) );
     h1err = l2err + math::sqrt(m_ev.integral( ( igrad(u_ex) - igrad(u_sol,G) ).sqNorm() * meas(G) ));	
     h2err = h1err + math::sqrt(m_ev.integral( ( ihess(u_ex) - ihess(u_sol,G) ).sqNorm() * meas(G) )); 
   };

   T numDofs() { return m_A.numDofs(); };

 private:

  void setMapperForBiharmonic(gsBoundaryConditions<> & bc, gsMappedBasis<2,real_t> & bb2, gsDofMapper & mapper)
  {
    mapper.setIdentity(bb2.nPatches(), bb2.size(), 1);

    gsMatrix<index_t> bnd;
    for (typename gsBoundaryConditions<real_t>::const_iterator
      it = bc.begin("Dirichlet"); it != bc.end("Dirichlet"); ++it)
    {
      bnd = bb2.basis(it->ps.patch).boundary(it->ps.side());
      mapper.markBoundary(it->ps.patch, bnd, 0);
    }

    for (typename gsBoundaryConditions<real_t>::const_iterator
     it = bc.begin("Neumann"); it != bc.end("Neumann"); ++it)
    {
      bnd = bb2.basis(it->ps.patch).boundaryOffset(it->ps.side(),1);
      mapper.markBoundary(it->ps.patch, bnd, 0);
    }
    mapper.finalize();
  };

  void gsDirichletNeumannValuesL2Projection(gsMultiPatch<> & mp,
   gsMultiBasis<> & dbasis,
   gsBoundaryConditions<> & bc,
   gsMappedBasis<2,real_t> & bb2,
   const expr::gsFeSpace<real_t> & u)
  {
    const gsDofMapper & mapper = u.mapper();

    gsMatrix<index_t> bnd = mapper.findFree(mapper.numPatches()-1);
    gsDofMapper mapperBdy;
      mapperBdy.setIdentity(bb2.nPatches(), bb2.size(), 1);  // bb2.nPatches() == 1
      mapperBdy.markBoundary(0, bnd, 0);
      mapperBdy.finalize();

      gsExprAssembler<real_t> A(1,1);
      A.setIntegrationElements(dbasis);

      auto G = A.getMap(mp);
      auto uu = A.getSpace(bb2);
      auto g_bdy = A.getBdrFunction(G);

      uu.setupMapper(mapperBdy);
      gsMatrix<real_t> & fixedDofs_A = const_cast<expr::gsFeSpace<real_t>&>(uu).fixedPart();
      fixedDofs_A.setZero( uu.mapper().boundarySize(), 1 );

      real_t lambda = 1e-5;

      A.initSystem();
      A.assembleBdr(bc.get("Dirichlet"), uu * uu.tr() * meas(G));
      A.assembleBdr(bc.get("Dirichlet"), uu * g_bdy * meas(G));
      A.assembleBdr(bc.get("Neumann"),
        lambda * (igrad(uu, G) * nv(G).normalized()) *
        (igrad(uu, G) * nv(G).normalized()).tr() * meas(G));
      A.assembleBdr(bc.get("Neumann"),
        lambda *  (igrad(uu, G) * nv(G).normalized()) * (g_bdy.tr() * nv(G).normalized()) * meas(G));

      gsSparseSolver<real_t>::SimplicialLDLT solver;
      solver.compute( A.matrix() );
      gsMatrix<real_t> & fixedDofs = const_cast<expr::gsFeSpace<real_t>& >(u).fixedPart();
      fixedDofs = solver.solve(A.rhs());
    };

  protected:

    gsMultiPatch<T> & m_mp;
    gsMultiBasis<T> & m_mb;
    
    gsBoundaryConditions<T> m_bc;
    gsFunctionExpr<T> m_f;
    gsFunctionExpr<T> m_ms;

    gsExprAssembler<T> m_A;
    gsExprEvaluator<T> m_ev;

    gsMatrix<T> m_solVector;
    
  }; // class gsBiharmonicAssembler2
} // namespace gismo
