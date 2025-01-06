/** @file gsBiharmonic_test.cpp

    @brief Testing the unstructured spline constructions with the biharmonic problem

    This is an example unit test that doesn't really do anything useful.
    It is here as a reference for you when creating additional unit tests.
    For additional reference information, see the "test.h" header.

    == BASIC REFERENCE ==
         - TEST(NAME_OF_TEST) { body_of_test }
         - TEST_FIXTURE(NAME_OF_FIXTURE,NAME_OF_TEST){ body_of_test }

    == CHECK MACRO REFERENCE ==
         - CHECK(EXPR);
         - CHECK_EQUAL(EXPECTED,ACTUAL);
         - CHECK_CLOSE(EXPECTED,ACTUAL,EPSILON);
         - CHECK_ARRAY_EQUAL(EXPECTED,ACTUAL,LENGTH);
         - CHECK_ARRAY_CLOSE(EXPECTED,ACTUAL,LENGTH,EPSILON);
         - CHECK_ARRAY2D_EQUAL(EXPECTED,ACTUAL,ROWCOUNT,COLCOUNT);
         - CHECK_ARRAY2D_CLOSE(EXPECTED,ACTUAL,ROWCOUNT,COLCOUNT,EPSILON);
         - CHECK_THROW(EXPR,EXCEPTION_TYPE_EXPECTED);

    == TIME CONSTRAINTS ==
         - UNITTEST_TIME_CONSTRAINT(TIME_IN_MILLISECONDS);
         - UNITTEST_TIME_CONSTRAINT_EXEMPT();

    == MORE INFO ==
         See: https://unittest-cpp.github.io/

    Author(s): P. Weinmueller
 **/

#include "gismo_unittest.h"       // Brings in G+Smo and the UnitTest++ framework

#include <gsUnstructuredSplines/src/gsApproxC1Spline.h>

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
}

void gsDirichletNeumannValuesL2Projection(gsMultiPatch<> & mp, gsMultiBasis<> & dbasis, gsBoundaryConditions<> & bc,
                                           gsMappedBasis<2,real_t> & bb2, const expr::gsFeSpace<real_t> & u)
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
                  lambda * (igrad(uu, G) * nv(G).normalized()) * (igrad(uu, G) * nv(G).normalized()).tr() * meas(G));
    A.assembleBdr(bc.get("Neumann"),
                  lambda *  (igrad(uu, G) * nv(G).normalized()) * (g_bdy.tr() * nv(G).normalized()) * meas(G));

    gsSparseSolver<real_t>::SimplicialLDLT solver;
    solver.compute( A.matrix() );
    gsMatrix<real_t> & fixedDofs = const_cast<expr::gsFeSpace<real_t>& >(u).fixedPart();
    fixedDofs = solver.solve(A.rhs());
}

void runBiharmonicTest ()
{
  gsInfo << "Test loaded successful\n";

  real_t convratetolh2 = 1.85; // Convergence rate for H2 should be around 2
  real_t convratetolh1 = 2.85; // Convergence rate for H1 should be around 3
  real_t convratetoll2 = 3.75;  // Convergence rate for L2 should be around 4
  
  bool second = false;

  index_t numRefine  = 3;
  index_t degree = 3;
  index_t smoothness = 2;
  
  gsMultiPatch<real_t> mp;
  gsBoundaryConditions<> bc;
  gsFunctionExpr<real_t> f, ms;
    
  
  std::string string_geo = "planar/1p_square.xml";
  gsInfo << "Filedata: " << string_geo << "\n";
  gsReadFile<>(string_geo, mp);
  mp.clearTopology();
  mp.computeTopology();

  gsFunctionExpr<>source("256*pi*pi*pi*pi*(4*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))",2);
  f.swap(source);
  gsInfo << "Source function " << f << "\n";
  
  gsFunctionExpr<> solution("(cos(4*pi*x) - 1) * (cos(4*pi*y) - 1)",2);
  ms.swap(solution);
  gsInfo << "Exact function " << ms << "\n";

  for (gsMultiPatch<>::const_biterator bit = mp.bBegin(); bit != mp.bEnd(); ++bit) {
    // Laplace
    gsFunctionExpr<> laplace ("-16*pi*pi*(2*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))",2);
    
    // Neumann
    gsFunctionExpr<> sol1der("-4*pi*(cos(4*pi*y) - 1)*sin(4*pi*x)",
			     "-4*pi*(cos(4*pi*x) - 1)*sin(4*pi*y)", 2);
    
    bc.addCondition(*bit, condition_type::dirichlet, ms);
    if (second)
	bc.addCondition(*bit, condition_type::laplace, laplace);
    else
      bc.addCondition(*bit, condition_type::neumann, sol1der);  
  }
  bc.setGeoMap(mp);
  gsInfo << "Boundary conditions:\n" << bc << "\n";

  gsMultiBasis<real_t> dbasis(mp, false);//true: poly-splines (not NURBS)
  dbasis.setDegree( degree); // preserve smoothness

  if (dbasis.basis(0).numElements() < 4) 
    dbasis.uniformRefine(1, degree-smoothness);

  //! [Problem setup]
  gsExprAssembler<real_t> A(1,1);

  // Elements used for numerical integration
  A.setIntegrationElements(dbasis);
  gsExprEvaluator<real_t> ev(A);
  
  // Set the geometry map
  auto G = A.getMap(mp);
  
  // Set the source term
  auto ff = A.getCoeff(f, G); // Laplace example
  
  // Set the discretization space
  gsMappedBasis<2,real_t> bb2;
  auto u = A.getSpace(bb2);

  // Solution vector and solution variable
  gsMatrix<real_t> solVector;
  auto u_sol = A.getSolution(u, solVector);
  
  // Recover manufactured solution
  auto u_ex = ev.getVariable(ms, G);
  //! [Problem setup]

  //! [Solver loop]
  gsVector<real_t> l2err(numRefine+1), h1err(numRefine+1), h2err(numRefine+1),
    IFaceErr(numRefine+1), meshsize(numRefine+1), dofs(numRefine+1);
 
  gsInfo<< "(dot1=approxC1construction, dot2=assembled, dot3=solved, dot4=got_error)\n"
    "\nDoFs: ";
  double setup_time(0), ma_time(0), slv_time(0), err_time(0);
  gsStopwatch timer;
  for (int r=0; r<=numRefine; ++r) {
    dbasis.uniformRefine(1,degree -smoothness);
    meshsize[r] = dbasis.basis(0).getMinCellLength();
    
    // The approx. C1 space
    gsApproxC1Spline<2,real_t> approxC1(mp,dbasis);
    approxC1.update(bb2);
    gsInfo<< "." <<std::flush; // Approx C1 construction done
    
    // Setup the mapper
    gsDofMapper map;
    setMapperForBiharmonic(bc, bb2,map);
    
    // Setup the system
    u.setupMapper(map);
    gsDirichletNeumannValuesL2Projection(mp, dbasis, bc, bb2, u);
  
    // Initialize the system
    A.initSystem();
    setup_time += timer.stop();
    
    dofs[r] = A.numDofs();
    gsInfo<< A.numDofs() <<std::flush;
    
    timer.restart();
    // Compute the system matrix and right-hand side
    A.assemble(ilapl(u, G) * ilapl(u, G).tr() * meas(G),u * ff * meas(G));
    
    // Enforce Laplace conditions to right-hand side
    auto g_L = A.getBdrFunction(G); // Set the laplace bdy value
    //auto g_L = A.getCoeff(laplace, G);
    A.assembleBdr(bc.get("Laplace"), (igrad(u, G) * nv(G)) * g_L.tr() );
    
    ma_time += timer.stop();
    gsInfo<< "." <<std::flush;// Assemblying done
    
    timer.restart();
    gsSparseSolver<real_t>::SimplicialLDLT solver;
    solver.compute( A.matrix() );
    solVector = solver.solve(A.rhs());
    
    slv_time += timer.stop();
    gsInfo<< "." <<std::flush; // Linear solving done
    
    timer.restart();
    //linferr[r] = ev.max( f-s ) / ev.max(f);
    
    l2err[r]= math::sqrt( ev.integral( (u_ex - u_sol).sqNorm() * meas(G) ) );
    h1err[r]= l2err[r] +
      math::sqrt(ev.integral( ( igrad(u_ex) - igrad(u_sol,G) ).sqNorm() * meas(G) ));
    
    h2err[r]= h1err[r] +
      math::sqrt(ev.integral( ( ihess(u_ex) - ihess(u_sol,G) ).sqNorm() * meas(G) )); 

    gsMatrix<real_t> solFull;
    u_sol.extractFull(solFull);
    gsMappedSpline<2, real_t> mappedSpline(bb2, solFull);
    
    auto ms_sol = A.getCoeff(mappedSpline);
    IFaceErr[r] = math::sqrt(ev.integralInterface(((igrad(ms_sol.left(), G.left()) -
						    igrad(ms_sol.right(), G.right())) *
						   nv(G).normalized()).sqNorm() * meas(G)));
 
    gsInfo<< ". " <<std::flush; // Error computations done
  } //for loop
    //! [Solver loop]

  timer.stop();
  gsInfo<<"\n\nTotal time: "<< setup_time+ma_time+slv_time+err_time <<"\n";
  gsInfo<<"     Setup: "<< setup_time <<"\n";
  gsInfo<<"  Assembly: "<< ma_time    <<"\n";
  gsInfo<<"   Solving: "<< slv_time   <<"\n";
  gsInfo<<"     Norms: "<< err_time   <<"\n";
  
    //! [Error and convergence rates]
  real_t convratelast = math::log(h2err[h2err.size()-2]/h2err[h2err.size()-1])/ std::log(2.0);
  CHECK( convratelast > convratetolh2 );

  convratelast = math::log(h1err[h1err.size()-2]/h1err[h1err.size()-1])/ std::log(2.0);
  CHECK( convratelast > convratetolh1 );

  convratelast = math::log(l2err[l2err.size()-2]/l2err[l2err.size()-1])/ std::log(2.0);
  CHECK( convratelast > convratetoll2 );
  
  /*
  gsInfo<< "\nL2 error: "<<std::scientific<<std::setprecision(3)<<l2err.transpose()<<"\n";
  gsInfo<< "H1 error: "<<std::scientific<<h1err.transpose()<<"\n";
  gsInfo<< "H2 error: "<<std::scientific<<h2err.transpose()<<"\n";
  gsInfo<< "Deriv Interface error: "<<std::scientific<<IFaceErr.transpose()<<"\n";
  
  if (numRefine>0) {
    gsInfo<< "EoC (L2): " << std::fixed<<std::setprecision(2)
	  <<  ( l2err.head(numRefine).array()  /
		  l2err.tail(numRefine).array() ).log().transpose() / std::log(2.0)
	  <<"\n";
    
    gsInfo<<   "EoC (H1): "<< std::fixed<<std::setprecision(2)
	  <<( h1err.head(numRefine).array() /
	      h1err.tail(numRefine).array() ).log().transpose() / std::log(2.0) <<"\n";
    
    gsInfo<<   "EoC (H2): "<< std::fixed<<std::setprecision(2)
	    <<( h2err.head(numRefine).array() /
		h2err.tail(numRefine).array() ).log().transpose() / std::log(2.0) <<"\n";
    
    gsInfo<<   "EoC (Iface): "<< std::fixed<<std::setprecision(2)
	  <<( IFaceErr.head(numRefine).array() /
	      IFaceErr.tail(numRefine).array() ).log().transpose() / std::log(2.0) <<"\n";
  }
  */
  //! [Error and convergence rates]
}


SUITE(gsBiharmonic_test)          // The suite should have the same name as the file
{

    TEST(approxC1)               // Declares a test named "gsBiharmonic_test:approxC1"
    {
        CHECK(true);              // We certainly hope that true is true
        CHECK_EQUAL(2,1+1);       // The value 1+1 should equal 2

        //CHECK_EQUAL(3,1+1);     // The value 1+1 should NOT equal 3

        int x[] = {1,2,3};
        int y[] = {1,2,3};
        CHECK_ARRAY_EQUAL(x,y,3); // These arrays of length 3 are equal

	      runBiharmonicTest();
	
        double a = 1.51;
        double b = 1.52;
        CHECK_CLOSE(a,b,0.1);     // These equal within 0.1
    }

    TEST(ASG1)
    {
      // CHECK_EQUAL(3,1+1);     // The value 1+1 should NOT equal 3
    }

    // ...

}
