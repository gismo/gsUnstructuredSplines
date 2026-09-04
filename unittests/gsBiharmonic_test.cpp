/** @file gsBiharmonic_test.cpp

    @brief Testing the unstructured spline constructions with the biharmonic problem

    Solves the biharmonic problem with a manufactured solution on unstructured
    spline spaces (Approx-C1 and AS-G1) and checks both the discretization-error
    convergence rates and the (approximate) C1/G1 continuity across patch
    interfaces, which is the whole point of these constructions.

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
#include <gsUnstructuredSplines/src/gsC1SurfSpline.h>

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

/// Which unstructured spline construction to build the biharmonic space from.
enum class BiharmonicMethod
{
    ApproxC1, // Approximate C1 (Farahat et al.): globally C1, works on any topology
    ASG1      // Analysis-suitable G1 (gsC1SurfSpline): only approximately G1 by
              // construction; requires interior regularity r = p-2
};

/// Parameters for one convergence run of the biharmonic test problem on an
/// unstructured spline space.
struct BiharmonicTestParams
{
    BiharmonicMethod method;
    std::string      geoFile;
    index_t          degree;
    index_t          smoothness;
    index_t          numRefine;

    // One-sided lower bounds on the observed EoC (estimated from the last two
    // refinement levels) for the discretization error in each norm.
    real_t tolL2;
    real_t tolH1;
    real_t tolH2;

    // Number of patch interfaces this fixture must produce once its topology
    // is (re)computed from the geometry; 0 for a single patch. Also gates the
    // interface-error check below: with no interfaces IFaceErr is trivially 0,
    // so asserting on it would be vacuous rather than a real continuity check.
    index_t expectedInterfaces;

    // Upper bound on the final-level (finest mesh) interface derivative jump
    // ||[[grad u]]||_{L2(interfaces)}. Ignored when expectedInterfaces == 0.
    real_t tolIFace;
};

void runBiharmonicTest (const BiharmonicTestParams & prm)
{
  gsInfo << "Test loaded successful\n";

  bool second = false;

  gsMultiPatch<real_t> mp;
  gsBoundaryConditions<> bc;
  gsFunctionExpr<real_t> f, ms;


  gsInfo << "Filedata: " << prm.geoFile << "\n";
  gsReadFile<>(prm.geoFile, mp);
  mp.clearTopology();
  mp.fixOrientation();
  mp.computeTopology();
  // Guard against a silently empty interface set: with zero interfaces the
  // IFaceErr check below would pass vacuously (integralInterface() == 0)
  // regardless of the construction's actual continuity.
  CHECK_EQUAL( prm.expectedInterfaces, static_cast<index_t>(mp.nInterfaces()) );

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
  dbasis.setDegree( prm.degree); // preserve smoothness

  // AS-G1 rebuilds the basis from the (possibly h-refined) geometry every
  // level -- see the method dispatch inside the loop below -- so the geometry
  // itself must first be raised to the target degree.
  if (prm.method == BiharmonicMethod::ASG1)
      mp.degreeElevate(prm.degree - mp.patch(0).degree(0));

  if (dbasis.basis(0).numElements() < 4)
  {
    dbasis.uniformRefine(1, prm.degree - prm.smoothness);
    if (prm.method == BiharmonicMethod::ASG1)
        mp.uniformRefine(1, prm.degree - prm.smoothness);
  }

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
  gsVector<real_t> l2err(prm.numRefine+1), h1err(prm.numRefine+1), h2err(prm.numRefine+1),
    IFaceErr(prm.numRefine+1), meshsize(prm.numRefine+1), dofs(prm.numRefine+1);

  gsInfo<< "(dot1=construction, dot2=assembled, dot3=solved, dot4=got_error)\n"
    "\nDoFs: ";
  double setup_time(0), ma_time(0), slv_time(0), err_time(0);
  gsStopwatch timer;
  for (int r=0; r<=prm.numRefine; ++r) {

    if (prm.method == BiharmonicMethod::ApproxC1)
    {
        dbasis.uniformRefine(1, prm.degree - prm.smoothness);
        meshsize[r] = dbasis.basis(0).getMinCellLength();

        // The approx. C1 space
        gsApproxC1Spline<2,real_t> approxC1(mp,dbasis);
        approxC1.update(bb2);
    }
    else // BiharmonicMethod::ASG1
    {
        // AS-G1 refines the geometry and rebuilds the basis from it directly
        // (mirrors biharmonic_planar_example.cpp's SURFASG1 branch), rather
        // than refining an independently-tracked basis.
        mp.uniformRefine(1, prm.degree - prm.smoothness);
        dbasis = gsMultiBasis<real_t>(mp);
        meshsize[r] = dbasis.basis(0).getMinCellLength();

        // The domain gsExprAssembler caches on setIntegrationElements() holds
        // non-owning pointers into the multibasis's knot vectors. dbasis was
        // just replaced wholesale above (unlike the Approx-C1 branch, which
        // refines the same basis objects in place), so those pointers now
        // target freed memory and the domain must be re-taken here.
        A.setIntegrationElements(dbasis);

        gsC1SurfSpline<2,real_t> smoothC1(mp,dbasis);
        smoothC1.init();
        smoothC1.compute();

        gsSparseMatrix<real_t> global2local = smoothC1.getSystem();
        global2local = global2local.transpose();
        gsMultiBasis<real_t> basis_temp;
        smoothC1.getMultiBasis(basis_temp);
        bb2.init(basis_temp, global2local);
    }
    gsInfo<< "." <<std::flush; // Construction done

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

    l2err[r]= math::sqrt( ev.integral( (u_ex - u_sol).sqNorm() * meas(G) ) );
    h1err[r]= l2err[r] +
      math::sqrt(ev.integral( ( igrad(u_ex) - igrad(u_sol,G) ).sqNorm() * meas(G) ));

    h2err[r]= h1err[r] +
      math::sqrt(ev.integral( ( ihess(u_ex) - ihess(u_sol,G) ).sqNorm() * meas(G) ));

    gsMatrix<real_t> solFull;
    u_sol.extractFull(solFull);
    gsMappedSpline<2, real_t> mappedSpline(bb2, solFull);

    // Jump of the gradient of the discrete solution across patch interfaces:
    // the direct measure of the construction's continuity. Zero for a
    // single-patch fixture (no interfaces to integrate over).
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

  gsInfo<< "\nL2 error: "<<std::scientific<<std::setprecision(3)<<l2err.transpose()<<"\n";
  gsInfo<< "H1 error: "<<std::scientific<<h1err.transpose()<<"\n";
  gsInfo<< "H2 error: "<<std::scientific<<h2err.transpose()<<"\n";
  gsInfo<< "Deriv Interface error: "<<std::scientific<<IFaceErr.transpose()<<"\n";

    //! [Error and convergence rates]
  real_t convratelast = math::log(h2err[h2err.size()-2]/h2err[h2err.size()-1])/ std::log(2.0);
  gsInfo<< "EoC (H2, last): "<<convratelast<<" (tol "<<prm.tolH2<<")\n";
  CHECK( convratelast > prm.tolH2 );

  convratelast = math::log(h1err[h1err.size()-2]/h1err[h1err.size()-1])/ std::log(2.0);
  gsInfo<< "EoC (H1, last): "<<convratelast<<" (tol "<<prm.tolH1<<")\n";
  CHECK( convratelast > prm.tolH1 );

  convratelast = math::log(l2err[l2err.size()-2]/l2err[l2err.size()-1])/ std::log(2.0);
  gsInfo<< "EoC (L2, last): "<<convratelast<<" (tol "<<prm.tolL2<<")\n";
  CHECK( convratelast > prm.tolL2 );

  gsInfo<< "Final IFaceErr: "<<IFaceErr[IFaceErr.size()-1]<<" (tol "<<prm.tolIFace<<")\n";
  if (prm.expectedInterfaces > 0)
      CHECK( IFaceErr[IFaceErr.size()-1] < prm.tolIFace );
  //! [Error and convergence rates]
}


SUITE(gsBiharmonic_test)          // The suite should have the same name as the file
{

    TEST(approxC1)               // Declares a test named "gsBiharmonic_test:approxC1"
    {
        // Single patch: no interfaces, so gsApproxC1Spline never enters the
        // interface-vertex branches of gsApproxC1Vertex/gsApproxC1Utils, but
        // this is the original, long-standing gate on the interior C1 space
        // and its arguments must stay exactly as they were.
        BiharmonicTestParams prm;
        prm.method             = BiharmonicMethod::ApproxC1;
        prm.geoFile             = "planar/1p_square.xml";
        prm.degree              = 3;
        prm.smoothness           = 2;
        prm.numRefine            = 3;
        prm.tolL2                = 3.75; // Convergence rate for L2 should be around 4
        prm.tolH1                = 2.85; // Convergence rate for H1 should be around 3
        prm.tolH2                = 1.85; // Convergence rate for H2 should be around 2
        prm.expectedInterfaces   = 0;
        prm.tolIFace             = 0.0;  // unused: expectedInterfaces==0 skips the check

        runBiharmonicTest(prm);
    }

    TEST(approxC1_multipatch)
    {
        // Three patches meeting at a single valence-3 extraordinary vertex:
        // this is what actually exercises the isInterface[dir] branches of
        // gsApproxC1Vertex/gsApproxC1Utils, which a single-patch fixture
        // cannot reach. Approx-C1 is exactly C1 by construction, so the
        // interface derivative jump should reach the double-precision floor.
        BiharmonicTestParams prm;
        prm.method             = BiharmonicMethod::ApproxC1;
        prm.geoFile              = "planar/3p_hexagon.xml";
        prm.degree               = 3;
        prm.smoothness            = 2;
        prm.numRefine             = 3;
        prm.tolL2                 = 3.75;
        prm.tolH1                 = 2.85;
        prm.tolH2                 = 1.85;
        prm.expectedInterfaces    = 3;
        prm.tolIFace              = 1e-10;

        runBiharmonicTest(prm);
    }

    TEST(ASG1)
    {
        // AS-G1 (gsC1SurfSpline) on the same valence-3 fixture, exercising
        // the edge and vertex constructions. AS-G1 requires interior knot
        // multiplicity 2 (regularity r = p-2, i.e. degree - smoothness == 2);
        // with multiplicity 1 the construction emits out-of-range sparse
        // matrix indices and aborts in set_from_triplets.
        BiharmonicTestParams prm;
        prm.method             = BiharmonicMethod::ASG1;
        prm.geoFile              = "planar/3p_hexagon.xml";
        prm.degree               = 3;
        prm.smoothness            = 1; // degree - smoothness == 2, as AS-G1 requires
        prm.numRefine             = 3;
        // Observed orders for AS-G1 with p=3 are ~4.0 (L2), ~3.0 (H1), ~2.0 (H2),
        // same as the (exactly C1) Approx-C1 cases above.
        prm.tolL2                 = 3.75;
        prm.tolH1                 = 2.85;
        prm.tolH2                 = 1.85;
        prm.expectedInterfaces    = 3;
        // AS-G1 is only approximately G1 by construction (the coupling across
        // an interface is enforced through a reduced gluing-data space, not
        // exactly); its interface error sits on a plateau rather than
        // converging to zero as the mesh is refined, so this checks that the
        // plateau is small, not that it decreases.
        prm.tolIFace              = 1e-5;

        runBiharmonicTest(prm);
    }

    // ...

}
