/** @file biharmonic2_example.cpp

    @brief Tutorial on how to use expression assembler and the (approx.) C1 basis function
                to solve the Biharmonic equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/

//! [Include namespace]
#include <gismo.h>

#include <gsUnstructuredSplines/gsApproxC1Spline.h>
#include <gsUnstructuredSplines/gsDPatch.h>
#include <gsUnstructuredSplines/gsAlmostC1.h>

using namespace gismo;
//! [Include namespace]

/**
 * Smoothing method:
 * - s 0 == Approx C1 method
 * - s 1 == Nitsche's method
 * - s 2 == D-Patch's method
 */
enum MethodFlags
{
    APPROXC1       = 0 << 0, // Approx C1 Method
    DPATCH         = 1 << 0, // D-Patch
    ALMOSTC1       = 1 << 1, //
    //????      = 1 << 2, // ????
    //????      = 1 << 3, // ????
    // Add more [...]
};

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


int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    bool mesh = false;

    index_t method = 0;

    index_t numRefine  = 3;
    index_t degree = 3;
    index_t smoothness = 2;

    index_t gluingDataDegree = -1;
    index_t gluingDataSmoothness = -1;

    bool last = false;
    bool info = false;
    bool second = false;

    bool interpolation = false;

    std::string fn;
    std::string output;
    std::string geometry = "g1000";

    gsCmdLine cmd("Tutorial on solving a Biharmonic problem with different spaces.");
    // Flags related to the method (default: Approx C1 method)
    cmd.addInt( "m", "method", "The chosen method for the biharmonic problem", method );

    // Flags related to the problem (default: first biharmonic problem)
    cmd.addSwitch("second", "Solve the second biharmonic problem", second);

    // Flags related to the input/geometry
    cmd.addString( "f", "file", "Input geometry file from path (with .xml)", fn );
    cmd.addString( "g", "geometry", "Input geometry file",  geometry );

    // Flags related to the discrete settings
    cmd.addInt( "p", "degree", "Set the polynomial degree of the basis.", degree );
    cmd.addInt( "s", "smoothness", "Set the smoothness of the basis.",  smoothness );
    cmd.addInt( "r", "numRefine", "Number of refinement-loops.",  numRefine );

    // Flags related to the output
    cmd.addSwitch("last", "Solve solely for the last level of h-refinement", last);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("mesh", "Plot the mesh", mesh);

    cmd.addString("o", "output", "Output in xml (for python)", output);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Initialize data]
    gsMultiPatch<real_t> mp;
    gsBoundaryConditions<> bc;
    gsFunctionExpr<real_t> f, ms;
    gsOptionList optionList;
    //! [Initialize data]

    //! [Read Argument inputs]
    std::string string_geo;
    if (fn.empty())
        string_geo = "planar/geometries/" + geometry + ".xml";
    else
        string_geo = fn;

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

    //! [Boundary condition]
    for (gsMultiPatch<>::const_biterator bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
    {
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
    //! [Boundary condition]

    optionList = cmd;
    gsInfo << "OptionList: " << optionList << "\n";
    gsInfo << "Finished\n";
    //! [Read Argument inputs]

    //! [Refinement]
    gsMultiBasis<real_t> dbasis(mp, true);//true: poly-splines (not NURBS)

    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    dbasis.setDegree( degree); // preserve smoothness
    //dbasis.degreeElevate(degree- mp.patch(0).degree(0));

    if (method == MethodFlags::DPATCH || method == MethodFlags::ALMOSTC1)
        mp.degreeElevate(degree-mp.patch(0).degree(0));

    // h-refine each basis
    if (last)
    {
        for (int r =0; r < numRefine; ++r)
        {
            dbasis.uniformRefine(1, degree-smoothness);
            mp.uniformRefine(1, degree-smoothness);
        }
        numRefine = 0;
    }

    // Assume that the condition holds for each patch TODO
    // Refine once
    if (dbasis.basis(0).numElements() < 4)
    {
        dbasis.uniformRefine(1, degree-smoothness);
        if (method == MethodFlags::DPATCH || method == MethodFlags::ALMOSTC1)
            mp.uniformRefine(1, degree-smoothness);
    }

    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
#ifdef _OPENMP
    gsInfo<< "Available threads: "<< omp_get_max_threads() <<"\n";
#endif

    //! [Refinement]

    //! [Problem setup]
    gsExprAssembler<real_t> A(1,1);
    //gsInfo<<"Active options:\n"<< A.options() <<"\n";

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

    // The approx. C1 space
    gsApproxC1Spline<2,real_t> approxC1(mp,dbasis);
    approxC1.options().setSwitch("info",info);
    approxC1.options().setSwitch("plot",plot);
    approxC1.options().setSwitch("interpolation",interpolation);
    approxC1.options().setSwitch("second",second);
    approxC1.options().setInt("gluingDataDegree",gluingDataDegree);
    approxC1.options().setInt("gluingDataSmoothness",gluingDataSmoothness);

    // D-Patch
    gsMultiPatch<> geom = mp;

    // Solution vector and solution variable
    gsMatrix<real_t> solVector;
    auto u_sol = A.getSolution(u, solVector);

    // Recover manufactured solution
    auto u_ex = ev.getVariable(ms, G);
    //! [Problem setup]

    // For Nitsche
    gsMatrix<real_t> mu_interfaces(mp.nInterfaces(),1);

    //! [Solver loop]
    gsVector<real_t> l2err(numRefine+1), h1err(numRefine+1), h2err(numRefine+1),
            IFaceErr(numRefine+1), meshsize(numRefine+1), dofs(numRefine+1),
            cond_num(numRefine+1);
    gsMatrix<real_t> penalty(numRefine+1, mp.nInterfaces());
    gsInfo<< "(dot1=approxC1construction, dot2=assembled, dot3=solved, dot4=got_error)\n"
        "\nDoFs: ";
    double setup_time(0), ma_time(0), slv_time(0), err_time(0);
    gsStopwatch timer;
    for (int r=0; r<=numRefine; ++r)
    {
        if (method == MethodFlags::APPROXC1)
        {
            dbasis.uniformRefine(1,degree -smoothness);
            meshsize[r] = dbasis.basis(0).getMinCellLength();
            approxC1.update(bb2);
        }
        else if (method == MethodFlags::DPATCH)
        {
            geom.uniformRefine(1,degree-smoothness);
            dbasis.uniformRefine(1,degree-smoothness);

            meshsize[r] = dbasis.basis(0).getMinCellLength();

            gsSparseMatrix<real_t> global2local;
            gsDPatch<2,real_t> dpatch(geom);
            dpatch.matrix_into(global2local);
            global2local = global2local.transpose();
            mp = dpatch.exportToPatches();
            dbasis = dpatch.localBasis();
            bb2.init(dbasis,global2local);

            gsFileData<> fd;
            fd << mp;
            fd.save("mp_filedata.xml");
            fd.clear();
            fd << geom;
            fd.save("geom_filedata.xml");
        }
        else if (method == MethodFlags::ALMOSTC1)
        {
            geom.uniformRefine(1,degree-smoothness);
            dbasis.uniformRefine(1,degree-smoothness);

            meshsize[r] = dbasis.basis(0).getMinCellLength();

            gsSparseMatrix<real_t> global2local;
            gsAlmostC1<2,real_t> almostC1(geom);
            almostC1.matrix_into(global2local);
            global2local = global2local.transpose();
            mp = almostC1.exportToPatches();
            dbasis = almostC1.localBasis();
            bb2.init(dbasis,global2local);
        }
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

        l2err[r]= math::sqrt( ev.integral( (u_ex - u_sol).sqNorm() * meas(G) ) ); // / ev.integral(f.sqNorm()*meas(G)) );
        h1err[r]= l2err[r] +
            math::sqrt(ev.integral( ( igrad(u_ex) - igrad(u_sol,G) ).sqNorm() * meas(G) )); // /ev.integral( igrad(f).sqNorm()*meas(G) ) );

        h2err[r]= h1err[r] +
                 math::sqrt(ev.integral( ( ihess(u_ex) - ihess(u_sol,G) ).sqNorm() * meas(G) )); // /ev.integral( ihess(f).sqNorm()*meas(G) )


        gsMatrix<real_t> solFull;
        u_sol.extractFull(solFull);
        gsMappedSpline<2, real_t> mappedSpline(bb2, solFull);

        auto ms_sol = A.getCoeff(mappedSpline);
        IFaceErr[r] = math::sqrt(ev.integralInterface(((igrad(ms_sol.left(), G.left()) -
                                                            igrad(ms_sol.right(), G.right())) *
                                                           nv(G).normalized()).sqNorm() * meas(G)));

        err_time += timer.stop();
        gsInfo<< ". " <<std::flush; // Error computations done
    } //for loop
    //! [Solver loop]


    timer.stop();
    gsInfo<<"\n\nTotal time: "<< setup_time+ma_time+slv_time+err_time <<"\n";
    gsInfo<<"     Setup: "<< setup_time <<"\n";
    gsInfo<<"  Assembly: "<< ma_time    <<"\n";
    gsInfo<<"   Solving: "<< slv_time   <<"\n";
    gsInfo<<"     Norms: "<< err_time   <<"\n";

    gsInfo<< "\nMesh-size: " << meshsize.transpose() << "\n";

    //! [Error and convergence rates]
    gsInfo<< "\nL2 error: "<<std::scientific<<std::setprecision(3)<<l2err.transpose()<<"\n";
    gsInfo<< "H1 error: "<<std::scientific<<h1err.transpose()<<"\n";
    gsInfo<< "H2 error: "<<std::scientific<<h2err.transpose()<<"\n";
    gsInfo<< "Deriv Interface error: "<<std::scientific<<IFaceErr.transpose()<<"\n";

    if (!last && numRefine>0)
    {
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
    //! [Error and convergence rates]

    //! [Export visualization in ParaView]
    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview( mp, "geom",1000,true);
        ev.options().setSwitch("plot.elements", mesh);
        ev.options().setInt   ("plot.npts"    , 1000);
        ev.writeParaview( u_sol   , G, "solution");
        //ev.writeParaview( u_ex    , G, "solution_ex");
        //ev.writeParaview( grad(s), G, "solution_grad");
        //ev.writeParaview( grad(f), G, "solution_ex_grad");
        ev.writeParaview( (u_ex-u_sol), G, "error_pointwise");
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";
    //! [Export visualization in ParaView]

        //! [Export data to xml]
    if (!output.empty())
    {
        index_t cols = 7;
        gsMatrix<real_t> error_collection(l2err.rows(), cols);
        error_collection.col(0) = meshsize;
        error_collection.col(1) = dofs;
        error_collection.col(2) = l2err;
        error_collection.col(3) = h1err;
        error_collection.col(4) = h2err;
        error_collection.col(5) = IFaceErr;
        error_collection.col(6) = cond_num;

        gsFileData<real_t> xml_out;
        xml_out << error_collection;
        // Add solution
        // [...]
        xml_out.save(output + ".xml");
        gsInfo << "XML saved to " + output + ".xml" << "\n";
    }
    //! [Export data to xml]


    return EXIT_SUCCESS;
}// end main
