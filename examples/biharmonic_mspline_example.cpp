/** @file biharmonic_mspline_example.cpp

    @brief Tutorial on how to use expression assembler and 
           gsMappedSpline to solve the Biharmonic equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/

//! [Include namespace]
#include <gismo.h>

#include <gsAssembler/gsBiharmonicExprAssembler.h>
using namespace gismo;
//! [Include namespace]

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot  = false;
    bool mesh  = false;

    bool last = false;
    bool second = false;

    std::string xml;

    gsCmdLine cmd("Tutorial on solving a Biharmonic problem with different spaces.");
    // Flags related to the problem (default: first biharmonic problem)
    cmd.addSwitch("second", "Solve the second biharmonic problem", second);

    // Flags related to the input basis via mspline
    cmd.addString("x", "xml", "XML File for the basis", xml);

    // Flags related to the output
    cmd.addSwitch("last", "Solve solely for the last level of h-refinement", last);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("mesh", "Plot the mesh", mesh);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Initialize data]
    gsMultiPatch<real_t> mp;
    gsSparseMatrix<real_t> cf;
    gsBoundaryConditions<> bc;
    gsFunctionExpr<real_t> f, ms;
    gsOptionList optionList, Aopt;
    gsMultiBasis<real_t> dbasis;
    gsFileData<> fd;
    //! [Initialize data]

    //! [Read Argument inputs]
    if (!xml.empty())
    {
        //! [Read geometry]
        fd.read(xml);

	// GEOMETRY
	fd.getFirst(mp);
        mp.clearTopology();
        mp.fixOrientation();
        mp.computeTopology();

	fd.getFirst(cf); // Coefficient matrix
	fd.getFirst(dbasis); // Standard basis functions (bspline basis functions)

        gsFunctionExpr<>source("256*pi*pi*pi*pi*(4*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))",2);
        // gsFunctionExpr<>source("1*pi*pi*pi*pi*(4*cos(1*pi*x)*cos(1*pi*y) - cos(1*pi*x) - cos(1*pi*y))",2);
        f.swap(source);
        gsInfo << "Source function " << f << "\n";

        gsFunctionExpr<> solution("(cos(4*pi*x) - 1) * (cos(4*pi*y) - 1)",2);
        // gsFunctionExpr<> solution("(cos(1*pi*x) - 1) * (cos(1*pi*y) - 1)",2);
        ms.swap(solution);
        gsInfo << "Exact function " << ms << "\n";

        //! [Boundary condition]
        for (gsMultiPatch<>::const_biterator bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
        {
            // Laplace
            gsFunctionExpr<> laplace ("-16*pi*pi*(2*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))",2);
            // gsFunctionExpr<> laplace ("-1*pi*pi*(2*cos(1*pi*x)*cos(1*pi*y) - cos(4*pi*x) - cos(4*pi*y))",2);

            // Neumann
            gsFunctionExpr<> sol1der("-4*pi*(cos(4*pi*y) - 1)*sin(4*pi*x)",
                                     "-4*pi*(cos(4*pi*x) - 1)*sin(4*pi*y)", 2);
            // gsFunctionExpr<> sol1der("-1*pi*(cos(1*pi*y) - 1)*sin(1*pi*x)",
                                     // "-1*pi*(cos(1*pi*x) - 1)*sin(1*pi*y)", 2);


            bc.addCondition(*bit, condition_type::dirichlet, ms);
            if (second)
                bc.addCondition(*bit, condition_type::laplace, laplace);
            else
                bc.addCondition(*bit, condition_type::neumann, sol1der);

        }
        gsInfo << "Boundary conditions:\n" << bc << "\n";
        //! [Boundary condition]

        optionList = cmd;
        //gsInfo << "OptionList: " << optionList << "\n";

        gsInfo << "Finished\n";
    }
    else
        GISMO_ERROR("No XML file provided provided!");
    //! [Read XML file]

    // Set the discretization space
    gsMappedBasis<2,real_t> bb2;
    bb2.init(dbasis,cf);
    gsInfo << "The MappedBasis has " << bb2.size() << " basis functions for all patches! \n";

    //! [Solver loop]
    gsVector<real_t> l2err(1), h1err(1), h2err(1),
      IFaceErr(1), meshsize(1), dofs(1);
    
    gsInfo<< "(dot1=approxC1construction, dot2=assembled, dot3=solved, dot4=got_error)\n"
        "\nDoFs: ";
    double setup_time(0), ma_time(0), slv_time(0), err_time(0);
    gsStopwatch timer;
    gsFunctionSet<>::Ptr solution;
    gsBiharmonicExprAssembler<real_t> BA;

    int r = 0;      
    gsInfo<< "." <<std::flush; // Approx C1 construction done
    bc.setGeoMap(mp);

    std::vector<gsBasis<real_t> *> bases = bb2.getBasesCopy();
    gsMultiBasis<> intBasis;
    intBasis = gsMultiBasis<>(bases,mp.topology());

    BA = gsBiharmonicExprAssembler<real_t>(mp,intBasis,f,bc);
    gsOptionList exprAssemblerOpts = Aopt.wrapIntoGroup("ExprAssembler");
    BA.setOptions(exprAssemblerOpts);
    BA.setSpaceBasis(bb2);

    //! [Problem setup]
    setup_time += timer.stop();

    timer.restart();
    // Compute the system matrix and right-hand side
    BA.assemble();
    ma_time += timer.stop();

    dofs[r] = BA.numDofs();
    gsInfo<< BA.numDofs() <<std::flush;

    gsInfo<< "." <<std::flush;// Assemblying done

    timer.restart();
    gsSparseSolver<real_t>::SimplicialLDLT solver;

    solver.compute( BA.matrix() );
    gsMatrix<> solVector = solver.solve(BA.rhs());
    BA.constructSolution(solVector);
    solution = give(BA.getSolution());

    slv_time += timer.stop();
    gsInfo<< "." <<std::flush; // Linear solving done

    timer.restart();
    //linferr[r] = ev.max( f-s ) / ev.max(f);

    std::tie(l2err[r],h1err[r],h2err[r]) = BA.errors(solVector,ms);
    IFaceErr[r] = BA.interfaceError(solVector,ms);
	
    err_time += timer.stop();
    gsInfo<< ". " <<std::flush; // Error computations done
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
    //! [Error and convergence rates]

    //! [Export visualization in ParaView]
    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
	gsExprEvaluator<> ev;
	ev.setIntegrationElements(dbasis);
	auto u_sol = ev.getVariable(*solution);
	auto u_ex  = ev.getVariable(ms);
	auto G = ev.getMap(mp);
	
	ev.options().setSwitch("plot.elements", mesh);
	ev.options().setInt   ("plot.npts"    , 1000);
	ev.writeParaview( u_sol   , G, "solution");
	ev.writeParaview( u_ex    , G, "solution_ex");
	
        //gsWriteParaview( geom, "geom",1000,true);
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";
    //! [Export visualization in ParaView]
    
    return EXIT_SUCCESS;
}// end main
