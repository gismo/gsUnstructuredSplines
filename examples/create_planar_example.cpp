/** @file biharmonic_example.cpp

    @brief A Biharmonic example for a single patch.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/

# include <gismo.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    bool mesh = false;
    bool save = false;

    index_t numRefine  = 0;
    index_t degree = 3;
    index_t smoothness = 2;

    std::string fn;
    std::string geometry = "g1000";

    gsCmdLine cmd("Example for creating surface geometries with L2 projection.");
    // Flags related to the discrete settings
    cmd.addInt( "p", "degree", "Set the polynomial degree of the basis.", degree );
    cmd.addInt( "s", "smoothness", "Set the smoothness of the basis.",  smoothness );
    cmd.addInt( "r", "numRefine", "Number of refinement-loops.",  numRefine );

    cmd.addString("f", "file", "Input geometry file (with .xml)", fn);
    cmd.addString( "g", "geometry", "Input geometry file",  geometry );

    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("mesh", "Plot the mesh", mesh);

    cmd.addSwitch("save", "Save the solution", save);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Read Argument inputs]
    gsMultiPatch<real_t> mp, mp_init;
    std::string string_geo;
    if (fn.empty())
        string_geo = "planar/geometries/" + geometry + ".xml";
    else
        string_geo = fn;

    gsInfo << "Filedata: " << string_geo << "\n";
    gsReadFile<>(string_geo, mp);
    mp.clearTopology();
    mp.computeTopology();

    gsReadFile<>(string_geo, mp_init);
    mp_init.clearTopology();
    mp_init.computeTopology();
    //! [Read geometry]

    gsMultiBasis<> basis(mp);
    basis.degreeIncrease(degree-mp.patch(0).degree(0));
    for (int r =0; r < numRefine; ++r)
        basis.uniformRefine(1, degree-smoothness);


    gsInfo << "basis: " << basis.basis(0) << "\n";
    //! [Problem setup]
    gsExprAssembler<real_t> A(1,1);
    //gsInfo<<"Active options:\n"<< A.options() <<"\n";

    // Elements used for numerical integration
    A.setIntegrationElements(basis);
    gsExprEvaluator<real_t> ev(A);

    // Set the geometry map
    auto G = A.getMap(mp_init);

    // Set the source term
    auto ff = A.getCoeff(mp_init); // Laplace example

    // Set the discretization space
    auto u = A.getSpace(basis,2);

    // Solution vector and solution variable
    gsMatrix<real_t> solVector;
    auto u_sol = A.getSolution(u, solVector);
    gsMultiPatch<real_t> sol_coarse;

    // Recover manufactured solution
    auto u_ex = ev.getVariable(mp_init, G);
    //! [Problem setup]

#ifdef _OPENMP
    gsInfo << "Available threads: "<< omp_get_max_threads() <<"\n";
#endif

    //! [Solver loop]
    gsSparseSolver<real_t>::SimplicialLDLT solver;

    // Setup the system
    u.setup(0);

    // Initialize the system
    A.initSystem();

    // Compute the system matrix and right-hand side
    A.assemble(u * u.tr() * meas(G), // + 1e-5 * igrad(u, G) * igrad(u, G).tr() * meas(G) + 1e-5 * ilapl(u, G) * ilapl(u, G).tr() * meas(G),
               u * ff * meas(G));

    gsInfo << "." << std::flush;// Assemblying done

    solver.compute( A.matrix() );
    solVector = solver.solve(A.rhs());
    gsDebugVar(solVector);

    gsMultiPatch<> surface;
    u_sol.extract(surface);

    gsWriteParaview(surface, "solution_planar",2000);
    gsWriteParaview(mp_init, "initial_planar",2000);

    gsFileData<>fd;
    fd << surface;
    fd.save("new_geometry.xml");


//
//    real_t l2err = math::sqrt( ev.integral( (u_ex - u_sol).sqNorm() * meas(G) ) ); // / ev.integral(f.sqNorm()*meas(G)) );
//    gsInfo<< "\nL2 error: "<<std::scientific<<std::setprecision(3)<<l2err<<"\n";
//
//    gsMultiPatch<> surface;
//    u_sol.extract(surface);
//    auto ms_sol = A.getCoeff(surface);
//    real_t IFaceErr = math::sqrt(ev.integralInterface((( igrad(ms_sol.left(), G.left()) -
//                                                    igrad(ms_sol.right(), G.right())) *
//                                                            nv(G).normalized()).sqNorm() * meas(G)));
//    gsInfo<< "\nIFaceErr error: "<<std::scientific<<std::setprecision(3)<<IFaceErr<<"\n";
//
//    if (save)
//    {
//        gsMultiPatch<> mp_surf(mp);
//        mp_surf.embed(3);
//        for (size_t np = 0; np < mp_surf.nPatches(); np++)
//        {
//            gsMatrix<> coefs_patch;
//            u_sol.extract(coefs_patch, np);
//            mp_surf.patch(np).coefs().col(2) = coefs_patch;
//        }
//
//        gsFileData<> fd;
//        fd << mp_surf;
//        fd.save("surface");
//
//        gsMatrix<> points;
//        points.setZero(2,2);
//        points(0,0) = 1.0;
//        gsInfo << "Jac left: " << mp_surf.patch(0).jacobian(points.col(0)) << "\n";
//        gsInfo << "Jac right: " << mp_surf.patch(1).jacobian(points.col(1)) << "\n";
//
//        for (gsBoxTopology::const_iiterator it = mp_surf.interfaces().begin();
//             it != mp_surf.interfaces().end(); ++it )
//        {
//            const boundaryInterface &iFace = *it;
//            const index_t patch1 = iFace.first().patch;
//            const index_t patch2 = iFace.second().patch;
//
//            gsAffineFunction<real_t> interfaceMap(iFace.dirMap(), iFace.dirOrientation(),
//                                             mp_surf.basis(patch1).support(),
//                                             mp_surf.basis(patch2).support());
//
//            const index_t dir1 = iFace.first().side().index() < 3 ? 1 : 0;
//            const index_t dir2 = iFace.second().side().index() < 3 ? 1 : 0;
//
//            index_t N = 5;
//            gsMatrix<> pointsL, pointsR;
//            if (iFace.first().side().index() == 1 || iFace.first().side().index() == 3)
//                pointsL.setZero(2,N);
//            else
//                pointsL.setOnes(2,N);
//            gsVector<> pp;
//            pp.setLinSpaced(N,0,1);
//            pointsL.row(dir1) = pp;
//            interfaceMap.eval_into(pointsL, pointsR);
//
//            gsDebugVar(pointsL);
//            gsDebugVar(pointsR);
//
//            gsDebugVar(dir2);
//            gsDebugVar(iFace.second().side().index());
//
//            for (index_t i = 0; i < pointsL.cols(); i++)
//            {
//                gsMatrix<> mat_det(3,3);
//                mat_det.col(0) = mp_surf.patch(patch1).jacobian(pointsL.col(i)).col(1-dir1);
//                mat_det.col(1) = mp_surf.patch(patch2).jacobian(pointsR.col(i)).col(1-dir2);
//                mat_det.col(2) = mp_surf.patch(patch2).jacobian(pointsR.col(i)).col(dir2);
//                gsInfo << "Det: " << mat_det.determinant() << "\n";
//                //gsInfo << "Jac right: " << mp_surf.patch(patch2).jacobian(pointsR.col(i)) << "\n";
//            }
//        }
//
//        gsWriteParaview(mp_surf, "surface", 2000);
//    }

    //! [Export visualization in ParaView]
    if (plot)
    {
        // Write approximate and exact solution to paraview files
        gsInfo << "Plotting in Paraview...\n";
        ev.options().setSwitch("plot.elements", mesh);
        ev.writeParaview( u_sol   , G, "solution");
        ev.writeParaview( u_ex - u_sol   , G, "solution_pointwise");
        gsInfo << "Saved with solution.pvd \n";
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";
    //! [Export visualization in ParaView]

    return  EXIT_SUCCESS;
}
