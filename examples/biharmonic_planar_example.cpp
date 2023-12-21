/** @file biharmonic2_example.cpp

    @brief Tutorial on how to use expression assembler and the (approx.) C1 basis function
                to solve the Biharmonic equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller, H.M. Verhelst
*/

//! [Include namespace]
#include <gismo.h>

#include <gsUnstructuredSplines/src/gsApproxC1Spline.h>
#include <gsUnstructuredSplines/src/gsDPatch.h>
#include <gsUnstructuredSplines/src/gsAlmostC1.h>
#include <gsUnstructuredSplines/src/gsC1SurfSpline.h>
#include <gsAssembler/gsBiharmonicExprAssembler.h>

#ifdef gsSpectra_ENABLED
#include <gsSpectra/gsSpectra.h>
#endif
#include <gsUtils/gsL2Projection.h>
using namespace gismo;
//! [Include namespace]

/**
 * Smoothing method:
 * - m 0 == Approx C1 method
 * - m 1 == D-Patch method
 * - m 2 == Almost C1 method
 * - m 3 == Nitsche's method
 */
enum MethodFlags
{
    DPATCH         = 1, // D-Patch
    APPROXC1       = 2, // Approx C1 Method
    SURFASG1       = 3, // AS-G1
    ALMOSTC1       = 4, // Almost C1
    NITSCHE        = 5, // Nitsche
    // Add more [...]
};

void writeToCSVfile(std::string name, gsMatrix<> matrix)
{
  std::ofstream file(name.c_str());

  for(int  i = 0; i < matrix.rows(); i++){
      for(int j = 0; j < matrix.cols(); j++){
         if(j+1 == matrix.cols()){
             file<<std::to_string(matrix(i,j));
         }else{
             file<<std::to_string(matrix(i,j))<<',';
         }
      }
      file<<'\n';
  }
}

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot  = false;
    bool nested  = false;
    bool plotApproxC1 = false;
    bool mesh  = false;

    index_t method = 0;

    index_t numRefine  = 3;
    index_t numRefineIni  = 1;
    index_t degree = 3;
    index_t smoothness = 2;

    index_t gluingDataDegree = -1;
    index_t gluingDataSmoothness = -1;

    bool last = false;
    bool info = false;
    bool second = false;
    bool cond = false;
    bool writeMatrix= false;
    bool interpolation = false;

    bool project    = false;

    index_t PiMat = 1;

    real_t penalty_init = -1.0;
    std::string xml;
    std::string output;
    std::string write;
    std::string assemberOptionsFile("options/assembler_options.xml");

    std::string fn;

    real_t scaling = 1.0;
    real_t lambda = 1e-5;
    real_t beta = 0.4;
    gsCmdLine cmd("Tutorial on solving a Biharmonic problem with different spaces.");
    // Flags related to the method (default: Approx C1 method)
    cmd.addInt( "m", "method", "The chosen method for the biharmonic problem", method );

    // Flags related to the problem (default: first biharmonic problem)
    cmd.addSwitch("second", "Solve the second biharmonic problem", second);

    // Perform nested refinement for D-Patch or Almost-C1
    cmd.addSwitch("nested", "Perform nested refinement for D-Patch or Almost-C1", nested);

    // Flags related to the input/geometry
    cmd.addString( "f", "file", "Input geometry file from path (with .xml)", fn );
    cmd.addString("x", "xml", "Use the input from the xml file", xml);

    // Flags related to the discrete settings
    cmd.addInt( "p", "degree", "Set the polynomial degree of the basis.", degree );
    cmd.addInt( "s", "smoothness", "Set the smoothness of the basis.",  smoothness );
    cmd.addInt( "r", "numRefine", "Number of refinement-loops.",  numRefine );
    cmd.addInt( "i", "numRefineIni", "Number initial of refinements.",  numRefineIni );
    cmd.addInt( "", "PiMat", "Pi matrix to use (0: Idempotent, 1: Non-negative",  PiMat );

    // Flags related to the approximate C1 method
    cmd.addInt( "P", "gluingDataDegree","Set the polynomial degree for the gluing data", gluingDataDegree );
    cmd.addInt( "R", "gluingDataSmoothness", "Set the smoothness for the gluing data",  gluingDataSmoothness );
    cmd.addSwitch("interpolation", "Compute the basis constructions with interpolation", interpolation);
    cmd.addSwitch("info", "Getting the information inside of Approximate C1 basis functions", info);
    cmd.addSwitch("plotApproxC1", "Plot the approximate C1 basis functions", plotApproxC1);

    cmd.addReal( "S", "scaling", "2D geometry scaling",  scaling);
    cmd.addReal( "L", "lambda", "Lambda for BC projection",  lambda);

    cmd.addReal( "B", "beta", "Beta for D-Patch", beta);


    // Flags related to Nitsche's method
    cmd.addReal( "y", "penalty", "Fixed Penalty value for Nitsche's method",  penalty_init);

    // Flags related to the output
    cmd.addSwitch("last", "Solve solely for the last level of h-refinement", last);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("mesh", "Plot the mesh", mesh);
    cmd.addSwitch("cond", "Estimate condition number (slow!)", cond);

    cmd.addString("O", "Aopt", "Assembler options file", assemberOptionsFile);
    cmd.addString("o", "output", "Output in xml (for python)", output);
    cmd.addString("w", "write", "Write to csv", write);

    cmd.addSwitch("writeMat", "Write projection matrix",writeMatrix);
    cmd.addSwitch( "project", "Project the geometry on the initial basis (D-Patch and almost-C1)", project );

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Initialize data]
    gsMultiPatch<real_t> mp;
    gsBoundaryConditions<> bc;
    gsFunctionExpr<real_t> f, ms;
    gsOptionList optionList, Aopt;
    //! [Initialize data]

    //! [Read Argument inputs]
    if (xml.empty() && !fn.empty())
    {
        //! [Read geometry]
        gsInfo << "Filedata: " << fn << "\n";
        gsReadFile<>(fn, mp);
        mp.clearTopology();
        mp.fixOrientation();
        mp.computeTopology();

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
        gsInfo << "OptionList: " << optionList << "\n";

        gsFileData<> fd(assemberOptionsFile); // "planar/biharmonic_pde/bvp1.xml"
        fd.getAnyFirst(Aopt);

        gsInfo<<"Assembler options:\n"<< Aopt <<"\n";

        gsInfo << "Finished\n";
    }
    else if (!xml.empty())
    {
        // id=0 Boundary
        // id=1 Source function
        // id=2 Optionlist
        // id=3 Exact solution
        // id=X Geometry (should be last!)
        gsFileData<> fd(xml); // "planar/biharmonic_pde/bvp1.xml"

        // Geometry
        fd.getAnyFirst(mp);
        mp.computeTopology();
        gsInfo << "Multipatch " << mp << "\n";

        // Functions
        fd.getId(1, f); // Source solution
        gsInfo << "Source function " << f << "\n";

        fd.getId(3, ms); // Exact solution
        gsInfo << "Exact function " << ms << "\n";

        // Boundary condition
        fd.getId(0, bc); // id=2: boundary conditions
        gsInfo << "Boundary conditions:\n" << bc << "\n";

        // Option list
        fd.getId(2, optionList); // id=100: assembler options
        gsInfo << "OptionList: " << optionList << "\n";

        degree = optionList.getInt("degree");
        smoothness = optionList.getInt("smoothness");
        numRefine = optionList.getInt("numRefine");

        gluingDataDegree = optionList.getInt("gluingDataDegree");
        gluingDataSmoothness = optionList.getInt("gluingDataSmoothness");

        method = optionList.getInt("method");

        penalty_init = optionList.getReal("penalty");

        //cond = optionList.getSwitch("cond");
        plot = optionList.getSwitch("plot");
        mesh = optionList.getSwitch("mesh");
        interpolation = optionList.getSwitch("interpolation");
    }
    else
        GISMO_ERROR("No XML file provided and no geometry file provided!");

    //! [Read XML file]

//    gsMatrix<> coefs;
//    for (index_t i = 0; i < mp.nPatches(); i++)
//    {
//        coefs = 0.25 * mp.patch(i).coefs();
//        mp.patch(i).setCoefs(coefs);
//    }
//    gsFileData<> fd;
//    fd << mp;
//    fd.save("ScaledGeometry");

    //! [Refinement]
    gsMultiBasis<real_t> dbasis(mp, false);//true: poly-splines (not NURBS)

    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    dbasis.setDegree( degree); // preserve smoothness
    //dbasis.degreeElevate(degree- mp.patch(0).degree(0));

    if (method == MethodFlags::DPATCH || method == MethodFlags::ALMOSTC1 || method == MethodFlags::SURFASG1)
        mp.degreeElevate(degree-mp.patch(0).degree(0));

    // h-refine each basis
    if (last)
    {
        for (int r =0; r < numRefine; ++r)
        {
            dbasis.uniformRefine(1, degree-smoothness);
            if (method == MethodFlags::DPATCH || method == MethodFlags::ALMOSTC1 || method == MethodFlags::SURFASG1)
                mp.uniformRefine(1, degree-smoothness);
        }
        numRefine = 0;
    }

    for (int r =0; r < numRefineIni; ++r)
    {
        dbasis.uniformRefine(1, degree-smoothness);
        if (method == MethodFlags::DPATCH || method == MethodFlags::ALMOSTC1 || method == MethodFlags::SURFASG1)
            mp.uniformRefine(1, degree-smoothness);

        beta /= 2;
    }
    // numRefine -= numRefineIni;

    // Assume that the condition holds for each patch TODO
    // Refine once
    // if (dbasis.basis(0).numElements() < 4)
    // {
    //     dbasis.uniformRefine(1, degree-smoothness);
    //     if (method == MethodFlags::DPATCH || method == MethodFlags::ALMOSTC1 || method == MethodFlags::SURFASG1)
    //         mp.uniformRefine(1, degree-smoothness);
    //     beta /= 2;
    // }

    for (size_t p=0; p!=mp.nPatches(); p++)
        for (index_t d=0; d!=2; d++)
            mp.patch(p).coefs().col(d) *= scaling;

    // Set the discretization space
    gsMappedBasis<2,real_t> bb2;
    // gsMappedSpline<2,real_t> geom;
    gsMultiPatch<> geom0, geom;
    geom = geom0 = mp;

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
    gsFunctionSet<>::Ptr solution;
    gsBiharmonicExprAssembler<real_t> BA;
    for (int r=0; r<=numRefine; ++r)
    {
        if (method == MethodFlags::APPROXC1)
        {
            dbasis.uniformRefine(1,degree -smoothness);
            meshsize[r] = dbasis.basis(0).getMinCellLength();

            // The approx. C1 space
            gsApproxC1Spline<2,real_t> approxC1(geom,dbasis);
            approxC1.options().setSwitch("info",info);
            approxC1.options().setSwitch("plot",plotApproxC1);
            approxC1.options().setSwitch("interpolation",interpolation);
            approxC1.options().setSwitch("second",second);
            approxC1.options().setInt("gluingDataDegree",gluingDataDegree);
            approxC1.options().setInt("gluingDataSmoothness",gluingDataSmoothness);
            approxC1.update(bb2);
        }
        else if (method == MethodFlags::NITSCHE)
        {
            dbasis.uniformRefine(1,degree-smoothness);
            meshsize[r] = dbasis.basis(0).getMinCellLength();
        }
        else if (method == MethodFlags::DPATCH)
        {
            mp.uniformRefine(1,degree-smoothness);
            dbasis = gsMultiBasis<>(mp);
            // dbasis.uniformRefine(1,degree-smoothness);
            // geom.uniformRefine(1,degree-smoothness);
            // dbasis = gsMultiBasis<>(geom);
            if (gsHTensorBasis<2,real_t> * test = dynamic_cast<gsHTensorBasis<2,real_t>*>(&dbasis.basis(0)))
                meshsize[r] = test->tensorLevel(0).getMinCellLength();
            else if (gsTensorBasis<2,real_t> * test = dynamic_cast<gsTensorBasis<2,real_t>*>(&dbasis.basis(0)))
                meshsize[r] = test->getMinCellLength();

            // Construct the D-Patch on mp
            // gsSparseMatrix<real_t> global2local;

            gsDPatch<2,real_t> dpatch(dbasis);
            // gsDPatch<2,real_t> dpatch(geom);
            dpatch.options().setInt("RefLevel",r);
            dpatch.options().setReal("Beta",beta);
            dpatch.options().setInt("Pi",PiMat);
            dpatch.options().setSwitch("SharpCorners",false);
            dpatch.compute();
            dpatch.update(bb2);

            // dpatch.matrix_into(global2local);
            // global2local = global2local.transpose();
            gsMultiBasis<> localbasis = dpatch.localBasis();
            // bb2.init(dbasis,global2local);

            if (r==0)
                geom0 = geom = dpatch.exportToPatches(mp);
            else
            {
                gsDofMapper mapper(localbasis);
                mapper.finalize();
                gsMatrix<> coefs;
                // First project the geometry geom0 onto bb2 and make a mapped spline
                gsInfo<<"L2-Projection error of geom0 on bb2 = "<<gsL2Projection<real_t>::projectGeometry(localbasis,bb2,geom0,coefs)<<"\n";
                coefs.resize(coefs.rows()/geom0.geoDim(),geom0.geoDim());
                gsMappedSpline<2,real_t> mspline;
                mspline.init(bb2,coefs);
                if (plot) gsWriteParaview( mspline, "mspline");

                // Then project onto localbasis so that geom represents the mapped geometry
                gsInfo<<"L2-Projection error of geom0 on dbasis = "<<gsL2Projection<real_t>::projectGeometry(localbasis,mspline,coefs)<<"\n";
                coefs.resize(coefs.rows()/mp.geoDim(),mp.geoDim());

                index_t offset = 0;
                for (index_t p = 0; p != geom.nPatches(); p++)
                {
                    geom.patch(p) = give(*localbasis.basis(p).makeGeometry((coefs.block(offset,0,mapper.patchSize(p),mp.geoDim()))));
                    offset += mapper.patchSize(p);
                }

            }
            if (plot) gsWriteParaview( geom, "geom",1000,true,false);

            // gsMatrix<> coefs;
            // gsL2Projection<real_t>::projectGeometry(dbasis,bb2,geom0,coefs);
        }
        else if (method == MethodFlags::ALMOSTC1)
        {
            gsMultiPatch<> tmp;
            // Project the previous geometry on the finest tensor level:
            for (size_t p=0; p!=geom.nPatches(); p++)
            {
                gsTHBSplineBasis<2,real_t> * thbsplineBasis;
                gsTensorBSplineBasis<2,real_t> * tbsplineBasis;
                if ((thbsplineBasis = dynamic_cast<gsTHBSplineBasis<2,real_t> *>(&geom.basis(p)) ))
                {
                    gsMatrix<> coefs;
                    // First project the geometry geom0 onto bb2 and make a mapped spline
                    gsTHBSplineBasis<2,real_t> tbasis = thbsplineBasis->tensorLevel(thbsplineBasis->maxLevel());
                    gsInfo<<"L2-Projection error of geom patch"<<p<<" on bb2 = "<<gsL2Projection<real_t>::projectGeometry(tbasis,geom.patch(p),coefs)<<"\n";
                    coefs.resize(coefs.rows()/geom0.geoDim(),geom0.geoDim());
                    tmp.addPatch(tbasis.makeGeometry(coefs));
                }
                else if ((tbsplineBasis = dynamic_cast<gsTensorBSplineBasis<2,real_t> *>(&geom.basis(p)) ))
                {
                    geom.patch(p).uniformRefine(1,degree-smoothness);
                    tmp.addPatch(geom.patch(p));
                }
                else
                    GISMO_ERROR("Geometry type not understood");
            }
            tmp.computeTopology();
            geom.swap(tmp);
            dbasis = gsMultiBasis<>(geom);

            if (plot) gsWriteParaview( tmp, "geom_ini",1000,true,false);

            if (gsHTensorBasis<2,real_t> * test = dynamic_cast<gsHTensorBasis<2,real_t>*>(&geom.basis(0)))
                meshsize[r] = test->tensorLevel(0).getMinCellLength();
            else if (gsTensorBasis<2,real_t> * test = dynamic_cast<gsTensorBasis<2,real_t>*>(&geom.basis(0)))
                meshsize[r] = test->getMinCellLength();

            // Construct the D-Patch on mp
            gsSparseMatrix<real_t> global2local;
            gsAlmostC1<2,real_t> almostC1(geom);
            almostC1.compute();
            almostC1.matrix_into(global2local);
            global2local = global2local.transpose();
            dbasis = almostC1.localBasis();
            bb2.init(dbasis,global2local);

            if (r==0)
            {
                geom0 = geom = almostC1.exportToPatches();
            }
            else
            {
            //     geom = almostC1.exportToPatches();
            // }
                gsMatrix<> targetCoefs, sourceCoefs;
                gsL2Projection<real_t>::projectGeometry(dbasis,bb2,geom0,targetCoefs);
                targetCoefs.resize(targetCoefs.rows()/2,2);
                bb2.getMapper().mapToSourceCoefs(targetCoefs,sourceCoefs);
                gsDofMapper mapper(dbasis);
                index_t offset = 0;
                for (index_t p = 0; p != geom0.nPatches(); p++)
                {
                    geom.patch(p) = give(*dbasis.basis(p).makeGeometry((sourceCoefs.block(offset,0,mapper.patchSize(p),mp.geoDim()))));
                    offset += mapper.patchSize(p);
                }
            }
        }
        else if (method == MethodFlags::SURFASG1) // Andrea
        {
            mp.uniformRefine(1,degree-smoothness);
            dbasis = gsMultiBasis<>(mp);

            meshsize[r] = dbasis.basis(0).getMinCellLength();

            gsC1SurfSpline<2,real_t> smoothC1(mp,dbasis);
            smoothC1.init();
            smoothC1.compute();

            gsSparseMatrix<real_t> global2local;
            global2local = smoothC1.getSystem();
            global2local = global2local.transpose();
            gsMultiBasis<> basis_temp;
            smoothC1.getMultiBasis(basis_temp);
            bb2.init(basis_temp,global2local);
        }
        gsInfo<< "." <<std::flush; // Approx C1 construction done
        bc.setGeoMap(geom);

        // Elements used for numerical integration
        std::vector<gsBasis<real_t> *> bases = bb2.getBasesCopy();
        gsMultiBasis<> intBasis;
        if (method != MethodFlags::NITSCHE)
            intBasis = gsMultiBasis<>(bases,mp.topology());
        else
            intBasis = dbasis;

        BA = gsBiharmonicExprAssembler<real_t>(geom,intBasis,f,bc);
        gsOptionList exprAssemblerOpts = Aopt.wrapIntoGroup("ExprAssembler");
        BA.setOptions(exprAssemblerOpts);
        if (method!=MethodFlags::NITSCHE)
            BA.setSpaceBasis(bb2);
        else
        {
            BA.options().setInt("Continuity",1);
            BA.options().setReal("PenaltyIfc",penalty_init);
        }

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

        // Compute the condition-number for the matrix (Slow)
        if (cond)
        {
#ifdef gsSpectra_ENABLED
            gsDebug<<"Computing the condition number using Spectra\n";
            real_t minev, maxev;
            index_t sz = BA.matrix().cols();
            gsStopwatch ev_time;
            gsSparseMatrix<real_t> I(BA.matrix().rows(),BA.matrix().cols());
            I.setIdentity();

            ev_time.restart();
            gsSpectraGenSymSolver<gsSparseMatrix<real_t>,Spectra::GEigsMode::Cholesky> evsolver_upp(BA.matrix(),I,1, sz);
            evsolver_upp.init();
            evsolver_upp.compute(Spectra::SortRule::LargestMagn,100,1e-10,Spectra::SortRule::LargestMagn);
            gsDebug<<"Largest eigenvalue computation finished in "<<ev_time.stop()<<" seconds\n";

            maxev = evsolver_upp.eigenvalues()(0);
            if (evsolver_upp.info()==Spectra::CompInfo::Successful)         { gsDebug<<"Spectra converged in "<<evsolver_upp.num_iterations()<<" iterations and with "<<evsolver_upp.num_operations()<<"operations. \n"; }
            else if (evsolver_upp.info()==Spectra::CompInfo::NumericalIssue){ GISMO_ERROR("Spectra did not converge! Error code: NumericalIssue"); }
            else if (evsolver_upp.info()==Spectra::CompInfo::NotConverging) { GISMO_ERROR("Spectra did not converge! Error code: NotConverging"); }
            else if (evsolver_upp.info()==Spectra::CompInfo::NotComputed)   { GISMO_ERROR("Spectra did not converge! Error code: NotComputed");   }
            else                                                      { GISMO_ERROR("No error code known"); }
            gsDebug << "Eigenvalues A*x=lambda*x:\n" << evsolver_upp.eigenvalues().transpose() <<"\n\n";

            ev_time.restart();
            gsSpectraGenSymShiftSolver<gsSparseMatrix<real_t>,Spectra::GEigsMode::ShiftInvert> evsolver_low(BA.matrix(),I,1, sz,0);
            evsolver_low.init();
            evsolver_low.compute(Spectra::SortRule::LargestMagn,100,1e-10,Spectra::SortRule::LargestMagn);
            gsDebug<<"Smallest eigenvalue computation finished in "<<ev_time.stop()<<" seconds\n";
            if (evsolver_low.info()==Spectra::CompInfo::Successful)         { gsDebug<<"Spectra converged in "<<evsolver_low.num_iterations()<<" iterations and with "<<evsolver_low.num_operations()<<"operations. \n"; }
            else if (evsolver_low.info()==Spectra::CompInfo::NumericalIssue){ GISMO_ERROR("Spectra did not converge! Error code: NumericalIssue"); }
            else if (evsolver_low.info()==Spectra::CompInfo::NotConverging) { GISMO_ERROR("Spectra did not converge! Error code: NotConverging"); }
            else if (evsolver_low.info()==Spectra::CompInfo::NotComputed)   { GISMO_ERROR("Spectra did not converge! Error code: NotComputed");   }
            else                                                      { GISMO_ERROR("No error code known"); }
            gsDebug << "Eigenvalues A*x=lambda*x:\n" << evsolver_low.eigenvalues().transpose() <<"\n\n";
            minev = evsolver_low.eigenvalues()(0);


            gsDebug << "Cond Number: " <<maxev/minev<< "\n";
            cond_num[r] = maxev/minev;
#else
            gsConjugateGradient<> cg(BA.matrix());

            cg.setCalcEigenvalues(true);
            //cg.setTolerance(1e-15);
            cg.setMaxIterations(100000);

            gsMatrix<> rhs, result;
            rhs.setRandom( BA.matrix().rows(), 1 );
            result.setRandom( BA.matrix().rows(), 1 );

            cg.solve(rhs,result);

            gsInfo << "Tol: " << cg.error() << "\n";
            gsInfo << "Max it: " << cg.iterations() << "\n";

            gsMatrix<real_t> eigenvalues;
            cg.getEigenvalues(eigenvalues);

            gsInfo << "Cond Number: " << eigenvalues.bottomRows(1)(0,0)/ eigenvalues(0,0) << "\n";
            cond_num[r] = eigenvalues.bottomRows(1)(0,0)/ eigenvalues(0,0);
#endif
        }
        else
            cond_num[r] = 0;

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
    if (cond)
        gsInfo<< "\nCondition-number: " << cond_num.transpose() << "\n";
    if (method == MethodFlags::NITSCHE)
        gsInfo<< "\nStabilization: " << penalty.transpose() << "\n";

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

        if (cond)
            gsInfo<<   "EoC (Cnum): "<< std::fixed<<std::setprecision(2)
                  <<( cond_num.tail(numRefine).array() /
                          cond_num.head(numRefine).array() ).log().transpose() / std::log(2.0) <<"\n";
    }
    //! [Error and convergence rates]

    //! [Export visualization in ParaView]
    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview( geom, "geom",1000,true);

        gsExprEvaluator<> ev;
        ev.setIntegrationElements(dbasis);
        auto u_sol = ev.getVariable(*solution);
        auto u_ex  = ev.getVariable(ms);
        auto G = ev.getMap(geom);

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
        index_t cols = method == MethodFlags::NITSCHE ? 7+penalty.cols() : 7;
        gsMatrix<real_t> error_collection(l2err.rows(), cols);
        error_collection.col(0) = meshsize;
        error_collection.col(1) = dofs;
        error_collection.col(2) = l2err;
        error_collection.col(3) = h1err;
        error_collection.col(4) = h2err;
        error_collection.col(5) = IFaceErr;
        error_collection.col(6) = cond_num;
        if (method == MethodFlags::NITSCHE)
            error_collection.block(0,7,penalty.rows(),penalty.cols()) = penalty;

        gsFileData<real_t> xml_out;
        xml_out << error_collection;
        xml_out.addString("Meshsize, dofs, l2err, h1err, h2err, iFaceErr, cond_num, (penalty)","Label");
        xml_out.addString(std::to_string(degree),"Degree");
        xml_out.addString(std::to_string(smoothness),"Regularity");
        xml_out.addString(std::to_string(numRefine),"NumRefine");
        xml_out.addString(std::to_string(method),"Method");
        // Add solution
        // [...]
        xml_out.save(output);
        gsInfo << "XML saved to " + output << "\n";
    }
    //! [Export data to xml]

    if (!write.empty())
    {
        std::ofstream file(write.c_str());
        file<<"Meshsize, dofs, l2err, h1err, h2err, iFaceErr, cond"<<"\n";
        for (index_t k=0; k<meshsize.size(); ++k)
        {
            file<<meshsize[k]<<","<<dofs[k]<<","<<l2err[k]<<","<<h1err[k]<<","<<h2err[k]<<","<<IFaceErr[k]<<cond_num[k]<<"\n";
        }
        file.close();
    }

    return EXIT_SUCCESS;
}// end main
