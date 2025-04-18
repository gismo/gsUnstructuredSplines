/** @file strongCoupling_example.cpp

    @brief Multi-patch shell analysis via unstructured splines

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst (2019 - ..., TU Delft)
               A. Farahat   (2019 - 2023, RICAM Linz)
*/

#include <gismo.h>

#include <gsUnstructuredSplines/src/gsMPBESBasis.h>
#include <gsUnstructuredSplines/src/gsMPBESSpline.h>
#include <gsUnstructuredSplines/src/gsDPatch.h>
#include <gsUnstructuredSplines/src/gsAlmostC1.h>
#include <gsUnstructuredSplines/src/gsApproxC1Spline.h>
#include <gsUnstructuredSplines/src/gsC1SurfSpline.h>

#ifdef gsKLShell_ENABLED
#include <gsKLShell/src/gsThinShellAssembler.h>
#include <gsKLShell/src/gsMaterialMatrixLinear.h>
#include <gsKLShell/src/gsFunctionSum.h>
#endif

using namespace gismo;

#ifdef gsKLShell_ENABLED
int main(int argc, char *argv[])
{
    bool plot       = false;
    bool plotGeo    = false;
    bool mesh       = false;
    bool stress     = false;
    bool write      = false;
    bool nonlinear  = false;
    bool homogeneous = false;

    index_t nmodes = 10;
    index_t mode   = 0;
    std::string input;

    real_t bcDirichlet = 1e3;
    real_t bcClamped = 1e3;

    real_t shift = 0.01;

    std::string geomFileName = "multipatches/2p_cornerSimo_AlmostC1_geom.xml";
    std::string basisFileName= "multipatches/2p_cornerSimo_AlmostC1_basis.xml";
    std::string bvpFileName  = "pde/shell/2p_cornerSimo_bvp.xml";
    std::string optFileName = "options/solver_options.xml";
    std::string out = "ModalResults";

    gsCmdLine cmd("Composite basis tests.");
    cmd.addReal( "D", "Dir", "Dirichlet BC penalty scalar",  bcDirichlet );
    cmd.addReal( "C", "Cla", "Clamped BC penalty scalar",  bcClamped );
    cmd.addString( "G", "geom","File containing the geometry",  geomFileName );
    cmd.addString( "b", "bas", "File containing the basis (dbasis and global2local)",  basisFileName );
    cmd.addString( "B", "bvp", "File containing the Boundary Value Problem (BVP)",  bvpFileName );
    cmd.addString( "O", "opt", "File containing solver options",  optFileName );
    cmd.addString( "o", "out", "Output directory",  out );
    cmd.addInt( "N", "nmodes", "Number of modes", nmodes );
    cmd.addInt( "M", "mode", "Mode number", mode );
    cmd.addReal( "S", "shift", "Set the shift of the solver.",  shift );
    cmd.addSwitch("plot", "plot",plot);
    cmd.addSwitch("plotGeo", "plotGeo",plotGeo);
    cmd.addSwitch("mesh", "mesh",mesh);
    cmd.addSwitch("stress", "stress",stress);
    cmd.addSwitch("write", "write",write);
    cmd.addSwitch( "nl", "Print information", nonlinear );
    cmd.addSwitch("homogeneous", "homogeneous dirichlet BCs",homogeneous);
    // to do:
    // smoothing method add nitsche @Pascal

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    GISMO_ENSURE(!(geomFileName.empty()) && !(basisFileName.empty()) && !(bvpFileName.empty()),"Not all filenames have been provided.\n"
                        <<"geomFileName  = "<<geomFileName<<"\n"
                        <<"basisFileName = "<<basisFileName<<"\n"
                        <<"bvpFileName   = "<<bvpFileName);

    gsFileData<> fd;

    gsMultiPatch<> geom;
    gsSparseMatrix<> global2local;
    gsMultiBasis<> dbasis;

    gsMappedBasis<2,real_t> bb2;
    gsBoundaryConditions<> bc;

    gsInfo<<"Reading geometry from "<<geomFileName<<"..."<<std::flush;
    gsReadFile<>(geomFileName, geom);
    if (geom.nInterfaces()==0 && geom.nBoundary()==0)
    {
        gsInfo<<"No topology found. Computing it...";
        geom.computeTopology();
    }
    gsInfo<<"Finished\n";

    gsInfo<<"Reading mapped basis from "<<basisFileName<<"..."<<std::flush;
    fd.read(basisFileName);
    GISMO_ENSURE(fd.template getFirst<gsMultiBasis<>>(dbasis)        ,"No multibasis was read");
    GISMO_ENSURE(fd.template getFirst<gsSparseMatrix<>>(global2local),"No sparse matrix was read");
    gsInfo<<"Finished\n";
    bb2.init(dbasis,global2local);

    fd.read(bvpFileName);
    if (homogeneous)
    {
        for (gsMultiPatch<>::const_biterator bit = geom.bBegin(); bit != geom.bEnd(); ++bit)
        {
            bc.addCondition(*bit, condition_type::dirichlet, 0, 0, false, 0);
            bc.addCondition(*bit, condition_type::dirichlet, 0, 0, false, 1);
            bc.addCondition(*bit, condition_type::dirichlet, 0, 0, false, 2);

            bc.addCondition(*bit, condition_type::weak_clamped, 0, 0, false, 0);
            bc.addCondition(*bit, condition_type::weak_clamped, 0, 0, false, 1);
            bc.addCondition(*bit, condition_type::weak_clamped, 0, 0, false, 2);
        }
    }
    else
    {
        index_t num = 0;
        gsInfo<<"Reading BCs from "<<bvpFileName<<"...";
        num = fd.template count<gsBoundaryConditions<>>();
        GISMO_ENSURE(num==1,"Number of boundary condition objects in XML should be 1, but is "<<num);
        fd.template getFirst<gsBoundaryConditions<>>(bc); // Multipatch domain
        gsInfo<<"Finished\n";

    }

    gsInfo<<bc<<"\n";
    bc.setGeoMap(geom);

    // Loads
    gsFunctionExpr<> force, pressure;
    gsInfo<<"Reading force function from "<<bvpFileName<<" (ID=21) ...";
    fd.getId(21, force); // id=1: source function
    gsInfo<<"Finished\n";
    // fd.getId(22, pressure); // id=1: source function ------- TO DO!
    // gsInfo<<"Pressure function "<< force << "\n";

    // Loads
    gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();
    gsMatrix<> points,loads;
    gsMatrix<index_t> pid_ploads;
    if ( fd.hasId(30) )
        fd.getId(30,points);
    if ( fd.hasId(31) )
        fd.getId(31,loads);

    if ( fd.hasId(32) )
        fd.getId(32,pid_ploads);
    else
        pid_ploads = gsMatrix<index_t>::Zero(1,points.cols());

    for (index_t k =0; k!=points.cols(); k++)
        pLoads.addLoad(points.col(k), loads.col(k), pid_ploads.at(k) ); // in parametric domain!

    gsInfo<<pLoads;

    // Reference points
    gsMatrix<index_t> refPatches;
    gsMatrix<> refPoints, refPars, refValue; // todo: add refValue..
    gsInfo<<"Reading reference point locations from "<<bvpFileName<<" (ID=50) ...";
    if ( fd.hasId(50) )
        fd.getId(50,refPoints);
    if (refPoints.rows()==2)
    {
        refPars = refPoints;
        gsInfo<<"Reference points are provided in parametric coordinates.\n";
    }
    else if (refPoints.rows()==3)
        gsInfo<<"Reference points are provided in physical coordinates.\n";
    else
        gsInfo<<"No reference points are provided.\n";
    gsInfo<<"Finished\n";
    gsInfo<<"Reading reference patches from "<<bvpFileName<<" (ID=51) ...";
    if ( fd.hasId(51) )
        fd.getId(51,refPatches);
    gsInfo<<"Finished\n";
    gsInfo<<"Reading reference values from "<<bvpFileName<<" (ID=52) ...";
    if ( fd.hasId(52) )
        fd.getId(52,refValue);
    else
        refValue = gsMatrix<>::Zero(geom.geoDim(),refPoints.cols());
    gsInfo<<"Finished\n";
    GISMO_ENSURE(refPatches.cols()==refPoints.cols(),"Number of reference points and patches do not match");


    // Material properties
    gsFunctionExpr<> t,E,nu,rho;
    gsInfo<<"\tReading thickness from "<<bvpFileName<<" (ID=10) ..."<<std::flush;
    fd.getId(10,t);
    gsInfo<<"Finished\n";

    gsInfo<<"\tReading Young's Modulus from "<<bvpFileName<<" (ID=11) ..."<<std::flush;
    fd.getId(11,E);
    gsInfo<<"Finished\n";

    gsInfo<<"\tReading Poisson ratio from "<<bvpFileName<<" (ID=12) ..."<<std::flush;
    fd.getId(12,nu);
    gsInfo<<"Finished\n";

    gsInfo<<"\tReading density from "<<bvpFileName<<" (ID=13) ..."<<std::flush;
    fd.getId(13,rho);
    gsInfo<<"Finished\n";
    gsInfo<<"Finished\n";

    if (plot || write)
        gsFileManager::mkdir(out);

    if (plotGeo) gsWriteParaview(geom,"geom",1000,true,false);

    std::vector<gsFunctionSet<>*> parameters(2);
    parameters[0] = &E;
    parameters[1] = &nu;

    gsMaterialMatrixLinear<3,real_t> materialMatrix(geom,t,parameters,rho);

    gsThinShellAssembler<3, real_t, true> assembler;

    //! [Solver loop]
    // #ifdef GISMO_WITH_PARDISO
    // gsSparseSolver<>::PardisoLU solver;
    // #else
    gsSparseSolver<>::CGDiagonal solver;
    // #endif
    gsVector<> solVector;

    assembler = gsThinShellAssembler<3, real_t, true>(geom,dbasis,bc,force,&materialMatrix);
    assembler.options().setInt("Continuity",-1);
    assembler.options().setReal("WeakDirichlet",bcDirichlet);
    assembler.options().setReal("WeakClamped",bcClamped);
    assembler.setSpaceBasis(bb2);
    assembler.setPointLoads(pLoads);
    // gsOptionList options = assembler.options();
    // options.setInt("Continuity",1);
    // assembler.setOptions(options);

    // Initialize the system

    gsInfo<<"Assembling stiffness matrix and rhs vector..."<<std::flush;
    assembler.assemble();
    gsSparseMatrix<> matrix = assembler.matrix();
    gsVector<> vector = assembler.rhs();
    gsInfo<<"Finished\n";

    // Nonlinear
    // Function for the Jacobian
    typedef std::function<gsSparseMatrix<real_t> (gsVector<real_t> const &)>    Jacobian_t;
    typedef std::function<gsVector<real_t> (gsVector<real_t> const &) >         Residual_t;
    Jacobian_t Jacobian = [&geom,&bb2,&assembler](gsVector<real_t> const &x)
    {
        gsMatrix<real_t> solFull = assembler.fullSolutionVector(x);
        GISMO_ASSERT(solFull.rows() % 3==0,"Rows of the solution vector does not match the number of control points");
        solFull.resize(solFull.rows()/3,3);
        gsMappedSpline<2,real_t> mspline(bb2,solFull);
        gsFunctionSum<real_t> def(&geom,&mspline);

        assembler.assembleMatrix(def);
        // gsSparseMatrix<real_t> m =
        return assembler.matrix();
    };
    // Function for the Residual
    Residual_t Residual = [&geom,&bb2,&assembler](gsVector<real_t> const &x)
    {
        gsMatrix<real_t> solFull = assembler.fullSolutionVector(x);
        GISMO_ASSERT(solFull.rows() % 3==0,"Rows of the solution vector does not match the number of control points");
        solFull.resize(solFull.rows()/3,3);

        gsMappedSpline<2,real_t> mspline(bb2,solFull);
        gsFunctionSum<real_t> def(&geom,&mspline);

        assembler.assembleVector(def);
        return assembler.rhs();
    };

    gsInfo<<"Solving system"<<std::flush;
    // Linear solve
    solver.compute( matrix );
    solVector = solver.solve(vector);
    gsInfo<<"Finished.\n";
    if (nonlinear)
    {
        real_t residual = vector.norm();
        real_t residual0 = residual;
        real_t residualOld = residual;
        gsVector<real_t> updateVector = solVector;
        gsVector<real_t> resVec = Residual(solVector);
        gsSparseMatrix<real_t> jacMat;
        for (index_t it = 0; it != 100; ++it)
        {
            jacMat = Jacobian(solVector);
            solver.compute(jacMat);
            updateVector = solver.solve(resVec); // this is the UPDATE
            solVector += updateVector;

            resVec = Residual(solVector);
            residual = resVec.norm();

            gsInfo<<"Iteration: "<< it
               <<", residue: "<< residual
               <<", update norm: "<<updateVector.norm()
               <<", log(Ri/R0): "<< math::log10(residualOld/residual0)
               <<", log(Ri+1/R0): "<< math::log10(residual/residual0)
               <<"\n";

            residualOld = residual;

            if (updateVector.norm() < 1e-6)
                break;
            else if (it+1 == it)
                gsWarn<<"Maximum iterations reached!\n";
        }
    }
    // ! [Solver loop]

    /// Make a gsMappedSpline to represent the solution
    // 1. Get all the coefficients (including the ones from the eliminated BCs.)
    gsMatrix<real_t> solFull = assembler.fullSolutionVector(solVector);
    gsMatrix<real_t> solZero = solFull;
    solZero.setZero();

    // 2. Reshape all the coefficients to a Nx3 matrix
    GISMO_ASSERT(solFull.rows() % 3==0,"Rows of the solution vector does not match the number of control points");
    solZero.resize(solZero.rows()/3,3);
    solFull.resize(solFull.rows()/3,3);

    // 3. Make the mapped spline
    gsMappedSpline<2,real_t> mspline(bb2,solFull);

    gsField<> solField(geom, mspline,true);

    if (refPoints.cols()!=0)
    {
        gsMatrix<> result;
        if (refPoints.rows()==3) // then they are provided in the physical domain and should be mapped to the parametric domain
        {
            refPars.resize(2,refPoints.cols());
            for (index_t p = 0; p!=refPoints.cols(); p++)
            {
                geom.patch(refPatches(0,p)).invertPoints(refPoints.col(p),result,1e-10);
                if (result.at(0)==std::numeric_limits<real_t>::infinity()) // if failed
                    gsWarn<<"Point inversion failed\n";
                refPars.col(p) = result;
            }
        }

        gsInfo<<"Physical coordinates of points\n";
        for (index_t p=0; p!=refPars.cols(); p++)
        {
            gsInfo<<",x"<<std::to_string(p)<<",y"<<std::to_string(p)<<",z"<<std::to_string(p);
        }
        gsInfo<<"\n";

        for (index_t p=0; p!=refPars.cols(); ++p)
        {
            geom.patch(refPatches(0,p)).eval_into(refPars.col(p),result);
            gsInfo<<result.row(0)<<","<<result.row(1)<<","<<result.row(2)<<",";
        }
        gsInfo<<"\n";
        gsMatrix<> refs(1,geom.geoDim()*refPoints.cols());
        for (index_t p=0; p!=refPars.cols(); p++)
            refs.block(0,p*geom.geoDim(),1,geom.geoDim()) = geom.piece(refPatches(0,p)).eval(refPars.col(p)).transpose();
        gsInfo<<"Reference point coordinates\n";
        for (index_t p=0; p!=refPars.cols(); ++p)
            gsInfo<<"x"<<std::to_string(p)<<"\ty"<<std::to_string(p)<<"\tz"<<std::to_string(p)<<"\t";
        gsInfo<<"\n";
        for (index_t p=0; p!=refPars.cols(); ++p)
            gsInfo<<refs(0,geom.geoDim()*p)<<"\t"<<refs(0,geom.geoDim()*p+1)<<"\t"<<refs(0,geom.geoDim()*p+2)<<"\t";
        gsInfo<<"\n";

        for (index_t p=0; p!=refPars.cols(); p++)
            refs.block(0,p*geom.geoDim(),1,geom.geoDim()) = solField.value(refPars.col(p),refPatches(0,p)).transpose();
        gsInfo<<"Computed values\n";
        for (index_t p=0; p!=refPars.cols(); ++p)
            gsInfo<<"x"<<std::to_string(p)<<"\ty"<<std::to_string(p)<<"\tz"<<std::to_string(p)<<"\t";
        gsInfo<<"\n";
        for (index_t p=0; p!=refPars.cols(); ++p)
            gsInfo<<refs(0,geom.geoDim()*p)<<"\t"<<refs(0,geom.geoDim()*p+1)<<"\t"<<refs(0,geom.geoDim()*p+2)<<"\t";
        gsInfo<<"\n";

        gsInfo<<"Reference values\n"; // provided as geom.geoDim() x points.cols() matrix
        for (index_t p=0; p!=refValue.cols(); ++p)
            gsInfo<<"x"<<std::to_string(p)<<"\ty"<<std::to_string(p)<<"\tz"<<std::to_string(p)<<"\t";
        gsInfo<<"\n";
        for (index_t p=0; p!=refValue.cols(); ++p)
            for (index_t d=0; d!=geom.geoDim(); d++)
                gsInfo<<refValue(d,p)<<"\t";
        gsInfo<<"\n";
    }

    gsInfo << "Number of Dofs: " << assembler.numDofs() << "\n";

    //! [Export visualization in ParaView]
    if (plot)
    {
        // 4. Plot the mapped spline on the original geometry
        gsField<> solField(geom, mspline,true);
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview<>( solField, "Deformation", 1000, mesh);

        // // 4. Plot the mapped spline on the original geometry
        // gsField<> solField2(mp_def, def,true);
        // gsInfo<<"Plotting in Paraview...\n";
        // gsWriteParaview<>( solField2, "Deformation_", 1000, true);


        // // 5. Plot the mapped spline on the deformed geometry
        // gsField<> defField(def, def,true);
        // gsInfo<<"Plotting in Paraview...\n";
        // gsWriteParaview<>( defField, "mp_def", 1000, true);

        // gsMultiPatch<> mpatches = mbasis.exportToPatches(tmp);

        // gsField<> solfield(geom,def,true);
        // gsWriteParaview(solfield,"solfield",1000,true);
    }

    gsMappedSpline<2,real_t> mspline_ori(bb2,solZero);
    gsMappedSpline<2,real_t> mspline_def(bb2,solFull);

    gsFunctionSum<real_t> ori(&geom,&mspline_ori);
    gsFunctionSum<real_t> def(&geom,&mspline_def);

    if (stress)
    {
        gsPiecewiseFunction<> membraneStresses;
        gsDebugVar("MembraneStress construction");
        assembler.constructStress(ori,def,membraneStresses,stress_type::membrane);
        gsWriteParaview(ori,membraneStresses,"MembraneStress",5000);

        gsPiecewiseFunction<> membraneStressesVM;
        gsDebugVar("MembraneStress (VM) construction");
        assembler.constructStress(ori,def,membraneStressesVM,stress_type::von_mises_membrane);
        gsWriteParaview(ori,membraneStressesVM,"MembraneStressVM",5000);

        gsPiecewiseFunction<> flexuralStresses;
        gsDebugVar("FlexuralStress construction");
        assembler.constructStress(ori,def,flexuralStresses,stress_type::flexural);
        gsWriteParaview(geom,flexuralStresses,"FlexuralStress",5000);
    }
    // if (write)
    // {
    //     std::vector<std::string> headers = {"u_x","u_y","u_z"};
    //     gsStaticOutput<real_t> ptsWriter("pointcoordinates.csv",refPoints);
    //     ptsWriter.init(headers);
    //     gsMatrix<> pointResults(geom.geoDim(),refPoints.cols());
    //     for (index_t p=0; p!=refPoints.cols(); p++)
    //         pointResults.col(p) = geom.piece(refPatches(0,p)).eval(refPoints.col(p));
    //     ptsWriter.add(pointResults);

    //     gsStaticOutput<real_t> numWriter("numerical.csv",refPoints);
    //     numWriter.init(headers);
    //     for (index_t p=0; p!=refPoints.cols(); p++)
    //         pointResults.col(p) = solField.value(refPoints.col(p),refPatches(0,p));
    //     numWriter.add(pointResults);

    //     gsStaticOutput<real_t> refWriter("reference.csv",refPoints);
    //     refWriter.init(headers);
    //     refWriter.add(refValue);
    // }
    //! [Export visualization in ParaView]

    return EXIT_SUCCESS;
}
#else
int main(int argc, char *argv[])
{
    GISMO_ERROR("G+Smo is not compiled with the gsKLShell module.");
    return EXIT_FAILURE;
}
#endif
