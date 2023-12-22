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

#include <gsUnstructuredSplines/src/gsApproxC1Spline.h>
#include <gsUnstructuredSplines/src/gsDPatch.h>
#include <gsUnstructuredSplines/src/gsAlmostC1.h>
#include <gsUnstructuredSplines/src/gsC1SurfSpline.h>

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

void computeStabilityParameter(gsMultiPatch<> mp, gsMultiBasis<> dbasis, gsMatrix<real_t> & mu_interfaces)
{
    mu_interfaces.setZero();

    index_t i = 0;
    for ( typename gsMultiPatch<real_t>::const_iiterator it = mp.iBegin(); it != mp.iEnd(); ++it, ++i)
    {
        gsMultiPatch<> mp_temp;
        mp_temp.addPatch(mp.patch(it->first().patch));
        mp_temp.addPatch(mp.patch(it->second().patch));
        mp_temp.computeTopology();

        gsMultiBasis<> dbasis_temp;
        dbasis_temp.addBasis(dbasis.basis(it->first().patch).clone().release());
        dbasis_temp.addBasis(dbasis.basis(it->second().patch).clone().release());

        gsBoundaryConditions<> bc;

//        patchSide pS1 = mp_temp.interfaces()[0].first();
//        patchSide pS2 = mp_temp.interfaces()[0].second();
//
//
//        index_t side = pS1.index() < 3 ? (pS1.index() == 1 ? 2 : 1) : (pS1.index() == 3 ? 4 : 3);
//        bc.addCondition(patchSide(pS1.patchIndex(), side), condition_type::dirichlet, 0);
//
//        side = pS2.index() < 3 ? (pS2.index() == 1 ? 2 : 1) : (pS2.index() == 3 ? 4 : 3);
//        bc.addCondition(patchSide(pS2.patchIndex(), side), condition_type::dirichlet, 0);

        // Make the Eigenvalue problem to a homogeneous one
        for (gsMultiPatch<>::const_biterator bit = mp_temp.bBegin(); bit != mp_temp.bEnd(); ++bit)
            bc.addCondition(*bit, condition_type::dirichlet, 0);

        gsExprAssembler<real_t> A2(1, 1), B2(1, 1);

        // Elements used for numerical integration
        A2.setIntegrationElements(dbasis_temp);
        B2.setIntegrationElements(dbasis_temp);

        // Set the geometry map
        auto GA = A2.getMap(mp_temp);
        auto GB = B2.getMap(mp_temp);

        // Set the discretization space
        auto uA = A2.getSpace(dbasis_temp);
        auto uB = B2.getSpace(dbasis_temp);

        uA.setup(bc, dirichlet::homogeneous, 0);
        uB.setup(bc, dirichlet::homogeneous,0);
        //uA.setup(0);
        //uB.setup(0);

        A2.initSystem();
        B2.initSystem();

        real_t c = 0.25;
        A2.assembleIfc(mp_temp.interfaces(),
                       c * ilapl(uA.left(), GA.left()) * ilapl(uA.left(), GA.left()).tr() * nv(GA.left()).norm(),
                       c * ilapl(uA.left(), GA.left()) * ilapl(uA.right(), GA.right()).tr() * nv(GA.left()).norm(),
                       c * ilapl(uA.right(), GA.right()) * ilapl(uA.left(), GA.left()).tr() * nv(GA.left()).norm(),
                       c * ilapl(uA.right(), GA.right()) * ilapl(uA.right(), GA.right()).tr() * nv(GA.left()).norm());

        B2.assemble(ilapl(uB, GB) * ilapl(uB, GB).tr() * meas(GB));

        // TODO INSTABLE && SLOW
        gsMatrix<> AA = A2.matrix().toDense().cast<real_t>();
        gsMatrix<> BB = B2.matrix().toDense().cast<real_t>();
        gsEigen::GeneralizedSelfAdjointEigenSolver<gsEigen::MatrixXd> ges(AA, BB);

        real_t m_h      = dbasis_temp.basis(0).getMinCellLength(); //*dbasis.basis(0).getMinCellLength();
        mu_interfaces(i,0) = 16.0 * m_h * ges.eigenvalues().array().maxCoeff();
/*
        gsSparseSolver<>::SimplicialLDLT sol;
        sol.compute(B2.matrix());
        gsSparseMatrix<> R = sol.matrixU();
        gsSparseMatrix<> RT = sol.matrixL();
        gsMatrix<> AAA = RT.toDense().inverse() * AA * R.toDense().inverse();

        gsConjugateGradient<> cg(AAA);

        cg.setCalcEigenvalues(true);
        cg.setTolerance(1e-15);
        cg.setMaxIterations(100000);

        gsMatrix<> rhs, result;
        rhs.setRandom( AAA.rows(), 1 );
        result.setRandom( AAA.rows(), 1 );

        cg.solve(rhs,result);

        gsInfo << "Tol: " << cg.error() << "\n";
        gsInfo << "Max it: " << cg.iterations() << "\n";

        gsMatrix<real_t> eigenvalues;
        cg.getEigenvalues(eigenvalues);

        gsInfo << "Cond Number: " << eigenvalues.bottomRows(1)(0,0) << "\n";
*/
    }
}

void writeToCSVfile(std::string name, gsMatrix<> matrix)
{
    std::ofstream file(name.c_str());
    for(int  i = 0; i < matrix.rows(); i++){
        for(int j = 0; j < matrix.cols(); j++){
           std::string str = std::to_string(matrix(i,j));
           if(j+1 == matrix.cols()){
               file<<std::setprecision(10)<<str;
           }else{
               file<<std::setprecision(10)<<str<<',';
           }
        }
        file<<'\n';
    }
}

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot  = false;
    bool plotApproxC1 = false;
    bool mesh  = false;
    bool first      = false;
    bool write      = false;

    index_t method = 0;

    index_t numRefine  = 3;
    index_t degree = 3;
    index_t smoothness = 2;

    index_t gluingDataDegree = -1;
    index_t gluingDataSmoothness = -1;

    bool last = false;
    bool info = false;
    bool cond = false;
    bool interpolation = false;

    real_t penalty_init = -1.0;
    std::string output;

    std::string fn = "planar/1p_square.xml";

    gsCmdLine cmd("Tutorial on solving a Biharmonic problem with different spaces.");
    // Flags related to the method (default: Approx C1 method)
    cmd.addInt( "m", "method", "The chosen method for the biharmonic problem", method );

    // Flags related to the input/geometry
    cmd.addString( "f", "file", "Input geometry file from path (with .xml)", fn );

    // Flags related to the discrete settings
    cmd.addInt( "p", "degree", "Set the polynomial degree of the basis.", degree );
    cmd.addInt( "s", "smoothness", "Set the smoothness of the basis.",  smoothness );
    cmd.addInt( "r", "numRefine", "Number of refinement-loops.",  numRefine );

    // Flags related to the approximate C1 method
    cmd.addInt( "P", "gluingDataDegree","Set the polynomial degree for the gluing data", gluingDataDegree );
    cmd.addInt( "R", "gluingDataSmoothness", "Set the smoothness for the gluing data",  gluingDataSmoothness );
    cmd.addSwitch("interpolation", "Compute the basis constructions with interpolation", interpolation);
    cmd.addSwitch("info", "Getting the information inside of Approximate C1 basis functions", info);
    cmd.addSwitch("plotApproxC1", "Plot the approximate C1 basis functions", plotApproxC1);

    // Flags related to Nitsche's method
    cmd.addReal( "y", "penalty", "Fixed Penalty value for Nitsche's method",  penalty_init);

    // Flags related to the output
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("write", "Write convergence data to file", write);
    cmd.addSwitch("mesh", "Plot the mesh", mesh);
    cmd.addSwitch("first", "Plot only first", first);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Initialize data]
    gsMultiPatch<real_t> mp;
    gsBoundaryConditions<> bc;
    gsOptionList optionList;
    //! [Initialize data]

    real_t thickness = 0.01;
    real_t E_modulus = 1e5;
    real_t Density = 1e5;
    real_t PoissonRatio = 0.3;
    real_t D = E_modulus*math::pow(thickness,3)/(12*(1-math::pow(PoissonRatio,2)));

    //! [Read Argument inputs]
    //! [Read geometry]
    std::string string_geo;
    GISMO_ENSURE(!fn.empty(),"No XML file provided and no geometry file provided!");

    gsInfo << "Filedata: " << fn << "\n";
    gsReadFile<>(fn, mp);
    mp.clearTopology();
    mp.computeTopology();
    gsMultiBasis<real_t> dbasis(mp, false);//true: poly-splines (not NURBS)
    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    dbasis.setDegree( degree); // preserve smoothness

    if (method == MethodFlags::DPATCH || method == MethodFlags::ALMOSTC1 || method == MethodFlags::SURFASG1)
        mp.degreeElevate(degree-mp.patch(0).degree(0));

    for (int r =0; r < numRefine; ++r)
    {
        dbasis.uniformRefine(1, degree-smoothness);
        if (method == MethodFlags::DPATCH || method == MethodFlags::ALMOSTC1 || method == MethodFlags::SURFASG1)
            mp.uniformRefine(1, degree-smoothness);
    }

    // // Assume that the condition holds for each patch TODO
    // // Refine once
    // if (mp.basis(0).numElements() < 4)
    // {
    //     dbasis.uniformRefine(1, degree-smoothness);
    //     if (method == MethodFlags::DPATCH || method == MethodFlags::ALMOSTC1 || method == MethodFlags::SURFASG1)
    //         mp.uniformRefine(1, degree-smoothness);
    // }

    //! [Refinement]

    for (size_t p=0; p!=mp.nPatches(); p++)
        gsInfo<<"Basis "<<p<<": "<<mp.basis(p)<<"\n";

    // //! [Boundary condition]
    for (gsMultiPatch<>::const_biterator bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
        bc.addCondition(*bit, condition_type::dirichlet, 0);

    bc.setGeoMap(mp);
    gsInfo << "Boundary conditions:\n" << bc << "\n";
    //! [Boundary condition]

    optionList = cmd;
    gsInfo << "OptionList: " << optionList << "\n";
    gsInfo << "Finished\n";
    //! [Read Argument inputs]


    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
#ifdef _OPENMP
    gsInfo<< "Available threads: "<< omp_get_max_threads() <<"\n";
#endif

//    gsWriteParaview(mp, "geom", 2000);
//
//    gsVector<> vec;
//    vec.setLinSpaced(5,0,1);
//    gsMatrix<> points;
//    points.setZero(2,5);
//    points.row(1) = vec;
//    gsInfo << mp.patch(0).eval(points) << "\n";
//
//    points.setOnes(2,5);
//    points.row(1) = vec;
//    gsInfo << mp.patch(1).eval(points) << "\n";
//
//    mp.patch(0).degreeElevate(2);
//    mp.patch(1).degreeElevate(1);
//    //mp.patch(0).uniformRefine(1);
//
//    gsFileData<> fd;
//    fd << mp;
//    fd.save("geometry");
    //! [Refinement]

    //! [Problem setup]

    gsMappedBasis<2,real_t> bb2;
    if (method == MethodFlags::APPROXC1)
    {
        // The approx. C1 space
        gsApproxC1Spline<2,real_t> approxC1(mp,dbasis);
        approxC1.options().setSwitch("info",info);
        approxC1.options().setSwitch("plot",plotApproxC1);
        approxC1.options().setSwitch("interpolation",interpolation);
        approxC1.options().setInt("gluingDataDegree",gluingDataDegree);
        approxC1.options().setInt("gluingDataSmoothness",gluingDataSmoothness);
        approxC1.update(bb2);
    }
    else if (method == MethodFlags::DPATCH)
    {
        gsSparseMatrix<real_t> global2local;
        gsDPatch<2,real_t> dpatch(mp);
        dpatch.options().setInt("Pi",0);
        dpatch.options().setSwitch("SharpCorners",false);
        dpatch.compute();
        dpatch.matrix_into(global2local);
        global2local = global2local.transpose();
        mp = dpatch.exportToPatches();
        dbasis = dpatch.localBasis();
        bb2.init(dbasis,global2local);
    }
    else if (method == MethodFlags::ALMOSTC1)
    {
        gsSparseMatrix<real_t> global2local;
        gsAlmostC1<2,real_t> almostC1(mp);
        almostC1.compute();
        almostC1.matrix_into(global2local);
        global2local = global2local.transpose();
        mp = almostC1.exportToPatches();
        dbasis = almostC1.localBasis();
        bb2.init(dbasis,global2local);
    }
    else if (method == MethodFlags::SURFASG1) // Andrea
    {
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


    //! [Problem setup]
    gsExprAssembler<real_t> A(1,1);
    //gsInfo<<"Active options:\n"<< A.options() <<"\n";

    // Elements used for numerical integration
    A.setIntegrationElements(dbasis);
    gsExprEvaluator<real_t> ev(A);

    // Set the geometry map
    auto G = A.getMap(mp);

    // Set the discretization space
    auto u = method == MethodFlags::NITSCHE ? A.getSpace(dbasis) : A.getSpace(bb2);

    // Solution vector and solution variable
    gsMatrix<real_t> solVector;
    auto u_sol = A.getSolution(u, solVector);

    if (method == MethodFlags::NITSCHE)
        u.setup(bc,dirichlet::l2Projection,0);
    else
        u.setup(bc,dirichlet::l2Projection,-1);

    // Initialize the system
    A.initSystem();
    gsInfo<< A.numDofs() <<std::flush;

    // Compute the system matrix and right-hand side
    A.assemble(D * ilapl(u, G) * ilapl(u, G).tr() * meas(G));

    if (method == MethodFlags::NITSCHE)
    {
        // For Nitsche
        gsMatrix<real_t> mu_interfaces(mp.nInterfaces(),1);
        if (penalty_init == -1.0)
            computeStabilityParameter(mp, dbasis, mu_interfaces);

        index_t i = 0;
        for ( typename gsMultiPatch<real_t>::const_iiterator it = mp.iBegin(); it != mp.iEnd(); ++it, ++i)
        {
            real_t stab     = 4 * ( dbasis.maxCwiseDegree() + dbasis.dim() ) * ( dbasis.maxCwiseDegree() + 1 );
            real_t m_h      = dbasis.basis(0).getMinCellLength(); //*dbasis.basis(0).getMinCellLength();
            real_t mu       = 2 * stab / m_h;
            real_t alpha = 1;

            //mu = penalty_init == -1.0 ? mu : penalty_init / m_h;
            if (penalty_init == -1.0)
                mu = mu_interfaces(i,0) / m_h;
            else
                mu = penalty_init / m_h;

            std::vector<boundaryInterface> iFace;
            iFace.push_back(*it);
            A.assembleIfc(iFace,
                    //B11
                          -alpha * 0.5 * igrad(u.left(), G) * nv(G.left()).normalized() *
                          (ilapl(u.left(), G)).tr() * nv(G.left()).norm(),
                          -alpha * 0.5 *
                          (igrad(u.left(), G) * nv(G.left()).normalized() * (ilapl(u.left(), G)).tr()).tr() *
                          nv(G.left()).norm(),
                    //B12
                          -alpha * 0.5 * igrad(u.left(), G.left()) * nv(G.left()).normalized() *
                          (ilapl(u.right(), G.right())).tr() * nv(G.left()).norm(),
                          -alpha * 0.5 * (igrad(u.left(), G.left()) * nv(G.left()).normalized() *
                                          (ilapl(u.right(), G.right())).tr()).tr() * nv(G.left()).norm(),
                    //B21
                          alpha * 0.5 * igrad(u.right(), G.right()) * nv(G.left()).normalized() *
                          (ilapl(u.left(), G.left())).tr() * nv(G.left()).norm(),
                          alpha * 0.5 * (igrad(u.right(), G.right()) * nv(G.left()).normalized() *
                                         (ilapl(u.left(), G.left())).tr()).tr() * nv(G.left()).norm(),
                    //B22
                          alpha * 0.5 * igrad(u.right(), G.right()) * nv(G.left()).normalized() *
                          (ilapl(u.right(), G.right())).tr() * nv(G.left()).norm(),
                          alpha * 0.5 * (igrad(u.right(), G.right()) * nv(G.left()).normalized() *
                                         (ilapl(u.right(), G.right())).tr()).tr() * nv(G.left()).norm(),

                    // E11
                          mu * igrad(u.left(), G.left()) * nv(G.left()).normalized() *
                          (igrad(u.left(), G.left()) * nv(G.left()).normalized()).tr() * nv(G.left()).norm(),
                    //-E12
                          -mu * (igrad(u.left(), G.left()) * nv(G.left()).normalized()) *
                          (igrad(u.right(), G.right()) * nv(G.left()).normalized()).tr() * nv(G.left()).norm(),
                    //-E21
                          -mu * (igrad(u.right(), G.right()) * nv(G.left()).normalized()) *
                          (igrad(u.left(), G.left()) * nv(G.left()).normalized()).tr() * nv(G.left()).norm(),
                    // E22
                          mu * igrad(u.right(), G.right()) * nv(G.left()).normalized() *
                          (igrad(u.right(), G.right()) * nv(G.left()).normalized()).tr() * nv(G.left()).norm()
            );
        }
    }


    gsSparseMatrix<> K = A.matrix();

    A.initSystem();
    A.assemble(Density * thickness * u * u.tr() * meas(G));
    // if (method == MethodFlags::NITSCHE)
    // {
    //     // For Nitsche
    //     gsMatrix<real_t> mu_interfaces(mp.nInterfaces(),1);
    //     if (penalty_init == -1.0)
    //         computeStabilityParameter(mp, dbasis, mu_interfaces);

    //     index_t i = 0;
    //     for ( typename gsMultiPatch<real_t>::const_iiterator it = mp.iBegin(); it != mp.iEnd(); ++it, ++i)
    //     {
    //         real_t stab     = 4 * ( dbasis.maxCwiseDegree() + dbasis.dim() ) * ( dbasis.maxCwiseDegree() + 1 );
    //         real_t m_h      = dbasis.basis(0).getMinCellLength(); //*dbasis.basis(0).getMinCellLength();
    //         real_t mu       = 2 * stab / m_h;
    //         real_t alpha = 1;

    //         //mu = penalty_init == -1.0 ? mu : penalty_init / m_h;
    //         if (penalty_init == -1.0)
    //             mu = mu_interfaces(i,0) / m_h;
    //         else
    //             mu = penalty_init / m_h;

    //         std::vector<boundaryInterface> iFace;
    //         iFace.push_back(*it);
    //         A.assembleIfc(iFace,
    //                 //B11
    //                       -alpha * 0.5 * igrad(u.left(), G) * nv(G.left()).normalized() *
    //                       (ilapl(u.left(), G)).tr() * nv(G.left()).norm(),
    //                       -alpha * 0.5 *
    //                       (igrad(u.left(), G) * nv(G.left()).normalized() * (ilapl(u.left(), G)).tr()).tr() *
    //                       nv(G.left()).norm(),
    //                 //B12
    //                       -alpha * 0.5 * igrad(u.left(), G.left()) * nv(G.left()).normalized() *
    //                       (ilapl(u.right(), G.right())).tr() * nv(G.left()).norm(),
    //                       -alpha * 0.5 * (igrad(u.left(), G.left()) * nv(G.left()).normalized() *
    //                                       (ilapl(u.right(), G.right())).tr()).tr() * nv(G.left()).norm(),
    //                 //B21
    //                       alpha * 0.5 * igrad(u.right(), G.right()) * nv(G.left()).normalized() *
    //                       (ilapl(u.left(), G.left())).tr() * nv(G.left()).norm(),
    //                       alpha * 0.5 * (igrad(u.right(), G.right()) * nv(G.left()).normalized() *
    //                                      (ilapl(u.left(), G.left())).tr()).tr() * nv(G.left()).norm(),
    //                 //B22
    //                       alpha * 0.5 * igrad(u.right(), G.right()) * nv(G.left()).normalized() *
    //                       (ilapl(u.right(), G.right())).tr() * nv(G.left()).norm(),
    //                       alpha * 0.5 * (igrad(u.right(), G.right()) * nv(G.left()).normalized() *
    //                                      (ilapl(u.right(), G.right())).tr()).tr() * nv(G.left()).norm(),

    //                 // E11
    //                       mu * igrad(u.left(), G.left()) * nv(G.left()).normalized() *
    //                       (igrad(u.left(), G.left()) * nv(G.left()).normalized()).tr() * nv(G.left()).norm(),
    //                 //-E12
    //                       -mu * (igrad(u.left(), G.left()) * nv(G.left()).normalized()) *
    //                       (igrad(u.right(), G.right()) * nv(G.left()).normalized()).tr() * nv(G.left()).norm(),
    //                 //-E21
    //                       -mu * (igrad(u.right(), G.right()) * nv(G.left()).normalized()) *
    //                       (igrad(u.left(), G.left()) * nv(G.left()).normalized()).tr() * nv(G.left()).norm(),
    //                 // E22
    //                       mu * igrad(u.right(), G.right()) * nv(G.left()).normalized() *
    //                       (igrad(u.right(), G.right()) * nv(G.left()).normalized()).tr() * nv(G.left()).norm()
    //         );
    //     }
    // }

    gsSparseMatrix<> M = A.matrix();

    gsInfo<< "." <<std::flush;// Assemblying done

    gsEigen::GeneralizedSelfAdjointEigenSolver< gsMatrix<real_t>::Base >  eigSolver;
    eigSolver.compute(K,M);
    gsMatrix<> values  = eigSolver.eigenvalues();
    gsMatrix<> vectors = eigSolver.eigenvectors();
    gsInfo<< "." <<std::flush; // Linear solving done
    values = values.cwiseSqrt();

    gsInfo<<"Finished.\n";


    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
        int systemRet = system("mkdir -p ModalResults");
        GISMO_ASSERT(systemRet!=-1,"Something went wrong with calling the system argument");

        gsParaviewCollection collection("ModalResults/modes");

        int N = 1;
        if (!first)
          N = vectors.cols();
        for (index_t m=0; m<N; m++)
        {
            // Compute solution based on eigenmode with number 'mode'
            solVector = vectors.col(m);
            std::string fileName = "ModalResults/modes" + util::to_string(m);
            ev.writeParaview( u_sol, G, fileName);
            for (index_t p = 0; p!=mp.nPatches(); p++)
            {
                fileName = "modes" + util::to_string(m);
                collection.addPart(fileName + ".vts",m,"",p);
            }
        }
        collection.save();
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Analytical part
    ////////////////////////////////////////////////////////////////////////////////

    std::vector<real_t> omegas;
    for (index_t m=1; m!=100; m++)
      for (index_t n=1; n!=100; n++)
        omegas.push_back((math::pow(m/1.0,2)+math::pow(n/1.0,2))*math::pow(3.1415926535,2)*math::sqrt(D / (Density * thickness)));

    std::sort(omegas.begin(),omegas.end());
    GISMO_ENSURE(omegas.size()>=values.rows(),"Too few analytical eigenvalues");
    omegas.resize(values.rows());
    gsAsVector<> analytical(omegas);



    if (write)
    {
        int systemRet = system("mkdir -p ModalResults");
        GISMO_ASSERT(systemRet!=-1,"Something went wrong with calling the system argument");

        gsMatrix<> ratio = values;
        ratio.array() /= analytical.array();

        gsMatrix<> data(values.rows(),4);
        data.col(0) = values;
        data.col(1) = analytical;
        data.col(2) = ratio;
        data.col(3) = gsVector<>::LinSpaced(values.rows(),0,1);
        gsDebugVar(data);


        std::string wnM = "ModalResults/eigenvalues.txt";
        writeToCSVfile(wnM,data);
    }

    return EXIT_SUCCESS;
}// end main
