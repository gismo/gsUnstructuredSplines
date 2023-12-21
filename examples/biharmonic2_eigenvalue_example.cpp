/** @file biharmonic_example.cpp

    @brief A Biharmonic example for a single patch.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/

# include <gismo.h>
#include <gsSpectra/gsSpectra.h>

using namespace gismo;

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


void setMapperForBiharmonic(gsBoundaryConditions<> & bc, gsMultiBasis<> & basis, gsDofMapper & mapper)
{
    mapper.init(basis);

    for (gsBoxTopology::const_iiterator it = basis.topology().iBegin();
         it != basis.topology().iEnd(); ++it) // C^0 at the interface
    {
        basis.matchInterface(*it, mapper);
    }

    gsMatrix<index_t> bnd;
    for (typename gsBoundaryConditions<real_t>::const_iterator
                 it = bc.begin("Dirichlet"); it != bc.end("Dirichlet"); ++it)
    {
        bnd = basis.basis(it->ps.patch).boundary(it->ps.side());
        mapper.markBoundary(it->ps.patch, bnd, 0);
    }

    for (typename gsBoundaryConditions<real_t>::const_iterator
                 it = bc.begin("Neumann"); it != bc.end("Neumann"); ++it)
    {
        bnd = basis.basis(it->ps.patch).boundaryOffset(it->ps.side(),1);
        mapper.markBoundary(it->ps.patch, bnd, 0);
    }
    mapper.finalize();
}

void setDirichletNeumannValuesL2Projection(gsMultiPatch<> & mp, gsMultiBasis<> & basis, gsBoundaryConditions<> & bc, const expr::gsFeSpace<real_t> & u)
{
    const gsDofMapper & mapper = u.mapper();
    gsDofMapper mapperBdy(basis, u.dim());
    for (gsBoxTopology::const_iiterator it = basis.topology().iBegin();
         it != basis.topology().iEnd(); ++it) // C^0 at the interface
    {
        basis.matchInterface(*it, mapperBdy);
    }
    for (size_t np = 0; np < mp.nPatches(); np++)
    {
        gsMatrix<index_t> bnd = mapper.findFree(np);
        mapperBdy.markBoundary(np, bnd, 0);
    }
    mapperBdy.finalize();

    gsExprAssembler<real_t> A(1,1);
    A.setIntegrationElements(basis);

    auto G = A.getMap(mp);
    auto uu = A.getSpace(basis);
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
    gsMatrix<real_t> fixedDofs_temp = solver.solve(A.rhs());

    // Reordering the dofs of the boundary
    fixedDofs.setZero(mapper.boundarySize(),1);
    index_t sz = 0;
    for (size_t np = 0; np < mp.nPatches(); np++)
    {
        gsMatrix<index_t> bnd = mapperBdy.findFree(np);
        bnd.array() += sz;
        for (index_t i = 0; i < bnd.rows(); i++)
        {
            index_t ii = mapperBdy.asVector()(bnd(i,0));
            fixedDofs(mapper.global_to_bindex(mapper.asVector()(bnd(i,0))),0) = fixedDofs_temp(ii,0);
        }
        sz += mapperBdy.patchSize(np,0);
    }
}


int main(int argc, char *argv[])
{
    // Input options
    index_t degree      = 2;
    index_t smoothness  = 1;
    index_t numHref     = 4;
    bool plot       = false;
    bool write      = false;
    bool first      = false;

    gsCmdLine cmd("Modal analysis for thin shells.");
    cmd.addInt("r","hRefine",
               "Number of dyadic h-refinement (bisection) steps to perform before solving",
               numHref);
    cmd.addInt("d","degree","Set the degree",degree);
    cmd.addInt("s","smoothness","Set the smoothness",smoothness);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addSwitch("write", "Write convergence data to file", write);
    cmd.addSwitch("first", "Plot only first", first);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //! [Parse command line]

    //! [Read geometry]
    gsMultiPatch<> mp;

    real_t thickness = 0.01;
    real_t E_modulus = 1e5;
    real_t Density = 1e5;
    real_t PoissonRatio = 0.3;
    real_t D = E_modulus*math::pow(thickness,3)/(12*(1-math::pow(PoissonRatio,2)));
    // real_t D = 1;
    // real_t Density = 1;
    // real_t thickness = 1;


    mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // degree
    mp.addAutoBoundaries();

    if (mp.nPatches() != 1)
    {
        gsInfo << "The geometry has more than one patch. Run the code with a single patch!\n";
        return EXIT_FAILURE;
    }
    //! [Read geometry]

    //! [Refinement]
    gsMultiBasis<> basis(mp, false);//true: poly-splines (not NURBS)

    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    basis.setDegree(degree); // preserve smoothness

    for (index_t r =0; r < numHref; ++r)
        basis.uniformRefine(1, degree-smoothness);

    gsInfo<<"Basis has "<<basis.totalElements()<<" elements\n";
    gsInfo<<"Basis is: \n"<<basis.basis(0)<<"\n";
    //! [Boundary condition]
    gsBoundaryConditions<> bc;
    for (gsMultiPatch<>::const_biterator bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
        bc.addCondition(*bit, condition_type::dirichlet, 0);

    bc.setGeoMap(mp);
    gsInfo << "Boundary conditions:\n" << bc << "\n";
    //! [Boundary condition]

    //! [Problem setup]
    gsExprAssembler<real_t> A(1,1);
    //gsInfo<<"Active options:\n"<< A.options() <<"\n";

    // Elements used for numerical integration
    A.setIntegrationElements(basis);
    gsExprEvaluator<real_t> ev(A);

    // Set the geometry map
    auto G = A.getMap(mp);

    // Set the discretization space
    auto u = A.getSpace(basis);

    // Solution vector and solution variable
    gsMatrix<real_t> solVector;
    auto u_sol = A.getSolution(u, solVector);

    // Setup the system
    u.setup(bc,dirichlet::l2Projection,0);

    // Initialize the system
    A.initSystem();

    gsInfo<<"Solving system with "<< A.numDofs() <<" degrees of freedom"<<std::flush;

    // Compute the system matrix and right-hand side
    A.assemble(D * ilapl(u, G) * ilapl(u, G).tr() * meas(G));
    gsSparseMatrix<> K = A.matrix();

    A.initSystem();
    A.assemble(Density * thickness * u * u.tr() * meas(G));
    gsSparseMatrix<> M = A.matrix();


    // gsSpectraGenSymShiftSolver<gsSparseMatrix<>,Spectra::GEigsMode::ShiftInvert> solver(K,M,A.numDofs()-1,2*(A.numDofs()-1),0.01);
    // solver.init();
    // solver.compute(Spectra::SortRule::SmallestAlge,1000,1e-6,Spectra::SortRule::SmallestAlge);
    // gsMatrix<> values  = solver.eigenvalues();
    // gsMatrix<> vectors = solver.eigenvectors();

    gsEigen::GeneralizedSelfAdjointEigenSolver< typename gsMatrix<>::Base >  eigSolver;
    eigSolver.compute(K,M);
    gsMatrix<> values  = eigSolver.eigenvalues();
    gsMatrix<> vectors = eigSolver.eigenvectors();

    gsInfo<<"Finished.\n";

    values = values.cwiseSqrt();

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
            fileName = "modes" + util::to_string(m) + "0";
            collection.addPart(fileName + ".vts",m);
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
    GISMO_ENSURE((omegas.size()>=values.rows()),"Too few analytical eigenvalues");
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
    return  EXIT_SUCCESS;
}
