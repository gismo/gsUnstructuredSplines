/** @file gsKLShell/tutorials/linear_shell.cpp

    @brief Tutorial for assembling the Kirchhoff-Love shell

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst
*/

#include <gismo.h>

#include <gsUnstructuredSplines/src/gsDPatch.h>

#ifdef gsKLShell_ENABLED
#include <gsKLShell/src/gsThinShellAssembler.h>
#include <gsKLShell/src/getMaterialMatrix.h>
#include <gsKLShell/src/gsFunctionSum.h>
#endif

using namespace gismo;

// Choose among various shell examples, default = Thin Plate
int main(int argc, char *argv[])
{
#ifdef gsKLShell_ENABLED
    //! [Parse command line]
    bool plot  = false;
    index_t numRefine  = 1;
    index_t numElevate = 1;

    gsCmdLine cmd("Linear shell tutorial.");
    cmd.addInt( "e", "degreeElevation","Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Read geometry]
    // Initialize [ori]ginal and [def]ormed geometry
    gsMultiPatch<> ori, tmp;
    std::string fileName = "surfaces/6p_hyperboloid.xml";
    gsReadFile<>(fileName, ori);
    ori.computeTopology();
    //! [Read geometry]

    //! [Initialize geometry]
    // p-refine
    if (numElevate!=0)
        ori.degreeElevate(numElevate);
    // h-refine
    for (int r =0; r < numRefine; ++r)
        ori.uniformRefine();
    //! [Initialize geometry]

    //! [Construct basis]
    // Construct a multi-basis object
    gsMultiBasis<> bases(ori);
    // Compute the smooth D-Patch construction
    gsDPatch<2,real_t> dpatch(bases);
    dpatch.options().setSwitch("SharpCorners",false);
    dpatch.compute();

    // Store the basis mapping matrix
    gsSparseMatrix<> global2local;
    dpatch.matrix_into(global2local);
    global2local = global2local.transpose();
    // Overwrite the local basis (might change due to D-Patch)
    bases = dpatch.localBasis();

    // Store the modified, smooth geometry
    tmp = ori; // temporary storage
    ori = dpatch.exportToPatches(tmp);

    // Construct a mapped basis, using the local basis and the mapper
    gsMappedBasis<2,real_t> mbasis;
    mbasis.init(bases,global2local);
    //! [Construct basis]

    //! [Set boundary conditions]
    // Define the boundary conditions object
    gsBoundaryConditions<> bc;
    // Set the geometry map for computation of the Dirichlet BCs
    bc.setGeoMap(ori);

    // Set the boundary conditions
    bc.addCondition(0,boundary::west, condition_type::dirichlet, 0, 0 ,false, -1 );
    bc.addCondition(0,boundary::west, condition_type::clamped  , 0, 0 ,false, -1 );
    bc.addCondition(1,boundary::west, condition_type::dirichlet, 0, 0 ,false, -1 );
    bc.addCondition(1,boundary::west, condition_type::clamped  , 0, 0 ,false, -1 );
    //! [Set boundary conditions]

    //! [Set surface force]
    // The surface force is defined in the physical space, i.e. 3D
    // The load is -8000*thickness
    gsFunctionExpr<> force("0","0","-80",3);
    //! [Set surface force]

    //! [Define the material matrix class]
    // Define material parameters
    // The material parameters are defined in the physical domain as well (!)
    gsConstantFunction<> E (2E11,3);
    gsConstantFunction<> nu(0.3 ,3);
    gsConstantFunction<> t (0.01,3);

    // Define a linear material, see \ref gsMaterialMatrixLinear.h
    // The first parameter is the physical domain dimension

    // gsMaterialMatrixLinear<3,real_t> materialMatrix(ori,t,E,nu);

    std::vector<gsFunctionSet<>*> parameters{&E,&nu};
    gsOptionList options;
    options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
    options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
    gsMaterialMatrixBase<real_t> * materialMatrix = getMaterialMatrix<3,real_t>(ori,t,parameters,options);
    //! [Define the material matrix class]

    //! [Define the assembler]
    // Define the assembler.
    // The first parameter is the physical domain dimension
    // The second parameter is the real type
    // The third parameter controls the bending stiffness
    gsThinShellAssembler<3, real_t, true > assembler(ori,bases,bc,force,materialMatrix);
    assembler.setSpaceBasis(mbasis);
    //! [Define the assembler]

    //! [Assemble linear part]
    assembler.assemble();
    gsSparseMatrix<> matrix = assembler.matrix();
    gsVector<> vector = assembler.rhs();
    //! [Assemble linear part]

    //! [Solve linear problem]
    gsInfo<<"Solving system with "<<assembler.numDofs()<<" DoFs\n";
    gsVector<> solVector;
    gsSparseSolver<>::CGDiagonal solver;
    solver.compute( matrix );
    solVector = solver.solve(vector);
    //! [Solve linear problem]

    //! [Construct solution and deformed geometry]
    // We will construct a mapped spline using all coefficients
    gsMatrix<real_t> solFull = assembler.fullSolutionVector(solVector);

    // 2. Reshape all the coefficients to a Nx3 matrix
    GISMO_ASSERT(solFull.rows() % 3==0,"Rows of the solution vector does not match the number of control points");
    solFull.resize(solFull.rows()/3,3);

    // 3. Make the mapped spline
    gsMappedSpline<2,real_t> mspline(mbasis,solFull);

    // Deformation is, on every parametric point, the original geometry + the mapped spline
    gsFunctionSum<real_t> def(&ori,&mspline);
    //! [Construct solution and deformed geometry]


    //! [Evaluate solution]
    // Evaluate the solution on a reference point (parametric coordinates)
    // The reference points are stored column-wise
    gsMatrix<> refPoint(2,1);
    index_t refPatch;
    // In this case, the reference point is on the east boundary of patch 4
    refPoint<<1.0,0.5;
    refPatch = 4;
    // Compute the values
    gsVector<> physpoint = def  .piece(refPatch).eval(refPoint);
    gsVector<> refVals = mspline.piece(refPatch).eval(refPoint);
    // gsInfo << "Displacement at reference point: "<<numVal<<"\n";
    gsInfo  << "Displacement at reference point ("
            <<physpoint.transpose()<<"): "
            <<": "<<refVals.transpose()<<"\n";
    gsInfo  <<"The elastic energy is: "<<0.5*solVector.transpose()*matrix*solVector<<"\n";
    //! [Evaluate solution]

    // ! [Export visualization in ParaView]
    if (plot)
    {
        // Plot the displacements on the deformed geometry
        gsField<> solField(ori, mspline,true);
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview<>( solField, "Deformation", 1000, true);

        // Plot the membrane Von-Mises stresses on the geometry
        gsPiecewiseFunction<> VMm;
        assembler.constructStress(ori,def,VMm,stress_type::von_mises_membrane);
        gsField<> membraneStress(ori,VMm, true);
        gsWriteParaview(membraneStress,"MembraneStress",1000);
    }
    // ! [Export visualization in ParaView]

    return EXIT_SUCCESS;
#else
    GISMO_ERROR("The tutorial needs to be compiled with gsKLShell enabled");
    return EXIT_FAILURE;
#endif
}// end main
