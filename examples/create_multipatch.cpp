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

#include <gsUnstructuredSplines/src/gsMPBESBasis.h>
#include <gsUnstructuredSplines/src/gsMPBESSpline.h>
#include <gsUnstructuredSplines/src/gsApproxC1Spline.h>
#include <gsUnstructuredSplines/src/gsDPatch.h>
#include <gsUnstructuredSplines/src/gsAlmostC1.h>
#include <gsUnstructuredSplines/src/gsC1SurfSpline.h>
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

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot  = false;
    bool mesh  = false;
    index_t method = 0;

    index_t numRefine  = 3;
    index_t degree = 3;
    index_t smoothness = 2;

    std::string fn, basisOutput, geoOutput;

    gsCmdLine cmd("Tutorial on solving a Biharmonic problem with different spaces.");
    // Flags related to the method (default: Approx C1 method)
    cmd.addInt( "m", "method", "The chosen method for the biharmonic problem", method );

    // Flags related to the input/geometry
    cmd.addString( "f", "file", "Input geometry file from path (with .xml)", fn );

    // Flags related to the discrete settings
    cmd.addInt( "p", "degree", "Set the polynomial degree of the basis.", degree );
    cmd.addInt( "s", "smoothness", "Set the smoothness of the basis.",  smoothness );
    cmd.addInt( "r", "numRefine", "Number of refinement-loops.",  numRefine );

    // Flags related to the output
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("mesh", "Plot the mesh", mesh);

    cmd.addString("S", "basisOutput", "Output in xml", basisOutput);
    cmd.addString("G", "geoOutput", "Output in xml", geoOutput);

    // cmd.addString("w", "write", "Write to csv", write);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    if (fn.empty())
        GISMO_ERROR("No file specified");
    else
    {

        GISMO_ENSURE(degree>smoothness,"Degree must be larger than the smoothness!");
        GISMO_ENSURE(smoothness>=0,"Degree must be larger than the smoothness!");
        if (method==3)
            GISMO_ENSURE(smoothness>=1 || smoothness <= degree-2,"Exact C1 method only works for smoothness <= p-2, but smoothness="<<smoothness<<" and p-2="<<degree-2);
        if (method==2 || method==3)
            GISMO_ENSURE(degree > 2,"Degree must be larger than 2 for the approx and exact C1 methods, but it is "<<degree);


        gsMultiPatch<> mp, geom;
        gsMultiBasis<> dbasis;
        gsSparseMatrix<> global2local;
        gsFileData<> fd(fn);
        gsFileData<> out;
        // Take a multipatch, elevate and refine to preferred p and s and construct new matrix and basis
        if (fd.getAnyFirst<gsMultiPatch<>>(mp))
        {
            if (plot) gsWriteParaview(mp,"mp",400,true);

            // Check if any of the patches is a THB patch. If so, we transfer to tensor-bspline on the finest level
            // Also get the max level
            unsigned maxLvl = 0;
            bool THBcheck = false;
            gsHTensorBasis<2,real_t> * hbasis;
            for (index_t p=0; p!=mp.nPatches(); p++)
            {
                if (hbasis = dynamic_cast<gsHTensorBasis<2,real_t> *>(&mp.basis(p)))
                {
                    maxLvl = std::max(maxLvl,hbasis->maxLevel());
                    THBcheck = true;
                }
            }
            gsDebugVar(maxLvl);
            if (THBcheck)
            {
                //This is already one step of the nested refinement; hence
                numRefine--;

                //
                gsTHBSpline<2,real_t> * THBspline;
                gsTensorBSpline<2,real_t> bspline;
                for (index_t p=0; p!=mp.nPatches(); p++)
                {
                    if ((THBspline = dynamic_cast<gsTHBSpline<2,real_t> *>(&mp.patch(p))))
                    {
                        THBspline->convertToBSpline(bspline);
                        if (THBspline->basis().maxLevel()==0)
                            for (unsigned l=0; l!=maxLvl; l++)
                            {
                                bspline.uniformRefine();
                                geom.addPatch(bspline);
                            }
                        else if (THBspline->basis().maxLevel()!=maxLvl)
                            GISMO_ERROR("Something went wrong");
                        else
                        {
                                geom.addPatch(bspline);
                        }
                    }
                }
                geom.computeTopology();
            }
            else
            {
                geom = mp;
                geom.computeTopology();
            }

            gsDebugVar(geom);
            std::vector<std::vector<patchCorner> > corners;
            geom.getEVs(corners);
            for (std::vector<std::vector<patchCorner> >::iterator it = corners.begin(); it!=corners.end(); it++)
                gsDebug<<it[0][0].patch<<","<<it[0][0].corner()<<":\t"<<"valence: "<<it->size()<<"\n";
            geom.getEVs(corners,true);
            for (std::vector<std::vector<patchCorner> >::iterator it = corners.begin(); it!=corners.end(); it++)
                gsDebug<<it[0][0].patch<<","<<it[0][0].corner()<<":\t"<<"valence: "<<it->size()<<"\n";

            gsInfo<<"Refining and elevating geometry..."<<std::flush;
            geom.degreeIncrease(degree-geom.patch(0).degree(0));
            for (int r =0; r < numRefine; ++r)
                geom.uniformRefine(1, degree-smoothness);
            gsInfo<<"Finished.\n";

            dbasis = gsMultiBasis<>(geom);

            gsInfo<<"Constructing spline ..."<<std::flush;
            if (method==-1)
            {
                // identity map
                global2local.resize(dbasis.totalSize(),dbasis.totalSize());
                for (size_t k=0; k!=dbasis.totalSize(); ++k)
                    global2local.coeffRef(k,k) = 1;
                dbasis = gsMultiBasis<>(geom);
            }
            else if (method==0)
            {
                gsMPBESSpline<2,real_t> cgeom(geom,3);
                gsMappedBasis<2,real_t> basis = cgeom.getMappedBasis();
                auto container = basis.getBasesCopy();

                global2local = basis.getMapper().asMatrix();
                geom = cgeom.exportToPatches();
                dbasis = gsMultiBasis<>(container,geom.topology());
            }
            else if (method==1)
            {
                gsDPatch<2,real_t> dpatch(dbasis);
                dpatch.compute();
                dpatch.matrix_into(global2local);

                global2local = global2local.transpose();
                geom = dpatch.exportToPatches(geom);
                dbasis = dpatch.localBasis();
            }
            else if (method==2) // Pascal
            {
                // The approx. C1 space
                gsApproxC1Spline<2,real_t> approxC1(geom,dbasis);
                approxC1.options().setSwitch("interpolation",true);
                approxC1.options().setInt("gluingDataDegree",-1);
                approxC1.options().setInt("gluingDataSmoothness",-1);

                global2local = approxC1.getSystem();
                global2local = global2local.transpose();
                approxC1.getMultiBasis(dbasis);
            }
            else if (method==3) // Andrea
            {
                gsC1SurfSpline<2,real_t> smoothC1(geom,dbasis);
                smoothC1.init();
                smoothC1.compute();

                global2local = smoothC1.getSystem();
                global2local = global2local.transpose();
                smoothC1.getMultiBasis(dbasis);
            }
            else if (method==4)
            {
                gsAlmostC1<2,real_t> almostC1(geom);
                almostC1.options().setSwitch("SharpCorners",false);
                almostC1.compute();
                almostC1.matrix_into(global2local);

                global2local = global2local.transpose();
                geom = almostC1.exportToPatches();
                dbasis = almostC1.localBasis();
            }

            if (plot) gsWriteParaview(geom,"geom",400,true);
        }
        else if (fd.has<gsMultiBasis<>>() && fd.has<gsSparseMatrix<>>())
        {

        }
        gsInfo<<"Finished"<<std::flush;


        gsMappedBasis<2,real_t> mbasis(dbasis,global2local);
        gsDebugVar(mbasis.globalSize());

        if (!basisOutput.empty())
        {
            gsInfo<<"Writing data to "<<basisOutput<<std::flush;
            out.add(global2local);
            out.add(dbasis);
            out.save(basisOutput);
            gsInfo<<"Finished";
        }
        if (!geoOutput.empty())
        {
            gsInfo<<"Writing geometry to "<<geoOutput<<std::flush;
            gsWrite(geom,geoOutput);
            gsInfo<<"Finished";
        }


    }

    return EXIT_SUCCESS;
}// end main
