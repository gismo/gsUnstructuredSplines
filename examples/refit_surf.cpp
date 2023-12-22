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
    bool noTHB  = false;

    index_t method = 0;

    index_t numRefine  = 3;
    index_t degree = 3;
    index_t smoothness = 2;

    std::string fn, GeomOut, SysOut;

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

    cmd.addString("G", "GeomOut", "Geometry output in xml", GeomOut);
    cmd.addString("S", "SysOut", "System (basis and matrix) output in xml", SysOut);

    cmd.addSwitch("noTHB", "No THB in export", noTHB);


    // cmd.addString("w", "write", "Write to csv", write);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Initialize data]
    gsMultiPatch<real_t> mp;

    gsReadFile<>(fn,mp);

    gsInfo<<"Refining and elevating geometry..."<<std::flush;
    mp.degreeIncrease(degree-mp.patch(0).degree(0));
    for (int r =0; r < numRefine; ++r)
        mp.uniformRefine(1, degree-smoothness);
    gsWriteParaview(mp,"mp",400,true);
    gsInfo<<"Finished.\n";

    for (index_t p=0; p!=mp.nPatches(); p++)
        gsDebug<<"Basis "<<p<<": "<<mp.basis(p)<<"\n";

    gsInfo<<"Computing topology..."<<std::flush;
    mp.computeTopology();
    mp.fixOrientation();
    gsInfo<<"Finished.\n";

    // STEP 1: Get curve network with merged linear interfaces
    gsInfo<<"Exporting curve network..."<<std::flush;
    mp.constructInterfaceRep();
    mp.constructBoundaryRep();
    auto & irep = mp.interfaceRep();
    auto & brep = mp.boundaryRep();
    gsDebug <<" irep "<< irep.size() <<" \n" ;
    gsDebug <<" brep "<< brep.size() <<" \n" ;


    std::vector<std::vector<patchCorner> > cornerLists;
    mp.getEVs(cornerLists,true);
    for (index_t k=0; k!=cornerLists.size(); k++)
    {
        gsMatrix<> points(3,cornerLists[k].size());
        for (index_t l = 0; l!=cornerLists[k].size(); l++)
        {
            patchCorner corner = cornerLists[k].at(l);
            gsVector<bool> pars;
            corner.corner().parameters_into(mp.parDim(),pars); // get the parametric coordinates of the corner
            gsMatrix<> supp = mp.basis(corner.patch).support();
            gsVector<> vec(supp.rows());
            for (index_t d = 0; d!=supp.rows(); d++)
                vec(d) = supp(d,pars(d));

            gsMatrix<> tmp;
            mp.patch(corner.patch).eval_into(vec,tmp);
            points.col(l) = tmp;
        }
        gsWriteParaviewPoints(points,"points_" + std::to_string(k));
    }

    // outputing...
    gsMultiPatch<> crv_net, iface_net, bnd_net;
    for (auto it = irep.begin(); it!=irep.end(); ++it)
    {
        iface_net.addPatch((*it->second));
        crv_net.addPatch((*it->second));
    }
    for (auto it = brep.begin(); it!=brep.end(); ++it)
    {
        bnd_net.addPatch((*it->second));
        crv_net.addPatch((*it->second));
    }

    if (plot) gsWriteParaview(iface_net,"iface_net",100);
    if (plot) gsWriteParaview(bnd_net,"bnd_net",100);
    if (plot) gsWriteParaview(crv_net,"crv_net",100);
    gsInfo<<"Finished\n";

    gsSparseMatrix<> global2local;
    gsMultiPatch<> geom;
    gsMultiBasis<> dbasis;
    gsMappedBasis<2,real_t> bb2;

    if (method==0)
    {
        gsInfo<<"Computing almost C1..."<<std::flush;
        gsDPatch<2,real_t> dpatch(mp);
        dpatch.options().setSwitch("SharpCorners",false);
        dpatch.compute();
        gsInfo<<"Finished.\n";
        dpatch.matrix_into(global2local);
        global2local = global2local.transpose();
        geom = dpatch.exportToPatches();
        dbasis = dpatch.localBasis();
        bb2.init(dbasis,global2local);
    }
    else if (method==1)
    {
        gsInfo<<"Computing almost C1..."<<std::flush;
        gsAlmostC1<2,real_t> almostC1(mp);
        almostC1.options().setSwitch("SharpCorners",false);
        almostC1.compute();
        gsInfo<<"Finished.\n";
        almostC1.matrix_into(global2local);
        global2local = global2local.transpose();
        geom = almostC1.exportToPatches();
        dbasis = almostC1.localBasis();
        bb2.init(dbasis,global2local);
    }
    else
        GISMO_ERROR("Method " + std::to_string(method) + " unknown");


    // project geometry on bb2
    gsMatrix<> coefs;
    gsL2Projection<real_t>::projectGeometry(dbasis,bb2,mp,coefs);
    gsMatrix<> allCoefs = global2local*coefs.reshape(global2local.cols(),mp.geoDim());

    // Substitute all coefficients of geom with the newly obtained ones.
    gsMultiBasis<> geombasis(geom);
    gsDofMapper mapper(geombasis);
    mapper.finalize();
    for (index_t p = 0; p != geom.nPatches(); p++)
    {
        for (index_t k=0; k!=mapper.patchSize(p); k++)
        {
            geom.patch(p).coefs().row(k) = allCoefs.row(mapper.index(k,p));
        }
    }


    if (noTHB)
    {
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
    }
    if (plot) gsWriteParaview(geom,"geom",1000,mesh);

    if (!GeomOut.empty())
        gsWrite<>(geom,GeomOut);
    if (!SysOut.empty())
    {
        gsFileData<> fd;
        fd.add(global2local);
        fd.add(dbasis);
        fd.save(SysOut);
    }


    return EXIT_SUCCESS;
}// end main
