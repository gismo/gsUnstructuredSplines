/** @file gsBiharmonicAssembler-test.cpp

    @brief Tutorial on how to use expression assembler and the (approx.) C1 basis function
                to solve the Biharmonic equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/

#include <gismo.h>
#include <gsUnstructuredSplines/src/gsApproxC1Spline.h>
#include <gsUnstructuredSplines/src/gsBiharmonicAssembler2.h>

using namespace gismo;

int main()
{
  gsInfo << "Test loaded successful\n";
  
  bool second = false;

  index_t numRefine  = 3;
  index_t degree = 3;
  index_t smoothness = 2;
  
  gsMultiPatch<real_t> mp;
  gsBoundaryConditions<> bc;
  gsFunctionExpr<real_t> f, ms;
  // Set the discretization space
  gsMappedBasis<2,real_t> bb2;    
  
  std::string string_geo = "planar/geometries/g1021.xml";
  gsInfo << "Filedata: " << string_geo << "\n";
  gsReadFile<>(string_geo, mp);
  mp.clearTopology();
  mp.computeTopology();

  gsMultiBasis<real_t> dbasis(mp, false);//true: poly-splines (not NURBS)
  dbasis.setDegree(degree); // preserve smoothness
  if (dbasis.basis(0).numElements() < 4) 
    dbasis.uniformRefine(1, degree-smoothness);

  //! [Problem setup]
  gsBiharmonicAssembler2<real_t> biharmonicAssembler(mp, dbasis, second);
  //! [Problem setup]

  //! [Solver loop]
  gsVector<real_t> l2err(numRefine+1), h1err(numRefine+1), h2err(numRefine+1),
    IFaceErr(numRefine+1), meshsize(numRefine+1), dofs(numRefine+1);
 
  gsInfo<< "(dot1=approxC1construction, dot2=assembled, dot3=solved, dot4=got_error)\n"
    "\nDoFs: ";
  double setup_time(0), ma_time(0), slv_time(0), err_time(0);
  gsStopwatch timer;
  for (int r=0; r<=numRefine; ++r) 
  {
    dbasis.uniformRefine(1, degree-smoothness);

    meshsize[r] = dbasis.basis(0).getMinCellLength();

    biharmonicAssembler.initialize();
    
    // // The approx. C1 space
    // gsApproxC1Spline<2,real_t> approxC1(mp,dbasis);
    // approxC1.update(bb2);
    // gsInfo<< "." <<std::flush; // Approx C1 construction done
    
    // DPatch
    gsDPatch<2,real_t> dpatch(dbasis);
    dpatch.options().setInt("RefLevel",r);
    dpatch.options().setInt("Pi",0);
    dpatch.options().setSwitch("SharpCorners",false);
    dpatch.compute();
    dpatch.matrix_into(global2local);
    gsSparseMatrix<real_t> global2local = global2local.transpose();
    geom = dpatch.exportToPatches();
    gsMultiBasis<> localbasis = dpatch.localBasis();
    bb2.init(localbasis,global2local);

    // Setup the system
    biharmonicAssembler.assemble(bb2);
    
    ma_time += timer.stop();
    gsInfo<< "." <<std::flush;// Assemblying done

    dofs[r] = biharmonicAssembler.numDofs();
    gsInfo<<dofs[r];
    timer.restart();    
    biharmonicAssembler.solve();
    slv_time += timer.stop();
    gsInfo<< "." <<std::flush; // Linear solving done
    
    timer.restart();
    biharmonicAssembler.computeError(bb2, l2err[r], h1err[r], h2err[r]);    
    gsInfo<< ". " <<std::flush; // Error computations done
  } //for loop
    //! [Solver loop]

  timer.stop();
  gsInfo<<"\n\nTotal time: "<< setup_time+ma_time+slv_time+err_time <<"\n";
  gsInfo<<"     Setup: "<< setup_time <<"\n";
  gsInfo<<"  Assembly: "<< ma_time    <<"\n";
  gsInfo<<"   Solving: "<< slv_time   <<"\n";
  gsInfo<<"     Norms: "<< err_time   <<"\n";

  gsInfo<< "\nMeshsize: "<<meshsize.transpose()<<"\n";
    gsInfo<< "Dofs: "<<dofs.transpose()<<"\n";
  gsInfo<< "L2 error: "<<std::scientific<<std::setprecision(3)<<l2err.transpose()<<"\n";
  gsInfo<< "H1 error: "<<std::scientific<<h1err.transpose()<<"\n";
  gsInfo<< "H2 error: "<<std::scientific<<h2err.transpose()<<"\n";
  
  if (numRefine>0) {
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
   }  
}

