/** @file biharmonic_approxC1_example.cpp

    @brief Solving the biharmonic equation on a multi-patch using Approx C1 Basis

    Author(s): P. Weinmueller
 **/

//! [Include namespace]
#include <gismo.h>

#include <gsUnstructuredSplines/src/gsApproxC1Spline.h>
#include <gsUnstructuredSplines/src/gsBiharmonicAssembler2.h>
using namespace gismo;
//! [Include namespace]


int main()
{
  //! [Problem input]
  bool second = false; // Computing the first biharmonic equation (Dirichlet+Neumann BCs)
  index_t numRefine = 3, degree = 3, smoothness = 2;
  std::string string_geo = "planar/square_6patches_bicubic.xml";
  //! [Problem input]
  
  gsMultiPatch<real_t> mp;
  gsBoundaryConditions<> bc;
  gsFunctionExpr<real_t> f, ms;
  // Set the discretization space
  gsMappedBasis<2,real_t> bb2;    
  
  gsInfo << "Filedata: " << string_geo << "\n";
  gsReadFile<>(string_geo, mp);
  mp.clearTopology();
  mp.computeTopology();

  gsMultiBasis<real_t> dbasis(mp, false); //true: poly-splines (not NURBS)
  dbasis.setDegree( degree); // preserve smoothness
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
  for (int r=0; r<=numRefine; ++r) {
    dbasis.uniformRefine(1, degree-smoothness);
    meshsize[r] = dbasis.basis(0).getMinCellLength();

    biharmonicAssembler.initialize();
    
    // The approx. C1 space
    timer.restart();    
    gsApproxC1Spline<2,real_t> approxC1(mp,dbasis);
    approxC1.update(bb2);
    setup_time += timer.stop();
    gsInfo<< "." <<std::flush; // Approx C1 construction done

    // Setup the system
    biharmonicAssembler.assemble(bb2);
    
    ma_time += timer.stop();
    gsInfo<< "." <<std::flush; // Assemblying done

    dofs[r] = biharmonicAssembler.numDofs();
    gsInfo << dofs[r] <<std::flush;

    // Solving the system
    timer.restart();    
    biharmonicAssembler.solve();
    slv_time += timer.stop();
    gsInfo<< "." <<std::flush; // Linear solving done
    
    timer.restart();
    biharmonicAssembler.computeError(bb2, l2err[r], h1err[r], h2err[r], IFaceErr[r]);
    err_time += timer.stop();
    gsInfo<< ". " <<std::flush; // Error computations done
  } //for loop
  //! [Solver loop]

  timer.stop();
  gsInfo << "\n\nTotal time: " << setup_time+ma_time+slv_time+err_time << "\n";
  gsInfo << "     Setup: "     << setup_time << "\n";
  gsInfo << "  Assembly: "     << ma_time    << "\n";
  gsInfo << "   Solving: "     << slv_time   << "\n";
  gsInfo << "     Norms: "     << err_time   << "\n";
  
  //! [Error and convergence rates]  
  gsInfo<< "\nL2 error: "            <<std::scientific<<std::setprecision(3)<<l2err.transpose()<< "\n";
  gsInfo<< "H1 error: "              <<std::scientific<<h1err.transpose()<< "\n";
  gsInfo<< "H2 error: "              <<std::scientific<<h2err.transpose()<< "\n";
  gsInfo<< "Deriv Interface error: " <<std::scientific<<IFaceErr.transpose()<< "\n";
  
  if (numRefine>0) {
    gsInfo<< "EoC (L2): " << std::fixed<<std::setprecision(2)
	  <<  ( l2err.head(numRefine).array() /
		  l2err.tail(numRefine).array() ).log().transpose() / std::log(2.0)
	  << "\n";
    gsInfo<< "EoC (H1): "<< std::fixed<<std::setprecision(2)
	  <<( h1err.head(numRefine).array() /
	      h1err.tail(numRefine).array() ).log().transpose() / std::log(2.0) << "\n";
    gsInfo<< "EoC (H2): "<< std::fixed<<std::setprecision(2)
	    <<( h2err.head(numRefine).array() /
		h2err.tail(numRefine).array() ).log().transpose() / std::log(2.0) << "\n";
    gsInfo<< "EoC (Iface): "<< std::fixed<<std::setprecision(2)
	  <<( IFaceErr.head(numRefine).array() /
	      IFaceErr.tail(numRefine).array() ).log().transpose() / std::log(2.0) << "\n";
  }
  //! [Error and convergence rates]
}
