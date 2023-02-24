/** @file gsGluingData.h

    @brief Compute the gluing data for one interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/

#pragma once

#include <gsUnstructuredSplines/src/gsApproxC1Utils.h>

namespace gismo
{

template<short_t d, class T>
class gsApproxC1GluingData
{
private:
    typedef typename std::vector<gsPatchReparameterized<d,T>> C1AuxPatchContainer;

    /// Shared pointer for gsApproxGluingData
    typedef memory::shared_ptr<gsApproxC1GluingData> Ptr;

    /// Unique pointer for gsApproxGluingData
    typedef memory::unique_ptr<gsApproxC1GluingData> uPtr;

public:
    gsApproxC1GluingData()
    { }


    gsApproxC1GluingData(C1AuxPatchContainer const & auxPatchContainer,
                       gsOptionList const & optionList,
                       std::vector<patchSide> sidesContainer,
                       std::vector<bool> isInterface = std::vector<bool>{},
                       gsTensorBSplineBasis<d, T> basis = gsTensorBSplineBasis<d, T>())
        : m_auxPatches(auxPatchContainer), m_optionList(optionList)
    {
        alphaSContainer.resize(2);
        betaSContainer.resize(2);
        if (m_auxPatches.size() == 2) // Interface
        {
            setGlobalGluingData(1,0); // u
            setGlobalGluingData(0,1); // v
        }
        else if (m_auxPatches.size() == 1 && sidesContainer.size() == 2) // Vertex
        {
            for (size_t dir = 0; dir < sidesContainer.size(); dir++)
            {
                // Map global side to local side
                index_t localSide = auxPatchContainer[0].getMapIndex(sidesContainer[dir].index());
                //gsInfo << "Global: " << sidesContainer[dir] << " : " << localSide << "\n";
                index_t localDir = localSide < 3 ? 1 : 0;

                createGluingDataSpace(m_auxPatches[0].getPatchRotated(), basis,
                                      localDir, bsp_gD, m_optionList.getInt("gluingDataDegree"), m_optionList.getInt("gluingDataSmoothness"));

                if(isInterface[dir]) // West
                {
                    setGlobalGluingData(0, localDir);
                }

                else
                {
                    // empty
                }
            }
        }
        //else
        //    gsInfo << "I am here \n";

    }

    // Computed the gluing data globally
    void setGlobalGluingData(index_t patchID = 0,  index_t dir = 1);

    gsBSpline<T> & alphaS(index_t patchID) { return alphaSContainer[patchID]; }
    gsBSpline<T> & betaS(index_t patchID) { return betaSContainer[patchID]; }

protected:

    // Spline space for the gluing data + multiPatch
    C1AuxPatchContainer m_auxPatches;

    gsBSplineBasis<T> bsp_gD;

    const gsOptionList m_optionList;

    // Result
    std::vector<gsBSpline<T>> alphaSContainer, betaSContainer;

}; // class gsGluingData


template<short_t d, class T>
void gsApproxC1GluingData<d, T>::setGlobalGluingData(index_t patchID, index_t dir)
{
    // Interpolate boundary yes or no //
    bool interpolate_boundary = false;
    // Interpolate boundary yes or no //

    gsTensorBSplineBasis<d, T> basis = dynamic_cast<const gsTensorBSplineBasis<d, T> &>(m_auxPatches[patchID].getBasisRotated().piece(0));;
    if (m_auxPatches.size() == 2) // Interface
    {
        gsTensorBSplineBasis<d, T> basis2 = dynamic_cast<const gsTensorBSplineBasis<d, T> &>(m_auxPatches[1-patchID].getBasisRotated().piece(0));
        if (basis.component(dir).numElements() > basis2.component(1-dir).numElements())
            basis.component(dir) = basis2.component(1-dir);

        // ======== Space for gluing data : S^(p_tilde, r_tilde) _k ========
        createGluingDataSpace(m_auxPatches[patchID].getPatchRotated(), basis,
                              dir, bsp_gD, m_optionList.getInt("gluingDataDegree"), m_optionList.getInt("gluingDataSmoothness"));
    }



    //! [Problem setup]
    typename gsSparseSolver<T>::SimplicialLDLT solver;
    gsExprAssembler<T> A(1,1);

    // Elements used for numerical integration
    gsMultiBasis<T> BsplineSpace(bsp_gD);
    A.setIntegrationElements(BsplineSpace);
    gsExprEvaluator<T> ev(A);

    gsAlpha<T> alpha(m_auxPatches[patchID].getPatchRotated(), dir);
    auto aa = A.getCoeff(alpha);

    // Set the discretization space
    auto u = A.getSpace(BsplineSpace);

    // Create Mapper
    gsDofMapper map(BsplineSpace);
    gsMatrix<index_t> act(2,1);
    act(0,0) = 0;
    act(1,0) = BsplineSpace[0].size()-1; // First and last
    if (interpolate_boundary)
        map.markBoundary(0, act); // Patch 0
    map.finalize();

    u.setupMapper(map);

    gsMatrix<T> & fixedDofs = const_cast<expr::gsFeSpace<T>&>(u).fixedPart();
    fixedDofs.setZero( u.mapper().boundarySize(), 1 );

    // For the boundary
    gsMatrix<T> points_bdy(1,2);
    points_bdy << 0.0, 1.0;
    if (interpolate_boundary)
        fixedDofs = alpha.eval(points_bdy).transpose();

    A.initSystem();
    A.assemble(u * u.tr(), u * aa);

    solver.compute( A.matrix() );
    gsMatrix<T> solVector = solver.solve(A.rhs());

    auto u_sol = A.getSolution(u, solVector);
    gsMatrix<T> sol;
    u_sol.extract(sol);

    typename gsGeometry<T>::uPtr tilde_temp;
    tilde_temp = bsp_gD.makeGeometry(sol);
    alphaSContainer[dir] = dynamic_cast<gsBSpline<T> &> (*tilde_temp);

    gsBeta<T> beta(m_auxPatches[patchID].getPatchRotated(), dir);
    auto bb = A.getCoeff(beta);

    // For the boundary
    if (interpolate_boundary)
        fixedDofs = beta.eval(points_bdy).transpose();

    A.initSystem();
    A.assemble(u * u.tr(), u * bb);

    solver.compute( A.matrix() );
    solVector = solver.solve(A.rhs());

    auto u_sol2 = A.getSolution(u, solVector);
    u_sol2.extract(sol);

    tilde_temp = bsp_gD.makeGeometry(sol);
    betaSContainer[dir] = dynamic_cast<gsBSpline<T> &> (*tilde_temp);

} // setGlobalGluingData


} // namespace gismo
