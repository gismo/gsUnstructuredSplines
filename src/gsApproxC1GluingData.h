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

#include <gsAssembler/gsDofMapperCreator.h>
#include <gsUtils/gsStopwatch.h>

namespace gismo
{

/** @brief Profiling accumulator for gsApproxC1Spline.

    A single process-wide instance (Meyers singleton, see gsApproxC1GetProfile())
    shared by gsApproxC1Spline, gsApproxC1Edge, gsApproxC1Vertex and
    gsApproxC1GluingData. Every field is an accumulated total across all
    edges/vertices/interfaces visited during one gsApproxC1Spline::update()
    call. Printing (gsApproxC1Profile::print()) is only ever invoked when the
    "info" gsOptionList switch is set, so the accumulation itself is gated
    the same way to keep ordinary runs at zero overhead.
*/
struct gsApproxC1Profile
{
    real_t t_init            = 0; // gsApproxC1Spline::init()
    real_t t_gluingData      = 0; // gsApproxC1GluingData::setGlobalGluingData
    real_t t_edgeSetup       = 0; // one-off mass matrix assemble + factorize, per edge side
    real_t t_edgeBfAssemble  = 0; // per-basis-function A.assemble(u*aa) (the RHS quadrature pass)
    real_t t_edgeBfSolve     = 0; // per-basis-function solver.solve(...)
    real_t t_vertex          = 0; // gsApproxC1Vertex construction (per patch-corner, plus the
                                   // once-per-vertex computeKernel() step), excluding nested
                                   // gluing-data time
    real_t t_finalAssembly   = 0; // gsApproxC1Spline::compute() global matrix insertion

    // Sub-buckets of t_vertex (already summed into it, so NOT added again in total()):
    real_t t_vertexSetup      = 0; // one-off mass matrix assemble + factorize, per patch-corner
    real_t t_vertexBfAssemble = 0; // batched A.initVector(6) + getCoeff + assemble(u * aa.tr())
    real_t t_vertexBfSolve    = 0; // the multi-column solver.solve(...)
    real_t t_vertexKernel     = 0; // computeKernel(), once per vertex (0 for internal vertices)

    index_t n_gluingData = 0;
    index_t n_edgeSides  = 0;
    index_t n_edgeBf     = 0;
    index_t n_vertex     = 0;

    void reset() { *this = gsApproxC1Profile(); }

    real_t total() const
    {
        return t_init + t_gluingData + t_edgeSetup + t_edgeBfAssemble
             + t_edgeBfSolve + t_vertex + t_finalAssembly;
    }

    void print(std::ostream & os) const
    {
        const real_t tot = total();
        const real_t pct = (tot > 0 ? (real_t)100.0/tot : (real_t)0.0);
        os << "=== gsApproxC1Spline attributed profile ===\n";
        os << "  init()                    : " << t_init           << " s  (" << t_init*pct           << "%)\n";
        os << "  gluing data (solve)       : " << t_gluingData      << " s  (" << t_gluingData*pct      << "%)  count=" << n_gluingData << "\n";
        os << "  edge mass assemble+factor : " << t_edgeSetup       << " s  (" << t_edgeSetup*pct       << "%)  count=" << n_edgeSides  << " sides\n";
        os << "  edge per-bf RHS assemble  : " << t_edgeBfAssemble  << " s  (" << t_edgeBfAssemble*pct  << "%)  count=" << n_edgeBf     << " bf"
           << (n_edgeBf>0 ? "  (" + util::to_string(t_edgeBfAssemble/n_edgeBf*1000) + " ms/bf)" : "") << "\n";
        os << "  edge per-bf solve         : " << t_edgeBfSolve     << " s  (" << t_edgeBfSolve*pct     << "%)  count=" << n_edgeBf     << " bf\n";
        os << "  vertex construction       : " << t_vertex          << " s  (" << t_vertex*pct          << "%)  count=" << n_vertex     << " corners\n";
        os << "    of which mass setup     : " << t_vertexSetup      << " s\n";
        os << "    of which bf RHS assemble: " << t_vertexBfAssemble << " s\n";
        os << "    of which bf solve       : " << t_vertexBfSolve    << " s\n";
        os << "    of which kernel         : " << t_vertexKernel     << " s  (0 for internal vertices)\n";
        os << "    (breakdown is partial: reparametrisation, mapper/space setup and\n"
              "     parametrizeBasisBack are not separately timed and make up the rest)\n";
        os << "  final matrix assembly     : " << t_finalAssembly   << " s  (" << t_finalAssembly*pct   << "%)\n";
        os << "  TOTAL (accounted)         : " << tot << " s\n";
    }
};

/// Meyers singleton accessor: one gsApproxC1Profile shared by every TU that
/// includes this header, regardless of how many object files instantiate it.
inline gsApproxC1Profile & gsApproxC1GetProfile()
{
    static gsApproxC1Profile profile;
    return profile;
}

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
    const bool profile = m_optionList.getSwitch("info");
    gsStopwatch clock;
    if (profile) clock.restart();

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
    gsDofMapper map = createMapper(BsplineSpace, 1, false);
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

    if (profile)
    {
        gsApproxC1Profile & prof = gsApproxC1GetProfile();
        prof.t_gluingData += clock.stop();
        ++prof.n_gluingData;
    }
} // setGlobalGluingData


} // namespace gismo
