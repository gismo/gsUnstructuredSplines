/** @file gsGluingData.h

    @brief Compute the gluing data for one interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Farahat
*/

#pragma once

#include <gsUnstructuredSplines/src/gsC1SurfGD.h>
#include <gsUnstructuredSplines/src/gsC1SurfGluingDataAssembler.h>
#include <gsUnstructuredSplines/src/gsC1SurfGluingDataVisitor.h>
#include <gsUnstructuredSplines/src/gsG1AuxiliaryPatch.h>

#include <gsUtils/gsStopwatch.h>
#include <cstdlib>

namespace gismo
{

/** @brief Profiling accumulator for gsC1SurfSpline / AS-G1.

    Mirrors gsApproxC1Profile (see gsApproxC1GluingData.h). Unlike the
    Approx-C1 classes, none of gsC1SurfGluingData/gsC1SurfGluingDataAssembler/
    gsC1SurfBasisEdge/gsC1SurfBasisVertex take a gsOptionList, so there is no
    existing "info" switch to gate on here. Printing (and therefore the
    accumulation, kept zero-cost otherwise) is gated on the environment
    variable GISMO_C1SURF_PROFILE instead -- set it to any non-empty value to
    turn profiling on. This is a deliberate deviation from the "info" switch
    used for Approx-C1, needed because AS-G1
    has no equivalent option infrastructure to hook into without changing
    out-of-scope files (gsC1SurfSpline.h/.hpp).

    gsC1SurfSpline::compute() is out of scope for this task, so there is no
    single point after the last edge/vertex where a one-shot final report can
    be printed. Instead gsC1SurfBasisEdge::setG1BasisEdge() (once per edge/
    boundary) and gsC1SurfBasisVertex::setG1BasisVertex() (once per vertex)
    each reprint the *cumulative* report; the last one printed during a run
    is the complete total.
*/
struct gsC1SurfProfile
{
    real_t t_gluingData   = 0; // gsC1SurfGluingData construction (7x7 + 4x4 solves, once per interface/boundary)
    real_t t_edgeApply    = 0; // gsC1SurfBasisEdge::apply(): the FULL 2-D beginAll()/endAll() quadrature pass, once per patch side
    real_t t_edgeSolve    = 0; // gsC1SurfBasisEdge::setG1BasisEdge(): one factorisation + multi-column solve, once per patch side
    real_t t_vertexApply  = 0; // gsC1SurfBasisVertex::apply(): ONE full 2-D pass builds the single six-column system
    real_t t_vertexSolve  = 0; // gsC1SurfBasisVertex::solve(): one factorisation + six-column solve per vertex

    index_t n_gluingData      = 0;
    index_t n_edgeApplyCalls  = 0; // == number of edge (patch-side) apply() calls
    index_t n_edgeElemVisits  = 0; // sum, over all edgeApply calls, of elements visited per call
    index_t n_vertexApplyCalls = 0;
    index_t n_vertexElemVisits = 0;

    void reset() { *this = gsC1SurfProfile(); }

    real_t total() const { return t_gluingData + t_edgeApply + t_edgeSolve + t_vertexApply + t_vertexSolve; }

    void print(std::ostream & os) const
    {
        const real_t tot = total();
        const real_t pct = (tot > 0 ? (real_t)100.0/tot : (real_t)0.0);
        os << "=== gsC1SurfSpline (AS-G1) attributed profile (cumulative) ===\n";
        os << "  gluing data solves        : " << t_gluingData  << " s  (" << t_gluingData*pct << "%)  count=" << n_gluingData << "\n";
        os << "  edge apply() [full 2D pass, per patch side] : " << t_edgeApply << " s  (" << t_edgeApply*pct << "%)  calls="
           << n_edgeApplyCalls << "  total-elem-visits=" << n_edgeElemVisits
           << (n_edgeApplyCalls>0 ? "  (" + util::to_string((double)n_edgeElemVisits/(double)n_edgeApplyCalls) + " elem/call)" : "") << "\n";
        os << "  edge solve (compute+solve): " << t_edgeSolve << " s  (" << t_edgeSolve*pct << "%)  calls=" << n_edgeApplyCalls << "\n";
        os << "  vertex apply() [1 full 2D pass, all 6 bf together] : " << t_vertexApply << " s  (" << t_vertexApply*pct << "%)  calls="
           << n_vertexApplyCalls << "  total-elem-visits=" << n_vertexElemVisits << "\n";
        os << "  vertex solve (compute+six-column solve) : " << t_vertexSolve << " s  (" << t_vertexSolve*pct << "%)  calls=" << n_vertexApplyCalls << "\n";
        os << "  TOTAL (accounted)         : " << tot << " s\n";
    }
};

inline gsC1SurfProfile & gsC1SurfGetProfile()
{
    static gsC1SurfProfile profile;
    return profile;
}

/// Gate: profiling (accumulation + printing) is enabled iff GISMO_C1SURF_PROFILE
/// is set in the environment. Checked once (function-local static).
inline bool gsC1SurfProfilingEnabled()
{
    static const bool enabled = (std::getenv("GISMO_C1SURF_PROFILE") != nullptr);
    return enabled;
}

template<class T, class Visitor = gsC1SurfGluingDataVisitor<T>>
class gsC1SurfGluingData : public gsC1SurfGD<T>
{

public:
    gsC1SurfGluingData()
    {
        setGDEdge();
    }

    gsC1SurfGluingData(gsMultiPatch<T> const & mp,
    gsMultiBasis<T> & mb)
    :  gsC1SurfGD<T>(mp, mb)
    {
        const bool profile = gsC1SurfProfilingEnabled();
        gsStopwatch clock;
        if (profile) clock.restart();

        // Solve the system for alpha_L and alpha_R and beta
        refresh();
        assemble();
        solve();

        // Solve the system for beta_L and beta_R
        refreshBeta();
        assembleBeta();
        solveBeta();

        if (profile)
        {
            gsC1SurfProfile & prof = gsC1SurfGetProfile();
            prof.t_gluingData += clock.stop();
            ++prof.n_gluingData;
        }
    }

    gsMatrix<T> evalAlpha_R(gsMatrix<T> points)
    {
        gsMatrix<T> ones(1, points.cols());
        ones.setOnes();
        return sol.row(0) * ( ones - points ) + sol.row(1) * points;
    }

    gsMatrix<T> evalAlpha_L(gsMatrix<T> points)
    {
        gsMatrix<T> ones(1, points.cols());
        ones.setOnes();
        return sol.row(2) * ( ones - points ) + sol.row(3) * points;
    }

    gsMatrix<T> evalBeta_R(gsMatrix<T> points)
    {
        gsMatrix<T> ones(1, points.cols());
        ones.setOnes();
        return solBeta.row(0) * ( ones - points ) + solBeta.row(1) * points;
    }

    gsMatrix<T> evalBeta_L(gsMatrix<T> points)
    {
        gsMatrix<T> ones(1, points.cols());
        ones.setOnes();
        return solBeta.row(2) * ( ones - points ) + solBeta.row(3) * points;
    }

    gsMatrix<T> evalBeta(gsMatrix<T> points)
    {
        gsMatrix<T> ones(1, points.cols());
        ones.setOnes();
        return sol.row(4) * ( ones - points ).cwiseProduct( ones - points)
            + sol.row(5) * ( ones - points ).cwiseProduct(points) + sol.row(6) * points;
    }

    gsMatrix<T> getSol(){ return sol; }

    gsMatrix<T> getSolBeta(){ return solBeta; }



protected:

    gsSparseSystem<T> mSys;
    gsSparseSystem<T> mSysBeta;
    gsMatrix<T> dirichletDofs;
    gsMatrix<T> dirichletDofsBeta;

    gsMatrix<T> sol; // In order, it contains: alpha_0L, alpha_1L, alpha_0R, alpha_1R, beta_0, beta_1, beta_2
                    // to construct the linear combination of the GD:
                    // alpha_L = ( 1 - t ) * alpha_0L + alpha_1L * t
                    // alpha_R = ( 1 - t ) * alpha_0R + alpha_1R * t
                    //beta = ( 1 - t )^2 * beta_0 + 2 * t * ( 1 - t ) * beta_1 + t^2 * beta_2

    gsMatrix<T> solBeta;

    void refresh()
    {
        gsVector<index_t> size(1);
        size << 7;

        gsDofMapper map(size);
        map.finalize();

        gsSparseSystem<T> sys(map);
        mSys = sys;
    }

    void refreshBeta()
    {
        gsVector<index_t> size(1);
        size << 4;

        gsDofMapper mapBeta(size);
        mapBeta.finalize();

        gsSparseSystem<T> sysBeta(mapBeta);
        mSysBeta = sysBeta;
    }

    void assemble()
    {
        mSys.reserve(49, 1); // Reserve for the matrix 7x7 values

        dirichletDofs.setZero(mSys.colMapper(0).boundarySize(),1);

        // Assemble volume integrals
        Visitor visitor;
        apply(visitor);

        mSys.matrix().makeCompressed();
    }

    void assembleBeta()
    {
        mSysBeta.reserve(16, 1); // Reserve for the matrix 4x4 values

        dirichletDofsBeta.setZero(mSysBeta.colMapper(0).boundarySize(), 1);

        // Assemble volume integrals
        Visitor visitorBeta;
        applyBeta(visitorBeta);

        mSysBeta.matrix().makeCompressed();
    }

    void apply(Visitor visitor)
    {
    #pragma omp parallel
        {
            Visitor
    #ifdef _OPENMP
            // Create thread-private visitor
            visitor_(visitor);
            const int tid = omp_get_thread_num();
            const int nt  = omp_get_num_threads();
    #else
                    &visitor_ = visitor;
    #endif
            gsQuadRule<T> quRule ; // Quadrature rule
            gsMatrix<T> quNodes  ; // Temp variable for mapped nodes
            gsVector<T> quWeights; // Temp variable for mapped weights

            const gsBasis<T> & basis = this->m_mb[0].basis(0).component(1); // = 0

            // Initialize reference quadrature rule and visitor data
            visitor_.initialize(basis,quRule);

            // Initialize domain element iterator
            typename gsBasis<T>::domainIter domIt    = basis.domain()->beginAll();
            typename gsBasis<T>::domainIter domItEnd = basis.domain()->endAll();

#           ifdef _OPENMP
            domIt += tid;
            for ( ; domIt<domItEnd; domIt+=nt )
#           else
            for (; domIt<domItEnd; ++domIt )
#           endif
            {
                // Map the Quadrature rule to the element
                quRule.mapTo( domIt.lowerCorner(), domIt.upperCorner(), quNodes, quWeights );

                // Perform required evaluations on the quadrature nodes
                visitor_.evaluate(quNodes, this->m_mp);

                // Assemble on element
                visitor_.assemble(domIt, quWeights);

                // Push to global matrix and right-hand side vector
    #pragma omp critical(localToGlobal)
                visitor_.localToGlobal( dirichletDofs, mSys); // omp_locks inside
            }
        }//omp parallel
    }

    void applyBeta(Visitor visitorBeta)
    {
    #pragma omp parallel
        {
            Visitor
    #ifdef _OPENMP
            // Create thread-private visitor
            visitor_Beta(visitorBeta);
            const int tid = omp_get_thread_num();
            const int nt  = omp_get_num_threads();
    #else
                    &visitor_Beta = visitorBeta;
    #endif
            gsQuadRule<T> quRule ; // Quadrature rule
            gsMatrix<T> quNodes  ; // Temp variable for mapped nodes
            gsVector<T> quWeights; // Temp variable for mapped weights

            const gsBasis<T> & basis = this->m_mb[0].basis(0).component(1); // = 0

            // Initialize reference quadrature rule and visitor data
            visitor_Beta.initialize(basis,quRule);

            // Initialize domain element iterator
            typename gsBasis<T>::domainIter domIt    = basis.domain()->beginAll();
            typename gsBasis<T>::domainIter domItEnd = basis.domain()->endAll();

#           ifdef _OPENMP
            domIt += tid;
            for ( ; domIt<domItEnd; domIt+=nt )
#           else
            for (; domIt<domItEnd; ++domIt )
#           endif
            {
                // Map the Quadrature rule to the element
                quRule.mapTo( domIt.lowerCorner(), domIt.upperCorner(), quNodes, quWeights );

                // Perform required evaluations on the quadrature nodes
                visitor_Beta.evaluateBeta(quNodes, this->m_mp, sol);

                // Assemble on element
                visitor_Beta.assembleBeta(domIt, quWeights);

                // Push to global matrix and right-hand side vector
    #pragma omp critical(localToGlobal)
                visitor_Beta.localToGlobalBeta( dirichletDofsBeta, mSysBeta); // omp_locks inside
            }
        }//omp parallel
    }



    void solve()
    {
        typename gsSparseSolver<T>::SimplicialLDLT solver;

        solver.compute(mSys.matrix());
        sol = solver.solve(mSys.rhs()); // My solution
    }

    void solveBeta()
    {
        typename gsSparseSolver<T>::SimplicialLDLT solver;

        solver.compute(mSysBeta.matrix());
        solBeta = solver.solve(mSysBeta.rhs()); // My solution
    }

    void setGDEdge()
    {
        gsMatrix<T> solTMP(7, 1);
        gsMatrix<T> solBetaTMP(4, 1);

        solTMP(0, 0) = 1;
        solTMP(1, 0) = 1;
        solTMP(2, 0) = 1;
        solTMP(3, 0) = 1;
        solTMP(4, 0) = 0;
        solTMP(5, 0) = 0;
        solTMP(6, 0) = 0;

        solBetaTMP(0, 0) = 0;
        solBetaTMP(1, 0) = 0;
        solBetaTMP(2, 0) = 0;
        solBetaTMP(3, 0) = 0;

        sol = solTMP;
        solBeta = solBetaTMP;
    }
}; // class gsGluingData



} // namespace gismo

