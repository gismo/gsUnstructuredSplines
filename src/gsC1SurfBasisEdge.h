/** @file gsC1SurfEdge.h

    @brief Creates the (approx) C1 Edge space.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Farahat & P. Weinmueller
*/

#pragma once

#include <gsUnstructuredSplines/src/gsC1SurfGluingData.h>
#include <gsUnstructuredSplines/src/gsC1SurfVisitorBasisEdge.h>

#include <gsAssembler/gsDofMapperCreator.h>

namespace gismo
{
    template<class T, class bhVisitor = gsC1SurfVisitorBasisEdge<T>>
    class gsC1SurfBasisEdge : public gsAssembler<T>
    {
    public:
        typedef gsAssembler<T> Base;

    public:
        gsC1SurfBasisEdge(const gsMultiPatch<T> & mp, // single patch
                          const gsMultiBasis<T> & basis, // single basis
                        index_t uv, // !!! 0 == u; 1 == v !!!
                        bool isBoundary,
                        gsC1SurfGluingData<T> gluingD)
                : m_mp(mp), m_basis(basis), m_uv(uv), m_isBoundary(isBoundary), m_gD(gluingD)
        {

            // Computing the G1 - basis function at the edge
            gsBSplineBasis<T> basis_edge = dynamic_cast<gsBSplineBasis<T> &>(m_basis.basis(0).component(m_uv)); // 0 -> v, 1 -> u

            gsBSplineBasis<T> basis_plus(basis_edge);
            basis_plus.elevateContinuity(1);
            m_basis_plus = basis_plus;
            //gsDebugVar(basis_plus.knots().asMatrix());


            gsBSplineBasis<T> basis_minus(basis_edge);
            basis_minus.degreeReduce(1);
            m_basis_minus = basis_minus;
            //gsDebugVar(basis_minus.knots().asMatrix());

            // Basis for the G1 basis
            m_basis_g1 = m_basis.basis(0);

//        gsTensorBSplineBasis<2, T> basis_edge_ab = dynamic_cast<gsTensorBSplineBasis<2, T> &>(m_basis_g1.basis(0));
//        gsInfo << "Basis edge 0: " << basis_edge_ab.component(0).knots().asMatrix() << "\n";
//        gsInfo << "Basis edge 1: " << basis_edge_ab.component(1).knots().asMatrix() << "\n";
        }

        // Computed the gluing data globally
        void setG1BasisEdge(gsMultiPatch<T> & result);

        void refresh();
        void assemble(index_t plus_lo, index_t plus_hi, index_t minus_lo, index_t minus_hi);
        inline void apply(bhVisitor & visitor, index_t plus_lo, index_t plus_hi,
                          index_t minus_lo, index_t minus_hi);
        void solve();

        using Base::constructSolution;
        void constructSolution(const gsMatrix<T> & solVector,
                               gsMultiPatch<T> & result, short_t unk = 0) const;

    private:
        // Avoid hidden overloads w.r.t. gsAssembler
        void assemble()
        { GISMO_NO_IMPLEMENTATION; }

        using gsAssembler<T>::assemble;
        void assemble(const gsMultiPatch<T> & /* curSolution */)
        { GISMO_NO_IMPLEMENTATION; }

    protected:

        // Input
        gsMultiPatch<T> m_mp;
        gsMultiBasis<T> m_basis;
        index_t m_uv;
        bool m_isBoundary;

        // Gluing data
        gsC1SurfGluingData<T> m_gD;

        // Basis for getting the G1 Basis
        gsBSplineBasis<T> m_basis_plus;
        gsBSplineBasis<T> m_basis_minus;

        // Basis for the G1 Basis
        gsMultiBasis<T> m_basis_g1;

        // Basis for Integration
        gsMultiBasis<T> m_geo;

        // For Dirichlet boundary
        using Base::m_ddof;
        using Base::m_system;


    }; // class gsG1BasisEdge

    template <class T, class bhVisitor>
    void gsC1SurfBasisEdge<T,bhVisitor>::setG1BasisEdge(gsMultiPatch<T> & result)
    {
        result.clear();

        const index_t n_plus  = m_basis_plus.size();
        const index_t n_minus = m_basis_minus.size();

        // first 3 and last 3 "plus", first 2 and last 2 "minus" basis functions are eliminated
        const index_t plus_lo = 3, plus_hi = n_plus - 3;
        const index_t minus_lo = 2, minus_hi = n_minus - 2;
        const index_t n_p = math::max(plus_hi  - plus_lo,  (index_t)0);
        const index_t n_m = math::max(minus_hi - minus_lo, (index_t)0);

        gsMultiPatch<T> g1EdgeBasis;
        if (n_p + n_m > 0)
        {
            m_geo = m_basis_g1; // Basis for Integration
            refresh();
            assemble(plus_lo, plus_lo + n_p, minus_lo, minus_lo + n_m);

            typename gsSparseSolver<T>::SimplicialLDLT solver;
            solver.compute(m_system.matrix());
            gsMatrix<T> sol = solver.solve(m_system.rhs());   // one column per basis function

            constructSolution(sol, g1EdgeBasis);
        }
        result = g1EdgeBasis;
    } // setG1BasisEdge

    template <class T, class bhVisitor>
    void gsC1SurfBasisEdge<T,bhVisitor>::constructSolution(const gsMatrix<T> & solVector, gsMultiPatch<T> & result, short_t unk) const
    {
        GISMO_UNUSED(unk);

        // solVector holds one column per basis function of this patch side; each
        // column becomes its own single-coefficient-column patch, appended in
        // column order (all "plus" ascending, then all "minus" ascending) --
        // gsC1SurfSpline::compute() consumes these patches by that fixed index.
        const gsDofMapper & mapper = m_system.colMapper(0); // unknown = 0
        const index_t sz = m_basis.basis(0).size();

        for (index_t col = 0; col != solVector.cols(); ++col)
        {
            gsMatrix<T> coeffs(sz, 1);
            for (index_t i = 0; i < sz; ++i)
            {
                if (mapper.is_free(i, 0)) // DoF value is in the solVector // 0 = unitPatch
                    coeffs(i, 0) = solVector(mapper.index(i, 0), col);
                else // eliminated DoF: fill with Dirichlet data
                    coeffs(i, 0) = m_ddof[0](mapper.bindex(i, 0), col); // = 0
            }
            result.addPatch(m_basis_g1.basis(0).makeGeometry(give(coeffs)));
        }
    }

    template <class T, class bhVisitor>
    void gsC1SurfBasisEdge<T,bhVisitor>::refresh()
    {
        // 1. Obtain a map from basis functions to matrix columns and rows
        gsDofMapper map = createMapper(m_basis.basis(0), 1, false);

        gsMatrix<index_t> act;

        for (index_t i = 2; i < m_basis.basis(0).component(1-m_uv).size(); i++) // only the first two u/v-columns are Dofs (0/1)
        {
            act = m_basis.basis(0).boundaryOffset(m_uv == 0 ? 3 : 1, i); // WEST
            map.markBoundary(0, act); // Patch 0
        }

        map.finalize();

        // 2. Create the sparse system
        m_system = gsSparseSystem<T>(map);

    } // refresh()

    template <class T, class bhVisitor>
    void gsC1SurfBasisEdge<T,bhVisitor>::assemble(index_t plus_lo, index_t plus_hi, index_t minus_lo, index_t minus_hi)
    {
        const index_t M = (plus_hi - plus_lo) + (minus_hi - minus_lo);

        // Reserve sparse system
        const index_t nz = gsAssemblerOptions::numColNz(m_basis[0],2,1,0.333333);
        m_system.reserve(nz, M);

        if(m_ddof.size()==0)
            m_ddof.resize(1); // 0,1

        const gsDofMapper & map = m_system.colMapper(0); // Map same for every functions

        m_ddof[0].setZero(map.boundarySize(), M );

        // Assemble volume integrals
        bhVisitor visitor;
        apply(visitor, plus_lo, plus_hi, minus_lo, minus_hi);

        m_system.matrix().makeCompressed();

    } // assemble()

    template <class T, class bhVisitor>
    void gsC1SurfBasisEdge<T,bhVisitor>::apply(bhVisitor & visitor, index_t plus_lo, index_t plus_hi,
                                               index_t minus_lo, index_t minus_hi)
    {
#pragma omp parallel
        {

            gsQuadRule<T> quRule ; // Quadrature rule
            gsMatrix<T> quNodes  ; // Temp variable for mapped nodes
            gsVector<T> quWeights; // Temp variable for mapped weights

            bhVisitor
#ifdef _OPENMP
            // Create thread-private visitor
    visitor_(visitor);
    const int tid = omp_get_thread_num();
    const int nt  = omp_get_num_threads();
#else
                    &visitor_ = visitor;
#endif

            gsBasis<T> & basis_g1 = m_basis_g1.basis(0); // basis for construction

            // Same for all patches
            gsBasis<T> & basis_geo = m_basis.basis(0).component(1-m_uv);
            gsBasis<T> & basis_plus = m_basis_plus;
            gsBasis<T> & basis_minus = m_basis_minus;

            // Initialize reference quadrature rule and visitor data
            visitor_.initialize(basis_g1, quRule);

            const gsGeometry<T> & patch = m_mp.patch(0);

            // refresh() eliminates every DoF whose transverse (1-m_uv) tensor
            // index is >= 2, so the only free DoFs are transverse functions 0
            // and 1. A 2-D tensor basis function is active on an element only
            // if both its 1-D factors are, so a free DoF is active on an
            // element only if the element's transverse index is within
            // function 1's support -- a condition on the transverse element
            // coordinate alone. basis_geo and the domain iterator below both
            // come from m_basis.basis(0) (m_geo == m_basis_g1 ==
            // m_basis.basis(0)), so cutoff and element boundaries are values
            // from the same knot vector.
            const index_t transDir = 1 - m_uv;
            const bool canSkip = basis_geo.size() > 1;
            const T cutoff = basis_geo.support(canSkip ? 1 : 0)(0, 1);
            // Element boundaries and cutoff are both knot values; the
            // tolerance only has to absorb floating-point roundoff at exact
            // equality. A 1e-12 relative fraction of the parameter-domain
            // width stays orders of magnitude below the smallest element
            // width at any refinement level.
            const gsMatrix<T> domainSupport = basis_geo.support();
            const T tol = (domainSupport(0, 1) - domainSupport(0, 0)) * (T)1e-12;

            // Initialize domain element iterator
            typename gsBasis<T>::domainIter domIt    = m_geo.basis(0).domain()->beginAll();
            typename gsBasis<T>::domainIter domItEnd = m_geo.basis(0).domain()->endAll();

#           ifdef _OPENMP
            domIt += tid;
            for ( ; domIt<domItEnd; domIt+=nt )
#           else
            for (; domIt<domItEnd; ++domIt )
#           endif
            {
                // No free DoF is active on this element: rhsVals and every
                // localMat(i,j) contribution for it are already discarded by
                // push() below, so skipping it changes no assembled value.
                if (canSkip && domIt.lowerCorner()(transDir) >= cutoff - tol) continue;

                // Map the Quadrature rule to the element
                quRule.mapTo( domIt.lowerCorner(), domIt.upperCorner(), quNodes, quWeights );

                // Perform required evaluations on the quadrature nodes
                visitor_.evaluate(plus_lo, plus_hi, minus_lo, minus_hi, basis_g1, basis_geo, basis_plus, basis_minus, patch, quNodes, m_uv, m_gD, m_isBoundary);

                // Assemble on element
                visitor_.assemble(domIt, quWeights);

                // Push to global matrix and right-hand side vector
#pragma omp critical(localToGlobal)
                visitor_.localToGlobal(0, m_ddof, m_system); // omp_locks inside // patchIndex == 0
            }
        }//omp parallel
    } // apply

} // namespace gismo
