/** @file gsC1SurfVertex.h

    @brief Creates the (approx.) C1 Vertex space.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Farahat & P. Weinmueller
*/

#pragma once


#include <gsUnstructuredSplines/src/gsC1SurfGluingData.h>
#include <gsUnstructuredSplines/src/gsC1SurfVisitorBasisVertex.h>

#include <gsAssembler/gsDofMapperCreator.h>



namespace gismo
{
    template<class T, class bhVisitor = gsG1ASVisitorBasisVertex<T>>
    class gsC1SurfBasisVertex : public gsAssembler<T>
    {
    public:
        typedef gsAssembler<T> Base;

    public:
        gsC1SurfBasisVertex(gsMultiPatch<T> mp, // Single Patch
                          gsMultiBasis<T> basis, // Single Basis
                          std::vector<bool> isBoundary,
                          gsMatrix<T> &Phi,
                          gsMatrix<T> gluingD)
                : m_mp(mp), m_basis(basis), m_isBoundary(isBoundary), m_Phi(Phi), m_gD(gluingD)
        {

            for (index_t dir = 0; dir < m_mp.parDim(); dir++) // For the TWO directions
            {
                // Computing the G1 - basis function at the edge
                // Spaces for computing the g1 basis
                gsBSplineBasis<T> basis_edge = dynamic_cast<gsBSplineBasis<T> &>(m_basis.basis(0).component(dir)); // 0 -> u, 1 -> v

                gsBSplineBasis<T> basis_plus(basis_edge);
                basis_plus.elevateContinuity(1);
                m_basis_plus.push_back(basis_plus);

                gsBSplineBasis<T> basis_minus(basis_edge);
                basis_minus.degreeReduce(1);
                m_basis_minus.push_back(basis_minus);
            }

            // Basis for the G1 basis
            m_basis_g1 = m_basis.basis(0);


        }


        void refresh();
        using Base::assemble;
        void assemble();
        inline void apply(bhVisitor & visitor, index_t patchIndex);
        void solve();

        void constructSolution(gsMultiPatch<T> & result);

        void setG1BasisVertex(gsMultiPatch<T> & result)
        {
            m_geo = m_basis_g1;
            refresh();
            assemble();
            solve();

            constructSolution(result);

            if (gsC1SurfProfilingEnabled())
                gsC1SurfGetProfile().print(gsInfo); // cumulative report after this vertex
        }

    private:
        // Avoid hidden overloads w.r.t. gsAssembler
        void constructSolution(const gsMatrix<T>& /* solVector */,
                               gsMultiPatch<T>& /* result */, short_t /* unk */) const
        { GISMO_NO_IMPLEMENTATION; }

        void constructSolution(const gsMatrix<T>& /* solVector */,
                               gsMultiPatch<T>& /* result */,
                               const gsVector<index_t>  & /* unknowns */) const
        { GISMO_NO_IMPLEMENTATION; }

    protected:

        // Input
        gsMultiPatch<T> m_mp;
        gsMultiBasis<T> m_basis;
        std::vector<bool> m_isBoundary;
        gsMatrix<T> m_Phi;

        // Gluing data
        gsMatrix<T> m_gD;

        // Basis for getting the G1 Basis
        std::vector<gsBSplineBasis<T>> m_basis_plus;
        std::vector<gsBSplineBasis<T>> m_basis_minus;

        // Basis for the G1 Basis
        gsMultiBasis<T> m_basis_g1;

        // Basis for Integration
        gsMultiBasis<T> m_geo;

        // System
        using Base::m_system;

        // For Dirichlet boundary
        using Base::m_ddof;

        gsMatrix<T> m_solMat;

    }; // class gsG1BasisEdge


    template <class T, class bhVisitor>
    void gsC1SurfBasisVertex<T,bhVisitor>::constructSolution(gsMultiPatch<T> & result)
    {

        result.clear();

        const gsDofMapper & mapper = m_system.colMapper(0); // unknown = 0, same mapper for all six patches

        // Reconstruct solution coefficients on patch p
        const index_t sz = m_basis.basis(0).size();

        gsMatrix<T> coeffs;
        for (index_t p = 0; p < 6; ++p)
        {
            coeffs.resize(sz, 1);

            for (index_t i = 0; i < sz; ++i)
            {
                if (mapper.is_free(i, 0)) // DoF value is in the solVector // 0 = unitPatch
                    coeffs(i, 0) = m_solMat(mapper.index(i, 0), p);
                else // eliminated DoF: fill with Dirichlet data
                    coeffs(i, 0) = m_ddof[0](mapper.bindex(i, 0), p); // = 0
            }
            result.addPatch(m_basis_g1.basis(0).makeGeometry(give(coeffs)));
        }
    }

    template <class T, class bhVisitor>
    void gsC1SurfBasisVertex<T,bhVisitor>::refresh()
    {
        // 1. Obtain a map from basis functions to matrix columns and rows
        gsDofMapper map = createMapper(m_basis.basis(0), 1, false);

        // SET THE DOFS

        map.finalize();


        // 2. Create the sparse system
        m_system = gsSparseSystem<T>(map);

    } // refresh()

    template <class T, class bhVisitor>
    void gsC1SurfBasisVertex<T,bhVisitor>::assemble()
    {
        // Reserve sparse system
        const index_t nz = gsAssemblerOptions::numColNz(m_basis[0],2,1,0.333333);
        m_system.reserve(nz, 6);

        if(m_ddof.size()==0)
            m_ddof.resize(1); // 0,1

        const gsDofMapper & map = m_system.colMapper(0); // Map same for every

        m_ddof[0].setZero(map.boundarySize(), 6 ); // plus

        // Assemble volume integrals
        bhVisitor visitor;
        apply(visitor,0); // patch 0

        m_system.matrix().makeCompressed();

    } // assemble()

    template <class T, class bhVisitor>
    void gsC1SurfBasisVertex<T,bhVisitor>::apply(bhVisitor & visitor, index_t patchIndex)
    {
        const bool profile = gsC1SurfProfilingEnabled();
        gsStopwatch applyClock;
        if (profile) applyClock.restart();

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
            gsBasis<T> & basis_geo = m_basis.basis(0); // 2D

            // Initialize reference quadrature rule and visitor data
            visitor_.initialize(basis_g1, quRule);

            const gsGeometry<T> & patch = m_mp.patch(0);

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
                // Map the Quadrature rule to the element
                quRule.mapTo( domIt.lowerCorner(), domIt.upperCorner(), quNodes, quWeights );
#pragma omp critical(evaluate)
                // Perform required evaluations on the quadrature nodes
                visitor_.evaluate(basis_g1, basis_geo, m_basis_plus, m_basis_minus, patch, quNodes, m_gD, m_isBoundary, m_Phi);

                // Assemble on element
                visitor_.assemble(domIt, quWeights);

                // Push to global matrix and right-hand side vector
#pragma omp critical(localToGlobal)
                visitor_.localToGlobal(patchIndex, m_ddof, m_system); // omp_locks inside // patchIndex == 0
            }
        }//omp parallel

        if (profile)
        {
            gsC1SurfProfile & prof = gsC1SurfGetProfile();
            prof.t_vertexApply += applyClock.stop();
            ++prof.n_vertexApplyCalls;
            // Unlike gsC1SurfBasisEdge, this is ONE full 2D pass that builds the
            // single six-column system for all six vertex basis functions together.
            prof.n_vertexElemVisits += (index_t)m_geo.basis(0).numElements();
        }
    } // apply

    template <class T, class bhVisitor>
    void gsC1SurfBasisVertex<T,bhVisitor>::solve()
    {
        const bool profile = gsC1SurfProfilingEnabled();
        gsStopwatch clock;
        if (profile) clock.restart();

        typename gsSparseSolver<T>::SimplicialLDLT solver;
//    typename gsSparseSolver<T>::LU solver;

        solver.compute(m_system.matrix());
        m_solMat = solver.solve(m_system.rhs());   // (free dofs) x 6

        if (profile)
            gsC1SurfGetProfile().t_vertexSolve += clock.stop();
    } // solve

} // namespace gismo
