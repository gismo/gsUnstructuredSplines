/** @file gsApproxC1VertexBasisProjection.h

    @brief Provides assembler for a G1 Basis for multiPatch.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmüller
*/

#pragma once

#include <gsUnstructuredSplines/gsApproxC1VertexBasisProjectionVisitor.h>
#include <gsUnstructuredSplines/gsApproxGluingData.h>

namespace gismo
{

template<short_t d, class T, class bhVisitor = gsApproxC1VertexBasisProjectionVisitor<d, T>>
class gsApproxC1VertexBasisProjection : public gsAssembler<T>
{
private:
    typedef typename std::vector<gsPatchReparameterized<d,T>> C1AuxPatchContainer;

public:
    typedef gsAssembler<T> Base;

public:
    gsApproxC1VertexBasisProjection(C1AuxPatchContainer & auxPatches,
                                     gsApproxGluingData<d, T> & approxGluingData,
                                     const index_t corner,
                                     const std::vector<index_t> & sideContainer,
                                     std::vector<bool> & isInterface,
                                     const real_t sigma,
                                     const gsOptionList & optionList)
            : m_auxPatches(auxPatches), m_approxGluingData(approxGluingData), m_corner(corner), m_sigma(sigma), m_optionList(optionList)
    {
        // Collect the needed basis
        m_basis_g1 =  dynamic_cast<gsTensorBSplineBasis<d, T>&>(m_auxPatches[0].getBasisRotated().getBasis(corner+4));

        m_basis_plus.resize(2);
        m_basis_minus.resize(2);
        m_basis_geo.resize(2);

        kindOfEdge.resize(2);

        for (size_t dir = 0; dir < sideContainer.size(); ++dir)
        {
            index_t localdir = m_auxPatches[0].getMapIndex(sideContainer[dir]) < 3 ? 1 : 0;

            m_basis_plus[localdir] = dynamic_cast<gsBSplineBasis<T>&>(m_auxPatches[0].getBasisRotated().getHelperBasis(sideContainer[dir]-1, 0));
            m_basis_minus[localdir] = dynamic_cast<gsBSplineBasis<T>&>(m_auxPatches[0].getBasisRotated().getHelperBasis(sideContainer[dir]-1, 1));
            m_basis_geo[localdir] = dynamic_cast<gsBSplineBasis<T>&>(m_auxPatches[0].getBasisRotated().getHelperBasis(sideContainer[dir]-1, 2));

            //m_basis_geo[localdir].reduceContinuity(1);
            kindOfEdge[localdir] = isInterface[dir];
            //kindOfEdge[localdir] = m_auxPatches[0].getBasisRotated().isInterface(sideContainer[dir]);

        }
/*
        gsInfo << "Basis plus: " << m_basis_plus[0] << "\n";
        gsInfo << "Basis minus: " << m_basis_minus[0] << "\n";
        gsInfo << "Basis geo: " << m_basis_geo[0] << "\n";
        gsInfo << "Basis plus 1: " << m_basis_plus[1] << "\n";
        gsInfo << "Basis minus 1: " << m_basis_minus[1] << "\n";
        gsInfo << "Basis geo 1: " << m_basis_geo[1] << "\n";

        gsInfo << "basis g1: " << m_basis_g1 << "\n";
        gsInfo << "kindOfEdge " << kindOfEdge[0] << " : " << kindOfEdge[1] << " : " << true << "\n";
*/

    }


    void refresh();
    void assemble();
    inline void apply(bhVisitor & visitor, int patchIndex);
    void solve();

    void constructSolution(gsMultiPatch<T> & result);

    void setBasisVertex(gsMultiPatch<T> & result)
    {
        refresh();
        assemble();
        solve();

        constructSolution(result);
    }


protected:

    // Input
    C1AuxPatchContainer m_auxPatches;
    gsApproxGluingData<d, T> m_approxGluingData;
    index_t m_corner;
    real_t m_sigma;
    const gsOptionList m_optionList;


    // Basis for the G1 Basis
    gsTensorBSplineBasis<d, T> m_basis_g1;

    // Basis need for construction
    std::vector<gsBSplineBasis<T>> m_basis_plus;
    std::vector<gsBSplineBasis<T>> m_basis_minus;
    std::vector<gsBSplineBasis<T>> m_basis_geo;

    // System
    std::vector<gsSparseSystem<T> > m_systemContainer;

    // For gluing data
    std::vector<bool> kindOfEdge;

    // For Dirichlet boundary
    using Base::m_ddof;

    std::vector<gsMatrix<>> solVec;
    }; // class gsG1BasisEdge


template <short_t d, class T, class bhVisitor>
void gsApproxC1VertexBasisProjection<d, T, bhVisitor>::constructSolution(gsMultiPatch<T> & result)
{

    result.clear();

    // Dim is the same for all basis functions
    const index_t dim = ( 0!=solVec.at(0).cols() ? solVec.at(0).cols() :  m_ddof[0].cols() );

    gsMatrix<T> coeffs;
    for (index_t p = 0; p < 6; ++p)
    {

        const gsDofMapper & mapper = m_systemContainer[p].colMapper(0); // unknown = 0

        // Reconstruct solution coefficients on patch p
        index_t sz;
        sz = m_basis_g1.basis(0).size();

        coeffs.resize(sz, dim);

        for (index_t i = 0; i < sz; ++i)
        {
            if (mapper.is_free(i, 0)) // DoF value is in the solVector // 0 = unitPatch
            {
                coeffs.row(i) = solVec.at(p).row(mapper.index(i, 0));
            }
            else // eliminated DoF: fill with Dirichlet data
            {
                //gsInfo << "mapper index dirichlet: " << m_ddof[unk].row( mapper.bindex(i, p) ).head(dim) << "\n";
                coeffs.row(i) = m_ddof[p].row( mapper.bindex(i, 0) ).head(dim); // = 0
            }
        }
        result.addPatch(m_basis_g1.basis(0).makeGeometry(give(coeffs)));
    }
}

template <short_t d, class T, class bhVisitor>
void gsApproxC1VertexBasisProjection<d, T, bhVisitor>::refresh()
{
    // 1. Obtain a map from basis functions to matrix columns and rows
    gsDofMapper map(m_basis_g1);

    index_t p_1 = m_basis_g1.degree(0);
    index_t p_2 = m_basis_g1.degree(1);

    gsMatrix<index_t> act;
    if (3*p_1+1 < m_basis_g1.component(0).size())
        for (index_t i = 3*p_1+1; i < m_basis_g1.component(0).size(); i++) // only the first 2*p+1
        {
            act = m_basis_g1.boundaryOffset(1, i); // WEST
            map.markBoundary(0, act); // Patch 0
        }
    if (3*p_2+1 < m_basis_g1.component(1).size())
        for (index_t i = 3*p_2+1; i < m_basis_g1.component(1).size(); i++) // only the first 2*p+1
        {
            act = m_basis_g1.boundaryOffset(3, i); // WEST
            map.markBoundary(0, act); // Patch 0
        }

    map.finalize();
    //gsInfo << "map : " << map.asVector() << "\n";

    // 2. Create the sparse system
    m_systemContainer.clear();
    gsSparseSystem<T> m_system = gsSparseSystem<T>(map);
    for (index_t i = 0; i < 6; i++)
        m_systemContainer.push_back(m_system);

} // refresh()

template <short_t d, class T, class bhVisitor>
void gsApproxC1VertexBasisProjection<d, T, bhVisitor>::assemble()
{
    // Reserve sparse system
    const index_t nz = gsAssemblerOptions::numColNz(m_basis_g1,2,1,0.333333);
    for (unsigned i = 0; i < m_systemContainer.size(); i++)
        m_systemContainer[i].reserve(nz, 1);


    const gsDofMapper & map = m_systemContainer[0].colMapper(0); // Map same for every

    m_ddof.resize(6); // 0,1
    for (size_t i = 0; i < m_ddof.size(); i++)
        m_ddof[i].setZero(map.boundarySize(), 1 ); // plus


    // Assemble volume integrals
    bhVisitor visitor;
    apply(visitor,0); // patch 0

    for (unsigned i = 0; i < m_systemContainer.size(); i++)
        m_systemContainer[i].matrix().makeCompressed();

} // assemble()

template <short_t d, class T, class bhVisitor>
void gsApproxC1VertexBasisProjection<d, T, bhVisitor>::apply(bhVisitor & visitor, int patchIndex)
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

        gsBasis<T> & basis_g1 = m_basis_g1; // basis for construction

        // Initialize reference quadrature rule and visitor data
        visitor_.initialize(basis_g1, quRule);

        const gsGeometry<T> & patch = m_auxPatches[patchIndex].getPatchRotated(); // patchIndex = 0

        // Initialize domain element iterator
        typename gsBasis<T>::domainIter domIt = basis_g1.makeDomainIterator(boundary::none);

#ifdef _OPENMP
        for ( domIt->next(tid); domIt->good(); domIt->next(nt) )
#else
        for (; domIt->good(); domIt->next() )
#endif
        {
            // Map the Quadrature rule to the element
            quRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights );
#pragma omp critical(evaluate)
            // Perform required evaluations on the quadrature nodes
            visitor_.evaluate(patch, m_basis_g1, m_basis_plus, m_basis_minus, m_basis_geo, m_approxGluingData,
                              quNodes, m_sigma, kindOfEdge, m_optionList);

            // Assemble on element
            visitor_.assemble(*domIt, quWeights);

            // Push to global matrix and right-hand side vector
#pragma omp critical(localToGlobal)
            visitor_.localToGlobal(patchIndex, m_ddof, m_systemContainer); // omp_locks inside // patchIndex == 0
        }
    }//omp parallel
} // apply

template <short_t d, class T, class bhVisitor>
void gsApproxC1VertexBasisProjection<d, T, bhVisitor>::solve()
{
    gsSparseSolver<real_t>::LU solver;

    for (index_t i = 0; i < 6; i++) // Tilde
    {
        solver.compute(m_systemContainer.at(i).matrix());
        solVec.push_back(solver.solve(m_systemContainer.at(i).rhs()));
    }
} // solve

} // namespace gismo