/** @file gsApproxC1Edge.h

    @brief Creates the (approx) C1 Edge space.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/

#pragma once

#include <gsUnstructuredSplines/src/gsContainerBasis.h>
#include <gsUnstructuredSplines/src/gsPatchReparameterized.h>

#include <gsUnstructuredSplines/src/gsApproxC1GluingData.h>

#include <gsUnstructuredSplines/src/gsApproxC1Utils.h>

namespace gismo
{

/** @brief Batched sibling of gsTraceBasis (gsApproxC1Utils.h).

    gsApproxC1Edge::compute() needs the RHS of `u * aa` for every trace-basis
    function bfID in a range; the test space `u`, the quadrature nodes and the
    mass-matrix factorization are the same for all of them, so a separate
    quadrature pass per bfID is redundant work.

    This class evaluates the WHOLE range [bfID_lo, bfID_hi) of trace-basis
    functions in one call, one row per basis function, so gsExprAssembler can
    assemble a multi-column RHS (`A.initVector(M); A.assemble(u * aa.tr())`) in a
    single quadrature pass instead of M. `beta`, `N_0` and `N_1` (the u.row(1-uv)
    evaluations) are independent of bfID and are hoisted out of the per-row loop;
    only `m_basis_plus.evalSingle_into`/`derivSingle_into` (cheap 1-D evals) are
    repeated per row. The arithmetic per row is identical to gsTraceBasis::eval_into.
*/
template <class T>
class gsTraceBasisBatch : public gismo::gsFunction<T>
{
protected:
    gsGeometry<T> & _geo;

    gsBSpline<T>       m_basis_beta;
    gsBSplineBasis<T>  m_basis_plus;

    gsBasis<T> & m_basis;

    bool m_isboundary;
    const index_t m_bfID_lo, m_bfID_hi, m_uv;

public:
    gsTraceBasisBatch(gsGeometry<T> & geo,
                       gsBSpline<T> basis_beta,
                       gsBSplineBasis<T> basis_plus,
                       gsBasis<T> & basis,
                       bool isboundary,
                       index_t bfID_lo,
                       index_t bfID_hi,
                       index_t uv) :
            _geo(geo), m_basis_beta(basis_beta), m_basis_plus(basis_plus), m_basis(basis),
            m_isboundary(isboundary), m_bfID_lo(bfID_lo), m_bfID_hi(bfID_hi), m_uv(uv),
            _piece(nullptr)
    { }

    ~gsTraceBasisBatch() { delete _piece; }

    GISMO_CLONE_FUNCTION(gsTraceBasisBatch)

    short_t domainDim() const {return 2;}

    short_t targetDim() const {return m_bfID_hi - m_bfID_lo;}

    mutable gsTraceBasisBatch<T> * _piece;

    const gsFunction<T> & piece(const index_t k) const
    {
        GISMO_UNUSED(k);
        _piece = new gsTraceBasisBatch(*this);
        return *_piece;
    }

    // Input is parametric coordinates of 2-D \a mp; one output row per basis function
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize( targetDim(), u.cols() );

        // tau/p (bfID-independent)
        gsBSplineBasis<T> bsp_temp = dynamic_cast<gsBSplineBasis<T> & >(m_basis.component(1-m_uv));
        T p = bsp_temp.degree();
        T tau_1 = bsp_temp.knots().at(p + 1); // p + 2

        // beta, N_0, N_1: bfID-independent, hoisted out of the per-row loop
        gsMatrix<T> beta, N_0, N_1;
        if (!m_isboundary)
            m_basis_beta.eval_into(u.row(m_uv),beta); // 1-dir == PatchID
        else
            beta.setZero(1, u.cols());

        m_basis.component(1-m_uv).evalSingle_into(0,u.row(1-m_uv),N_0); // u
        m_basis.component(1-m_uv).evalSingle_into(1,u.row(1-m_uv),N_1); // u

        gsMatrix<T> N_i_plus, der_N_i_plus, temp;
        for (index_t bfID = m_bfID_lo; bfID < m_bfID_hi; ++bfID)
        {
            m_basis_plus.evalSingle_into(bfID,u.row(m_uv),N_i_plus); // v
            m_basis_plus.derivSingle_into(bfID,u.row(m_uv),der_N_i_plus);

            temp = beta.cwiseProduct(der_N_i_plus);
            result.row(bfID - m_bfID_lo) = (N_i_plus.cwiseProduct(N_0 + N_1) - temp.cwiseProduct(N_1) * tau_1 / p).row(0);
        }
    }
}; // class gsTraceBasisBatch

/// Batched sibling of gsNormalDerivBasis (gsApproxC1Utils.h) -- see gsTraceBasisBatch.
template <class T>
class gsNormalDerivBasisBatch : public gismo::gsFunction<T>
{
protected:
    gsGeometry<T> & _geo;

    gsBSpline<T>       m_basis_alpha;
    gsBSplineBasis<T>  m_basis_minus;

    gsBasis<T> & m_basis;

    bool m_isboundary;
    const index_t m_bfID_lo, m_bfID_hi, m_uv;

public:
    gsNormalDerivBasisBatch(gsGeometry<T> & geo,
                             gsBSpline<T> basis_alpha,
                             gsBSplineBasis<T> basis_minus,
                             gsBasis<T> & basis,
                             bool isboundary,
                             index_t bfID_lo,
                             index_t bfID_hi,
                             index_t uv) :
            _geo(geo), m_basis_alpha(basis_alpha), m_basis_minus(basis_minus), m_basis(basis),
            m_isboundary(isboundary), m_bfID_lo(bfID_lo), m_bfID_hi(bfID_hi), m_uv(uv),
            _piece(nullptr)
    { }

    ~gsNormalDerivBasisBatch() { delete _piece; }

    GISMO_CLONE_FUNCTION(gsNormalDerivBasisBatch)

    short_t domainDim() const {return 2;}

    short_t targetDim() const {return m_bfID_hi - m_bfID_lo;}

    mutable gsNormalDerivBasisBatch<T> * _piece;

    const gsFunction<T> & piece(const index_t k) const
    {
        GISMO_UNUSED(k);
        _piece = new gsNormalDerivBasisBatch(*this);
        return *_piece;
    }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize( targetDim(), u.cols() );

        gsBSplineBasis<T> bsp_temp = dynamic_cast<gsBSplineBasis<T> & >(m_basis.component(1-m_uv));
        T p = bsp_temp.degree();
        T tau_1 = bsp_temp.knots().at(p + 1); // p + 2

        gsMatrix<T> alpha, N_1;
        if (!m_isboundary)
            m_basis_alpha.eval_into(u.row(m_uv),alpha); // 1-dir == PatchID
        else
            alpha.setOnes(1, u.cols());

        m_basis.component(1-m_uv).evalSingle_into(1,u.row(1-m_uv),N_1); // u

        gsMatrix<T> N_j_minus;
        for (index_t bfID = m_bfID_lo; bfID < m_bfID_hi; ++bfID)
        {
            m_basis_minus.evalSingle_into(bfID,u.row(m_uv),N_j_minus); // v

            if (!m_isboundary)
                result.row(bfID - m_bfID_lo) = ((m_uv == 0 ? T(-1.0) : T(1.0)) * alpha.cwiseProduct(N_j_minus.cwiseProduct(N_1)) * tau_1 / p).row(0);
            else
                result.row(bfID - m_bfID_lo) = ((m_uv == 0 ? T(-1.0) : T(1.0)) * alpha.cwiseProduct(N_j_minus.cwiseProduct(N_1))).row(0);
        }
    }
}; // class gsNormalDerivBasisBatch


template<short_t d, class T>
class gsApproxC1Edge
{

private:
    typedef gsContainerBasis<d, T> Basis;
    typedef typename std::vector<Basis> BasisContainer;
    typedef typename std::vector<gsPatchReparameterized<d,T>> C1AuxPatchContainer;

    /// Shared pointer for gsApproxC1Edge
    typedef memory::shared_ptr<gsApproxC1Edge> Ptr;

    /// Unique pointer for gsApproxC1Edge
    typedef memory::unique_ptr<gsApproxC1Edge> uPtr;


public:
    /// Empty constructor
    ~gsApproxC1Edge() { }


    gsApproxC1Edge(gsMultiPatch<T> const & mp,
                   BasisContainer & bases,
                const boundaryInterface & item,
                size_t & numInt,
                const gsOptionList & optionList)
                : m_mp(mp), m_bases(bases), m_optionList(optionList)
    {
        GISMO_UNUSED(numInt);

        m_auxPatches.clear();
        m_auxPatches.push_back(gsPatchReparameterized<d,T>(m_mp.patch(item.first().patch), m_bases[item.first().patch]));
        m_auxPatches.push_back(gsPatchReparameterized<d,T>(m_mp.patch(item.second().patch), m_bases[item.second().patch]));

        std::vector<patchSide> sidesContainer(2);
        sidesContainer[0] = item.first();
        sidesContainer[1] = item.second();

        reparametrizeInterfacePatches(sidesContainer);

        compute(sidesContainer);
/*
        if (m_optionList.getSwitch("plot"))
        {
            std::string fileName;
            std::string basename = "InterfaceBasisFunctions" + util::to_string(numInt);
            gsParaviewCollection collection(basename);

            for (size_t i = 0; i< basisEdgeResult[0].nPatches(); i++)
            {
                // First Interface Side
                fileName = basename + "_0_" + util::to_string(i);
                gsField<T> temp_field(m_mp.patch(item.first().patch), basisEdgeResult[0].patch(i));
                gsWriteParaview(temp_field, fileName, 5000);
                collection.addTimestep(fileName, i, "0.vts");
                // Second Interface Side
                fileName = basename + "_1_" + util::to_string(i);
                gsField<T> temp_field_1(m_mp.patch(item.second().patch), basisEdgeResult[1].patch(i));
                gsWriteParaview(temp_field_1, fileName, 5000);
                collection.addTimestep(fileName, i, "0.vts");
            }
            collection.save();

            //gsWriteParaview(basisEdgeResult[0], "interface_basis", 20000);
        }
*/
    }

    gsApproxC1Edge(gsMultiPatch<T> const & mp,
                   BasisContainer & bases,
                const patchSide & item,
                size_t & numBdy,
                const gsOptionList & optionList)
                : m_mp(mp), m_bases(bases), m_optionList(optionList)
    {
        GISMO_UNUSED(numBdy);

        m_auxPatches.clear();
        m_auxPatches.push_back(gsPatchReparameterized<d,T>(m_mp.patch(item.patch), m_bases[item.patch]));

        std::vector<patchSide> sidesContainer(1);
        sidesContainer[0] = item;

        reparametrizeSinglePatch(item.side().index());

        compute(sidesContainer);
/*
        if (m_optionList.getSwitch("plot")) {
            std::string fileName;
            std::string basename = "BoundaryBasisFunctions" + util::to_string(numBdy);
            gsParaviewCollection collection(basename);

            for (size_t i = 0; i < basisEdgeResult[0].nPatches(); i++) {
                // First Interface Side
                fileName = basename + "_0_" + util::to_string(i);
                gsField<T> temp_field(m_mp.patch(item.patch), basisEdgeResult[0].patch(i));
                gsWriteParaview(temp_field, fileName, 5000);
                collection.newTimeStep(fileName, i, "0.vts");
            }
            collection.save();
        }
*/
    }

    std::vector<gsMultiPatch<T>> getEdgeBasis() { return basisEdgeResult; };

protected:

    // Input
    gsMultiPatch<T> const & m_mp;
    BasisContainer & m_bases;

    const gsOptionList & m_optionList;

    // Need for rotation, etc.
    C1AuxPatchContainer m_auxPatches;

    // Store temp solution
    std::vector<gsMultiPatch<T>> basisEdgeResult;

private:

    // Compute topology
    // After computeTopology() the patches will have the same patch-index as the position-index in auxGeom
    // EXAMPLE: global patch-index-order inside auxGeom: [2, 3, 4, 1, 0]
    //          in auxTop: 2->0, 3->1, 4->2, 1->3, 0->4
    void computeAuxTopology();

    void reparametrizeInterfacePatches(std::vector<patchSide> & sidesContainer);

    void reparametrizeSinglePatch(index_t side);

    void compute(std::vector<patchSide> & sidesContainer);

}; // Class gsApproxC1Edge

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsApproxC1Edge.hpp)
#endif
