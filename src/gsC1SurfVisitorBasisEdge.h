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

namespace gismo
{
    template <class T>
    class gsC1SurfVisitorBasisEdge
    {
    public:

        gsC1SurfVisitorBasisEdge()
        {
        }

        void initialize(const gsBasis<T>       & basis, //
                        gsQuadRule<T>    & rule)
        {
            gsVector<index_t> numQuadNodes( basis.dim() );
            for (index_t i = 0; i < basis.dim(); ++i) // to do: improve
                numQuadNodes[i] = basis.degree(i) + 1;

            // Setup Quadrature
            rule = gsGaussRule<T>(numQuadNodes);// NB!

            // Set Geometry evaluation flags
            md.flags = NEED_MEASURE ;
        }

        // Evaluate on element.
        // Batched over all "plus" (rows 0..plus_hi-plus_lo-1) and "minus" basis
        // functions of this patch side (the remaining rows), in that fixed
        // order: only rhsVals depends on the basis function, everything else
        // computed here (alpha/beta/N_0/N_1, and localMat in assemble()) is the
        // same mass matrix / gluing data for every bf on this side.
        inline void evaluate(const index_t plus_lo, const index_t plus_hi,
                             const index_t minus_lo, const index_t minus_hi,
                             gsBasis<T>       & basis, //
                             gsBasis<T>       & basis_geo,
                             gsBasis<T>       & basis_plus,
                             gsBasis<T>       & basis_minus,
                             const gsGeometry<T>    & geo, // patch
                             gsMatrix<T>            & quNodes,
                             index_t & uv,
                             gsC1SurfGluingData<T>  & gluingData,
                             bool & isBoundary)
        {
            md.points = quNodes;

            // Compute the active basis functions
            // Assumes actives are the same for all quadrature points on the elements
            basis.active_into(md.points.col(0), actives);

            // Evaluate basis functions on element
            basis.eval_into(md.points, basisData);

            // Compute geometry related values
            geo.computeMap(md);

            numActive = actives.rows();

            // tau/p
            const gsBSplineBasis<T> & bsp_temp = dynamic_cast<const gsBSplineBasis<T> & >(basis_geo);

            T p = basis_geo.maxDegree();
            T tau_1 = bsp_temp.knots().at(p + 2);

            gsMatrix<T> alpha, beta,
                    N_0, N_1,
                    N_j_minus, N_i_plus,
                    der_N_i_plus, temp;

            gsMatrix<T> pts1d;

            if (uv == 1) // edge is in v-direction
            {
                alpha = gluingData.evalAlpha_R(md.points.bottomRows(1));
                beta = gluingData.evalBeta_R(md.points.bottomRows(1));

                basis_geo.evalSingle_into(0,md.points.topRows(1),N_0); // u
                basis_geo.evalSingle_into(1,md.points.topRows(1),N_1); // u

                pts1d = md.points.bottomRows(1);
            } // Patch 0
            else // uv == 0, edge is in u-direction
            {
                alpha = gluingData.evalAlpha_L(md.points.topRows(1));
                beta = gluingData.evalBeta_L(md.points.topRows(1));

                basis_geo.evalSingle_into(0,md.points.bottomRows(1),N_0); // v
                basis_geo.evalSingle_into(1,md.points.bottomRows(1),N_1); // v

                pts1d = md.points.topRows(1);
            } // Patch 1

            if (isBoundary) beta.setZero();   // only the "plus" rows read beta
            if (isBoundary) alpha.setOnes();  // only the "minus" rows read alpha

            const index_t n_p = plus_hi  - plus_lo;
            const index_t n_m = minus_hi - minus_lo;

            rhsVals.resize(n_p + n_m, md.points.cols());

            // Neither operand depends on bfID: both are the same for every "plus" row.
            temp = beta.cwiseProduct(N_1);
            const gsMatrix<T> N_0_plus_N_1 = N_0 + N_1;

            index_t row = 0;
            for (index_t bfID = plus_lo; bfID < plus_hi; ++bfID, ++row)
            {
                basis_plus.evalSingle_into(bfID, pts1d, N_i_plus);
                basis_plus.derivSingle_into(bfID, pts1d, der_N_i_plus);

                rhsVals.row(row) = (N_i_plus.cwiseProduct(N_0_plus_N_1) - temp.cwiseProduct(der_N_i_plus) * tau_1 / p).row(0);
            } // n_plus

            if (uv == 1)
            {
                for (index_t bfID = minus_lo; bfID < minus_hi; ++bfID, ++row)
                {
                    basis_minus.evalSingle_into(bfID, pts1d, N_j_minus);
                    rhsVals.row(row) = (alpha.cwiseProduct(N_j_minus.cwiseProduct(N_1))).row(0);
                } // n_minus
            }
            else // uv == 0
            {
                for (index_t bfID = minus_lo; bfID < minus_hi; ++bfID, ++row)
                {
                    basis_minus.evalSingle_into(bfID, pts1d, N_j_minus);
                    rhsVals.row(row) = (- alpha.cwiseProduct(N_j_minus.cwiseProduct(N_1))).row(0);
                } // n_bar
            }

            localMat.setZero(numActive, numActive);
            localRhs.setZero(numActive, rhsVals.rows()); // multiple right-hand sides
        } // evaluate

        inline void assemble(gsDomainIteratorWrapper<T>    & element,
                             const gsVector<T>      & quWeights)
        {
            GISMO_UNUSED(element);

            gsMatrix<T> & basisVals  = basisData;

            // ( u, v)
            localMat.noalias() =
                    basisData * quWeights.asDiagonal() *
                    md.measures.asDiagonal() * basisData.transpose();

            for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
            {
                // Multiply weight by the geometry measure
                const T weight = quWeights[k] * md.measure(k);
                localRhs.noalias() += weight * (basisVals.col(k) * rhsVals.col(k).transpose());
            }

        }

        inline void localToGlobal(const index_t patchIndex,
                                  const std::vector<gsMatrix<T> >    & eliminatedDofs,
                                  gsSparseSystem<T>      & system)
        {
            // Map patch-local DoFs to global DoFs
            system.mapColIndices(actives, patchIndex, actives);
            // Add contributions to the system matrix and right-hand side
            system.push(localMat, localRhs, actives, eliminatedDofs[0], 0, 0);
        }

    protected:
        gsMatrix<index_t> actives;
        gsMatrix<T> basisData;
        index_t numActive;

    protected:
        // Local values of the right hand side
        gsMatrix<T>  rhsVals;

    protected:
        // Local matrices
        gsMatrix<T> localMat;
        gsMatrix<T> localRhs;

        gsMapData<T> md;

    }; // class gsVisitorG1BasisEdge
} // namespace gismo
