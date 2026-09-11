/** @file gsApproxC1Edge.hpp

    @brief Creates the approx C1 space.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/

#pragma once

#include <gsUnstructuredSplines/src/gsApproxC1Edge.h>

#include <gsUnstructuredSplines/src/gsPatchReparameterized.h>

#include <gsUnstructuredSplines/src/gsApproxC1GluingData.h>

#include <gsAssembler/gsDofMapperCreator.h>

namespace gismo
{
    template<short_t d,class T>
    void gsApproxC1Edge<d,T>::compute(std::vector<patchSide> & sidesContainer) {

        // Compute GLuing data
        gsApproxC1GluingData<d, T> approxGluingData(m_auxPatches, m_optionList, sidesContainer);

        gsTensorBSplineBasis<d, T> basis;
        index_t dir_pm, patch;
        if (sidesContainer.size() == 2) {
            if (m_auxPatches[0].getBasisRotated().piece(0).component(1).numElements() >
                m_auxPatches[1].getBasisRotated().piece(0).component(0).numElements()) {
                basis = dynamic_cast<const gsTensorBSplineBasis<d, T> &>(m_auxPatches[1].getBasisRotated().piece(0));
                dir_pm = 0;
                patch = 1;
            } else {
                basis = dynamic_cast<const gsTensorBSplineBasis<d, T> &>(m_auxPatches[0].getBasisRotated().piece(0));
                dir_pm = 1;
                patch = 0;
            }
        }
        else
        {
            basis = dynamic_cast<const gsTensorBSplineBasis<d, T> &>(m_auxPatches[0].getBasisRotated().piece(0));
            dir_pm = 1;
            patch = 0;
        }

        //! [Problem setup]
        basisEdgeResult.clear();
        for (size_t patchID = 0; patchID < sidesContainer.size(); patchID++)
        {
            gsMultiPatch<T> result;

            index_t dir = patchID == 0 ? 1 : 0;

            gsBSplineBasis<T> basis_plus, basis_minus;

            gsMultiBasis<T> initSpace(m_auxPatches[patchID].getBasisRotated().piece(0));
            createPlusSpace(m_auxPatches[patch].getPatchRotated(), basis, dir_pm, basis_plus);
            createMinusSpace(m_auxPatches[patch].getPatchRotated(), basis, dir_pm, basis_minus);

            gsGeometry<T> &geo = m_auxPatches[patchID].getPatchRotated();

            gsBSpline<T> beta, alpha;
            bool bdy = true;
            if (sidesContainer.size() == 2)
            {
                bdy = false;
                beta = approxGluingData.betaS(dir);
                alpha = approxGluingData.alphaS(dir);
            }

            // [!The same setup for each bf!]
            typename gsSparseSolver<T>::SimplicialLDLT solver;
            gsExprAssembler<T> A(1, 1);

            // Elements used for numerical integration
            gsMultiBasis<T> edgeSpace(
                    m_auxPatches[patchID].getBasisRotated().piece(sidesContainer[patchID]));

            A.setIntegrationElements(edgeSpace);
            gsExprEvaluator<T> ev(A);

            // Set the discretization space
            auto u = A.getSpace(edgeSpace);

            // Create Mapper
            gsDofMapper map = createMapper(edgeSpace, 1, false);
            if (!m_optionList.getSwitch("interpolation"))
            {
                gsMatrix<index_t> act;
                for (index_t i = 2; i < edgeSpace[0].component(1 - dir).size();
                     i++) // only the first two u/v-columns are Dofs (0/1)
                {
                    act = edgeSpace[0].boundaryOffset(dir == 0 ? 3 : 1, i); // WEST
                    map.markBoundary(0, act); // Patch 0
                }
                map.finalize();

                u.setupMapper(map);

                gsMatrix<T> &fixedDofs = const_cast<expr::gsFeSpace<T> &>(u).fixedPart();
                fixedDofs.setZero(u.mapper().boundarySize(), 1);

                A.initSystem();
                A.assemble(u * u.tr()); // The Matrix is the same for each bf
                solver.compute(A.matrix());
                // [!The same setup for each bf!]
            }

            index_t n_plus = basis_plus.size();
            index_t n_minus = basis_minus.size();

            const bool interpolation = m_optionList.getSwitch("interpolation");
            const index_t sz = edgeSpace.basis(0).size();

            index_t bfID_init = 3;
            if (interpolation)
            {
                for (index_t bfID = bfID_init; bfID < n_plus - bfID_init; bfID++) // first 3 and last 3 bf are eliminated
                {
                    gsTraceBasis<T> traceBasis(geo, beta, basis_plus, initSpace.basis(0), bdy, bfID, dir);
                    gsMatrix<T> anchors = edgeSpace.basis(0).anchors();
                    gsMatrix<T> values = traceBasis.eval(anchors);
                    result.addPatch(edgeSpace.basis(0).interpolateAtAnchors(give(values)));
                }
            }
            else if (n_plus - 2*bfID_init > 0)
            {
                // All (n_plus - 2*bfID_init) trace-basis right-hand sides are batched into
                // ONE quadrature pass (gsTraceBasisBatch) and ONE multi-column solve
                // against the mass-matrix factorization computed above.
                const index_t M = n_plus - 2*bfID_init;
                gsTraceBasisBatch<T> traceBasisBatch(geo, beta, basis_plus, initSpace.basis(0), bdy,
                                                      bfID_init, n_plus - bfID_init, dir);

                A.initVector(M);
                auto aa = A.getCoeff(traceBasisBatch);
                A.assemble(u * aa.tr());

                gsMatrix<T> solMat = solver.solve(A.rhs()); // sz_free x M

                // Build coefficients per column directly from the mapper, instead of
                // gsFeSpace::getCoeffs' multi-column path (which only fills the
                // eliminated-DoF value into column 0 of a multi-column result).
                // Eliminated rows read from fixedDofs, exactly as constructSolution()
                // does elsewhere in this submodule (e.g. gsC1SurfBasisEdge.h).
                const gsDofMapper & mapper = u.mapper();
                const gsMatrix<T> & fixedDofsRHS = u.fixedPart();
                for (index_t col = 0; col < M; ++col)
                {
                    gsMatrix<T> coefs(sz, 1);
                    for (index_t i = 0; i < sz; ++i)
                        coefs(i, 0) = mapper.is_free(i, 0)
                                    ? solMat(mapper.index(i, 0), col)
                                    : fixedDofsRHS(mapper.bindex(i, 0), 0);
                    result.addPatch(edgeSpace.basis(0).makeGeometry(give(coefs)));
                }
            }

            bfID_init = 2;
            if (interpolation)
            {
                for (index_t bfID = bfID_init; bfID < n_minus - bfID_init; bfID++)  // first 2 and last 2 bf are eliminated
                {
                    gsNormalDerivBasis<T> normalDerivBasis(geo, alpha, basis_minus, initSpace.basis(0), bdy, bfID,
                                                                dir);
                    gsMatrix<T> anchors = edgeSpace.basis(0).anchors();
                    gsMatrix<T> values = normalDerivBasis.eval(anchors);
                    result.addPatch(edgeSpace.basis(0).interpolateAtAnchors(give(values)));
                }
            }
            else if (n_minus - 2*bfID_init > 0)
            {
                const index_t M = n_minus - 2*bfID_init;
                gsNormalDerivBasisBatch<T> normalDerivBasisBatch(geo, alpha, basis_minus, initSpace.basis(0), bdy,
                                                                  bfID_init, n_minus - bfID_init, dir);

                A.initVector(M);
                auto aa = A.getCoeff(normalDerivBasisBatch);
                A.assemble(u * aa.tr());

                gsMatrix<T> solMat = solver.solve(A.rhs());

                const gsDofMapper & mapper = u.mapper();
                const gsMatrix<T> & fixedDofsRHS = u.fixedPart();
                for (index_t col = 0; col < M; ++col)
                {
                    gsMatrix<T> coefs(sz, 1);
                    for (index_t i = 0; i < sz; ++i)
                        coefs(i, 0) = mapper.is_free(i, 0)
                                    ? solMat(mapper.index(i, 0), col)
                                    : fixedDofsRHS(mapper.bindex(i, 0), 0);
                    result.addPatch(edgeSpace.basis(0).makeGeometry(give(coefs)));
                }
            }

            // parametrizeBasisBack
            m_auxPatches[patchID].parametrizeBasisBack(result);

            basisEdgeResult.push_back(result);
        }
    }


    template<short_t d,class T>
    void gsApproxC1Edge<d,T>::computeAuxTopology()
    {
        for(size_t i = 0; i <  m_auxPatches.size(); i++)
        {
            if(m_auxPatches[i].getPatchRotated().orientation() == -1)
                m_auxPatches[i].swapAxis();
        }
    }


    template<short_t d,class T>
    void gsApproxC1Edge<d,T>::reparametrizeInterfacePatches(std::vector<patchSide> & sidesContainer)
    {
        computeAuxTopology();

//        gsMultiPatch<T> temp_mp;
//        for(index_t i = 0; i <  m_auxPatches.size(); i++)
//            temp_mp.addPatch(m_auxPatches[i].getPatchRotated());
//
//        temp_mp.computeTopology();

        // Right patch along the interface. Patch 0 -> v coordinate. Edge west along interface
        switch (sidesContainer[0].side().index())
        {
            case 1:
                break;
            case 4: m_auxPatches[0].rotateParamClock();
                break;
            case 3: m_auxPatches[0].rotateParamAntiClock();
                break;
            case 2: m_auxPatches[0].rotateParamAntiClockTwice();
                break;
            default:
                break;
        }

        // Left patch along the interface. Patch 1 -> u coordinate. Edge south along interface
        switch (sidesContainer[1].side().index())
        {
            case 3:
                break;
            case 4: m_auxPatches[1].rotateParamAntiClockTwice();
                break;
            case 2: m_auxPatches[1].rotateParamAntiClock();
                break;
            case 1: m_auxPatches[1].rotateParamClock();
                break;
            default:
                break;
        }
    } // reparametrizeInterfacePatches


    template<short_t d,class T>
    void gsApproxC1Edge<d,T>::reparametrizeSinglePatch(index_t side)
    {
        computeAuxTopology();

        if(m_auxPatches[0].getOrient())
        {
            switch (side)
            {
                case 3:
                    break;
                case 2:
                    m_auxPatches[0].rotateParamClock();
                    break;
                case 4:
                    m_auxPatches[0].rotateParamAntiClockTwice();
                    break;
                case 1:
                    m_auxPatches[0].rotateParamAntiClock();
                    break;
            }
        }
        else
        {
            switch (side)
            {
                case 1:
                    break;
                case 4:
                    m_auxPatches[0].rotateParamClock();
                    break;
                case 2:
                    m_auxPatches[0].rotateParamAntiClockTwice();
                    break;
                case 3:
                    m_auxPatches[0].rotateParamAntiClock();
                    break;
            }
        }
    } // reparametrizeSinglePatch

} // namespace gismo
