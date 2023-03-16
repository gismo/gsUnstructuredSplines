/** @file gsApproxC1Vertex.h

    @brief Creates the (approx.) C1 Vertex space.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller & A. Farahat
*/

#pragma once

#include <gsUnstructuredSplines/src/gsApproxC1Utils.h>

#include <gsUnstructuredSplines/src/gsContainerBasis.h>
#include <gsUnstructuredSplines/src/gsPatchReparameterized.h>
#include <gsUnstructuredSplines/src/gsApproxC1GluingData.h>


namespace gismo {
template<short_t d, class T>
class gsApproxC1Vertex
{

private:
    typedef gsContainerBasis<d, T> Basis;
    typedef typename std::vector<Basis> BasisContainer;
    typedef typename std::vector<gsPatchReparameterized<d,T>> C1AuxPatchContainer;

    /// Shared pointer for gsApproxC1Vertex
    typedef memory::shared_ptr<gsApproxC1Vertex> Ptr;

    /// Unique pointer for gsApproxC1Vertex
    typedef memory::unique_ptr<gsApproxC1Vertex> uPtr;


public:
    /// Empty constructor
    ~gsApproxC1Vertex() { }


    gsApproxC1Vertex(gsMultiPatch<T> & mp,
                BasisContainer & bases,
                const std::vector<size_t> & patchesAroundVertex,
                const std::vector<size_t> & vertexIndices,
                const index_t & numVer,
                const gsOptionList & optionList)
                : m_mp(mp), m_bases(bases), m_patchesAroundVertex(patchesAroundVertex),
                m_vertexIndices(vertexIndices), m_optionList(optionList)
    {
        m_auxPatches.clear();
        basisVertexResult.clear();

        for(size_t i = 0; i < m_patchesAroundVertex.size(); i++)
        {
            index_t patch_1 = m_patchesAroundVertex[i];

            m_auxPatches.push_back(gsPatchReparameterized<d,T>(m_mp.patch(patch_1), m_bases[patch_1]));
        }

        reparametrizeVertexPatches();

        // Compute Sigma
        T sigma = computeSigma(m_vertexIndices);
        gsMatrix<T> Phi(6, 6);
        Phi.setIdentity();

        Phi.col(1) *= sigma;
        Phi.col(2) *= sigma;
        Phi.col(3) *= sigma * sigma;
        Phi.col(4) *= sigma * sigma;
        Phi.col(5) *= sigma * sigma;

        gsMultiPatch<T> rotPatch;
        if (m_auxPatches[0].getPatchRotated().parDim() + 1 == m_auxPatches[0].getPatchRotated().targetDim()) // Surface
        {
            gsMatrix<T> zero;
            zero.setZero(2, 1);
            gsMatrix<T> Jk = m_auxPatches[0].getPatchRotated().jacobian(zero);
            gsMatrix<T> G = Jk.transpose() * Jk; // Symmetric
            gsMatrix<T> G_inv = G.cramerInverse(); // Symmetric

            gsMatrix<T> geoMapDeriv1 = m_auxPatches[0].getPatchRotated()
                    .deriv(zero); // First derivative of the geometric mapping with respect to the parameter coordinates
            gsMatrix<T> geoMapDeriv2 = m_auxPatches[0].getPatchRotated()
                    .deriv2(zero); // Second derivative of the geometric mapping with respect to the parameter coordinates

            //Computing the normal vector to the tangent plane along the boundary curve
            gsVector<T> n(3);
            n.setZero();
            n(0) = Jk(1,0)*Jk(2,1)-Jk(2,0)*Jk(1,1);
            n(1) = Jk(2,0)*Jk(0,1)-Jk(0,0)*Jk(2,1);
            n(2) = Jk(0,0)*Jk(1,1)-Jk(1,0)*Jk(0,1);

            gsVector<T> z(3);
            z.setZero();
            z(2) = 1.0;

            gsVector<T> rotVec(3);
            rotVec.setZero(3);
            rotVec(0) = n(1,0)*z(2,0)-n(2,0)*z(1,0);
            rotVec(1) = n(2,0)*z(0,0)-n(0,0)*z(2,0);
            rotVec(2) = n(0,0)*z(1,0)-n(1,0)*z(0,0);

            T cos_t = (n.dot(z))/ (n.norm() * z.norm());
            T sin_t = rotVec.norm() / (n.norm() * z.norm());

//                Rotation matrix
            gsMatrix<T> R(3, 3);
            R.setZero();
//                Row 0
            R(0, 0) = cos_t + rotVec.x() * rotVec.x() * (1 - cos_t);
            R(0, 1) = rotVec.x() * rotVec.y() * (1 - cos_t) - rotVec.z() * sin_t;
            R(0, 2) = rotVec.x() * rotVec.z() * (1 - cos_t) + rotVec.y() * sin_t;
//                Row 1
            R(1, 0) = rotVec.x() * rotVec.y() * (1 - cos_t) + rotVec.z() * sin_t;
            R(1, 1) = cos_t + rotVec.y() * rotVec.y() * (1 - cos_t);
            R(1, 2) = rotVec.y() * rotVec.z() * (1 - cos_t) - rotVec.x() * sin_t;
//                Row 2
            R(2, 0) = rotVec.x() * rotVec.z() * (1 - cos_t) - rotVec.y() * sin_t;
            R(2, 1) = rotVec.y() * rotVec.z() * (1 - cos_t) + rotVec.x() * sin_t;
            R(2, 2) = cos_t + rotVec.z() * rotVec.z() * (1 - cos_t);

            for (size_t np = 0; np < m_auxPatches.size(); np++)
            {
                gsMatrix<T> coeffPatch = m_auxPatches[np].getPatchRotated().coefs();

                for (index_t i = 0; i < coeffPatch.rows(); i++)
                {
                    coeffPatch.row(i) =
                            (coeffPatch.row(i) - coeffPatch.row(0)) * R.transpose() + coeffPatch.row(0);
                }

                rotPatch.addPatch(m_auxPatches[np].getPatchRotated());
                rotPatch.patch(np).setCoefs(coeffPatch);
            }

            Phi.resize(13, 6);
            Phi.setZero();
            Phi(0, 0) = 1;
            Phi(1, 1) = sigma;
            Phi(2, 2) = sigma;
            Phi(4, 3) = sigma * sigma;
            Phi(5, 4) = sigma * sigma;
            Phi(7, 4) = sigma * sigma;
            Phi(8, 5) = sigma * sigma;
        }
        else
        {
            for (size_t np = 0; np < m_auxPatches.size(); np++)
                rotPatch.addPatch(m_auxPatches[np].getPatchRotated());
        }

        // TODO fix the gluing Data stuff
        for(size_t i = 0; i < m_patchesAroundVertex.size(); i++)
        {
            C1AuxPatchContainer auxPatchSingle;
            auxPatchSingle.push_back(m_auxPatches[i]);

            std::vector<patchSide> containingSides;  // global sides plus index
            patchCorner pC(m_patchesAroundVertex[i], m_vertexIndices[i]);
            pC.getContainingSides(d, containingSides);

//            std::vector<bool> isInterface(2);
//            isInterface[0] = m_mp.isInterface(patchSide(m_patchesAroundVertex[i], containingSides.at(0).side()));
//            isInterface[1] = m_mp.isInterface(patchSide(m_patchesAroundVertex[i], containingSides.at(1).side()));

            if (containingSides.at(0).side() < 3) // If isInterface_1 == v, then switch
            {
                patchSide side_temp = containingSides[0];
                containingSides[0] = containingSides[1];
                containingSides[1] = side_temp;

//                bool isInterface_tmp = isInterface[0];
//                isInterface[0] = isInterface[1];
//                isInterface[1] = isInterface_tmp;
            }

            std::vector<bool> isInterface(2);
            isInterface[0] = m_mp.isInterface(patchSide(m_patchesAroundVertex[i], containingSides.at(0).side()));  // global interface at u
            isInterface[1] = m_mp.isInterface(patchSide(m_patchesAroundVertex[i], containingSides.at(1).side()));  // global interface at v

            //Problem setup
            std::vector<gsBSpline<T>> alpha, beta;
            std::vector<gsBSplineBasis<T>> basis_plus, basis_minus;

            std::vector<bool> kindOfEdge;

            alpha.resize(2); beta.resize(2); basis_plus.resize(2); basis_minus.resize(2);
            kindOfEdge.resize(2);

            gsTensorBSplineBasis<d, T> basis = dynamic_cast<const gsTensorBSplineBasis<d, T> &>(auxPatchSingle[0].getBasisRotated().piece(0));
            gsTensorBSplineBasis<d, T> basis_pm = dynamic_cast<const gsTensorBSplineBasis<d, T> &>(auxPatchSingle[0].getBasisRotated().piece(0));
            for (size_t dir = 0; dir < containingSides.size(); ++dir)
            {
                //index_t localdir = auxPatchSingle[0].getMapIndex(containingSides[dir].index()) < 3 ? 1 : 0;
                if (isInterface[dir])
                {
                    patchSide result;
                    m_mp.getNeighbour(containingSides[dir], result);

                    index_t patch2 = -1;
                    for (size_t i_tmp = 0; i_tmp < m_patchesAroundVertex.size(); i_tmp++)
                        if (m_patchesAroundVertex[i_tmp] == size_t(result.patch))
                            patch2 = i_tmp;

                    GISMO_ASSERT(patch2 > -1, "Something went wrong");

                    gsTensorBSplineBasis<d, T> basis2 = dynamic_cast<const gsTensorBSplineBasis<d, T> &>(m_auxPatches[patch2].getBasisRotated().piece(
                            0));
                    index_t dir_1 = auxPatchSingle[0].getMapIndex(containingSides[dir].side()) < 3 ? 1 : 0;
                    index_t dir_2 = m_auxPatches[patch2].getMapIndex(result.side().index()) < 3 ? 1 : 0;
                    if (basis.component(dir_1).numElements() > basis2.component(dir_2).numElements())
                        basis_pm.component(dir_1) = basis2.component(dir_2);

                }
            }

            // Compute Gluing data
            // Stored locally
            gsApproxC1GluingData<d, T> approxGluingData(auxPatchSingle, m_optionList, containingSides, isInterface, basis_pm);

            //gsGeometry<T> & geo = auxPatchSingle[0].getPatchRotated();
            gsGeometry<T> & geo = rotPatch.patch(i);
            gsMultiBasis<T> initSpace(auxPatchSingle[0].getBasisRotated().piece(0));
            for (size_t dir = 0; dir < containingSides.size(); ++dir)
            {
                index_t localdir = auxPatchSingle[0].getMapIndex(containingSides[dir].index()) < 3 ? 1 : 0;

                if (isInterface[dir])
                {
                    alpha[localdir] = approxGluingData.alphaS(localdir);
                    beta[localdir] = approxGluingData.betaS(localdir);
                }
                // Store the isInterface locally
                kindOfEdge[localdir] = isInterface[dir];

                gsBSplineBasis<T> b_plus, b_minus;
                if (isInterface[dir])
                {
                    // TODO FIX
                    createPlusSpace(geo, basis_pm, localdir, b_plus);
                    createMinusSpace(geo, basis_pm, localdir, b_minus);
                    //createPlusSpace(geo, initSpace.basis(0), dir, b_plus);
                    //createMinusSpace(geo, initSpace.basis(0), dir, b_minus);
                }
                else
                {
                    createPlusSpace(geo, initSpace.basis(0), localdir, b_plus);
                    createMinusSpace(geo, initSpace.basis(0), localdir, b_minus);
                }

//                gsDebugVar(b_plus);
//                gsDebugVar(b_minus);
//                gsDebugVar(dir);
//                gsDebugVar(isInterface[localdir]);
                basis_plus[localdir] = b_plus;
                basis_minus[localdir] = b_minus;
            }

            typename gsSparseSolver<T>::SimplicialLDLT solver;
            gsExprAssembler<T> A(1, 1);

            // Elements used for numerical integration
            gsMultiBasis<T> vertexSpace(auxPatchSingle[0].getBasisRotated().piece(m_vertexIndices[i] + 4));
            A.setIntegrationElements(vertexSpace);
            gsExprEvaluator<T> ev(A);

            // Set the discretization space
            auto u = A.getSpace(vertexSpace);

            // Create Mapper
            gsDofMapper map(vertexSpace);
            if (!m_optionList.getSwitch("interpolation"))
            {
                gsMatrix<index_t> act;
                for (index_t dir = 0; dir < vertexSpace.basis(0).domainDim(); dir++)
                    for (index_t i = 3 * vertexSpace.basis(0).degree(dir) + 1; i < vertexSpace.basis(0).component(
                            1 - dir).size(); i++) // only the first two u/v-columns are Dofs (0/1)
                    {
                        act = vertexSpace.basis(0).boundaryOffset(dir == 0 ? 3 : 1, i); // WEST
                        //map.markBoundary(0, act); // Patch 0
                    }
                map.finalize();

                u.setupMapper(map);

                gsMatrix<T> &fixedDofs = const_cast<expr::gsFeSpace<T> &>(u).fixedPart();
                fixedDofs.setZero(u.mapper().boundarySize(), 1);

                A.initSystem();
                A.assemble(u * u.tr());
                solver.compute(A.matrix());
            }

            // Create Basis functions
            gsMultiPatch<T> result_1;
            for (index_t bfID = 0; bfID < 6; bfID++)
            {
                gsVertexBasis<T> vertexBasis(geo, basis_pm, alpha, beta, basis_plus, basis_minus, Phi,
                                             kindOfEdge, bfID);
                if (m_optionList.getSwitch("interpolation"))
                {
                    //gsQuasiInterpolate<T>::Schoenberg(edgeSpace.basis(0), traceBasis, sol);
                    //result.addPatch(edgeSpace.basis(0).interpolateAtAnchors(give(values)));
                    gsMatrix<T> anchors = vertexSpace.basis(0).anchors();
                    gsMatrix<T> values = vertexBasis.eval(anchors);
                    result_1.addPatch(vertexSpace.basis(0).interpolateAtAnchors(give(values)));
                }
                else
                {
                    A.initVector();

                    auto aa = A.getCoeff(vertexBasis);
                    A.assemble(u * aa);

                    gsMatrix<T> solVector = solver.solve(A.rhs());

                    auto u_sol = A.getSolution(u, solVector);
                    gsMatrix<T> sol;
                    u_sol.extract(sol);

                    result_1.addPatch(vertexSpace.basis(0).makeGeometry(give(sol)));
                }
            }
            //Problem setup end

            // Store temporary
            basisVertexResult.push_back(result_1);
        }

        gsMultiPatch<T> temp_mp;
        for (size_t j = 0; j < m_patchesAroundVertex.size(); j++)
            temp_mp.addPatch(m_mp.patch(m_patchesAroundVertex[j]));
        temp_mp.computeTopology();

        if (m_patchesAroundVertex.size() != temp_mp.interfaces().size()) // No internal vertex
        {
            computeKernel();

            for(size_t i = 0; i < m_patchesAroundVertex.size(); i++)
                m_auxPatches[i].parametrizeBasisBack(basisVertexResult[i]); // parametrizeBasisBack
        }
        else // Internal vertex
        {
            for(size_t i = 0; i < m_patchesAroundVertex.size(); i++)
                m_auxPatches[i].parametrizeBasisBack(basisVertexResult[i]); // parametrizeBasisBack
        }
/*
        if (m_optionList.getSwitch("plot"))
        {
            std::string fileName;
            std::string basename = "VerticesBasisFunctions" + util::to_string(numVer);
            gsParaviewCollection collection(basename);

            for (size_t np = 0; np < m_patchesAroundVertex.size(); ++np)
            {
                if (basisVertexResult.size() != 0)
                    for (size_t i = 0; i < basisVertexResult[np].nPatches(); ++i)
                    {
                        fileName = basename + "_" + util::to_string(np) + "_" + util::to_string(i);
                        gsField<T> temp_field(m_mp.patch(m_patchesAroundVertex[np]), basisVertexResult[np].patch(i));
                        gsWriteParaview(temp_field, fileName, 5000);
                        collection.addTimestep(fileName, i, "0.vts");

                    }
            }
            collection.save();
            //if (m_patchesAroundVertex.size() == 2)
            //    gsWriteParaview(basisVertexResult[0], "vertex_basis", 20000);
        }
*/
    }

    void reparametrizeVertexPatches();

    void checkOrientation(size_t i);

    T computeSigma(const std::vector<size_t> & vertexIndices);

    void computeKernel();

    std::vector<gsMultiPatch<T>> getVertexBasis() { return basisVertexResult; }


protected:

    // Input
    gsMultiPatch<T> & m_mp;
    BasisContainer & m_bases;

    const std::vector<size_t> & m_patchesAroundVertex;
    const std::vector<size_t> & m_vertexIndices;

    const gsOptionList & m_optionList;

    // Need for rotation, etc.
    C1AuxPatchContainer m_auxPatches;

    // Store temp solution
    std::vector<gsMultiPatch<T>> basisVertexResult;

}; // gsApproxC1Vertex

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsApproxC1Vertex.hpp)
#endif
