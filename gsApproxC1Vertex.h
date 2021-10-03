/** @file gsApproxC1Vertex.h

    @brief Creates the (approx.) C1 Vertex space.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller & A. Farahat
*/

#pragma once

#include <gsUnstructuredSplines/gsC1AuxiliaryPatch.h>
#include <gsUnstructuredSplines/gsApproxC1VertexBasisProjection.h>



namespace gismo {
template<short_t d, class T>
class gsApproxC1Vertex
{

private:
    typedef gsC1Basis<d, T> Basis;
    typedef typename std::vector<Basis> C1BasisContainer;
    typedef typename std::vector<gsC1AuxiliaryPatch<d,T>> C1AuxPatchContainer;

    /// Shared pointer for gsApproxC1Vertex
    typedef memory::shared_ptr<gsApproxC1Vertex> Ptr;

    /// Unique pointer for gsApproxC1Vertex
    typedef memory::unique_ptr<gsApproxC1Vertex> uPtr;


public:
    /// Empty constructor
    ~gsApproxC1Vertex() { }


    gsApproxC1Vertex(gsMultiPatch<T> const & mp,
                C1BasisContainer & bases,
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
            index_t vertex_1 = m_vertexIndices[i];
            index_t patch_1 = m_patchesAroundVertex[i];

            m_auxPatches.push_back(gsC1AuxiliaryPatch<d,T>(m_mp.patch(patch_1), m_bases[patch_1], vertex_1));
        }

        reparametrizeVertexPatches();

        // Compute Sigma
        real_t sigma = computeSigma(m_vertexIndices);

        std::vector<gsApproxGluingData<d, T>> gD; // delete later

        for(size_t i = 0; i < m_patchesAroundVertex.size(); i++)
        {
            C1AuxPatchContainer auxPatchSingle;
            auxPatchSingle.push_back(m_auxPatches[i]);

            std::vector<index_t> sideContainer;
            if (auxPatchSingle[0].getOrient() == 0) // not rotated
            {
                switch (auxPatchSingle[0].side()) // corner
                {
                    case 1:
                        sideContainer.push_back(3); // u
                        sideContainer.push_back(1); // v
                        break;
                    case 2:
                        sideContainer.push_back(3); // u
                        sideContainer.push_back(2); // v
                        break;
                    case 3:
                        sideContainer.push_back(4); // u
                        sideContainer.push_back(1); // v
                        break;
                    case 4:
                        sideContainer.push_back(4); // u
                        sideContainer.push_back(2); // v
                        break;
                    default:
                        gsInfo << "Something went wrong\n";
                        break;
                }
            }
            else if (auxPatchSingle[0].getOrient() == 1) // rotated
            {
                switch (auxPatchSingle[0].side()) // corner
                {
                    case 1:
                        sideContainer.push_back(1); // u
                        sideContainer.push_back(3); // v
                        break;
                    case 2:
                        sideContainer.push_back(3); // u
                        sideContainer.push_back(2); // v
                        break;
                    case 3:
                        sideContainer.push_back(4); // u
                        sideContainer.push_back(1); // v
                        break;
                    case 4:
                        sideContainer.push_back(2); // u
                        sideContainer.push_back(4); // v
                        break;
                    default:
                        gsInfo << "Something went wrong\n";
                        break;
                }
            }

            //gsInfo << "Patch: " << m_patchesAroundVertex[i] << " with index: " << m_vertexIndices[i] << "\n";

            // Compute Gluing data
            gsApproxGluingData<d, T> approxGluingData(auxPatchSingle, m_optionList, sideContainer);

            // Create Basis functions
            gsMultiPatch<> result_1;
            gsApproxC1VertexBasisProjection<d, T> approxC1VertexBasis(auxPatchSingle, approxGluingData,
                                                                                m_vertexIndices[i], sideContainer, sigma, m_optionList);
            approxC1VertexBasis.setBasisVertex(result_1);


            // Store temporary
            basisVertexResult.push_back(result_1);

            gD.push_back(approxGluingData); // delete later
        }

        if (m_auxPatches[0].getC1BasisRotated().getKindOfVertex(m_vertexIndices[0]) != 0) // No internal vertex
        {
/*
            if (m_patchesAroundVertex.size() == 2 && m_vertexIndices[0] == 2)
            {
                for (size_t np = 0; np<basisVertexResult[0].nPatches(); np++) {
                    gsMatrix<> points_u, points_v;
                    points_u.setZero(2, 10);
                    points_v.setZero(2, 10);
                    gsVector<> lin;
                    lin.setLinSpaced(10, 0, 1);
                    points_u.row(0) = lin.transpose();
                    points_v.row(1) = lin.transpose();
*/
                    /*
                    gsInfo << "points: " << points_u << "\n";

                    gsInfo << "alphaS: " << gD[0].alphaS(0).eval(points_v.row(1)) << "\n";
                    gsInfo << "alphaS: " << gD[1].alphaS(1).eval(points_v.row(1)) << "\n";
                    gsInfo << "betaS: " << gD[0].betaS(0).eval(points_v.row(1)) << "\n";
                    gsInfo << "betaS: " << gD[1].betaS(1).eval(points_v.row(1)) << "\n";

                    gsInfo << "GLUINGDATA: " << gD[0].alphaS(0).eval(points_v.row(1)).cwiseProduct(
                            gD[1].betaS(1).eval(points_v.row(1))) +
                                                gD[1].alphaS(1).eval(points_v.row(1)).cwiseProduct(
                                                        gD[0].betaS(0).eval(points_v.row(1))) << "\n";

                    gsInfo << "DERIV: " << basisVertexResult[0].patch(np).deriv(points_u) << "\n";
                    gsInfo << "DERIV 2: " << basisVertexResult[1].patch(np).deriv(points_v) << "\n";
                    gsInfo << "part1: " << gD[0].alphaS(0).eval(points_v.row(1)).cwiseProduct(
                            basisVertexResult[1].patch(np).deriv(points_v).row(0)) << "\n";
                    gsInfo << "part2: " << gD[1].alphaS(1).eval(points_v.row(1)).cwiseProduct(
                            basisVertexResult[0].patch(np).deriv(points_u).row(1)) << "\n";
                    */
/*
                    gsInfo << "np: " << np << "\n";


                    gsInfo << "DERIV: " << basisVertexResult[0].patch(np).deriv(points_u) << "\n";
                    gsInfo << "DERIV 2: " << basisVertexResult[1].patch(np).deriv(points_v) << "\n";


                    gsInfo << "alphaS: " << gD[0].alphaS(0).eval(points_v.row(1)) << "\n";
                    gsInfo << "alphaS: " << gD[1].alphaS(1).eval(points_v.row(1)) << "\n";

                    gsInfo << "GLUINGDATA: " << gD[0].alphaS(0).eval(points_v.row(1)).cwiseProduct(
                            gD[1].betaS(1).eval(points_v.row(1))) +
                                                gD[1].alphaS(1).eval(points_v.row(1)).cwiseProduct(
                                                        gD[0].betaS(0).eval(points_v.row(1))) << "\n";

                    gsInfo << "test: " << gD[0].alphaS(0).eval(points_v.row(1)).cwiseProduct(
                            basisVertexResult[1].patch(np).deriv(points_v).row(0))
                                          + gD[1].alphaS(1).eval(points_v.row(1)).cwiseProduct(
                            basisVertexResult[0].patch(np).deriv(points_u).row(1)) << "\n";
                }
           }
*/
            computeKernel();
/*
            gsMatrix<> points;
            points.setZero(2,2);
            points(0,1) = 1.0;

            if (m_patchesAroundVertex.size() == 2 && !m_optionList.getSwitch("noVertex"))
                for (size_t np = 0; np < m_patchesAroundVertex.size(); ++np)
                    for (size_t i = 0; i < basisVertexResult[np].nPatches(); ++i)
                        gsInfo << i << " : " << basisVertexResult[np].patch(i).deriv(points.col(0)) << "\n\n";
*/



            for(size_t i = 0; i < m_patchesAroundVertex.size(); i++)
                m_auxPatches[i].parametrizeBasisBack(basisVertexResult[i]); // parametrizeBasisBack

/*
            if (m_patchesAroundVertex.size() == 2 && m_vertexIndices[0] == 2)
            {
                gsWriteParaview(basisVertexResult[0].patch(0),"LeftPatch",2000);
                gsWriteParaview(basisVertexResult[1].patch(0),"RightPatch",2000);

                std::string fileName;
                std::string basename = "VerticesBasisFunctions" + util::to_string(numVer);
                gsParaviewCollection collection(basename);

                for (size_t np = 0; np < m_patchesAroundVertex.size(); ++np)
                {
                    if (basisVertexResult.size() != 0)
                        for (size_t i = 0; i < basisVertexResult[np].nPatches(); ++i)
                        {
                            fileName = basename + "_" + util::to_string(np) + "_" + util::to_string(i);
                            gsField<> temp_field(m_auxPatches[m_patchesAroundVertex[np]].getPatch(), basisVertexResult[np].patch(i));
                            gsWriteParaview(temp_field, fileName, 5000);
                            collection.addTimestep(fileName, i, "0.vts");

                        }
                }
                collection.save();

            }
*/

/*
            gsInfo << "\n";

            if (m_patchesAroundVertex.size() == 2 && !m_optionList.getSwitch("noVertex"))
                for (size_t np = 0; np < m_patchesAroundVertex.size(); ++np)
                    for (size_t i = 0; i < basisVertexResult[np].nPatches(); ++i)
                        gsInfo << i << " : " << basisVertexResult[np].patch(i).deriv(points.col(1-np)) << "\n\n";
*/

        }
        else // Internal vertex
        {
            for(size_t i = 0; i < m_patchesAroundVertex.size(); i++)
                m_auxPatches[i].parametrizeBasisBack(basisVertexResult[i]); // parametrizeBasisBack
        }




        /*
        gsMatrix<> points;
        points.setZero(2,2);
        points(0,1) = 1.0;


        if (m_patchesAroundVertex.size() == 2 && !m_optionList.getSwitch("noVertex"))
            for (size_t np = 0; np < m_patchesAroundVertex.size(); ++np)
                for (size_t i = 0; i < basisVertexResult[np].nPatches(); ++i)
                    gsInfo << basisVertexResult[np].patch(i).deriv(points.col(1-np)) << "\n\n";
        */

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
                        gsField<> temp_field(m_mp.patch(m_patchesAroundVertex[np]), basisVertexResult[np].patch(i));
                        gsWriteParaview(temp_field, fileName, 5000);
                        collection.addTimestep(fileName, i, "0.vts");

                    }
            }
            collection.save();
        }
    }

    gsApproxC1Vertex(gsMultiPatch<T> const & mp,
                      C1BasisContainer & bases,
                      const std::vector<size_t> & patchesAroundVertex,
                      const std::vector<size_t> & vertexIndices,
                      const index_t & numVer,
                      std::vector<std::vector<gsMultiPatch<T>>> & vertex_bf,
                      const gsOptionList & optionList)
            : m_mp(mp), m_bases(bases), m_patchesAroundVertex(patchesAroundVertex),
              m_vertexIndices(vertexIndices), m_optionList(optionList)
    {
        m_auxPatches.clear();
        basisVertexResult.clear();

        for(size_t i = 0; i < m_patchesAroundVertex.size(); i++)
        {
            index_t vertex_1 = m_vertexIndices[i];
            index_t patch_1 = m_patchesAroundVertex[i];

            m_auxPatches.push_back(gsC1AuxiliaryPatch<d,T>(m_mp.patch(patch_1), m_bases[patch_1], vertex_1));
        }

        // 2 == u,v
        std::vector<std::vector<gsMultiPatch<T>>> vertexAuxiliary_bf(patchesAroundVertex.size(), std::vector<gsMultiPatch<T>>(2));
        reparametrizeVertexBasis(vertex_bf, vertexAuxiliary_bf);

        gsMultiPatch<T> mp_vertex;
        for(size_t i = 0; i < patchesAroundVertex.size(); i++)
            mp_vertex.addPatch(m_mp.patch(patchesAroundVertex[i]));
        mp_vertex.computeTopology();

        gsMatrix<T> coeffs_mat;
        coeffs_mat.setZero(mp_vertex.nPatches()*4, mp_vertex.nInterfaces()*5 ); // *8, + mp_vertex.nPatches()*4

        gsMatrix<> points(2, 1);
        points.setZero();

        for(size_t i = 0; i < patchesAroundVertex.size(); i++)
        {
            for(index_t dir = 0; dir < 2; dir++)
            {
                index_t i_shift = 0;
                for(size_t numInt = 0; numInt < mp_vertex.nInterfaces(); numInt++)
                {
                    if (mp_vertex.interfaces()[numInt].first().patch == (index_t) i)
                    {
                        if (mp_vertex.interfaces()[numInt].first().side().index() > 2 && dir == 0)
                            i_shift = numInt*5;
                        else if (mp_vertex.interfaces()[numInt].first().side().index() < 3 && dir == 1)
                            i_shift = numInt*5;
                    }
                    if (mp_vertex.interfaces()[numInt].second().patch == (index_t) i)
                    {
                        if (mp_vertex.interfaces()[numInt].second().side().index() > 2 && dir == 0)
                            i_shift = numInt*5;
                        else if (mp_vertex.interfaces()[numInt].second().side().index() < 3 && dir == 1)
                            i_shift = numInt*5;
                    }
                }

                for (size_t ii = 0; ii < vertexAuxiliary_bf[i][dir].nPatches(); ii++) {
                    gsMatrix<> basisVal = vertexAuxiliary_bf[i][dir].patch(ii).eval(points);
                    gsMatrix<> basisDer = vertexAuxiliary_bf[i][dir].patch(ii).deriv(points);
                    gsMatrix<> basisDer2 = vertexAuxiliary_bf[i][dir].patch(ii).deriv2(points);
                    if (basisVal(0, 0) * basisVal(0, 0) > 1e-25)
                        coeffs_mat(i * 4 + 0, i_shift + ii) = (dir == 1 ? -1 : 1) * basisVal(0, 0);
                    if (basisDer(0, 0) * basisDer(0, 0) > 1e-25)
                        coeffs_mat(i * 4 + 1, i_shift + ii) = (dir == 1 ? -1 : 1) * basisDer(0, 0); // u
                    if (basisDer(1, 0) * basisDer(1, 0) > 1e-25)
                        coeffs_mat(i * 4 + 2, i_shift + ii) = (dir == 1 ? -1 : 1) * basisDer(1, 0); // v
                    if (basisDer2(2, 0) * basisDer2(2, 0) > 1e-25)
                        coeffs_mat(i * 4 + 3, i_shift + ii) = (dir == 1 ? -1 : 1) * basisDer2(2, 0); // uv
                }
            }
        }
        // Correction term
/*
        std::vector<size_t> interface_indices;
        index_t patchID_next = 0;
        for(size_t i = 0; i < patchesAroundVertex.size(); i++)
        {
            index_t patchID = patchID_next;
            index_t i_shift = 0, dir = 0;
            for(size_t numInt = 0; numInt < mp_vertex.nInterfaces(); numInt++)
            {
                bool intexist = false;
                for (size_t i_indizes = 0; i_indizes < interface_indices.size(); i_indizes++)
                    if (interface_indices[i_indizes] == numInt)
                        intexist = true;

                if (mp_vertex.interfaces()[numInt].first().patch == patchID && !intexist)
                {
                    i_shift = numInt;
                    dir = mp_vertex.interfaces()[numInt].first().side().index() > 2 ? 0 : 1;
                    patchID_next = mp_vertex.interfaces()[numInt].second().patch;
                    interface_indices.push_back(numInt);
                    break;
                }
                if (mp_vertex.interfaces()[numInt].second().patch == patchID  && !intexist)
                {
                    i_shift = numInt;
                    dir = mp_vertex.interfaces()[numInt].second().side().index() > 2 ? 0 : 1;
                    patchID_next = mp_vertex.interfaces()[numInt].first().patch;
                    interface_indices.push_back(numInt);
                    break;
                }
            }

            for (size_t ii = 0; ii < vertexAuxiliary_bf[patchID][dir].nPatches(); ii++) {
                gsMatrix<> basisVal = vertexAuxiliary_bf[patchID][dir].patch(ii).eval(points);
                gsMatrix<> basisDer = vertexAuxiliary_bf[patchID][dir].patch(ii).deriv(points);
                gsMatrix<> basisDer2 = vertexAuxiliary_bf[patchID][dir].patch(ii).deriv2(points);
                if (basisVal(0, 0) * basisVal(0, 0) > 1e-25)
                    coeffs_mat(mp_vertex.nPatches() * 4 + i * 4 + 0, i_shift*5 + ii) =
                            (dir == 1 ? -1 : 1) * basisVal(0, 0);
                if (basisDer(0, 0) * basisDer(0, 0) > 1e-25)
                    coeffs_mat(mp_vertex.nPatches() * 4 + i * 4 + 1, i_shift*5 + ii) =
                            (dir == 1 ? -1 : 1) * basisDer(0, 0); // u
                if (basisDer(1, 0) * basisDer(1, 0) > 1e-25)
                    coeffs_mat(mp_vertex.nPatches() * 4 + i * 4 + 2, i_shift*5 + ii) =
                            (dir == 1 ? -1 : 1) * basisDer(1, 0); // v
                if (basisDer2(2, 0) * basisDer2(2, 0) > 1e-25)
                    coeffs_mat(mp_vertex.nPatches() * 4 + i * 4 + 3, i_shift*5 + ii) =
                            (dir == 1 ? -1 : 1) * basisDer2(2, 0); // uv
            }

            gsBasis<> &basis_temp = vertexAuxiliary_bf[patchID][dir].patch(0).basis();
            for (size_t ii = 0; ii < 4; ii++) {
                coeffs_mat(mp_vertex.nPatches() * 4 + i * 4 + ii, i_shift*4+ mp_vertex.nInterfaces()*5 + 0) = (dir == 1 ? 1 : -1) * basis_temp.evalSingle(0,points)(0,0);
                coeffs_mat(mp_vertex.nPatches() * 4 + i * 4 + ii, i_shift*4+ mp_vertex.nInterfaces()*5 + 1) = (dir == 1 ? 1 : -1) * basis_temp.derivSingle(1,points)(0,0);
                coeffs_mat(mp_vertex.nPatches() * 4 + i * 4 + ii, i_shift*4+ mp_vertex.nInterfaces()*5 + 2) = (dir == 1 ? 1 : -1) * basis_temp.derivSingle(basis_temp.component(0).size(),points)(1,0);
                coeffs_mat(mp_vertex.nPatches() * 4 + i * 4 + ii, i_shift*4+ mp_vertex.nInterfaces()*5 + 3) = (dir == 1 ? 1 : -1) * basis_temp.deriv2Single(basis_temp.component(0).size()+1,points)(2,0);
            }


        }
*/


        //gsInfo << "coeffs_mat " << coeffs_mat << "\n";
        real_t threshold = 1e-8;
        Eigen::FullPivLU<gsMatrix<>> KernelVertex(coeffs_mat);
        KernelVertex.setThreshold(threshold);
        //gsInfo << "Coefs: " << coefs_corner << "\n";
        while (KernelVertex.dimensionOfKernel() < (index_t) mp_vertex.nInterfaces()+3)
        {
            threshold += 1e-8;
            KernelVertex.setThreshold(threshold);
        }
        gsMatrix<T> kernel = KernelVertex.kernel();

        //gsInfo << "Kernel: " << kernel << "\n";
        //gsInfo << "Kernel dim: " << KernelVertex.dimensionOfKernel() << "\n";

        basisVertexResult.resize(patchesAroundVertex.size());
        for(size_t i = 0; i < patchesAroundVertex.size(); i++) {
            for (index_t kernel_dim = 0; kernel_dim < kernel.cols(); kernel_dim++) {
                gsMatrix<> coef_bf, coef_bf_singleU, coef_bf_singleV;
                coef_bf.setZero(vertexAuxiliary_bf[i][0].patch(0).coefs().rows(), 1);
                coef_bf_singleU.setZero(vertexAuxiliary_bf[i][0].patch(0).coefs().rows(), 1);
                coef_bf_singleV.setZero(vertexAuxiliary_bf[i][0].patch(0).coefs().rows(), 1);

                for (index_t dir = 0; dir < 2; dir++) {
                    index_t i_shift = 0;
                    for (size_t numInt = 0; numInt < mp_vertex.nInterfaces(); numInt++) {
                        if (mp_vertex.interfaces()[numInt].first().patch == (index_t) i) {
                            if (mp_vertex.interfaces()[numInt].first().side().index() > 2 && dir == 0)
                                i_shift = numInt * 5;
                            else if (mp_vertex.interfaces()[numInt].first().side().index() < 3 && dir == 1)
                                i_shift = numInt * 5;
                        }
                        if (mp_vertex.interfaces()[numInt].second().patch == (index_t) i) {
                            if (mp_vertex.interfaces()[numInt].second().side().index() > 2 && dir == 0)
                                i_shift = numInt * 5;
                            else if (mp_vertex.interfaces()[numInt].second().side().index() < 3 && dir == 1)
                                i_shift = numInt * 5;
                        }
                    }
                    for (size_t ii = 0; ii < vertexAuxiliary_bf[i][dir].nPatches(); ii++) {
                        coef_bf += vertexAuxiliary_bf[i][dir].patch(ii).coefs() * kernel(i_shift + ii, kernel_dim);
                        if (dir == 0)
                            coef_bf_singleU += vertexAuxiliary_bf[i][dir].patch(ii).coefs() * kernel(i_shift + ii, kernel_dim);
                        if (dir == 1)
                            coef_bf_singleV += vertexAuxiliary_bf[i][dir].patch(ii).coefs() * kernel(i_shift + ii, kernel_dim);
                    }
                }
                gsBasis<> &basis_temp = vertexAuxiliary_bf[i][0].patch(0).basis();
                gsGeometry<>::uPtr geo_tempU = basis_temp.makeGeometry(coef_bf_singleU);
                gsGeometry<>::uPtr geo_tempV = basis_temp.makeGeometry(coef_bf_singleV);
                gsMatrix<> coef_bf_corr;
                coef_bf_corr.setZero(vertexAuxiliary_bf[i][0].patch(0).coefs().rows(), 1);
                coef_bf_corr(0,0) = coef_bf_singleU(0,0);
                coef_bf_corr(1,0) = coef_bf_singleU(1,0);
                coef_bf_corr(basis_temp.component(0).size(),0) = coef_bf_singleU(basis_temp.component(0).size(),0);
                coef_bf_corr(basis_temp.component(0).size()+1,0) = coef_bf_singleU(basis_temp.component(0).size()+1,0);

                //gsInfo << "coef_bf_corr: " << coef_bf_corr << "\n";
                //gsInfo << "coef_bf: " << coef_bf << "\n";
                coef_bf -= coef_bf_corr;

/*
                gsGeometry<>::uPtr geo_temp = basis_temp.makeGeometry(coef_bf);
                gsInfo << "Patch: " << i << "\n";
                gsInfo << "eval: " << geo_temp->eval(points)(0,0) << "\n";
                gsInfo << "deriv: " << geo_temp->deriv(points)(0,0) << "\n";
                gsInfo << "deriv: " << geo_temp->deriv(points)(1,0) << "\n";
                gsInfo << "deriv2: " << geo_temp->deriv2(points)(2,0) << "\n";

                coef_bf_corr.setZero(vertexAuxiliary_bf[i][0].patch(0).coefs().rows(), 1);
                coef_bf_corr(0,0) = geo_tempU->eval(points)(0,0);
                coef_bf -= coef_bf_corr;

                gsInfo << "eval: " << geo_tempU->eval(points)(0,0) << " = " << geo_tempV->eval(points)(0,0) << "\n";
                gsInfo << "deriv: " << geo_tempU->deriv(points)(0,0) << " = " << geo_tempV->deriv(points)(0,0) << "\n";
                gsInfo << "deriv: " << geo_tempU->deriv(points)(1,0) << " = " << geo_tempV->deriv(points)(1,0) << "\n";
                gsInfo << "deriv2: " << geo_tempU->deriv2(points)(2,0) << " = " << geo_tempV->deriv2(points)(2,0) << "\n";

                coef_bf_corr.setZero();
                coef_bf_corr(1,0) = geo_tempU->deriv(points)(0,0) / basis_temp.derivSingle(1,points)(0,0);
                coef_bf -= coef_bf_corr;

                coef_bf_corr.setZero();
                coef_bf_corr(basis_temp.component(0).size(),0) = geo_tempU->deriv(points)(1,0) /
                        basis_temp.derivSingle(basis_temp.component(0).size(), points)(1,0);
                coef_bf -= coef_bf_corr;

                coef_bf_corr.setZero();
                coef_bf_corr(basis_temp.component(0).size()+1,0) = geo_tempU->deriv2(points)(2,0)/
                        basis_temp.deriv2Single(basis_temp.component(0).size()+1,points)(2,0);
                coef_bf -= coef_bf_corr;
*/
                basisVertexResult[i].addPatch(basis_temp.makeGeometry(coef_bf));
            }
        }

        for(size_t i = 0; i < m_patchesAroundVertex.size(); i++)
            m_auxPatches[i].parametrizeBasisBack(basisVertexResult[i]); // parametrizeBasisBack

        if (m_optionList.getSwitch("plot"))
        {
            std::string fileName;
            std::string basename = "NoVerticesBasisFunctions" + util::to_string(numVer);
            gsParaviewCollection collection(basename);

            for (size_t np = 0; np < m_patchesAroundVertex.size(); ++np)
            {
                if (basisVertexResult.size() != 0)
                    for (size_t i = 0; i < basisVertexResult[np].nPatches(); ++i)
                    {
                        fileName = basename + "_" + util::to_string(np) + "_" + util::to_string(i);
                        gsField<> temp_field(m_mp.patch(m_patchesAroundVertex[np]), basisVertexResult[np].patch(i));
                        gsWriteParaview(temp_field, fileName, 5000);
                        collection.addTimestep(fileName, i, "0.vts");

                    }
            }
            collection.save();
        }
    }

    void saveBasisVertex(gsSparseMatrix<T> & system)
    {
        for (size_t np = 0; np < m_patchesAroundVertex.size(); ++np)
        {
            index_t patch_1 = m_patchesAroundVertex[np];
            index_t corner = m_vertexIndices[np];

            index_t shift_row = 0, shift_col = 0;
            for (index_t np = 0; np < patch_1; ++np)
            {
                shift_row += m_bases[np].size_rows();
                shift_col += m_bases[np].size_cols();
            }

            index_t ii = 0;
            for (index_t i = m_bases[patch_1].rowBegin(corner+4); i < m_bases[patch_1].rowEnd(corner+4); ++i, ++ii)
            {
                index_t jj = 0;
                for (index_t j = m_bases[patch_1].colBegin(corner+4); j < m_bases[patch_1].colEnd(corner+4); ++j, ++jj)
                    if (basisVertexResult[np].patch(ii).coef(jj,0)*basisVertexResult[np].patch(ii).coef(jj,0)>1e-25)
                        system.insert(shift_row+i,shift_col+j) = basisVertexResult[np].patch(ii).coef(jj,0);
            }
        }
    }

    void reparametrizeVertexPatches()
    {
        for(size_t i = 0; i < m_auxPatches.size(); i++)
        {
            checkOrientation(i); // Check if the orientation is correct. If not, modifies vertex and edge vectors

            if(m_auxPatches[i].getOrient() == 0) // not switched
                switch (m_auxPatches[i].side()) // == vertex
                {
                    case 1:
                        //gsInfo << "Patch: " << m_patchesAroundVertex[i] << " with side " << m_vertexIndices[i]  << " not rotated\n";
                        break;
                    case 4:
                        m_auxPatches[i].rotateParamAntiClockTwice();
                        //gsInfo << "Patch: " << m_patchesAroundVertex[i] << " with side " << m_vertexIndices[i]  << " rotated twice anticlockwise\n";
                        break;
                    case 2:
                        m_auxPatches[i].rotateParamAntiClock();
                        //gsInfo << "Patch: " << m_patchesAroundVertex[i] << " with side " << m_vertexIndices[i]  << " rotated anticlockwise\n";
                        break;
                    case 3:
                        m_auxPatches[i].rotateParamClock();
                        //gsInfo << "Patch: " << m_patchesAroundVertex[i] << " with side " << m_vertexIndices[i]  << " rotated clockwise\n";
                        break;
                }
            else if (m_auxPatches[i].getOrient() == 1) // switched
                switch (m_auxPatches[i].side()) // == vertex
                {
                    case 1:
                        //gsInfo << "Patch: " << m_patchesAroundVertex[i] << " with side " << m_vertexIndices[i]  << " not rotated\n";
                        break;
                    case 4:
                        m_auxPatches[i].rotateParamAntiClockTwice();
                        //gsInfo << "Patch: " << m_patchesAroundVertex[i] << " with side " << m_vertexIndices[i]  << " rotated twice anticlockwise\n";
                        break;
                    case 3:
                        m_auxPatches[i].rotateParamAntiClock();
                        //gsInfo << "Patch: " << m_patchesAroundVertex[i] << " with side " << m_vertexIndices[i]  << " rotated anticlockwise\n";
                        break;
                    case 2:
                        m_auxPatches[i].rotateParamClock();
                        //gsInfo << "Patch: " << m_patchesAroundVertex[i] << " with side " << m_vertexIndices[i]  << " rotated clockwise\n";
                        break;
                }
        }
    }

    void reparametrizeVertexBasis(std::vector<std::vector<gsMultiPatch<T>>> & vertex_bf, std::vector<std::vector<gsMultiPatch<T>>> & result) // only for noVertex
    {
        for(size_t i = 0; i < m_auxPatches.size(); i++)
        {
            size_t patch = m_patchesAroundVertex[i];
            size_t vertex = m_vertexIndices[i];
            index_t dir_u = -1, dir_v = -1;
            switch (vertex) // == vertex
            {
                case 1:
                    dir_u = 4;
                    dir_v = 0;
                    break;
                case 2:
                    dir_u = 5;
                    dir_v = 2;
                    break;
                case 3:
                    dir_u = 6;
                    dir_v = 1;
                    break;
                case 4:
                    dir_u = 7;
                    dir_v = 3;
                    break;
                default:
                    break;
            }

            gsMultiPatch<> basisVertex_u = vertex_bf[patch][dir_u];
            gsMultiPatch<> basisVertex_v = vertex_bf[patch][dir_v];

            if (m_auxPatches[i].getPatch().orientation() == -1)
            {
                m_auxPatches[i].swapBasisAxis(basisVertex_u);
                m_auxPatches[i].swapBasisAxis(basisVertex_v);
                gsInfo << "Changed axis on patch: " << m_patchesAroundVertex[i] << " with side " << m_vertexIndices[i] << "\n";
            }

            if(m_auxPatches[i].getOrient() == 0) // not switched
                switch (m_auxPatches[i].side()) // == vertex
                {
                    case 1:
                        //gsInfo << "Patch: " << m_patchesAroundVertex[i] << " with side " << m_vertexIndices[i]  << " not rotated\n";
                        break;
                    case 4:
                        m_auxPatches[i].rotateBasisAntiClockTwice(basisVertex_u);
                        m_auxPatches[i].rotateBasisAntiClockTwice(basisVertex_v);
                        m_auxPatches[i].setNumberOfRotation(2);
                        //gsInfo << "Patch: " << m_patchesAroundVertex[i] << " with side " << m_vertexIndices[i]  << " rotated twice anticlockwise\n";
                        break;
                    case 2:
                        m_auxPatches[i].rotateBasisAntiClock(basisVertex_u);
                        m_auxPatches[i].rotateBasisAntiClock(basisVertex_v);
                        m_auxPatches[i].setNumberOfRotation(1);
                        //gsInfo << "Patch: " << m_patchesAroundVertex[i] << " with side " << m_vertexIndices[i]  << " rotated anticlockwise\n";
                        break;
                    case 3:
                        m_auxPatches[i].rotateBasisClock(basisVertex_u);
                        m_auxPatches[i].rotateBasisClock(basisVertex_v);
                        m_auxPatches[i].setNumberOfRotation(-1);
                        //gsInfo << "Patch: " << m_patchesAroundVertex[i] << " with side " << m_vertexIndices[i]  << " rotated clockwise\n";
                        break;
                }
            else if (m_auxPatches[i].getOrient() == 1) // switched
                switch (m_auxPatches[i].side()) // == vertex
                {
                    case 1:
                        //gsInfo << "Patch: " << m_patchesAroundVertex[i] << " with side " << m_vertexIndices[i]  << " not rotated\n";
                        break;
                    case 4:
                        m_auxPatches[i].rotateBasisAntiClockTwice(basisVertex_u);
                        m_auxPatches[i].rotateBasisAntiClockTwice(basisVertex_v);
                        m_auxPatches[i].setNumberOfRotation(2);
                        //gsInfo << "Patch: " << m_patchesAroundVertex[i] << " with side " << m_vertexIndices[i]  << " rotated twice anticlockwise\n";
                        break;
                    case 3:
                        m_auxPatches[i].rotateBasisAntiClock(basisVertex_u);
                        m_auxPatches[i].rotateBasisAntiClock(basisVertex_v);
                        m_auxPatches[i].setNumberOfRotation(1);
                        //gsInfo << "Patch: " << m_patchesAroundVertex[i] << " with side " << m_vertexIndices[i]  << " rotated anticlockwise\n";
                        break;
                    case 2:
                        m_auxPatches[i].rotateBasisClock(basisVertex_u);
                        m_auxPatches[i].rotateBasisClock(basisVertex_v);
                        m_auxPatches[i].setNumberOfRotation(-1);
                        //gsInfo << "Patch: " << m_patchesAroundVertex[i] << " with side " << m_vertexIndices[i]  << " rotated clockwise\n";
                        break;
                }

            result[i][0] = basisVertex_u;
            result[i][1] = basisVertex_v;

            //gsWriteParaview(basisVertex_u,"basisUVertexRot"+util::to_string(i),2000);
            //gsWriteParaview(basisVertex_v,"basisVVertexRot"+util::to_string(i),2000);
        }
    }

    void checkOrientation(size_t i)
    {
        if (m_auxPatches[i].getPatch().orientation() == -1)
        {
            m_auxPatches[i].swapAxis();
            //gsInfo << "Changed axis on patch: " << m_patchesAroundVertex[i] << " with side " << m_vertexIndices[i] << "\n";
        }
    }

    void interpolateBasisVertex(size_t i, gsApproxGluingData<d, T> & approxGluingData, index_t corner,
                                std::vector<index_t> & sideContainer, real_t & sigma, gsMultiPatch<> & result_1)
    {

        // Collect the needed basis
        gsTensorBSplineBasis<d, T> basis_vertex = m_auxPatches[i].getC1BasisRotated().getVertexBasis(corner);

        std::vector<gsBSplineBasis<T>> basis_plus(2);
        std::vector<gsBSplineBasis<T>> basis_minus(2);
        std::vector<gsBSplineBasis<T>> basis_geo(2);

        std::vector<bool> kindOfEdge(2);

        for (size_t dir = 0; dir < sideContainer.size(); ++dir)
        {
            index_t localdir = m_auxPatches[i].getMapIndex(sideContainer[dir]) < 3 ? 1 : 0;

            basis_plus[localdir] = m_auxPatches[i].getC1BasisRotated().getBasisPlus(sideContainer[dir]);
            basis_minus[localdir] = m_auxPatches[i].getC1BasisRotated().getBasisMinus(sideContainer[dir]);
            basis_geo[localdir] = m_auxPatches[i].getC1BasisRotated().getBasisGeo(sideContainer[dir]);

            kindOfEdge[localdir] = m_auxPatches[i].getC1BasisRotated().isInterface(sideContainer[dir]);

        }

        gsGeometry<T> & geo = m_auxPatches[i].getPatch();

        // Points to interpolate at (Greville points):
        gsMatrix<> points = basis_vertex.anchors();

        // Computing the basis functions at the vertex
        gsMatrix<> Phi(6,6);
        Phi.setIdentity();

        Phi.row(1) *= sigma;
        Phi.row(2) *= sigma;
        Phi.row(3) *= sigma * sigma;
        Phi.row(4) *= sigma * sigma;
        Phi.row(5) *= sigma * sigma;

        // Computing c, c+ and c-
        // Point zero
        gsMatrix<> zero;
        zero.setZero(2,1);

        std::vector<gsMatrix<>> c_0, c_1;
        std::vector<gsMatrix < >> c_0_plus, c_1_plus, c_2_plus;
        std::vector<gsMatrix < >> c_0_plus_deriv, c_1_plus_deriv, c_2_plus_deriv;
        std::vector<gsMatrix < >> c_0_minus, c_1_minus;
        for (index_t i = 0; i < 2; i++) // i == 0 == u , i == 1 == v
        {
            gsMatrix<> b_0, b_1;
            gsMatrix<> b_0_plus, b_1_plus, b_2_plus;
            gsMatrix<> b_0_plus_deriv, b_1_plus_deriv, b_2_plus_deriv;
            gsMatrix<> b_0_minus, b_1_minus;

            gsBSplineBasis<T> bsp_temp = dynamic_cast<gsBSplineBasis<> & >(basis_geo[1-i]);
            real_t p = bsp_temp.degree();
            real_t h_geo = bsp_temp.knots().at(p + 1);

            basis_geo[1-i].evalSingle_into(0, points.row(i),b_0); // first
            basis_geo[1-i].evalSingle_into(1, points.row(i),b_1); // second

            basis_plus[i].evalSingle_into(0, points.row(i),b_0_plus);
            basis_plus[i].evalSingle_into(1, points.row(i),b_1_plus);
            basis_plus[i].evalSingle_into(2, points.row(i),b_2_plus);

            basis_plus[i].derivSingle_into(0, points.row(i),b_0_plus_deriv);
            basis_plus[i].derivSingle_into(1, points.row(i),b_1_plus_deriv);
            basis_plus[i].derivSingle_into(2, points.row(i),b_2_plus_deriv);

            basis_minus[i].evalSingle_into(0, points.row(i),b_0_minus);
            basis_minus[i].evalSingle_into(1, points.row(i),b_1_minus);

            c_0.push_back(b_0 + b_1);
            c_1.push_back((h_geo / p) * b_1);

            c_0_minus.push_back(b_0_minus + b_1_minus);
            c_1_minus.push_back(h_geo/ (p-1) * b_1_minus);

            gsMatrix<> der_b_1_plus_0, der2_b_1_plus_0, der2_b_2_plus_0;
            basis_plus[i].derivSingle_into(1, zero.row(i), der_b_1_plus_0);
            basis_plus[i].deriv2Single_into(1, zero.row(i), der2_b_1_plus_0);
            basis_plus[i].deriv2Single_into(2, zero.row(i), der2_b_2_plus_0);

            real_t factor_c_1_plus = 1/der_b_1_plus_0(0,0);
            real_t factor2_c_1_plus = -der2_b_1_plus_0(0,0)/(der_b_1_plus_0(0,0)*der2_b_2_plus_0(0,0));
            real_t factor_c_2_plus = 1/der2_b_2_plus_0(0,0);

            c_0_plus.push_back(b_0_plus + b_1_plus + b_2_plus);
            c_1_plus.push_back(factor_c_1_plus * b_1_plus + factor2_c_1_plus * b_2_plus);
            c_2_plus.push_back(factor_c_2_plus * b_2_plus );

            c_0_plus_deriv.push_back(b_0_plus_deriv + b_1_plus_deriv + b_2_plus_deriv);
            c_1_plus_deriv.push_back(factor_c_1_plus * b_1_plus_deriv + factor2_c_1_plus * b_2_plus_deriv);
            c_2_plus_deriv.push_back(factor_c_2_plus * b_2_plus_deriv);
        }

        std::vector<gsMatrix<>> alpha, beta, alpha_0, beta_0, alpha_deriv, beta_deriv;

        gsMatrix < T > temp_mat;
        if (kindOfEdge[0])
        {
            approxGluingData.alphaS(0).eval_into(points.row(0),temp_mat); // 1-dir == PatchID
            alpha.push_back(temp_mat); // u

            approxGluingData.alphaS(0).eval_into(zero.row(0),temp_mat); // 1-dir == PatchID
            alpha_0.push_back(temp_mat); // u

            approxGluingData.alphaS(0).deriv_into(zero.row(0),temp_mat); // 1-dir == PatchID
            alpha_deriv.push_back(temp_mat); // u

            approxGluingData.betaS(0).eval_into(points.row(0),temp_mat); // 1-dir == PatchID
            beta.push_back(temp_mat); // u

            approxGluingData.betaS(0).eval_into(zero.row(0),temp_mat); // 1-dir == PatchID
            beta_0.push_back(temp_mat); // u

            approxGluingData.betaS(0).deriv_into(zero.row(0),temp_mat); // 1-dir == PatchID
            beta_deriv.push_back(temp_mat); // u
        }
        else
        {
            temp_mat.setOnes(1, points.cols());
            alpha.push_back(temp_mat); // u

            temp_mat.setOnes(1, zero.cols());
            alpha_0.push_back(temp_mat); // u

            temp_mat.setZero(1, zero.cols());
            alpha_deriv.push_back(temp_mat); // u

            temp_mat.setZero(1, points.cols());
            beta.push_back(temp_mat); // u

            temp_mat.setZero(1, zero.cols());
            beta_0.push_back(temp_mat); // u

            temp_mat.setZero(1, zero.cols());
            beta_deriv.push_back(temp_mat); // u
        }



        if (kindOfEdge[1]) {
            approxGluingData.alphaS(1).eval_into(points.row(1), temp_mat); // 1-dir == PatchID
            alpha.push_back(temp_mat); // v

            approxGluingData.alphaS(1).eval_into(zero.row(0), temp_mat); // 1-dir == PatchID
            alpha_0.push_back(temp_mat); // v

            approxGluingData.alphaS(1).deriv_into(zero.row(0), temp_mat); // 1-dir == PatchID
            alpha_deriv.push_back(temp_mat); // v

            approxGluingData.betaS(1).eval_into(points.row(1), temp_mat); // 1-dir == PatchID
            beta.push_back(temp_mat); // v

            approxGluingData.betaS(1).eval_into(zero.row(0), temp_mat); // 1-dir == PatchID
            beta_0.push_back(temp_mat); // v

            approxGluingData.betaS(1).deriv_into(zero.row(0), temp_mat); // 1-dir == PatchID
            beta_deriv.push_back(temp_mat); // v
        }
        else
        {
            temp_mat.setOnes(1, points.cols());
            alpha.push_back(temp_mat); // u

            temp_mat.setOnes(1, zero.cols());
            alpha_0.push_back(temp_mat); // u

            temp_mat.setZero(1, zero.cols());
            alpha_deriv.push_back(temp_mat); // u

            temp_mat.setZero(1, points.cols());
            beta.push_back(temp_mat); // u

            temp_mat.setZero(1, zero.cols());
            beta_0.push_back(temp_mat); // u

            temp_mat.setZero(1, zero.cols());
            beta_deriv.push_back(temp_mat); // u
        }

        // Geo data:
        gsMatrix<T> geo_jac = geo.jacobian(zero);
        gsMatrix<T> geo_der2 = geo.deriv2(zero);

        // Compute dd^^(i_k) and dd^^(i_k-1)
        gsMatrix<> dd_ik_plus, dd_ik_minus;
        gsMatrix<> dd_ik_minus_deriv, dd_ik_plus_deriv;
        dd_ik_minus = -1/(alpha_0[0](0,0)) * (geo_jac.col(1) +
                                              beta_0[0](0,0) * geo_jac.col(0));

        dd_ik_plus = 1/(alpha_0[1](0,0)) * (geo_jac.col(0) +
                                            beta_0[1](0,0) * geo_jac.col(1));

        gsMatrix<> geo_deriv2_12(2,1), geo_deriv2_11(2,1), geo_deriv2_22(2,1);
        geo_deriv2_12.row(0) = geo_der2.row(2);
        geo_deriv2_12.row(1) = geo_der2.row(5);
        geo_deriv2_11.row(0) = geo_der2.row(0);
        geo_deriv2_11.row(1) = geo_der2.row(3);
        geo_deriv2_22.row(0) = geo_der2.row(1);
        geo_deriv2_22.row(1) = geo_der2.row(4);
        gsMatrix<> alpha_squared_u = alpha_0[0]*alpha_0[0];
        gsMatrix<> alpha_squared_v = alpha_0[1]*alpha_0[1];

        dd_ik_minus_deriv = -1/(alpha_squared_u(0,0)) * // N^2
                            ((geo_deriv2_12 + (beta_deriv[0](0,0) * geo_jac.col(0) +
                                               beta_0[0](0,0) * geo_deriv2_11))*alpha_0[0](0,0) -
                             (geo_jac.col(1) + beta_0[0](0,0) * geo_jac.col(0)) *
                             alpha_deriv[0](0,0));

        dd_ik_plus_deriv = 1/(alpha_squared_v(0,0)) *
                           ((geo_deriv2_12 + (beta_deriv[1](0,0) * geo_jac.col(1) +
                                              beta_0[1](0,0) * geo_deriv2_22))*alpha_0[1](0,0) -
                            (geo_jac.col(0) + beta_0[1](0,0) * geo_jac.col(1)) *
                            alpha_deriv[1](0,0));

/*
        gsInfo << "Transversal\n";
        if (kindOfEdge[1])
            gsInfo << "dd_ik_minus_deriv " << dd_ik_minus_deriv << "\n";
        if (kindOfEdge[0])
            gsInfo << "dd_ik_minus_deriv " << dd_ik_plus_deriv << "\n";

        gsInfo << geo_jac.col(0) << " : " << geo_jac.col(1) << "\n";
*/
        // Comupute d_(0,0)^(i_k), d_(1,0)^(i_k), d_(0,1)^(i_k), d_(1,1)^(i_k) ; i_k == 2
        std::vector<gsMatrix<>> d_ik;
        d_ik.push_back(Phi.col(0));
        d_ik.push_back(Phi.block(0,1,6,2) * geo_jac.col(0) ); // deriv into u
        d_ik.push_back(Phi.block(0,1,6,2) * geo_jac.col(1) ); // deriv into v
        d_ik.push_back((geo_jac(0,0) * Phi.col(3) + geo_jac(1,0) * Phi.col(4))*geo_jac(0,1) +
                       (geo_jac(0,0) * Phi.col(4) + geo_jac(1,0) * Phi.col(5))*geo_jac(1,1) +
                       Phi.block(0,1,6,1) * geo_der2.row(2) +
                       Phi.block(0,2,6,1) * geo_der2.row(5)); // Hessian

        // Compute d_(*,*)^(il,ik)
        std::vector<gsMatrix<>> d_ilik_minus, d_ilik_plus;
        d_ilik_minus.push_back(Phi.col(0));
        d_ilik_minus.push_back(Phi.block(0,1,6,2) * geo_jac.col(0));
        d_ilik_minus.push_back((geo_jac(0,0) * Phi.col(3) + geo_jac(1,0) * Phi.col(4))*geo_jac(0,0) +
                               (geo_jac(0,0) * Phi.col(4) + geo_jac(1,0) * Phi.col(5))*geo_jac(1,0) +
                               Phi.block(0,1,6,1) * geo_der2.row(0) +
                               Phi.block(0,2,6,1) * geo_der2.row(3));
        d_ilik_minus.push_back(Phi.block(0,1,6,2) * dd_ik_minus);
        d_ilik_minus.push_back((geo_jac(0,0) * Phi.col(3) + geo_jac(1,0) * Phi.col(4))*dd_ik_minus(0,0) +
                               (geo_jac(0,0) * Phi.col(4) + geo_jac(1,0) * Phi.col(5))*dd_ik_minus(1,0) +
                               Phi.block(0,1,6,1) * dd_ik_minus_deriv.row(0) +
                               Phi.block(0,2,6,1) * dd_ik_minus_deriv.row(1));

        d_ilik_plus.push_back(Phi.col(0));
        d_ilik_plus.push_back(Phi.block(0,1,6,2) * geo_jac.col(1));
        d_ilik_plus.push_back((geo_jac(0,1) * Phi.col(3) + geo_jac(1,1) * Phi.col(4))*geo_jac(0,1) +
                              (geo_jac(0,1) * Phi.col(4) + geo_jac(1,1) * Phi.col(5))*geo_jac(1,1) +
                              Phi.block(0,1,6,1) * geo_der2.row(1) +
                              Phi.block(0,2,6,1) * geo_der2.row(4));
        d_ilik_plus.push_back(Phi.block(0,1,6,2) * dd_ik_plus);
        d_ilik_plus.push_back((geo_jac(0,1) * Phi.col(3) + geo_jac(1,1) * Phi.col(4))*dd_ik_plus(0,0) +
                              (geo_jac(0,1) * Phi.col(4) + geo_jac(1,1) * Phi.col(5))*dd_ik_plus(1,0) +
                              Phi.block(0,1,6,1) * dd_ik_plus_deriv.row(0) +
                              Phi.block(0,2,6,1) * dd_ik_plus_deriv.row(1));

        for (index_t i = 0; i < 6; i++)
        {

            gsMatrix<> fValues = d_ilik_minus.at(0)(i,0) * (c_0_plus.at(0).cwiseProduct(c_0.at(1)) -
                                                       beta[0].cwiseProduct(c_0_plus_deriv.at(0).cwiseProduct(c_1.at(1)))) +
                            d_ilik_minus.at(1)(i,0) * (c_1_plus.at(0).cwiseProduct(c_0.at(1)) -
                                                       beta[0].cwiseProduct(c_1_plus_deriv.at(0).cwiseProduct(c_1.at(1)))) +
                            d_ilik_minus.at(2)(i,0) * (c_2_plus.at(0).cwiseProduct(c_0.at(1)) -
                                                       beta[0].cwiseProduct(c_2_plus_deriv.at(0).cwiseProduct(c_1.at(1)))) -
                            d_ilik_minus.at(3)(i,0) * alpha[0].cwiseProduct(c_0_minus.at(0).cwiseProduct(c_1.at(1))) -
                            d_ilik_minus.at(4)(i,0) * alpha[0].cwiseProduct(c_1_minus.at(0).cwiseProduct(c_1.at(1))); // f*_(ik-1,ik)

            //if (kindOfEdge[0])
            //rhsVals.at(i).setZero();

            //if (!kindOfEdge[1])
            fValues += d_ilik_plus.at(0)(i,0) * (c_0_plus.at(1).cwiseProduct(c_0.at(0)) -
                                                       beta[1].cwiseProduct(c_0_plus_deriv.at(1).cwiseProduct(c_1.at(0)))) +
                             d_ilik_plus.at(1)(i,0) * (c_1_plus.at(1).cwiseProduct(c_0.at(0)) -
                                                       beta[1].cwiseProduct(c_1_plus_deriv.at(1).cwiseProduct(c_1.at(0)))) +
                             d_ilik_plus.at(2)(i,0) * (c_2_plus.at(1).cwiseProduct(c_0.at(0)) -
                                                       beta[1].cwiseProduct(c_2_plus_deriv.at(1).cwiseProduct(c_1.at(0)))) +
                             d_ilik_plus.at(3)(i,0) * alpha[1].cwiseProduct(c_0_minus.at(1).cwiseProduct(c_1.at(0))) +
                             d_ilik_plus.at(4)(i,0) * alpha[1].cwiseProduct(c_1_minus.at(1).cwiseProduct(c_1.at(0))); // f*_(ik+1,ik)

            fValues -= d_ik.at(0)(i,0) * c_0.at(0).cwiseProduct(c_0.at(1)) + d_ik.at(2)(i,0) * c_0.at(0).cwiseProduct(c_1.at(1)) +
                             d_ik.at(1)(i,0) * c_1.at(0).cwiseProduct(c_0.at(1)) + d_ik.at(3)(i,0) * c_1.at(0).cwiseProduct(c_1.at(1)); // f*_(ik)

            // Returns a geometry with basis = tBasis
            // and coefficients being
            // computed as the interpolant of \a funct
            gsGeometry<>::uPtr interpolant = basis_vertex.interpolateAtAnchors(fValues);
            result_1.addPatch(*interpolant);
        }

    } // interpolateBasisVertex

    real_t computeSigma(const std::vector<size_t> & vertexIndices)
    {
        real_t sigma = 0;

        real_t p = 0;
        real_t h_geo = 1;
        for(size_t i = 0; i < m_auxPatches.size(); i++)
        {
            gsTensorBSplineBasis<2, real_t> bsp_temp = m_auxPatches[0].getC1BasisRotated().getVertexBasis(vertexIndices[i]);

            real_t p_temp = math::max(bsp_temp.degree(0), bsp_temp.degree(1));

            p = (p < p_temp ? p_temp : p);

            for(index_t j = 0; j < m_auxPatches[i].getPatch().parDim(); j++)
            {
                real_t h_geo_temp = bsp_temp.knot(j,p + 1);
                h_geo = (h_geo > h_geo_temp ? h_geo_temp : h_geo);
            }
        }

        gsMatrix<> zero;
        zero.setZero(2,1);
        for (size_t i = 0; i < m_auxPatches.size(); i++)
            sigma += gsMatrix<> (m_auxPatches[i].getPatch().deriv(zero)).lpNorm<Eigen::Infinity>();
        sigma *= h_geo/(m_auxPatches.size()*p);

        return (1 / sigma);
    }

    void computeKernel()
    {

        // TODO Boundary vertex with valence > 2
        gsMultiPatch<T> mp_vertex;
        for(size_t i = 0; i < m_patchesAroundVertex.size(); i++)
        {
            mp_vertex.addPatch(m_auxPatches[i].getPatch());
        }
        mp_vertex.computeTopology();

        index_t dim_mat = 0;
        std::vector<index_t> dim_u, dim_v, side, patchID;
        gsMatrix<> matrix_det(2,2), points(2,1);
        points.setZero();
        for(size_t np = 0; np < mp_vertex.nPatches(); np++)
        {
            if (mp_vertex.isBoundary(np,3)) // u
            {
                side.push_back(3);
                patchID.push_back(np);
                dim_u.push_back(basisVertexResult[np].basis(0).component(0).size());
                dim_v.push_back(basisVertexResult[np].basis(0).component(1).size());
                dim_mat += basisVertexResult[np].basis(0).component(0).size();

                matrix_det.col(0) = m_auxPatches[np].getPatch().jacobian(points).col(0); // u
            }
            if (mp_vertex.isBoundary(np,1)) // v
            {
                side.push_back(1);
                patchID.push_back(np);
                dim_u.push_back(basisVertexResult[np].basis(0).component(0).size());
                dim_v.push_back(basisVertexResult[np].basis(0).component(1).size());
                dim_mat += basisVertexResult[np].basis(0).component(1).size();

                matrix_det.col(1) = m_auxPatches[np].getPatch().jacobian(points).col(1); // u
            }
        }
        if (patchID.size() != 2)
            gsInfo << "Something went wrong \n";

        index_t dofsCorner = 3;
        if (matrix_det.determinant()*matrix_det.determinant() > 1e-15) // There is (numerically) a kink
            dofsCorner = 1;

        for(size_t np = 0; np < mp_vertex.nPatches(); np++)
            m_bases[m_patchesAroundVertex[np]].setNumDofsVertex(dofsCorner, m_vertexIndices[np]);

        if (m_optionList.getSwitch("info"))
            gsInfo << "Det: " << matrix_det.determinant() << "\n";

        gsMatrix<T> coefs_corner(dim_mat, 6);
        coefs_corner.setZero();

        index_t shift_row = 0;
        for (size_t bdy_index = 0; bdy_index < patchID.size(); ++bdy_index)
        {
            if (side[bdy_index] < 3) // v
            {
                for (index_t i = 0; i < dim_v[bdy_index]; ++i)
                {
                    for (index_t j = 0; j < 6; ++j)
                    {
                        T coef_temp = basisVertexResult[patchID[bdy_index]].patch(j).coef(i*dim_u[bdy_index], 0); // v = 0
                        if (coef_temp * coef_temp > 1e-25)
                            coefs_corner(shift_row+i, j) = coef_temp;
                    }
                }
                shift_row += dim_v[bdy_index];
            }
            else // u
            {
                for (index_t i = 0; i < dim_u[bdy_index]; ++i)
                {
                    for (index_t j = 0; j < 6; ++j)
                    {
                        T coef_temp = basisVertexResult[patchID[bdy_index]].patch(j).coef(i, 0); // v = 0
                        if (coef_temp * coef_temp > 1e-25)
                            coefs_corner(shift_row+i, j) = coef_temp;
                    }
                }
                shift_row += dim_u[bdy_index];
            }
        }

        real_t threshold = 1e-10;
        Eigen::FullPivLU<gsMatrix<>> KernelCorner(coefs_corner);
        KernelCorner.setThreshold(threshold);
        //gsInfo << "Coefs: " << coefs_corner << "\n";
        while (KernelCorner.dimensionOfKernel() < dofsCorner)
        {
            threshold += 1e-8;
            KernelCorner.setThreshold(threshold);
        }
        if (m_optionList.getSwitch("info"))
            gsInfo << "Dimension of Kernel: " << KernelCorner.dimensionOfKernel() << " With " << threshold << "\n";

        gsMatrix<> vertBas;
        vertBas.setIdentity(6, 6);

        gsMatrix<T> kernel = KernelCorner.kernel();

        size_t count = 0;
        while (kernel.cols() < 6) {
            kernel.conservativeResize(kernel.rows(), kernel.cols() + 1);
            kernel.col(kernel.cols() - 1) = vertBas.col(count);

            Eigen::FullPivLU<gsMatrix<>> ker_temp(kernel);
            ker_temp.setThreshold(1e-6);
            if (ker_temp.dimensionOfKernel() != 0) {
                kernel = kernel.block(0, 0, kernel.rows(), kernel.cols() - 1);
            }
            count++;
        }
        if (m_optionList.getSwitch("info"))
            gsInfo << "Kernel: " << kernel << "\n";

        for(size_t np = 0; np < m_patchesAroundVertex.size(); np++)
        {
            gsMultiPatch<> temp_result_0 = basisVertexResult[np];

            for (size_t j = 0; j < 6; ++j)
            {
                index_t dim_uv = temp_result_0.basis(j).size();
                gsMatrix<> coef_bf;
                coef_bf.setZero(dim_uv, 1);
                for (index_t i = 0; i < 6; ++i)
                    if (kernel(i, j) * kernel(i, j) > 1e-25)
                        coef_bf += temp_result_0.patch(i).coefs() * kernel(i, j);

                basisVertexResult[np].patch(j).setCoefs(coef_bf);
            }
        }

    }


protected:

    // Input
    gsMultiPatch<T> const & m_mp;
    C1BasisContainer & m_bases;

    const std::vector<size_t> & m_patchesAroundVertex;
    const std::vector<size_t> & m_vertexIndices;

    const gsOptionList & m_optionList;

    // Need for rotation, etc.
    C1AuxPatchContainer m_auxPatches;

    // Store temp solution
    std::vector<gsMultiPatch<T>> basisVertexResult;

}; // gsApproxC1Vertex

} // namespace gismo
