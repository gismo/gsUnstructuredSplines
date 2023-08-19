/** @file gsC1SurfVertex.h

    @brief Creates the (approx.) C1 Vertex space.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Farahat & P. Weinmueller  
*/

#pragma once


#include <gsUnstructuredSplines/src/gsG1AuxiliaryPatch.h>
#include <gsUnstructuredSplines/src/gsC1SurfBasisVertex.h>


namespace gismo {
template<short_t d, class T>
class gsC1SurfVertex
{
    /// Shared pointer for gsC1SurfVertex
    typedef memory::shared_ptr<gsC1SurfVertex> Ptr;

    /// Unique pointer for gsC1SurfVertex
    typedef memory::unique_ptr<gsC1SurfVertex> uPtr;


public:
    /// Empty constructor
    ~gsC1SurfVertex() { }

    gsC1SurfVertex(gsMultiPatch<T> & mp, const std::vector<index_t > patchesAroundVertex, const std::vector<index_t > vertexIndices) : m_mp(mp)
    {
        for(size_t  i = 0; i < patchesAroundVertex.size(); i++)
        {
            auxGeom.push_back(gsG1AuxiliaryPatch<d,T>(mp.patch(patchesAroundVertex[i]), patchesAroundVertex[i]));
            auxVertexIndices.push_back(vertexIndices[i]);
            checkBoundary(mp, patchesAroundVertex[i], vertexIndices[i]);
        }
        sigma = 0.0;
    }

    gsMultiPatch<T> computeAuxTopology(){
        gsMultiPatch<T> auxTop;
        for(size_t  i = 0; i <  auxGeom.size(); i++)
        {
            auxTop.addPatch(auxGeom[i].getPatch());
        }
        auxTop.computeTopology();
        return auxTop;
    }

    void reparametrizeG1Vertex()
    {
        for(size_t  i = 0; i < auxGeom.size(); i++)
        {
            checkOrientation(i); // Check if the orientation is correct. If not, modifies vertex and edge vectors

            switch (auxVertexIndices[i])
            {
                case 1:
                    break;
                case 4:
                    auxGeom[i].rotateParamAntiClockTwice();
                    break;
                case 2:
                    auxGeom[i].rotateParamAntiClock();
                    break;
                case 3:
                    auxGeom[i].rotateParamClock();
                    break;
            }
        }
    }

    index_t kindOfVertex()
    {
        if(auxGeom.size() == 1)
            return -1; // Boundary vertex

        gsMultiPatch<T> top(computeAuxTopology());
        index_t  nInt = top.interfaces().size();
        if((index_t)auxGeom.size() == nInt)
            return 0; // Internal vertex
        else
            return 1; // Interface-Boundary vertex
    }

    void checkOrientation(index_t  i)
    {
        if (auxGeom[i].getPatch().orientation() == -1)
        {
            auxGeom[i].swapAxis();
            this->swapBdy(i); //Swap boundary edge bool-value

            // Swap vertices index after swapping axis
            if(auxVertexIndices[i] == 2)
                auxVertexIndices[i] = 3;
            else
            if(auxVertexIndices[i] == 3)
                auxVertexIndices[i] = 2;
        }
    }

    void computeSigma()
    {
        T p = 0;
        T h_geo = 0;
        for(size_t  i = 0; i < auxGeom.size(); i++)
        {
            gsTensorBSplineBasis<d, T> & bsp_temp = dynamic_cast<gsTensorBSplineBasis<d, T> & >(auxGeom[i].getPatch().basis());
            T p_temp = bsp_temp.maxDegree();
            p = (p < p_temp ? p_temp : p);

            for(index_t j = 0; j < auxGeom[i].getPatch().parDim(); j++)
            {
                T h_geo_temp = bsp_temp.component(j).knots().at(p + 2);
                h_geo = (h_geo < h_geo_temp ? h_geo_temp : h_geo);
            }
        }
        T val = auxGeom.size();

        gsMatrix<T> zero;
        zero.setZero(2,1);
        for (index_t i = 0; i < val; i++)
            sigma += auxGeom[i].getPatch().deriv(zero).template lpNorm<gsEigen::Infinity>();
        sigma *= h_geo/(val*p);
        sigma = 1 / sigma;
    }

    void checkBoundary(gsMultiPatch<T> & mpTmp, index_t   patchInd, index_t  sideInd)
    {
        std::vector<bool> tmp;
        switch (sideInd)
        {
            case 1: tmp.push_back(mpTmp.isBoundary(patchInd,3));
                tmp.push_back(mpTmp.isBoundary(patchInd,1));
                break;
            case 2: tmp.push_back(mpTmp.isBoundary(patchInd, 2));
                tmp.push_back(mpTmp.isBoundary(patchInd, 3));
                break;
            case 3: tmp.push_back(mpTmp.isBoundary(patchInd, 1));
                tmp.push_back(mpTmp.isBoundary(patchInd, 4));
                break;
            case 4: tmp.push_back(mpTmp.isBoundary(patchInd, 4));
                tmp.push_back(mpTmp.isBoundary(patchInd, 2));
                break;
            default:
                break;
        }
        isBdy.push_back(tmp);
    }

    void swapBdy(index_t  i)
    {
        bool tmp = isBdy[i][0];
        isBdy[i][0] = isBdy[i][1];
        isBdy[i][1] = tmp;
    }

    gsMatrix<T> computeBigSystemMatrix( index_t np)
    {
        gsMultiBasis<T> bas(auxGeom[np].getPatch());
        gsTensorBSplineBasis<d, T> & temp_L = dynamic_cast<gsTensorBSplineBasis<d, T> &>(bas.basis(0));
        index_t  dimU = temp_L.size(0);
        index_t  dimV = temp_L.size(1);

        gsMatrix<T> BigMatrix;
        BigMatrix.setZero( 2 * (dimU + dimV - 2),auxGeom[np].getG1Basis().nPatches());

        for(index_t  bf = 0; bf < auxGeom[np].getG1Basis().nPatches(); bf++)
        {
            for (index_t  i = 0; i < 2 * dimU; i++)
            {
                if (auxGeom[np].getG1BasisCoefs(bf).at(i) * auxGeom[np].getG1BasisCoefs(bf).at(i) > m_zero*m_zero)
                    BigMatrix(i, bf) = auxGeom[np].getG1BasisCoefs(bf).at(i);
            }

            for (index_t  i = 1; i < dimV - 1; i++)
            {
                for(index_t  j = i; j < i + 2; j++)
                {
                    if (auxGeom[np].getG1BasisCoefs(bf).at((i + 1) * dimU + j - i) * auxGeom[np].getG1BasisCoefs(bf).at((i + 1) * dimU + j - i) > m_zero*m_zero)
                        BigMatrix(i + j + (2 * dimU ) - 2, bf) = auxGeom[np].getG1BasisCoefs(bf).at((i + 1) * dimU + j - i);
                }
            }
        }
        return BigMatrix;
    }

    gsMatrix<T> computeSmallSystemMatrix( index_t np)
    {
        gsMultiBasis<T> bas(auxGeom[np].getPatch());
        gsTensorBSplineBasis<d, T> & temp_L = dynamic_cast<gsTensorBSplineBasis<d, T> &>(bas.basis(0));
        index_t  dimU = temp_L.size(0);
        index_t  dimV = temp_L.size(1);

        gsMatrix<T> SmallMatrix;
        SmallMatrix.setZero((dimU + dimV - 1),auxGeom[np].getG1Basis().nPatches());

        for(index_t  bf = 0; bf < auxGeom[np].getG1Basis().nPatches(); bf++)
        {
            for (index_t  i = 0; i < dimU; i++)
            {
                if (auxGeom[np].getG1BasisCoefs(bf).at(i) * auxGeom[np].getG1BasisCoefs(bf).at(i) > m_zero*m_zero)
                    SmallMatrix(i, bf) = auxGeom[np].getG1BasisCoefs(bf).at(i);
            }

            for (index_t  i = 1; i < dimV; i++)
            {
                if (auxGeom[np].getG1BasisCoefs(bf).at(i * dimU) * auxGeom[np].getG1BasisCoefs(bf).at(i * dimU) > m_zero*m_zero)
                    SmallMatrix(i + dimU -1, bf) = auxGeom[np].getG1BasisCoefs(bf).at(i * dimU);
            }
        }
        return SmallMatrix;
    }

    gsMatrix<T> leftBoundaryBigSystem(index_t np)
    {
        gsMultiBasis<T> bas(auxGeom[np].getPatch());
        gsTensorBSplineBasis<d, T> & temp_L = dynamic_cast<gsTensorBSplineBasis<d, T> &>(bas.basis(0));
        index_t  dimU = temp_L.size(0);
        index_t  dimV = temp_L.size(1);

        gsMatrix<T> BigMatrix;
        BigMatrix.setZero( 2 * dimV,6);

        for(index_t  bf = 0; bf < 6; bf++)
        {
            for (index_t  i = 0; i < dimV ; i++)
            {
                for(index_t  j = i; j < i + 2; j++)
                {
                    if (auxGeom[np].getG1BasisCoefs(bf).at(i * dimU + j - i) * auxGeom[np].getG1BasisCoefs(bf).at(i * dimU + j - i) > m_zero*m_zero)
                        BigMatrix(i + j, bf) = auxGeom[np].getG1BasisCoefs(bf).at(i  * dimU + j - i);
                    else
                        auxGeom[np].getG1BasisCoefs(bf).at(i) *= 0;
                }
            }
        }
        return BigMatrix;
    }

    gsMatrix<T> rightBoundaryBigSystem( index_t np)
    {
        gsMultiBasis<T> bas(auxGeom[np].getPatch());
        gsTensorBSplineBasis<d, T> & temp_L = dynamic_cast<gsTensorBSplineBasis<d, T> &>(bas.basis(0));
        index_t  dimU = temp_L.size(0);

        gsMatrix<T> BigMatrix;
        BigMatrix.setZero( 2 * dimU ,6);

        for(index_t  bf = 0; bf < 6; bf++)
        {
            for (index_t  i = 0; i < 2 * dimU; i++)
            {
                if (auxGeom[np].getG1BasisCoefs(bf).at(i) * auxGeom[np].getG1BasisCoefs(bf).at(i) > m_zero*m_zero)
                    BigMatrix(i, bf) = auxGeom[np].getG1BasisCoefs(bf).at(i);
                else
                    auxGeom[np].getG1BasisCoefs(bf).at(i) *= 0;
            }
        }
        return BigMatrix;
    }

    gsMatrix<T> leftBoundarySmallSystem( index_t np)
    {
        gsMultiBasis<T> bas(auxGeom[np].getPatch());
        gsTensorBSplineBasis<d, T> & temp_L = dynamic_cast<gsTensorBSplineBasis<d, T> &>(bas.basis(0));
        index_t  dimU = temp_L.size(0);
        index_t  dimV = temp_L.size(1);

        gsMatrix<T> SmallMatrix;
        SmallMatrix.setZero(dimV,6);

        for(index_t  bf = 0; bf < 6; bf++)
        {
            for (index_t  i = 0; i < dimV; i++)
            {
                if (auxGeom[np].getG1BasisCoefs(bf).at(i * dimU) * auxGeom[np].getG1BasisCoefs(bf).at(i * dimU) > m_zero*m_zero)
                    SmallMatrix(i, bf) = auxGeom[np].getG1BasisCoefs(bf).at(i * dimU);
                else
                    auxGeom[np].getG1BasisCoefs(bf).at(i * dimU) *= 0;
            }
        }
        return SmallMatrix;
    }

    gsMatrix<T> rightBoundarySmallSystem( index_t np)
    {
        gsMultiBasis<T> bas(auxGeom[np].getPatch());
        gsTensorBSplineBasis<d, T> & temp_L = dynamic_cast<gsTensorBSplineBasis<d, T> &>(bas.basis(0));
        index_t  dimU = temp_L.size(0);

        gsMatrix<T> SmallMatrix;
        SmallMatrix.setZero( dimU, 6);

        for(index_t  bf = 0; bf < 6; bf++)
        {
            for (index_t  i = 0; i < dimU; i++)
            {
                if (auxGeom[np].getG1BasisCoefs(bf).at(i ) * auxGeom[np].getG1BasisCoefs(bf).at(i ) > m_zero*m_zero)
                    SmallMatrix(i, bf) = auxGeom[np].getG1BasisCoefs(bf).at(i);
                else
                    auxGeom[np].getG1BasisCoefs(bf).at(i) *= 0;
            }
        }
        return SmallMatrix;
    }

    gsMatrix<T> bigInternalBoundaryPatchSystem( index_t np)
    {
        gsMultiBasis<T> bas(auxGeom[np].getPatch());
        gsTensorBSplineBasis<d, T> & temp_L = dynamic_cast<gsTensorBSplineBasis<d, T> &>(bas.basis(0));
        index_t  dimU = temp_L.size(0);

        gsMatrix<T> Matrix;
        Matrix.setZero( 3 ,6);

        for(index_t  bf = 0; bf < 6; bf++)
        {
            Matrix(0, bf) = auxGeom[np].getG1BasisCoefs(bf).at(0);
            Matrix(1, bf) = auxGeom[np].getG1BasisCoefs(bf).at(1);
            Matrix(2, bf) = auxGeom[np].getG1BasisCoefs(bf).at(dimU);

        }
        return Matrix;
    }

    gsMatrix<T> smallInternalBoundaryPatchSystem( index_t np)
    {
        gsMatrix<T> Matrix;
        Matrix.setZero( 1 ,6);

        for(index_t  bf = 0; bf < 6; bf++)
        {
            Matrix(0, bf) = auxGeom[np].getG1BasisCoefs(bf).at(0);
        }
        return Matrix;
    }

    std::pair<gsMatrix<T>, gsMatrix<T>> createSinglePatchSystem(index_t np)
    {
        if(isBdy[np][1] == 1)
            return std::make_pair(leftBoundaryBigSystem(np), leftBoundarySmallSystem(np));
        else
        {
            if (isBdy[np][0] == 1)
                return std::make_pair(rightBoundaryBigSystem(np), rightBoundarySmallSystem(np));
            else
                return std::make_pair(bigInternalBoundaryPatchSystem(np), smallInternalBoundaryPatchSystem(np));
        }
    }

    void checkValues(gsMatrix<T> & mat)
    {
        for(index_t bk = 0; bk < mat.cols(); bk++ )
        {
            for(index_t r=0; r < mat.rows(); r++)
            {
                if( abs( mat(r, bk) * mat(r, bk) )   < (m_zero*m_zero) )
                    mat(r, bk) = 0;
            }
        }
    }

    void addVertexBasis(gsMatrix<T> & basisV)
    {
        gsMatrix<T> vertBas;
        vertBas.setIdentity(auxGeom[0].getG1Basis().nPatches(), auxGeom[0].getG1Basis().nPatches());
        index_t  count = 0;
        index_t numBF = auxGeom[0].getG1Basis().nPatches();
        while (basisV.cols() < numBF)
        {
            basisV.conservativeResize(basisV.rows(), basisV.cols() + 1);
            basisV.col(basisV.cols() - 1) = vertBas.col(count);

            gsEigen::FullPivLU<gsMatrix<T>> ker(basisV);
            ker.setThreshold(1e-10);
            if (ker.dimensionOfKernel() != 0)
            {
                basisV = basisV.block(0, 0, basisV.rows(), basisV.cols() - 1);
            }
            count++;
        }
    }

    void addSmallKerBasis(gsMatrix<T> & basisV, gsMatrix<T> & smallK, index_t smallKDim)
    {
        for(index_t i=0; i < smallKDim; i++)
        {
            basisV.conservativeResize(basisV.rows(), basisV.cols() + 1);
            basisV.col(basisV.cols()-1) = smallK.col(i);

            gsEigen::FullPivLU<gsMatrix<T>> ker(basisV);
            ker.setThreshold(1e-10);
            if(ker.dimensionOfKernel() != 0)
            {
                basisV = basisV.block(0, 0, basisV.rows(), basisV.cols()-1);
            }
        }

    }

    std::pair<gsMatrix<T>, std::vector<index_t>> selectVertexBoundaryBasisFunction(gsMatrix<T> bigKernel, index_t bigKerDim, gsMatrix<T> smallKernel, index_t smallKerDim)
    {
        gsMatrix<T> basisVect;
        std::vector<index_t> numberPerType;

        numberPerType.push_back(bigKerDim); // Number of basis which has to be moved to the internal
        numberPerType.push_back(smallKerDim - bigKerDim); // Number of basis which are boundary function of FIRST TYPE
        numberPerType.push_back(bigKernel.cols() - smallKerDim); // Number of basis which are boundary function of SECOND TYPE

        if(bigKerDim != 0)
        {
            checkValues(bigKernel);
            checkValues(smallKernel);

            basisVect = bigKernel;
            addSmallKerBasis(basisVect, smallKernel, smallKerDim);
            addVertexBasis(basisVect);
        }
        else
        {
            if (smallKerDim != 0)
            {
                checkValues(smallKernel);

                basisVect = smallKernel;
                addVertexBasis(basisVect);
            }
            else
            {
                gsMatrix<T> vertBas;
                vertBas.setIdentity(auxGeom[0].getG1Basis().nPatches(), auxGeom[0].getG1Basis().nPatches());
                basisVect = vertBas;
            }
        }
        return std::make_pair(basisVect, numberPerType);
    }

    gsMatrix<T> selectGD(index_t i)
    {
        gsMatrix<T> coefs(4, 2);

        if( kindOfVertex() == 1 ) // If the boundary it´s along u and along v there is an interface (Right Patch) or viceversa
        {
            gsMultiPatch<T> tmp(this->computeAuxTopology());
            for(auto iter : tmp.interfaces())
            {
                if( i == iter.first().patch || i == iter.second().patch )
                {
                    gsMultiPatch<T> aux;
                    if( iter.first().index() == 1 )
                    {
                        aux.addPatch(tmp.patch(iter.first().patch));
                        aux.addPatch(tmp.patch(iter.second().patch));
                    }
                    else
                    {
                        aux.addPatch(tmp.patch(iter.second().patch));
                        aux.addPatch(tmp.patch(iter.first().patch));
                    }

                    aux.computeTopology();
                    gsMultiBasis<T> auxB(aux);
                    gsC1SurfGluingData<T> ret(aux, auxB);
                    gsMatrix<T> sol = ret.getSol();
                    gsMatrix<T> solBeta = ret.getSolBeta();

                    if( (isBdy[i][0] == 0) && (isBdy[i][1] == 0))
                    {
                        if ( (i == iter.first().patch && iter.first().index() == 3)
                             || (i == iter.second().patch && iter.second().index() == 3) )
                        {
                            coefs(0, 0) = sol(2, 0);
                            coefs(1, 0) = sol(3, 0);
                            coefs(2, 0) = solBeta(2, 0);
                            coefs(3, 0) = solBeta(3, 0);
                        }
                        else
                        if ( (i == iter.first().patch && iter.first().index() == 1)
                             || (i == iter.second().patch && iter.second().index() == 1))
                        {
                            coefs(0, 1) = sol(0, 0);
                            coefs(1, 1) = sol(1, 0);
                            coefs(2, 1) = solBeta(0, 0);
                            coefs(3, 1) = solBeta(1, 0);
                        }
                    }
                    else
                    {
                        if ((i == iter.first().patch && iter.first().index() == 3)
                            || (i == iter.second().patch && iter.second().index() == 3))
                        {
                            coefs(0, 0) = sol(2, 0);
                            coefs(0, 1) = 1;
                            coefs(1, 0) = sol(3, 0);
                            coefs(1, 1) = 1;
                            coefs(2, 0) = solBeta(2, 0);
                            coefs(2, 1) = 0;
                            coefs(3, 0) = solBeta(3, 0);
                            coefs(3, 1) = 0;
                        }
                        else if ((i == iter.first().patch && iter.first().index() == 1)
                                 || (i == iter.second().patch && iter.second().index() == 1))
                        {
                            coefs(0, 0) = 1;
                            coefs(0, 1) = sol(0, 0);
                            coefs(1, 0) = 1;
                            coefs(1, 1) = sol(1, 0);
                            coefs(2, 0) = 0;
                            coefs(2, 1) = solBeta(0, 0);
                            coefs(3, 0) = 0;
                            coefs(3, 1) = solBeta(1, 0);
                        }
                    }
                }
            }
        }
        else
        if( kindOfVertex() == -1 ) // Single patch corner
        {
            coefs.setZero();

            coefs(0, 0) = 1;
            coefs(0, 1) = 1;
            coefs(1, 0) = 1;
            coefs(1, 1) = 1;
        }
        else
        if( kindOfVertex() == 0 ) // Internal vertex -> Two interfaces
        {
            gsMultiPatch<T> tmp(this->computeAuxTopology());
            for(auto iter : tmp.interfaces())
            {
                if( (i == iter.first().patch) || (i == iter.second().patch) )
                {
                    gsMultiPatch<T> aux;
                    if( iter.first().index() == 1 )
                    {
                        aux.addPatch(tmp.patch(iter.first().patch));
                        aux.addPatch(tmp.patch(iter.second().patch));
                    }
                    else
                    {
                        aux.addPatch(tmp.patch(iter.second().patch));
                        aux.addPatch(tmp.patch(iter.first().patch));
                    }

                    aux.computeTopology();
                    gsMultiBasis<T> auxB(aux);
                    gsC1SurfGluingData<T> ret(aux, auxB);
                    gsMatrix<T> sol = ret.getSol();
                    gsMatrix<T> solBeta = ret.getSolBeta();

                    if ( (i == iter.first().patch && iter.first().index() == 3)
                         || (i == iter.second().patch && iter.second().index() == 3) )
                    {
                        coefs(0, 0) = sol(2, 0);
                        coefs(1, 0) = sol(3, 0);
                        coefs(2, 0) = solBeta(2, 0);
                        coefs(3, 0) = solBeta(3, 0);
                    }
                    else
                    if ( (i == iter.first().patch && iter.first().index() == 1)
                         || (i == iter.second().patch && iter.second().index() == 1))
                    {
                        coefs(0, 1) = sol(0, 0);
                        coefs(1, 1) = sol(1, 0);
                        coefs(2, 1) = solBeta(0, 0);
                        coefs(3, 1) = solBeta(1, 0);
                    }
                }
            }
        }
        return coefs;
    }

    void computeG1InternalVertexBasis()
    {

        m_zero = 1e-10;

        this->reparametrizeG1Vertex();
        this->computeSigma();

        std::vector<gsMultiPatch<T>> g1BasisVector;
        std::pair<gsMatrix<T>, std::vector<index_t>> vertexBoundaryBasis;

        gsMatrix<T> Phi(6, 6);
        Phi.setIdentity();

        Phi.col(1) *= sigma;
        Phi.col(2) *= sigma;
        Phi.col(3) *= sigma * sigma;
        Phi.col(4) *= sigma * sigma;
        Phi.col(5) *= sigma * sigma;

        gsMultiPatch<T> rotPatch;

        if (auxGeom[0].getPatch().parDim() + 1 == auxGeom[0].getPatch().targetDim())
        {
            gsMatrix<T> zero;
            zero.setZero(2, 1);
            gsMatrix<T> Jk = auxGeom[0].getPatch().jacobian(zero);
            gsMatrix<T> G = Jk.transpose() * Jk; // Symmetric
            gsMatrix<T> G_inv = G.cramerInverse(); // Symmetric


            gsMatrix<T> geoMapDeriv1 = auxGeom[0].getPatch()
                    .deriv(zero); // First derivative of the geometric mapping with respect to the parameter coordinates
            gsMatrix<T> geoMapDeriv2 = auxGeom[0].getPatch()
                    .deriv2(zero); // Second derivative of the geometric mapping with respect to the parameter coordinates

            //Computing the normal vector to the tangent plane along the boundary curve
//            gsVector<T,3> t1 = Jk.col(0);
//            gsVector<T,3> t2 = Jk.col(1);
//
//            gsVector<T,3> n = t1.cross(t2);
//
//            gsVector<T> normal = n.normalized();
//            n = n.normalized();
//            gsVector<T,3> z(0, 0, 1);
//
//            gsVector<T,3> rotVec = n.cross(z);
//            rotVec = rotVec.normalized();
//
//            T cos_t = n.dot(z) / (n.norm() * z.norm());
//            T sin_t = (n.cross(z)).norm() / (n.norm() * z.norm());

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

            for (size_t  np = 0; np < auxGeom.size(); np++)
            {
                gsMatrix<T> coeffPatch = auxGeom[np].getPatch().coefs();

                for (index_t i = 0; i < coeffPatch.rows(); i++)
                {
                    coeffPatch.row(i) =
                            (coeffPatch.row(i) - coeffPatch.row(0)) * R.transpose() + coeffPatch.row(0);
                }

                rotPatch.addPatch(auxGeom[np].getPatch());
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

            for (size_t  i = 0; i < auxGeom.size(); i++)
            {
                gsMatrix<T> gdCoefs(selectGD(i));
                gsMultiPatch<T> g1Basis;

                gsC1SurfBasisVertex<T> g1BasisVertex_0(rotPatch.patch(i), rotPatch.patch(i).basis(), isBdy[i], Phi, gdCoefs);

                g1BasisVertex_0.setG1BasisVertex(g1Basis);

                g1BasisVector.push_back(g1Basis);
                auxGeom[i].setG1Basis(g1Basis);
            }
        }
        else
        {
            for (size_t  i = 0; i < auxGeom.size(); i++)
            {
                gsMatrix<T> gdCoefs(selectGD(i));
                gsMultiPatch<T> g1Basis;

                gsC1SurfBasisVertex<T> g1BasisVertex_0(auxGeom[i].getPatch(), auxGeom[i].getPatch().basis(), isBdy[i], Phi, gdCoefs);

                g1BasisVertex_0.setG1BasisVertex(g1Basis);

                g1BasisVector.push_back(g1Basis);
                auxGeom[i].setG1Basis(g1Basis);
            }
        }

        // OLD
//        if (this->kindOfVertex() == 1) // Interface-Boundary vertex
//        {
//            gsMatrix<T> bigMatrix(0,0);
//            gsMatrix<T> smallMatrix(0,0);
//            for (index_t  i = 0; i < auxGeom.size(); i++)
//            {
//                std::pair<gsMatrix<T>, gsMatrix<T>> tmp;
//                tmp = createSinglePatchSystem(i);
//                index_t  row_bigMatrix = bigMatrix.rows();
//                index_t  row_smallMatrix = smallMatrix.rows();
//
//                bigMatrix.conservativeResize(bigMatrix.rows() + tmp.first.rows(), 6);
//                smallMatrix.conservativeResize(smallMatrix.rows() + tmp.second.rows(), 6);
//
//                bigMatrix.block(row_bigMatrix, 0, tmp.first.rows(), 6) = tmp.first;
//                smallMatrix.block(row_smallMatrix, 0, tmp.second.rows(), 6) = tmp.second;
//            }
//
//            gsEigen::FullPivLU<gsMatrix<T>> BigLU(bigMatrix);
//            gsEigen::FullPivLU<gsMatrix<T>> SmallLU(smallMatrix);
//            SmallLU.setThreshold(1e-10);
//            BigLU.setThreshold(1e-10);
//
////            if (!g1OptionList.getSwitch("neumann"))
////                dim_kernel = SmallLU.dimensionOfKernel();
////            else
////                dim_kernel = BigLU.dimensionOfKernel();
//
//            vertexBoundaryBasis = selectVertexBoundaryBasisFunction(BigLU.kernel(), BigLU.dimensionOfKernel(), SmallLU.kernel(), SmallLU.dimensionOfKernel());
//
//        }
//        else if(this->kindOfVertex() == -1) // Boundary vertex
//        {
//            gsEigen::FullPivLU<gsMatrix<T>> BigLU(computeBigSystemMatrix(0));
//            gsEigen::FullPivLU<gsMatrix<T>> SmallLU(computeSmallSystemMatrix(0));
//            SmallLU.setThreshold(1e-10);
//            BigLU.setThreshold(1e-10);
//
////            if (!g1OptionList.getSwitch("neumann"))
////                dim_kernel = SmallLU.dimensionOfKernel();
////            else
////                dim_kernel = BigLU.dimensionOfKernel();
//
//            vertexBoundaryBasis = selectVertexBoundaryBasisFunction(BigLU.kernel(), BigLU.dimensionOfKernel(), SmallLU.kernel(), SmallLU.dimensionOfKernel());
//
//        }
//
//        if (this->kindOfVertex() != 0)
//            for (index_t  i = 0; i < auxGeom.size(); i++)
//            {
//                gsMultiPatch<T> temp_mp_g1 = g1BasisVector[i];
//                for (index_t  bf = 0; bf < temp_mp_g1.nPatches(); bf++)
//                {
////                    gsDebug << "coeffbf: " << temp_mp_g1.patch(bf).coefs().transpose() << "\n";
//                    gsMatrix<T> coef_bf;
//                    coef_bf.setZero(temp_mp_g1.patch(bf).coefs().dim().first,1);
//                    for (index_t  lambda = 0; lambda < temp_mp_g1.nPatches(); lambda++)
//                        coef_bf += temp_mp_g1.patch(lambda).coefs() * vertexBoundaryBasis.first(lambda,bf);
//
//                    g1BasisVector[i].patch(bf).setCoefs(coef_bf);
//                }
//                auxGeom[i].parametrizeBasisBack(g1BasisVector[i]);
//            }
//        else
//            for (index_t  i = 0; i < auxGeom.size(); i++)
//                auxGeom[i].parametrizeBasisBack(g1BasisVector[i]);
        // END

        // NEW
        if (this->kindOfVertex() != 0)
            computeKernel(g1BasisVector);

        for (size_t  i = 0; i < auxGeom.size(); i++)
            auxGeom[i].parametrizeBasisBack(g1BasisVector[i]);
        // END

        for (size_t  i = 0; i < auxGeom.size(); i++)
        {
            for (size_t  ii = 0; ii < auxGeom[i].getG1Basis().nPatches(); ii++)
            {
                gsMatrix<T> coefs_temp;
                coefs_temp.setZero(auxGeom[i].getG1Basis().patch(ii).coefs().rows(),1);
                for (index_t j = 0; j < auxGeom[i].getG1Basis().patch(ii).coefs().rows(); j++)
                {
                    if (!math::almostEqual(auxGeom[i].getG1Basis().patch(ii).coefs().at(j)*auxGeom[i].getG1Basis().patch(ii).coefs().at(j),T(0)))
                        coefs_temp(j,0) = auxGeom[i].getG1Basis().patch(ii).coefs()(j,0);
                }
                auxGeom[i].getG1Basis().patch(ii).setCoefs(coefs_temp);
            }
            basisVertexResult.push_back(auxGeom[i].getG1Basis());
        }


        // just for plotting
//        std::string fileName;
//        std::string basename = "VerticesBasisFunctions" + util::to_string(auxGeom.size());
//        gsParaviewCollection collection(basename);
//
//        for (index_t  np = 0; np < auxGeom.size(); ++np)
//        {
//            if (basisVertexResult.size() != 0)
//                for (index_t  i = 0; i < basisVertexResult[np].nPatches(); ++i)
//                {
//                    fileName = basename + "_" + util::to_string(np) + "_" + util::to_string(i);
//                    gsField<T> temp_field(m_mp.patch(auxGeom[np].getGlobalPatchIndex()), basisVertexResult[np].patch(i));
//                    gsWriteParaview(temp_field, fileName, 5000);
//                    collection.addTimestep(fileName, i, "0.vts");
//
//                }
//        }
//        collection.save();
    }

    gsG1AuxiliaryPatch<d,T> & getSinglePatch(const index_t  i){ return auxGeom[i]; }

    std::vector<gsMultiPatch<T>> getBasis()
    { return basisVertexResult; }




    void computeKernel(std::vector<gsMultiPatch<T>> & g1BasisVector)
    {

        gsMultiPatch<T> mp_vertex;
        for(size_t  i = 0; i < auxGeom.size(); i++)
            mp_vertex.addPatch(auxGeom[i].getPatch());

        mp_vertex.computeTopology();

        index_t dim_mat = 0;
        std::vector<index_t> dim_u, dim_v, side;
        std::vector<index_t> dim_u_iFace;
        std::vector<size_t> patchID, patchID_iFace;
        gsMatrix<T> matrix_det(m_mp.targetDim(), m_mp.targetDim()), points(m_mp.parDim(),1);
        points.setZero();
        for(size_t np = 0; np < mp_vertex.nPatches(); np++)
        {
            if (mp_vertex.isBoundary(np,3)) // u
            {
                side.push_back(3);
                patchID.push_back(np);
                dim_u.push_back(auxGeom[np].getG1Basis().basis(0).component(0).size());
                dim_v.push_back(auxGeom[np].getG1Basis().basis(0).component(1).size());
                dim_mat += auxGeom[np].getG1Basis().basis(0).component(0).size();
                dim_mat += auxGeom[np].getG1Basis().basis(0).component(0).size();

                if(m_mp.parDim() == m_mp.targetDim()) // Planar
                    matrix_det.col(0) = auxGeom[np].getPatch().jacobian(points).col(0); // u
                else if(m_mp.parDim() + 1 == m_mp.targetDim()) // Surface
                {
                    gsMatrix<T> N, ev;
                    auxGeom[np].getPatch().jacobian_into(points, ev);
                    N.setZero(3,1);
                    N(0,0) = ev(1,0)*ev(2,1)-ev(2,0)*ev(1,1);
                    N(1,0) = ev(2,0)*ev(0,1)-ev(0,0)*ev(2,1);
                    N(2,0) = ev(0,0)*ev(1,1)-ev(1,0)*ev(0,1);
                    matrix_det.col(0) = ev.col(0);
                    matrix_det.col(2) = N;
                }
            }
            if (mp_vertex.isBoundary(np,1)) // v
            {
                side.push_back(1);
                patchID.push_back(np);
                dim_u.push_back(auxGeom[np].getG1Basis().basis(0).component(0).size());
                dim_v.push_back(auxGeom[np].getG1Basis().basis(0).component(1).size());
                dim_mat += auxGeom[np].getG1Basis().basis(0).component(1).size();
                dim_mat += auxGeom[np].getG1Basis().basis(0).component(1).size();

                if(m_mp.parDim() == m_mp.targetDim()) // Planar
                    matrix_det.col(1) = auxGeom[np].getPatch().jacobian(points).col(1); // u
                else if(m_mp.parDim() + 1 == m_mp.targetDim()) // Surface
                {
                    gsMatrix<T> N, ev;
                    auxGeom[np].getPatch().jacobian_into(points, ev);
                    N.setZero(3,1);
                    N(0,0) = ev(1,0)*ev(2,1)-ev(2,0)*ev(1,1);
                    N(1,0) = ev(2,0)*ev(0,1)-ev(0,0)*ev(2,1);
                    N(2,0) = ev(0,0)*ev(1,1)-ev(1,0)*ev(0,1);
                    matrix_det.col(1) = ev.col(1);
                }
            }
            if(mp_vertex.isInterface(patchSide(np,1)) && mp_vertex.isInterface(patchSide(np,3)))
            {
                dim_mat += 4;

                patchID_iFace.push_back(np);
                dim_u_iFace.push_back(auxGeom[np].getG1Basis().basis(0).component(0).size());
            }

        }
        GISMO_ASSERT(patchID.size()==2,"Something went wrong");

        index_t dofsCorner = 1;
        if (matrix_det.determinant()*matrix_det.determinant() > 1e-15) // There is (numerically) a kink
        {
            dofsCorner = 0;  // With Neumann
        }

        // gsDebug << "Det: " << matrix_det.determinant() << "\n";

        gsMatrix<T> coefs_corner(dim_mat, 6);
        coefs_corner.setZero();

        index_t shift_row = 0;
        for (size_t  bdy_index = 0; bdy_index < patchID.size(); ++bdy_index)
        {
            if (side[bdy_index] < 3) // v
            {
                index_t shift_row_neumann = dim_v[bdy_index];
                for (index_t i = 0; i < dim_v[bdy_index]; ++i)
                {
                    for (index_t j = 0; j < 6; ++j)
                    {
                        T coef_temp = auxGeom[patchID[bdy_index]].getG1Basis().patch(j).coef(i*dim_u[bdy_index], 0); // v = 0
                        if (!math::almostEqual(coef_temp * coef_temp,T(0)))
                            coefs_corner(shift_row+i, j) = coef_temp;

                        coef_temp = auxGeom[patchID[bdy_index]].getG1Basis().patch(j).coef(i*dim_u[bdy_index] +1, 0); // v = 0
                        if (!math::almostEqual(coef_temp * coef_temp,T(0)))
                            coefs_corner(shift_row_neumann + shift_row+i, j) = coef_temp;

                    }
                }
                shift_row += dim_v[bdy_index];
                shift_row += dim_v[bdy_index];
            }
            else // u
            {
                index_t shift_row_neumann = dim_u[bdy_index];
                for (index_t i = 0; i < dim_u[bdy_index]; ++i)
                {
                    for (index_t j = 0; j < 6; ++j)
                    {
                        T coef_temp = auxGeom[patchID[bdy_index]].getG1Basis().patch(j).coef(i, 0); // v = 0
                        if (!math::almostEqual(coef_temp * coef_temp,T(0)))
                            coefs_corner(shift_row+i, j) = coef_temp;

                        coef_temp = auxGeom[patchID[bdy_index]].getG1Basis().patch(j).coef(i+dim_u[bdy_index], 0); // v = 0
                        if (!math::almostEqual(coef_temp * coef_temp,T(0)))
                            coefs_corner(shift_row_neumann + shift_row+i, j) = coef_temp;

                    }
                }
                shift_row += dim_u[bdy_index];
                shift_row += dim_u[bdy_index];
            }
        }
        for (size_t iFace_index = 0; iFace_index < patchID_iFace.size(); ++iFace_index)
        {
            for (index_t i = 0; i < 2; ++i) // Only the first two
            {
                for (index_t j = 0; j < 6; ++j)
                {
                    T coef_temp = auxGeom[patchID_iFace[iFace_index]].getG1Basis().patch(j).coef(i, 0); // v = 0
                    if (!math::almostEqual(coef_temp * coef_temp,T(0)))
                        coefs_corner(shift_row + i, j) = coef_temp;
                    coef_temp = auxGeom[patchID_iFace[iFace_index]].getG1Basis().patch(j).coef(i + dim_u_iFace[iFace_index],
                                                                                      0); // v = 0
                    if (!math::almostEqual(coef_temp * coef_temp,T(0)))
                        coefs_corner(shift_row + 2 + i, j) = coef_temp; //  +2 bcs of the previous adding
                }
            }
            shift_row += 4;
        }

        gsMatrix<T> kernel;
        if (dofsCorner > 0)
        {
            T threshold = 1e-10;
            gsEigen::FullPivLU<gsMatrix<T>> KernelCorner(coefs_corner);
            KernelCorner.setThreshold(threshold);
            //gsDebug << "Coefs: " << coefs_corner << "\n";
            while (KernelCorner.dimensionOfKernel() < dofsCorner) {
                threshold += 1e-8;
                KernelCorner.setThreshold(threshold);
            }
            // gsDebug << "Dimension of Kernel: " << KernelCorner.dimensionOfKernel() << " With " << threshold << "\n";

            gsMatrix<T> vertBas;
            vertBas.setIdentity(6, 6);

            kernel = KernelCorner.kernel();

            index_t  count = 0;
            while (kernel.cols() < 6) {
                kernel.conservativeResize(kernel.rows(), kernel.cols() + 1);
                kernel.col(kernel.cols() - 1) = vertBas.col(count);

                gsEigen::FullPivLU<gsMatrix<T>> ker_temp(kernel);
                ker_temp.setThreshold(1e-6);
                if (ker_temp.dimensionOfKernel() != 0) {
                    kernel = kernel.block(0, 0, kernel.rows(), kernel.cols() - 1);
                }
                count++;
            }
        }
        else
            kernel.setIdentity(6, 6);

        // gsDebug << "NumDofs: " << dofsCorner << " with Kernel: \n" << kernel << "\n";

        for(size_t  np = 0; np < auxGeom.size(); np++)
        {
            gsMultiPatch<T> temp_result_0 = auxGeom[np].getG1Basis();

            for (index_t  j = 0; j < 6; ++j)
            {
                index_t dim_uv = temp_result_0.basis(j).size();
                gsMatrix<T> coef_bf;
                coef_bf.setZero(dim_uv, 1);
                for (index_t i = 0; i < 6; ++i)
                    if (!math::almostEqual(kernel(i, j) * kernel(i, j),T(0)))
                        coef_bf += temp_result_0.patch(i).coefs() * kernel(i, j);


                g1BasisVector[np].patch(j).setCoefs(coef_bf);
            }
        }

    }


protected:

    std::vector<gsG1AuxiliaryPatch<d,T>> auxGeom;
    std::vector<index_t > auxVertexIndices;
    std::vector< std::vector<bool>> isBdy;

    T sigma;
    index_t  dim_kernel;

    T m_zero;

    // Store temp solution
    std::vector<gsMultiPatch<T>> basisVertexResult;

    // only for plotting
    gsMultiPatch<T> & m_mp;

}; // gsC1SurfVertex

} // namespace gismo
