/** @file gsC1Basis.hpp

    @brief Provides declaration of Basis abstract interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#pragma once

#include <gsUnstructuredSplines/gsC1Basis.h>
#include <gsIO/gsFileData.h>
#include <gsCore/gsMultiPatch.h>
#include <gsNurbs/gsTensorBSplineBasis.h>

namespace gismo
{

template<short_t d,class T>
void gsC1Basis<d,T>::defaultOptions()
{
    m_options.addInt("test","a test option",0);
}

template<short_t d, class T>
void gsC1Basis<d,T>::print_spaces()
{
    // Some tests:
    for (index_t i = 0; i < 9; i++)
        if(basisG1Container[i].getMinCellLength() != basisG1Container[i].getMaxCellLength())
            gsInfo << "Different mesh-sizes is not implemented! \n";


    gsInfo << "-------------------------- Spaces for patch " << m_patchID << " --------------------------\n";
    gsInfo << "Interior space: S_1(" << basisG1Container[0].degree(0) << ", [";
    std::vector<index_t> kv_mult = basisG1Container[0].knots(0).multiplicities();
    for (size_t j = 1; j < kv_mult.size()-1; j++)
        gsInfo << kv_mult[j] << " ";
    gsInfo << "], " << basisG1Container[0].getMinCellLength() <<") ";
    gsInfo << "x S_2(" << basisG1Container[0].degree(1) << ", [";
    std::vector<index_t> kv_mult2 = basisG1Container[0].knots(1).multiplicities();
    for (size_t j = 1; j < kv_mult2.size()-1; j++)
        gsInfo << kv_mult2[j] << " ";
    gsInfo << "], " << basisG1Container[0].getMinCellLength() <<")\n";
    gsInfo << "\n------ Edge space:\n";
    for (index_t i = 1; i < 5; i++)
    {
        gsInfo << (kindOfEdge[i-1] ? "Interface-edge" : "Boundary-edge") << " space: S_1(" << basisG1Container[i].degree(0) << ", [";
        std::vector<index_t> kv_mult = basisG1Container[i].knots(0).multiplicities();
        for (size_t j = 1; j < kv_mult.size()-1; j++)
            gsInfo << kv_mult[j] << " ";
        gsInfo << "], " << basisG1Container[i].getMinCellLength() <<") ";
        gsInfo << "x S_2(" << basisG1Container[i].degree(1) << ", [";
        std::vector<index_t> kv_mult2 = basisG1Container[i].knots(1).multiplicities();
        for (size_t j = 1; j < kv_mult2.size()-1; j++)
            gsInfo << kv_mult2[j] << " ";
        gsInfo << "], " << basisG1Container[i].getMinCellLength() <<")\n";
    }
    gsInfo << "\n------ Vertex space:\n";
    for (index_t i = 5; i < 9; i++)
    {
        gsInfo << (kindOfVertex[i-5] == -1 ? "Boundary-vertex" : (kindOfVertex[i-5] == 0 ? "Internal-vertex" : "Interface-vertex"))
               << " space: S_1(" << basisG1Container[i].degree(0) << ", [";
        std::vector<index_t> kv_mult = basisG1Container[i].knots(0).multiplicities();
        for (size_t j = 1; j < kv_mult.size()-1; j++)
            gsInfo << kv_mult[j] << " ";
        gsInfo << "], " << basisG1Container[i].getMinCellLength() <<") ";
        gsInfo << "x S_2(" << basisG1Container[i].degree(1) << ", [";
        std::vector<index_t> kv_mult2 = basisG1Container[i].knots(1).multiplicities();
        for (size_t j = 1; j < kv_mult2.size()-1; j++)
            gsInfo << kv_mult2[j] << " ";
        gsInfo << "], " << basisG1Container[i].getMinCellLength() <<")\n";
    }
    gsInfo << "\n------ Plus/Minus space:\n";
    for (index_t i = 0; i < 4; i++)
    {
        gsInfo << "Plus space: S_1(" << basisPlusContainer[i].degree() << ", [";
        std::vector<index_t> kv_mult = basisPlusContainer[i].knots().multiplicities();
        for (size_t j = 1; j < kv_mult.size()-1; j++)
            gsInfo << kv_mult[j] << " ";
        gsInfo << "], " << basisPlusContainer[i].getMinCellLength() <<") ";
        gsInfo << "Minus space: S_1(" << basisMinusContainer[i].degree() << ", [";
        std::vector<index_t> kv_mult2 = basisMinusContainer[i].knots().multiplicities();
        for (size_t j = 1; j < kv_mult2.size()-1; j++)
            gsInfo << kv_mult2[j] << " ";
        gsInfo << "], " << basisMinusContainer[i].getMinCellLength() <<")\n";
    }
    gsInfo << "\n------ Gluing data/Geo space:\n";
    for (index_t i = 0; i < 4; i++)
    {
        gsInfo << "Gluing data space: S_1(" << basisGluingDataContainer[i].degree() << ", [";
        std::vector<index_t> kv_mult = basisGluingDataContainer[i].knots().multiplicities();
        for (size_t j = 1; j < kv_mult.size()-1; j++)
            gsInfo << kv_mult[j] << " ";
        gsInfo << "], " << basisGluingDataContainer[i].getMinCellLength() <<") ";
        gsInfo << "Geo space: S_1(" << basisGeoContainer[i].degree() << ", [";
        std::vector<index_t> kv_mult2 = basisGeoContainer[i].knots().multiplicities();
        for (size_t j = 1; j < kv_mult2.size()-1; j++)
            gsInfo << kv_mult2[j] << " ";
        gsInfo << "], " << basisGeoContainer[i].getMinCellLength() <<")\n";
    }
}

template<short_t d, class T>
void gsC1Basis<d,T>::uniformRefine()
{
    for (size_t i=0; i< basisG1Container.size(); ++i)
        basisG1Container[i].uniformRefine();

    for (size_t i=0; i< basisPlusContainer.size(); ++i)
        basisPlusContainer[i].uniformRefine();
    for (size_t i=0; i< basisMinusContainer.size(); ++i)
        basisMinusContainer[i].uniformRefine();
    for (size_t i=0; i< basisGeoContainer.size(); ++i)
        basisGeoContainer[i].uniformRefine();
    for (size_t i=0; i< basisGluingDataContainer.size(); ++i)
        basisGluingDataContainer[i].uniformRefine();
}

template<short_t d, class T>
void gsC1Basis<d,T>::swapAxis()
{
    for (size_t i=0; i< basisG1Container.size(); ++i)
    {
        gsTensorBSplineBasis<d, T> newTensorBasis(basisG1Container[i].knots(1),basisG1Container[i].knots(0));
        basisG1Container[i].swap(newTensorBasis);
    }
}

template<short_t d, class T>
void gsC1Basis<d,T>::init()
{
    // Col == number of coefs
    // Row == number of basis functions

    // Cols:
    for (size_t i = 0; i<basisG1Container.size(); ++i)
        if (basisG1Container[i].size() == 1)
            colContainer[i] = 0;
        else
            colContainer[i] = basisG1Container[i].size();

    // Inner basis functions
    {
        index_t dim_u = basisG1Container[0].component(0).size();
        index_t dim_v = basisG1Container[0].component(1).size();
        rowContainer[0] = (dim_u - 4) * (dim_v - 4);
    }

    // Interface basis functions
    for (index_t i = 0; i<4; ++i)
    {
        index_t dim_u = basisG1Container[i+1].component(0).size();
        index_t dim_v = basisG1Container[i+1].component(1).size();

        if (m_patches.isBoundary(m_patchID,i+1)) // +1 of side index
        {
            rowContainer[1+i] = math::max(basisPlusContainer[i].size()+basisMinusContainer[i].size() - 10, 0);
            kindOfEdge[i] = false; // bdy
        }
        else
        {
            rowContainer[1+i] = math::max(basisPlusContainer[i].size()+basisMinusContainer[i].size() - 10, 0);
            kindOfEdge[i] = true; // interface
        }

    }

    // Vertex basis functions
    for (index_t i = 0; i<4; ++i)
    {
        rowContainer[1+4+i] = valenceOfVertex[i];
    }

    if (m_options.getSwitch("info"))
    {
        gsInfo << "Patch: " << m_patchID << "\n";
        for (size_t i = 0; i < colContainer.size(); ++i)
            gsInfo << colContainer[i] << ", ";
        gsInfo << "\n";
        for (size_t i = 0; i < rowContainer.size(); ++i)
            gsInfo << rowContainer[i] << ", ";
        gsInfo << "\n";
        gsInfo << "Kind of Vertex\n";
        for (size_t i = 0; i < kindOfVertex.size(); ++i)
            gsInfo << kindOfVertex[i] << ", ";
        gsInfo << "\n";

    }
}

template<short_t d, class T>
void gsC1Basis<d,T>::matchWith(const boundaryInterface &bi, const gsBasis<T> &other, gsMatrix<int> &bndThis,
                               gsMatrix<int> &bndOther) const
{
    bndThis = this->boundaryOffset(bi.first().side(), 0);
    bndOther = other.boundaryOffset(bi.second().side(), 0);
}


template<short_t d, class T>
gsMatrix<int> gsC1Basis<d,T>::boundaryOffset(const boxSide &side, int offset) const
{
    if (offset > 1)
        gsInfo << "Offset > 1 is not implemented! \n";

    std::vector<boxCorner> containedCorners;
    side.getContainedCorners(d, containedCorners);

    short_t side_id = side.index();
    index_t num = 0;

    if (offset == 0)
        num = basisPlusContainer[side_id - 1].size() - 6; // -6 is shifting, the same for bdy and interface
    else if (offset == 1)
        num = basisMinusContainer[side_id - 1].size() - 4; // Interface
    else
        gsInfo << "Offset > 1 is not implemented! \n";

    num = num < 0 ? 0 : num; // if there are no bf at the interface

    gsMatrix<index_t> indizes(num , 1);

    index_t start = rowBegin(side_id); // The first num basis functions
    if (offset == 1)
        start += basisPlusContainer[side_id - 1].size() - 6; // second row shift

        index_t ii = 0;
    for (index_t i = start; i < start + num; i++, ii++) // Single basis function
        indizes(ii, 0) = i;

    return indizes;

    /*
    if (side.index() < 5) // Edge
    {
        short_t side_id = side.index();
        index_t num = 0;
        if (offset == 0)
        {
            index_t bdy_shift = 6;
            if (!kindOfEdge[side_id - 1])
                num = basisPlusContainer[side_id - 1].size() - bdy_shift; // Boundary
            else
                num = basisPlusContainer[side_id - 1].size() - 6; // Interface
        }
        else if (offset == 1)
        {
            index_t bdy_shift = 6;
            if (!kindOfEdge[side_id - 1])
                num = basisMinusContainer[side_id - 1].size() - bdy_shift; // Boundary might not used and wrong
            else
                num = basisMinusContainer[side_id - 1].size() - 4; // Interface
        } else
            gsInfo << "Offset > 1 is not implemented! \n";

        num = num < 0 ? 0 : num; // if there are no bf at the interface

        index_t ii = 0;
        gsMatrix<index_t> indizes(num , 1);
        //for(index_t of = 0;of<=offset;++of)
        {
            index_t start = rowBegin(side_id); // The first num basis functions

            if (offset == 1) {
                index_t bdy_shift = 6;
                if (!kindOfEdge[side_id - 1])
                    start += basisPlusContainer[side_id - 1].size() - bdy_shift; // Boundary
                else
                    start += basisPlusContainer[side_id - 1].size() - 6; // Interface
            }

            for (index_t i = start; i < start + num; i++, ii++) // Single basis function
                indizes(ii, 0) = i;
        }
        return indizes;
    }
    else if (side.index() > 4)
    {
        index_t corner_id = side.index(); // + 4 already included!
        if (offset == 0 && rows(corner_id ) != 0) {

            if (kindOfVertex[corner_id - 4 - 1] != 0)
            {
                index_t ii = 0;
                gsMatrix<index_t> indizes(valenceOfVertex[corner_id - 4 - 1] - numDofsVertex[corner_id - 4 - 1], 1);
                index_t corner_id = side.index(); // + 4 already included!
                index_t start = rowBegin(corner_id); // The first 3 basis functions
                for (index_t i = start + numDofsVertex[corner_id - 4 - 1];
                     i < start + valenceOfVertex[corner_id - 4 - 1]; i++, ii++) // Single basis function
                    indizes(ii, 0) = i;
                return indizes;
            }
            else
            {
                gsMatrix<index_t> null(1, 1);
                null(0, 0) = -1;
                return null;
            }
        }
        else if (offset == 1)
        {
            index_t ii = 0;
            gsMatrix<index_t> indizes(valenceOfVertex[corner_id - 4 - 1], 1);
            index_t start = rowBegin(corner_id); // The first 3 basis functions
            for (index_t i = start;
                 i < start + valenceOfVertex[corner_id - 4 - 1]; i++, ii++) // Single basis function
                indizes(ii, 0) = i;
            return indizes;
        }
    }
    gsMatrix<index_t> null(1, 1);
    null(0, 0) = -1;

    return null;
     */
}

template<short_t d,class T>
gsC1Basis<d,T>::gsC1Basis(  gsMultiPatch<T> & mp,
                                    index_t patchID )
:
m_patches(mp), m_patchID(patchID)
{
    //defaultOptions();
    //getOptions();

    // 9 Subspaces for the single patch
    basisG1Container.resize(9);

    // For each side:
    basisMinusContainer.resize(4);
    basisPlusContainer.resize(4);
    basisGeoContainer.resize(4);
    basisGluingDataContainer.resize(4);

    // For size
    rowContainer.resize(9);
    colContainer.resize(9);

    // For topology
    kindOfEdge.resize(4);
    kindOfVertex.resize(4);

    // For C1 Basis at vertex
    valenceOfVertex.resize(4);
    for (size_t i = 0; i<valenceOfVertex.size(); i++)
        valenceOfVertex[i] = 6; // to get for boundary vertex == 6

    // For boundary
    numDofsVertex.resize(4);
    for (size_t i = 0; i<numDofsVertex.size(); i++)
        numDofsVertex[i] = 1; // to get for boundary vertex == 1
}

// template<short_t d,class T>
// gsC1Basis<d,T>::gsC1Basis( gsC1Basis<d,T> other )
// :
// m_patches(other.m_patches),
// m_patchID(other.m_patchID),
// m_options(other.m_options)
// {

// }




} // namespace gismo
