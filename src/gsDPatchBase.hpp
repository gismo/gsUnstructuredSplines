/** @file gsDPatchBase.hpp

    @brief Creates the D-Patch smoothing matrix.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gsNurbs/gsTensorNurbsBasis.h>
#include <gsHSplines/gsHTensorBasis.h>
#include <gsHSplines/gsTHBSpline.h>

#include <gsAssembler/gsExprHelper.h>
#include <gsAssembler/gsExprEvaluator.h>
#include <gsAssembler/gsAssembler.h>

namespace gismo
{

    template<short_t d,class T>
    void gsDPatchBase<d,T>::compute()
    {
        _initialize();
        _computeMapper();
        _computeSmoothMatrix();
        GISMO_ASSERT(this->_checkMatrix(m_matrix),"Mapper does not have column sum equal to 1");
        _makeTHB();
        _computeEVs();
        GISMO_ASSERT(this->_checkMatrix(m_matrix),"Mapper does not have column sum equal to 1");
        m_computed=true;
    }

    template<short_t d,class T>
    void gsDPatchBase<d,T>::defaultOptions()
    {
        m_options.addSwitch("SharpCorners","Reproduce sharp corners",true);
        m_options.addReal("SharpCornerTolerance","Sharp corner tolerance",1e-2);
        m_options.addSwitch("Verbose","Verbose output",false);
    }

    template<short_t d,class T>
    std::vector<bool> gsDPatchBase<d,T>::getSharpCorners(T tol) const
    {
        std::fill(m_vertCheck.begin(), m_vertCheck.end(), false);
        std::vector<bool> result(m_vertCheck.size());

        if (m_patches.empty())
        {
            gsWarn<<"Sharp corners should be detected, but no multi-patch object is assigned to this class\n";
            return result;
        }

        std::fill(result.begin(), result.end(), false);
        index_t cidx;
        patchCorner pcorner;
        for (size_t p=0; p<m_Bbases.nBases(); p++)
        {
            for (index_t c=1; c<5; c++)
            {
                cidx = _vertIndex(p,c);
                if (m_vertCheck.at(cidx))
                    continue;

                // Get the patchCorners which are on the boundary
                pcorner = patchCorner(p,c);
                std::vector<patchCorner> pcorners;
                m_topology.getCornerList(pcorner,pcorners);
                std::vector<std::pair<patchCorner,patchSide>> boundaries;
                for (std::vector<patchCorner>::iterator it=pcorners.begin(); it!=pcorners.end(); it++)
                {
                    std::vector<patchSide> psides;
                    it->getContainingSides(d,psides);
                    for (size_t s=0; s!=psides.size(); s++)
                    {
                        // the 0,k (k=0,1) DoF should be eliminated
                        if (!m_topology.isInterface(psides[s]))
                        {
                            boundaries.push_back(std::make_pair(*it,psides[s]));
                            continue;
                        }
                    }

                    // Mark the corner as passed
                    m_vertCheck[ _vertIndex(it->patch, it->corner()) ] = true;
                }

                if (boundaries.size()==0)
                    continue;
                else if (boundaries.size()==2)
                {
                    std::vector<gsVector<T>> normals;
                    for (std::vector<std::pair<patchCorner,patchSide>>::iterator it=boundaries.begin(); it!=boundaries.end(); it++)
                    {
                        gsMapData<T> md;
                        md.flags = NEED_OUTER_NORMAL;

                        // Get the coordinates of the corner
                        gsVector<bool> pars;
                        it->first.parameters_into(m_Bbases.domainDim(),pars); // get the parametric coordinates of the corner
                        gsMatrix<T> supp = m_Bbases.basis(it->first.patch).support();
                        gsVector<T> vec(supp.rows());
                        for (index_t r = 0; r!=supp.rows(); r++)
                            vec(r) = supp(r,pars(r));
                        // Assign the corner coordinates to the mapper
                        md.points = vec;
                        md.side = it->second;
                        m_patches.patch(it->first.patch).computeMap(md);

                        normals.push_back(md.outNormal(0).normalized());
                    }

                    if ( (std::abs(normals[0].transpose() * normals[1] - 1)) > tol )
                        for (std::vector<patchCorner>::iterator it=pcorners.begin(); it!=pcorners.end(); it++)
                            result[ _vertIndex(it->patch, it->corner()) ] = true;
                }
                else
                    GISMO_ERROR("Size of the stored boundary corners and sides must be 0 or 2, but is "<<boundaries.size());
            }

        }
        return result;
    }

    /*=====================================================================================
                                    Information functions
    =====================================================================================*/

    template<short_t d,class T>
    void gsDPatchBase<d,T>::mapperInfo() const
    {
        for (size_t p=0; p!=m_bases.nBases(); p++)
        {
            gsInfo<<"----------------------------------\n";
            gsInfo<<"patch "<<p<<"\n";
            index_t size = m_mapModified.patchSize(p);
            for (index_t k=0; k!=size; k++)
            {
                bool free = m_mapModified.is_free(k,p);
                std::string str = free ? "free" : "eliminated";
                gsInfo<<"DoF "<<k<<" is "<<str<<"\n";
            }
        }
    }

    template<short_t d,class T>
    void gsDPatchBase<d,T>::vertexInfo(patchCorner corner) const
    {
        std::pair<index_t,bool> data = this->_vertexData(corner);
        gsInfo<<"Patch "<<corner.patch<<", corner "<<corner<<" has valence "<<data.first<<" and is "<<(data.second ? "an interior vertex" : "a boundary vertex")<<"\n";

    }

    template<short_t d,class T>
    void gsDPatchBase<d,T>::sideInfo(patchSide side) const
    {
        gsInfo<<"Patch "<<side.patch<<", side "<<side<<" is "<<(m_topology.isBoundary(side) ? "a boundary side" : "an interface")<<"\n";
    }

    template<short_t d,class T>
    void gsDPatchBase<d,T>::sideInfo() const
    {
        gsInfo<<"**D-Patch Side info**\n";
        for(size_t i = 0;i<m_bases.nBases();++i)
            for(index_t j=1;j<=4;++j)
                sideInfo(patchSide(i,j));
    }

    template<short_t d,class T>
    void gsDPatchBase<d,T>::cornerInfo() const
    {
        gsInfo<<"**D-Patch Corner info**\n";
        for(size_t i = 0;i<m_bases.nBases();++i)
            for(index_t j=1;j<=4;++j)
                vertexInfo(patchCorner(i,j));
    }

    /*=====================================================================================
                                    Output functions
    =====================================================================================*/
    template<short_t d,class T>
    gsGeometry<T>* gsDPatchBase<d,T>::exportPatch(index_t patch, bool computeCoefs)
    {
        GISMO_ASSERT(m_computed,"The method has not been computed! Call compute().");
        ////////////////////////////////////////////////
        // This can be done more efficient!!
        // Do it once instead of for every patch
        ////////////////////////////////////////////////
        if (computeCoefs || m_coefs.rows()==0)
        {
            m_coefs = this->_preCoefficients(); // gets coefficients of the modified size
            m_coefs = m_matrix.transpose() * m_coefs; // maps to local size
        }

        ////////////////////////////////////////////////
        index_t size,offset = 0;
        for (index_t p=0; p!=patch; p++)
            offset += m_mapOriginal.patchSize(p);

        size = m_mapOriginal.patchSize(patch);
        gsMatrix<T> local = m_coefs.block(offset,0,size,m_coefs.cols());
        return m_bases[patch].makeGeometry( give(local) ).release();
    }

    // template<short_t d,class T>
    // gsMultiPatch<T> gsDPatchBase<d,T>::exportToPatches()
    // {
    //     m_coefs = this->_preCoefficients();
    //     m_coefs = m_matrix.transpose() * m_coefs;

    //     std::vector<gsGeometry<T> *> patches(m_bases.nBases());
    //     for (size_t p=0; p!=m_bases.nBases(); p++)
    //         patches[p]= this->exportPatch(p,false);

    //     return gsMultiPatch<T>(patches,m_topology.boundaries(),m_topology.interfaces());
    // }

    /*=====================================================================================
                                    Output functions
    =====================================================================================*/

    template<short_t d,class T>
    const index_t gsDPatchBase<d,T>::_indexFromSides(index_t index1, const patchSide side1, index_t index2, const patchSide side2)
    {
        /*
            Finds the index index1 away from side1 and index2 from side2
            index 1 is the index parallel to side 1
            index 2 is the index parallel to side 2
        */
        GISMO_ASSERT(side1.patch==side2.patch,"Sides must be from the same patch");
        GISMO_ASSERT(side1.side().direction()!=side2.side().direction(),"Sides must have different direction");
        index_t index;

        gsBasis<T> * basis = &m_bases.basis(side1.patch);

        gsVector<index_t> indices1 = static_cast<gsVector<index_t>>(basis->boundaryOffset(side1.side(),index2));

        index_t n = indices1.rows();
        if (side1.side()==1) //west
            if (side2.side().parameter()==0) //south
                index = indices1(index1);
            else
                index = indices1(n-1-index1); //north
        else if (side1.side()==2) //east
            if (side2.side().parameter()==0) //south
                index = indices1(index1);
            else
                index = indices1(n-1-index1); //north
        else if (side1.side()==3) //south
            if (side2.side().parameter()==0) //west
                index = indices1(index1);
            else
                index = indices1(n-1-index1); //east
        else if (side1.side()==4) //north
            if (side2.side().parameter()==0) //west
                index = indices1(index1);
            else
                index = indices1(n-1-index1); //east
        else
            GISMO_ERROR("Side unknown. index = "<<side1.side());
        return index;
    }

    template<short_t d,class T>
    const index_t gsDPatchBase<d,T>::_indexFromVert(const index_t index, const patchCorner corner, const patchSide side, const index_t offset) const
    {
        return _indexFromVert(&m_bases.basis(corner.patch),index, corner, side, offset);
    }

    template<short_t d,class T>
    const index_t gsDPatchBase<d,T>::_indexFromVert(const gsMultiBasis<T> & bases, const index_t index, const patchCorner corner, const patchSide side, const index_t offset) const
    {
        return _indexFromVert(&bases.basis(corner.patch),index, corner, side, offset);
    }

    template<short_t d,class T>
    const index_t gsDPatchBase<d,T>::_indexFromVert(const gsBasis<T> * basis, const index_t index, const patchCorner corner, const patchSide side, const index_t offset) const
    {
        const gsTensorBSplineBasis<d,T> * tbbasis;
        const gsTensorNurbsBasis<d,T> * tnbasis;
        const gsHTensorBasis<d,T> * thbasis;

        if ((index==0) && (offset==0))
        {
            // gsHTensorBasis<d,T> *basis = dynamic_cast<gsHTensorBasis<d,T>*>(&bases.basis(corner.patch));
            index_t idx = basis->functionAtCorner(corner.corner());
            return idx;
        }
        else
        {
            if ((tbbasis = dynamic_cast<const gsTensorBSplineBasis<d,T> * >(basis)))
                return _indexFromVert_impl(tbbasis,index,corner,side,offset);
            else if ((tnbasis = dynamic_cast<const gsTensorNurbsBasis<d,T> * >(basis)))
                return _indexFromVert_impl(&tnbasis->source(),index,corner,side,offset);
            else if ((thbasis = dynamic_cast<const gsHTensorBasis<d,T> * >(basis)))
                return _indexFromVert_impl(thbasis,index,corner,side,offset);
            else
                GISMO_ERROR("Basis type unknown");
        }
    }

    template<short_t d,class T>
    template<class U>
    typename util::enable_if<util::is_same<U,   const gsTensorBSplineBasis<d,T> *>::value,
                                                const index_t>::type
    gsDPatchBase<d,T>::_indexFromVert_impl(U basis, const index_t index, const patchCorner corner, const patchSide side, const index_t offset) const
    {
        /*
            Finds indices i in the direction of side away from the vertex
            if index = 0, the corner index is requested
        */

        index_t result = -1;

        gsVector<index_t,2> sizes;
        basis->size_cwise(sizes);
        if (side.side()==1) //west
        {
            GISMO_ASSERT(offset<sizes[0],"Offset is out of bounds");
            GISMO_ASSERT(index<sizes[1],"Index is out of bounds");
            if (corner.corner()==1)//southwest
                result = basis->index(offset,index);
            else if (corner.corner()==3) //northwest
                result = basis->index(offset,sizes[1]-1-index);
            else
                GISMO_ERROR(corner.corner() << " is not adjacent to side "<<side.side()<<"!");
        }
        else if (side.side()==2) //east
        {
            GISMO_ASSERT(offset<sizes[0],"Offset is out of bounds");
            GISMO_ASSERT(index<sizes[1],"Index is out of bounds");
            if (corner.corner()==2)//southeast
                result = basis->index(sizes[0]-1-offset,index);
            else if (corner.corner()==4) //northeast
                result = basis->index(sizes[0]-1-offset,sizes[1]-1-index);
            else
                GISMO_ERROR(corner.corner() << " is not adjacent to side "<<side.side()<<"!");
        }
        else if (side.side()==3) //south
        {
            GISMO_ASSERT(offset<sizes[1],"Offset is out of bounds");
            GISMO_ASSERT(index<sizes[0],"Index is out of bounds");
            if (corner.corner()==1)//southwest
                result = basis->index(index,offset);
            else if (corner.corner()==2) //southeast
                result = basis->index(sizes[0]-1-index,offset);
            else
                GISMO_ERROR(corner.corner() << " is not adjacent to side "<<side.side()<<"!");
        }
        else if (side.side()==4) //north
        {
            GISMO_ASSERT(offset<sizes[1],"Offset is out of bounds");
            GISMO_ASSERT(index<sizes[0],"Index is out of bounds");
            if (corner.corner()==3)//northwest
                result = basis->index(index,sizes[1]-1-offset);
            else if (corner.corner()==4) //northeast
                result = basis->index(sizes[0]-1-index,sizes[1]-1-offset);
            else
                GISMO_ERROR(corner.corner() << " is not adjacent to side "<<side.side()<<"!");

        }

        return result;
    }

    template<short_t d,class T>
    template<class U>
    typename util::enable_if<util::is_same<U,   const gsHTensorBasis<d,T> *>::value,
                                                const index_t>::type
    gsDPatchBase<d,T>::_indexFromVert_impl(U basis, const index_t index, const patchCorner corner, const patchSide side, const index_t offset) const
    {
        /*
            Finds indices i in the direction of side away from the vertex
            if index = 0, the corner index is requested
        */

        index_t result = -1;
        index_t level = basis->levelAtCorner(corner.corner());
        gsTensorBSplineBasis<d,T> tbasis;
        tbasis = basis->tensorLevel(level);
        result = _indexFromVert(&tbasis,index,corner,side,offset);
        result = basis->flatTensorIndexToHierachicalIndex(result,level);
        if (result!=-1) return result;

        GISMO_ASSERT(basis->maxLevel()!=0,"Basis must have more than one level!");
        GISMO_ASSERT(level!=0,"The level in the corner should not be 0!");

        // Find last active index
        index_t nActive_f = 0;
        result = 0;
        while (result!=-1)
        {
            result = _indexFromVert(&tbasis,nActive_f++,corner,side,offset);
            result = basis->flatTensorIndexToHierachicalIndex(result,level);
        }
        nActive_f--;
        // nActive_f is the number of functions from the corner that are active in the finest level

        // Now, we find the number of functions that are not active in the basis one level coarser
        index_t nInactive_c = 0;
        level--;
        result = -1;
        tbasis = basis->tensorLevel(level);
        while (result==-1) //&& nInactive_c < MAX!!
        {
            result = _indexFromVert(&tbasis,nInactive_c++,corner,side,offset);
            result = basis->flatTensorIndexToHierachicalIndex(result,level);
        }
        nInactive_c--;
        // nInactive_c number of functions from the corner that are inactive one level coarser

        // Now, the distance of the function with nInactive_c is  with respect to the corner, is nActive_f
        index_t tmpIdx = index - nActive_f + nInactive_c;
        result = _indexFromVert(&tbasis,tmpIdx,corner,side,offset);
        result = basis->flatTensorIndexToHierachicalIndex(result,level);
        GISMO_ASSERT(result!=-1,"Something went wrong");
        return result;
    }

    template<short_t d,class T>
    const std::pair<index_t,bool> gsDPatchBase<d,T>::_vertexData(const patchCorner corner) const
    {
        std::vector<patchCorner> corners;
        std::pair<index_t,bool> output;
        output.second = m_topology.getCornerList(corner,corners); // bool is true if interior vertex
        output.first = corners.size();
        return output;
    }

    template<short_t d,class T>
    void gsDPatchBase<d,T>::_getLowestCorners(std::vector<patchCorner> & pcorners, index_t n) const
    {
        struct {
            bool operator()(patchCorner a, patchCorner b) const { return a.patch < b.patch; }
        } customLess;

        // Sort
        std::sort(pcorners.begin(), pcorners.end(),customLess);
        // Resize; this vector are the indices we want to keep the DoFs of
        pcorners.resize(n);
    }

    template<short_t d,class T>
    void gsDPatchBase<d,T>::_removeLowestCorners(std::vector<patchCorner> & pcorners, index_t n) const
    {
        index_t size = pcorners.size();
        GISMO_ASSERT(n<=size,"You cannot remove more corners than there are actually stored in the container. Container size = "<<size<<" no. corners to be removed = "<<n);
        struct {
            bool operator()(patchCorner a, patchCorner b) const { return a.patch > b.patch; }
        } customGreater;

        // Sort
        std::sort(pcorners.begin(), pcorners.end(),customGreater);
        // Resize; this vector are the indices we want to keep the DoFs of
        pcorners.resize(size-n);
    }

    template<short_t d,class T>
    void gsDPatchBase<d,T>::_getLowestIndices(std::vector<std::pair<index_t,index_t>> & indices, index_t n) const
    {
        // indices are pairs with first=patchID, second=localID
        struct {
            bool operator()(std::pair<index_t,index_t> a, std::pair<index_t,index_t> b) const
            {
                return
                             (a.first < b.first)
                            ||
                            ((a.first == b.first) &&
                             (a.second < b.second)     );
            }
        } customLess;

        // Sort
        std::sort(indices.begin(), indices.end(),customLess);
        // Resize; this vector are the indices we want to keep the DoFs of
        indices.resize(n);
    }

    template<short_t d,class T>
    void gsDPatchBase<d,T>::_removeLowestIndices(std::vector<std::pair<index_t,index_t>> & indices, index_t n) const
    {
        index_t size = indices.size();
        GISMO_ASSERT(n<=size,"You cannot remove more corners than there are actually stored in the container. Container size = "<<size<<" no. corners to be removed = "<<n);
        // indices are pairs with first=patchID, second=localID
        struct {
            bool operator()(std::pair<index_t,index_t> a, std::pair<index_t,index_t> b) const
            {
                return
                             (a.first > b.first)
                            ||
                            ((a.first == b.first) &&
                             (a.second > b.second)     );
            }
        } customGreater;

        // Sort
        std::sort(indices.begin(), indices.end(),customGreater);
        // Resize; this vector are the indices we want to keep the DoFs of
        indices.resize(size-n);
    }

    // gets all interface DoFs on \a depth from the corner
    template<short_t d,class T>
    std::vector<std::pair<index_t,index_t>> gsDPatchBase<d,T>::_getInterfaceIndices(patchCorner pcorner, index_t depth, const gsMultiBasis<T> & mbasis) const
    {
        std::vector<std::pair<index_t,index_t>> result;
        std::vector<patchSide> csides(d);
        index_t index, patch;

        pcorner.getContainingSides(d,csides);
        patch = pcorner.patch;
        const gsBasis<T> * basis = &mbasis.basis(patch);

        if (depth==0)
        {
            // add the 0,0 one
            index = basis->functionAtCorner(pcorner);
            result.push_back(std::make_pair(patch,index));
        }
        else
        {
            for (index_t i = 0; i!=d; i++)
            {
                if ( m_topology.isBoundary(csides[i]) ) continue;
                index = _indexFromVert(mbasis,depth,pcorner,csides[i],0);
                result.push_back(std::make_pair(patch,index));
            }
        }
        return result;
    }

    // gets all interface DoFs on \a depth from the all corners surrounding pcorner
    template<short_t d,class T>
    std::vector<std::pair<index_t,index_t>> gsDPatchBase<d,T>::_getAllInterfaceIndices(patchCorner pcorner, index_t depth, const gsMultiBasis<T> & mbasis) const
    {
        std::vector<std::pair<index_t,index_t>> result, tmp;
        std::vector<patchCorner> pcorners;

        m_topology.getCornerList(pcorner,pcorners); // bool is true if interior vertex
        for (std::vector<patchCorner>::const_iterator it=pcorners.begin(); it!=pcorners.end(); it++)
        {
            tmp = this->_getInterfaceIndices(*it,depth,mbasis);
            result.insert(result.end(), tmp.begin(), tmp.end());
        }
        return result;
    }

    template<short_t d,class T>
    bool gsDPatchBase<d,T>::_checkMatrix(const gsSparseMatrix<T> & matrix) const // ! makes a deep copy (otherwise the contents of m_matrix get destroyed somehow...)
    {
        GISMO_ASSERT(matrix.cols()==matrix.outerSize(),"is the matrix ColMajor?");
        gsVector<T> colSums(matrix.cols());
        colSums.setZero();
        for (index_t i = 0; i<matrix.outerSize(); ++i)
            for (typename gsSparseMatrix<T>::iterator it = matrix.begin(i); it; ++it)
                colSums.at(i) += it.value();

        // TODO: Make tolerance depending on real_t
        colSums.array() -= 1;
        return (colSums.array() < 1e-8).all() && (colSums.array() > -1e-8).all(); // Lower?
    }

//     /*=====================================================================================
//                                     Construction functions
//     =====================================================================================*/


//     template<short_t d,class T>
//     void gsDPatchBase<d,T>::_makeTHB()
//     {
//         // prepare the geometry
//         std::vector<std::vector<patchCorner> > cornerLists;
//         // m_RefPatches.getEVs(cornerLists,true);

//         // get the corners that need refinement
//         std::vector<patchCorner> cornerList;
//         patchCorner pcorner;
//         index_t cidx;
//         for(size_t p = 0;p<m_RefPatches.nPatches();++p)
//         {
//             for(index_t c=1;c<=4;++c)
//             {
//                 pcorner=patchCorner(p,c);
//                 cidx = _vertIndex(p,c);
//                 bool C0 = m_C0s[cidx];
//                 bool isCycle = m_RefPatches.getCornerList(pcorner,cornerList);
//                 bool alreadyReached = false;
//                 for(size_t k = 0;k<cornerList.size();++k)
//                     if((size_t)cornerList[k].patch<p)
//                         alreadyReached = true;

//                 // add if
//                 // interior vertex with valence!=4
//                 // or
//                 // boundary vertex with valence > 2 (unless C0, then valence > 1)
//                 if(((isCycle&&cornerList.size()!=4)||((!isCycle)&&cornerList.size()>2-C0))&&!alreadyReached)
//                     cornerLists.push_back(cornerList);
//             }
//         }

//         if (cornerLists.size()!=0)
//         {
//             /// Change the coefficients
//             gsMatrix<T> coefs = this->freeCoefficients(); // gets coefficients of the modified size
//             coefs = m_matrix.transpose() * coefs; // maps to local size

//             this->setCoefficients(coefs,m_RefPatches);

//             /// Handle the EVs
//             std::vector< std::vector<index_t> > elVec(m_RefPatches.nPatches());
//             for (size_t v =0; v!=cornerLists.size(); v++)
//                 for (size_t c = 0; c!=cornerLists[v].size(); c++)
//                 {
//                     patchCorner corner = cornerLists[v].at(c);
//                     gsVector<bool> pars;
//                     corner.corner().parameters_into(m_RefPatches.parDim(),pars); // get the parametric coordinates of the corner
//                     gsMatrix<T> mat = pars.template cast<T>(); // cast to real coordinates

//                     gsMatrix<T> boxes(m_RefPatches.parDim(),2);
//                     boxes.col(0) << mat;
//                     boxes.col(1) << mat;

//                     gsHTensorBasis<d,T> *basis = dynamic_cast<gsHTensorBasis<d,T>*>(&m_RefPatches.basis(corner.patch));
//                     std::vector<index_t> elements = basis->asElements(boxes,0); // 0-ring

//                     elVec.at(corner.patch).insert(elVec.at(corner.patch).end(), elements.begin(), elements.end());

//                     // gsHTensorBasis<d,T> *basis = dynamic_cast<gsHTensorBasis<d,T>*>(&m_RefPatches.basis(corner.patch));

//                     // basis->refineElements(elements, m_tMatrix);
//                 }

//             gsSparseMatrix<T> tmp;
//             index_t rows = 0, cols = 0;
//             std::vector<gsEigen::Triplet<T,index_t>> tripletList;
//             for (size_t p=0; p!=m_RefPatches.nPatches(); p++)
//             {
//                 gsHTensorBasis<d,T> *basis = dynamic_cast<gsHTensorBasis<d,T>*>(&m_RefPatches.basis(p));
//                 std::vector< gsSortedVector< index_t > > xmat = basis->getXmatrix();

//                 m_RefPatches.patch(p).refineElements(elVec[p]);

//                 basis->transfer(xmat,tmp);

//                 for (index_t i = 0; i<tmp.outerSize(); ++i)
//                     for (typename gsSparseMatrix<T>::iterator it(tmp,i); it; ++it)
//                         tripletList.push_back(gsEigen::Triplet<T,index_t>(it.row()+rows,it.col()+cols,it.value()));

//                 rows += tmp.rows();
//                 cols += tmp.cols();
//             }

//             m_tMatrix.resize(rows,cols);
//             m_tMatrix.setFromTriplets(tripletList.begin(), tripletList.end());

//             m_tMatrix.makeCompressed();
//             m_bases = gsMultiBasis<T>(m_RefPatches);
//         }

//         // redefine the mappers
//         m_mapOriginal = gsDofMapper(m_bases);
//         m_mapOriginal.finalize();

//         // gsWriteParaview<T>(m_RefPatches,"mp_ref",1000,true);
//     }

//     template<short_t d,class T>
//     void gsDPatchBase<d,T>::_computeEVs()
//     {
//         /*
//             Our goal is to create three vectors c11, c12, c21 which all contain the
//             c11, c12 and c21 coefficients of the patches around the EV in the right order
//             (counter)-clockwise.
//         */

//         std::vector<std::vector<patchCorner> > cornerLists;
//         // m_topology.getEVs(cornerLists);

//         // get the corners that need refinement
//         std::vector<patchCorner> cornerList;
//         patchCorner pcorner;
//         index_t cidx;
//         for(size_t p = 0;p<m_RefPatches.nPatches();++p)
//         {
//             for(index_t c=1;c<=4;++c)
//             {
//                 pcorner=patchCorner(p,c);
//                 cidx = _vertIndex(p,c);
//                 bool C0 = m_C0s[cidx];
//                 bool isCycle = m_RefPatches.getCornerList(pcorner,cornerList);
//                 bool alreadyReached = false;
//                 for(size_t k = 0;k<cornerList.size();++k)
//                     if((size_t)cornerList[k].patch<p)
//                         alreadyReached = true;

//                 // add if
//                 // interior vertex with valence!=4
//                 // or
//                 // boundary vertex with valence > 2 (unless C0, then valence > 1)
//                 if(((isCycle&&cornerList.size()!=4)||((!isCycle)&&cornerList.size()>2-C0))&&!alreadyReached)
//                     cornerLists.push_back(cornerList);
//             }
//         }

//         std::fill(m_vertCheck.begin(), m_vertCheck.end(), false);
//         if (cornerLists.size()!=0)
//         {
//             m_matrix = m_matrix * m_tMatrix.transpose();

//             std::vector<patchCorner> pcorners;
//             patchCorner pcorner;
//             gsMatrix<T> Cg;         // coefficients
//             gsMatrix<T> ub;         // baricentric coordinates
//             gsMatrix<index_t> uind; // column indices of baricentric coordinates
//             index_t cidx;

//             for (std::vector<std::vector<patchCorner> >::iterator it=cornerLists.begin(); it!=cornerLists.end(); it++)
//             {

//                 std::vector<patchCorner> pcorners = *it;
//                 pcorner = it->at(0);
//                 cidx = _vertIndex(pcorner.patch,pcorner.corner());
//                 if (m_vertCheck.at(cidx))
//                     continue;

//                 std::pair<index_t,bool> vdata = _vertexData(pcorner); // corner c

//                 // get the triangle
//                 gsMatrix<T> Cg;
//                 std::tie(Cg,ub,uind) = _makeTriangle(pcorner);

//                 // The barycentric coordinates will be attached to the matrix rows corresponding to the 0,0 corners (and the three lowest patch index corners whenever valence > 3)
//                 // We use _getLowestIndices such that the corners are assigned to increasing patch corners
//                 // We need the index on the old basis (the unrefined basis), because we plug it into the mapModified (which maps the local DoFs to the global ones)
//                 std::vector<std::pair<index_t,index_t>> indices, tmp;
//                 if (vdata.first==2)
//                 {
//                     // These are two indices
//                     indices  = _getAllInterfaceIndices(pcorner,0,m_Bbases);
//                     tmp      = _getAllInterfaceIndices(pcorner,1,m_Bbases);
//                     _getLowestIndices(tmp,1);
//                     indices.push_back(tmp[0]);
//                 }
//                 else
//                 {
//                     indices  = _getAllInterfaceIndices(pcorner,0,m_Bbases);
//                     _getLowestIndices(indices,3);
//                 }


//                 std::vector<index_t> rowIndices;
//                 rowIndices.reserve(3);
//                 for (std::vector<std::pair<index_t,index_t>>::iterator it=indices.begin(); it!=indices.end(); it++)
//                 {
//                     // We need the index on the old basis (the unrefined basis), because we plug it into the mapModified (which maps the local DoFs to the global ones)
//                     GISMO_ASSERT(m_mapModified.is_free(it->second,it->first),"This DoF must be free! patch = "<<it->first<<"; index = "<<it->first);
//                     rowIndices.push_back(m_mapModified.index(it->second,it->first));
//                 }

//                 index_t rowIdx,colIdx;
//                 // set the colums related to the barycentric columns equal to zero
//                 for (index_t j=0; j!=ub.cols(); j++)
//                 {
//                     colIdx = uind(0,j);
//                     m_matrix.prune(
//                                     [&colIdx](index_t i, index_t j, T)
//                                     { return j!=colIdx; }
//                                     );
//                 }

//                 for (index_t i=0; i!=ub.rows(); i++)
//                     for (index_t j=0; j!=ub.cols(); j++)
//                     {
//                         rowIdx = rowIndices[i];
//                         colIdx = uind(0,j);
//                         m_matrix(rowIdx,colIdx) = ub(i,j);
//                     }

//                 for (size_t k = 0; k!=pcorners.size(); k++)
//                     m_vertCheck[ _vertIndex(pcorners[k].patch, pcorners[k].corner()) ] = true;
//             }
//             m_matrix.makeCompressed();
//         }
//     }
//

    template<short_t d,class T>
    void gsDPatchBase<d,T>::_initialize() // also initialize the mappers!
    {
        _initChecks();

        m_C0s.resize(m_nVerts);
        if (m_options.getSwitch("SharpCorners"))
            m_C0s = getSharpCorners(m_options.getReal("SharpCornerTolerance"));
        else
            std::fill(m_C0s.begin(), m_C0s.end(), false);

        _initTHB();
        _initBasis();
        _countDoFs();
        _initMappers();
        _initMatrix();
    }

    template<short_t d,class T>
    void gsDPatchBase<d,T>::_initChecks()
    {
        m_topology.checkConsistency();
        m_nSides = 2*m_topology.nInterfaces() + m_topology.nBoundary();
        m_nVerts = 4*m_Bbases.nBases();

        m_sideCheck.resize(m_nSides);
        std::fill(m_sideCheck.begin(), m_sideCheck.end(), false);
        m_vertCheck.resize(m_nVerts);
        std::fill(m_vertCheck.begin(), m_vertCheck.end(), false);
    }

    template<short_t d,class T>
    void gsDPatchBase<d,T>::_initTHB()
    {
        // Cast all patches of the mp object to THB splines
        for (size_t k=0; k!=m_Bbases.nBases(); ++k)
        {
            if ( (dynamic_cast<const gsTensorBSplineBasis<d,T> * > (&m_Bbases.basis(k))) )
                m_bases.addBasis(new gsTHBSplineBasis<d,T>(m_Bbases.basis(k)));
            else if ((dynamic_cast<const gsTHBSplineBasis<d,T> * > (&m_Bbases.basis(k))))
                m_bases.addBasis(m_Bbases.basis(k).clone());
            else
                gsWarn<<"No THB basis was constructed";
        }
        m_bases.setTopology(m_topology);
    }

    template<short_t d,class T>
    void gsDPatchBase<d,T>::_initBasis()
    {
    }

    template<short_t d,class T>
    void gsDPatchBase<d,T>::_initMappers()
    {
        m_mapModified = gsDofMapper(m_bases);
        m_mapOriginal = gsDofMapper(m_bases);
    }

    template<short_t d,class T>
    void gsDPatchBase<d,T>::_initMatrix()
    {
        m_matrix.resize(m_size,m_bases.totalSize());
    }

    template<short_t d,class T>
    void gsDPatchBase<d,T>::_whichHandled()
    {
        for (size_t p=0; p!=m_bases.nBases(); p++)
        {
            gsInfo<<"----------------------------------\n";
            gsInfo<<"patch "<<p<<"\n";
            for (size_t b=0; b!=m_mapOriginal.patchSize(p); b++)
            {
                index_t idx = m_mapModified.index(b,p);
                if (m_mapModified.is_free_index(idx))
                    gsInfo<<"basis function "<<b<<" (check="<<m_basisCheck[idx]<<") on patch "<<p<<" is "<<(m_basisCheck[idx] ? "":"not ")<<"handled\n";
                else
                    gsInfo<<"basis function "<<b<<" on patch "<<p<<" is "<<"eliminated\n";
            }
        }
    }

    template<short_t d,class T>
    void gsDPatchBase<d,T>::_performChecks(bool basis)
    {
        bool checkSides = std::all_of(m_sideCheck.begin(), m_sideCheck.end(), [](bool m_sideCheck) { return m_sideCheck; });
        GISMO_ENSURE(checkSides,"Not all sides are checked");
        bool checkVerts = std::all_of(m_vertCheck.begin(), m_vertCheck.end(), [](bool m_vertCheck) { return m_vertCheck; });
        GISMO_ENSURE(checkVerts,"Not all vertices are checked");

        if (!basis)
            return;

        bool checkBasis = std::all_of(m_basisCheck.begin(), m_basisCheck.end(), [](bool m_basisCheck) { return m_basisCheck; });
        GISMO_ENSURE(checkBasis,"Not all basis functions are checked");
    }

    template<short_t d,class T>
    void gsDPatchBase<d,T>::_resetChecks(bool basis)
    {
        m_sideCheck.resize(m_nSides);
        std::fill(m_sideCheck.begin(), m_sideCheck.end(), false);
        m_vertCheck.resize(m_nVerts);
        std::fill(m_vertCheck.begin(), m_vertCheck.end(), false);

        if (!basis)
            return;

        m_basisCheck.resize(m_size);
        std::fill(m_basisCheck.begin(), m_basisCheck.end(), false);

    }

    template<short_t d,class T>
    void gsDPatchBase<d,T>::_handleInteriorVertex(patchCorner pcorner, index_t valence)
    {
        std::vector<patchSide> psides;
        std::vector<patchCorner> corners;
        sparseEntry_t entries;

        boundaryInterface iface;
        patchSide otherSide;

        index_t colIdx, rowIdx;
        T weight;

        pcorner.getContainingSides(d,psides);

        index_t b11_p1 = _indexFromVert(1,pcorner,psides[0],1); // point 1,1 (does not matter which reference side is taken)
        rowIdx = m_mapModified.index(b11_p1,pcorner.patch);
        colIdx = m_mapOriginal.index(b11_p1,pcorner.patch);
        // Influence of 1,1 to itself
        weight = 1.;
        entries.push_back(std::make_tuple(rowIdx,colIdx,weight));

        m_topology.getCornerList(pcorner,corners);

        // The 1,1 corners of each patch will be given 0.5 weight in the interface handling, but in addition they will receive a 1/v weight from the 0,0 DoFs on each patch
        gsBasis<T> * tmpBasis;
        index_t index;
        for (std::vector<patchCorner>::iterator corn = corners.begin(); corn != corners.end(); ++corn)
        {
            tmpBasis = &m_bases.basis(corn->patch);
            index = tmpBasis->functionAtCorner(corn->corner());
            colIdx = m_mapOriginal.index(index,corn->patch);
            weight = 1./valence;
            entries.push_back(std::make_tuple(rowIdx,colIdx,weight));

            // m_matrix(rowIdx,colIdx) = 1./valence;
        }

        // Interface contributions
        for (std::vector<patchSide>::iterator side = psides.begin(); side != psides.end(); ++side)
        {
            GISMO_ENSURE(m_topology.getInterface(*side,iface),"Side must be an interface!");
            m_topology.getNeighbour(*side,otherSide);
            patchCorner otherCorner = iface.mapCorner(pcorner);

            index_t b10_p1 = _indexFromVert(1,pcorner,*side,0); // index from vertex pcorners[c] along side psides[0] with offset 0.
            index_t b10_p2 = _indexFromVert(1,otherCorner,otherSide,0); // point 0,1

            weight = 0.5;
            colIdx = m_mapOriginal.index(b10_p1,side->patch);
            entries.push_back(std::make_tuple(rowIdx,colIdx,weight));
            colIdx = m_mapOriginal.index(b10_p2,otherSide.patch);
            entries.push_back(std::make_tuple(rowIdx,colIdx,weight));
        }


        #pragma omp critical (handle_interior_vertex)
        {
            _pushAndCheck(entries);
            m_vertCheck[ _vertIndex(pcorner.patch, pcorner.corner()) ] = true;
        }
    }


    template<short_t d,class T>
    void gsDPatchBase<d,T>::_handleVertex(patchCorner pcorner)
    {

        if (m_vertCheck[ _vertIndex(pcorner.patch, pcorner.corner()) ])
        {
            // gsDebug<<"corner "<<pcorner.corner()<<" ("<<pcorner.patch<<") skipped!\n";
            return;
        }

        std::pair<index_t,bool> vdata = _vertexData(pcorner); // corner c

        bool C0 = m_C0s[_vertIndex(pcorner.patch,pcorner.corner())];
        if (!vdata.second) // boundary vertices
        {
            if (vdata.first==1)
                _handleRegularCorner(pcorner);
            else if (vdata.first==2 && !C0)
                _handleRegularBoundaryVertexSmooth(pcorner,vdata.first);
            else if (vdata.first==2 &&  C0)
                _handleRegularBoundaryVertexNonSmooth(pcorner,vdata.first);
            else if (vdata.first> 2 && !C0)
                _handleIrregularBoundaryVertexSmooth(pcorner,vdata.first);
            else if (vdata.first> 2 &&  C0)
                _handleIrregularBoundaryVertexNonSmooth(pcorner,vdata.first);
            else
                GISMO_ERROR("Something went wrong");
        } // end boundary vertices
        else // interior vertices
            _handleInteriorVertex(pcorner,vdata.first);
    }

    template<short_t d,class T>
    void gsDPatchBase<d,T>::_handleInterface(boundaryInterface iface)
    {
        if (m_sideCheck[ _sideIndex(iface.first().patch, iface.first().side()) ] || m_sideCheck[ _sideIndex(iface.second().patch, iface.second().side()) ])
        {
            // gsDebug<<"sides "<<iface.first().side()<<" ("<<iface.first().patch<<") and "<<iface.second().side()<<" ("<<iface.second().patch<<") skipped!\n";
            return;
        }

        std::vector<patchCorner> pcorners;

        std::vector<std::vector<index_t>> selectedIndices(2);
        std::vector<std::vector<index_t>> selectedOIndices(2);

        std::vector<gsBasis<T> *> basis(2);
        std::vector<gsMatrix<index_t>> indices(2); // interface indices
        std::vector<gsMatrix<index_t>> oindices(2); // interface indices

        sparseEntry_t entries;

        // get the bases belonging to both patches
        for (index_t p =0; p!=2; p++)
            basis[p] = &m_bases.basis(iface[p].patch);

        // this assumes the directions are handled correctly in matchWith (indices has the same direction as oindices)
        basis[0]->matchWith(iface,*basis[1],indices[0],indices[1],0);
        basis[0]->matchWith(iface,*basis[1],oindices[0],oindices[1],1);

        index_t np;
        for (index_t p =0; p!=2; p++)
        {
            np = 1-p; // not index p;

            iface[p].getContainedCorners(d,pcorners);
            for (index_t c =0; c!=2; c++)
            {
                selectedIndices[p].push_back(_indexFromVert(0,pcorners[c],iface[p],0)); // index from vertex pcorners[c] along side psides[0] with offset 0.
                selectedIndices[p].push_back(_indexFromVert(1,pcorners[c],iface[p],0)); // index from vertex pcorners[c] along side psides[0] with offset 0.

                selectedOIndices[p].push_back(_indexFromVert(0,pcorners[c],iface[p],1)); // index from vertex pcorners[c] along side psides[0] with offset 0.
                selectedOIndices[p].push_back(_indexFromVert(1,pcorners[c],iface[p],1)); // index from vertex pcorners[c] along side psides[0] with offset 0.
            }
            std::vector<index_t> allIndices(indices[p].data(), indices[p].data() + indices[p].rows() * indices[p].cols());
            std::vector<index_t> result;
            std::copy_if(allIndices.begin(), allIndices.end(), std::back_inserter(result),
                [&selectedIndices,&p] (index_t entry)
                {
                    std::vector<index_t>::const_iterator res = std::find(selectedIndices[p].begin(), selectedIndices[p].end(), entry);
                    return (res == selectedIndices[p].end());
                });
            indices[p] = gsAsMatrix<index_t>(result);

            std::vector<index_t> allOIndices(oindices[p].data(), oindices[p].data() + oindices[p].rows() * oindices[p].cols());
            result.clear();
            std::copy_if(allOIndices.begin(), allOIndices.end(), std::back_inserter(result),
                [&selectedOIndices,&p] (index_t entry)
                {
                    std::vector<index_t>::const_iterator res = std::find(selectedOIndices[p].begin(), selectedOIndices[p].end(), entry);
                    return (res == selectedOIndices[p].end());
                });
            oindices[p] = gsAsMatrix<index_t>(result);
        }

        GISMO_ASSERT(indices[0].size()==indices[1].size(),"Indices do not have the right size, indices[0].size()="<<indices[0].size()<<",indices[1].size()="<<indices[1].size());
        GISMO_ASSERT(oindices[0].size()==oindices[1].size(),"Offset indices do not have the right size, oindices[0].size()="<<oindices[0].size()<<",oindices[1].size()="<<oindices[1].size());

        index_t rowIdx,colIdx;
        T weight;
        // loop over adjacent patches and couple the DoFs.
        for (index_t p =0; p!= 2; p++)
        {
            np = 1-p; // not index p;
            for (index_t k=0; k!= indices[p].size(); k++ )
            {
                GISMO_ASSERT(m_mapModified.is_free(oindices[p].at(k),iface[p].patch),"Index "<<oindices[p].at(k)<<" on patch "<<iface[p].patch<<" is eliminated. Something went wrong?");
                rowIdx = m_mapModified.index(oindices[p].at(k),iface[p].patch);
                // rowIdx1 = m_mapOriginal.index(oindices[p].at(k),patches[p]);
                GISMO_ASSERT(m_mapOriginal.is_free(oindices[p].at(k),iface[p].patch),"Index is eliminated. Something went wrong?");
                colIdx = m_mapOriginal.index(oindices[p].at(k),iface[p].patch);
                // m_matrix(rowIdx,colIdx) = 1.0;
                weight = 1.0;
                entries.push_back(std::make_tuple(rowIdx,colIdx,weight));

                GISMO_ASSERT(m_mapOriginal.is_free(indices[p].at(k),iface[p].patch),"Index is eliminated. Something went wrong?");
                colIdx = m_mapOriginal.index(indices[p].at(k),iface[p].patch);
                // m_matrix(rowIdx,colIdx) = 0.5;
                weight = 0.5;
                entries.push_back(std::make_tuple(rowIdx,colIdx,weight));

                GISMO_ASSERT(m_mapOriginal.is_free(indices[np].at(k),iface[np].patch),"Index is eliminated. Something went wrong?");
                colIdx = m_mapOriginal.index(indices[np].at(k),iface[np].patch);
                weight = 0.5;
                // m_matrix(rowIdx,colIdx) = 0.5;
                entries.push_back(std::make_tuple(rowIdx,colIdx,weight));

                // m_basisCheck[rowIdx] = true;
            }
            // m_sideCheck[ _sideIndex(iface[p].patch, iface[p].side()) ] = true; // side finished
        }

        #pragma omp critical (handle_interface)
        {
            _pushAndCheck(entries);

            for (index_t p =0; p!= 2; p++)
                m_sideCheck[ _sideIndex(iface[p].patch, iface[p].side()) ] = true; // side finished
        }

    }

    template<short_t d,class T>
    void gsDPatchBase<d,T>::_handleBoundary(patchSide side)
    {
        std::vector<patchCorner> pcorners;
        std::vector<index_t> selectedIndices;
        gsBasis<T> * basis = &m_bases.basis(side.patch);
        sparseEntry_t entries;

        gsMatrix<index_t> indices = basis->boundaryOffset(side.side(),0);
        side.getContainedCorners(d,pcorners);
        for (index_t c =0; c!=2; c++)
        {
            selectedIndices.push_back(_indexFromVert(0,pcorners[c],side,0)); // index from vertex pcorners[c] along side psides[0] with offset 0.
            selectedIndices.push_back(_indexFromVert(1,pcorners[c],side,0)); // index from vertex pcorners[c] along side psides[0] with offset 0.
        }

        std::sort(selectedIndices.begin(),selectedIndices.end());
        std::vector<index_t> allIndices(indices.data(), indices.data() + indices.rows() * indices.cols());
        std::vector<index_t> result(allIndices.size());
        std::vector<index_t>::iterator it=std::set_difference (allIndices.begin(), allIndices.end(), selectedIndices.begin(), selectedIndices.end(), result.begin());
        result.resize(it-result.begin());

        index_t rowIdx,colIdx;
        T weight;
        for (std::vector<index_t>::iterator it = result.begin(); it!=result.end(); ++it)
        {
            rowIdx = m_mapModified.index(*it,side.patch);
            colIdx = m_mapOriginal.index(*it,side.patch);
            weight = 1.0;
            entries.push_back(std::make_tuple(rowIdx,colIdx,weight));

            // m_matrix(rowIdx,colIdx) = 1.0;
            // m_basisCheck[rowIdx] = true;
        }

        #pragma omp critical (handle_boundary)
        {
            _pushAndCheck(entries);

            m_sideCheck.at( _sideIndex(side.patch,side.side()) ) = true;
        }
    }

    template<short_t d,class T>
    void gsDPatchBase<d,T>::_handleInterior()
    {
        #pragma omp critical (handle_interior)
        {
        index_t rowIdx,colIdx;
        for(size_t p=0; p!=m_bases.nBases(); p++)
            for (index_t b=0; b!=m_bases.basis(p).size(); b++)
            {
                rowIdx = m_mapModified.index(b,p);
                // rowIdx = m_mapOriginal.index(b,p);
                if ( (!m_mapModified.is_free(b,p)) || (m_basisCheck[rowIdx]) )
                // if ( (m_basisCheck[rowIdx]) )
                    continue;
                colIdx = m_mapOriginal.index(b,p);
                m_matrix(rowIdx,colIdx) = 1;
                m_basisCheck[rowIdx] = true;
                // gsInfo<<"Basis function "<<rowIdx<<"(patch: "<<p<<"; fun: "<<b<<") is "<< (m_basisCheck[rowIdx] ? "" : "not ")<<"processed\n";
            }
        }
    }

    template<short_t d,class T>
    void gsDPatchBase<d,T>::_computeSmoothMatrix()
    {
        GISMO_ASSERT((size_t)m_mapModified.freeSize()==m_size,"Size does not match predicted size, m_mapModified.freeSize()="<<m_mapModified.freeSize()<<"; m_size="<<m_size);

        _resetChecks(true);

        // iterate over the vertices
// #pragma omp parallel
// {
        // #pragma omp parallel for collapse(2)
        for (size_t p=0; p<m_bases.nBases(); p++)
            for (index_t c=1; c<5; c++)
                _handleVertex(patchCorner(p,c));

        // #pragma omp parallel for
        for(gsBoxTopology::const_iiterator iit = m_topology.iBegin(); iit< m_topology.iEnd(); iit++)
            _handleInterface(*iit);

        // boundaries
        // #pragma omp parallel for
        for(gsBoxTopology::const_biterator bit = m_topology.bBegin(); bit< m_topology.bEnd(); bit++)
            _handleBoundary(*bit);

        _handleInterior();
// }
        if (m_options.getSwitch("Verbose")) { _whichHandled(); }

        _performChecks(true);
        m_matrix.makeCompressed();
    }

    template<short_t d,class T>
    void gsDPatchBase<d,T>::_computeMapper() // also initialize the mappers!
    {
        // interfaces
        _resetChecks(false);

        patchCorner pcorner;

        // For the interfaces, we eliminate all DoFs located on the interface, except the ones coinciding with the end vertices
// #pragma omp parallel
// {
//         #pragma omp parallel for
        for(gsBoxTopology::const_iiterator iit = m_topology.iBegin(); iit!= m_topology.iEnd(); iit++)
            _computeInterfaceMapper(*iit);
// }
        // On the boundaries, we don't do anything
// #pragma omp parallelz
// {
        // #pragma omp parallel for
        for(gsBoxTopology::const_biterator bit = m_topology.bBegin(); bit!= m_topology.bEnd(); bit++)
            _computeBoundaryMapper(*bit);
// }

// m_mapModified.finalize();
// m_mapOriginal.finalize();


        // For the vertices, we eliminate as follows (where v is the valence):
        // - No elimination when v==1
        // - One on each side when v==2
        // - All but three when v>2
        for (size_t p=0; p!=m_bases.nBases(); p++)
        {
            for (index_t c=1; c<5; c++)
            {
                pcorner = patchCorner(p,c);
                /// NOT PARALLEL YET
                _computeVertexMapper(pcorner);
            }
        }
        m_mapModified.finalize();
        m_mapOriginal.finalize();

        GISMO_ASSERT((size_t)m_mapModified.freeSize()==m_size,"Size does not match predicted size, m_mapModified.freeSize()="<<m_mapModified.freeSize()<<"; m_size="<<m_size);
        m_matrix.resize( m_size, m_mapOriginal.freeSize() );

        // gsDebugVar(m_mapModified.coupledSize());
        // gsDebugVar(m_mapModified.boundarySize());

    }

    // Same for DPatch and AlmostC1
    template<short_t d,class T>
    void gsDPatchBase<d,T>::_computeInterfaceMapper(boundaryInterface iface)
    {
        index_t sidx = _sideIndex( iface.second().patch,iface.second().side());
        if (m_sideCheck.at(sidx))
            return;

        gsBasis<T> * basis;
        std::pair<index_t,bool> vdata1, vdata2;
        std::vector<index_t> patches(2);
        std::vector<patchSide> psides(2);
        gsVector<index_t> indices;
        std::vector<patchCorner> pcorners;

        patches[0] = iface.first().patch;
        patches[1] = iface.second().patch;
        psides[0] = patchSide(iface.first().patch,iface.first().side()); // the interface on the first patch
        psides[1] = patchSide(iface.second().patch,iface.second().side()); // the interface on the second patch

        for (index_t p = 0; p != 2; p++)
        {
            sidx = _sideIndex( patches[p] ,psides[p] );
            if (m_sideCheck.at(sidx))
                continue;

            /*
                Eliminates the interior nodes on the interfaces

                o o o @ X | X @ o o o               X: eliminated DoFs by interface/boundary rule
                o o o @ X | X @ o o o               @: modified DoFs by interface rule
                o o o @ X | X @ o o o               o: preserved DoFs (interior)
                o o o x x | x x @ @ @               ?: Depends on the vertex (X and @ if not (interior vertex & valence = 4))
                o o o x x | x x X X X               x: handled in vertex rule
                ----------|----------
                          | x x X X X
                          | x x @ @ @
                          | o o o o o
                          | o o o o o
                          | o o o o o

            */

            basis = &m_bases.basis(patches[p]);
            indices = static_cast<gsVector<index_t>>( basis->boundary(psides[p]) );

            patchSide(patches[p],psides[p]).getContainedCorners(d,pcorners);
            vdata1 = this->_vertexData(pcorners[0]);
            vdata2 = this->_vertexData(pcorners[1]);

            // cast indices to an std::vector
            std::vector<index_t> allIndices(indices.data(), indices.data() + indices.rows() * indices.cols());
            // for(size_t i=0; i < allIndices.size(); i++)
            //     std::cout << allIndices.at(i) << ' ';

            std::vector<index_t> selectedIndices;
            // for both vertices of the side, add the indices at the vertex and one inside
            for (index_t c =0; c!=2; c++)
            {
                selectedIndices.push_back(_indexFromVert(0,pcorners[c],psides[p],0)); // index from vertex pcorners[c] along side psides[0] with offset 0.
                selectedIndices.push_back(_indexFromVert(1,pcorners[c],psides[p],0)); // index from vertex pcorners[c] along side psides[0] with offset 0.
            }

            std::sort(selectedIndices.begin(),selectedIndices.end());
            // for(size_t i=0; i < selectedIndices.size(); i++)
            //     std::cout << selectedIndices.at(i) << ' ';

            std::vector<index_t> result(allIndices.size());
            std::vector<index_t>::iterator it=std::set_difference (allIndices.begin(), allIndices.end(), selectedIndices.begin(), selectedIndices.end(), result.begin());
            result.resize(it-result.begin());

            gsAsMatrix<index_t> indices(result,result.size(),1);

            #pragma omp critical (side_interface)
            {
                m_mapModified.markBoundary(patches[p], indices);
                m_sideCheck.at(sidx) = true;
            }
        }
    }

    template<short_t d,class T>
    void gsDPatchBase<d,T>::_computeBoundaryMapper(patchSide boundary)
    {
        index_t sidx = _sideIndex(boundary.patch,boundary.side());
        m_sideCheck.at(sidx) = true;
    }

    template<short_t d,class T>
    void gsDPatchBase<d,T>::_computeVertexMapper(patchCorner pcorner)
    {
        index_t cidx = _vertIndex(pcorner.patch,pcorner.corner());
        if (m_vertCheck.at(cidx))
            return;

        bool C0 = m_C0s[cidx];
        bool interior;
        index_t valence;

        std::tie(valence,interior) = _vertexData(pcorner); // corner c
        if (!interior && valence==1) //valence = 1
            _computeMapperRegularCorner_v1(pcorner,valence);
        else if (!interior && valence==2 && C0)
            _computeMapperRegularBoundaryVertexNonSmooth_v2(pcorner,valence);
        else if (!interior && valence==2 && !C0)
            _computeMapperRegularBoundaryVertexSmooth_v2(pcorner,valence);
        else if (!interior && valence >2 && C0)
            _computeMapperIrregularBoundaryVertexNonSmooth_v(pcorner,valence);
        else if (!interior && valence >2 && !C0)
            _computeMapperIrregularBoundaryVertexSmooth_v(pcorner,valence);
        else if (interior)
            _computeMapperInteriorVertex_v(pcorner,valence);
        else
            GISMO_ERROR("Something went terribly wrong, interior="<<interior<<"; valence="<<valence);

        // label vertex as processed
        m_vertCheck[ cidx ] = true;
    }

    template<short_t d,class T>
    void gsDPatchBase<d,T>::_computeMapperRegularCorner_v1(patchCorner pcorner, index_t valence)
    {
        // do nothing,
    }

    // Same for DPatch and AlmostC1
    template<short_t d,class T>
    void gsDPatchBase<d,T>::_computeMapperRegularBoundaryVertexSmooth_v2(patchCorner pcorner, index_t valence)
    {
        std::vector<patchSide> psides(2);

        /*
            o o o @ X |i| X @ o o o                 x: eliminated DoFs by vertex rule
            o o o @ X |n| X @ o o o                 X: eliminated DoFs by interface/boundary rule
            o o o @ X |t| X @ o o o                 o: preserved DoFs (interior)
            o o o * x |e| x * o o o                 @: modified DoFs by interface rule
            o o o * x |r| x * o o o                 *: modified DoFs by vertex rule (unique DoFs)
            -----------------------                 %: modified DoFs by vertex rule (matched DoFs)
            -boundary-| | -boundary-
            -----------------------

            v = 4 (almost C1)
                o o o @ X |i| X @ o o o                 x: eliminated DoFs by vertex rule
                o o o @ X |n| X @ o o o                 X: eliminated DoFs by interface/boundary rule
                o o o @ X |t| X @ o o o                 o: preserved DoFs (interior)
                @ @ @ * x |e| x * @ @ @                 @: modified DoFs by interface rule
                X X X x x |r| x x X X X                 *: modified DoFs by vertex rule (unique DoFs)
                -----------------------
                interface-| |-interface
                -----------------------
                X X X x x |i| x x X X X
                @ @ @ * x |n| x * @ @ @
                o o o @ X |t| X @ o o o
                o o o @ X |e| X @ o o o
                o o o @ X |r| X @ o o o

        */

        // we mark the nodes belonging to the interface
        pcorner.getContainingSides(d,psides);
        for (size_t p=0; p!=psides.size(); p++)
        {
            // the 0,k (k=0,1) DoF should be eliminated
            if (m_topology.isInterface(psides[p]))
            {
                for (index_t k=0; k!=2; k++)
                    m_mapModified.eliminateDof(_indexFromVert(k,pcorner,psides[p],0),pcorner.patch);
            }
        }

        // // we mark the nodes belonging to the interface
        // pcorner.getContainingSides(d,psides);
        // for (size_t p=0; p!=psides.size(); p++)
        // {
        //     if (m_topology.isInterface(psides[p]))
        //     {
        //         // the 0,0 vertex should be eliminated
        //         m_mapModified.eliminateDof(basis->functionAtCorner(pcorner),pcorner.patch);
        //     }
        // }
    }

    // Different for DPatch and AlmostC1
    template<short_t d,class T>
    void gsDPatchBase<d,T>::_computeMapperRegularBoundaryVertexNonSmooth_v2(patchCorner pcorner, index_t valence)
    {
        std::vector<patchSide> psides(2);

        /*
            o o o @ X |i| X @ o o o                 x: eliminated DoFs by vertex rule
            o o o @ X |n| X @ o o o                 X: eliminated DoFs by interface/boundary rule
            o o o @ X |t| X @ o o o                 o: preserved DoFs (interior)
            o o o * x |e| x * o o o                 @: modified DoFs by interface rule
            o o o o o |r| o o o o o                 *: modified DoFs by vertex rule (unique DoFs)
            -----------------------                 %: modified DoFs by vertex rule (matched DoFs)
            -boundary-| | -boundary-
            -----------------------
        */
        // we mark the nodes belonging to the interface
        pcorner.getContainingSides(d,psides);
        for (size_t p=0; p!=psides.size(); p++)
        {
            if (m_topology.isInterface(psides[p]))
            {
                m_mapModified.eliminateDof(this->_indexFromVert(1,pcorner,psides[p],0),pcorner.patch);
            }
        }
    }

    template<short_t d,class T>
    void gsDPatchBase<d,T>::_computeMapperIrregularBoundaryVertexSmooth_v(patchCorner pcorner, index_t valence)
    {
        // std::vector<std::pair<index_t,index_t>> indices0  = _getAllInterfaceIndices(pcorner,0,m_bases);
        // std::vector<std::pair<index_t,index_t>> indices1  = _getAllInterfaceIndices(pcorner,1,m_bases);
        // std::vector<patchCorner> pcorners;
        // m_topology.getCornerList(pcorner,pcorners);
        // for (std::vector<patchCorner>::iterator it=pcorners.begin(); it!=pcorners.end(); it++)
        // {
        //     // mark the vertex as passed
        //     m_vertCheck[ _vertIndex(it->patch, it->corner()) ] = true;
        // }

        // // Eliminate the 1,0 and 0,1s
        // for (std::vector<std::pair<index_t,index_t>>::iterator it=indices1.begin(); it!=indices1.end(); it++)
        //     m_mapModified.eliminateDof(it->second,it->first);

        // _removeLowestIndices(indices0,3);
        // for (std::vector<std::pair<index_t,index_t>>::iterator it=indices0.begin(); it!=indices0.end(); it++)
        //     m_mapModified.eliminateDof(it->second,it->first);
        gsWarn<<"C0 handling for boundary corners with valence >2 has not yet been implemented. Using the default approach\n";
        this->_computeMapperIrregularBoundaryVertexNonSmooth_v(pcorner,valence);
    }

    template<short_t d,class T>
    void gsDPatchBase<d,T>::_computeMapperIrregularBoundaryVertexNonSmooth_v(patchCorner pcorner, index_t valence)
    {
        // Get all the 0,1 or 1,0 DoFs on the interfaces
        // The 0,0 DoFs are kept but deleted later from the matrix
        std::vector<patchSide> psides(2);
        pcorner.getContainingSides(d,psides);
        for (size_t p=0; p!=psides.size(); p++)
        {
            // the 0,k (k=0,1) DoF should be eliminated
            if (m_topology.isInterface(psides[p]))
                m_mapModified.eliminateDof(_indexFromVert(1,pcorner,psides[p],0),pcorner.patch);
        }
    }

    template<short_t d,class T>
    void gsDPatchBase<d,T>::_computeMapperInteriorVertex_v(patchCorner pcorner, index_t valence)
    {
        std::vector<patchSide> psides(2);
        /*
            o o o @ X |i| X @ o o o                 x: eliminated DoFs by vertex rule
            o o o @ X |n| X @ o o o                 X: eliminated DoFs by interface/boundary rule
            o o o @ X |t| X @ o o o                 o: preserved DoFs (interior)
            @ @ @ * x |e| x * @ @ @                 @: modified DoFs by interface rule
            X X X x x |r| x x X X X                 *: modified DoFs by vertex rule (unique DoFs)
            -----------------------
            interface-| |-interface
            -----------------------
            X X X x x |i| x x X X X
            @ @ @ * x |n| x * @ @ @
            o o o @ X |t| X @ o o o
            o o o @ X |e| X @ o o o
            o o o @ X |r| X @ o o o

        */
        // we mark the nodes belonging to the interfaces (both sides bordering the vertex)
        pcorner.getContainingSides(d,psides);
        for (size_t p=0; p!=psides.size(); p++)
        {
            m_mapModified.eliminateDof(this->_indexFromVert(0,pcorner,psides[p],0),pcorner.patch);
            m_mapModified.eliminateDof(this->_indexFromVert(1,pcorner,psides[p],0),pcorner.patch);
        }
    }

    template<short_t d,class T>
    void gsDPatchBase<d,T>::_handleRegularCorner(patchCorner pcorner)
    {
        std::vector<patchSide> psides;
        std::vector<index_t> indices(4);
        sparseEntry_t entries;

        pcorner.getContainingSides(d,psides);
        indices[0] = _indexFromVert(0,pcorner,psides[0],0); // b00
        indices[1] = _indexFromVert(1,pcorner,psides[0],0); // b01
        indices[2] = _indexFromVert(1,pcorner,psides[1],0); // b10
        indices[3] = _indexFromVert(1,pcorner,psides[1],1); // b11

        T weight = 1.0;
        index_t colIdx, rowIdx;
        for (std::vector<index_t>::iterator it = indices.begin(); it!=indices.end(); ++it)
        {
            rowIdx = m_mapModified.index(*it,pcorner.patch);
            colIdx = m_mapOriginal.index(*it,pcorner.patch);
            entries.push_back(std::make_tuple(rowIdx,colIdx,weight));
        }

        #pragma omp critical (handle_boundary_vertex_tt)
        {
            _pushAndCheck(entries);
            m_vertCheck[ _vertIndex(pcorner.patch, pcorner.corner()) ] = true;
        }
        // gsInfo<<"patch = "<<pcorner.patch<<", corner = "<<pcorner.corner()<<"\n";
        return;
    }

    template<short_t d,class T>
    void gsDPatchBase<d,T>::_handleRegularBoundaryVertexSmooth(patchCorner pcorner, index_t valence)
    {
        std::vector<patchSide> psides;
        std::vector<index_t> indices(3);

        sparseEntry_t entries;

        index_t colIdx, rowIdx;
        T weight;
        boundaryInterface iface;

        // Get the sides joining at the corner.
        pcorner.getContainingSides(d,psides);

        // 1. find the interface
        index_t iindex = m_topology.isInterface(psides[0]) ? 0 : 1;

        GISMO_ENSURE(m_topology.getInterface(psides[iindex],iface),"Must be an interface");

        // 2. collect indices
        // If we want C0 at this vertex, we only handle the row k=1.
        patchSide otherSide = iface.other(psides[iindex]);
        patchCorner otherCorner = iface.mapCorner(pcorner);
        for (index_t k = 0; k!=2; k++) // index of point over the interface
        {
            indices[0] = _indexFromVert(k,pcorner,psides[iindex],1); // bk1 on patch of iface
            indices[1] = _indexFromVert(k,pcorner,psides[iindex],0); // bk0 on patch of iface
            indices[2] = _indexFromVert(k,otherCorner,otherSide,0); // bk0 on other patch

            rowIdx = m_mapModified.index(indices[0],pcorner.patch);
            colIdx = m_mapOriginal.index(indices[0],pcorner.patch);
            weight = 1.0;
            entries.push_back(std::make_tuple(rowIdx,colIdx,weight));
            colIdx = m_mapOriginal.index(indices[1],psides[iindex].patch);
            weight = 0.5;
            entries.push_back(std::make_tuple(rowIdx,colIdx,weight));
            colIdx = m_mapOriginal.index(indices[2],otherSide.patch);
            weight = 0.5;
            entries.push_back(std::make_tuple(rowIdx,colIdx,weight));
        }

        #pragma omp critical (handle_boundary_vertex_tt)
        {
            _pushAndCheck(entries);
            m_vertCheck[ _vertIndex(pcorner.patch, pcorner.corner()) ] = true;
        }
    }

    template<short_t d,class T>
    void gsDPatchBase<d,T>::_handleRegularBoundaryVertexNonSmooth(patchCorner pcorner, index_t valence)
    {
        std::vector<patchSide> psides;
        std::vector<index_t> indices(3);

        sparseEntry_t entries;

        index_t colIdx, rowIdx;
        T weight;
        boundaryInterface iface;

        // Get the sides joining at the corner.
        pcorner.getContainingSides(d,psides);

        // 1. find the interface
        index_t iindex = m_topology.isInterface(psides[0]) ? 0 : 1;

        GISMO_ENSURE(m_topology.getInterface(psides[iindex],iface),"Must be an interface");

        // 2. collect indices
        // If we want C0 at this vertex, we only handle the row k=1.
        patchSide otherSide = iface.other(psides[iindex]);
        patchCorner otherCorner = iface.mapCorner(pcorner);
        for (index_t k = 1; k!=2; k++) // index of point over the interface
        {
            indices[0] = _indexFromVert(k,pcorner,psides[iindex],1); // bk1 on patch of iface
            indices[1] = _indexFromVert(k,pcorner,psides[iindex],0); // bk0 on patch of iface
            indices[2] = _indexFromVert(k,otherCorner,otherSide,0); // bk0 on other patch

            rowIdx = m_mapModified.index(indices[0],pcorner.patch);
            colIdx = m_mapOriginal.index(indices[0],pcorner.patch);
            weight = 1.0;
            entries.push_back(std::make_tuple(rowIdx,colIdx,weight));
            colIdx = m_mapOriginal.index(indices[1],psides[iindex].patch);
            weight = 0.5;
            entries.push_back(std::make_tuple(rowIdx,colIdx,weight));
            colIdx = m_mapOriginal.index(indices[2],otherSide.patch);
            weight = 0.5;
            entries.push_back(std::make_tuple(rowIdx,colIdx,weight));
        }

        #pragma omp critical (handle_boundary_vertex_tt)
        {
            _pushAndCheck(entries);
            m_vertCheck[ _vertIndex(pcorner.patch, pcorner.corner()) ] = true;
        }
    }

    template<short_t d,class T>
    void gsDPatchBase<d,T>::_handleIrregularBoundaryVertexSmooth(patchCorner pcorner, index_t valence)
    {
        gsWarn<<"C1 handling for boundary corners with valence >2 has not yet been implemented. Using the default approach\n";
        this->_handleIrregularBoundaryVertexNonSmooth(pcorner,valence);
    }

    template<short_t d,class T>
    void gsDPatchBase<d,T>::_handleIrregularBoundaryVertexNonSmooth(patchCorner pcorner, index_t valence)
    {
        std::vector<patchSide> psides;
        std::vector<patchCorner> corners;
        std::vector<index_t> indices;
        sparseEntry_t entries;

        boundaryInterface iface;
        patchSide otherSide;

        index_t colIdx, rowIdx;
        T weight;

        pcorner.getContainingSides(d,psides);

        std::vector<index_t> rowIndices, colIndices, patchIndices;

        // pcorner is the current corner
        m_topology.getCornerList(pcorner,corners);

        ////////////////////////////////////////////////////////////////////////////////
        // Influence of 1,1 to itself
        index_t b11_p1 = _indexFromVert(1,pcorner,psides[0],1); // point 1,1 (does not matter which reference side is taken)
        rowIdx = m_mapModified.index(b11_p1,pcorner.patch);
        colIdx = m_mapOriginal.index(b11_p1,pcorner.patch);

        weight = 1.;
        entries.push_back(std::make_tuple(rowIdx,colIdx,weight));

        ////////////////////////////////////////////////////////////////////////////////
        index_t idx;
        for (index_t k = 0; k!=2; k++)
        {
            // Check if one of the adjacent interface is a boundary;
            // if so, add weight 1.0 to itself
            if (!m_topology.getInterface(psides[k],iface)) // check if the side is NOT an interface
            {
                idx = _indexFromVert(1,pcorner,psides[k],0);
                rowIdx = m_mapModified.index(idx,pcorner.patch); //1,0 corner (on the boundary)
                colIdx = m_mapOriginal.index(idx,pcorner.patch); //1,0 corner (on the boundary)
                weight = 1.0;
                entries.push_back(std::make_tuple(rowIdx,colIdx,weight));
            }
            // else, add weight of 0.5 from the 0,1 or 1,0 vertices across the interface
            else
            {
                weight = 0.5;
                rowIdx = m_mapModified.index(b11_p1,pcorner.patch);

                patchSide otherSide = iface.other(psides[k]);
                patchCorner otherCorner = iface.mapCorner(pcorner);

                idx = _indexFromVert(1,pcorner,psides[k],0); // bk0 on patch of iface
                colIdx = m_mapOriginal.index(idx,pcorner.patch);
                entries.push_back(std::make_tuple(rowIdx,colIdx,weight));

                idx = _indexFromVert(1,otherCorner,otherSide,0); // bk0 on other patch
                colIdx = m_mapOriginal.index(idx,otherCorner.patch);
                entries.push_back(std::make_tuple(rowIdx,colIdx,weight));
            }
        }

        // Lastly, give the 0,0 a weight 1 to itself
        index_t b00_p1 = _indexFromVert(0,pcorner,psides[0],0); // point 0,0 (does not matter which reference side is taken)
        rowIdx = m_mapModified.index(b00_p1,pcorner.patch);
        colIdx = m_mapOriginal.index(b00_p1,pcorner.patch);

        weight = 1.;
        entries.push_back(std::make_tuple(rowIdx,colIdx,weight));

        #pragma omp critical (handle_boundary_vertex_ff)
        {
            _pushAndCheck(entries);
            m_vertCheck[ _vertIndex(pcorner.patch, pcorner.corner()) ] = true;
        }
    }

} // namespace gismo
