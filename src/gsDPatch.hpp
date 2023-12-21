/** @file gsDPatch.hpp

    @brief Creates the D-Patch smoothing matrix.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gsHSplines/gsHTensorBasis.h>
#include <gsHSplines/gsTHBSplineBasis.h>
#include <gsHSplines/gsTHBSpline.h>
#include <gsSolver/gsBlockOp.h>
#include <gsSolver/gsMatrixOp.h>

#include <gsIO/gsWriteParaview.h>

// #define PI 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647

namespace gismo
{
    // Constructors
    template<short_t d,class T>
    gsDPatch<d,T>::gsDPatch(const gsMultiPatch<T> & patches)
    :
    Base(patches)
    {
        this->defaultOptions();
    }

    // Constructors
    template<short_t d,class T>
    gsDPatch<d,T>::gsDPatch(const gsMultiBasis<T> & bases)
    :
    Base(bases)
    {
        this->defaultOptions();
    }

    template<short_t d,class T>
    gsDPatch<d,T>::~gsDPatch()
    {
        // freeAll(m_bases);
    }

    template<short_t d,class T>
    void gsDPatch<d,T>::defaultOptions()
    {
        Base::defaultOptions();
        m_options.addInt("Pi","Pi matrix to be applied, 0: Non-negative, 1: Idempotent",0);
        m_options.addInt("RefLevel","Refinement level",0);
        m_options.addReal("Beta","Beta parameter",0.4);
    }

    /*=====================================================================================
                                    Coefficients
    =====================================================================================*/
    template<short_t d,class T>
    gsMatrix<T> gsDPatch<d,T>::_preCoefficients(const gsMultiPatch<T> & patches)
    {
        GISMO_ASSERT(m_mapModified.isFinalized(),"Mapper is not finalized");

        gsMatrix<T> coefs(m_mapModified.freeSize(),patches.geoDim());

        index_t size;
        for (size_t p=0; p!=m_bases0.nBases(); p++) // patches
        {
            gsMatrix<T> tmpCoefs;
            if (m_tMatrices[p].rows()!=0 && m_tMatrices[p].cols()!=0)
                tmpCoefs = m_tMatrices[p]*patches.patch(p).coefs();
            else
                tmpCoefs = patches.patch(p).coefs();

            size = m_mapModified.patchSize(p);
            for (index_t k=0; k!=size; k++)
            {
                if (m_mapModified.is_free(k,p,0))
                    coefs.row(m_mapModified.index(k,p,0)) = tmpCoefs.row(k);
            }
        }

        // Correct the v=3 boundary vertices:
        std::fill(m_vertCheck.begin(), m_vertCheck.end(), false);
        for (size_t p=0; p!=m_bases0.nBases(); p++)
        {
            for (index_t c=1; c<5; c++)
            {
                index_t idx = this->_vertIndex(p,c);
                if(m_vertCheck[ idx] )
                    continue;

                patchCorner pcorner(p,c);
                std::pair<index_t,bool> vdata = this->_vertexData(pcorner); // corner c

                if (std::count(m_C0s.begin(), m_C0s.end(), pcorner))
                    continue;
                if (vdata.first==3 && !vdata.second)
                {
                    std::vector<patchCorner> otherCorners;
                    std::vector<patchSide> csides;

                    m_topology.getCornerList(pcorner,otherCorners);
                    index_t b00 = -1, b11i = -1;
                    std::vector<index_t> b11b;
                    for (std::vector<patchCorner>::iterator corner = otherCorners.begin(); corner != otherCorners.end(); corner++)
                    {
                        corner->getContainingSides(d,csides);
                        if ( m_topology.isBoundary(csides[0]) || m_topology.isBoundary(csides[1]) ) //
                            b11b.push_back(m_mapModified.index( this->_indexFromVert(m_bases0,1,*corner,csides[0],1) , corner->patch) );
                        else if (b11i==-1)
                            b11i = m_mapModified.index( this->_indexFromVert(m_bases0,1,*corner,csides[0],1) , corner->patch);
                        else
                            GISMO_ERROR("b11i is already assigned?");

                        if (corner==otherCorners.begin())
                        {
                            const gsBasis <T> * basis = &m_bases0.basis(corner->patch);
                            b00 = m_mapModified.index( basis->functionAtCorner(corner->corner()), corner->patch );
                        }

                        idx = this->_vertIndex(corner->patch,corner->corner());
                        m_vertCheck[ idx ] = true;
                    }
                    coefs.row(b00) = coefs.row(b11b[0]) + coefs.row(b11b[1]) - coefs.row(b11i);
                }
                else
                    m_vertCheck[ idx ] = true;

            }
        }

        return coefs;
    }

    /*=====================================================================================
                                    Construction functions
    =====================================================================================*/

    template<short_t d,class T>
    void gsDPatch<d,T>::_initTHB()
    {
        // Cast all patches of the mp object to THB splines
        for (size_t k=0; k!=m_Bbases.nBases(); ++k)
        {
            // Check the basis and make the new level(s)
            std::vector<gsKnotVector<T>> KVs(d);
            index_t degree;
            if ( const gsTensorBSplineBasis<d,T> * tbasis0 = dynamic_cast<const gsTensorBSplineBasis<d,T> * > (&m_Bbases.basis(k)) )
            {
                gsTHBSplineBasis<d,T> thbBasis(*tbasis0,true);
                for (short_t dim=0; dim!=d; dim++)
                {
                    GISMO_ENSURE(tbasis0->degree(dim)>=2,"Degree of the basis must be larger than or equal to 2, but is "<<tbasis0->degree(dim)<<" (component "<<d<<")");
                    const gsKnotVector<T> & KV = tbasis0->knots(dim);
                    KVs[dim] = tbasis0->knots(dim);
                    degree = KVs[dim].degree();

                    // Every knot needs multiplicity p-1. For degree 2, this gives the same basis!
                    for (typename gsKnotVector<T>::uiterator knot = std::next(KV.ubegin()); knot!=std::prev(KV.uend()); knot++)
                        KVs[dim].insert(*knot,degree-KV.multiplicity(*knot)-1);

                }
                // Create level 1
                gsTensorBSplineBasis<d,T> tbasis1(KVs);
                thbBasis.addLevel(tbasis1);

                // Create level 2
                gsTensorBSplineBasis<d,T> tbasis2 = tbasis1;
                for (short_t dim = 0; dim!=d; dim++)
                    tbasis2.uniformRefine(1,tbasis2.degree(dim)-1,dim);
                thbBasis.addLevel(tbasis2);

                m_bases0.addBasis(thbBasis.clone());
            }
            else
                GISMO_ERROR("Basis can only be constructed on gsTensorBSplineBasis");
        }
        m_bases0.setTopology(m_topology);
        m_bases = m_bases0;
    }

    template <short_t d, class T>
    void gsDPatch<d,T>::_refBoxes(std::vector<std::vector<index_t>> & patchBoxes)
    {
        patchBoxes.clear();
        patchBoxes.resize(m_bases0.size());

        // prepare the geometry
        gsMatrix<index_t> box(d,2);
        std::vector<index_t> boxes;
        gsVector<bool> pars;
        index_t nelements;
        index_t degree;
        patchCorner corner;
        std::vector<patchCorner> cornerList;
        std::vector<std::vector<patchCorner> > cornerLists;

        m_topology.getEVs(cornerLists);

        index_t N = cornerLists.size();

        // Make a mask of corners per patch to track which ones have been handled
        gsMatrix<bool> mask(m_bases0.nBases(),math::pow(2,d));
        mask.setConstant(false);

        for (index_t v =0; v<N; v++)
        {// Loop over EVs
            for (size_t c = 0; c<cornerLists[v].size(); c++)
            {// Loop over corners per EV
                corner = cornerLists[v].at(c);
                gsHTensorBasis<d,T> * basis = dynamic_cast<gsHTensorBasis<d,T>*>(&m_bases0.basis(corner.patch));

                if (mask(corner.patch,corner.corner()-1))
                    continue;

                corner.parameters_into(d,pars);
                box.setZero();
                for (short_t dim = 0; dim!=d; dim++)
                {
                    const gsKnotVector<T> & KV = basis->tensorLevel(0).knots(dim);
                    degree = KV.degree();
                    nelements = (degree < 4) ? 2 : 1;
                    nelements *= std::pow(2,m_options.getInt("RefLevel"));

                    GISMO_ASSERT(nelements<=(index_t)KV.numElements(),"Need more elements than available for refinement around corner "<<corner.corner()<<" of patch "<<corner.patch<<".\n"<<"nelements = "<<nelements<<"; KV.numElements() = "<<KV.numElements()<<"\n");

                    box.row(dim).setConstant(pars(dim)*(KV.uSize()-1));
                    box(dim,!pars(dim)) += ( 1 - 2*pars(dim) ) * nelements; // subtracts from box(d,0) if pars(d)==1, adds to box(d,1) if pars(d)==0

                    // If all elements in this direction are refined, we need to add the other corner of this side to the list of corners to be refined
                    if ((index_t)KV.numElements()==nelements)
                    {
                        // Get the patch side in the direction of the knot vector KV
                        GISMO_ASSERT(d==2,"This does not work for d!=2!");
                        boxCorner otherCorner;
                        if      ((corner.m_index==1 && dim==0) || (corner.m_index==4 && dim==1)) // the other corner is south east
                            otherCorner = 2;
                        else if ((corner.m_index==1 && dim==1) || (corner.m_index==4 && dim==0)) // the other corner is north west
                            otherCorner = 3;
                        else if ((corner.m_index==2 && dim==0) || (corner.m_index==3 && dim==1)) // the other corner is south west
                            otherCorner = 1;
                        else if ((corner.m_index==2 && dim==1) || (corner.m_index==3 && dim==0)) // the other corner is north east
                            otherCorner = 4;
                        else
                            GISMO_ERROR("Combination unknown...");

                        patchCorner otherPCorner(corner.patch,otherCorner.m_index);

                        cornerList.clear();
                        m_topology.getCornerList(otherPCorner,cornerList);
                        cornerLists.push_back(cornerList);
                        N++;
                    }
                }
                boxes.clear();
                // Assign boxes. This is the box on the current level, level 0.
                boxes.push_back(0);
                boxes.insert(boxes.end(), box.data(), box.data()+box.size());
                patchBoxes.at(corner.patch).insert(patchBoxes.at(corner.patch).end(), boxes.begin(), boxes.end());

                mask(corner.patch,corner.corner()-1) = true;
            }// Loop over corners per EV
        }// Loop over EVs
    }

    template<short_t d,class T>
    void gsDPatch<d,T>::_initBasis()
    {
        std::vector< std::vector<index_t> > elVec;
        this->_refBoxes(elVec);
        m_tMatrices.resize(m_bases0.nBases());
        for (size_t p=0; p!=m_bases0.nBases(); p++)
        {
            // Transform using gsAsMatrix
            gsAsMatrix<index_t> boxMat(elVec[p],2*d+1,elVec[p].size()/(2*d+1));
            boxMat.row(0).array() += 1;
            boxMat.block(1,0,boxMat.rows()-1,boxMat.cols()).array() *= 2;

            // Refine elements
            gsHTensorBasis<d,T> *basis = dynamic_cast<gsHTensorBasis<d,T>*>(&m_bases0.basis(p));
            std::vector< gsSortedVector< index_t > > xmat = basis->getXmatrix();
            basis->refineElements_withTransfer(elVec[p],m_tMatrices[p]);
        }
        m_bases = m_bases0;
    }

    template<short_t d,class T>
    void gsDPatch<d,T>::_makeTHB() //IMPLEMENT THIS
    {
        gsMultiBasis<> refBases = m_bases;
        std::vector< std::vector<index_t> > elVec;
        this->_refBoxes(elVec);

        gsSparseMatrix<T> tmp;
        index_t rows = 0, cols = 0;
        std::vector<gsEigen::Triplet<T,index_t>> tripletList;
        for (size_t p=0; p!=m_bases0.nBases(); p++)
        {
            // Transform using gsAsMatrix
            gsAsMatrix<index_t> boxMat(elVec[p],2*d+1,elVec[p].size()/(2*d+1));
            boxMat.row(0).array() += 2;
            boxMat.block(1,0,boxMat.rows()-1,boxMat.cols()).array() *= 4;

            gsHTensorBasis<d,T> *basis = dynamic_cast<gsHTensorBasis<d,T>*>(&refBases.basis(p));
            std::vector< gsSortedVector< index_t > > xmat = basis->getXmatrix();
            basis->refineElements_withTransfer(elVec[p],tmp);
            for (index_t i = 0; i<tmp.outerSize(); ++i)
                for (typename gsSparseMatrix<T>::iterator it(tmp,i); it; ++it)
                    tripletList.push_back(gsEigen::Triplet<T,index_t>(it.row()+rows,it.col()+cols,it.value()));

            rows += tmp.rows();
            cols += tmp.cols();
        }

        m_tMatrix.resize(rows,cols);
        m_tMatrix.setFromTriplets(tripletList.begin(), tripletList.end());

        m_tMatrix.makeCompressed();
        m_bases = refBases;

        // redefine the mappers
        // m_mapModified = gsDofMapper(m_bases);
        m_mapOriginal = gsDofMapper(m_bases);
        m_mapOriginal.finalize();
    }

    template<short_t d,class T>
    void gsDPatch<d,T>::_computeEVs()
    {
        /*
            Our goal is to create three vectors cij[0], cij[1], cij[2] which all contain the
            cij[0], cij[1] and cij[2] coefficients of the patches around the EV in the right order
            (counter)-clockwise.
        */

        std::vector<std::vector<patchCorner> > cornerLists;
        m_topology.getEVs(cornerLists);

        if (cornerLists.size()!=0)
        {
            // gsWriteCSV(m_matrix.toDense(),"matrix0.csv");
            m_matrix = m_matrix * m_tMatrix.transpose();
            // gsWriteCSV(m_tMatrix.toDense(),"matrixTHB.csv");
            // gsWriteCSV(m_matrix.toDense(),"matrix1.csv");

            gsSparseMatrix<T> pi;
            std::vector<patchSide> sides(2);
            std::vector<patchSide> allSides;
            std::vector<std::vector<patchCorner> > cornerLists;
            std::vector<patchCorner> corners;
            std::vector<index_t> allPatches;
            std::map<index_t,index_t> patches;
            std::vector<boundaryInterface> interfaces;
            m_topology.getEVs(cornerLists);

            sparseEntry_t entries;

            if (cornerLists.size()!=0)
            {
                for (size_t v =0; v!=cornerLists.size(); v++) // over EVs
                {
                    patches.clear();
                    index_t N = cornerLists[v].size();

                    allPatches.resize(m_bases.nBases());
                    corners.resize(N);
                    interfaces.resize(N);

                    std::vector<gsMatrix<index_t>> cij(3);
                    for (index_t k=0; k!=3; k++)
                        cij.at(k) = gsMatrix<index_t>::Zero(N,1);

                    std::vector<gsMatrix<index_t>> cijo = cij;

                    /*
                        First, we loop over all the interfaces to construct a (counter)clock-wise map for our coefficients
                        Looping clock-wise or counter clock-wise does not matter, as long as the coefficients are neighboring
                    */
                    // Loop over all sides such that we can fill cij[1] and cij[2]. We just start from side 0 of corner 0
                    patchCorner corner = cornerLists[v][0];
                    patchCorner otherCorner = patchCorner(0,0);
                    corner.getContainingSides(d,sides);
                    patchSide side = sides[0];
                    patchSide otherSide;
                    std::vector<patchCorner> pcorners(2);

                    // Initialize patch map
                    for (index_t i = 0; i!=N; i++) // over interfaces
                    {
                        patches.insert(std::make_pair(side.patch,i));
                        corners[i] = corner;
                        bool isInterface = m_topology.getInterface(side,interfaces[i]);
                        GISMO_ENSURE(isInterface,"Side must be an interface!");

                        std::vector<boxCorner> adjcorners;
                        m_topology.getNeighbour(side,otherSide);
                        otherSide.getContainedCorners(d,adjcorners);
                        for (index_t k=0; k!=N; k++)
                        {
                            if (cornerLists[v][k] == patchCorner(otherSide.patch,adjcorners[0]))
                                otherCorner = patchCorner(otherSide.patch,adjcorners[0]);
                            else if (cornerLists[v][k] == patchCorner(otherSide.patch,adjcorners[1]))
                                otherCorner = patchCorner(otherSide.patch,adjcorners[1]);
                            else continue;
                        }
                        GISMO_ENSURE(otherCorner!=patchCorner(0,0),"Error");

                        // interfaces[i] = boundaryInterface(side,otherSide,d);
                        // GISMO_ASSERT(corners[i].patch==interfaces[i].first().patch,"Must be true");

                        // get the NEXT side
                        otherCorner.getContainingSides(d,sides);
                        if (otherSide == sides[0])
                            otherSide = sides[1];
                        else if (otherSide == sides[1])
                            otherSide = sides[0];
                        else
                            GISMO_ERROR("An error occurred.");

                        corner = otherCorner;
                        side = otherSide;
                    }

                    for (index_t i = 0; i!=N; i++) // over corners in EVs
                    {

                        otherCorner = interfaces[i].mapCorner(corners[i]);

                        if (interfaces[i].first().patch==corners[i].patch)
                        {
                            otherSide = interfaces[i].second();
                            side = interfaces[i].first();
                        }
                        else
                        {
                            otherSide = interfaces[i].first();
                            side = interfaces[i].second();
                        }

                        // C11 coefficients
                        cij[0](patches[side.patch],0) = this->_indexFromVert(m_bases,1,corners[i],side,1);
                        // C21 coefficients
                        cij[1](patches[otherSide.patch],0) = this->_indexFromVert(m_bases,2,otherCorner,otherSide,1);
                        // C12 coefficients
                        cij[2](patches[side.patch],0) = this->_indexFromVert(m_bases,2,corners[i],side,1);

                        // C11 coefficients
                        cijo[0](patches[side.patch],0) = this->_indexFromVert(m_bases0,1,corners[i],side,1);
                        // C21 coefficients
                        cijo[1](patches[otherSide.patch],0) = this->_indexFromVert(m_bases0,2,otherCorner,otherSide,1);
                        // C12 coefficients
                        cijo[2](patches[side.patch],0) = this->_indexFromVert(m_bases0,2,corners[i],side,1);
                    }

                    std::vector<gsMatrix<index_t>> rowIndices(3);
                    for (index_t k=0; k!=3; k++)
                        rowIndices.at(k) = gsMatrix<index_t>::Zero(N,1);

                    for (index_t i = 0; i!=N; i++)
                    {
                        corner = corners[i];
                        corner.getContainingSides(d,sides);
                        // we look for the 1,1 index so it does not matter which side we use
                        for (index_t k=0; k!=3; k++)
                        {
                            GISMO_ASSERT(m_mapModified.is_free(cijo[k](i,0),corner.patch),"Something went wrong in the indexing of the sparse matrix for EVs.\n corner = "<<corner.corner()<<"\n patch = "<<corner.patch<<"\n k = "<<k<<"\n i = "<<i<<"\n cijo[k](i,0) = "<<cijo[k](i,0)<<"\n cijo[k] = "<<cijo[k]<<"\n");
                            rowIndices[k](i,0) = m_mapModified.index(cijo[k](i,0),corner.patch);
                            GISMO_ASSERT(m_mapOriginal.is_free(cij[k](i,0),corner.patch),"Something went wrong in the indexing of the sparse matrix for EVs");
                            cij[k](i,0) = m_mapOriginal.index(cij[k](i,0),corner.patch);
                        }
                    }

                    gsMatrix<T> Pi = _makePi(N);
                    gsVector<T> c(3*N);
                    index_t colIdx = 0, idx = 0;
                    for (index_t i = 0; i!=N; i++) // for all involved corners
                    {
                        for (index_t k = 0; k!=3; k++) // loop over ij=11,12,21
                        {
                            for (index_t j=0; j!=N; j++) // loop over the connected corners
                                for (index_t l = 0; l!=3; l++) // loop over ij=11,12,21
                                    c.at(j+l*N) = m_matrix.coeff(rowIndices[k](i,0),cij[l](j,0));
                            c = Pi * c;
                            for (index_t j=0; j!=N; j++) // loop over the connected corners
                                for (index_t l = 0; l!=3; l++) // loop over ij=11,12,21
                                    entries.push_back(std::make_tuple(rowIndices[k](i,0),cij[l](j,0),c.at(j+l*N)));
                        }
                    }

                    #pragma omp critical (_computeEV1)
                        _push(entries);

                    entries.clear();

                    // smoothing center, i.e all 0,0, 1,0 and 0,1 basis functions get the same value as the 1,1
                    for (index_t i = 0; i!=N; i++) // for all involved corners
                    {
                        for (index_t j=0; j!=N; j++) // loop over the connected corners
                        {
                            for (index_t k = 0; k!=3; k++) // loop over ij=11,12,21
                            {
                                corners[j].getContainingSides(d,sides);
                                /////////////////////////////////////////////////////////////////
                                colIdx = this->_indexFromVert(m_bases,0,corners[j],sides[0],0); // 0,0
                                GISMO_ASSERT(m_mapOriginal.is_free(colIdx,corners[j].patch),"Something went wrong in the indexing of the sparse matrix for EVs");
                                colIdx = m_mapOriginal.index(colIdx,corners[j].patch);
                                entries.push_back(std::make_tuple(rowIndices[k](i,0),colIdx,m_matrix.coeff(rowIndices[k](i,0),cij[0](j,0))));

                                colIdx = this->_indexFromVert(m_bases,1,corners[j],sides[0],0); // 1,0
                                GISMO_ASSERT(m_mapOriginal.is_free(colIdx,corners[j].patch),"Something went wrong in the indexing of the sparse matrix for EVs");
                                colIdx = m_mapOriginal.index(colIdx,corners[j].patch);
                                entries.push_back(std::make_tuple(rowIndices[k](i,0),colIdx,m_matrix.coeff(rowIndices[k](i,0),cij[0](j,0))));

                                colIdx = this->_indexFromVert(m_bases,1,corners[j],sides[1],0); // 0,1
                                GISMO_ASSERT(m_mapOriginal.is_free(colIdx,corners[j].patch),"Something went wrong in the indexing of the sparse matrix for EVs");
                                colIdx = m_mapOriginal.index(colIdx,corners[j].patch);
                                entries.push_back(std::make_tuple(rowIndices[k](i,0),colIdx,m_matrix.coeff(rowIndices[k](i,0),cij[0](j,0))));

                            }
                        }
                    }

                    #pragma omp critical (_computeEV2)
                        _push(entries);

                    entries.clear();

                    // interface smoothing
                    for (index_t i = 0; i!=N; i++) // for all involved interfaces
                    {
                        patchCorner corner = corners[patches[interfaces[i][0].patch]];
                        patchCorner otherCorner = corners[patches[interfaces[i][1].patch]];
                        patchSide side = interfaces[i][0];
                        patchSide otherSide = interfaces[i][1];
                        for (index_t k = 2; k!=4 ; k++)// std::max(basis1->maxDegree(),basis2->maxDegree())+
                        {
                            idx = this->_indexFromVert(m_bases,k,corner,side,0);
                            GISMO_ASSERT(m_mapOriginal.is_free(idx,side.patch),"Something went wrong in the indexing of the sparse matrix for EVs");
                            index_t j0k = m_mapOriginal.index(idx,side.patch);

                            idx = this->_indexFromVert(m_bases,k,otherCorner,otherSide,0);
                            GISMO_ASSERT(m_mapOriginal.is_free(idx,otherSide.patch),"Something went wrong in the indexing of the sparse matrix for EVs");
                            index_t jk0 = m_mapOriginal.index(idx,otherSide.patch);

                            idx = this->_indexFromVert(m_bases,k,corner,side,1);         // point (k,0)
                            GISMO_ASSERT(m_mapOriginal.is_free(idx,side.patch),"Something went wrong in the indexing of the sparse matrix for EVs");
                            index_t jk1 = m_mapOriginal.index(idx,side.patch); // point (k,0)

                            idx = this->_indexFromVert(m_bases,k,otherCorner,otherSide,1);         // point (k,0)
                            GISMO_ASSERT(m_mapOriginal.is_free(idx,otherSide.patch),"Something went wrong in the indexing of the sparse matrix for EVs");
                            index_t j1k = m_mapOriginal.index(idx,otherSide.patch); // point (k,0)

                            index_t row;
                            for (index_t l = 0; l!=3; l++) // loop over ij=11,12,21
                            {
                                for (index_t r=0; r!=N; r++)
                                {
                                    row = rowIndices[l](r,0);
                                    entries.push_back(std::make_tuple(rowIndices[l](r,0),j0k,0.5 * ( m_matrix.coeff(row,jk1) + m_matrix.coeff(row,j1k) )));
                                    entries.push_back(std::make_tuple(rowIndices[l](r,0),jk0,0.5 * ( m_matrix.coeff(row,jk1) + m_matrix.coeff(row,j1k) )));
                                }
                            }
                        }
                    }

                    #pragma omp critical (_computeEV3)
                        _push(entries);
                }
            }
        }

        m_matrix.makeCompressed();
        // gsDebugVar(m_matrix.toDense());
    }

    template<short_t d,class T>
    gsMatrix<T> gsDPatch<d,T>::_makePi(index_t valence)
    {
        gsMatrix<T> Pi(3*valence,3*valence);
        gsMatrix<T> P(valence,9);
        P.setZero();

        T pi = 4*std::atan(1);
        T phi = 2*pi / valence;
        std::complex<T> I(0,1);
        T beta = m_options.getReal("Beta") * std::pow(0.5,m_options.getInt("RefLevel"));
        T psi = std::arg( std::complex<T>((T(1.0)+I*beta*T(math::sin(phi)))*math::exp( -I*phi / T(2.0) ) ));
        if (m_options.getInt("Pi")==0)//idempotent
        {
            for (index_t j=0; j!=valence; j++)
            {
                P(j,0) = P(j,1) = P(j,2) = P(j,3) = P(j,6) = 1.0 / (3.0 * valence);;
                P(j,4) = P(j,8) = 1.0 / (3.0 * valence) * ( 1.0 + 3.0*math::cos( j * phi ) );
                P(j,5) = 1.0 / (3.0 * valence) * ( 1.0 + 3.0*math::cos( 2.0 * psi + j * phi ) );
                P(j,7) = 1.0 / (3.0 * valence) * ( 1.0 + 3.0*math::cos( 2.0 * psi - j * phi ) );
            }
        }
        else if (m_options.getInt("Pi")==1) //non-negative entries
        {
            for (index_t j=0; j!=valence; j++)
            {
                P(j,0) = P(j,3) = P(j,6) = 0;
                P(j,1) = P(j,2) = 1.0 / (2.0 * valence);
                P(j,4) = P(j,8) = 1.0 / (2.0 * valence) * ( 1.0 + math::cos( j * phi ) );
                P(j,5) = 1.0 / (2.0 * valence) * ( 1.0 + math::cos( 2.0 * psi + j * phi ) );
                P(j,7) = 1.0 / (2.0 * valence) * ( 1.0 + math::cos( 2.0 * psi - j * phi ) );
            }
        }
        else
            GISMO_ERROR("Pi option unknown");

        index_t offsetI, offsetJ = 0;
        gsMatrix<T> tmp(valence,valence);
        for (index_t i=0; i!=9; i++)
        {
            offsetI = (i / 3)*valence;// std::floor(i/3)
            offsetJ = (i % 3)*valence;
            for (index_t j=0; j!=valence; j++ )
                for (index_t k=0; k!=valence; k++ )
                {
                    index_t c = (j-k) % valence;
                    if (c < 0)
                        c += valence;
                    tmp(j,k) = P( c, i);
                }

            // gsDebugVar(offsetI);
            // gsDebugVar(offsetJ);
            // tmp.setOnes();
            // tmp *= i;
            Pi.block(offsetI,offsetJ,valence, valence) = tmp;
        }

        return Pi;
    }

    template<short_t d,class T>
    void gsDPatch<d,T>::_countDoFs() // also initialize the mappers!
    {
        size_t tmp;
        m_size = tmp = 0;
        // number of interior basis functions
        for (size_t p=0; p!=m_bases.nBases(); p++)
        {
            tmp += m_bases.basis(p).size();
            for (index_t k=0; k!=2; k++)
            {
                tmp -= m_bases.basis(p).boundaryOffset(boxSide(1),k).size();
                tmp -= m_bases.basis(p).boundaryOffset(boxSide(2),k).size();
                tmp -= m_bases.basis(p).boundaryOffset(boxSide(3),k).size()-4;
                tmp -= m_bases.basis(p).boundaryOffset(boxSide(4),k).size()-4;
            }
        }

        // gsDebug<<"Number of interior DoFs: "<<tmp<<"\n";
        m_size += tmp;

        // interfaces
        gsBasis<T> * basis1;
        gsBasis<T> * basis2;
        gsVector<index_t> indices1,indices2;
        tmp = 0;
        for(gsBoxTopology::const_iiterator iit = m_topology.iBegin(); iit!= m_topology.iEnd(); iit++)
        {
            basis1 = &m_bases.basis(iit->first().patch);
            basis2 = &m_bases.basis(iit->second().patch);
            tmp += basis1->boundary(iit->first().side()).size() - 4;
            tmp += basis2->boundary(iit->second().side()).size() - 4;
        }
        // gsDebug<<"Number of interface DoFs: "<<tmp<<"\n";
        m_size += tmp;

        // boundaries
        tmp = 0;
        for(gsBoxTopology::const_biterator bit = m_topology.bBegin(); bit!= m_topology.bEnd(); bit++)
        {
            basis1 = &m_bases.basis(bit->patch);
            tmp += (basis1->boundaryOffset(bit->side(),0).size() - 4);
            tmp += (basis1->boundaryOffset(bit->side(),1).size() - 4);
        }
        // gsDebug<<"Number of boundary DoFs: "<<tmp<<"\n";
        m_size += tmp;

        // vertices
        tmp = 0;
        std::vector<bool> passed(m_bases.nBases()*4);
        std::fill(passed.begin(), passed.end(), false);

        std::vector<patchCorner> corners;
        // index_t corn = 0;
        for (size_t p=0; p!=m_bases.nBases(); p++)
            for (index_t c=1; c<5; c++)
            {
                index_t idx = this->_vertIndex(p,c);
                if (!passed.at(idx))
                {
                    m_topology.getCornerList(patchCorner(p,c),corners);

                    for (size_t k=0; k!=corners.size(); k++)
                        passed.at(this->_vertIndex(corners[k].patch,corners[k])) = true;

                    std::pair<index_t,bool> vdata = this->_vertexData(patchCorner(p,c)); // corner c
                    bool C0 = m_C0s[idx];
                    if ((!vdata.second) && vdata.first==1) // valence = 1, must be a boundary vertex
                        tmp += 4;
                    else if ((!vdata.second) && vdata.first==2 && !C0)
                        tmp += 4;
                    else if ((!vdata.second) && vdata.first>2 && !C0)
                        tmp += 4;
                    else if ((!vdata.second) && vdata.first==2 && C0)
                        tmp += 6;
                    else if ((!vdata.second) && vdata.first>2 && C0)
                        tmp += 4;
                    else
                        tmp += vdata.first; // valence;

                    // corn +=1;
                }
            }
        // gsDebug<<"Number of unique corners: "<<corn<<"\n";

        // gsDebug<<"Number of vertex DoFs: "<<tmp<<"\n";

        m_size += tmp;
    }

    template<short_t d,class T>
    void gsDPatch<d,T>::_handleIrregularBoundaryVertexSmooth(patchCorner pcorner, index_t valence)
    {
        // // 2. make container for the interfaces
        // std::vector<boundaryInterface> ifaces;
        // boundaryInterface iface;
        // std::vector<patchSide> boundaries;
        // index_t extraRow;
        // std::vector<patchCorner> temp_corners, corners;
        // std::vector<patchSide> psides;
        // std::vector<index_t> colIndices,rowIndices, patches;

        // sparseEntry_t entries;

        // m_topology.getCornerList(pcorner,corners);

        // /*
        //     Warning:
        //     This case handles all patchCorners related to this vertex at once!
        //     */

        // // 1. get all adjacent vertices

        // // find corner with the lowest patch index
        // patchCorner lowest = corners[0];
        // for (size_t k=0; k!=corners.size(); k++)
        //     lowest = (corners[k].patch < lowest.patch ) ? corners[k] : lowest;

        // lowest.getContainingSides(d,psides);
        // // get (0,0) index from the corner with lowest patch number and store the row index.
        // extraRow = this->_indexFromVert(0,lowest,psides[0],0);
        // extraRow = m_mapModified.index(extraRow,lowest.patch);

        // // 2. loop over the adjacent vertices
        // colIndices.clear();
        // rowIndices.clear();
        // patches.clear();
        // for (std::vector<patchCorner>::iterator it = corners.begin(); it!=corners.end(); ++it)
        // {
        //     // *it is a patchCorner
        //     // 3. determine if one of the contained sides is a boundary
        //     it->getContainingSides(d,psides);

        //     // store the 0,0 indices
        //     colIndices.push_back( this->_indexFromVert(0,*it,psides[0],0) );
        //     rowIndices.push_back( this->_indexFromVert(1,*it,psides[0],1) );
        //     patches.push_back(it->patch);

        //     for (index_t k = 0; k!=2; k++) // index of point over the interface
        //     {
        //         if ( m_topology.getInterface(psides[k],iface) ) // if it is an interface, store it
        //         {
        //             ifaces.push_back(iface);
        //         }
        //         else                                            // if not, then store the side
        //         {
        //             //
        //             index_t colIdx = m_mapOriginal.index(this->_indexFromVert(1,*it,psides[k],0),it->patch);
        //             index_t rowIdx = m_mapModified.index(this->_indexFromVert(1,*it,psides[k],1),it->patch);
        //             entries.push_back(std::make_tuple(rowIdx,colIdx,0.5));
        //             entries.push_back(std::make_tuple(extraRow,colIdx,0.5));
        //             // m_matrix(rowIdx,colIdx) = 0.5;
        //             // m_matrix(extraRow,colIdx) = 0.5;
        //             // m_basisCheck[rowIdx] = true;
        //         }
        //     }
        // }

        // // GISMO_ASSERT(boundaries.size()==2,"There must be two boundaries that are not an interface!");
        // // ifaces.push_back(boundaryInterface(boundaries[0],boundaries[1],d));

        // // the extra (0,0) node gets 0.25 for all (0,0) entries
        // for (size_t k = 0; k!=colIndices.size(); k++)
        // {
        //     index_t colIdx = m_mapOriginal.index(colIndices[k],patches[k]);
        //     entries.push_back(std::make_tuple(extraRow,colIdx,0.25));
        //     // m_matrix(extraRow,colIdx) = 0.25;
        // }

        // // Fill the matrix entries related to the 1,1 coefs (stored in indices) with 1 for itself and with 0.25 for the others
        // for (size_t k = 0; k!=rowIndices.size(); k++)
        // {
        //     index_t rowIdx = m_mapModified.index(rowIndices[k],patches[k]);
        //     index_t colIdx = m_mapOriginal.index(rowIndices[k],patches[k]);
        //     // Fill the matrix entries related to itself (1,1) with a 1.0
        //     entries.push_back(std::make_tuple(rowIdx,colIdx,1.0));
        //     // m_matrix(rowIdx,colIdx) = 1.0;

        //     for (size_t l = 0; l!=colIndices.size(); l++)
        //     {
        //         // Fill the matrix entries related to the 0,0 coefs (stored in indices) 0.25 for all corners
        //         colIdx = m_mapOriginal.index(colIndices[l],patches[l]);
        //         entries.push_back(std::make_tuple(extraRow,colIdx,0.25));
        //         // m_matrix(rowIdx,colIdx) = 0.25;
        //     }

        //     // m_basisCheck[rowIdx] = true;
        // }

        // rowIndices.resize(2);
        // colIndices.resize(2);
        // patches.resize(2);

        // // extra point handling
        // colIndices.resize(2);
        // for (std::vector<patchSide>::iterator it = boundaries.begin(); it!=boundaries.end(); ++it)
        // {
        //     // find which corner of the interface
        //     it->getContainedCorners(d,temp_corners);
        //     for (std::vector<patchCorner>::iterator corn = temp_corners.begin(); corn!=temp_corners.end(); ++corn)
        //     {
        //         if ( std::find(corners.begin(), corners.end(), *corn) == corners.end() ) // the contained corner is not in corners
        //             continue;

        //         index_t colIdx = this->_indexFromVert(1,*corn,*it,0);
        //         colIdx = m_mapOriginal.index(colIdx,corn->patch);
        //         entries.push_back(std::make_tuple(extraRow,colIdx,0.5));
        //         // m_matrix(extraRow,colIdx) = 0.5;

        //     }
        // }

        // // m_basisCheck[extraRow] = true;
        // // }
        // // Interface handling
        // for (std::vector<boundaryInterface>::iterator it = ifaces.begin(); it!=ifaces.end(); ++it)
        // {
        //     // if (!m_topology.isInterface(it->first()))
        //     //     continue;

        //     // find which corner of the interface
        //     it->first().getContainedCorners(d,temp_corners);

        //     std::vector<patchSide> isides(2);
        //     std::vector<patchCorner> icorners(2);

        //     for (std::vector<patchCorner>::iterator corn = temp_corners.begin(); corn!=temp_corners.end(); ++corn)
        //     {
        //         if ( std::find(corners.begin(), corners.end(), *corn) == corners.end() ) // the contained corner is not in corners
        //         {
        //             continue;
        //         }

        //         // Now we need to fill the matrix for the rows corresponding to the (1,1) coefficients of the corners
        //         icorners[0] = *corn;
        //         icorners[1] = it->mapCorner(*corn);
        //         isides[0] = it->first();
        //         isides[1] = it->second();

        //         // get rowIndices
        //         for (size_t k=0; k!=icorners.size(); k++)
        //         {
        //             rowIndices[k] = this->_indexFromVert(1,icorners[k],isides[k],1);
        //             colIndices[k] = this->_indexFromVert(1,icorners[k],isides[k],0);
        //             patches[k] = icorners[k].patch;
        //         }

        //         for (size_t k = 0; k!=rowIndices.size(); k++)
        //         {
        //             index_t rowIdx = m_mapModified.index(rowIndices[k],patches[k]);
        //             // Fill the matrix entries related to the 0,0 coefs (stored in indices) 0.25 for all corners
        //             for (size_t l = 0; l!=colIndices.size(); l++)
        //             {
        //                 index_t colIdx = m_mapOriginal.index(colIndices[l],patches[l]);
        //                 entries.push_back(std::make_tuple(rowIdx,colIdx,0.5));
        //                 // m_matrix(rowIdx,colIdx) = 0.5;
        //             }
        //             // m_basisCheck[rowIdx] = true;
        //         }
        //     }
        // }

        // #pragma omp critical (handle_boundary_vertex_ff)
        // {
        //     _pushAndCheck(entries);

        //     // Furthermore, all adjacent vertices are checked
        //     for (std::vector<patchCorner>::iterator it = corners.begin(); it!=corners.end(); ++it)
        //         m_vertCheck[ this->_vertIndex(it->patch, it->corner()) ] = true;
        // }

        // 2. make container for the interfaces
        std::vector<boundaryInterface> ifaces;
        boundaryInterface iface;
        std::vector<patchSide> boundaries;
        index_t extraRow;
        std::vector<patchCorner> temp_corners, corners;
        std::vector<patchSide> psides;
        std::vector<index_t> colIndices,rowIndices, patches;

        sparseEntry_t entries;

        m_topology.getCornerList(pcorner,corners);

        // if (!(std::count(m_C0s.begin(), m_C0s.end(), pcorner)))
        // {
        if ((std::count(m_C0s.begin(), m_C0s.end(), pcorner)))
            gsWarn<<"C0 handling for boundary corners with valence 3 has not yet been implemented\n";


        /*
            Warning:
            This case handles all patchCorners related to this vertex at once!
            */

        // 1. get all adjacent vertices

        // find corner with the lowest patch index
        patchCorner lowest = corners[0];
        for (size_t k=0; k!=corners.size(); k++)
            lowest = (corners[k].patch < lowest.patch ) ? corners[k] : lowest;

        lowest.getContainingSides(d,psides);
        // get (0,0) index from the corner with lowest patch number and store the row index.
        extraRow = this->_indexFromVert(0,lowest,psides[0],0);
        extraRow = m_mapModified.index(extraRow,lowest.patch);

        // 2. loop over the adjacent vertices
        colIndices.clear();
        rowIndices.clear();
        patches.clear();
        for (std::vector<patchCorner>::iterator it = corners.begin(); it!=corners.end(); ++it)
        {
            // *it is a patchCorner
            // 3. determine if one of the contained sides is a boundary
            it->getContainingSides(d,psides);

            // store the 0,0 indices
            colIndices.push_back( this->_indexFromVert(0,*it,psides[0],0) );
            rowIndices.push_back( this->_indexFromVert(1,*it,psides[0],1) );
            patches.push_back(it->patch);

            for (index_t k = 0; k!=2; k++) // index of point over the interface
            {
                if ( m_topology.getInterface(psides[k],iface) ) // if it is an interface, store it
                {
                    ifaces.push_back(iface);
                }
                else                                            // if not, then store the side
                {
                    //
                    index_t colIdx = m_mapOriginal.index(this->_indexFromVert(1,*it,psides[k],0),it->patch);
                    index_t rowIdx = m_mapModified.index(this->_indexFromVert(1,*it,psides[k],1),it->patch);
                    // m_matrix(rowIdx,colIdx) = 0.5;
                    // m_matrix(extraRow,colIdx) = 0.5;
                    // m_basisCheck[rowIdx] = true;
                    entries.push_back(std::make_tuple(rowIdx,colIdx,0.5));
                    entries.push_back(std::make_tuple(extraRow,colIdx,0.5));
                    // boundaries.push_back(psides[k]);
                }
            }
        }

        // GISMO_ASSERT(boundaries.size()==2,"There must be two boundaries that are not an interface!");
        // ifaces.push_back(boundaryInterface(boundaries[0],boundaries[1],d));




        // the extra (0,0) node gets 0.25 for all (0,0) entries
        for (size_t k = 0; k!=colIndices.size(); k++)
        {
            index_t colIdx = m_mapOriginal.index(colIndices[k],patches[k]);
            // m_matrix(extraRow,colIdx) = 0.25;
            entries.push_back(std::make_tuple(extraRow,colIdx,0.25));

        }
        // Fill the matrix entries related to the 1,1 coefs (stored in indices) with 1 for itself and with 0.25 for the others
        for (size_t k = 0; k!=rowIndices.size(); k++)
        {
            index_t rowIdx = m_mapModified.index(rowIndices[k],patches[k]);
            index_t colIdx = m_mapOriginal.index(rowIndices[k],patches[k]);
            // Fill the matrix entries related to itself (1,1) with a 1.0
            entries.push_back(std::make_tuple(rowIdx,colIdx,1.0));
            // m_matrix(rowIdx,colIdx) = 1.0;

            for (size_t l = 0; l!=colIndices.size(); l++)
            {
                // Fill the matrix entries related to the 0,0 coefs (stored in indices) 0.25 for all corners
                colIdx = m_mapOriginal.index(colIndices[l],patches[l]);
                // m_matrix(rowIdx,colIdx) = 0.25;
                entries.push_back(std::make_tuple(rowIdx,colIdx,0.25));
            }

            // m_basisCheck[rowIdx] = true;
        }

        rowIndices.resize(2);
        colIndices.resize(2);
        patches.resize(2);

        // extra point handling
        colIndices.resize(2);
        for (std::vector<patchSide>::iterator it = boundaries.begin(); it!=boundaries.end(); ++it)
        {
            // find which corner of the interface
            it->getContainedCorners(d,temp_corners);
            for (std::vector<patchCorner>::iterator corn = temp_corners.begin(); corn!=temp_corners.end(); ++corn)
            {
                if ( std::find(corners.begin(), corners.end(), *corn) == corners.end() ) // the contained corner is not in corners
                    continue;

                index_t colIdx = this->_indexFromVert(1,*corn,*it,0);
                colIdx = m_mapOriginal.index(colIdx,corn->patch);
                entries.push_back(std::make_tuple(extraRow,colIdx,0.5));
                // m_matrix(extraRow,colIdx) = 0.5;

            }
        }

        // m_basisCheck[extraRow] = true;
        // }
        // Interface handling
        for (std::vector<boundaryInterface>::iterator it = ifaces.begin(); it!=ifaces.end(); ++it)
        {
            // if (!m_topology.isInterface(it->first()))
            //     continue;

            // find which corner of the interface
            it->first().getContainedCorners(d,temp_corners);

            std::vector<patchSide> isides(2);
            std::vector<patchCorner> icorners(2);

            for (std::vector<patchCorner>::iterator corn = temp_corners.begin(); corn!=temp_corners.end(); ++corn)
            {
                if ( std::find(corners.begin(), corners.end(), *corn) == corners.end() ) // the contained corner is not in corners
                {
                    continue;
                }

                // Now we need to fill the matrix for the rows corresponding to the (1,1) coefficients of the corners
                icorners[0] = *corn;
                icorners[1] = it->mapCorner(*corn);
                isides[0] = it->first();
                isides[1] = it->second();

                // get rowIndices
                for (size_t k=0; k!=icorners.size(); k++)
                {
                    rowIndices[k] = this->_indexFromVert(1,icorners[k],isides[k],1);
                    colIndices[k] = this->_indexFromVert(1,icorners[k],isides[k],0);
                    patches[k] = icorners[k].patch;
                }

                for (size_t k = 0; k!=rowIndices.size(); k++)
                {
                    index_t rowIdx = m_mapModified.index(rowIndices[k],patches[k]);
                    // Fill the matrix entries related to the 0,0 coefs (stored in indices) 0.25 for all corners
                    for (size_t l = 0; l!=colIndices.size(); l++)
                    {
                        index_t colIdx = m_mapOriginal.index(colIndices[l],patches[l]);
                        entries.push_back(std::make_tuple(rowIdx,colIdx,0.5));
                        // m_matrix(rowIdx,colIdx) = 0.5;
                    }
                    // m_basisCheck[rowIdx] = true;
                }
            }

        }

        #pragma omp critical (handle_boundary_vertex_ff)
        {
            _pushAndCheck(entries);

            // Furthermore, all adjacent vertices are checked
            for (std::vector<patchCorner>::iterator it = corners.begin(); it!=corners.end(); ++it)
                m_vertCheck[ this->_vertIndex(it->patch, it->corner()) ] = true;
        }
    }

    template<short_t d,class T>
    void gsDPatch<d,T>::_handleIrregularBoundaryVertexNonSmooth(patchCorner pcorner, index_t valence)
    {
        gsWarn<<"C0 handling for boundary corners with valence 3 has not yet been implemented. Using the default approach\n";
        this->_handleIrregularBoundaryVertexSmooth(pcorner,valence);
    }

    template<short_t d,class T>
    void gsDPatch<d,T>::_computeVertexMapper(patchCorner pcorner)
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
        else if (!interior && valence==3 && C0)
            _computeMapperIrregularBoundaryVertexNonSmooth_v3(pcorner,valence);
        else if (!interior && valence==3 && !C0)
            _computeMapperIrregularBoundaryVertexSmooth_v3(pcorner,valence);
        else if (!interior && valence >3 && C0)
            _computeMapperIrregularBoundaryVertexNonSmooth_v(pcorner,valence);
        else if (!interior && valence >3 && !C0)
            _computeMapperIrregularBoundaryVertexSmooth_v(pcorner,valence);
        else if (interior)
            _computeMapperInteriorVertex_v(pcorner,valence);
        else
            GISMO_ERROR("Something went terribly wrong, interior="<<interior<<"; valence="<<valence);

        // label vertex as processed
        m_vertCheck[ cidx ] = true;
    }

    // Different for DPatch and AlmostC1
    template<short_t d,class T>
    void gsDPatch<d,T>::_computeMapperIrregularBoundaryVertexSmooth_v3(patchCorner pcorner, index_t valence)
    {
        std::vector<patchSide> psides(2);
        std::vector<patchCorner> pcorners;
        /*
            o o o @ X |i| X @ o o o                 x: eliminated DoFs by vertex rule
            o o o @ X |n| X @ o o o                 X: eliminated DoFs by interface/boundary rule
            o o o @ X |t| X @ o o o                 o: preserved DoFs (interior)
            @ @ @ * x |e| x * @ @ @                 @: modified DoFs by interface rule
            X X X x % |r| % x X X X                 *: modified DoFs by vertex rule (unique DoFs)
            -----------------------                 %: modified DoFs by vertex rule (matched DoFs)
            -boundary-| |-interface
            -----------------------
                     b| | % x X X X
                     o| | x * @ @ @
                     u| | X @ o o o
                     n| | X @ o o o
                     d| | X @ o o o

        */
        // if (!(std::count(m_C0s.begin(), m_C0s.end(), pcorner)))
        // {

        // We handle all corners associated to pcorner
        m_topology.getCornerList(pcorner,pcorners);

        for (size_t c=0; c!=pcorners.size(); c++)
        {
            // Eliminate their 0,1 and 1,0 vertices
            pcorners[c].getContainingSides(d,psides);
            m_mapModified.eliminateDof(this->_indexFromVert(1,pcorners[c],psides[0]),pcorners[c].patch);
            m_mapModified.eliminateDof(this->_indexFromVert(1,pcorners[c],psides[1]),pcorners[c].patch);

            // And match the 0,0 vertex (i.e. the corner) to the corner that is first in the list pcorners.
            patchSide pseudo; // this side does not contribute since we use index = 0 in this->_indexFromVert
            if (c!=0)
                m_mapModified.matchDof(pcorners[0].patch,this->_indexFromVert(0,pcorners[0],pseudo,0),pcorners[c].patch,this->_indexFromVert(0,pcorners[c],pseudo,0));
            // mark the vertex as passed
            m_vertCheck[ this->_vertIndex(pcorners[c].patch, pcorners[c].corner()) ] = true;
        }
        // }
        // else
        // {
            // pcorner.getContainingSides(d,psides);
            // for (size_t p=0; p!=psides.size(); p++)
            // {
            //     if (m_topology.isInterface(psides[p]))
            //     {
            //             m_mapModified.eliminateDof(this->_indexFromVert(1,pcorner,psides[p],0),pcorner.patch);
            //     }
            // }
        // }
    }

    // Different for DPatch and AlmostC1
    template<short_t d,class T>
    void gsDPatch<d,T>::_computeMapperIrregularBoundaryVertexNonSmooth_v3(patchCorner pcorner, index_t valence)
    {
        gsWarn<<"C0 handling for boundary corners with valence 3 has not yet been implemented. Using the default approach\n";
        this->_computeMapperIrregularBoundaryVertexSmooth_v3(pcorner,valence);
    }

    template<short_t d,class T>
    void gsDPatch<d,T>::_computeMapperIrregularBoundaryVertexSmooth_v(patchCorner pcorner, index_t valence)
    {
        GISMO_ERROR("Boundary vertex on patch"<<pcorner.patch<<" with index "<<pcorner.corner()<<" with valence = "<<valence<<" has no implementation");
    }

    template<short_t d,class T>
    void gsDPatch<d,T>::_computeMapperIrregularBoundaryVertexNonSmooth_v(patchCorner pcorner, index_t valence)
    {
        GISMO_ERROR("Boundary vertex on patch"<<pcorner.patch<<" with index "<<pcorner.corner()<<" with valence = "<<valence<<" has no implementation");
    }

} // namespace gismo
