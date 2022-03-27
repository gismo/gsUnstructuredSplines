/** @file gsDPatch.hpp

    @brief Creates the D-Patch smoothing matrix.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gsHSplines/gsHTensorBasis.h>
#include <gsHSplines/gsTHBSpline.h>
#include <gsSolver/gsBlockOp.h>
#include <gsSolver/gsMatrixOp.h>
#include <typeinfo>

#define PI 3.141592653589793

namespace gismo
{
    // Constructors
    template<short_t d,class T>
    gsDPatch<d,T>::gsDPatch(const gsMultiPatch<T> & patches)
    :
    Base(patches)
    {
        defaultOptions();
        compute();
    }

    template<short_t d,class T>
    gsDPatch<d,T>::~gsDPatch()
    {
        freeAll(m_bases);
    }

    template<short_t d,class T>
    void gsDPatch<d,T>::defaultOptions()
    {
        m_options.addSwitch("Verbose","Verbose output",false);
        m_options.addInt("Pi","Pi matrix to be applied, 0: Non-negative, 1: Idempoten",0);
    }

    /*=====================================================================================
                                    Coefficients
    =====================================================================================*/
    template<short_t d,class T>
    gsMatrix<T> gsDPatch<d,T>::_preCoefficients()
    {
        GISMO_ASSERT(m_mapModified.isFinalized(),"Mapper is not finalized, run XXXX first");

        gsMatrix<T> coefs(m_mapModified.freeSize(),m_patches.geoDim());

        index_t size;
        for (size_t p=0; p!=m_patches.nPatches(); p++) // patches
        {
            size = m_mapModified.patchSize(p);
            for (index_t k=0; k!=size; k++)
            {
                if (m_mapModified.is_free(k,p))
                    coefs.row(m_mapModified.index(k,p,0)) = m_patches.patch(p).coefs().row(k);
            }
        }

        // Correct the v=3 boundary vertices:
        std::fill(m_vertCheck.begin(), m_vertCheck.end(), false);
        for (size_t p=0; p!=m_patches.nPatches(); p++)
        {
            for (index_t c=1; c<5; c++)
            {
                index_t idx = _vertIndex(p,c);
                if(m_vertCheck[ idx] )
                    continue;

                patchCorner pcorner(p,c);
                std::pair<index_t,bool> vdata = _vertexData(pcorner); // corner c

                if (std::count(m_C0s.begin(), m_C0s.end(), pcorner))
                    continue;
                if (vdata.first==3 && !vdata.second)
                {
                    std::vector<patchCorner> otherCorners;
                    std::vector<patchSide> csides;

                    m_patches.getCornerList(pcorner,otherCorners);
                    index_t b00 = -1, b11i = -1;
                    std::vector<index_t> b11b;
                    for (std::vector<patchCorner>::iterator corner = otherCorners.begin(); corner != otherCorners.end(); corner++)
                    {
                        corner->getContainingSides(d,csides);
                        if ( m_patches.isBoundary(csides[0]) || m_patches.isBoundary(csides[1]) ) //
                            b11b.push_back(m_mapModified.index( _indexFromVert(m_Bbases,1,*corner,csides[0],1) , corner->patch) );
                        else if (b11i==-1)
                            b11i = m_mapModified.index( _indexFromVert(m_Bbases,1,*corner,csides[0],1) , corner->patch);
                        else
                            GISMO_ERROR("b11i is already assigned?");

                        if (corner==otherCorners.begin())
                        {
                            gsBasis <> * basis = &m_Bbases.basis(corner->patch);
                            b00 = m_mapModified.index( basis->functionAtCorner(corner->corner()), corner->patch );
                        }

                        idx = _vertIndex(corner->patch,corner->corner());
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
    void gsDPatch<d,T>::_makeTHB()
    {
        m_RefPatches = gsMultiPatch<T>(m_patches);
        // gsMultiPatch<T> refPatches(m_patches);
        // prepare the geometry
        std::vector<std::vector<patchCorner> > cornerLists;
        m_RefPatches.getEVs(cornerLists);

        if (cornerLists.size()!=0)
        {
            std::vector< std::vector<index_t> > elVec(m_RefPatches.nPatches());
            for (size_t v =0; v!=cornerLists.size(); v++)
                for (size_t c = 0; c!=cornerLists[v].size(); c++)
                {
                    patchCorner corner = cornerLists[v].at(c);
                    gsVector<bool> pars;
                    corner.corner().parameters_into(m_RefPatches.parDim(),pars); // get the parametric coordinates of the corner
                    gsMatrix<T> mat = pars.template cast<T>(); // cast to real coordinates

                    gsMatrix<T> boxes(m_RefPatches.parDim(),2);
                    boxes.col(0) << mat;
                    boxes.col(1) << mat;

                    gsHTensorBasis<2,T> *basis = dynamic_cast<gsHTensorBasis<2,T>*>(&m_RefPatches.basis(corner.patch));
                    std::vector<index_t> elements = basis->asElements(boxes,1);

                    elVec.at(corner.patch).insert(elVec.at(corner.patch).end(), elements.begin(), elements.end());

                    // gsHTensorBasis<2,T> *basis = dynamic_cast<gsHTensorBasis<2,T>*>(&m_RefPatches.basis(corner.patch));

                    // basis->refineElements(elements, m_tMatrix);
                    // gsDebugVar(m_tMatrix.toDense());
                }

            gsSparseMatrix<T> tmp;
            index_t rows = 0, cols = 0;
            std::vector<Eigen::Triplet<T,index_t>> tripletList;
            for (size_t p=0; p!=m_RefPatches.nPatches(); p++)
            {
                gsHTensorBasis<2,T> *basis = dynamic_cast<gsHTensorBasis<2,T>*>(&m_RefPatches.basis(p));
                std::vector< gsSortedVector< index_t > > xmat = basis->getXmatrix();

                m_RefPatches.patch(p).refineElements(elVec[p]);

                basis->transfer(xmat,tmp);

                for (index_t i = 0; i<tmp.outerSize(); ++i)
                    for (typename gsSparseMatrix<T>::iterator it(tmp,i); it; ++it)
                        tripletList.push_back(Eigen::Triplet<T,index_t>(it.row()+rows,it.col()+cols,it.value()));

                rows += tmp.rows();
                cols += tmp.cols();
            }

            m_tMatrix.resize(rows,cols);
            m_tMatrix.setFromTriplets(tripletList.begin(), tripletList.end());

            m_tMatrix.makeCompressed();
            m_bases = gsMultiBasis<T>(m_RefPatches);
        }

        // redefine the mappers
        // m_mapModified = gsDofMapper(m_bases);
        m_mapOriginal = gsDofMapper(m_bases);
        m_mapOriginal.finalize();
    }

    template<short_t d,class T>
    void gsDPatch<d,T>::_computeEVs()
    {
        /*
            Our goal is to create three vectors c11, c12, c21 which all contain the
            c11, c12 and c21 coefficients of the patches around the EV in the right order
            (counter)-clockwise.
        */

        std::vector<std::vector<patchCorner> > cornerLists;
        m_patches.getEVs(cornerLists);

        if (cornerLists.size()!=0)
        {
            // gsWriteCSV(m_matrix.toDense(),"matrix0.csv");
            m_matrix = m_matrix * m_tMatrix.transpose();
            // gsWriteCSV(m_tMatrix.toDense(),"matrixTHB.csv");
            // gsWriteCSV(m_matrix.toDense(),"matrix1.csv");

            gsSparseMatrix<> pi;
            std::vector<patchSide> sides(2);
            std::vector<patchSide> allSides;
            std::vector<std::vector<patchCorner> > cornerLists;
            std::vector<patchCorner> corners;
            std::vector<index_t> allPatches;
            std::map<index_t,index_t> patches;
            std::vector<boundaryInterface> interfaces;
            m_patches.getEVs(cornerLists);

            if (cornerLists.size()!=0)
            {
                for (size_t v =0; v!=cornerLists.size(); v++) // over EVs
                {
                    patches.clear();
                    index_t N = cornerLists[v].size();

                    allPatches.resize(m_patches.nPatches());
                    corners.resize(N);
                    interfaces.resize(N);

                    gsMatrix<index_t> c11(N,1);
                    gsMatrix<index_t> c12(N,1);
                    gsMatrix<index_t> c21(N,1);
                    gsMatrix<index_t> c11o(N,1);
                    gsMatrix<index_t> c12o(N,1);
                    gsMatrix<index_t> c21o(N,1);


                    /*
                        First, we loop over all the interfaces to construct a (counter)clock-wise map for our coefficients
                        Looping clock-wise or counter clock-wise does not matter, as long as the coefficients are neighboring
                    */
                    // Loop over all sides such that we can fill c12 and c21. We just start from side 0 of corner 0
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
                        GISMO_ASSERT(m_patches.getInterface(side,interfaces[i]),"Side must be an interface!");

                        std::vector<boxCorner> adjcorners;
                        m_patches.getNeighbour(side,otherSide);
                        otherSide.getContainedCorners(d,adjcorners);
                        for (index_t k=0; k!=N; k++)
                        {
                            if (cornerLists[v][k] == patchCorner(otherSide.patch,adjcorners[0]))
                                otherCorner = patchCorner(otherSide.patch,adjcorners[0]);
                            else if (cornerLists[v][k] == patchCorner(otherSide.patch,adjcorners[1]))
                                otherCorner = patchCorner(otherSide.patch,adjcorners[1]);
                            else continue;
                        }
                        GISMO_ASSERT(otherCorner!=patchCorner(0,0),"Error");

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

                        c11(patches[side.patch],0) = _indexFromVert(1,corners[i],side,1,0);
                        c12(patches[side.patch],0) = _indexFromVert(2,corners[i],side,1,0);
                        c21(patches[otherSide.patch],0) = _indexFromVert(2,otherCorner,otherSide,1,0);

                        c11o(patches[side.patch],0) = _indexFromVert(m_Bbases,1,corners[i],side,1);
                        c12o(patches[side.patch],0) = _indexFromVert(m_Bbases,2,corners[i],side,1);
                        c21o(patches[otherSide.patch],0) = _indexFromVert(m_Bbases,2,otherCorner,otherSide,1);

                    }

                    // to do: integrate this loop
                    gsMatrix<index_t> rowIndices(3*N,1);
                    for (index_t i = 0; i!=N; i++)
                    {
                        corner = corners[i];
                        corner.getContainingSides(d,sides);
                        // we look for the 1,1 index so it does not matter which side we use
                        // rowIndices(i,0) = m_mapModified.index(_indexFromVert(1,corner,sides[0],1,-1),corner.patch);
                        rowIndices(i,0)     = m_mapModified.index(c11o(i,0),corners[i].patch);
                        rowIndices(i+N,0)   = m_mapModified.index(c12o(i,0),corners[i].patch);
                        rowIndices(i+2*N,0) = m_mapModified.index(c21o(i,0),corners[i].patch);

                        c11(i,0) = m_mapOriginal.index(c11(i,0),corners[i].patch);
                        c12(i,0) = m_mapOriginal.index(c12(i,0),corners[i].patch);
                        c21(i,0) = m_mapOriginal.index(c21(i,0),corners[i].patch);
                    }

                    gsMatrix<T> Pi = _makePi(N);

                    gsVector<T> c(3*N);
                    index_t colIdx = 0, idx = 0;
                    for (index_t i = 0; i!=N; i++) // for all involved corners
                    {
                        for (index_t k=0; k!=3; k++)// c11, c12, c21
                        {
                            for (index_t j=0; j!=N; j++) // loop over the connected corners
                            {
                                c.at(j) = m_matrix(rowIndices(i+k*N,0),c11(j,0));
                                c.at(j+N) = m_matrix(rowIndices(i+k*N,0),c21(j,0));
                                c.at(j+2*N) = m_matrix(rowIndices(i+k*N,0),c12(j,0));
                            }
                            c = Pi * c;
                            for (index_t j=0; j!=N; j++) // loop over the connected corners
                            {
                                m_matrix(rowIndices(i+k*N,0),c11(j,0)) = c.at(j);
                                m_matrix(rowIndices(i+k*N,0),c21(j,0)) = c.at(j+N);
                                m_matrix(rowIndices(i+k*N,0),c12(j,0)) = c.at(j+2*N);
                            }
                        }
                    }

                    // smoothing center
                    for (index_t i = 0; i!=N; i++) // for all involved corners
                    {
                        for (index_t k=0; k!=3; k++)
                        {
                            for (index_t j=0; j!=N; j++) // loop over the connected corners
                            {
                                corners[j].getContainingSides(d,sides);
                                colIdx = _indexFromVert(0,corners[j],sides[0],0,0); // 0,0
                                colIdx = m_mapOriginal.index(colIdx,corners[j].patch);
                                m_matrix(rowIndices(i+k*N,0),colIdx) = m_matrix.coeff(rowIndices(i+k*N,0),c11(i,0)); // by construction >>>> PROBLEMATIC

                                colIdx = _indexFromVert(1,corners[j],sides[0],0,0); // 1,0
                                colIdx = m_mapOriginal.index(colIdx,corners[j].patch);
                                m_matrix(rowIndices(i+k*N,0),colIdx) = m_matrix.coeff(rowIndices(i+k*N,0),c11(i,0)); // by construction >>>> PROBLEMATIC

                                colIdx = _indexFromVert(1,corners[j],sides[1],0,0); // 0,1
                                colIdx = m_mapOriginal.index(colIdx,corners[j].patch);
                                m_matrix(rowIndices(i+k*N,0),colIdx) = m_matrix.coeff(rowIndices(i+k*N,0),c11(i,0)); // by construction >>>> PROBLEMATIC
                            }
                        }
                    }

                    // interface smoothing
                    for (index_t i = 0; i!=N; i++) // for all involved interfaces
                    {
                        for (index_t k = 2; k!=4 ; k++)// std::max(basis1->maxDegree(),basis2->maxDegree())+
                        {
                            patchCorner corner = corners[patches[interfaces[i][0].patch]];
                            patchCorner otherCorner = corners[patches[interfaces[i][1].patch]];
                            patchSide side = interfaces[i][0];
                            patchSide otherSide = interfaces[i][1];

                            idx = _indexFromVert(k,corner,side,0,0);
                            index_t j0k = m_mapOriginal.index(idx,side.patch);

                            idx = _indexFromVert(k,otherCorner,otherSide,0,0);
                            index_t jk0 = m_mapOriginal.index(idx,otherSide.patch);

                            idx = _indexFromVert(k,corner,side,1,0);         // point (k,0)
                            index_t jk1 = m_mapOriginal.index(idx,side.patch); // point (k,0)

                            idx = _indexFromVert(k,otherCorner,otherSide,1,0);         // point (k,0)
                            index_t j1k = m_mapOriginal.index(idx,otherSide.patch); // point (k,0)

                            index_t row;
                            for (index_t r=0; r!=rowIndices.rows(); r++)
                            {
                                row = rowIndices(r,0);
                                m_matrix(row,j0k) = 0.5 * ( m_matrix.coeff(row,jk1) + m_matrix.coeff(row,j1k) );
                                m_matrix(row,jk0) = 0.5 * ( m_matrix.coeff(row,jk1) + m_matrix.coeff(row,j1k) ); //  >>>> PROBLEMATIC
                            }

                        }
                    }
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

        T phi = 2*PI / valence;
        std::complex<T> I(1,1);
        T beta = 0.4;
        T psi = std::arg( (T(1.0)+I*beta*T(math::sin(phi)) ) * std::complex<T>(math::exp( -I*phi / T(2.0) ) ));

        if (m_options.getInt("Pi")==0)
        {
            for (index_t j=0; j!=valence; j++)
            {
                P(j,0) = P(j,1) = P(j,2) = P(j,3) = P(j,6) = 1.0 / (3.0 * valence);;
                P(j,4) = P(j,8) = 1.0 / (3.0 * valence) * ( 1.0 + 3.0*math::cos( j * phi ) );
                P(j,5) = 1.0 / (3.0 * valence) * ( 1.0 + 3.0*math::cos( 2.0 * psi + j * phi ) );
                P(j,7) = 1.0 / (3.0 * valence) * ( 1.0 + 3.0*math::cos( 2.0 * psi - j * phi ) );
            }
        }
        else if (m_options.getInt("Pi")==1)
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
            GISMO_ERROR("Pi ption unknown");

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
        for (size_t p=0; p!=m_patches.nPatches(); p++)
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
        for(gsBoxTopology::const_iiterator iit = m_patches.iBegin(); iit!= m_patches.iEnd(); iit++)
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
        for(gsBoxTopology::const_biterator bit = m_patches.bBegin(); bit!= m_patches.bEnd(); bit++)
        {
            basis1 = &m_bases.basis(bit->patch);
            tmp += (basis1->boundaryOffset(bit->side(),0).size() - 4);
            tmp += (basis1->boundaryOffset(bit->side(),1).size() - 4);
        }
        // gsDebug<<"Number of boundary DoFs: "<<tmp<<"\n";
        m_size += tmp;

        // vertices
        tmp = 0;
        std::vector<bool> passed(m_patches.nPatches()*4);
        std::fill(passed.begin(), passed.end(), false);

        std::vector<patchCorner> corners;
        // index_t corn = 0;
        for (size_t p=0; p!=m_patches.nPatches(); p++)
            for (index_t c=1; c<5; c++)
            {
                index_t idx = _vertIndex(p,c);
                if (!passed.at(idx))
                {
                    m_patches.getCornerList(patchCorner(p,c),corners);

                    for (size_t k=0; k!=corners.size(); k++)
                        passed.at(_vertIndex(corners[k].patch,corners[k])) = true;

                    std::pair<index_t,bool> vdata = _vertexData(patchCorner(p,c)); // corner c
                    if (!vdata.second) // boundary vertex
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
    void gsDPatch<d,T>::_computeMapper() // also initialize the mappers!
    {
        // interfaces
        _resetChecks(false);

        gsBasis<T> * basis;

        std::vector<index_t> patches(2);
        std::vector<patchSide> psides(2);
        gsVector<index_t> indices;
        std::vector<patchCorner> pcorners;
        patchCorner pcorner;
        index_t cidx, sidx;
        std::pair<index_t,bool> vdata1, vdata2, vdata;

        for(gsBoxTopology::const_iiterator iit = m_patches.iBegin(); iit!= m_patches.iEnd(); iit++)
        {
            sidx = _sideIndex( iit->second().patch,iit->second().side());
            if (m_sideCheck.at(sidx))
                continue;

            patches[0] = iit->first().patch;
            patches[1] = iit->second().patch;
            psides[0] = patchSide(iit->first().patch,iit->first().side()); // the interface on the first patch
            psides[1] = patchSide(iit->second().patch,iit->second().side()); // the interface on the second patch

            for (index_t p = 0; p != 2; p++)
            {
                sidx = _sideIndex( patches[p] ,psides[p] );
                if (m_sideCheck.at(sidx))
                    continue;

                /*
                    Eliminates the interior nodes on the interfaces

                    o o o @ X | X @ o o o               X: eliminated DoFs by interface/boundary rule
                    o o o @ X | X @ o o o               @: modified DoFs by interface rule
                    o o o ? ? | ? ? o o o               o: preserved DoFs (interior)
                    o o ? x x | x x ? @ @               ?: Depends on the vertex (X and @ if not (interior vertex & valence = 4))
                    o o ? x x | x x ? X X               x: handled in vertex rule
                    ----------|----------
                              | x x ? X X
                              | x x ? @ @
                              | ? ? o o o
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
                    selectedIndices.push_back(_indexFromVert(0,pcorners[c],psides[p],0,0)); // index from vertex pcorners[c] along side psides[0] with offset 0.
                    selectedIndices.push_back(_indexFromVert(1,pcorners[c],psides[p],0,0)); // index from vertex pcorners[c] along side psides[0] with offset 0.
                }

                std::sort(selectedIndices.begin(),selectedIndices.end());
                // for(size_t i=0; i < selectedIndices.size(); i++)
                //     std::cout << selectedIndices.at(i) << ' ';

                std::vector<index_t> result(allIndices.size());
                std::vector<index_t>::iterator it=std::set_difference (allIndices.begin(), allIndices.end(), selectedIndices.begin(), selectedIndices.end(), result.begin());
                result.resize(it-result.begin());

                gsAsMatrix<index_t> indices(result,result.size(),1);
                m_mapModified.markBoundary(patches[p], indices);

                // gsDebug<<"Eliminated "<<indices.transpose()<<" of basis "<<patches[p]<<"\n";
                m_sideCheck.at(sidx) = true;
            }
        }

        // boundaries
        for(gsBoxTopology::const_biterator bit = m_patches.bBegin(); bit!= m_patches.bEnd(); bit++)
        {
            sidx = _sideIndex(bit->patch,bit->side());
            m_sideCheck.at(sidx) = true;
        }

        // vertices
        for (size_t p=0; p!=m_patches.nPatches(); p++)
        {
            for (index_t c=1; c<5; c++)
            {
                cidx = _vertIndex(p,c);
                if (m_vertCheck.at(cidx))
                    continue;

                pcorner = patchCorner(p,c);

                basis = &m_bases.basis(p);

                // get valence and vertex info of corner
                std::pair<index_t,bool> vdata = _vertexData(pcorner); // corner c

                if (!vdata.second) // boundary vertex
                {
                    if (vdata.first==2) //valence = 2
                    {
                        /*
                            o o o @ X |i| X @ o o o                 x: eliminated DoFs by vertex rule
                            o o o @ X |n| X @ o o o                 X: eliminated DoFs by interface/boundary rule
                            o o o @ X |t| X @ o o o                 o: preserved DoFs (interior)
                            o o o * x |e| x * o o o                 @: modified DoFs by interface rule
                            o o o * x |r| x * o o o                 *: modified DoFs by vertex rule (unique DoFs)
                            -----------------------                 %: modified DoFs by vertex rule (matched DoFs)
                            -boundary-| | -boundary-
                            -----------------------
                        */
                        // we mark the nodes belonging to the interface
                        pcorner.getContainingSides(d,psides);
                        for (size_t p=0; p!=psides.size(); p++)
                        {
                            if (m_patches.isInterface(psides[p]))
                            {
                                if (std::count(m_C0s.begin(), m_C0s.end(), pcorner))
                                {
                                    m_mapModified.eliminateDof(_indexFromVert(1,pcorner,psides[p],0),pcorner.patch);
                                }
                                else
                                {
                                    m_mapModified.eliminateDof(_indexFromVert(0,pcorner,psides[p],0),pcorner.patch);
                                    m_mapModified.eliminateDof(_indexFromVert(1,pcorner,psides[p],0),pcorner.patch);
                                    // gsDebug<<"Eliminated "<<_indexFromVert(1,pcorner,psides[p],0)<<" and "<<_indexFromVert(1,pcorner,psides[p],0)<<" of basis "<<pcorner.patch<<"\n";
                                }
                            }
                        }
                    }
                    else if (vdata.first==3) //valence = 3
                    {
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
                        if ((std::count(m_C0s.begin(), m_C0s.end(), pcorner)))
                            gsWarn<<"C0 handling for boundary corners with valence 3 has not yet been implemented\n";

                        // We handle all corners associated to pcorner
                        m_patches.getCornerList(pcorner,pcorners);

                        for (size_t c=0; c!=pcorners.size(); c++)
                        {

                            // Eliminate their 0,1 and 1,0 vertices
                            pcorners[c].getContainingSides(d,psides);
                            m_mapModified.eliminateDof(_indexFromVert(1,pcorners[c],psides[0]),pcorners[c].patch);
                            m_mapModified.eliminateDof(_indexFromVert(1,pcorners[c],psides[1]),pcorners[c].patch);

                            // And match the 0,0 vertex (i.e. the corner) to the corner that is first in the list pcorners.
                            if (c!=0)
                            {
                                patchSide pseudo = patchSide(pcorner.patch,1); // this side does not contribute since we use index = 0 in _indexFromVert
                                // m_mapModified.eliminateDof(_indexFromVert(0,pcorners[c],pseudo),pcorners[c].patch);
                                m_mapModified.matchDof(pcorners[0].patch,_indexFromVert(0,pcorners[0],pseudo),pcorners[c].patch,_indexFromVert(0,pcorners[c],pseudo));
                                // gsDebug<<"Matched "<<_indexFromVert(0,pcorners[0],pseudo)<<" of basis "<<pcorners[0].patch<<" with "<<_indexFromVert(0,pcorners[c],pseudo)<<" of basis "<<pcorners[c].patch<<"\n";

                            }
                            // mark the vertex as passed
                            m_vertCheck[ _vertIndex(pcorners[c].patch, pcorners[c].corner()) ] = true;
                        }
                        // }
                        // else
                        // {
                            // pcorner.getContainingSides(d,psides);
                            // for (size_t p=0; p!=psides.size(); p++)
                            // {
                            //     if (m_patches.isInterface(psides[p]))
                            //     {
                            //             m_mapModified.eliminateDof(_indexFromVert(1,pcorner,psides[p],0),pcorner.patch);
                            //     }
                            // }
                        // }
                    }
                    else if (vdata.first==1) //valence = 1
                    {  } // do nothing
                    else
                    {
                        // GISMO_ERROR("Boundary vertex with valence = "<<vdata.first<<" has no implementation");
                    }
                }
                else // interior vertex
                {
                    /*
                        o o o @ X |i| X @ o o o                 x: eliminated DoFs by vertex rule
                        o o o @ X |n| X @ o o o                 X: eliminated DoFs by interface/boundary rule
                        o o o @ X |t| X @ o o o                 o: preserved DoFs (interior)
                        @ @ @ * x |e| x * @ @ @                 @: modified DoFs by interface rule
                        X X X x x |r| x x X X X                 *: modified DoFs by vertex rule (unique DoFs)
                        -----------------------
                        -boundary-| |-interface
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
                        m_mapModified.eliminateDof(_indexFromVert(0,pcorner,psides[p],0),pcorner.patch);
                        m_mapModified.eliminateDof(_indexFromVert(1,pcorner,psides[p],0),pcorner.patch);
                    }
                }

                // label vertex as processed
                m_vertCheck[ cidx ] = true;
            }
        }
        m_mapModified.finalize();
        m_mapOriginal.finalize();

        m_matrix.resize( m_mapModified.freeSize(), m_mapOriginal.freeSize() );

        // gsDebugVar(m_mapModified.coupledSize());
        // gsDebugVar(m_mapModified.boundarySize());
    }

    template<short_t d,class T>
    void gsDPatch<d,T>::_handleVertex(patchCorner pcorner)
    {

        if (m_vertCheck[ _vertIndex(pcorner.patch, pcorner.corner()) ])
        {
            // gsDebug<<"corner "<<pcorner.corner()<<" ("<<pcorner.patch<<") skipped!\n";
            return;
        }
        // if (std::count(m_C0s.begin(), m_C0s.end(), pcorner))
        // {
        //     gsInfo<<"patch = "<<pcorner.patch<<", corner = "<<pcorner.corner()<<"\n";
        //     m_vertCheck[ _vertIndex(pcorner.patch, pcorner.corner()) ] = true;
        //     return;
        // }

        std::vector<patchSide> psides;
        std::vector<patchCorner> corners;
        std::vector<index_t> indices, rowIndices, colIndices, patches;
        std::vector<index_t> weights;


        std::pair<index_t,bool> vdata = _vertexData(pcorner); // corner c
        pcorner.getContainingSides(d,psides);

        if (!vdata.second) // boundary vertices
        {
            // Correct
            if (vdata.first==1)
            {
                GISMO_ASSERT(psides.size()==2,"Must have 2 adjacent sides");
                indices.resize(4);
                indices[0] = _indexFromVert(0,pcorner,psides[0],0); // b00
                indices[1] = _indexFromVert(1,pcorner,psides[0],0); // b01
                indices[2] = _indexFromVert(1,pcorner,psides[1],0); // b10
                indices[3] = _indexFromVert(1,pcorner,psides[1],1); // b11

                for (std::vector<index_t>::iterator it = indices.begin(); it!=indices.end(); ++it)
                {
                    index_t rowIdx = m_mapModified.index(*it,pcorner.patch);
                    index_t colIdx = m_mapOriginal.index(*it,pcorner.patch);
                    m_matrix(rowIdx,colIdx) = 1.0;

                    m_basisCheck[rowIdx] = true;
                }

                m_vertCheck[ _vertIndex(pcorner.patch, pcorner.corner()) ] = true;
                // gsInfo<<"patch = "<<pcorner.patch<<", corner = "<<pcorner.corner()<<"\n";
                return;
            }
            // Correct
            else if (vdata.first==2)
            {
                // 1. find the interface
                index_t iindex = m_patches.isInterface(psides[0]) ? 0 : 1;

                boundaryInterface iface;
                GISMO_ENSURE(m_patches.getInterface(psides[iindex],iface),"Must be an interface");
                // 2. collect indices
                // If we want C0 at this vertex, we only handle the row k=1.
                bool C0 = (std::count(m_C0s.begin(), m_C0s.end(), pcorner)) ? true : false;
                indices.resize(3);
                patchSide otherSide = iface.other(psides[iindex]);
                patchCorner otherCorner = iface.mapCorner(pcorner);
                for (index_t k = (C0 ? 1 : 0); k!=2; k++) // index of point over the interface
                {
                    indices[0] = _indexFromVert(k,pcorner,psides[iindex],1); // bk1 on patch of iface
                    indices[1] = _indexFromVert(k,pcorner,psides[iindex],0); // bk0 on patch of iface
                    indices[2] = _indexFromVert(k,otherCorner,otherSide,0); // bk0 on other patch

                    index_t rowIdx = m_mapModified.index(indices[0],pcorner.patch);
                    index_t colIdx = m_mapOriginal.index(indices[0],pcorner.patch);

                    m_matrix(rowIdx,colIdx) = 1.0;
                    colIdx = m_mapOriginal.index(indices[1],psides[iindex].patch);
                    m_matrix(rowIdx,colIdx) = 0.5;
                    colIdx = m_mapOriginal.index(indices[2],otherSide.patch);
                    m_matrix(rowIdx,colIdx) = 0.5;

                    m_basisCheck[rowIdx] = true;
                }


                m_vertCheck[ _vertIndex(pcorner.patch, pcorner.corner()) ] = true;
            }
            // Correct
            else if (vdata.first==3)
            {
                // 2. make container for the interfaces
                std::vector<boundaryInterface> ifaces;
                boundaryInterface iface;
                std::vector<patchSide> boundaries;
                index_t extraRow;
                std::vector<patchCorner> temp_corners;

                m_patches.getCornerList(pcorner,corners);

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
                extraRow = _indexFromVert(0,lowest,psides[0],0);
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
                    colIndices.push_back( _indexFromVert(0,*it,psides[0],0) );
                    rowIndices.push_back( _indexFromVert(1,*it,psides[0],1) );
                    patches.push_back(it->patch);

                    for (index_t k = 0; k!=2; k++) // index of point over the interface
                    {
                        if ( m_patches.getInterface(psides[k],iface) ) // if it is an interface, store it
                        {
                            ifaces.push_back(iface);
                        }
                        else                                            // if not, then store the side
                        {
                            //
                            index_t colIdx = m_mapOriginal.index(_indexFromVert(1,*it,psides[k],0),it->patch);
                            index_t rowIdx = m_mapModified.index(_indexFromVert(1,*it,psides[k],1),it->patch);
                            m_matrix(rowIdx,colIdx) = 0.5;
                            m_matrix(extraRow,colIdx) = 0.5;
                            m_basisCheck[rowIdx] = true;
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
                    m_matrix(extraRow,colIdx) = 0.25;
                }

                // Fill the matrix entries related to the 1,1 coefs (stored in indices) with 1 for itself and with 0.25 for the others
                for (size_t k = 0; k!=rowIndices.size(); k++)
                {
                    index_t rowIdx = m_mapModified.index(rowIndices[k],patches[k]);
                    index_t colIdx = m_mapOriginal.index(rowIndices[k],patches[k]);
                    // Fill the matrix entries related to itself (1,1) with a 1.0
                    m_matrix(rowIdx,colIdx) = 1.0;

                    for (size_t l = 0; l!=colIndices.size(); l++)
                    {
                        // Fill the matrix entries related to the 0,0 coefs (stored in indices) 0.25 for all corners
                        colIdx = m_mapOriginal.index(colIndices[l],patches[l]);
                        m_matrix(rowIdx,colIdx) = 0.25;
                    }

                    m_basisCheck[rowIdx] = true;
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

                        index_t colIdx = _indexFromVert(1,*corn,*it,0);
                        colIdx = m_mapOriginal.index(colIdx,corn->patch);
                        m_matrix(extraRow,colIdx) = 0.5;

                    }
                }

                m_basisCheck[extraRow] = true;
                // }
                // Interface handling
                for (std::vector<boundaryInterface>::iterator it = ifaces.begin(); it!=ifaces.end(); ++it)
                {
                    // if (!m_patches.isInterface(it->first()))
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
                            rowIndices[k] = _indexFromVert(1,icorners[k],isides[k],1);
                            colIndices[k] = _indexFromVert(1,icorners[k],isides[k],0);
                            patches[k] = icorners[k].patch;
                        }

                        for (size_t k = 0; k!=rowIndices.size(); k++)
                        {
                            index_t rowIdx = m_mapModified.index(rowIndices[k],patches[k]);
                            // Fill the matrix entries related to the 0,0 coefs (stored in indices) 0.25 for all corners
                            for (size_t l = 0; l!=colIndices.size(); l++)
                            {
                                index_t colIdx = m_mapOriginal.index(colIndices[l],patches[l]);
                                m_matrix(rowIdx,colIdx) = 0.5;
                            }
                            m_basisCheck[rowIdx] = true;
                        }
                    }
                }



                for (std::vector<patchCorner>::iterator it = corners.begin(); it!=corners.end(); ++it)
                    m_vertCheck[ _vertIndex(it->patch, it->corner()) ] = true;
            }
            else
                GISMO_ERROR("Boundaries with valence higher than 3 cannot be handled.");
        }
        else // interior vertices
        {
            m_patches.getCornerList(pcorner,corners);

            gsBasis<T> * tmpBasis;
            // find the patch corner which shares the interface

            // 1. give the basis function a weight 1 from itself
            index_t b11_p1 = _indexFromVert(1,pcorner,psides[0],1); // point 1,1 (does not matter which reference side is taken)
            index_t rowIdx = m_mapModified.index(b11_p1,pcorner.patch);
            index_t colIdx = m_mapOriginal.index(b11_p1,pcorner.patch);
            m_matrix(rowIdx,colIdx) = 1.0;

            for (std::vector<patchSide>::iterator side = psides.begin(); side != psides.end(); ++side)
            {
                boundaryInterface iface;
                GISMO_ASSERT(m_patches.getInterface(*side,iface),"Side must be an interface!");
                patchSide otherSide;
                m_patches.getNeighbour(*side,otherSide);
                patchCorner otherCorner = iface.mapCorner(pcorner);

                index_t b10_p1 = _indexFromVert(1,pcorner,*side,0); // index from vertex pcorners[c] along side psides[0] with offset 0.
                index_t b10_p2 = _indexFromVert(1,otherCorner,otherSide,0); // point 0,1

                // 2. give the basis function a weight 1/valence from the other (0,0) basis functions
                index_t index;
                for (std::vector<patchCorner>::iterator corn = corners.begin(); corn != corners.end(); ++corn)
                {
                    tmpBasis = &m_bases.basis(corn->patch);
                    index = tmpBasis->functionAtCorner(corn->corner());
                    colIdx = m_mapOriginal.index(index,corn->patch);
                    m_matrix(rowIdx,colIdx) = 1./vdata.first;
                }

                // 3. add weight 1/2 from the interface bcs.
                colIdx = m_mapOriginal.index(b10_p1,side->patch);
                m_matrix(rowIdx,colIdx) = 0.5;

                colIdx = m_mapOriginal.index(b10_p2,otherSide.patch);
                m_matrix(rowIdx,colIdx) = 0.5;

            }
            m_basisCheck[rowIdx] = true;
            m_vertCheck[ _vertIndex(pcorner.patch, pcorner.corner()) ] = true;
        }
    }

    template<short_t d,class T>
    void gsDPatch<d,T>::_handleInterface(boundaryInterface iface)
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
                selectedIndices[p].push_back(_indexFromVert(0,pcorners[c],iface[p],0,0)); // index from vertex pcorners[c] along side psides[0] with offset 0.
                selectedIndices[p].push_back(_indexFromVert(1,pcorners[c],iface[p],0,0)); // index from vertex pcorners[c] along side psides[0] with offset 0.

                selectedOIndices[p].push_back(_indexFromVert(0,pcorners[c],iface[p],1,0)); // index from vertex pcorners[c] along side psides[0] with offset 0.
                selectedOIndices[p].push_back(_indexFromVert(1,pcorners[c],iface[p],1,0)); // index from vertex pcorners[c] along side psides[0] with offset 0.
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
        // loop over adjacent patches and couple the DoFs.
        for (index_t p =0; p!= 2; p++)
        {
            np = 1-p; // not index p;
            for (index_t k=0; k!= indices[p].size(); k++ )
            {
                rowIdx = m_mapModified.index(oindices[p].at(k),iface[p].patch);
                // rowIdx1 = m_mapOriginal.index(oindices[p].at(k),patches[p]);
                colIdx = m_mapOriginal.index(oindices[p].at(k),iface[p].patch);
                m_matrix(rowIdx,colIdx) = 1.0;

                colIdx = m_mapOriginal.index(indices[p].at(k),iface[p].patch);
                m_matrix(rowIdx,colIdx) = 0.5;

                colIdx = m_mapOriginal.index(indices[np].at(k),iface[np].patch);
                m_matrix(rowIdx,colIdx) = 0.5;

                m_basisCheck[rowIdx] = true;
            }
            m_sideCheck[ _sideIndex(iface[p].patch, iface[p].side()) ] = true; // side finished
        }
    }

    template<short_t d,class T>
    void gsDPatch<d,T>::_handleBoundary(patchSide side)
    {
            std::vector<patchCorner> pcorners;
            std::vector<index_t> selectedIndices;
            gsBasis<T> * basis = &m_bases.basis(side.patch);
            gsMatrix<index_t> indices = basis->boundaryOffset(side.side(),0);
            side.getContainedCorners(d,pcorners);
            for (index_t c =0; c!=2; c++)
            {
                selectedIndices.push_back(_indexFromVert(0,pcorners[c],side,0,0)); // index from vertex pcorners[c] along side psides[0] with offset 0.
                selectedIndices.push_back(_indexFromVert(1,pcorners[c],side,0,0)); // index from vertex pcorners[c] along side psides[0] with offset 0.
            }

            std::sort(selectedIndices.begin(),selectedIndices.end());
            std::vector<index_t> allIndices(indices.data(), indices.data() + indices.rows() * indices.cols());
            std::vector<index_t> result(allIndices.size());
            std::vector<index_t>::iterator it=std::set_difference (allIndices.begin(), allIndices.end(), selectedIndices.begin(), selectedIndices.end(), result.begin());
            result.resize(it-result.begin());

            index_t rowIdx,colIdx;
            for (std::vector<index_t>::iterator it = result.begin(); it!=result.end(); ++it)
            {
                rowIdx = m_mapModified.index(*it,side.patch);
                colIdx = m_mapOriginal.index(*it,side.patch);
                m_matrix(rowIdx,colIdx) = 1.0;
                m_basisCheck[rowIdx] = true;
            }

            m_sideCheck.at( _sideIndex(side.patch,side.side()) ) = true;
    }

    template<short_t d,class T>
    void gsDPatch<d,T>::_handleInterior()
    {
        index_t rowIdx,colIdx;
        for(size_t p=0; p!=m_patches.nPatches(); p++)
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

    template<short_t d,class T>
    void gsDPatch<d,T>::_computeSmoothMatrix()
    {
        _resetChecks(true);

        // iterate over the vertices
        for (size_t p=0; p!=m_patches.nPatches(); p++)
            for (index_t c=1; c<5; c++)
                _handleVertex(patchCorner(p,c));

        for(gsBoxTopology::const_iiterator iit = m_patches.iBegin(); iit!= m_patches.iEnd(); iit++)
            _handleInterface(*iit);

        // boundaries
        for(gsBoxTopology::const_biterator bit = m_patches.bBegin(); bit!= m_patches.bEnd(); bit++)
            _handleBoundary(*bit);

        _handleInterior();

        if (m_options.getSwitch("Verbose")) { _whichHandled(); }

        _performChecks(true);

        m_matrix.makeCompressed();
    }


} // namespace gismo