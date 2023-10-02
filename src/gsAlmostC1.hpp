/** @file gsAlmostC1.hpp

    @brief Creates the D-Patch smoothing matrix.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gsIO/gsWriteParaview.h>
#include <gsHSplines/gsHTensorBasis.h>
#include <gsHSplines/gsTHBSpline.h>

#include <gsAssembler/gsExprHelper.h>
#include <gsAssembler/gsExprEvaluator.h>
#include <gsAssembler/gsAssembler.h>

namespace gismo
{
    // Constructors
    template<short_t d,class T>
    gsAlmostC1<d,T>::gsAlmostC1(const gsMultiPatch<T> & patches)
    :
    Base(patches)
    {
        this->defaultOptions();

        for (size_t p=0; p!=m_Bbases.nBases(); p++)
            for (short_t dim=0; dim!=d; dim++)
                GISMO_ENSURE(m_Bbases.basis(p).degree(dim)==2,"Degree of the basis ( dimension "<<dim<<" ) of patch "<<p<<" is "<<m_Bbases.basis(p).degree(dim)<<", but should be 2!");
    }

    template<short_t d,class T>
    gsAlmostC1<d,T>::~gsAlmostC1()
    {
        freeAll(m_bases);
    }

    // template<short_t d,class T>
    // void gsAlmostC1<d,T>::compute()
    // {
    //     // m_RefPatches = m_patches;
    //     _initialize();
    //     _computeMapper();
    //     _computeSmoothMatrix();
    //     GISMO_ASSERT(this->_checkMatrix(m_matrix),"Mapper does not have column sum equal to 1");
    //     _makeTHB();
    //     _computeEVs();
    //     GISMO_ASSERT(this->_checkMatrix(m_matrix),"Mapper does not have column sum equal to 1");
    // }


    /*=====================================================================================
                                    Special functions
    =====================================================================================*/
    template<short_t d,class T>
    gsMatrix<T> gsAlmostC1<d,T>::_getNormals(const std::vector<patchCorner> & corners) const
    {
        gsMatrix<T> normals(3,corners.size());

        gsVector<bool> pars;
        gsMatrix<T> mat;

        gsExprEvaluator<T> ev;
        typename gsExprEvaluator<T>::geometryMap Gm = ev.getMap(m_RefPatches);
        index_t k = 0;
        for (typename std::vector<patchCorner>::const_iterator it = corners.begin(); it!=corners.end(); it++, k++)
        {
            it->corner().parameters_into(m_RefPatches.parDim(),pars); // get the parametric coordinates of the corner
            gsMatrix<T> supp = m_RefPatches.basis(it->patch).support();
            gsVector<T> vec(supp.rows());
            for (index_t r = 0; r!=supp.rows(); r++)
                vec(r) = supp(r,pars(r));

            normals.col(k) = ev.eval(sn(Gm).normalized(),vec,it->patch);
        }
        return normals;
    }


    template<short_t d,class T>
    std::tuple<gsMatrix<T>,gsMatrix<T>,gsMatrix<index_t>> gsAlmostC1<d,T>::_makeTriangle(const patchCorner & corner) const
    {
        GISMO_ASSERT(m_RefPatches.nPatches()!=0,"Are the patches refined?");

        index_t tdim = m_RefPatches.targetDim();

        std::vector<patchCorner> corners;
        m_RefPatches.getCornerList(corner,corners);

        gsVector<bool> pars;
        gsMatrix<T> mat;
        // 1. Get the coordinates of the vertex and set its z coordinate to 0
        gsMatrix<T> um(3,1), midpoint;
        um.setZero();
        corner.corner().parameters_into(m_RefPatches.parDim(),pars); // get the parametric coordinates of the corner
        gsMatrix<T> supp = m_RefPatches.basis(corner.patch).support();
        gsVector<T> vec(supp.rows());
        for (index_t r = 0; r!=supp.rows(); r++)
            vec(r) = supp(r,pars(r));

        um.block(0,0,tdim,1) = m_RefPatches.patch(corner.patch).eval(vec);
        midpoint = um; // store the original midpoint

        // 2. Get the 0,0;0,1; 1,0; 1,1 coordinates
        gsMatrix<T> u(3,corners.size()*4);
        u.setZero();
        gsMatrix<index_t> uind(1,corners.size()*4);
        uind.setZero();

        std::vector<patchSide> csides;
        index_t idx;
        for (size_t c = 0; c!=corners.size(); c++)
        {
            corners[c].getContainingSides(d,csides);
            index_t k=0;
            for (index_t i=0; i!=2; i++)
                for (index_t j=0; j!=2; j++,k++)
                {
                    idx = _indexFromVert(i,corners[c],csides[0],j);
                    uind(0,4*c+k) = m_mapOriginal.index(idx,corners[c].patch);
                    u.block(0,4*c+k,m_RefPatches.targetDim(),1) = m_RefPatches.patch(corners[c].patch).coefs().row(idx).transpose();
                }
        }

        // 3. Translate all points to a coordinate system with origin um
        gsMatrix<T> up = u;
        for (index_t k=0; k!=up.cols(); k++)
            up.col(k) -= um;

        // 4. Rotate the points parallel the xy-plane and set their z-coordinates to 0
        gsMatrix<T,3,3> Rn, Rx;
        Rn.setIdentity();
        if (m_RefPatches.targetDim()==2)
        {
            // do nothing
        }
        else if(m_RefPatches.targetDim()==3)
        {
            // Get the average normal at the corner
            gsVector<T> avgnormal = _getNormals(corners).rowwise().mean();
            // If the norm is zero, the normals are likely to be opposite
            if (avgnormal.norm()==0) avgnormal = _getNormals(corners).col(0);
            // Find the rotation matrix that maps the average normal to the z axis
            gsVector<T,3> ez;
            ez<<0,0,1;
            Rn = _getRotationMatrix(avgnormal.normalized(),ez);

            for (index_t k=0; k!=up.cols(); k++)
                up.col(k).applyOnTheLeft(Rn);

            up.row(2).setZero(); // all points
            um.row(2).setZero();// midpoint
        }
        else
            GISMO_ERROR("Target dimension of the multipatch should be 2 or 3, but is "<<m_RefPatches.targetDim());

        // 5. Find the maximum distance from the midpoint to all points
        T distance, maxDistance = 0;
        gsMatrix<T> umax;
        for (index_t k = 0; k!=up.cols(); k++)
        {
            distance = (up.col(k)).norm();
            if (distance > maxDistance)
            {
                maxDistance = distance;
                umax = up.col(k);
            }
        }

        gsVector<T,3> ex;
        ex<<1,0,0;

        // 6. Rotate all points such that the maximum point is aligned with the x-axis
        Rx = _getRotationMatrix(umax.normalized(),ex);
        for (index_t k=0; k!=up.cols(); k++)
            up.col(k).applyOnTheLeft(Rx);

        // 7. Obtain the coordinates of the triangle that encloses the circle with radius maxDistance in the xy plane
        T r = maxDistance;
        T a = 1. / ( 1./6. * std::sqrt(3) ) * r;
        T rr = 1. / 3. * std::sqrt(3) * a;

        gsMatrix<T> Cp(2,3);
        Cp.col(0)<<rr,0;
        Cp.col(1)<<-r, 0.5*a;
        Cp.col(2)<<-r,-0.5*a;

        // 8. Get the barycentric coordinates of the points
        gsMatrix<T> ub = up;
        up.row(2).setOnes(); // project the parametric points to z=1
        gsMatrix<T> A(3,3);
        A.block(0,0,2,3) = Cp;
        A.row(2).setOnes();

        for (index_t k = 0; k!=ub.cols(); k++)
        {
            ub.col(k) = A.colPivHouseholderQr().solve(up.col(k));
            GISMO_ASSERT((Cp * ub.col(k)-up.col(k).head(2)).norm()<1e-12,"Something went wrong with the computation of the barycentric coordinates. (Cp * ub.col(k)-up.col(k).head(2)).norm() = "<<(Cp * ub.col(k)-up.col(k).head(2)).norm()<<"; Cp * ub.col(k) = "<<Cp * ub.col(k)<<"; up.col(k).head(2) = "<<up.col(k).head(2));
        }

        // 9. Move the corners of the triangle back to physical coordinates
        gsMatrix<T> Cg(3,3);
        Cg.setZero();
        Cg.block(0,0,2,3) = Cp;

        for (index_t k = 0; k!=Cg.cols(); k++)
        {
            Cg.col(k).applyOnTheLeft((Rx).transpose());
            Cg.col(k).applyOnTheLeft((Rn).transpose());
            Cg.col(k) += midpoint;
        }

        if (m_RefPatches.targetDim()==2)
            Cg.conservativeResize(2,gsEigen::NoChange);

        return std::make_tuple(Cg,ub,uind);
    }

    template<short_t d,class T>
    gsMatrix<T,3,3> gsAlmostC1<d,T>::_getRotationMatrix(const gsVector<T,3> & a, const gsVector<T,3> & b) const
    {
        GISMO_ASSERT(std::abs(a.norm()-1)<1e-14,"A must be a unit vector, a.norm() = "<<a.norm());
        GISMO_ASSERT(std::abs(b.norm()-1)<1e-14,"B must be a unit vector, b.norm() = "<<b.norm());

        gsVector<T,3> v = a.cross(b);
        v.normalize();
        T theta = std::acos( a.dot(b) / ( a.norm() * b.norm() ) );

        T s = std::sin(theta);
        T c = std::cos(theta);
        gsMatrix<T,3,3> R,vx,tmp, I;
        R.setZero();
        vx.setZero();

        vx.row(0)<<0,-v.at(2),v.at(1);
        vx.row(1)<<v.at(2),0,-v.at(0);
        vx.row(2)<<-v.at(1),v.at(0),0;

        I.setIdentity();
        R += I*c;
        R += vx * s;
        tmp = (v*v.transpose()) * (1-c);
        R += tmp;

        GISMO_ASSERT((R * a - b).norm() < 1e-12,"Rotation matrix is wrong, R*a = "<<R*a<<"; b = "<<b);
        return R;
    }


    /*=====================================================================================
                                    Coefficients
    =====================================================================================*/
    // ADD THE COEFFICIENTS OF THE TRIANGLES AS EXTRA COEFFICIENTS

    template<short_t d,class T>
    gsMatrix<T> gsAlmostC1<d,T>::freeCoefficients()
    {
        GISMO_ASSERT(m_mapModified.isFinalized(),"Mapper is not finalized, run XXXX first");

        GISMO_ASSERT((size_t)m_mapModified.freeSize()==m_size,"Size does not match predicted size, m_mapModified.freeSize()="<<m_mapModified.freeSize()<<"; m_size="<<m_size);
        gsMatrix<T> coefs(m_mapModified.freeSize(),m_patches.geoDim());

        index_t size;
        for (size_t p=0; p!=m_bases.nBases(); p++) // patches
        {
            size = m_mapModified.patchSize(p);
            for (index_t k=0; k!=size; k++)
            {
                if (m_mapModified.is_free(k,p))
                    coefs.row(m_mapModified.index(k,p,0)) = m_patches.patch(p).coefs().row(k);
            }
        }
        return coefs;
    }

    template<short_t d,class T>
    gsMatrix<T> gsAlmostC1<d,T>::_preCoefficients()
    {
        GISMO_ASSERT(m_mapModified.isFinalized(),"Mapper is not finalized, run XXXX first");

        gsMatrix<T> coefs = this->freeCoefficients();

        // Correct the EVs
        std::fill(m_vertCheck.begin(), m_vertCheck.end(), false);
        index_t cidx;
        std::vector<patchCorner> pcorners;
        patchCorner pcorner;
        for (size_t p=0; p!=m_bases.nBases(); p++)
        {
            for (index_t c=1; c<5; c++)
            {
                cidx = _vertIndex(p,c);
                if (m_vertCheck.at(cidx))
                    continue;

                bool C0 = m_C0s[cidx];
                pcorner = patchCorner(p,c);
                m_topology.getCornerList(pcorner,pcorners);
                std::pair<index_t,bool> vdata = _vertexData(pcorner); // corner c
                if (vdata.first > 2 && !(vdata.first==4 && vdata.second)) // valence must be 3 or larger, but it must not be an interior vertex with v=4
                {
                    // get the triangle
                    gsMatrix<T> Cg;
                    std::tie(Cg,std::ignore,std::ignore) = _makeTriangle(pcorner);

                    // The barycentric coordinates will be attached to the matrix rows corresponding to the 0,0 corners (and the three lowest patch index corners whenever valence > 3)
                    // We use _getLowestIndices such that the corners are assigned to increasing patch corners
                    std::vector<std::pair<index_t,index_t>> indices  = _getAllInterfaceIndices(pcorner,0,m_Bbases);
                    _getLowestIndices(indices,3);

                    std::vector<index_t> rowIndices;
                    for (std::vector<std::pair<index_t,index_t>>::iterator it=indices.begin(); it!=indices.end(); it++)
                    {
                        GISMO_ASSERT(m_mapModified.is_free(it->second,it->first),"This DoF must be free! patch = "<<it->first<<"; index = "<<it->first);
                        rowIndices.push_back(m_mapModified.index(it->second,it->first));
                    }

                    index_t rowIdx;
                    for (index_t j=0; j!=Cg.cols(); j++)
                    {
                        rowIdx = rowIndices[j];
                        coefs.row(rowIdx) = Cg.col(j).transpose();
                    }

                    for (size_t k = 0; k!=pcorners.size(); k++)
                        m_vertCheck[ _vertIndex(pcorners[k].patch, pcorners[k].corner()) ] = true;
                }
                else if (vdata.first == 2 && C0) // valence must be 3 or larger, but it must not be an interior vertex with v=4
                {
                    // get the triangle
                    gsMatrix<T> Cg;
                    std::tie(Cg,std::ignore,std::ignore) = _makeTriangle(pcorner);

                    // The barycentric coordinates will be attached to the matrix rows corresponding to the 0,0 corners (and the three lowest patch index corners whenever valence > 3)
                    // We use _getLowestIndices such that the corners are assigned to increasing patch corners
                    std::vector<std::pair<index_t,index_t>> indices0  = _getAllInterfaceIndices(pcorner,0,m_Bbases);
                    std::vector<std::pair<index_t,index_t>> indices1  = _getAllInterfaceIndices(pcorner,1,m_Bbases);
                    _getLowestIndices(indices1,1);
                    indices0.push_back(indices1[0]);

                    std::vector<index_t> rowIndices;
                    for (std::vector<std::pair<index_t,index_t>>::iterator it=indices0.begin(); it!=indices0.end(); it++)
                    {
                        GISMO_ASSERT(m_mapModified.is_free(it->second,it->first),"This DoF must be free! patch = "<<it->first<<"; index = "<<it->first);
                        rowIndices.push_back(m_mapModified.index(it->second,it->first));
                    }

                    index_t rowIdx;
                    for (index_t j=0; j!=Cg.cols(); j++)
                    {
                        rowIdx = rowIndices[j];
                        coefs.row(rowIdx) = Cg.col(j).transpose();
                    }

                    for (size_t k = 0; k!=pcorners.size(); k++)
                        m_vertCheck[ _vertIndex(pcorners[k].patch, pcorners[k].corner()) ] = true;
                }
                else
                {
                    m_vertCheck[ _vertIndex(pcorner.patch, pcorner.corner()) ] = true;
                    continue;
                }
            }
        }
        return coefs;
    }

    template<short_t d,class T>
    void gsAlmostC1<d,T>::setCoefficients(const gsMatrix<T> & coefs, gsMultiPatch<T> & mp) const
    {
        std::vector<index_t> sizes(mp.nPatches());
        index_t totalsize = 0;
        for (size_t p=0; p!=mp.nPatches(); p++) // patches
        {
            sizes.at(p) = mp.patch(p).coefs().rows();
            totalsize += sizes.at(p);
        }

        GISMO_ASSERT(totalsize==coefs.rows(),"Sizes do not agree");

        gsMultiBasis<T> basis(mp);
        gsDofMapper tmpMap(basis);
        tmpMap.finalize();

        index_t offset = 0;
        for (size_t p=0; p!=mp.nPatches(); p++) // patches
        {
            for (index_t k=0; k!=sizes.at(p); k++)
            {
                mp.patch(p).coefs().row(k) = coefs.row(tmpMap.index(k,p));
            }
            offset += sizes.at(p);
        }

    }

    /*=====================================================================================
                                    Construction functions
    =====================================================================================*/

    template<short_t d,class T>
    void gsAlmostC1<d,T>::_initBasis()
    {
        m_bases = gsMultiBasis<T>(m_patches);
    }

    template<short_t d,class T>
    void gsAlmostC1<d,T>::_initTHB()
    {
        m_RefPatches = m_patches;
        // Cast all patches of the mp object to THB splines
        gsTHBSpline<d,T> thb;
        gsTensorBSpline<d,T> * geo;
        for (size_t k=0; k!=m_RefPatches.nPatches(); ++k)
        {
            if ( (geo = dynamic_cast< gsTensorBSpline<d,T> * > (&m_RefPatches.patch(k))) )
            {
                thb = gsTHBSpline<d,T>(*geo);
                m_RefPatches.patch(k) = thb;
            }
            else if (dynamic_cast< gsTHBSpline<d,T> * > (&m_RefPatches.patch(k)))
            { }
            else
                gsWarn<<"No THB basis was constructed";
        }
    }


    template <short_t d, class T>
    void gsAlmostC1<d,T>::_refBoxes(std::vector<std::vector<index_t>> & patchBoxes)
    {
        patchBoxes.clear();
        patchBoxes.resize(m_RefPatches.nPatches());

        // prepare the geometry
        gsMatrix<index_t> box(d,2);
        std::vector<index_t> boxes;
        gsVector<bool> pars;
        index_t nelements;
        patchCorner corner;
        std::vector<patchCorner> cornerList;
        std::vector<std::vector<patchCorner> > cornerLists = _getSpecialCornerLists(m_RefPatches);

        index_t N = cornerLists.size();

        // Make a mask of corners per patch to track which ones have been handled
        gsMatrix<bool> mask(m_RefPatches.nPatches(),math::pow(2,d));
        mask.setConstant(false);

        for (index_t v =0; v<N; v++)
        {// Loop over EVs
            for (size_t c = 0; c<cornerLists[v].size(); c++)
            {// Loop over corners per EV
                corner = cornerLists[v].at(c);
                gsHTensorBasis<d,T> * basis = dynamic_cast<gsHTensorBasis<d,T>*>(&m_RefPatches.basis(corner.patch));

                if (mask(corner.patch,corner.corner()-1))
                    continue;

                corner.parameters_into(d,pars);
                box.setZero();
                for (short_t dim = 0; dim!=d; dim++)
                {
                    const gsKnotVector<T> & KV = basis->tensorLevel(0).knots(dim);
                    nelements = 1;
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
    void gsAlmostC1<d,T>::_makeTHB()
    {
        // prepare the geometry
        std::vector<std::vector<patchCorner> > cornerLists = _getSpecialCornerLists(m_RefPatches);

        if (cornerLists.size()!=0)
        {
            /// Change the coefficients
            gsMatrix<T> coefs = this->freeCoefficients(); // gets coefficients of the modified size
            coefs = m_matrix.transpose() * coefs; // maps to local size

            this->setCoefficients(coefs,m_RefPatches);

            /// Handle the EVs
            std::vector< std::vector<index_t> > elVec;
            this->_refBoxes(elVec);

            gsSparseMatrix<T> tmp;
            index_t rows = 0, cols = 0;
            std::vector<gsEigen::Triplet<T,index_t>> tripletList;
            for (size_t p=0; p!=m_RefPatches.nPatches(); p++)
            {
                // Transform using gsAsMatrix
                gsAsMatrix<index_t> boxMat(elVec[p],2*d+1,elVec[p].size()/(2*d+1));
                boxMat.row(0).array() += 1;
                boxMat.block(1,0,boxMat.rows()-1,boxMat.cols()).array() *= 2;

                gsHTensorBasis<2,T> *basis = dynamic_cast<gsHTensorBasis<2,T>*>(&m_RefPatches.basis(p));
                std::vector< gsSortedVector< index_t > > xmat = basis->getXmatrix();

                m_RefPatches.patch(p).refineElements(elVec[p]);

                basis->transfer(xmat,tmp);
                for (index_t i = 0; i<tmp.outerSize(); ++i)
                    for (typename gsSparseMatrix<T>::iterator it(tmp,i); it; ++it)
                        tripletList.push_back(gsEigen::Triplet<T,index_t>(it.row()+rows,it.col()+cols,it.value()));

                rows += tmp.rows();
                cols += tmp.cols();
            }

            m_tMatrix.resize(rows,cols);
            m_tMatrix.setFromTriplets(tripletList.begin(), tripletList.end());

            m_tMatrix.makeCompressed();
            m_bases = gsMultiBasis<T>(m_RefPatches);
        }

        // redefine the mappers
        m_mapOriginal = gsDofMapper(m_bases);
        m_mapOriginal.finalize();

        // gsWriteParaview<T>(m_RefPatches,"mp_ref",100,true);
    }

    template<short_t d, class T>
    std::vector<std::vector<patchCorner> > gsAlmostC1<d,T>::_getSpecialCornerLists(const gsMultiPatch<T> & patches)
    {
        std::vector<std::vector<patchCorner> > cornerLists;
        // get the corners that need refinement
        std::vector<patchCorner> cornerList;
        patchCorner pcorner;
        index_t cidx;
        for(size_t p = 0;p<patches.nPatches();++p)
        {
            for(index_t c=1;c<5;++c)
            {
                pcorner=patchCorner(p,c);
                cidx = _vertIndex(p,c);
                bool C0 = m_C0s[cidx];
                bool isCycle = patches.getCornerList(pcorner,cornerList);
                bool alreadyReached = false;
                for(size_t k = 0;k<cornerList.size();++k)
                    if((size_t)cornerList[k].patch<p)
                        alreadyReached = true;

                // add if
                // interior vertex with valence!=4
                // or
                // boundary vertex with valence > 2 (unless C0, then valence > 1)
                if(((isCycle&&cornerList.size()!=4)||((!isCycle)&&cornerList.size()>2-(size_t)C0))&&!alreadyReached)
                    cornerLists.push_back(cornerList);
            }
        }
        return cornerLists;
    }

    template<short_t d,class T>
    void gsAlmostC1<d,T>::_computeEVs()
    {
        /*
            Our goal is to create three vectors c11, c12, c21 which all contain the
            c11, c12 and c21 coefficients of the patches around the EV in the right order
            (counter)-clockwise.
        */

        std::vector<std::vector<patchCorner> > cornerLists = _getSpecialCornerLists(m_RefPatches);


        std::fill(m_vertCheck.begin(), m_vertCheck.end(), false);
        if (cornerLists.size()!=0)
        {
            m_matrix = m_matrix * m_tMatrix.transpose();

            std::vector<patchCorner> pcorners;
            patchCorner pcorner;
            gsMatrix<T> Cg;         // coefficients
            gsMatrix<T> ub;         // baricentric coordinates
            gsMatrix<index_t> uind; // column indices of baricentric coordinates
            index_t cidx;

            for (std::vector<std::vector<patchCorner> >::iterator it=cornerLists.begin(); it!=cornerLists.end(); it++)
            {

                std::vector<patchCorner> pcorners = *it;
                pcorner = it->at(0);
                cidx = _vertIndex(pcorner.patch,pcorner.corner());
                if (m_vertCheck.at(cidx))
                    continue;

                std::pair<index_t,bool> vdata = _vertexData(pcorner); // corner c

                // get the triangle
                gsMatrix<T> Cg;
                std::tie(Cg,ub,uind) = _makeTriangle(pcorner);

                // The barycentric coordinates will be attached to the matrix rows corresponding to the 0,0 corners (and the three lowest patch index corners whenever valence > 3)
                // We use _getLowestIndices such that the corners are assigned to increasing patch corners
                // We need the index on the old basis (the unrefined basis), because we plug it into the mapModified (which maps the local DoFs to the global ones)
                std::vector<std::pair<index_t,index_t>> indices, tmp;
                if (vdata.first==2)
                {
                    // These are two indices
                    indices  = _getAllInterfaceIndices(pcorner,0,m_Bbases);
                    tmp      = _getAllInterfaceIndices(pcorner,1,m_Bbases);
                    _getLowestIndices(tmp,1);
                    indices.push_back(tmp[0]);
                }
                else
                {
                    indices  = _getAllInterfaceIndices(pcorner,0,m_Bbases);
                    _getLowestIndices(indices,3);
                }


                std::vector<index_t> rowIndices;
                rowIndices.reserve(3);
                for (std::vector<std::pair<index_t,index_t>>::iterator it=indices.begin(); it!=indices.end(); it++)
                {
                    // We need the index on the old basis (the unrefined basis), because we plug it into the mapModified (which maps the local DoFs to the global ones)
                    GISMO_ASSERT(m_mapModified.is_free(it->second,it->first),"This DoF must be free! patch = "<<it->first<<"; index = "<<it->first);
                    rowIndices.push_back(m_mapModified.index(it->second,it->first));
                }

                index_t rowIdx,colIdx;
                // set the colums related to the barycentric columns equal to zero
                for (index_t j=0; j!=ub.cols(); j++)
                {
                    colIdx = uind(0,j);
                    m_matrix.prune(
                                    [&colIdx](index_t i, index_t j, T)
                                    { return j!=colIdx; }
                                    );
                }

                for (index_t i=0; i!=ub.rows(); i++)
                    for (index_t j=0; j!=ub.cols(); j++)
                    {
                        rowIdx = rowIndices[i];
                        colIdx = uind(0,j);
                        m_matrix(rowIdx,colIdx) = ub(i,j);
                    }

                for (size_t k = 0; k!=pcorners.size(); k++)
                    m_vertCheck[ _vertIndex(pcorners[k].patch, pcorners[k].corner()) ] = true;
            }
            m_matrix.makeCompressed();
        }
    }

    template<short_t d,class T>
    void gsAlmostC1<d,T>::_countDoFs() // also initialize the mappers!
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

        // add DoFs for the vertices (denoted by v) if
        // - part of a boundary vertex with valence 1
        // - valence >2 (interior or boundary vertex) [add 3 in total]

        // vertices (denoted by v)
        tmp = 0;
        std::vector<bool> passed(m_bases.nBases()*4);
        std::fill(passed.begin(), passed.end(), false);

        std::vector<patchCorner> corners;
        for (size_t p=0; p!=m_bases.nBases(); p++)
            for (index_t c=1; c<5; c++)
            {
                index_t idx = _vertIndex(p,c);
                if (!passed.at(idx))
                {
                    m_topology.getCornerList(patchCorner(p,c),corners);
                    for (size_t k=0; k!=corners.size(); k++)
                        passed.at(_vertIndex(corners[k].patch,corners[k])) = true;

                    std::pair<index_t,bool> vdata = _vertexData(patchCorner(p,c)); // corner c
                    bool C0 = m_C0s[idx];
                    // 1,1; 0,0; 0,1; 1,0 DoFs
                    if ((!vdata.second) && vdata.first==1) // valence = 1, must be a boundary vertex
                        tmp += 4;

                    // both 1,1 DoFs + 2 for the boundary 1,0 or 0,1 DoFs
                    else if ((!vdata.second) && vdata.first==2 && !C0) // valence = 1, must be a boundary vertex
                        tmp += 4;

                    // all 1,1 DoFs + 3 for the triangle + 2 for the boundary 1,0 or 0,1 DoFs
                    else if ((!vdata.second) && vdata.first>2 && !C0)
                        tmp += vdata.first+3+2;

                    // all 1,1 DoFs + 0,0 DoFs + 2 for the boundary 1,0 or 0,1 DoFs + 1 for the triangle
                    else if ((!vdata.second) && vdata.first==2 && C0)
                        tmp += 2*vdata.first+2+1;

                    // all 1,1 DoFs + 0,0 DoFs + 2 for the boundary 1,0 or 0,1 DoFs
                    else if ((!vdata.second) && vdata.first>2 && C0)
                        tmp += 2*vdata.first+2;

                    // all 1,1 DoFs
                    else if (( vdata.second) && vdata.first==4) // valence = 1, must be a boundary vertex
                        tmp += 4;

                    // all 1,1 DoFs + 3 for the triangle
                    else
                        tmp += vdata.first+3;
                }
            }

        // gsDebug<<"Number of vertex DoFs: "<<tmp<<"\n";
        m_size += tmp;
    }

    template<short_t d,class T>
    void gsAlmostC1<d,T>::_computeVertexMapper(patchCorner pcorner)
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
        else if (interior && valence==4)
            _computeMapperInteriorVertex_v4(pcorner,valence);
        else if (interior && (valence==3 || valence> 4) )
            _computeMapperInteriorVertex_v(pcorner,valence);
        else
            GISMO_ERROR("Something went terribly wrong, interior="<<interior<<"; valence="<<valence);

        // label vertex as processed
        m_vertCheck[ cidx ] = true;
    }

    template<short_t d,class T>
    void gsAlmostC1<d,T>::_computeMapperRegularBoundaryVertexNonSmooth_v2(patchCorner pcorner, index_t valence)
    {
        // for C0 vertices, the 1,0 and 0,1 DoFs on the interface need to be eliminated
        // However, we need a minimum of 3 DoFs around the vertex and we keep the 0,0s anyways
        // Using _removeLowestIndices(indices,3), only the 0,1 or 1,0 index with the highest is selected for removal
        // The 1,0 or 0,1s on the boundary are also kept

        // The 0,0s are kept
        // std::vector<std::pair<index_t,index_t>> indices0  = _getAllInterfaceIndices(pcorner,0,m_bases);
        // From the 1,0s we take the highest and eliminate it
        std::vector<std::pair<index_t,index_t>> indices1  = _getAllInterfaceIndices(pcorner,1,m_bases);

        _removeLowestIndices(indices1,1);
        for (std::vector<std::pair<index_t,index_t>>::iterator it=indices1.begin(); it!=indices1.end(); it++)
            m_mapModified.eliminateDof(it->second,it->first);

        std::vector<patchCorner> pcorners;
        m_topology.getCornerList(pcorner,pcorners);
        for (std::vector<patchCorner>::iterator it=pcorners.begin(); it!=pcorners.end(); it++)
        {
            // mark the vertex as passed
            m_vertCheck[ _vertIndex(it->patch, it->corner()) ] = true;
        }
    }

    // Reimplemented
    template<short_t d,class T>
    void gsAlmostC1<d,T>::_computeMapperIrregularBoundaryVertexSmooth_v(patchCorner pcorner, index_t valence)
    {
        std::vector<std::pair<index_t,index_t>> indices0  = _getAllInterfaceIndices(pcorner,0,m_bases);
        std::vector<std::pair<index_t,index_t>> indices1  = _getAllInterfaceIndices(pcorner,1,m_bases);
        std::vector<patchCorner> pcorners;
        m_topology.getCornerList(pcorner,pcorners);
        for (std::vector<patchCorner>::iterator it=pcorners.begin(); it!=pcorners.end(); it++)
        {
            // mark the vertex as passed
            m_vertCheck[ _vertIndex(it->patch, it->corner()) ] = true;
        }

        // Eliminate the 1,0 and 0,1s
        for (std::vector<std::pair<index_t,index_t>>::iterator it=indices1.begin(); it!=indices1.end(); it++)
            m_mapModified.eliminateDof(it->second,it->first);

        _removeLowestIndices(indices0,3);
        for (std::vector<std::pair<index_t,index_t>>::iterator it=indices0.begin(); it!=indices0.end(); it++)
            m_mapModified.eliminateDof(it->second,it->first);
    }

    // Extra compared to bass class
    template<short_t d,class T>
    void gsAlmostC1<d,T>::_computeMapperInteriorVertex_v4(patchCorner pcorner, index_t valence)
    {
        this->_computeMapperRegularBoundaryVertexSmooth_v2(pcorner,valence);
    }

    // Reimplemented
    template<short_t d,class T>
    void gsAlmostC1<d,T>::_computeMapperInteriorVertex_v(patchCorner pcorner, index_t valence)
    {
        std::vector<std::pair<index_t,index_t>> indices0  = _getAllInterfaceIndices(pcorner,0,m_bases);
        std::vector<std::pair<index_t,index_t>> indices1  = _getAllInterfaceIndices(pcorner,1,m_bases);
        std::vector<patchCorner> pcorners;
        m_topology.getCornerList(pcorner,pcorners);
        for (std::vector<patchCorner>::iterator it=pcorners.begin(); it!=pcorners.end(); it++)
        {
            // mark the vertex as passed
            m_vertCheck[ _vertIndex(it->patch, it->corner()) ] = true;
        }

        // Eliminate the left-over 0,0s
        _removeLowestIndices(indices0,3);
        for (std::vector<std::pair<index_t,index_t>>::iterator it=indices0.begin(); it!=indices0.end(); it++)
            m_mapModified.eliminateDof(it->second,it->first);
        // ... and the 1,0 and 0,1s
        for (std::vector<std::pair<index_t,index_t>>::iterator it=indices1.begin(); it!=indices1.end(); it++)
            m_mapModified.eliminateDof(it->second,it->first);
    }

    template<short_t d,class T>
    void gsAlmostC1<d,T>::_handleRegularBoundaryVertexNonSmooth(patchCorner pcorner, index_t valence)
    {
        this->_handleIrregularBoundaryVertexNonSmooth(pcorner,valence);
    }

    template<short_t d,class T>
    void gsAlmostC1<d,T>::_handleIrregularBoundaryVertexSmooth(patchCorner pcorner, index_t valence)
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

        gsBasis<T> * basis;

        // 2. make container for the interfaces
        std::vector<index_t> rowIndices, colIndices, patchIndices;

        // pcorner is the current corner
        m_topology.getCornerList(pcorner,corners);

        index_t b11_p1 = _indexFromVert(1,pcorner,psides[0],1); // point 1,1 (does not matter which reference side is taken)
        rowIdx = m_mapModified.index(b11_p1,pcorner.patch);
        colIdx = m_mapOriginal.index(b11_p1,pcorner.patch);
        // Influence of 1,1 to itself
        weight = 1.;
        entries.push_back(std::make_tuple(rowIdx,colIdx,weight));

        for (std::vector<patchSide>::iterator side = psides.begin(); side != psides.end(); ++side)
        {
            if (!m_topology.isInterface(*side))
                continue;

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

        // The 1,1 corners of each patch will be given 0.5 weight in the interface handling, but in addition they will receive a 1/(v+2) weight from the 0,0 DoFs on each patch
        pcorner.getContainingSides(d,psides);

        // colIndices stores the 0,0 corners (including the 0,0s of the boundary sides)
        for (typename std::vector<patchCorner>::iterator it = corners.begin(); it!=corners.end(); it++)
        {
            basis = &m_bases.basis(it->patch);
            colIndices.push_back(basis->functionAtCorner(*it));
            patchIndices.push_back(it->patch);
        }

        basis = &m_bases.basis(pcorner.patch);
        // Check if one of the adjacent interfaces is a boundary; if so, add weight 1.0 to itself and add it to the rowIndices
        index_t idx;

        for (index_t k = 0; k!=2; k++)
            if (!m_topology.getInterface(psides[k],iface)) // check if the side is NOT an interface
            {
                idx = _indexFromVert(1,pcorner,psides[k],0);
                rowIdx = m_mapModified.index(idx,pcorner.patch); //1,0 corner (on the boundary)
                colIdx = m_mapOriginal.index(idx,pcorner.patch); //1,0 corner (on the boundary)
                weight = 1.0;
                entries.push_back(std::make_tuple(rowIdx,colIdx,weight));
                rowIndices.push_back(rowIdx);
            }

        GISMO_ASSERT(rowIndices.size()<2,"Size issue, the boundary vertex is adjacent to two boundaries??" << rowIndices.size());

        if (rowIndices.size()==1)
        {
            rowIdx = rowIndices[0];
            for (size_t k=0; k!=colIndices.size(); k++)
            {
                colIdx = m_mapOriginal.index(colIndices.at(k),patchIndices.at(k));
                weight = 1. / 2.;
                entries.push_back(std::make_tuple(rowIdx,colIdx,weight));
            }
        }

        #pragma omp critical (handle_boundary_vertex_ff)
        {
            _pushAndCheck(entries);

            // Furthermore, if the corner is one of the three DoFs that is preserved, we mark the 0,0;0,1;1,0 DoF as handled (should be a zero-row)
            std::vector<std::pair<index_t,index_t>> indices = this->_getInterfaceIndices(pcorner,0,m_bases);
            for (std::vector<std::pair<index_t,index_t>>::iterator it=indices.begin(); it!=indices.end(); it++)
                if (m_mapModified.is_free(it->second,it->first))
                {
                    rowIdx = m_mapModified.index(it->second,it->first);
                    m_basisCheck[rowIdx] = true;
                }

            m_basisCheck[rowIdx] = true;
            m_vertCheck[ _vertIndex(pcorner.patch, pcorner.corner()) ] = true;
        }
    }

    template<short_t d,class T>
    void gsAlmostC1<d,T>::_handleIrregularBoundaryVertexNonSmooth(patchCorner pcorner, index_t valence)
    {
        Base::_handleIrregularBoundaryVertexNonSmooth(pcorner,valence);
        #pragma omp critical (handle_boundary_vertex_ff)
        {
            index_t rowIdx;
            // Furthermore, if the corner is one of the three DoFs that is preserved, we mark the 0,0;0,1;1,0 DoF as handled (should be a zero-row)
            std::vector<std::pair<index_t,index_t>> indices = this->_getInterfaceIndices(pcorner,1,m_bases);
            for (std::vector<std::pair<index_t,index_t>>::iterator it=indices.begin(); it!=indices.end(); it++)
                if (m_mapModified.is_free(it->second,it->first))
                {
                    rowIdx = m_mapModified.index(it->second,it->first);
                    m_basisCheck[rowIdx] = true;
                }
        }
    }

    template<short_t d,class T>
    void gsAlmostC1<d,T>::_handleInteriorVertex(patchCorner pcorner, index_t valence)
    {
        Base::_handleInteriorVertex(pcorner,valence);
        #pragma omp critical (handle_interior_vertex)
        {
            index_t rowIdx;
            // Furthermore, if the corner is one of the three DoFs that is preserved, we mark the 0,0;0,1;1,0 DoF as handled (should be a zero-row)
            std::vector<std::pair<index_t,index_t>> indices = this->_getInterfaceIndices(pcorner,0,m_bases);
            for (std::vector<std::pair<index_t,index_t>>::iterator it=indices.begin(); it!=indices.end(); it++)
                if (m_mapModified.is_free(it->second,it->first))
                {
                    rowIdx = m_mapModified.index(it->second,it->first);
                    m_basisCheck[rowIdx] = true;
                }
        }
    }


} // namespace gismo
