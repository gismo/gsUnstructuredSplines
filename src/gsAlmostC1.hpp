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

        for (size_t p=0; p!=m_patches.nPatches(); p++)
            for (short_t dim=0; dim!=d; dim++)
                GISMO_ENSURE(m_patches.basis(p).degree(dim)==2,"Degree of the basis ( dimension "<<dim<<" ) of patch "<<p<<" is "<<m_patches.basis(p).degree(dim)<<", but should be 2!");

        compute();
    }

    template<short_t d,class T>
    gsAlmostC1<d,T>::~gsAlmostC1()
    {
        freeAll(m_bases);
    }

    template<short_t d,class T>
    void gsAlmostC1<d,T>::defaultOptions()
    {
        m_options.addSwitch("SharpCorners","Reproduce sharp corners",true);
        m_options.addReal("SharpCornerTolerance","Sharp corner tolerance",1e-2);
        m_options.addSwitch("Verbose","Verbose output",false);
    }


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
        typename gsExprEvaluator<T>::geometryMap Gm = ev.getMap(m_patches);
        index_t k = 0;
        for (typename std::vector<patchCorner>::const_iterator it = corners.begin(); it!=corners.end(); it++, k++)
        {
            it->corner().parameters_into(m_patches.parDim(),pars); // get the parametric coordinates of the corner
            mat = pars.template cast<T>(); // cast to real coordinates
            normals.col(k) = ev.eval(sn(Gm).normalized(),mat,it->patch);
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
        mat = pars.template cast<T>(); // cast to real coordinates
        um.block(0,0,tdim,1) = m_RefPatches.patch(corner.patch).eval(mat);
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
            GISMO_ASSERT((Cp * ub.col(k)-up.col(k).head(2)).norm()<1e-14,"Something went wrong with the computation of the barycentric coordinates");
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
            Cg.conservativeResize(2,Eigen::NoChange);

        return std::make_tuple(Cg,ub,uind);
    }

    template<short_t d,class T>
    gsMatrix<T,3,3> gsAlmostC1<d,T>::_getRotationMatrix(const gsVector<T,3> & a, const gsVector<T,3> & b) const
    {
        GISMO_ASSERT(std::abs(a.norm()-1)<1e-14,"A must be a unit vector, a.norm() = "<<std::abs(a.norm()-1));
        GISMO_ASSERT(std::abs(b.norm()-1)<1e-14,"A must be a unit vector, b.norm() = "<<std::abs(b.norm()-1));

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
        for (size_t p=0; p!=m_patches.nPatches(); p++) // patches
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
        for (size_t p=0; p!=m_patches.nPatches(); p++)
        {
            for (index_t c=1; c<5; c++)
            {
                cidx = _vertIndex(p,c);
                if (m_vertCheck.at(cidx))
                    continue;

                bool C0 = m_C0s[cidx];
                pcorner = patchCorner(p,c);
                m_patches.getCornerList(pcorner,pcorners);
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
    void gsAlmostC1<d,T>::_makeTHB()
    {
        m_RefPatches = m_patches;
        // prepare the geometry
        std::vector<std::vector<patchCorner> > cornerLists = _getSpecialCornerLists(m_RefPatches);

        if (cornerLists.size()!=0)
        {
            /// Change the coefficients
            gsMatrix<T> coefs = this->freeCoefficients(); // gets coefficients of the modified size
            coefs = m_matrix.transpose() * coefs; // maps to local size

            this->setCoefficients(coefs,m_RefPatches);

            /// Handle the EVs
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
                    std::vector<index_t> elements = basis->asElements(boxes,0); // 0-ring

                    elVec.at(corner.patch).insert(elVec.at(corner.patch).end(), elements.begin(), elements.end());

                    // gsHTensorBasis<2,T> *basis = dynamic_cast<gsHTensorBasis<2,T>*>(&m_RefPatches.basis(corner.patch));

                    // basis->refineElements(elements, m_tMatrix);
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
        m_mapOriginal = gsDofMapper(m_bases);
        m_mapOriginal.finalize();

        // gsWriteParaview<>(m_RefPatches,"mp_ref",1000,true);
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
            for(int c=1;c<5;++c)
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
    void gsAlmostC1<d,T>::_initialize() // also initialize the mappers!
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
        _initCoefs();
    }

    template<short_t d,class T>
    void gsAlmostC1<d,T>::_countDoFs() // also initialize the mappers!
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

        // add DoFs for the vertices (denoted by v) if
        // - part of a boundary vertex with valence 1
        // - valence >2 (interior or boundary vertex) [add 3 in total]

        // vertices (denoted by v)
        tmp = 0;
        std::vector<bool> passed(m_patches.nPatches()*4);
        std::fill(passed.begin(), passed.end(), false);

        std::vector<patchCorner> corners;
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
    void gsAlmostC1<d,T>::_computeInterfaceMapper(boundaryInterface iface)
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

            #pragma omp critical (side_interface)
            {
                m_mapModified.markBoundary(patches[p], indices);
                m_sideCheck.at(sidx) = true;
            }
        }
    }

    template<short_t d,class T>
    void gsAlmostC1<d,T>::_computeBoundaryMapper(patchSide boundary)
    {
        index_t sidx = _sideIndex(boundary.patch,boundary.side());
        m_sideCheck.at(sidx) = true;
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
            _computeVertexMapper_impl<true,1>(pcorner,valence);
        else if (!interior && valence==2 && C0)
            _computeVertexMapper_impl<true,2,false>(pcorner,valence);
        else if (!interior && valence==2 && !C0)
            _computeVertexMapper_impl<true,2,true>(pcorner,valence);
        else if (!interior && valence >2 && C0)
            _computeVertexMapper_impl<true,-1,false>(pcorner,valence);
        else if (!interior && valence >2 && !C0)
            _computeVertexMapper_impl<true,-1,true>(pcorner,valence);
        else if (interior && valence==4)
            _computeVertexMapper_impl<false,4>(pcorner,valence);
        else if (interior && (valence==3 || valence> 4) )
            _computeVertexMapper_impl<false,-1>(pcorner,valence);
        else
            GISMO_ERROR("Something went terribly wrong, interior="<<interior<<"; valence="<<valence);

        // label vertex as processed
        m_vertCheck[ cidx ] = true;
    }

    template<short_t d,class T>
    template<bool _boundary, index_t _v>
    typename std::enable_if<  _boundary && _v==1, void>::type
    gsAlmostC1<d,T>::_computeVertexMapper_impl(patchCorner pcorner, index_t valence)
    {
        // do nothing,
    }

    template<short_t d,class T>
    template<bool _boundary, index_t _v, bool _smooth>
    typename std::enable_if< _boundary && _v==2 && _smooth, void>::type
    gsAlmostC1<d,T>::_computeVertexMapper_impl(patchCorner pcorner, index_t valence)
    {
        std::vector<patchSide> psides(2);

        /*
        v = 2
            o o o @ X |i| X @ o o o                 x: eliminated DoFs by vertex rule
            o o o @ X |n| X @ o o o                 X: eliminated DoFs by interface/boundary rule
            o o o @ X |t| X @ o o o                 o: preserved DoFs (interior)
            o o o @ x |e| x @ o o o                 @: modified DoFs by interface rule
            o o o * x |r| x * o o o                 *: modified DoFs by vertex rule (unique DoFs)
            -----------------------                 %: modified DoFs by vertex rule (matched DoFs)
            -boundary-| | -boundary-
            -----------------------

        v = 4
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
            if (m_patches.isInterface(psides[p]))
                for (index_t k=0; k!=2; k++)
                    m_mapModified.eliminateDof(_indexFromVert(k,pcorner,psides[p],0,0),pcorner.patch);
        }

        // // we mark the nodes belonging to the interface
        // pcorner.getContainingSides(d,psides);
        // for (size_t p=0; p!=psides.size(); p++)
        // {
        //     if (m_patches.isInterface(psides[p]))
        //     {
        //         // the 0,0 vertex should be eliminated
        //         m_mapModified.eliminateDof(basis->functionAtCorner(pcorner),pcorner.patch);
        //     }
        // }
    }

    template<short_t d,class T>
    template<bool _boundary, index_t _v, bool _smooth>
    typename std::enable_if<  _boundary && _v==2 && (!_smooth), void>::type
    gsAlmostC1<d,T>::_computeVertexMapper_impl(patchCorner pcorner, index_t valence)
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
        m_patches.getCornerList(pcorner,pcorners);
        for (std::vector<patchCorner>::iterator it=pcorners.begin(); it!=pcorners.end(); it++)
        {
            // mark the vertex as passed
            m_vertCheck[ _vertIndex(it->patch, it->corner()) ] = true;
        }
    }

    template<short_t d,class T>
    template<bool _boundary, index_t _v, bool _smooth>
    typename std::enable_if<  _boundary && _v==-1 && _smooth, void>::type
    gsAlmostC1<d,T>::_computeVertexMapper_impl(patchCorner pcorner, index_t valence)
    {
        std::vector<std::pair<index_t,index_t>> indices0  = _getAllInterfaceIndices(pcorner,0,m_bases);
        std::vector<std::pair<index_t,index_t>> indices1  = _getAllInterfaceIndices(pcorner,1,m_bases);
        std::vector<patchCorner> pcorners;
        m_patches.getCornerList(pcorner,pcorners);
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

    template<short_t d,class T>
    template<bool _boundary, index_t _v, bool _smooth>
    typename std::enable_if<  _boundary && _v==-1 && (!_smooth), void>::type
    gsAlmostC1<d,T>::_computeVertexMapper_impl(patchCorner pcorner, index_t valence)
    {
        // Get all the 0,1 or 1,0 DoFs on the interfaces
        // The 0,0 DoFs are kept but deleted later from the matrix
        std::vector<patchSide> psides(2);
        pcorner.getContainingSides(d,psides);
        for (size_t p=0; p!=psides.size(); p++)
        {
            // the 0,k (k=0,1) DoF should be eliminated
            if (m_patches.isInterface(psides[p]))
                m_mapModified.eliminateDof(_indexFromVert(1,pcorner,psides[p],0,0),pcorner.patch);
        }
    }

    template<short_t d,class T>
    template<bool _boundary, index_t _v>
    typename std::enable_if<  (!_boundary) && _v==4, void>::type
    gsAlmostC1<d,T>::_computeVertexMapper_impl(patchCorner pcorner, index_t valence)
    {
        this->_computeVertexMapper_impl<true,2,true>(pcorner,valence);
    }

    template<short_t d,class T>
    template<bool _boundary, index_t _v>
    typename std::enable_if<  (!_boundary) && _v==-1, void>::type
    gsAlmostC1<d,T>::_computeVertexMapper_impl(patchCorner pcorner, index_t valence)
    {
        std::vector<std::pair<index_t,index_t>> indices0  = _getAllInterfaceIndices(pcorner,0,m_bases);
        std::vector<std::pair<index_t,index_t>> indices1  = _getAllInterfaceIndices(pcorner,1,m_bases);
        std::vector<patchCorner> pcorners;
        m_patches.getCornerList(pcorner,pcorners);
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
    void gsAlmostC1<d,T>::_computeMapper() // also initialize the mappers!
    {
        // interfaces
        _resetChecks(false);

        std::vector<index_t> patches(2);
        std::vector<patchSide> psides(2);
        patchCorner pcorner;
        std::pair<index_t,bool> vdata1, vdata2, vdata;

        // For the interfaces, we eliminate all DoFs located on the interface, except the ones coinciding with the end vertices
#pragma omp parallel
{
        #pragma omp parallel for
        for(gsBoxTopology::const_iiterator iit = m_patches.iBegin(); iit!= m_patches.iEnd(); iit++)
            _computeInterfaceMapper(*iit);
}
        // On the boundaries, we don't do anything
#pragma omp parallel
{
        #pragma omp parallel for
        for(gsBoxTopology::const_biterator bit = m_patches.bBegin(); bit!= m_patches.bEnd(); bit++)
            _computeBoundaryMapper(*bit);
}

// m_mapModified.finalize();
// m_mapOriginal.finalize();


        // For the vertices, we eliminate as follows (where v is the valence):
        // - No elimination when v==1
        // - One on each side when v==2
        // - All but three when v>2
        for (size_t p=0; p!=m_patches.nPatches(); p++)
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

    template<short_t d,class T>
    void gsAlmostC1<d,T>::_handleRegularCorner(patchCorner pcorner)
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
    template<bool _regular, bool _smooth>
    typename std::enable_if<  _regular  &&   _smooth    , void>::type
    gsAlmostC1<d,T>::_handleBoundaryVertex(patchCorner pcorner, index_t valence)
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
        index_t iindex = m_patches.isInterface(psides[0]) ? 0 : 1;

        GISMO_ENSURE(m_patches.getInterface(psides[iindex],iface),"Must be an interface");

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
    template<bool _regular, bool _smooth>
    typename std::enable_if<  _regular  && (!_smooth)  , void>::type
    gsAlmostC1<d,T>::_handleBoundaryVertex(patchCorner pcorner, index_t valence)
    {
        this->_handleBoundaryVertex<false,false>(pcorner,valence);
    }

    template<short_t d,class T>
    template<bool _regular, bool _smooth>
    typename std::enable_if<(!_regular) &&   _smooth    , void>::type
    gsAlmostC1<d,T>::_handleBoundaryVertex(patchCorner pcorner, index_t valence)
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
        m_patches.getCornerList(pcorner,corners);

        index_t b11_p1 = _indexFromVert(1,pcorner,psides[0],1); // point 1,1 (does not matter which reference side is taken)
        rowIdx = m_mapModified.index(b11_p1,pcorner.patch);
        colIdx = m_mapOriginal.index(b11_p1,pcorner.patch);
        // Influence of 1,1 to itself
        weight = 1.;
        entries.push_back(std::make_tuple(rowIdx,colIdx,weight));

        for (std::vector<patchSide>::iterator side = psides.begin(); side != psides.end(); ++side)
        {
            if (!m_patches.isInterface(*side))
                continue;

            GISMO_ENSURE(m_patches.getInterface(*side,iface),"Side must be an interface!");

            m_patches.getNeighbour(*side,otherSide);
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
            if (!m_patches.getInterface(psides[k],iface)) // check if the side is NOT an interface
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
    template<bool _regular, bool _smooth>
    typename std::enable_if<(!_regular) && (!_smooth)   , void>::type
    gsAlmostC1<d,T>::_handleBoundaryVertex(patchCorner pcorner, index_t valence)
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
        m_patches.getCornerList(pcorner,corners);

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
            if (!m_patches.getInterface(psides[k],iface)) // check if the side is NOT an interface
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

            // Furthermore, if the corner is one of the three DoFs that is preserved, we mark the 0,0;0,1;1,0 DoF as handled (should be a zero-row)
            std::vector<std::pair<index_t,index_t>> indices = this->_getInterfaceIndices(pcorner,1,m_bases);
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
    void gsAlmostC1<d,T>::_handleInteriorVertex(patchCorner pcorner, index_t valence)
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

        m_patches.getCornerList(pcorner,corners);

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



        // Influence of 1,1 to itself
        for (std::vector<patchSide>::iterator side = psides.begin(); side != psides.end(); ++side)
        {
            GISMO_ENSURE(m_patches.getInterface(*side,iface),"Side must be an interface!");
            m_patches.getNeighbour(*side,otherSide);
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

            // Furthermore, if the corner is one of the three DoFs that is preserved, we mark the 0,0;0,1;1,0 DoF as handled (should be a zero-row)
            std::vector<std::pair<index_t,index_t>> indices = this->_getInterfaceIndices(pcorner,0,m_bases);
            for (std::vector<std::pair<index_t,index_t>>::iterator it=indices.begin(); it!=indices.end(); it++)
                if (m_mapModified.is_free(it->second,it->first))
                {
                    rowIdx = m_mapModified.index(it->second,it->first);
                    m_basisCheck[rowIdx] = true;
                }

            // Mark the vertex as processed
            m_vertCheck[ _vertIndex(pcorner.patch, pcorner.corner()) ] = true;
        }
    }

    template<short_t d,class T>
    void gsAlmostC1<d,T>::_handleVertex(patchCorner pcorner)
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
                _handleBoundaryVertex<true,true>(pcorner,vdata.first);
            else if (vdata.first==2 &&  C0)
                _handleBoundaryVertex<true,false>(pcorner,vdata.first);
            else if (vdata.first> 2 && !C0)
                _handleBoundaryVertex<false,true>(pcorner,vdata.first);
            else if (vdata.first> 2 &&  C0)
                _handleBoundaryVertex<false,false>(pcorner,vdata.first);
            else
                GISMO_ERROR("Something went wrong");
        } // end boundary vertices
        else // interior vertices
            _handleInteriorVertex(pcorner,vdata.first);
    }

    template<short_t d,class T>
    void gsAlmostC1<d,T>::_handleInterface(boundaryInterface iface)
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
        T weight;
        // loop over adjacent patches and couple the DoFs.
        for (index_t p =0; p!= 2; p++)
        {
            np = 1-p; // not index p;
            for (index_t k=0; k!= indices[p].size(); k++ )
            {
                rowIdx = m_mapModified.index(oindices[p].at(k),iface[p].patch);
                // rowIdx1 = m_mapOriginal.index(oindices[p].at(k),patches[p]);
                colIdx = m_mapOriginal.index(oindices[p].at(k),iface[p].patch);
                // m_matrix(rowIdx,colIdx) = 1.0;
                weight = 1.0;
                entries.push_back(std::make_tuple(rowIdx,colIdx,weight));

                colIdx = m_mapOriginal.index(indices[p].at(k),iface[p].patch);
                // m_matrix(rowIdx,colIdx) = 0.5;
                weight = 0.5;
                entries.push_back(std::make_tuple(rowIdx,colIdx,weight));

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
    void gsAlmostC1<d,T>::_handleBoundary(patchSide side)
    {
            std::vector<patchCorner> pcorners;
            std::vector<index_t> selectedIndices;
            gsBasis<T> * basis = &m_bases.basis(side.patch);
            sparseEntry_t entries;

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
    void gsAlmostC1<d,T>::_handleInterior()
    {
        #pragma omp critical (handle_interior)
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
    }

    template<short_t d,class T>
    void gsAlmostC1<d,T>::_computeSmoothMatrix()
    {
        GISMO_ASSERT((size_t)m_mapModified.freeSize()==m_size,"Size does not match predicted size, m_mapModified.freeSize()="<<m_mapModified.freeSize()<<"; m_size="<<m_size);

        _resetChecks(true);

        // iterate over the vertices
#pragma omp parallel
{
        #pragma omp parallel for collapse(2)
        for (size_t p=0; p!=m_patches.nPatches(); p++)
            for (index_t c=1; c<5; c++)
                _handleVertex(patchCorner(p,c));

        #pragma omp parallel for
        for(gsBoxTopology::const_iiterator iit = m_patches.iBegin(); iit!= m_patches.iEnd(); iit++)
            _handleInterface(*iit);

        // boundaries
        #pragma omp parallel for
        for(gsBoxTopology::const_biterator bit = m_patches.bBegin(); bit!= m_patches.bEnd(); bit++)
            _handleBoundary(*bit);

        _handleInterior();
}
        if (m_options.getSwitch("Verbose")) { _whichHandled(); }

        _performChecks(true);
    }
} // namespace gismo