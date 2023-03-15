/** @file gsDPatchBase.hpp

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
        std::fill(result.begin(), result.end(), false);
        index_t cidx;
        patchCorner pcorner;
        for (size_t p=0; p<m_patches.nPatches(); p++)
        {
            for (index_t c=1; c<5; c++)
            {
                cidx = _vertIndex(p,c);
                if (m_vertCheck.at(cidx))
                    continue;

                // Get the patchCorners which are on the boundary
                pcorner = patchCorner(p,c);
                std::vector<patchCorner> pcorners;
                m_patches.getCornerList(pcorner,pcorners);
                std::vector<std::pair<patchCorner,patchSide>> boundaries;
                for (std::vector<patchCorner>::iterator it=pcorners.begin(); it!=pcorners.end(); it++)
                {
                    std::vector<patchSide> psides;
                    it->getContainingSides(d,psides);
                    for (size_t s=0; s!=psides.size(); s++)
                    {
                        // the 0,k (k=0,1) DoF should be eliminated
                        if (!m_patches.isInterface(psides[s]))
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
                        gsVector<bool> pars;
                        it->first.parameters_into(m_patches.parDim(),pars); // get the parametric coordinates of the corner
                        md.points = pars.template cast<T>();
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


//     /*=====================================================================================
//                                     Special functions
//     =====================================================================================*/
//     template<short_t d,class T>
//     gsMatrix<T> gsDPatchBase<d,T>::_getNormals(const std::vector<patchCorner> & corners) const
//     {
//         gsMatrix<T> normals(3,corners.size());

//         gsVector<bool> pars;
//         gsMatrix<T> mat;

//         gsExprEvaluator<T> ev;
//         typename gsExprEvaluator<T>::geometryMap Gm = ev.getMap(m_patches);
//         index_t k = 0;
//         for (typename std::vector<patchCorner>::const_iterator it = corners.begin(); it!=corners.end(); it++, k++)
//         {
//             it->corner().parameters_into(m_patches.parDim(),pars); // get the parametric coordinates of the corner
//             mat = pars.template cast<T>(); // cast to real coordinates
//             normals.col(k) = ev.eval(sn(Gm).normalized(),mat,it->patch);
//         }
//         return normals;
//     }


//     template<short_t d,class T>
//     std::tuple<gsMatrix<T>,gsMatrix<T>,gsMatrix<index_t>> gsDPatchBase<d,T>::_makeTriangle(const patchCorner & corner) const
//     {
//         GISMO_ASSERT(m_RefPatches.nPatches()!=0,"Are the patches refined?");

//         index_t tdim = m_RefPatches.targetDim();

//         std::vector<patchCorner> corners;
//         m_RefPatches.getCornerList(corner,corners);

//         gsVector<bool> pars;
//         gsMatrix<T> mat;
//         // 1. Get the coordinates of the vertex and set its z coordinate to 0
//         gsMatrix<T> um(3,1), midpoint;
//         um.setZero();
//         corner.corner().parameters_into(m_RefPatches.parDim(),pars); // get the parametric coordinates of the corner
//         mat = pars.template cast<T>(); // cast to real coordinates
//         um.block(0,0,tdim,1) = m_RefPatches.patch(corner.patch).eval(mat);
//         midpoint = um; // store the original midpoint

//         // 2. Get the 0,0;0,1; 1,0; 1,1 coordinates
//         gsMatrix<T> u(3,corners.size()*4);
//         u.setZero();
//         gsMatrix<index_t> uind(1,corners.size()*4);
//         uind.setZero();

//         std::vector<patchSide> csides;
//         index_t idx;
//         for (size_t c = 0; c!=corners.size(); c++)
//         {
//             corners[c].getContainingSides(d,csides);
//             index_t k=0;
//             for (index_t i=0; i!=2; i++)
//                 for (index_t j=0; j!=2; j++,k++)
//                 {
//                     idx = _indexFromVert(i,corners[c],csides[0],j);
//                     uind(0,4*c+k) = m_mapOriginal.index(idx,corners[c].patch);
//                     u.block(0,4*c+k,m_RefPatches.targetDim(),1) = m_RefPatches.patch(corners[c].patch).coefs().row(idx).transpose();
//                 }
//         }

//         // 3. Translate all points to a coordinate system with origin um
//         gsMatrix<T> up = u;
//         for (index_t k=0; k!=up.cols(); k++)
//             up.col(k) -= um;

//         // 4. Rotate the points parallel the xy-plane and set their z-coordinates to 0
//         gsMatrix<T,3,3> Rn, Rx;
//         Rn.setIdentity();
//         if (m_RefPatches.targetDim()==2)
//         {
//             // do nothing
//         }
//         else if(m_RefPatches.targetDim()==3)
//         {
//             // Get the average normal at the corner
//             gsVector<T> avgnormal = _getNormals(corners).rowwise().mean();

//             // Find the rotation matrix that maps the average normal to the z axis
//             gsVector<T,3> ez;
//             ez<<0,0,1;
//             Rn = _getRotationMatrix(avgnormal.normalized(),ez);

//             for (index_t k=0; k!=up.cols(); k++)
//                 up.col(k).applyOnTheLeft(Rn);

//             up.row(2).setZero(); // all points
//             um.row(2).setZero();// midpoint
//         }
//         else
//             GISMO_ERROR("Target dimension of the multipatch should be 2 or 3, but is "<<m_RefPatches.targetDim());

//         // 5. Find the maximum distance from the midpoint to all points
//         T distance, maxDistance = 0;
//         gsMatrix<T> umax;
//         for (index_t k = 0; k!=up.cols(); k++)
//         {
//             distance = (up.col(k)).norm();
//             if (distance > maxDistance)
//             {
//                 maxDistance = distance;
//                 umax = up.col(k);
//             }
//         }

//         gsVector<T,3> ex;
//         ex<<1,0,0;

//         // 6. Rotate all points such that the maximum point is aligned with the x-axis
//         Rx = _getRotationMatrix(umax.normalized(),ex);
//         for (index_t k=0; k!=up.cols(); k++)
//             up.col(k).applyOnTheLeft(Rx);

//         // 7. Obtain the coordinates of the triangle that encloses the circle with radius maxDistance in the xy plane
//         T r = maxDistance;
//         T a = 1. / ( 1./6. * std::sqrt(3) ) * r;
//         T rr = 1. / 3. * std::sqrt(3) * a;

//         gsMatrix<T> Cp(2,3);
//         Cp.col(0)<<rr,0;
//         Cp.col(1)<<-r, 0.5*a;
//         Cp.col(2)<<-r,-0.5*a;

//         // 8. Get the barycentric coordinates of the points
//         gsMatrix<T> ub = up;
//         up.row(2).setOnes(); // project the parametric points to z=1
//         gsMatrix<T> A(3,3);
//         A.block(0,0,2,3) = Cp;
//         A.row(2).setOnes();

//         for (index_t k = 0; k!=ub.cols(); k++)
//         {
//             ub.col(k) = A.colPivHouseholderQr().solve(up.col(k));
//             GISMO_ASSERT((Cp * ub.col(k)-up.col(k).head(2)).norm()<1e-14,"Something went wrong with the computation of the barycentric coordinates");
//         }

//         // 9. Move the corners of the triangle back to physical coordinates
//         gsMatrix<T> Cg(3,3);
//         Cg.setZero();
//         Cg.block(0,0,2,3) = Cp;

//         for (index_t k = 0; k!=Cg.cols(); k++)
//         {
//             Cg.col(k).applyOnTheLeft((Rx).transpose());
//             Cg.col(k).applyOnTheLeft((Rn).transpose());
//             Cg.col(k) += midpoint;
//         }

//         if (m_RefPatches.targetDim()==2)
//             Cg.conservativeResize(2,Eigen::NoChange);

//         return std::make_tuple(Cg,ub,uind);
//     }

//     template<short_t d,class T>
//     void gsDPatchBase<d,T>::_toBarycentricCoordinates(const gsMatrix<T> & Cs, gsMatrix<T> & u) const
//     {

//     }

//     template<short_t d,class T>
//     gsMatrix<T,3,3> gsDPatchBase<d,T>::_getRotationMatrix(const gsVector<T,3> & a, const gsVector<T,3> & b) const
//     {
//         GISMO_ASSERT(std::abs(a.norm()-1)<1e-14,"A must be a unit vector, a.norm() = "<<std::abs(a.norm()-1));
//         GISMO_ASSERT(std::abs(b.norm()-1)<1e-14,"A must be a unit vector, b.norm() = "<<std::abs(b.norm()-1));

//         gsVector<T,3> v = a.cross(b);
//         v.normalize();
//         T theta = std::acos( a.dot(b) / ( a.norm() * b.norm() ) );

//         T s = std::sin(theta);
//         T c = std::cos(theta);
//         gsMatrix<T,3,3> R,vx,tmp, I;
//         R.setZero();
//         vx.setZero();

//         vx.row(0)<<0,-v.at(2),v.at(1);
//         vx.row(1)<<v.at(2),0,-v.at(0);
//         vx.row(2)<<-v.at(1),v.at(0),0;

//         I.setIdentity();
//         R += I*c;
//         R += vx * s;
//         tmp = (v*v.transpose()) * (1-c);
//         R += tmp;

//         GISMO_ASSERT((R * a - b).norm() < 1e-12,"Rotation matrix is wrong, R*a = "<<R*a<<"; b = "<<b);
//         return R;
//     }

    /*=====================================================================================
                                    Information functions
    =====================================================================================*/

    template<short_t d,class T>
    void gsDPatchBase<d,T>::mapperInfo() const
    {
        for (size_t p=0; p!=m_patches.nPatches(); p++)
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
        gsInfo<<"Patch "<<side.patch<<", side "<<side<<" is "<<(m_patches.isBoundary(side) ? "a boundary side" : "an interface")<<"\n";
    }

    template<short_t d,class T>
    void gsDPatchBase<d,T>::sideInfo() const
    {
        gsInfo<<"**D-Patch Side info**\n";
        for(size_t i = 0;i<m_patches.nPatches();++i)
            for(int j=1;j<=4;++j)
                sideInfo(patchSide(i,j));
    }

    template<short_t d,class T>
    void gsDPatchBase<d,T>::cornerInfo() const
    {
        gsInfo<<"**D-Patch Corner info**\n";
        for(size_t i = 0;i<m_patches.nPatches();++i)
            for(int j=1;j<=4;++j)
                vertexInfo(patchCorner(i,j));
    }

//     /*=====================================================================================
//                                     Coefficients
//     =====================================================================================*/
//     // ADD THE COEFFICIENTS OF THE TRIANGLES AS EXTRA COEFFICIENTS

//     template<short_t d,class T>
//     gsMatrix<T> gsDPatchBase<d,T>::freeCoefficients()
//     {
//         GISMO_ASSERT(m_mapModified.isFinalized(),"Mapper is not finalized, run XXXX first");

//         GISMO_ASSERT((size_t)m_mapModified.freeSize()==m_size,"Size does not match predicted size, m_mapModified.freeSize()="<<m_mapModified.freeSize()<<"; m_size="<<m_size);
//         gsMatrix<T> coefs(m_mapModified.freeSize(),m_patches.geoDim());

//         index_t size;
//         for (size_t p=0; p!=m_patches.nPatches(); p++) // patches
//         {
//             size = m_mapModified.patchSize(p);
//             for (index_t k=0; k!=size; k++)
//             {
//                 if (m_mapModified.is_free(k,p))
//                     coefs.row(m_mapModified.index(k,p,0)) = m_patches.patch(p).coefs().row(k);
//             }
//         }
//         return coefs;
//     }

//     template<short_t d,class T>
//     gsMatrix<T> gsDPatchBase<d,T>::_preCoefficients()
//     {
//         GISMO_ASSERT(m_mapModified.isFinalized(),"Mapper is not finalized, run XXXX first");

//         gsMatrix<T> coefs = this->freeCoefficients();

//         // Correct the EVs
//         std::fill(m_vertCheck.begin(), m_vertCheck.end(), false);
//         index_t cidx;
//         std::vector<patchCorner> pcorners;
//         patchCorner pcorner;
//         for (size_t p=0; p!=m_patches.nPatches(); p++)
//         {
//             for (index_t c=1; c<5; c++)
//             {
//                 cidx = _vertIndex(p,c);
//                 if (m_vertCheck.at(cidx))
//                     continue;

//                 bool C0 = m_C0s[cidx];
//                 pcorner = patchCorner(p,c);
//                 m_patches.getCornerList(pcorner,pcorners);
//                 std::pair<index_t,bool> vdata = _vertexData(pcorner); // corner c
//                 if (vdata.first > 2 && !(vdata.first==4 && vdata.second)) // valence must be 3 or larger, but it must not be an interior vertex with v=4
//                 {
//                     // get the triangle
//                     gsMatrix<T> Cg;
//                     std::tie(Cg,std::ignore,std::ignore) = _makeTriangle(pcorner);

//                     // The barycentric coordinates will be attached to the matrix rows corresponding to the 0,0 corners (and the three lowest patch index corners whenever valence > 3)
//                     // We use _getLowestIndices such that the corners are assigned to increasing patch corners
//                     std::vector<std::pair<index_t,index_t>> indices  = _getAllInterfaceIndices(pcorner,0,m_Bbases);
//                     _getLowestIndices(indices,3);

//                     std::vector<index_t> rowIndices;
//                     for (std::vector<std::pair<index_t,index_t>>::iterator it=indices.begin(); it!=indices.end(); it++)
//                     {
//                         GISMO_ASSERT(m_mapModified.is_free(it->second,it->first),"This DoF must be free! patch = "<<it->first<<"; index = "<<it->first);
//                         rowIndices.push_back(m_mapModified.index(it->second,it->first));
//                     }

//                     index_t rowIdx;
//                     for (index_t j=0; j!=Cg.cols(); j++)
//                     {
//                         rowIdx = rowIndices[j];
//                         coefs.row(rowIdx) = Cg.col(j).transpose();
//                     }

//                     for (size_t k = 0; k!=pcorners.size(); k++)
//                         m_vertCheck[ _vertIndex(pcorners[k].patch, pcorners[k].corner()) ] = true;
//                 }
//                 else if (vdata.first == 2 && C0) // valence must be 3 or larger, but it must not be an interior vertex with v=4
//                 {
//                     // get the triangle
//                     gsMatrix<T> Cg;
//                     std::tie(Cg,std::ignore,std::ignore) = _makeTriangle(pcorner);

//                     // The barycentric coordinates will be attached to the matrix rows corresponding to the 0,0 corners (and the three lowest patch index corners whenever valence > 3)
//                     // We use _getLowestIndices such that the corners are assigned to increasing patch corners
//                     std::vector<std::pair<index_t,index_t>> indices0  = _getAllInterfaceIndices(pcorner,0,m_Bbases);
//                     std::vector<std::pair<index_t,index_t>> indices1  = _getAllInterfaceIndices(pcorner,1,m_Bbases);
//                     _getLowestIndices(indices1,1);
//                     indices0.push_back(indices1[0]);

//                     std::vector<index_t> rowIndices;
//                     for (std::vector<std::pair<index_t,index_t>>::iterator it=indices0.begin(); it!=indices0.end(); it++)
//                     {
//                         GISMO_ASSERT(m_mapModified.is_free(it->second,it->first),"This DoF must be free! patch = "<<it->first<<"; index = "<<it->first);
//                         rowIndices.push_back(m_mapModified.index(it->second,it->first));
//                     }

//                     index_t rowIdx;
//                     for (index_t j=0; j!=Cg.cols(); j++)
//                     {
//                         rowIdx = rowIndices[j];
//                         coefs.row(rowIdx) = Cg.col(j).transpose();
//                     }

//                     for (size_t k = 0; k!=pcorners.size(); k++)
//                         m_vertCheck[ _vertIndex(pcorners[k].patch, pcorners[k].corner()) ] = true;
//                 }
//                 else
//                 {
//                     m_vertCheck[ _vertIndex(pcorner.patch, pcorner.corner()) ] = true;
//                     continue;
//                 }
//             }
//         }
//         return coefs;
//     }

    template<short_t d,class T>
    gsMatrix<T> gsDPatchBase<d,T>::allCoefficients() const
    {
        std::vector<index_t> sizes(m_patches.nPatches());
        index_t totalsize = 0;
        for (size_t p=0; p!=m_patches.nPatches(); p++) // patches
        {
            sizes.at(p) = m_patches.patch(p).coefs().rows();
            totalsize += sizes.at(p);
        }

        gsMultiBasis<T> basis(m_patches);
        gsDofMapper tmpMap(basis);
        tmpMap.finalize();

        gsMatrix<T> coefs(totalsize,m_patches.geoDim());
        index_t offset = 0;
        for (size_t p=0; p!=m_patches.nPatches(); p++) // patches
        {
            for (index_t k=0; k!=sizes.at(p); k++)
            {
                    coefs.row(tmpMap.index(k,p)) = m_patches.patch(p).coefs().row(k);
            }
            offset += sizes.at(p);
        }

        return coefs;
    }

//     template<short_t d,class T>
//     void gsDPatchBase<d,T>::setCoefficients(const gsMatrix<T> & coefs, gsMultiPatch<T> & mp) const
//     {
//         std::vector<index_t> sizes(mp.nPatches());
//         index_t totalsize = 0;
//         for (size_t p=0; p!=mp.nPatches(); p++) // patches
//         {
//             sizes.at(p) = mp.patch(p).coefs().rows();
//             totalsize += sizes.at(p);
//         }

//         GISMO_ASSERT(totalsize==coefs.rows(),"Sizes do not agree");

//         gsMultiBasis<T> basis(mp);
//         gsDofMapper tmpMap(basis);
//         tmpMap.finalize();

//         index_t offset = 0;
//         for (size_t p=0; p!=mp.nPatches(); p++) // patches
//         {
//             for (index_t k=0; k!=sizes.at(p); k++)
//             {
//                 mp.patch(p).coefs().row(k) = coefs.row(tmpMap.index(k,p));
//             }
//             offset += sizes.at(p);
//         }

//     }

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
        if (computeCoefs)
        {
            m_coefs = this->_preCoefficients(); // gets coefficients of the modified size
            m_coefs = m_matrix.transpose() * m_coefs; // maps to local size
        }

        ////////////////////////////////////////////////
        index_t size,offset = 0;
        for (index_t p=0; p!=patch; p++)
            offset += m_mapOriginal.patchSize(p);

        size = m_mapOriginal.patchSize(patch);
        gsMatrix<T> local = m_coefs.block(offset,0,size,m_patches.geoDim());
        return m_bases[patch].makeGeometry( give(local) ).release();
    }

    // template<short_t d,class T>
    // gsMultiPatch<T> gsDPatchBase<d,T>::exportToPatches()
    // {
    //     m_coefs = this->_preCoefficients();
    //     m_coefs = m_matrix.transpose() * m_coefs;

    //     std::vector<gsGeometry<T> *> patches(m_patches.nPatches());
    //     for (size_t p=0; p!=m_patches.nPatches(); p++)
    //         patches[p]= this->exportPatch(p,false);

    //     return gsMultiPatch<T>(patches,m_patches.boundaries(),m_patches.interfaces());
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
    const gsVector<index_t> gsDPatchBase<d,T>::_indicesFromVert(index_t index, const patchCorner corner, const patchSide side, index_t offset)
    {
        /*
            Finds indices i1,...,in in the direction of side away from the vertex
        */

        gsVector<index_t> result(index);

        gsBasis<T> * basis = &m_bases.basis(side.patch);

        gsVector<index_t> indices = static_cast<gsVector<index_t>>(basis->boundaryOffset(side.side(),offset));
        if (side.side()==1) //west
        {
            if (corner.corner()==1)//southwest
                result = indices.head(index);
            else if (corner.corner()==3) //northwest
                result = indices.tail(index);
            else
                GISMO_ERROR(corner.corner() << "is not on side "<<side.side()<<"!");
        }
        else if (side.side()==2) //east
        {
            if (corner.corner()==2)//southeast
                result = indices.head(index);
            else if (corner.corner()==4) //northeast
                result = indices.tail(index);
            else
                GISMO_ERROR(corner.corner() << "is not on side "<<side.side()<<"!");
        }
        else if (side.side()==3) //south
        {
            if (corner.corner()==1)//southwest
                result = indices.head(index);
            else if (corner.corner()==2) //southeast
                result = indices.tail(index);
            else
                GISMO_ERROR(corner.corner() << "is not adjacent to side "<<side.side()<<"!");
        }
        else if (side.side()==4) //north
        {
            if (corner.corner()==3)//northwest
                result = indices.head(index);
            else if (corner.corner()==4) //northeast
                result = indices.tail(index);
            else
                GISMO_ERROR(corner.corner() << "is not adjacent to side "<<side.side()<<"!");
        }
        return result;
    }

    template<short_t d,class T>
    const index_t gsDPatchBase<d,T>::_indexFromVert(index_t index, const patchCorner corner, const patchSide side, index_t offset, index_t levelOffset) const
    {
        return _indexFromVert(m_bases,index, corner, side, offset, levelOffset);
    }

    template<short_t d,class T>
    const index_t gsDPatchBase<d,T>::_indexFromVert(const gsMultiBasis<T> & bases, index_t index, const patchCorner corner, const patchSide side, index_t offset, index_t levelOffset) const
    {
        return _indexFromVert(&bases.basis(corner.patch),index, corner, side, offset, levelOffset);
    }

    template<short_t d,class T>
    const index_t gsDPatchBase<d,T>::_indexFromVert(const gsBasis<T> * basis, index_t index, const patchCorner corner, const patchSide side, index_t offset, index_t levelOffset) const
    {
        const gsTensorBSplineBasis<d,T> * tbbasis;
        const gsHTensorBasis<d,T> * thbasis;
        if ((tbbasis = dynamic_cast<const gsTensorBSplineBasis<d,T> * >(basis)))
            return _indexFromVert_impl(tbbasis,index,corner,side,offset,levelOffset);
        else if ((thbasis = dynamic_cast<const gsHTensorBasis<d,T> * >(basis)))
            return _indexFromVert_impl(thbasis,index,corner,side,offset,levelOffset);
        else
            GISMO_ERROR("Basis type unknown");
    }

    template<short_t d,class T>
    template<class U>
    typename util::enable_if<util::is_same<U,   const gsHTensorBasis<d,T> *>::value,
                                                const index_t>::type
    gsDPatchBase<d,T>::_indexFromVert_impl(U basis, index_t index, const patchCorner corner, const patchSide side, index_t offset, index_t levelOffset) const
    {
        if ((index==0) && (offset==0))
        {
            // gsHTensorBasis<d,T> *basis = dynamic_cast<gsHTensorBasis<d,T>*>(&bases.basis(corner.patch));
            index_t idx = basis->functionAtCorner(corner.corner());
            return idx;
        }
        else
        {
            std::vector<index_t> indices(1);
            indices.at(0) = index;
            std::vector<index_t> result = _indexFromVert(basis,indices,corner,side,offset,levelOffset);
            return result.at(0);
        }
    }

    template<short_t d,class T>
    template<class U>
    typename util::enable_if<util::is_same<U,   const gsTensorBSplineBasis<d,T> *>::value,
                                                const index_t>::type
    gsDPatchBase<d,T>::_indexFromVert_impl(U basis, index_t index, const patchCorner corner, const patchSide side, index_t offset, index_t levelOffset) const
    {
        if ((index==0) && (offset==0))
        {
            // gsHTensorBasis<d,T> *basis = dynamic_cast<gsHTensorBasis<d,T>*>(&bases.basis(corner.patch));
            index_t idx = basis->functionAtCorner(corner.corner());
            return idx;
        }
        else
        {
            std::vector<index_t> indices(1);
            indices.at(0) = index;
            std::vector<index_t> result = _indexFromVert(basis,indices,corner,side,offset,levelOffset);
            return result.at(0);
        }
    }

    template<short_t d,class T>
    const std::vector<index_t> gsDPatchBase<d,T>::_indexFromVert(const std::vector<index_t> & index, const patchCorner corner, const patchSide side, index_t offset, index_t levelOffset) const
    {
        return _indexFromVert(m_bases,index,corner,side,offset,levelOffset);
    }

    template<short_t d,class T>
    const std::vector<index_t> gsDPatchBase<d,T>::_indexFromVert(const gsMultiBasis<T> & bases, const std::vector<index_t> & index, const patchCorner corner, const patchSide side, index_t offset, index_t levelOffset) const
    {
        return _indexFromVert(bases,index, corner, side, offset, levelOffset);
    }


    template<short_t d,class T>
    const std::vector<index_t> gsDPatchBase<d,T>::_indexFromVert(const gsHTensorBasis<d,T> * basis, const std::vector<index_t> & index, const patchCorner corner, const patchSide side, index_t offset, index_t levelOffset) const
    {
        /*
            Finds indices i in the direction of side away from the vertex
            if index = 0, the corner index is requested
        */

        std::vector<index_t> result(index.size());
        for (size_t k=0; k!=index.size(); k++)
        {
            index_t level = basis->levelAtCorner(corner.corner()) + levelOffset;
            gsTensorBSplineBasis<d,T> tbasis = basis->tensorLevel(level);
            result = this->_indexFromVert(&tbasis,index,corner,side,offset);
            if (levelOffset==0) // if levelOffset!=0, the index in the original basis is requested
                result[k] = basis->flatTensorIndexToHierachicalIndex(result[k],level);
        }

        return result;
    }

    template<short_t d,class T>
    const std::vector<index_t> gsDPatchBase<d,T>::_indexFromVert(const gsTensorBSplineBasis<d,T> * tbasis, const std::vector<index_t> & index, const patchCorner corner, const patchSide side, index_t offset, index_t levelOffset) const
    {
        /*
            Finds indices i in the direction of side away from the vertex
            if index = 0, the corner index is requested
        */

        std::vector<index_t> result(index.size());

        gsVector<index_t,2> sizes;
        tbasis->size_cwise(sizes);
        for (size_t k=0; k!=index.size(); k++)
        {
            if (side.side()==1) //west
            {
                if (corner.corner()==1)//southwest
                    result[k] = tbasis->index(offset,index[k]);
                else if (corner.corner()==3) //northwest
                    result[k] = tbasis->index(offset,sizes[1]-1-index[k],0);
                else
                    GISMO_ERROR(corner.corner() << " is not adjacent to side "<<side.side()<<"!");
            }
            else if (side.side()==2) //east
            {
                if (corner.corner()==2)//southeast
                    result[k] = tbasis->index(sizes[0]-1-offset,index[k]);
                else if (corner.corner()==4) //northeast
                    result[k] = tbasis->index(sizes[0]-1-offset,sizes[1]-1-index[k]);
                else
                    GISMO_ERROR(corner.corner() << " is not adjacent to side "<<side.side()<<"!");
            }
            else if (side.side()==3) //south
            {
                if (corner.corner()==1)//southwest
                    result[k] = tbasis->index(index[k],offset);
                else if (corner.corner()==2) //southeast
                    result[k] = tbasis->index(sizes[0]-1-index[k],offset);
                else
                    GISMO_ERROR(corner.corner() << " is not adjacent to side "<<side.side()<<"!");
            }
            else if (side.side()==4) //north
            {
                if (corner.corner()==3)//northwest
                    result[k] = tbasis->index(index[k],sizes[1]-1-offset,0);
                else if (corner.corner()==4) //northeast
                    result[k] = tbasis->index(sizes[0]-1-index[k],sizes[1]-1-offset);
                else
                    GISMO_ERROR(corner.corner() << " is not adjacent to side "<<side.side()<<"!");

            }
        }

        return result;
    }

    template<short_t d,class T>
    const std::pair<index_t,bool> gsDPatchBase<d,T>::_vertexData(const patchCorner corner) const
    {
        std::vector<patchCorner> corners;
        std::pair<index_t,bool> output;
        output.second = m_patches.getCornerList(corner,corners); // bool is true if interior vertex
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
                if ( m_patches.isBoundary(csides[i]) ) continue;
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

        m_patches.getCornerList(pcorner,pcorners); // bool is true if interior vertex
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

        return (colSums.array() < 1+1e-8).any() && (colSums.array() > 1-1e-8).any();
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
//             for(int c=1;c<=4;++c)
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
//             std::vector<Eigen::Triplet<T,index_t>> tripletList;
//             for (size_t p=0; p!=m_RefPatches.nPatches(); p++)
//             {
//                 gsHTensorBasis<d,T> *basis = dynamic_cast<gsHTensorBasis<d,T>*>(&m_RefPatches.basis(p));
//                 std::vector< gsSortedVector< index_t > > xmat = basis->getXmatrix();

//                 m_RefPatches.patch(p).refineElements(elVec[p]);

//                 basis->transfer(xmat,tmp);

//                 for (index_t i = 0; i<tmp.outerSize(); ++i)
//                     for (typename gsSparseMatrix<T>::iterator it(tmp,i); it; ++it)
//                         tripletList.push_back(Eigen::Triplet<T,index_t>(it.row()+rows,it.col()+cols,it.value()));

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

//         // gsWriteParaview<>(m_RefPatches,"mp_ref",1000,true);
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
//         // m_patches.getEVs(cornerLists);

//         // get the corners that need refinement
//         std::vector<patchCorner> cornerList;
//         patchCorner pcorner;
//         index_t cidx;
//         for(size_t p = 0;p<m_RefPatches.nPatches();++p)
//         {
//             for(int c=1;c<=4;++c)
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
        _initCoefs();
    }

    template<short_t d,class T>
    void gsDPatchBase<d,T>::_initChecks()
    {
        m_patches.checkConsistency();
        m_nSides = 2*m_patches.nInterfaces() + m_patches.nBoundary();
        m_nVerts = 4*m_patches.nPatches();

        m_sideCheck.resize(m_nSides);
        std::fill(m_sideCheck.begin(), m_sideCheck.end(), false);
        m_vertCheck.resize(m_nVerts);
        std::fill(m_vertCheck.begin(), m_vertCheck.end(), false);
    }

    template<short_t d,class T>
    void gsDPatchBase<d,T>::_initTHB()
    {
        // Cast all patches of the mp object to THB splines
        gsTHBSpline<d,T> thb;
        gsTensorBSpline<d,T> * geo;
        for (size_t k=0; k!=m_patches.nPatches(); ++k)
        {
            if ( (geo = dynamic_cast< gsTensorBSpline<d,T> * > (&m_patches.patch(k))) )
            {
                thb = gsTHBSpline<d,T>(*geo);
                m_patches.patch(k) = thb;
            }
            else if (dynamic_cast< gsTHBSpline<d,T> * > (&m_patches.patch(k)))
            { }
            else
                gsWarn<<"No THB basis was constructed";
        }
    }

    template<short_t d,class T>
    void gsDPatchBase<d,T>::_initBasis()
    {
        m_Bbases = m_bases = gsMultiBasis<T>(m_patches);
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
    void gsDPatchBase<d,T>::_initCoefs()
    {
        m_coefs = this->allCoefficients();
    }

    template<short_t d,class T>
    void gsDPatchBase<d,T>::_whichHandled()
    {
        for (size_t p=0; p!=m_patches.nPatches(); p++)
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
        GISMO_ASSERT(checkSides,"Not all sides are checked");
        bool checkVerts = std::all_of(m_vertCheck.begin(), m_vertCheck.end(), [](bool m_vertCheck) { return m_vertCheck; });
        GISMO_ASSERT(checkVerts,"Not all vertices are checked");

        if (!basis)
            return;

        bool checkBasis = std::all_of(m_basisCheck.begin(), m_basisCheck.end(), [](bool m_basisCheck) { return m_basisCheck; });
        GISMO_ASSERT(checkBasis,"Not all basis functions are checked");
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
    void gsDPatchBase<d,T>::_handleInterior()
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
    void gsDPatchBase<d,T>::_computeSmoothMatrix()
    {
        GISMO_ASSERT((size_t)m_mapModified.freeSize()==m_size,"Size does not match predicted size, m_mapModified.freeSize()="<<m_mapModified.freeSize()<<"; m_size="<<m_size);

        _resetChecks(true);

        // iterate over the vertices
#pragma omp parallel
{
        #pragma omp parallel for collapse(2)
        for (size_t p=0; p<m_patches.nPatches(); p++)
            for (index_t c=1; c<5; c++)
                _handleVertex(patchCorner(p,c));

        #pragma omp parallel for
        for(gsBoxTopology::const_iiterator iit = m_patches.iBegin(); iit< m_patches.iEnd(); iit++)
            _handleInterface(*iit);

        // boundaries
        #pragma omp parallel for
        for(gsBoxTopology::const_biterator bit = m_patches.bBegin(); bit< m_patches.bEnd(); bit++)
            _handleBoundary(*bit);

        _handleInterior();
}
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
        for(gsBoxTopology::const_iiterator iit = m_patches.iBegin(); iit!= m_patches.iEnd(); iit++)
            _computeInterfaceMapper(*iit);
// }
        // On the boundaries, we don't do anything
// #pragma omp parallel
// {
        // #pragma omp parallel for
        for(gsBoxTopology::const_biterator bit = m_patches.bBegin(); bit!= m_patches.bEnd(); bit++)
            _computeBoundaryMapper(*bit);
// }

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
            if (m_patches.isInterface(psides[p]))
            {
                for (index_t k=0; k!=2; k++)
                    m_mapModified.eliminateDof(_indexFromVert(k,pcorner,psides[p],0,0),pcorner.patch);
            }
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
            if (m_patches.isInterface(psides[p]))
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
        // m_patches.getCornerList(pcorner,pcorners);
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
            if (m_patches.isInterface(psides[p]))
                m_mapModified.eliminateDof(_indexFromVert(1,pcorner,psides[p],0,0),pcorner.patch);
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
        index_t iindex = m_patches.isInterface(psides[0]) ? 0 : 1;

        GISMO_ENSURE(m_patches.getInterface(psides[iindex],iface),"Must be an interface");

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
            m_vertCheck[ _vertIndex(pcorner.patch, pcorner.corner()) ] = true;
        }
    }

} // namespace gismo
