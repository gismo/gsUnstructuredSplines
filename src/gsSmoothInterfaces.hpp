/** @file gsSmoothInterfaces.hpp

    @brief Creates the D-Patch smoothing matrix.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

// #define PI 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647

#include<gsHSplines/gsTHBSpline.h>

namespace gismo
{
    // Constructors
    template<short_t d,class T>
    gsSmoothInterfaces<d,T>::gsSmoothInterfaces(const gsMultiPatch<T> & patches)
    :
    Base(patches)
    {
        this->defaultOptions();
    }

    /*=====================================================================================
                                    Coefficients
    =====================================================================================*/
    template<short_t d,class T>
    gsMatrix<T> gsSmoothInterfaces<d,T>::_preCoefficients(const gsMultiPatch<T> & patches)
    {
        GISMO_ASSERT(m_mapModified.isFinalized(),"Mapper is not finalized, run XXXX first");
        GISMO_ASSERT(!patches.empty(),"The reference multipatch is empty!");

        gsMatrix<T> coefs(m_mapModified.freeSize(),patches.geoDim());

        index_t size;
        for (size_t p=0; p!=m_bases.nBases(); p++) // patches
        {
            size = m_mapModified.patchSize(p);
            for (index_t k=0; k!=size; k++)
            {
                if (m_mapModified.is_free(k,p))
                    coefs.row(m_mapModified.index(k,p,0)) = patches.patch(p).coefs().row(k);
            }
        }

        return coefs;
    }

    /*=====================================================================================
                                    Construction functions
    =====================================================================================*/

    // template<short_t d,class T>
    // gsGeometry<T>* gsSmoothInterfaces<d,T>::exportPatch(index_t patch, bool computeCoefs)
    // {
    //     GISMO_ASSERT(m_computed,"The method has not been computed! Call compute().");
    //     ////////////////////////////////////////////////
    //     // This can be done more efficient!!
    //     // Do it once instead of for every patch
    //     ////////////////////////////////////////////////
    //     if (computeCoefs)
    //     {
    //         m_coefs = this->_preCoefficients(); // gets coefficients of the modified size
    //         m_coefs = m_matrix.transpose() * m_coefs; // maps to local size
    //     }

    //     ////////////////////////////////////////////////
    //     index_t size,offset = 0;
    //     for (index_t p=0; p!=patch; p++)
    //         offset += m_mapOriginal.patchSize(p);

    //     size = m_mapOriginal.patchSize(patch);
    //     gsMatrix<T> local = m_coefs.block(offset,0,size,m_patches.geoDim());


    //     bool check = dynamic_cast< gsTensorBSplineBasis<d,T> * > (&m_Bbases.basis(patch));
    //     gsDebugVar(check);

    //     return m_Bbases[patch].makeGeometry( give(local) ).release();
    // }

    template<short_t d,class T>
    void gsSmoothInterfaces<d,T>::_initBasis()
    {
        m_bases = m_Bbases;
    }

    template<short_t d,class T>
    void gsSmoothInterfaces<d,T>::_makeTHB()
    {
    }

    template<short_t d,class T>
    void gsSmoothInterfaces<d,T>::_initTHB()
    {

    }

    template<short_t d,class T>
    void gsSmoothInterfaces<d,T>::_computeEVs()
    {
        m_matrix.makeCompressed();
    }

    template<short_t d,class T>
    void gsSmoothInterfaces<d,T>::_countDoFs() // also initialize the mappers!
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
                        tmp += 2*vdata.first+2;
                    else if ((!vdata.second) && vdata.first==2 && C0)
                        tmp += 6;
                    else if ((!vdata.second) && vdata.first>2 && C0)
                        tmp += 2*vdata.first+2;
                    else
                        tmp += vdata.first; // valence;

                    // corn +=1;
                }
            }
        // gsDebug<<"Number of unique corners: "<<corn<<"\n";

        // gsDebug<<"Number of vertex DoFs: "<<tmp<<"\n";

        m_size += tmp;
    }

} // namespace gismo