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
    DPatch(patches)
    {
        this->defaultOptions();
    }

    /*=====================================================================================
                                    Coefficients
    =====================================================================================*/
    template<short_t d,class T>
    gsMatrix<T> gsSmoothInterfaces<d,T>::_preCoefficients()
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
    void gsSmoothInterfaces<d,T>::_makeTHB()
    {
        m_RefPatches = gsMultiPatch<T>(m_patches);
        // m_bases = m_Bbases;

        // gsTHBSpline<d,T> thb;
        // gsTensorBSpline<d,T> * geo;
        // for (size_t k=0; k!=m_RefPatches.nPatches(); ++k)
        // {
        //     if ( (geo = dynamic_cast< gsTensorBSpline<d,T> * > (&m_RefPatches.patch(k))) )
        //     {
        //         gsDebugVar("BSpline");
        //     }
        //     else if (dynamic_cast< gsTHBSpline<d,T> * > (&m_RefPatches.patch(k)))
        //     {
        //         gsDebugVar("THBSpline");
        //     }
        //     else
        //         gsWarn<<"No THB basis was constructed";
        // }

        // gsTHBSplineBasis<d,T> * thb_basis;
        // gsTensorBSplineBasis<d,T> * geo_basis;
        // for (size_t k=0; k!=m_bases.nBases(); ++k)
        // {
        //     if ( (geo_basis = dynamic_cast< gsTensorBSplineBasis<d,T> * > (&m_bases.basis(k))) )
        //     {
        //         gsDebugVar("BSpline");
        //     }
        //     else if ( (thb_basis = dynamic_cast< gsTHBSplineBasis<d,T> * > (&m_bases.basis(k))) )
        //     {
        //         m_bases.addBasis(&thb_basis->tensorLevel(0));
        //         bool check = dynamic_cast< gsTHBSplineBasis<d,T> * > (&m_bases.basis(k));
        //         gsDebugVar(check);
        //         gsDebugVar("THBSpline");
        //     }
        //     else
        //         gsWarn<<"No THB basis was constructed";
        // }

    }

    template<short_t d,class T>
    void gsSmoothInterfaces<d,T>::_initTHB()
    {
    //     // Cast all patches of the mp object to B splines
    //     gsTHBSpline<d,T> * thb;
    //     gsTensorBSpline<d,T> geo;
    //     gsTensorBSplineBasis<d,T> basis;
    //     gsDebugVar("hi");
    //     for (size_t k=0; k!=m_RefPatches.nPatches(); ++k)
    //     {
    //         bool test1 = dynamic_cast< gsTHBSpline<d,T> * > (&m_RefPatches.patch(k));
    //         bool test2 = dynamic_cast< gsTensorBSpline<d,T> * > (&m_RefPatches.patch(k));
    //         gsDebugVar(test1);
    //         gsDebugVar(test2);
    //         if ( (thb = dynamic_cast< gsTHBSpline<d,T> * > (&m_RefPatches.patch(k))) )
    //         {
    //             gsDebugVar(thb->basis().maxLevel());
    //             basis = thb->basis().tensorLevel(0);
    //             geo = gsTensorBSpline<d,T>(basis,thb->coefs());
    //             m_RefPatches.patch(k) = geo;
    //         }
    //         else if (dynamic_cast< gsTensorBSpline<d,T> * > (&m_RefPatches.patch(k)))
    //         { }
    //         else
    //             gsWarn<<"No THB basis was constructed";
    //     }
    }

    template<short_t d,class T>
    void gsSmoothInterfaces<d,T>::_computeEVs()
    {
        m_matrix.makeCompressed();
        // gsDebugVar(m_matrix.toDense());
    }

} // namespace gismo