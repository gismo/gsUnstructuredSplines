/** @file gsSmoothInterfaces.h

    @brief Creates the D-Patch smoothing matrix.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#pragma once

#include <gsUnstructuredSplines/src/gsDPatchBase.h>

namespace gismo
{


/**
 * @brief      Constructs the D-Patch, from which the transformation matrix can be called
 *
 * @tparam     d     parametric dimension
 */
template<short_t d,class T>
class gsSmoothInterfaces :  public gsDPatchBase<d,T>
{

public:
    using Base = gsDPatchBase<d,T>;

    /// Shared pointer for gsSmoothInterfaces
    typedef memory::shared_ptr< gsSmoothInterfaces > Ptr;

    /// Unique pointer for gsSmoothInterfaces
    typedef memory::unique_ptr< gsSmoothInterfaces > uPtr;

    /// Empty constructor
    gsSmoothInterfaces() : Base()
    { }

    ~gsSmoothInterfaces() {}

//    using Base::compute;

    /**
     * @brief      Default constructor
     *
     * @param      mp    Multipatch of the geometry
     */
    gsSmoothInterfaces(gsMultiPatch<T> const & mp) ;

    GISMO_CLONE_FUNCTION(gsSmoothInterfaces)

    // ~gsSmoothInterfaces();

    using Base::defaultOptions;

protected:
    /**
     * @brief       Computes the C1 coefficients for pre-multiplication to make the multipatch
     *
     * Takes the coefficients which are tagged as "free" in the modified DoFMapper (m_mapModified) and when a boundary vertex with valence=3 is present, this one is shifted.
     *
     */
    gsMatrix<T> _preCoefficients(const gsMultiPatch<T> & patches);
//    using Base::_preCoefficients;

//    using Base::allCoefficients;

//    using Base::exportPatch;
    // gsGeometry<T> * exportPatch(index_t patch, bool computeCoefs);

protected:

//    using Base::_indexFromSides;

//    using Base::_indicesFromVert;

//    using Base::_indexFromVert;

//    using Base::_vertexData;

//    using Base::_sideIndex;

//    using Base::_vertIndex;

//    using Base::_getLowestCorners;

//    using Base::_removeLowestCorners;

//    using Base::_getLowestIndices;

//    using Base::_removeLowestIndices;

//    using Base::_getInterfaceIndices;

//    using Base::_getAllInterfaceIndices;

protected:

    void _countDoFs() override;

    void _initBasis() override;

    void _makeTHB() override;

    void _initTHB() override;

    void _computeEVs() override;

    /**
     * @brief      Makes the Pi matrix
     *
     * This matrix is used to transform the coefficients of the D-Patch smoothing matrix
     *
     * @param[in]  valence  The valence
     *
     * @return     Matrix for smoothing around an EV}
     */
//    using Base::_makePi;

protected:

//    using Base::_performChecks;
//    using Base::_resetChecks;

protected:

protected:
//    using Base::_whichHandled;

protected:
   using Base::m_patches;
   using Base::m_topology;
   using Base::m_computed;
   using Base::m_bases;
   using Base::m_Bbases;
   using Base::m_tMatrix;
   using Base::m_sideCheck;
   using Base::m_vertCheck;
   using Base::m_basisCheck;
   using Base::m_C0s;

   using Base::m_mapModified;
   using Base::m_mapOriginal;

   using Base::m_matrix;

   using Base::m_options;

   using Base::m_size;

   using Base::m_coefs;

   using Base::m_nSides;
   using Base::m_nVerts;

};



}

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsSmoothInterfaces.hpp)
#endif
