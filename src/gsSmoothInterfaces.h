/** @file gsSmoothInterfaces.h

    @brief Creates the D-Patch smoothing matrix.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#pragma once

// #include <gsUnstructuredSplines/src/gsDPatchBase.h>
#include <gsUnstructuredSplines/src/gsAlmostC1.h>
// #include <gsUnstructuredSplines/src/gsDPatchBase.hpp>

namespace gismo
{


/**
 * @brief      Constructs the D-Patch, from which the transformation matrix can be called
 *
 * @tparam     d     parametric dimension
 */
template<short_t d,class T>
class gsSmoothInterfaces :  public gsAlmostC1<d,T>
{

public:
    // using Base = gsDPatchBase<d,T>;
    using DPatch = gsAlmostC1<d,T>;

    /// Shared pointer for gsSmoothInterfaces
    typedef memory::shared_ptr< gsSmoothInterfaces > Ptr;

    /// Unique pointer for gsSmoothInterfaces
    typedef memory::unique_ptr< gsSmoothInterfaces > uPtr;

    /// Empty constructor
    gsSmoothInterfaces() : DPatch()
    { }

    ~gsSmoothInterfaces() {}

//    using DPatch::compute;

    /**
     * @brief      Default constructor
     *
     * @param      mp    Multipatch of the geometry
     */
    gsSmoothInterfaces(gsMultiPatch<T> const & mp) ;

    GISMO_CLONE_FUNCTION(gsSmoothInterfaces)

    // ~gsSmoothInterfaces();


//    using DPatch::exportToPatches;


protected:
    /**
     * @brief       Computes the C1 coefficients for pre-multiplication to make the multipatch
     *
     * Takes the coefficients which are tagged as "free" in the modified DoFMapper (m_mapModified) and when a boundary vertex with valence=3 is present, this one is shifted.
     *
     */
    gsMatrix<T> _preCoefficients();
//    using DPatch::_preCoefficients;

//    using DPatch::allCoefficients;

//    using DPatch::exportPatch;
    // gsGeometry<T> * exportPatch(index_t patch, bool computeCoefs);

protected:

//    using DPatch::_indexFromSides;

//    using DPatch::_indicesFromVert;

//    using DPatch::_indexFromVert;

//    using DPatch::_vertexData;

//    using DPatch::_sideIndex;

//    using DPatch::_vertIndex;

//    using DPatch::_getLowestCorners;

//    using DPatch::_removeLowestCorners;

//    using DPatch::_getLowestIndices;

//    using DPatch::_removeLowestIndices;

//    using DPatch::_getInterfaceIndices;

//    using DPatch::_getAllInterfaceIndices;

protected:

//    using DPatch::_countDoFs;

//    using DPatch::_computeMapper; // also initialize the mappers!

//    using DPatch::_computeSmoothMatrix;

    void _makeTHB();
    void _initTHB();

//    using DPatch::_makeTHB;

    void _computeEVs();
//    using DPatch::_computeEVs;

    /**
     * @brief      Makes the Pi matrix
     *
     * This matrix is used to transform the coefficients of the D-Patch smoothing matrix
     *
     * @param[in]  valence  The valence
     *
     * @return     Matrix for smoothing around an EV}
     */
//    using DPatch::_makePi;

protected:

//    using DPatch::_performChecks;
//    using DPatch::_resetChecks;

protected:

/*    *
     * @brief      Handles a vertex in the global matrix
     *
     * We use the following notation convention (per patch!):
     * b00 is the basis function at the vertex
     * b10 is the basis function next to the vertex along the first interface that connects to the vertex
     * b20 is the basis function next to b10 along the first interface that connects to the vertex
     * etc.
     *
     * b01 is the basis function next to the vertex along the second interface that connects to the vertex
     * b02 is the basis function next to b01 along the second interface that connects to the vertex
     * etc.
     *
     * b11 is the basis function with offset 1 from both interfaces and from the vertex itself
     * b22 is the basis function with offset 2 from both interfaces and from the vertex itself
     * etc.
     *
     * There are different options.
     * a) Boundary vertices
     *      i)  Valence 1: b00, b10, b01 and b00 all get weight 1.0 w.r.t the same basis function in the local basis
     *      ii) Valence 2: This case contains an interface between two patches. We use index k to denote the row basis functions along the interface. So k=0 corresponds to the basis functions on the boundary and k=1 corresponds to the basis functions with offset 1 from the boundaries. Using this convention, the functions bk1 in the local basis, are coupled to bk1 in the global basis with weight 1. The functions bk0 in the local basis (on the considered patch) are coupled to bk1 in the global basis with weight 0.5. The functions bk0 in the local basis (on the other patch) are coupled to bk1 (on the considered patch) in the global basis with weight 0.5.
     *      iii)Valence 3: In this case, the matched vertices on all the adjacent patches are treated in one go! Note that all the basis functions corresponding to the vertex (b00) on all patches are matched! We couple the b00 functions of all patches (in the local basis) with weight 1/4 to the b00 of the adjacent patch with the lowest number in the global basis. Then, the b11 on the considered patch is coupled with weight 1 to itself and with weight 0.25 to the b00s of the other patches. Then, we will handle the vertices where an interface and a boundary meet (there are two of these). For the patch corners that are on an interface, we find the b11 and b10 vertices (orthogonal to the interface) and we give all b10s weight 0.5 w.r.t. the b11s in the global basis (on both patches). Lastly, we add weight 0.5 for the b10s along the boundaries (so only for two patches) to the (matched) b00 basis function (all b00s refer to the same dof in the global basis).
     * b) Interior vertices (all valences):
     *      i)  b11 gets weight 1.0 w.r.t the same basis function in the local basis
     *      ii) all associated b00s in the local basis get weight 1/valence w.r.t. b11 in the global basis
     *      iii)the b10 and b01 in the local basis get weight 1/2 w.r.t. b11 in the global basis
     *
     * @param[in]  pcorner  The patchcorner
*/

//    using DPatch::_handleVertex;
    /**
     * @brief      Handles an interface in the global matrix
     *
     * Gives all the DoFs that have offset 1 (orthogonal) from the interface weight 1.0 w.r.t itself. All the DoFs ON the interface (on both patches) will have weight 0.5 to the DoF with offset 1.
     * This interface handling excludes the indices that are in the 0 and 1 ring around vertices.
     *
     * @param[in]  iface  The interface
     */
//    using DPatch::_handleInterface;
    /**
     * @brief      Handles a boundary in the global matrix
     *
     * Handles all DoFs on the boundary with unit-weight, except the ones in the 0 and 1 rings around the vertices.
     *
     * @param[in]  side  The boundary side
     */
//    using DPatch::_handleBoundary;
    /**
     * @brief      Handles the interior in the global matrix
     *
     * Gives all left-over DoFs, which are in the interior, weight 1 w.r.t. itself
     */
//    using DPatch::_handleInterior;
    /**
     * @brief      Prints which DoFs have been handled and which have been eliminated
     */

protected:
//    using DPatch::_whichHandled;

protected:
   using DPatch::m_patches;
   using DPatch::m_computed;
   using DPatch::m_RefPatches;
   using DPatch::m_bases;
   using DPatch::m_Bbases;
   using DPatch::m_tMatrix;
   using DPatch::m_sideCheck;
   using DPatch::m_vertCheck;
   using DPatch::m_basisCheck;
   using DPatch::m_C0s;

   using DPatch::m_mapModified;
   using DPatch::m_mapOriginal;

   using DPatch::m_matrix;

   using DPatch::m_options;

   using DPatch::m_size;

   using DPatch::m_coefs;

   using DPatch::m_nSides;
   using DPatch::m_nVerts;

};



}

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsSmoothInterfaces.hpp)
#endif