/** @file gsDPatch.h

    @brief Creates the D-Patch smoothing matrix.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#pragma once

#include <gsCore/gsBoxTopology.h>
#include <gsCore/gsMultiPatch.h>
#include <gsIO/gsOptionList.h>
#include <gsUnstructuredSplines/src/gsDPatchBase.h>
// #include <gsUnstructuredSplines/src/gsDPatchBase.hpp>
namespace gismo
{


/**
 * @brief      Constructs the D-Patch, from which the transformation matrix can be called
 *
 * @tparam     d     parametric dimension
 *
 * @ingroup    UnstructuredSplines
 */
template<short_t d,class T>
class gsDPatch : public gsDPatchBase<d,T>
{

typedef typename std::vector<std::tuple<index_t,index_t,T>> sparseEntry_t;

public:
    using Base = gsDPatchBase<d,T>;

    /// Shared pointer for gsDPatch
    typedef memory::shared_ptr< gsDPatch > Ptr;

    /// Unique pointer for gsDPatch
    typedef memory::unique_ptr< gsDPatch > uPtr;

    /// Empty constructor
    gsDPatch() : Base()
    { }

    // using Base::compute;

    /**
     * @brief      Default constructor
     *
     * @param      mp    Multipatch of the geometry
     */
    gsDPatch(const gsMultiBasis<T> & mb) ;


    /**
     * @brief      Default constructor
     *
     * @param      mp    Multipatch of the geometry
     */
    gsDPatch(const gsMultiPatch<T> & mp) ;

    GISMO_CLONE_FUNCTION(gsDPatch)

    ~gsDPatch();

    void defaultOptions();

    // using Base::exportToPatches;

protected:
    /**
     * @brief       Computes the C1 coefficients for pre-multiplication to make the multipatch
     *
     * Takes the coefficients which are tagged as "free" in the modified DoFMapper (m_mapModified) and when a boundary vertex with valence=3 is present, this one is shifted.
     *
     */
    gsMatrix<T> _preCoefficients(const gsMultiPatch<T> & patches) override;
    using Base::_preCoefficients;

    // using Base::exportPatch;

protected:

    // using Base::_indexFromSides;

    // using Base::_indicesFromVert;

    // using Base::_indexFromVert;

    // using Base::_vertexData;

    // using Base::_sideIndex;

    // using Base::_vertIndex;

    // using Base::_getLowestCorners;

    // using Base::_removeLowestCorners;

    // using Base::_getLowestIndices;

    // using Base::_removeLowestIndices;

    // using Base::_getInterfaceIndices;

    // using Base::_getAllInterfaceIndices;

protected:

    using Base::_performChecks;
    using Base::_resetChecks;

    void _countDoFs();

    // void _computeMapper(); // also initialize the mappers!
    using Base::_computeMapper;

    // void _computeSmoothMatrix();
    using Base::_computeSmoothMatrix;

    void _initTHB() override;

    void _refBoxes(std::vector<std::vector<index_t>> & patchBoxes);

    void _initBasis() override;

    void _makeTHB() override;

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
    gsMatrix<T> _makePi(index_t valence);

    using Base::getSharpCorners;
    using Base::_indexFromSides;
    using Base::_indexFromVert;
    using Base::_vertexData;
    using Base::_sideIndex;
    using Base::_vertIndex;
    using Base::_getLowestCorners;
    using Base::_removeLowestCorners;
    using Base::_getLowestIndices;
    using Base::_removeLowestIndices;

    using Base::_push;
    using Base::_pushAndCheck;



protected:
    // void _computeInterfaceMapper(boundaryInterface iface);
    using Base::_computeInterfaceMapper;

    // void _computeBoundaryMapper(patchSide boundary);
    using Base::_computeBoundaryMapper;

    void _computeVertexMapper(patchCorner pcorner);


private:
    // // Boundary vertex of valence 1
    // template<bool _boundary, index_t _v> // valence=2
    // typename std::enable_if<  _boundary && _v==1, void>::type
    // // SAME
    // _computeVertexMapperBoundary_v1(patchCorner pcorner, index_t valence);

    using Base::_computeMapperRegularCorner_v1;

    // // Boundary vertex of valence 2 with C1 smoothness
    // template<bool _boundary, index_t _v, bool _smooth> // valence=2
    // typename std::enable_if<  _boundary && _v==2 && _smooth, void>::type
    // // SAME
    // _computeVertexMapperBoundarySmooth_v2(patchCorner pcorner, index_t valence);
    using Base::_computeMapperRegularBoundaryVertexSmooth_v2;

    // // Boundary vertex of valence 2 with C0 smoothness
    // template<bool _boundary, index_t _v, bool _smooth> // valence=2
    // typename std::enable_if<  _boundary && _v==2 && (!_smooth), void>::type
    // // DIFFERENT
    // _computeVertexMapperBoundaryNonSmooth_v2(patchCorner pcorner, index_t valence);
    using Base::_computeMapperRegularBoundaryVertexNonSmooth_v2;

    // Boundary vertex of valence 3 with C1 smoothness
    // ONLY DPATCH
    void _computeMapperIrregularBoundaryVertexSmooth_v3(patchCorner pcorner, index_t valence);

    // Boundary vertex of valence 3 with C0 smoothness
    // template<bool _boundary, index_t _v, bool _smooth> // valence=2
    // typename std::enable_if<  _boundary && _v==3 && (!_smooth), void>::type
    // ONLY DPATCH
    void _computeMapperIrregularBoundaryVertexNonSmooth_v3(patchCorner pcorner, index_t valence);

    // Boundary vertex of valence !(1,2,3) with C1 smoothness
    // template<bool _boundary, index_t _v, bool _smooth>
    // typename std::enable_if<  _boundary && _v==-1 && _smooth, void>::type
    // DIFFERENT
    void _computeMapperIrregularBoundaryVertexSmooth_v(patchCorner pcorner, index_t valence);

    // Boundary vertex of valence !(1,2,3) with C0 smoothness
    // template<bool _boundary, index_t _v, bool _smooth>
    // typename std::enable_if<  _boundary && _v==-1 && (!_smooth), void>::type
    // DIFFERENT
    void _computeMapperIrregularBoundaryVertexNonSmooth_v(patchCorner pcorner, index_t valence);

    // // Interior vertex
    // template<bool _boundary, index_t _v>
    // typename std::enable_if<  (!_boundary) && _v==-1, void>::type
    // // DIFFERENT
    // _computeVertexMapperInterior_v(patchCorner pcorner, index_t valence);
    using Base::_computeMapperInteriorVertex_v;

protected:

    /**
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
    // void _handleVertex(patchCorner pcorner);
    // interior vertices
    // void _handleInteriorVertex(patchCorner pcorner, index_t valence);

    /**
     * @brief      Handles an interface in the global matrix
     *
     * Gives all the DoFs that have offset 1 (orthogonal) from the interface weight 1.0 w.r.t itself. All the DoFs ON the interface (on both patches) will have weight 0.5 to the DoF with offset 1.
     * This interface handling excludes the indices that are in the 0 and 1 ring around vertices.
     *
     * @param[in]  iface  The interface
     */
    // void _handleInterface(boundaryInterface iface);
    /**
     * @brief      Handles a boundary in the global matrix
     *
     * Handles all DoFs on the boundary with unit-weight, except the ones in the 0 and 1 rings around the vertices.
     *
     * @param[in]  side  The boundary side
     */
    // void _handleBoundary(patchSide side);
    /**
     * @brief      Handles the interior in the global matrix
     *
     * Gives all left-over DoFs, which are in the interior, weight 1 w.r.t. itself
     */
    // void _handleInterior();
    /**
     * @brief      Prints which DoFs have been handled and which have been eliminated
     */

private:
    /**
     * @brief      Handles a regular corner
     *
     * @param[in]  pcorner  The pcorner
     */
    // void _handleRegularCorner(patchCorner pcorner);


    // template<bool _regular, bool _smooth> // valence=2
    // typename std::enable_if<  _regular  &&   _smooth   , void>::type
    // _handleBoundaryVertex(patchCorner pcorner, index_t valence);

    // template<bool _regular, bool _smooth> // valence=2
    // typename std::enable_if<  _regular  && (!_smooth)   , void>::type
    // _handleBoundaryVertex(patchCorner pcorner, index_t valence);

    // template<bool _regular, bool _smooth> // valence > 2
    // typename std::enable_if<(!_regular) &&   _smooth    , void>::type
    // _handleBoundaryVertex(patchCorner pcorner, index_t valence);

    // template<bool _regular, bool _smooth> // valence > 1
    // typename std::enable_if<(!_regular) && (!_smooth)   , void>::type
    // _handleBoundaryVertex(patchCorner pcorner, index_t valence);

    // void _handleRegularBoundaryVertexSmooth(patchCorner pcorner, index_t valence);

    void _handleIrregularBoundaryVertexSmooth(patchCorner pcorner, index_t valence);

    void _handleIrregularBoundaryVertexNonSmooth(patchCorner pcorner, index_t valence);

protected:
    using Base::_whichHandled;

protected:
    using Base::m_computed;
    using Base::m_topology;

    using Base::m_bases;
    gsMultiBasis<T> m_bases0;
    using Base::m_Bbases;
    using Base::m_tMatrix;
    std::vector<gsSparseMatrix<T>> m_tMatrices;
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
#include GISMO_HPP_HEADER(gsDPatch.hpp)
#endif
