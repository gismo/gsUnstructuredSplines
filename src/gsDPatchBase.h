/** @file gsDPatchBase.h

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

// #include <gsUnstructuredSplines/src/gsDPatchBase.hpp>

namespace gismo
{


/**
 * @brief      Constructs the D-Patch, from which the transformation matrix can be called
 *
 * @tparam     d     parametric dimension
 */
template<short_t d,class T>
class gsDPatchBase
{
protected:
    typedef typename std::vector<std::tuple<index_t,index_t,T>> tuple_t;

public:

    /// Shared pointer for gsDPatchBase
    typedef memory::shared_ptr< gsDPatchBase > Ptr;

    /// Unique pointer for gsDPatchBase
    typedef memory::unique_ptr< gsDPatchBase > uPtr;

    /// Empty constructor
    gsDPatchBase() : m_patches(gsMultiPatch<T>())
    { }


    /**
     * @brief      Default constructor
     *
     * @param      mp    Multipatch of the geometry
     */
    gsDPatchBase(gsMultiPatch<T> const & mp)
    :
    m_patches(mp)
    {

    }

    ~gsDPatchBase()
    {
        freeAll(m_bases);
    }

public:
//----------------------------------------------------------------------------------------------------------------------------
    /**
     * @brief      Computes the construction
     */
    virtual void compute();

    /**
     * @brief      Sets the default options
     */
    virtual void defaultOptions()
    {
        GISMO_NO_IMPLEMENTATION;
    }
    // {
    //     GISMO_NO_IMPLEMENTATION;
    // }

    virtual gsOptionList & options() { return m_options; }


//----------------------------------------------------------------------------------------------------------------------------
// Getters for basis, geometry and map
    /**
     * @brief       Returns the basis that is used for the D-Patch. Could be THB refined.
     *
     */
    virtual gsMultiBasis<T> localBasis() const {return m_bases;}

    /**
     * @brief       Returns the multipatch that is used for the D-Patch
     *
     */
    virtual gsMultiPatch<T> getGeometry() const {return m_patches;}

    /**
     * @brief       Exports a single modified patch with index \a patch
     *
     * The patch is obtained by transforming the coefficients of the D-Patch to the original basis, such that the original basis functions can be used to plot the geometry (and the patch structure will remain intact).
     * To construct the geometry, the coefficients for the C1 basis are multiplied with the transpose of the transformation matrix. The C1 coefficients are obtained with \ref _preCoefficients().
     *
     */
    virtual gsGeometry<T>* exportPatch(index_t patch, bool computeCoefs=true);

    /**
     * @brief      Exports the modified geometry to a @a gsMultiPatch object
     *
     * @return     A multipatch with the geometry
     */
    virtual gsMultiPatch<T> exportToPatches()
    {
        m_coefs = this->_preCoefficients();
        m_coefs = m_matrix.transpose() * m_coefs;

        std::vector<gsGeometry<T> *> patches(m_patches.nPatches());
        for (size_t p=0; p!=m_patches.nPatches(); p++)
            patches[p]= this->exportPatch(p,false);

        return gsMultiPatch<T>(patches,m_patches.boundaries(),m_patches.interfaces());
    }

    /**
     * @brief      Returns the smoothing matrix into \a matrix
     *
     * @param      matrix  The matrix
     */
    virtual const void matrix_into(gsSparseMatrix<T> & matrix) const
    { matrix = m_matrix; }

    /**
     * @brief      Returns the smoothing matrix
     *
     * The number of columns of the matrix corresponds to the number of basis functions in the local basis; this is the sum of all the basis functions over all the patches.
     * The number of rows of the matrix corresponds to the number of global basis functions, i.e. the number of basis functions corresponding to the D-Patch.
     * Multiplying the basis with the local basis function values gives the values of the basis functions in the global basis.
     *
     * @return     A matrix \a result to transfer local to global coefficients
     */
    virtual const gsSparseMatrix<T> matrix() const
    {
        gsSparseMatrix<T> result; matrix_into(result);
        return result;
    }

//----------------------------------------------------------------------------------------------------------------------------
// Info functions
    /**
     * @brief       Returns for each basis function if it is free or eliminated
     *
     * Returns for each basis function if it is free or eliminated and checks if the internal mapper is defined correctly
     */
    virtual void mapperInfo() const;

    /**
     * @brief      Returns information about a vertex
     *
     * @param[in]  corner  The \ref patchCorner
     *
     * @return     Prints the patch number, the corner index, the valence and if the corner is an interior or a boundary vertex
     */
    virtual void vertexInfo(patchCorner corner) const;

    /**
     * @brief      Returns information about a vertex
     *
     * @param[in]  patch  The \ref patchSide
     *
     * @return     Prints the patch number, the side index, the valence and if the side is a boundary or an interface
     */
    virtual void sideInfo(patchSide side) const;

    /**
     * @brief       Returns information for all the sides in the topology.
     *
     * Returns for all the patches and for all sides (1 to 4) if it is a boundary or an interface.
     */
    virtual void sideInfo() const;

    /**
     * @brief       Returns information for all the corners in the topology.
     *
     * Returns for all the patches and for all corners (1 to 4) the valence and if it is an interior vertex or a boundary vertex.
     */
    virtual void cornerInfo() const;


protected:
//----------------------------------------------------------------------------------------------------------------------------
// Pure virtual functions (to be overloaded)
    /**
     * @brief       Computes the C1 coefficients for pre-multiplication to make the multipatch
     *
     * Takes the coefficients which are tagged as "free" in the modified DoFMapper (m_mapModified) and when a boundary vertex with valence=3 is present, this one is shifted.
     *
     */
    virtual gsMatrix<T> _preCoefficients() = 0;

    /**
     * @brief      Initializes the class:
     *              -
     *              -
     */
    virtual void _initialize();


    /**
     * @brief      Calculates the mapper.
     */
    virtual void _computeMapper() = 0;

    /**
     * @brief      Calculates the smooth matrix.
     */
    virtual void _computeSmoothMatrix() = 0;

    /**
     * @brief      Makes a THB basis.
     */
    virtual void _makeTHB() = 0;

    /**
     * @brief      Corrects the EVs
     */
    virtual void _computeEVs() = 0;

//----------------------------------------------------------------------------------------------------------------------------
// Virtual functions (could be overloaded)
    /**
     * @brief       Computes the local coefficients and puts them in one big matrix
     */
    virtual gsMatrix<T> allCoefficients() const;

protected:
//----------------------------------------------------------------------------------------------------------------------------
// Helper functions
    /**
     * @brief      Checks if corners are sharp or not
     *
     * @param[in]  tol   The tolerance
     *
     * @return     A vector with a boolean true if the corner is C0 (following in-class corner indexing convention)
     */
    virtual std::vector<bool> getSharpCorners(T tol = 1e-2) const;

    /**
     * @brief      Computes the index of a basis function using sides as reference
     *
     * @param[in]  index1  The index of the basis function parallel to the first side
     * @param[in]  side1   The first side
     * @param[in]  index2  The index of the basis function parallel to the second side
     * @param[in]  side2   The second side
     *
     * @return     Index that is \a index1 in direction of \a side1 and \a index2 in direction of \a side2
     */
    virtual const index_t _indexFromSides(index_t index1, const patchSide side1, index_t index2, const patchSide side2);


    /**
     * @brief      Computes the index of a basis function taking one corner and one side as reference
     *
     * @param[in]  index   Offset of the basis function parallel to the side \a side, measured from \a corner
     * @param[in]  corner  The corner to be measured from
     * @param[in]  side    The side which contains \a corner
     * @param[in]  offset  The offset from the side (orthogonal to the side)
     *
     * @return     Index of \a index places from \a corner along \a side, with offset \a offset
     */
    virtual const gsVector<index_t> _indicesFromVert(index_t index, const patchCorner corner, const patchSide side, index_t offset = 0);


    /**
     * @brief      Computes the index of a basis function taking one corner and one side as reference
     *
     * @param[in]  bases   (optional) Multibasis to evaluate the index on
     * @param[in]  index   Offset of the basis function parallel to the side \a side, measured from \a corner
     * @param[in]  corner  The corner to be measured from
     * @param[in]  side    The side which contains \a corner
     * @param[in]  offset  The offset from the side (orthogonal to the side)
     * @param[in]  levelOffset  The level to be computed from. \a levelOffset = 0 returns the deepest THB level, and any positive number will give one level coarser
     *
     * @return     Index of \a index places from \a corner along \a side, with offset \a offset and with offset \a levelOffset from the deepest level
     */
    virtual const index_t _indexFromVert(gsMultiBasis<T> bases, index_t index, const patchCorner corner, const patchSide side, index_t offset = 0, index_t levelOffset = 0) const;
    virtual const index_t _indexFromVert(index_t index, const patchCorner corner, const patchSide side, index_t offset = 0, index_t levelOffset = 0) const;


    /**
     * @brief      Computes the index of a basis function taking one corner and one side as reference (multiple indices)
     *
     * @param[in]  bases   (optional) Multibasis to evaluate the index on
     * @param[in]  index        Vector with offsets of the basis function parallel to the side \a side, measured from \a corner
     * @param[in]  corner       The corner
     * @param[in]  side         The side
     * @param[in]  offset       The offset
     * @param[in]  levelOffset  The level offset
     *
     * @return     { description_of_the_return_value }
     */
    virtual const std::vector<index_t> _indexFromVert(gsMultiBasis<T> bases, std::vector<index_t> index, const patchCorner corner, const patchSide side, index_t offset = 0, index_t levelOffset = 0) const;
    virtual const std::vector<index_t> _indexFromVert(std::vector<index_t> index, const patchCorner corner, const patchSide side, index_t offset = 0, index_t levelOffset = 0) const;

    /**
     * @brief      Returns the valence and whether a corner is interior or boundary
     *
     * @param[in]  corner  The \ref patchCorner
     *
     * @return     A pair with .first giving the valence and .second being true if the vertex is interior and false if the vertex is on a boundary
     */
    virtual const std::pair<index_t,bool> _vertexData(const patchCorner corner) const;

    /**
     * @brief      Gets the valence.
     *
     * @param[in]  corner  The corner
     *
     * @return     The valence.
     */
    virtual const index_t _getValence( patchCorner corner) const
    { return this->_vertexData(corner).first; }

    /**
     * @brief      Determines whether the specified corner is interior vertex.
     *
     * @param[in]  corner  The corner
     *
     * @return     True if the specified corner is interior vertex, False otherwise.
     */
    virtual const bool _isInteriorVertex( patchCorner corner) const
    { return this->_vertexData(corner).second; }

    /**
     * @brief      Computes global index of the side
     *
     * @param[in]  patch    The patch number
     * @param[in]  bside    The \ref boxSide
     *
     * @return     Returns a global index of the side
     */
    virtual const index_t _sideIndex( index_t patch,  boxSide bside)     const
    { return 4*patch + bside - 1; }
    /**
     * @brief      Computes global index of the side
     *
     * @param[in]  pside    The \ref patchSide
     *
     * @return     Returns a global index of the side
     */
    virtual const index_t _sideIndex( patchSide pside)     const
    { return _sideIndex( pside.patch , pside.side() ); }

    /**
     * @brief      Computes global index of the corner
     *
     * @param[in]  patch    The patch number
     * @param[in]  corner   The \ref boxCorner
     *
     * @return     Returns a global index of the corner
     */
    virtual const index_t _vertIndex( index_t patch,  boxCorner corner)  const
    { return 4*patch + corner -1; }

    /**
     * @brief      Computes global index of the corner
     *
     * @param[in]  pcorner   The \ref patchCorner
     *
     * @return     Returns a global index of the side
     */
    virtual const index_t _vertIndex( patchCorner pcorner)     const
    { return _vertIndex( pcorner.patch , pcorner.corner() ); }


    /**
     * @brief      From a list of \a patchCorners pcorners, get the lowest \a n corners
     *
     * @param      pcorners  The pcorners
     * @param[in]  n         The number of corners
     */
    virtual void _getLowestCorners(std::vector<patchCorner> & pcorners, index_t n = 3) const;

    /**
     * @brief      From a list of \a patchCorners pcorners, remove all but the lowest \a n corners
     *
     * @param      pcorners  The pcorners
     * @param[in]  n         The number of corners
     */
    virtual void _removeLowestCorners(std::vector<patchCorner> & pcorners, index_t n = 3) const;

    /**
     * @brief      From a list of tuples (patch,index), get the lowest \a n tuples
     *
     * @param      pcorners  The pcorners
     * @param[in]  n         The number of corners
     */
    virtual void _getLowestIndices(std::vector<std::pair<index_t,index_t>> & indices, index_t n = 3) const;

    /**
     * @brief      From a list of tuples (patch,index), remove all but the lowest \a n tuples
     *
     * @param      pcorners  The pcorners
     * @param[in]  n         The number of corners
     */
    virtual void _removeLowestIndices(std::vector<std::pair<index_t,index_t>> & indices, index_t n = 3) const;


    virtual std::vector<std::pair<index_t,index_t>> _getInterfaceIndices(patchCorner pcorner, index_t depth, const gsMultiBasis<T> & mbasis) const;
    virtual std::vector<std::pair<index_t,index_t>> _getAllInterfaceIndices(patchCorner pcorner, index_t depth, const gsMultiBasis<T> & mbasis) const;


    virtual bool _checkMatrix(const gsSparseMatrix<T> & matrix) const;

// protected:
//     /**
//      * @brief      Prepares the THB basis if needed.
//      *
//      * This function constructs THB refinements on the places where they are needed, i.e. around EVs. It also constructs the transfer matrix (m_tMatrix) forms the transformation between the original B-spline basis and the THB-Spline basis.
//      */
//     void _makeTHB();

//     /**
//      * @brief      Computes D-Patch smoothing
//      *
//      * Given a basis with THB refinement around the EVs, this function computes the D-Patch smoothing
//      */
//     void _computeEVs();

//     /**
//      * @brief      Makes the Pi matrix
//      *
//      * This matrix is used to transform the coefficients of the D-Patch smoothing matrix
//      *
//      * @param[in]  valence  The valence
//      *
//      * @return     Matrix for smoothing around an EV}
//      */
//     gsMatrix<T> _makePi(index_t valence);

    /**
     * @brief      Initializes the matrix, the basis and the mappers
     */
    virtual void _initChecks();

    /**
     * @brief      Initializes the matrix, the basis and the mappers
     */
    virtual void _initTHB();

    /**
     * @brief      Initializes the basis.
     */
    virtual void _initBasis();

    /**
     *
     * @brief      Initializes the matrix, the basis and the mappers
     */
    virtual void _initMappers();

    /**
     * @brief      Initializes the matrix, the basis and the mappers
     */
    virtual void _countDoFs()
    {
        GISMO_NO_IMPLEMENTATION;
    }

    /**
     * @brief      Initializes the matrix.
     */
    virtual void _initMatrix();

    /**
     * @brief      Initializes the coefficients.
     */
    virtual void _initCoefs();

    /**
     * @brief      Performs checks on sides, vertices and bases
     */
    virtual void _performChecks(bool basis);

    /**
     * @brief      Resets checks on sides, vertices and bases
     */
    virtual void _resetChecks(bool basis);


//     *
//      * @brief      Computes the modified mapper
//      *
//      * The modified mapper is computed based on the elimination of different functions with different conditions.
//      * 1) For interfaces, it eliminates all the nodes except the first two and the last two
//      * 2) For boundaries, there is no elimination
//      * 3) For vertices, there are few options
//      *  a) Boundary vertices
//      *      i)  Valence 1: No eliminations
//      *      ii) Valence 2: the two outer-most basis functions on the interface (the one at the vertex and the one next to it on the interface) are both eliminated
//      *      iii)Valence 3: On all the patches, the basis functions corresponding to the vertex are coupled to eachother. The basis functions next to this one (on an interface OR on a boundary) are eliminated
//      *  b) Interior vertices: all basis functions along the interface are eliminated if not done so

//     void _computeMapper(); // also initialize the mappers!

// protected:

//     void _pushToMatrix(tuple_t entries)
//     {
//         index_t rowIdx,colIdx;
//         T weight;
//         for (typename tuple_t::const_iterator it=entries.begin(); it!=entries.end(); it++)
//         {
//             std::tie(rowIdx,colIdx,weight) = *it;
//             m_matrix(rowIdx,colIdx) = weight;
//             m_basisCheck[rowIdx] = true;
//         }
//     }

//     void _computeInterfaceMapper(boundaryInterface iface);

//     void _computeBoundaryMapper(patchSide boundary);

//     void _computeVertexMapper(patchCorner pcorner);

//     // Boundary vertex of valence 1
//     template<bool _boundary, index_t _v> // valence=2
//     typename std::enable_if<  _boundary && _v==1, void>::type
//     _computeVertexMapper_impl(patchCorner pcorner, index_t valence);

//     // Boundary vertex of valence 2 with C1 smoothness
//     template<bool _boundary, index_t _v, bool _smooth> // valence=2
//     typename std::enable_if<  _boundary && _v==2 && _smooth, void>::type
//     _computeVertexMapper_impl(patchCorner pcorner, index_t valence);

//     // Boundary vertex of valence 2 with C0 smoothness
//     template<bool _boundary, index_t _v, bool _smooth> // valence=2
//     typename std::enable_if<  _boundary && _v==2 && (!_smooth), void>::type
//     _computeVertexMapper_impl(patchCorner pcorner, index_t valence);

//     // Boundary vertex of valence !(1,2,3) with C1 smoothness
//     template<bool _boundary, index_t _v, bool _smooth>
//     typename std::enable_if<  _boundary && _v==-1 && _smooth, void>::type
//     _computeVertexMapper_impl(patchCorner pcorner, index_t valence);

//     // Boundary vertex of valence !(1,2,3) with C0 smoothness
//     template<bool _boundary, index_t _v, bool _smooth>
//     typename std::enable_if<  _boundary && _v==-1 && (!_smooth), void>::type
//     _computeVertexMapper_impl(patchCorner pcorner, index_t valence);

//     // Ordinary interior vertex
//     template<bool _boundary, index_t _v> // valence=2
//     typename std::enable_if<  (!_boundary) && _v==4, void>::type
//     _computeVertexMapper_impl(patchCorner pcorner, index_t valence);

//     // Extraordinary interior vertex
//     template<bool _boundary, index_t _v>
//     typename std::enable_if<  (!_boundary) && _v==-1, void>::type
//     _computeVertexMapper_impl(patchCorner pcorner, index_t valence);

// public:

//     /**
//      * @brief      Handles a vertex in the global matrix
//      *
//      * We use the following notation convention (per patch!):
//      * b00 is the basis function at the vertex
//      * b10 is the basis function next to the vertex along the first interface that connects to the vertex
//      * b20 is the basis function next to b10 along the first interface that connects to the vertex
//      * etc.
//      *
//      * b01 is the basis function next to the vertex along the second interface that connects to the vertex
//      * b02 is the basis function next to b01 along the second interface that connects to the vertex
//      * etc.
//      *
//      * b11 is the basis function with offset 1 from both interfaces and from the vertex itself
//      * b22 is the basis function with offset 2 from both interfaces and from the vertex itself
//      * etc.
//      *
//      * There are different options.
//      * a) Boundary vertices
//      *      i)  Valence 1: b00, b10, b01 and b00 all get weight 1.0 w.r.t the same basis function in the local basis
//      *      ii) Valence 2: This case contains an interface between two patches. We use index k to denote the row basis functions along the interface. So k=0 corresponds to the basis functions on the boundary and k=1 corresponds to the basis functions with offset 1 from the boundaries. Using this convention, the functions bk1 in the local basis, are coupled to bk1 in the global basis with weight 1. The functions bk0 in the local basis (on the considered patch) are coupled to bk1 in the global basis with weight 0.5. The functions bk0 in the local basis (on the other patch) are coupled to bk1 (on the considered patch) in the global basis with weight 0.5.
//      *      iii)Valence 3: In this case, the matched vertices on all the adjacent patches are treated in one go! Note that all the basis functions corresponding to the vertex (b00) on all patches are matched! We couple the b00 functions of all patches (in the local basis) with weight 1/4 to the b00 of the adjacent patch with the lowest number in the global basis. Then, the b11 on the considered patch is coupled with weight 1 to itself and with weight 0.25 to the b00s of the other patches. Then, we will handle the vertices where an interface and a boundary meet (there are two of these). For the patch corners that are on an interface, we find the b11 and b10 vertices (orthogonal to the interface) and we give all b10s weight 0.5 w.r.t. the b11s in the global basis (on both patches). Lastly, we add weight 0.5 for the b10s along the boundaries (so only for two patches) to the (matched) b00 basis function (all b00s refer to the same dof in the global basis).
//      * b) Interior vertices (all valences):
//      *      i)  b11 gets weight 1.0 w.r.t the same basis function in the local basis
//      *      ii) all associated b00s in the local basis get weight 1/valence w.r.t. b11 in the global basis
//      *      iii)the b10 and b01 in the local basis get weight 1/2 w.r.t. b11 in the global basis
//      *
//      * @param[in]  pcorner  The patchcorner
//      */
//     void _handleVertex(patchCorner pcorner);

// protected:

//     void _handleRegularCorner(patchCorner pcorner);

//     template<bool _regular, bool _smooth> // valence=2
//     typename std::enable_if<  _regular  &&   _smooth   , void>::type
//     _handleBoundaryVertex(patchCorner pcorner, index_t valence);

//     template<bool _regular, bool _smooth> // valence=2
//     typename std::enable_if<  _regular  && (!_smooth)   , void>::type
//     _handleBoundaryVertex(patchCorner pcorner, index_t valence);

//     template<bool _regular, bool _smooth> // valence > 2
//     typename std::enable_if<(!_regular) &&   _smooth    , void>::type
//     _handleBoundaryVertex(patchCorner pcorner, index_t valence);

//     template<bool _regular, bool _smooth> // valence > 1
//     typename std::enable_if<(!_regular) && (!_smooth)   , void>::type
//     _handleBoundaryVertex(patchCorner pcorner, index_t valence);

//     // interior vertices
//     void _handleInteriorVertex(patchCorner pcorner, index_t valence);

// public:

//     /**
//      * @brief      Handles an interface in the global matrix
//      *
//      * Gives all the DoFs that have offset 1 (orthogonal) from the interface weight 1.0 w.r.t itself. All the DoFs ON the interface (on both patches) will have weight 0.5 to the DoF with offset 1.
//      * This interface handling excludes the indices that are in the 0 and 1 ring around vertices.
//      *
//      * @param[in]  iface  The interface
//      */
//     void _handleInterface(boundaryInterface iface);
//     /**
//      * @brief      Handles a boundary in the global matrix
//      *
//      * Handles all DoFs on the boundary with unit-weight, except the ones in the 0 and 1 rings around the vertices.
//      *
//      * @param[in]  side  The boundary side
//      */
//     void _handleBoundary(patchSide side);
//     /**
//      * @brief      Handles the interior in the global matrix
//      *
//      * Gives all left-over DoFs, which are in the interior, weight 1 w.r.t. itself
//      */
//     void _handleInterior();
    /**
     * @brief      Prints which DoFs have been handled and which have been eliminated
     */
    virtual void _whichHandled();



// protected:
//     bool _checkMatrix(const gsSparseMatrix<T> & matrix) const // ! makes a deep copy (otherwise the contents of m_matrix get destroyed somehow...)
//     {
//         GISMO_ASSERT(matrix.cols()==matrix.outerSize(),"is the matrix ColMajor?");
//         gsVector<T> colSums(matrix.cols());
//         colSums.setZero();
//         for (index_t i = 0; i<matrix.outerSize(); ++i)
//             for (typename gsSparseMatrix<T>::iterator it = matrix.begin(i); it; ++it)
//                 colSums.at(i) += it.value();

//         return (colSums.array() < 1+1e-8).any() && (colSums.array() > 1-1e-8).any();
//     }

protected:
    const gsMultiPatch<T> & m_patches;
    gsMultiPatch<T> m_RefPatches;
    gsMultiBasis<T> m_bases, m_Bbases;
    mutable gsSparseMatrix<T> m_tMatrix;
    mutable std::vector<bool> m_sideCheck;
    mutable std::vector<bool> m_vertCheck;
    mutable std::vector<bool> m_basisCheck;
    mutable std::vector<bool> m_C0s;

    mutable gsDofMapper m_mapModified,m_mapOriginal;

    mutable gsSparseMatrix<T> m_matrix;

    mutable size_t m_size;

    mutable gsMatrix<T> m_coefs;

    gsOptionList m_options;

    size_t m_nSides;
    size_t m_nVerts;

};



}

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsDPatchBase.hpp)
#endif
