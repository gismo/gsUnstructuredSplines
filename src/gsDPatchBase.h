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

#include <gsIO/gsOptionList.h>

#include <gsMSplines/gsMappedBasis.h>
#include <gsMSplines/gsMappedSpline.h>

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
    typedef typename std::vector<std::tuple<index_t,index_t,T>> sparseEntry_t;

public:

    /// Shared pointer for gsDPatchBase
    typedef memory::shared_ptr< gsDPatchBase > Ptr;

    /// Unique pointer for gsDPatchBase
    typedef memory::unique_ptr< gsDPatchBase > uPtr;

    /// Empty constructor
    gsDPatchBase()
    :
    gsDPatchBase(gsMultiBasis<T>(),gsMultiPatch<T>())
    { }

    /**
     * @brief      Default constructor
     *
     * @param      mn    MultiBasis
     */
    gsDPatchBase(const gsMultiBasis<T> & mb, const gsMultiPatch<T> & mp)
    :
    m_patches(mp),
    m_Bbases(mb),
    m_topology(m_Bbases.topology()),
    m_computed(false)             
    {
    }

    gsDPatchBase(const gsMultiBasis<T> & mb)
    :
    gsDPatchBase(mb, gsMultiPatch<T>())
    {
    }

    /**
     * @brief      Default constructor
     *
     * @param      mp    Multipatch of the geometry
     */
    gsDPatchBase(const gsMultiPatch<T> & mp)
    :
    gsDPatchBase(give(gsMultiBasis<T>(mp)), mp)
    {
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
    virtual void defaultOptions();
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
    virtual const gsMultiBasis<T> & localBasis() const
    {
        return m_bases;
    }

    virtual void localBasis_into(gsMultiBasis<T> & localBasis) const
    {
        localBasis = m_bases;
    }

    /**
     * @brief       Returns the modified geometry  corresponding to the local basis
     *
     */
    virtual void localGeometry_into(gsMultiPatch<T> & localGeometry)
    {
        localGeometry = exportToPatches(m_patches);
    }

    /**
     * @brief       Returns the multipatch that is used for the D-Patch
     *
     */
    virtual const gsMultiPatch<T> & getGeometry() const
    {
        return m_patches;
    }

    /**
     * @brief       Returns the basis that is used for the D-Patch. Could be THB refined.
     *
     */
    virtual void globalBasis_into(gsMappedBasis<d,T> & mbasis) const
    {
        GISMO_ASSERT(m_computed,"The method has not been computed! Call compute().");
        gsSparseMatrix<T> matrix = m_matrix.transpose();
        mbasis.init(m_bases,matrix);
        gsDebugVar("hi");
    }

    /**
     * @brief       Returns the multipatch that is used for the D-Patch
     *
     */
    virtual void globalGeometry_into(const gsMultiPatch<T> & patches, gsMappedSpline<d,T> & mspline)
    {
        GISMO_ASSERT(!patches.empty() && patches.nPatches()==m_bases.nBases(),"The reference multipatch is empty!");
        gsMappedBasis<d,T> mbasis;
        this->globalBasis_into(mbasis);
        gsMatrix<T> localCoefs = this->_preCoefficients(patches);
        mspline.init(mbasis,localCoefs);
    }

    virtual void globalGeometry_into(gsMappedSpline<d,T> & mspline)
    {
        this->globalGeometry_into(m_patches,mspline);
    }

    virtual void update( gsMappedBasis<d,T> & mbasis ) const
    {
        mbasis.init(m_bases,m_matrix.transpose());
    }

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
    virtual gsMultiPatch<T> exportToPatches(const gsMultiPatch<T> & patches)
    {
        GISMO_ASSERT(m_computed,"The method has not been computed! Call compute().");
        GISMO_ASSERT(!patches.empty(),"The reference multipatch is empty!");
        m_coefs = this->_preCoefficients(patches);
        m_coefs = m_matrix.transpose() * m_coefs;

        std::vector<gsGeometry<T> *> PatchContainer(patches.nPatches());
        for (size_t p=0; p!=patches.nPatches(); p++)
            PatchContainer[p]= this->exportPatch(p,false);

        return gsMultiPatch<T>(PatchContainer,m_topology.boundaries(),m_topology.interfaces());
    }

    virtual gsMultiPatch<T> exportToPatches()
    {
        GISMO_ASSERT(!m_patches.empty(),"The reference multipatch is empty!");
        return this->exportToPatches(m_patches);
    }


    /**
     * @brief      Returns the smoothing matrix into \a matrix
     *
     * @param      matrix  The matrix
     */
    virtual const void matrix_into(gsSparseMatrix<T> & matrix) const
    {
        matrix = this->matrix();
    }

    virtual const gsSparseMatrix<T> & Tmatrix() const
    {
        return m_tMatrix;
    }

    /**
     * @brief      Returns the smoothing matrix
     *
     * The number of columns of the matrix corresponds to the number of basis functions in the local basis; this is the sum of all the basis functions over all the patches.
     * The number of rows of the matrix corresponds to the number of global basis functions, i.e. the number of basis functions corresponding to the D-Patch.
     * Multiplying the basis with the local basis function values gives the values of the basis functions in the global basis.
     *
     * @return     A matrix \a result to transfer local to global coefficients
     */
    virtual const gsSparseMatrix<T> & matrix() const
    {
        GISMO_ASSERT(m_computed,"The method has not been computed! Call compute().");
        return m_matrix;
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


    gsMatrix<T> preCoefficients() { return _preCoefficients();};

protected:
//----------------------------------------------------------------------------------------------------------------------------
// Pure virtual functions (to be overloaded)
    /**
     * @brief       Computes the C1 coefficients for pre-multiplication to make the multipatch
     *
     * Takes the coefficients which are tagged as "free" in the modified DoFMapper (m_mapModified) and when a boundary vertex with valence=3 is present, this one is shifted.
     *
     */
    virtual gsMatrix<T> _preCoefficients(const gsMultiPatch<T> & patches) = 0;
    virtual gsMatrix<T> _preCoefficients()
    {
        return _preCoefficients(m_patches);
    }

    /**
     * @brief      Initializes the class:
     *              -
     *              -
     */
    virtual void _initialize();


    /**
     * @brief      Calculates the mapper.
     */
    virtual void _computeMapper();

    virtual void _computeInterfaceMapper(boundaryInterface iface);

    virtual void _computeBoundaryMapper(patchSide boundary);

    virtual void _computeVertexMapper(patchCorner pcorner);

    /**
     * @brief      Calculates the smooth matrix.
     */
    virtual void _computeSmoothMatrix();

    void _push(sparseEntry_t entries)
    {
        index_t rowIdx,colIdx;
        T weight;
        for (typename sparseEntry_t::const_iterator it=entries.begin(); it!=entries.end(); it++)
        {
            std::tie(rowIdx,colIdx,weight) = *it;
            m_matrix(rowIdx,colIdx) = weight;
        }
    }

    void _pushAndCheck(sparseEntry_t entries)
    {
        index_t rowIdx,colIdx;
        T weight;
        for (typename sparseEntry_t::const_iterator it=entries.begin(); it!=entries.end(); it++)
        {
            std::tie(rowIdx,colIdx,weight) = *it;
            m_matrix(rowIdx,colIdx) = weight;
            m_basisCheck[rowIdx] = true;
        }
    }

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
    virtual void _handleVertex(patchCorner pcorner);
    // interior vertices
    virtual void _handleInteriorVertex(patchCorner pcorner, index_t valence);

    /**
     * @brief      Handles an interface in the global matrix
     *
     * Gives all the DoFs that have offset 1 (orthogonal) from the interface weight 1.0 w.r.t itself. All the DoFs ON the interface (on both patches) will have weight 0.5 to the DoF with offset 1.
     * This interface handling excludes the indices that are in the 0 and 1 ring around vertices.
     *
     * @param[in]  iface  The interface
     */
    virtual void _handleInterface(boundaryInterface iface);
    /**
     * @brief      Handles a boundary in the global matrix
     *
     * Handles all DoFs on the boundary with unit-weight, except the ones in the 0 and 1 rings around the vertices.
     *
     * @param[in]  side  The boundary side
     */
    virtual void _handleBoundary(patchSide side);
    /**
     * @brief      Handles the interior in the global matrix
     *
     * Gives all left-over DoFs, which are in the interior, weight 1 w.r.t. itself
     */
    virtual void _handleInterior();

protected:
    /**
     * @brief      Handles a regular corner
     *
     * @param[in]  pcorner  The pcorner
     */
    virtual void _handleRegularCorner(patchCorner pcorner);

    virtual void _handleRegularBoundaryVertexSmooth(patchCorner pcorner, index_t valence);

    virtual void _handleRegularBoundaryVertexNonSmooth(patchCorner pcorner, index_t valence);

    virtual void _handleIrregularBoundaryVertexSmooth(patchCorner pcorner, index_t valence);

    virtual void _handleIrregularBoundaryVertexNonSmooth(patchCorner pcorner, index_t valence);

    /**
     * @brief      Prints which DoFs have been handled and which have been eliminated
     */

protected:

    /**
     * @brief      Makes a THB basis.
     */
    virtual void _makeTHB() = 0;

    /**
     * @brief      Corrects the EVs
     */
    virtual void _computeEVs() = 0;

protected:
    // Boundary vertex of valence 1
    // template<bool _boundary, index_t _v> // valence=2
    // virtual typename std::enable_if<  _boundary && _v==1, void>::type
    // SAME
    virtual void _computeMapperRegularCorner_v1(patchCorner pcorner, index_t valence);

    // Boundary vertex of valence 2 with C1 smoothness
    // template<bool _boundary, index_t _v, bool _smooth> // valence=2
    // virtual  typename std::enable_if<  _boundary && _v==2 && _smooth, void>::type
    // SAME
    virtual void _computeMapperRegularBoundaryVertexSmooth_v2(patchCorner pcorner, index_t valence);

    // Boundary vertex of valence 2 with C0 smoothness
    // template<bool _boundary, index_t _v, bool _smooth> // valence=2
    // virtual  typename std::enable_if<  _boundary && _v==2 && (!_smooth), void>::type
    // DIFFERENT
    virtual void _computeMapperRegularBoundaryVertexNonSmooth_v2(patchCorner pcorner, index_t valence);

    // Boundary vertex of valence !(1,2,3) with C1 smoothness
    // template<bool _boundary, index_t _v, bool _smooth>
    // virtual  typename std::enable_if<  _boundary && _v==-1 && _smooth, void>::type
    // DIFFERENT
    virtual void _computeMapperIrregularBoundaryVertexSmooth_v(patchCorner pcorner, index_t valence);

    // Boundary vertex of valence !(1,2,3) with C0 smoothness
    // template<bool _boundary, index_t _v, bool _smooth>
    // virtual  typename std::enable_if<  _boundary && _v==-1 && (!_smooth), void>::type
    // DIFFERENT
    virtual void _computeMapperIrregularBoundaryVertexNonSmooth_v(patchCorner pcorner, index_t valence);

    // Interior vertex
    // template<bool _boundary, index_t _v>
    // virtual  typename std::enable_if<  (!_boundary) && _v==-1, void>::type
    // DIFFERENT
    virtual void _computeMapperInteriorVertex_v(patchCorner pcorner, index_t valence);


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
     * @param[in]  bases   (optional) Multibasis to evaluate the index on
     * @param[in]  index   Offset of the basis function parallel to the side \a side, measured from \a corner
     * @param[in]  corner  The corner to be measured from
     * @param[in]  side    The side which contains \a corner
     * @param[in]  offset  The offset from the side (orthogonal to the side)
     *
     * @return     Index of \a index places from \a corner along \a side, with offset \a offset and with offset \a levelOffset from the deepest level
     */
    virtual const index_t _indexFromVert(const index_t index, const patchCorner corner, const patchSide side, const index_t offsets = 0) const;
    virtual const index_t _indexFromVert(const gsMultiBasis<T> & bases, const index_t index, const patchCorner corner, const patchSide side, const index_t offsets = 0) const;
    virtual const index_t _indexFromVert(const gsBasis<T> * basis, const index_t index, const patchCorner corner, const patchSide side, const index_t offsets = 0) const;
private:
    template<class U>
    typename util::enable_if<util::is_same<U, const gsHTensorBasis<d,T> *>::value,const index_t>::type
    _indexFromVert_impl(U basis, const index_t index, const patchCorner corner, const patchSide side, const index_t offsets = 0) const;

    template<class U>
    typename util::enable_if<util::is_same<U, const gsTensorBSplineBasis<d,T> *>::value,const index_t>::type
    _indexFromVert_impl(U basis, const index_t index, const patchCorner corner, const patchSide side, const index_t offsets = 0) const;
protected:

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
     * @brief      Performs checks on sides, vertices and bases
     */
    virtual void _performChecks(bool basis);

    /**
     * @brief      Resets checks on sides, vertices and bases
     */
    virtual void _resetChecks(bool basis);

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
    const gsMultiBasis<T> m_Bbases; // reference?
    
    gsMultiBasis<T> m_bases;

    gsBoxTopology m_topology;

    bool m_computed;

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
