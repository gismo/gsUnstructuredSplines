/** @file gsMPBESSpline_test.cpp

    @brief Regression tests for the gsMPBESSpline ownership fixes: the
    multipatch constructors must populate every piece, and assigning into an
    already-populated gsMPBESSpline must preserve the derived (AS-G1) basis
    type rather than slicing it to the gsMappedBasis base.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Verhelst
**/

#include "gismo_unittest.h"
#include <gsUnstructuredSplines/src/gsMPBESSpline.h>
#include <gsUnstructuredSplines/src/gsMPBESBasis.h>

using namespace gismo;

namespace
{

// Two unit squares side by side, degree-elevated to 3 so that
// incrSmoothness=3 (used throughout this file, matching
// create_multipatch.cpp's method==0 path) is admissible.
gsMultiPatch<real_t> twoPatchGrid()
{
    gsMultiPatch<real_t> mp = gsNurbsCreator<real_t>::BSplineSquareGrid(2, 1, 1.0);
    mp.degreeElevate(2);
    return mp;
}

}

SUITE(gsMPBESSpline_test)
{

    TEST(MultipatchCtorsPopulatePieces)
    {
        gsMultiPatch<real_t> mp = twoPatchGrid();

        gsMPBESSpline<2,real_t> cgeom(mp, 3);
        CHECK_EQUAL(2, cgeom.nPieces());
        gsMatrix<real_t> pt(2,1); pt << 0.5, 0.5;
        cgeom.piece(0).eval(pt);
        cgeom.piece(1).eval(pt);

        std::vector<patchCorner> c0List; // no forced C0 corners
        gsMPBESSpline<2,real_t> cgeomC0(mp, c0List, 3);
        CHECK_EQUAL(2, cgeomC0.nPieces());
        cgeomC0.piece(0).eval(pt);
        cgeomC0.piece(1).eval(pt);
    }

    // The end-to-end motivating case: gsMPBESBasis::clone() is a pure
    // virtual override (GISMO_UPTR_FUNCTION_PURE in gsMPBESBasis.h) so that
    // gsMappedSpline::operator=, calling `other.m_mbases->clone()`
    // polymorphically, produces another gsMPBESBasis. If that override were
    // ever dropped, clone() would resolve to gsMappedBasis's own
    // GISMO_CLONE_FUNCTION, which constructs a plain gsMappedBasis via its
    // own copy constructor regardless of the runtime type — silently
    // slicing the AS-G1 basis away.
    TEST(AssignmentAfterPopulatedCtor)
    {
        gsMultiPatch<real_t> mp = twoPatchGrid();
        gsMPBESSpline<2,real_t> sp(mp, 3);

        gsMultiPatch<real_t> mp2 = twoPatchGrid();
        gsMPBESSpline<2,real_t> cp(mp2, 3); // already-populated target

        cp = sp;

        CHECK_EQUAL(sp.nPieces(), cp.nPieces());
        gsMatrix<real_t> pt(2,1); pt << 0.3, 0.6;
        cp.piece(0).eval(pt);

        const gsMPBESBasis<2,real_t> * derived = dynamic_cast<const gsMPBESBasis<2,real_t>*>(&cp.getMappedBasis());
        CHECK(nullptr != derived);
    }

}
