/** @file biharmonic_example.cpp

    @brief A Biharmonic example for a single patch.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/

# include <gismo.h>

#include <gsUnstructuredSplines/src/gsApproxC1Spline.h>
#include <gsUnstructuredSplines/src/gsDPatch.h>
#include <gsUnstructuredSplines/src/gsAlmostC1.h>
#include <gsUnstructuredSplines/src/gsC1SurfSpline.h>

/**
 * Smoothing method:
 * - m 0 == Approx C1 method
 * - m 1 == D-Patch method
 * - m 2 == Almost C1 method
 * - m 3 == Nitsche's method
 */
enum MethodFlags
{
    APPROXC1       = 0, // Approx C1 Method
    DPATCH         = 1, // D-Patch
    ALMOSTC1       = 2, // Almost C1
    NITSCHE        = 3, // Nitsche
    SPLINE         = 4, // Spline (only for single patch)
    SURFASG1       = 5, // Only for AS-G1 geometries
    // Add more [...]
};

namespace gismo{
    namespace expr{
        template<class E>
        class deriv2_expr : public _expr<deriv2_expr<E> >
        {
            typename E::Nested_t _u;

        public:
            // enum {ColBlocks = E::rowSpan }; // ????
            enum{ Space = E::Space, ScalarValued= 0, ColBlocks = (1==E::Space?1:0) };

            typedef typename E::Scalar Scalar;

            deriv2_expr(const E & u) : _u(u) { }

            mutable gsMatrix<Scalar> res, tmp;

            const gsMatrix<Scalar> & eval(const index_t k) const {return eval_impl(_u,k); }

            index_t rows() const //(components)
            {
                return 3; // _u.dim() for space or targetDim() for geometry
            }

            index_t cols() const
            {
                return _u.source().domainDim() * ( _u.source().domainDim() + 1 ) / 2;
            }

            void parse(gsExprHelper<Scalar> & evList) const
            {
                _u.parse(evList);
                _u.data().flags |= NEED_DERIV2;
            }

            const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
            const gsFeSpace<Scalar> & colVar() const { return _u.colVar(); }

            index_t cardinality_impl() const { return _u.cardinality_impl(); }

            void print(std::ostream &os) const { os << "deriv2("; _u.print(os); os <<")"; }

            private:
                template<class U> inline
                typename util::enable_if< util::is_same<U,gsGeometryMap<Scalar> >::value, const gsMatrix<Scalar> & >::type
                eval_impl(const U & u, const index_t k)  const
                {
                    /*
                        Here, we compute the hessian of the geometry map.
                        The hessian of the geometry map c has the form: hess(c)
                        [d11 c1, d11 c2, d11 c3]
                        [d22 c1, d22 c2, d22 c3]
                        [d12 c1, d12 c2, d12 c3]
                        The geometry map has components c=[c1,c2,c3]
                    */
                    // evaluate the geometry map of U
                    res = _u.data().values[2].reshapeCol(k, cols(), _u.data().dim.second );
                    return res;
                }

                template<class U> inline
                typename util::enable_if< util::is_same<U,gsComposition<Scalar> >::value, const gsMatrix<Scalar> & >::type
                eval_impl(const U & u, const index_t k)  const
                {
                    /*
                        Here, we compute the hessian of the geometry map.
                        The hessian of the geometry map c has the form: hess(c)
                        [d11 c1, d11 c2, d11 c3]
                        [d22 c1, d22 c2, d22 c3]
                        [d12 c1, d12 c2, d12 c3]
                        The geometry map has components c=[c1,c2,c3]
                    */
                    // evaluate the geometry map of U
                    tmp =  _u.data().values[2];
                    gsDebugVar(tmp);
                    // Cast to Voight notation
                    tmp.resize(4,1);
                    std::swap( tmp(3,0), tmp(1,0) );

                    res = tmp;
                    return res;
                }

                template<class U> inline
                typename util::enable_if< util::is_same<U,gsFeSolution<Scalar> >::value, const gsMatrix<Scalar> & >::type
                eval_impl(const U & u, const index_t k)  const
                {
                    /*
                        Here, we compute the hessian of the geometry map.
                        The hessian of the geometry map c has the form: hess(c)
                        [d11 c1, d11 c2, d11 c3]
                        [d22 c1, d22 c2, d22 c3]
                        [d12 c1, d12 c2, d12 c3]
                        The geometry map has components c=[c1,c2,c3]
                    */
                    // evaluate the geometry map of U
                    hess_expr<gsFeSolution<Scalar>> sHess = hess_expr<gsFeSolution<Scalar>>(_u);
                    res = sHess.eval(k).transpose();
                    return res;
                }

                /// Spexialization for a space
                template<class U> inline
                typename util::enable_if<util::is_same<U,gsFeSpace<Scalar> >::value, const gsMatrix<Scalar> & >::type
                eval_impl(const U & u, const index_t k) const
                {
                    /*
                        Here, we compute the hessian of the basis with n actives.
                        The hessian of the basis u has the form: hess(u)
                            active 1                active 2                        active n = cardinality
                        [d11 u1, d11 u2, d11 u3] [d11 u1, d11 u2, d11 u3] ... [d11 u1, d11 u2, d11 u3]
                        [d22 u1, d22 u2, d22 u3] [d22 u1, d22 u2, d22 u3] ... [d22 u1, d22 u2, d22 u3]
                        [d12 u1, d12 u2, d12 u3] [d12 u1, d12 u2, d12 u3] ... [d12 u1, d12 u2, d12 u3]
                        Here, the basis function has components u = [u1,u2,u3]. Since they are evaluated for scalars
                        we use blockDiag to make copies for all components ui
                            active 1     active 2     active k = cardinality/dim   active 1           active 2k       active 1           active 2k
                        [d11 u, 0, 0] [d11 u, 0, 0] ... [d11 u, 0, 0]            [0, d11 u, 0]  ... [0, d11 u, 0]  [0, d11 u, 0]  ... [0, d11 u, 0]
                        [d22 u, 0, 0] [d22 u, 0, 0] ... [d22 u, 0, 0]            [0, d22 u, 0]  ... [0, d22 u, 0]  [0, d22 u, 0]  ... [0, d22 u, 0]
                        [d12 u, 0, 0] [d12 u, 0, 0] ... [d12 u, 0, 0]            [0, d12 u, 0]  ... [0, d12 u, 0]  [0, d12 u, 0]  ... [0, d12 u, 0]
                    */
                    const index_t numAct = u.data().values[0].rows();   // number of actives of a basis function
                    const index_t cardinality = u.cardinality();        // total number of actives (=3*numAct)

                    res.resize(rows(), _u.dim() *_u.cardinality()); // (3 x 3*cardinality)
                    res.setZero();

                    tmp = _u.data().values[2].reshapeCol(k, cols(), numAct );
                    for (index_t d = 0; d != cols(); ++d)
                    {
                        const index_t s = d*(cardinality + 1);
                        for (index_t i = 0; i != numAct; ++i)
                            res.col(s+i*_u.cols()) = tmp.col(i);
                    }

                    return res;
                }
        };

        // vector v should be a row vector
        template<class E1, class E2>
        class deriv2dot_expr : public _expr<deriv2dot_expr<E1, E2> >
        {
            typename E1::Nested_t _u;
            typename E2::Nested_t _v;

        public:
            enum{ Space = E1::Space, ScalarValued= 0, ColBlocks= 0 };
            // Note: what happens if E2 is a space? The following can fix it:
            // enum{ Space = (E1::Space == 1 || E2::Space == 1) ? 1 : 0, ScalarValued= 0, ColBlocks= 0 };

            typedef typename E1::Scalar Scalar;

            deriv2dot_expr(const E1 & u, const E2 & v) : _u(u), _v(v) { }

            mutable gsMatrix<Scalar> res,tmp, vEv;

            const gsMatrix<Scalar> & eval(const index_t k) const {return eval_impl(_u,k); }

            index_t rows() const
            {
                return 1; //since the product with another vector is computed
            }

            index_t cols() const
            {
                return 1;
            }

            void parse(gsExprHelper<Scalar> & evList) const
            {
                parse_impl<E1>(evList);
            }

            const gsFeSpace<Scalar> & rowVar() const
            {
                if      (E1::Space == 1 && E2::Space == 0)
                    return _u.rowVar();
                else if (E1::Space == 0 && E2::Space == 1)
                    return _v.rowVar();
                else
                    return gsNullExpr<Scalar>::get();
            }

            const gsFeSpace<Scalar> & colVar() const
            {
                if      (E1::Space == 1 && E2::Space == 0)
                    return _v.colVar();
                else if (E1::Space == 0 && E2::Space == 1)
                    return _u.colVar();
                else
                    return gsNullExpr<Scalar>::get();
            }

            void print(std::ostream &os) const { os << "deriv2("; _u.print(os); _v.print(os); os <<")"; }

        private:
            template<class U> inline
            typename util::enable_if< !util::is_same<U,gsFeSolution<Scalar> >::value,void>::type
            parse_impl(gsExprHelper<Scalar> & evList) const
            {
                evList.add(_u);   // We manage the flags of _u "manually" here (sets data)
                _u.data().flags |= NEED_DERIV2; // define flags

                _v.parse(evList); // We need to evaluate _v (_v.eval(.) is called)

                // Note: evList.parse(.) is called only in exprAssembler for the global expression
            }

            template<class U> inline
            typename util::enable_if< util::is_same<U,gsFeSolution<Scalar> >::value,void>::type
            parse_impl(gsExprHelper<Scalar> & evList) const
            {
                _u.parse(evList); //
                hess(_u).parse(evList); //

                // evList.add(_u);   // We manage the flags of _u "manually" here (sets data)
                _v.parse(evList); // We need to evaluate _v (_v.eval(.) is called)
            }

            template<class U> inline
            typename util::enable_if< util::is_same<U,gsGeometryMap<Scalar> >::value, const gsMatrix<Scalar> & >::type
            eval_impl(const U & u, const index_t k)  const
            {
                /*
                    Here, we multiply the hessian of the geometry map by a vector, which possibly has multiple actives.
                    The hessian of the geometry map c has the form: hess(c)
                    [d11 c1, d11 c2, d11 c3]
                    [d22 c1, d22 c2, d22 c3]
                    [d12 c1, d12 c2, d12 c3]
                    And we want to compute [d11 c .v; d22 c .v;  d12 c .v] ( . denotes a dot product and c and v are both vectors)
                    So we simply evaluate for every active basis function v_k the product hess(c).v_k
                */
                // evaluate the geometry map of U
                tmp =_u.data().values[2].reshapeCol(k, cols(), _u.data().dim.second );
                vEv = _v.eval(k);
                res = vEv * tmp.transpose();
                return res;
            }

            template<class U> inline
            typename util::enable_if<util::is_same<U,gsFeSpace<Scalar> >::value, const gsMatrix<Scalar> & >::type
            eval_impl(const U & u, const index_t k) const
            {
                /*
                    We assume that the basis has the form v*e_i where e_i is the unit vector with 1 on index i and 0 elsewhere
                    This implies that hess(v) = [hess(v_1), hess(v_2), hess(v_3)] only has nonzero entries in column i. Hence,
                    hess(v) . normal = hess(v_i) * n_i (vector-scalar multiplication. The result is then of the form
                    [hess(v_1)*n_1 .., hess(v_2)*n_2 .., hess(v_3)*n_3 ..]. Here, the dots .. represent the active basis functions.
                */
                GISMO_ASSERT(_u.source().domainDim()==2,"Domain dimension should be 2x2. Other dimensions are not yet implemented");
                GISMO_ASSERT(_v.rows()==2,"v must be 2x2, but is = "<<_v.rows()<<" x "<<_v.cols());
                GISMO_ASSERT(_v.cols()==2,"v must be 2x2, but is = "<<_v.rows()<<" x "<<_v.cols());
                const index_t numAct = u.data().values[0].rows();   // number of actives of a basis function
                const index_t cardinality = u.cardinality();        // total number of actives (=3*numAct)
                const index_t nDers = _u.source().domainDim() * ( _u.source().domainDim() + 1 ) / 2;
                res.resize(rows()*cardinality, cols() );
                // tmp is ordered as
                // [d11 d22 d12] (in 2D)
                tmp =_u.data().values[2].reshapeCol(k, nDers, numAct );

                // vEv is ordered as
                // [v11, v12; v21 v22] (in 2D)
                vEv = _v.eval(k);

                // Cast to Voight notation
                vEv.resize(4,1);
                std::swap( vEv(3,0), vEv(1,0) );
                vEv.conservativeResize(3,1);

                // [dxx u, dyy u, dxy u] [1 0 0] [v11]
                //                       [0 1 0] [v22]
                //                       [0 0 2] [v12]
                gsMatrix<Scalar,3,3> ones;
                ones.setIdentity();
                ones(2,2) = 2.0;
                for (index_t i = 0; i!=numAct; i++)
                    res.row(i) = tmp.col(i).transpose() * ones * vEv;

                return res;
            }

            template<class U> inline
            typename util::enable_if<util::is_same<U,gsFeSolution<Scalar> >::value, const gsMatrix<Scalar> & >::type
            eval_impl(const U & u, const index_t k) const
            {
                /*
                    We assume that the basis has the form v*e_i where e_i is the unit vector with 1 on index i and 0 elsewhere
                    This implies that hess(v) = [hess(v_1), hess(v_2), hess(v_3)] only has nonzero entries in column i. Hence,
                    hess(v) . normal = hess(v_i) * n_i (vector-scalar multiplication. The result is then of the form
                    [hess(v_1)*n_1 .., hess(v_2)*n_2 .., hess(v_3)*n_3 ..]. Here, the dots .. represent the active basis functions.
                */
                hess_expr<gsFeSolution<Scalar>> sHess = hess_expr<gsFeSolution<Scalar>>(_u); // NOTE: This does not parse automatically!
                tmp = sHess.eval(k);
                vEv = _v.eval(k);

                // Cast to Voight notation
                vEv.resize(4,1);
                std::swap( vEv(3,0), vEv(1,0) );
                vEv.conservativeResize(3,1);

                // Cast to Voight notation
                tmp.resize(4,1);
                std::swap( tmp(3,0), tmp(1,0) );
                tmp.conservativeResize(3,1);


                // [dxx u, dyy u, dxy u] [1 0 0] [v11]
                //                       [0 1 0] [v22]
                //                       [0 0 2] [v12]
                gsMatrix<Scalar,3,3> ones;
                ones.setIdentity();
                ones(2,2) = 2.0;
                res = tmp.transpose() * ones * vEv;
                return res;
            }
        };

        template<class E> EIGEN_STRONG_INLINE
        deriv2_expr<E> deriv2(const E & u) { return deriv2_expr<E>(u); }

        template<class E1, class E2> EIGEN_STRONG_INLINE
        deriv2dot_expr<E1, E2> deriv2(const E1 & u, const E2 & v) { return deriv2dot_expr<E1, E2>(u,v); }

    }
}

using namespace gismo;

void setMapperForBiharmonic(gsBoundaryConditions<> & bc, gsMappedBasis<2,real_t> & bb2, gsDofMapper & mapper)
{
    mapper.setIdentity(bb2.nPatches(), bb2.size(), 1);

    gsMatrix<index_t> bnd;
    for (typename gsBoundaryConditions<real_t>::const_iterator
                 it = bc.begin("Dirichlet"); it != bc.end("Dirichlet"); ++it)
    {
        bnd = bb2.basis(it->ps.patch).boundary(it->ps.side());
        mapper.markBoundary(it->ps.patch, bnd, 0);
    }

    for (typename gsBoundaryConditions<real_t>::const_iterator
             it = bc.begin("Neumann"); it != bc.end("Neumann"); ++it)
    {
        bnd = bb2.basis(it->ps.patch).boundaryOffset(it->ps.side(),1);
        mapper.markBoundary(it->ps.patch, bnd, 0);
    }
    mapper.finalize();
}

void gsDirichletNeumannValuesL2Projection(gsMultiPatch<> & mp, gsMultiBasis<> & dbasis, gsBoundaryConditions<> & bc,
                                           gsMappedBasis<2,real_t> & bb2, const expr::gsFeSpace<real_t> & u)
{
    const gsDofMapper & mapper = u.mapper();

    gsMatrix<index_t> bnd = mapper.findFree(mapper.numPatches()-1);
    gsDofMapper mapperBdy;
    mapperBdy.setIdentity(bb2.nPatches(), bb2.size(), 1);  // bb2.nPatches() == 1
    mapperBdy.markBoundary(0, bnd, 0);
    mapperBdy.finalize();

    gsExprAssembler<real_t> A(1,1);
    A.setIntegrationElements(dbasis);

    auto G = A.getMap(mp);
    auto uu = A.getSpace(bb2);
    auto gg = pow(fform(G).det(),0.5);
    auto g_bdy = A.getBdrFunction(G); // Bug?!?

//    gsFunctionExpr<> ms("(cos(4*pi*x) - 1) * (cos(4*pi*y) - 1) + 0*z",3);
//    // Neumann
//    gsFunctionExpr<> sol1der("-4*pi*(cos(4*pi*y) - 1)*sin(4*pi*x)",
//                             "-4*pi*(cos(4*pi*x) - 1)*sin(4*pi*y)",
//                             "0", 3);
//    auto f_dir = A.getCoeff(ms, G);
//    auto f_neu = A.getCoeff(sol1der, G);
//
//        gsExprEvaluator<real_t> ev(A);
//        index_t N = 5;
//        gsMatrix<> points;
//        points.setZero(2,N);
//        gsVector<> pp;
//        pp.setLinSpaced(N,0,1);
//        points.row(0) = pp;
//
//
//
//        for (index_t i = 0; i < points.cols(); i++)
//        {
//            gsDebugVar(f_dir.parDim());
//            gsDebugVar(f_dir.targetDim());
//            gsDebugVar(f_dir.rows());
//            gsDebugVar(points.rows());
//            gsDebugVar(ev.eval(f_dir, points.col(i)));
//        }


    uu.setupMapper(mapperBdy);
    gsMatrix<real_t> & fixedDofs_A = const_cast<expr::gsFeSpace<real_t>&>(uu).fixedPart();
    fixedDofs_A.setZero( uu.mapper().boundarySize(), 1 );

    real_t lambda = 1e-5;

    // Mapped Spline
    A.initSystem();
    A.assembleBdr(bc.get("Dirichlet"), uu * uu.tr() * gg);
    A.assembleBdr(bc.get("Dirichlet"), uu * g_bdy * gg);
    A.assembleBdr(bc.get("Neumann"),
                  lambda * ((jac(G) * fform(G).inv() * igrad(uu).tr()).tr() * nv(G).normalized())
                               * ((jac(G) * fform(G).inv() * igrad(uu).tr()).tr() * nv(G).normalized()).tr() * gg);
    A.assembleBdr(bc.get("Neumann"),
                  lambda *  ((jac(G) * fform(G).inv() * igrad(uu).tr()).tr() * nv(G).normalized()) * (g_bdy.tr() * nv(G).normalized()) * gg);

    gsSparseSolver<real_t>::SimplicialLDLT solver;
    solver.compute( A.matrix() );
    gsMatrix<real_t> & fixedDofs = const_cast<expr::gsFeSpace<real_t>& >(u).fixedPart();
    fixedDofs = solver.solve(A.rhs());
}

void setMapperForBiharmonic(gsBoundaryConditions<> & bc, gsMultiBasis<> & dbasis, gsDofMapper & mapper)
{
    mapper.init(dbasis);

    for (gsBoxTopology::const_iiterator it = dbasis.topology().iBegin();
         it != dbasis.topology().iEnd(); ++it) // C^0 at the interface
    {
        dbasis.matchInterface(*it, mapper);
    }

    gsMatrix<index_t> bnd;
    for (typename gsBoundaryConditions<real_t>::const_iterator
                 it = bc.begin("Dirichlet"); it != bc.end("Dirichlet"); ++it)
    {
        bnd = dbasis.basis(it->ps.patch).boundary(it->ps.side());
        mapper.markBoundary(it->ps.patch, bnd, 0);
    }

    for (typename gsBoundaryConditions<real_t>::const_iterator
                 it = bc.begin("Neumann"); it != bc.end("Neumann"); ++it)
    {
        bnd = dbasis.basis(it->ps.patch).boundaryOffset(it->ps.side(),1);
        mapper.markBoundary(it->ps.patch, bnd, 0);
    }
    mapper.finalize();
}

void gsDirichletNeumannValuesL2Projection(gsMultiPatch<> & mp, gsMultiBasis<> & dbasis, gsBoundaryConditions<> & bc, const expr::gsFeSpace<real_t> & u)
{
    gsDofMapper mapper = u.mapper();
    gsDofMapper mapperBdy(dbasis, u.dim());
    for (gsBoxTopology::const_iiterator it = dbasis.topology().iBegin();
         it != dbasis.topology().iEnd(); ++it) // C^0 at the interface
    {
        dbasis.matchInterface(*it, mapperBdy);
    }
    for (size_t np = 0; np < mp.nPatches(); np++)
    {
        gsMatrix<index_t> bnd = mapper.findFree(np);
        mapperBdy.markBoundary(np, bnd, 0);
    }
    mapperBdy.finalize();

    gsExprAssembler<real_t> A(1,1);
    A.setIntegrationElements(dbasis);

    auto G = A.getMap(mp);
    auto uu = A.getSpace(dbasis);
    auto gg = pow(fform(G).det(),0.5);
    auto g_bdy = A.getBdrFunction(G); // Bug?!?

    uu.setupMapper(mapperBdy);
    gsMatrix<real_t> & fixedDofs_A = const_cast<expr::gsFeSpace<real_t>&>(uu).fixedPart();
    fixedDofs_A.setZero( uu.mapper().boundarySize(), 1 );

    real_t lambda = 1e-5;

    // Standard Spline
    A.initSystem();
    A.assembleBdr(bc.get("Dirichlet"), uu * uu.tr() * gg);
    A.assembleBdr(bc.get("Dirichlet"), uu * g_bdy * gg);
    A.assembleBdr(bc.get("Neumann"),
                  lambda * ((jac(G) * fform(G).inv() * igrad(uu).tr()).tr() * nv(G).normalized())
                  * ((jac(G) * fform(G).inv() * igrad(uu).tr()).tr() * nv(G).normalized()).tr() * gg);
    A.assembleBdr(bc.get("Neumann"),
                  lambda *  ((jac(G) * fform(G).inv() * igrad(uu).tr()).tr() * nv(G).normalized()) * (g_bdy.tr()  * nv(G).normalized()) * gg);

    gsSparseSolver<real_t>::SimplicialLDLT solver;
    solver.compute( A.matrix() );
    gsMatrix<real_t> & fixedDofs = const_cast<expr::gsFeSpace<real_t>& >(u).fixedPart();
    gsMatrix<real_t> fixedDofs_temp = solver.solve(A.rhs());

    // Reordering the dofs of the boundary
    fixedDofs.setZero(mapper.boundarySize(),1);
    index_t sz = 0;
    for (size_t np = 0; np < mp.nPatches(); np++)
    {
        gsMatrix<index_t> bnd = mapperBdy.findFree(np);
        bnd.array() += sz;
        for (index_t i = 0; i < bnd.rows(); i++)
        {
            index_t ii = mapperBdy.asVector()(bnd(i,0));
            fixedDofs(mapper.global_to_bindex(mapper.asVector()(bnd(i,0))),0) = fixedDofs_temp(ii,0);
        }
        sz += mapperBdy.patchSize(np,0);
    }
}


int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    bool mesh = false;

    index_t method = 0;

    index_t numRefine  = 5;
    index_t degree = 3;
    index_t smoothness = 2;

    real_t penalty_init = -1;

    bool info = false;
    bool last = false;
    bool second = false;
    bool residual = false;

    std::string output;
    std::string fn = "surfaces/surface_roof.xml";
    std::string geometry;

    gsCmdLine cmd("Example for solving the biharmonic problem (single patch only).");
    cmd.addInt( "m", "method", "The chosen method for the biharmonic problem", method );

    cmd.addInt("p", "degree","Set discrete polynomial degree", degree);
    cmd.addInt("s", "smoothness", "Set discrete regularity",  smoothness);
    cmd.addInt("r", "refinementLoop", "Number of refinement steps", numRefine);

    cmd.addString("f", "file", "Input geometry file (with .xml)", fn);
    cmd.addString( "g", "geometry", "Input geometry file",  geometry );

    // Flags related to Nitsche's method
    cmd.addReal( "y", "penalty", "Fixed Penalty value for Nitsche's method",  penalty_init);

    cmd.addSwitch("info", "Plot the information of the Approx C1 Basis", info);
    cmd.addSwitch("last", "Solve problem only on the last level of h-refinement", last);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("mesh", "Plot the mesh", mesh);
    cmd.addSwitch("residual", "Compute the error with residual", residual);

    cmd.addSwitch("second", "Compute second biharmonic problem with u = g1 and Delta u = g2 "
                            "(default first biharmonic problem: u = g1 and partial_n u = g2)", second);

    cmd.addString("o", "output", "Output in xml (for python)", output);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Read Argument inputs]
    gsMultiPatch<real_t> mp;
    std::string string_geo;
    if (geometry.empty())
        string_geo = fn;
    else
        string_geo = "surfaces/geometries/" + geometry + ".xml";

    gsInfo << "Filedata: " << string_geo << "\n";
    gsReadFile<>(string_geo, mp);
    mp.clearTopology();
    mp.computeTopology();
    //! [Read geometry]

//    gsFunctionExpr<> laplace ("(1/((1 + 4 * x^2 + 4 * y^2)^2) ) * (-8 * cos(8 * x) * ( (25 + 144 * x^4 + 164 * y^2 + 256 * y^4 + \n"
//                              "        8 * x^2 * (17 + 50 * y^2) ) * cos(1 - 6 * y) + 12 * y * (1 + 2 * x^2 + 2 * y^2) * sin(1 - 6 * y)) + \n"
//                              "  128 * x * sin(8 * x) * ( (1 + 2 * x^2 + 2 * y^2) * cos(1 - 6 * y) + 6 * y * (1 + 4 * x^2 + 4 * y^2) * sin(1 - 6 * y)))",3);
//
//    gsFunctionExpr<> ms("2 * cos(8 * x) * cos(1 - 6 * y)",3);
//
//    gsFunctionExpr<>sol1der ("-((16 * ( (1 + 4 * y^2) * cos(1 - 6 * y) * sin(8 * x) + 3 * x * y * cos(8 * x) * sin(1 - 6 * y) ) ) / (1 + 4 * x^2 + 4 * y^2))",
//                             "(64 * x * y * cos(1 - 6 * y) * sin(8 * x) + 12 * (1 + 4 * x^2) * cos(8 * x) * sin(1 - 6 * y) ) / (1 + 4 * x^2 + 4 * y^2)",
//                             "(8 * (4 * x * cos(1 - 6 * y) * sin(8 * x) - 3 * y * cos(8 * x) * sin(1 - 6 * y) ) ) / (1 + 4 * x^2 + 4 * y^2)",3);
//
//
////    gsFunctionExpr<>sol2der ("-( ( 2 * x * (1 + x^2) * y ) / (1 + x^2 + y^2)^2)",
////                             "-( ( 2 * x * y * (1 + y^2) ) / (1 + x^2 + y^2)^2)",
////                             "(2 * x * y ) / (1 + x^2 + y^2)^2",
////                             "(1 + y^2 + x^2 * (1 + 2 * y^2) ) / (1 + x^2 + y^2)^2",
////                             "(x * (1 + x^2 - y^2) ) / (1 + x^2 + y^2)^2",
////                             "(y * (1 - x^2 + y^2) ) / (1 + x^2 + y^2)^2", 3);
//
//
//    gsFunctionExpr<> f  ("(1 / ( (1 + 4 * x^2 + 4 * y^2)^5) ) * 8 * ( (1 + 4 * x^2 + 4 * y^2) * (5005 + 55732 * y^2 + \n"
//                              "      4 * (38416 * x^8 + 67765 * y^4 + 24 * y^6 * (5815 + 4374 * y^2) + 8 * x^6 * (8677 + 57232 * y^2) + \n"
//                              "         x^4 * (43637 + 397208 * y^2 + 905440 * y^4) + x^2 * (10637 + 126254 * y^2 + 467352 * y^4 + 590976 * y^6) ) ) * cos(9 * x) * cos(1 - 7 * y) - \n"
//                              "   288 * x * (36 + 367 * x^2 + 1521 * x^4 + 2880 * x^6 + 2352 * x^8 + 4 * (67 + 513 * x^2 + 1368 * x^4 + 1560 * x^6) * y^2 + \n"
//                              "      9 * (59 + 256 * (x^2 + 2 * x^4) ) * y^4 - 96 * (3 + x^2) * y^6 - 816 * y^8) * cos(1 - 7 * y) * sin(9 * x) + \n"
//                              "   56 * y * (4 * (36 - 5424 * x^8 + 415 * y^2 + 2001 * y^4 + 4416 * y^6 + 3888 * y^8 - 288 * x^6 * (17 + 43 * y^2) + \n"
//                              "         4 * x^2 * (31 + 273 * y^2 + 984 * y^4 + 1560 * y^6) - 3 * x^4 * (303 + 256 * y^2 * (7 + 6 * y^2) ) ) * cos(9 * x) - \n"
//                              "      9 * x * (1 + 4 * x^2 + 4 * y^2) * (103 + 1568 * x^6 + 786 * y^2 + 2448 * y^4 + 2592 * y^6 + 16 * x^4 * (121 + 358 * y^2) + \n"
//                              "         x^2 * (722 + 4384 * y^2 + 6752 * y^4) ) * sin(9 * x) ) * sin(1 - 7 * y) )",3);


    //gsFunctionExpr<>f("256*pi*pi*pi*pi*(4*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))",3);
    //gsFunctionExpr<>f("256*pi*pi*pi*pi*(4*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))",2);
    gsFunctionExpr<>f("5",3);
    gsInfo << "Source function: " << f << "\n";

    //gsFunctionExpr<> ms("(cos(4*pi*x) - 1) * (cos(4*pi*y) - 1)",3);
    //gsFunctionExpr<> ms("(cos(4*pi*x) - 1) * (cos(4*pi*y) - 1)",2);
    gsFunctionExpr<> ms("0",3);
    gsInfo << "Exact function: " << ms << "\n";

    //! [Refinement]
    gsMultiBasis<real_t> basis(mp, false);//true: poly-splines (not NURBS)

    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    basis.setDegree(degree); // preserve smoothness

    if (method == MethodFlags::DPATCH || method == MethodFlags::ALMOSTC1 || method == MethodFlags::SURFASG1)
        mp.degreeElevate(degree-mp.patch(0).degree(0));

    // h-refine each basis
    if (last)
    {
        for (index_t r =0; r < numRefine; ++r)
            basis.uniformRefine(1, degree-smoothness);

        numRefine = 0;
    }

    // Assume that the condition holds for each patch TODO
    // Refine once
    if (basis.basis(0).numElements() < 4)
    {
        basis.uniformRefine(1, degree-smoothness);
        if (method == MethodFlags::DPATCH || method == MethodFlags::ALMOSTC1 || method == MethodFlags::SURFASG1)
            mp.uniformRefine(1, degree-smoothness);
    }
    //! [Refinement]

    gsInfo << "OptionList: " << cmd << "\n";

    //! [Boundary condition]
    // Laplace
    //gsFunctionExpr<> laplace ("-16*pi*pi*(2*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))",3);
    //gsFunctionExpr<> laplace ("-16*pi*pi*(2*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))",2);
    gsFunctionExpr<> laplace ("0",3);
    // Neumann
//    gsFunctionExpr<> sol1der("-4*pi*(cos(4*pi*y) - 1)*sin(4*pi*x)",
//                     "-4*pi*(cos(4*pi*x) - 1)*sin(4*pi*y)",
//                     "0", 3);
//    gsFunctionExpr<> sol1der("-4*pi*(cos(4*pi*y) - 1)*sin(4*pi*x)",
//                     "-4*pi*(cos(4*pi*x) - 1)*sin(4*pi*y)",2);
    gsFunctionExpr<> sol1der("0",
                     "0",
                     "0", 3);

    gsBoundaryConditions<> bc;
    for (gsMultiPatch<>::const_biterator bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
    {
        bc.addCondition(*bit, condition_type::dirichlet, ms);
        if (second)
            bc.addCondition(*bit, condition_type::laplace, laplace);
        else
            bc.addCondition(*bit, condition_type::neumann, sol1der);
    }
    bc.setGeoMap(mp);
    gsInfo << "Boundary conditions:\n" << bc << "\n";
    //! [Boundary condition]

    //! [Problem setup]
    gsExprAssembler<real_t> A(1,1);
    //gsInfo<<"Active options:\n"<< A.options() <<"\n";

    // Elements used for numerical integration
    A.setIntegrationElements(basis);
    gsExprEvaluator<real_t> ev(A);

    // Set the geometry map
    auto G = A.getMap(mp);

    // Set the source term
    auto ff = A.getCoeff(f, G); // Laplace example

    // Set the discretization space
    gsMappedBasis<2,real_t> bb2;
    auto u = method == MethodFlags::NITSCHE ? A.getSpace(basis) : A.getSpace(bb2);

    // Solution vector and solution variable
    gsMatrix<real_t> solVector, solVec2;
    auto u_sol = A.getSolution(u, solVector);
    gsMappedSpline<2, real_t> ms_coarse;
    gsMultiPatch<real_t> sol_nitsche_coarse;

    // Recover manufactured solution
    auto u_ex = ev.getVariable(ms, G);
    //! [Problem setup]

#ifdef _OPENMP
    gsInfo << "Available threads: "<< omp_get_max_threads() <<"\n";
#endif

    //! [Solver loop]
    gsSparseSolver<real_t>::SimplicialLDLT solver;

    gsVector<real_t> l2err(numRefine+1), h1err(numRefine+1), h2err(numRefine+1), IFaceErr(numRefine+1),
            dofs(numRefine+1), meshsize(numRefine+1);
    gsMatrix<real_t> penalty(numRefine+1, mp.nInterfaces());
    gsInfo<< "(dot1=assembled, dot2=solved, dot3=got_error)\n"
             "\nDoFs: ";
    double setup_time(0), ma_time(0), slv_time(0), err_time(0);
    gsStopwatch timer;
    for (index_t r=0; r<=numRefine; ++r)
    {

        if (method == MethodFlags::SPLINE)
        {
            if (mp.nPatches() != 1)
            {
                gsInfo << "The geometry has more than one patch. Run the code with a single patch!\n";
                return EXIT_FAILURE;
            }

            // Refine uniform once
            basis.uniformRefine(1,degree -smoothness);
            meshsize[r] = basis.basis(0).getMaxCellLength();

            gsSparseMatrix<real_t> global2local(basis.size(), basis.size());
            global2local.setIdentity();
            bb2.init(basis,global2local);
            gsInfo << "Spline basis created \n";
        }
        else if (method == MethodFlags::APPROXC1)
        {
            basis.uniformRefine(1,degree -smoothness);
            meshsize[r] = basis.basis(0).getMinCellLength();

            // The approx. C1 space
            gsApproxC1Spline<2,real_t> approxC1(mp,basis);
            approxC1.options().setSwitch("info",info);
            approxC1.options().setSwitch("plot",plot);
            approxC1.options().setSwitch("interpolation",true);
            approxC1.options().setSwitch("second",second);
            approxC1.options().setInt("gluingDataDegree",-1);
            approxC1.options().setInt("gluingDataSmoothness",-1);
            approxC1.update(bb2);
        }
        else if (method == MethodFlags::DPATCH)
        {
            mp.uniformRefine(1,degree-smoothness);
            basis.uniformRefine(1,degree-smoothness);

            if (gsHTensorBasis<2,real_t> * test = dynamic_cast<gsHTensorBasis<2,real_t>*>(&basis.basis(0)))
                meshsize[r] = test->tensorLevel(0).getMinCellLength();
            else if (gsTensorBasis<2,real_t> * test = dynamic_cast<gsTensorBasis<2,real_t>*>(&basis.basis(0)))
                meshsize[r] = test->getMinCellLength();

            gsSparseMatrix<real_t> global2local;
            gsDPatch<2,real_t> dpatch(mp);
            dpatch.matrix_into(global2local);
            global2local = global2local.transpose();
            mp = dpatch.exportToPatches();
            basis = dpatch.localBasis();
            bb2.init(basis,global2local);
            gsInfo << "DPATCH basis created \n";
        }
        else if (method == MethodFlags::ALMOSTC1)
        {
            mp.uniformRefine(1,degree-smoothness);
            basis.uniformRefine(1,degree-smoothness);

            if (gsHTensorBasis<2,real_t> * test = dynamic_cast<gsHTensorBasis<2,real_t>*>(&basis.basis(0)))
                meshsize[r] = test->tensorLevel(0).getMinCellLength();
            else if (gsTensorBasis<2,real_t> * test = dynamic_cast<gsTensorBasis<2,real_t>*>(&basis.basis(0)))
                meshsize[r] = test->getMinCellLength();

            gsSparseMatrix<real_t> global2local;
            gsAlmostC1<2,real_t> almostC1(mp);
            almostC1.matrix_into(global2local);
            global2local = global2local.transpose();
            mp = almostC1.exportToPatches();
            basis = almostC1.localBasis();
            bb2.init(basis,global2local);
            gsInfo << "ALMOSTC1 basis created \n";
        }
        else if (method == MethodFlags::NITSCHE)
        {
            basis.uniformRefine(1,degree-smoothness);
            meshsize[r] = basis.basis(0).getMinCellLength();
        }
        else if (method == MethodFlags::SURFASG1) // Andrea
        {
            mp.uniformRefine(1,degree-smoothness);
            basis.uniformRefine(1,degree-smoothness);

            meshsize[r] = basis.basis(0).getMinCellLength();

            gsC1SurfSpline<2,real_t> smoothC1(mp,basis);
            smoothC1.init();
            smoothC1.compute();

            gsSparseMatrix<real_t> global2local;
            global2local = smoothC1.getSystem();
            global2local = global2local.transpose();
            gsMultiBasis<> basis_temp;
            smoothC1.getMultiBasis(basis_temp);
            bb2.init(basis_temp,global2local);
        }

        // Setup the mapper
        if (method == MethodFlags::APPROXC1 || method == MethodFlags::DPATCH || method == MethodFlags::ALMOSTC1
                                                                                || method == MethodFlags::SURFASG1) // MappedBasis
        {
            gsDofMapper map;
            setMapperForBiharmonic(bc, bb2,map);

            // Setup the system
            u.setupMapper(map);
            gsDirichletNeumannValuesL2Projection(mp, basis, bc, bb2, u);
        }
        else if (method == MethodFlags::NITSCHE) // Nitsche
        {
            gsDofMapper map;
            setMapperForBiharmonic(bc, basis,map);

            // Setup the system
            u.setupMapper(map);
            gsDirichletNeumannValuesL2Projection(mp, basis, bc, u);
        }

        // Initialize the system
        A.initSystem();
        setup_time += timer.stop();

        gsInfo<< A.numDofs() <<std::flush;

        timer.restart();
        // Compute the system matrix and right-hand side
        auto gg = pow(fform(G).det(),0.5);
        auto G0 = 1.0/gg * fform(G);
        auto G0inv = G0.inv();
        A.assemble(( deriv2(u,G0inv) * (deriv2(u,G0inv)).tr()) * 1.0/gg,
                   u * ff * gg);
        //gsInfo << "Finished \n";
        // Enforce Laplace conditions to right-hand side
        // auto g_L = A.getBdrFunction(G); // Set the laplace bdy value  // Bug, doesnt work
        auto g_L = A.getCoeff(laplace, G);
        A.assembleBdr(bc.get("Laplace"), ((jac(G) * fform(G).inv() * igrad(u).tr()).tr() * nv(G).normalized()) * g_L.tr() * gg );

        if (method == MethodFlags::NITSCHE)
        {
            index_t i = 0;
            for ( typename gsMultiPatch<real_t>::const_iiterator it = mp.iBegin(); it != mp.iEnd(); ++it, ++i)
            {
                real_t stab     = 4 * ( basis.maxCwiseDegree() + basis.dim() ) * ( basis.maxCwiseDegree() + 1 );
                real_t m_h      = basis.basis(0).getMinCellLength(); //*dbasis.basis(0).getMinCellLength();
                real_t mu       = 2 * stab / m_h;
                real_t alpha = 1;

                if (penalty_init != -1.0)
                    mu = penalty_init / m_h;

                penalty(r,i) = mu;

                auto ggL = pow(fform(G.left()).det(),0.5);
                auto G0L = 1.0/ggL * fform(G.left());
                auto G0invL = G0L.inv();

                auto ggR = pow(fform(G.right()).det(),0.5);
                auto G0R = 1.0/ggR * fform(G.right());
                auto G0invR = G0R.inv();

                std::vector<boundaryInterface> iFace;
                iFace.push_back(*it);
                A.assembleIfc(iFace,
                        //B11
                              -alpha * 0.5 * ((jac(G.left()) * fform(G.left()).inv() * igrad(u.left()).tr()).tr() * nv(G.left()).normalized()) *
                              (deriv2(u.left(),G0invL)).tr() * gg,
                              -alpha * 0.5 *
                                      (((jac(G.left()) * fform(G.left()).inv() * igrad(u.left()).tr()).tr() * nv(G.left()).normalized())
                              * (deriv2(u.left(),G0invL)).tr()).tr() *
                              gg,
                        //B12
                              -alpha * 0.5 * ((jac(G.left()) * fform(G.left()).inv() * igrad(u.left()).tr()).tr() * nv(G.left()).normalized()) *
                              (deriv2(u.right(),G0invR)).tr() * gg,
                              -alpha * 0.5 * (((jac(G.left()) * fform(G.left()).inv() * igrad(u.left()).tr()).tr() * nv(G.left()).normalized()) *
                                              (deriv2(u.right(),G0invR)).tr()).tr() * gg,
                        //B21
                              alpha * 0.5 * ((jac(G.right()) * fform(G.right()).inv() * igrad(u.right()).tr()).tr() * nv(G.left()).normalized()) *
                              (deriv2(u.left(),G0invL)).tr() * gg,
                              alpha * 0.5 * (((jac(G.right()) * fform(G.right()).inv() * igrad(u.right()).tr()).tr() * nv(G.left()).normalized()) *
                                             (deriv2(u.left(),G0invL)).tr()).tr() * gg,
                        //B22
                              alpha * 0.5 * ((jac(G.right()) * fform(G.right()).inv() * igrad(u.right()).tr()).tr() * nv(G.left()).normalized()) *
                              (deriv2(u.right(),G0invR)).tr() * gg,
                              alpha * 0.5 * (((jac(G.right()) * fform(G.right()).inv() * igrad(u.right()).tr()).tr() * nv(G.left()).normalized()) *
                                             (deriv2(u.right(),G0invR)).tr()).tr() * gg,

                        // E11
                              mu * ((jac(G.left()) * fform(G.left()).inv() * igrad(u.left()).tr()).tr() * nv(G.left()).normalized()) *
                              (((jac(G.left()) * fform(G.left()).inv() * igrad(u.left()).tr()).tr() * nv(G.left()).normalized())).tr() * gg,
                        //-E12
                              -mu * (((jac(G.left()) * fform(G.left()).inv() * igrad(u.left()).tr()).tr() * nv(G.left()).normalized())) *
                              (((jac(G.right()) * fform(G.right()).inv() * igrad(u.right()).tr()).tr() * nv(G.left()).normalized())).tr() * gg,
                        //-E21
                              -mu * (((jac(G.right()) * fform(G.right()).inv() * igrad(u.right()).tr()).tr() * nv(G.left()).normalized())) *
                              (((jac(G.left()) * fform(G.left()).inv() * igrad(u.left()).tr()).tr() * nv(G.left()).normalized())).tr() * gg,
                        // E22
                              mu * ((jac(G.right()) * fform(G.right()).inv() * igrad(u.right()).tr()).tr() * nv(G.left()).normalized()) *
                              (((jac(G.right()) * fform(G.right()).inv() * igrad(u.right()).tr()).tr() * nv(G.left()).normalized())).tr() * gg
                );
            }
        }

        dofs[r] = A.numDofs();
        ma_time += timer.stop();
        gsInfo << "." << std::flush;// Assemblying done

        timer.restart();
        solver.compute( A.matrix() );
        solVector = solver.solve(A.rhs());

        slv_time += timer.stop();
        gsInfo << "." << std::flush; // Linear solving done

        timer.restart();

//        index_t N = 5;
//        gsMatrix<> points;
//        points.setZero(2,N);
//        gsVector<> pp;
//        pp.setLinSpaced(N,0,1);
//        points.row(0) = pp;
//
//        for (index_t i = 0; i < points.cols(); i++)
//        {
//            gsDebugVar(ev.eval(u_ex, points.col(i)));
//            gsDebugVar( ev.eval(u_sol, points.col(i) ));
//        }

        //linferr[r] = ev.max( f-s ) / ev.max(f);
        if (residual)
        {
            if (method != MethodFlags::NITSCHE) {
                if (r != 0) {
                    auto u_coarse = A.getCoeff(ms_coarse);
                    l2err[r] = math::sqrt(
                            ev.integral((u_coarse - u_sol).sqNorm() * gg)); // / ev.integral(ff.sqNorm()*meas(G)) );
                    h1err[r] = l2err[r] +
                               math::sqrt(ev.integral((igrad(u_coarse) - igrad(u_sol)).sqNorm() *
                                                      gg)); // /ev.integral( igrad(f).sqNorm()*meas(G) ) );

                    h2err[r] = h1err[r] +
                               math::sqrt(ev.integral((ihess(u_coarse) - ihess(u_sol)).sqNorm() *
                                                      gg)); // /ev.integral( ihess(f).sqNorm()*meas(G) )
                } else {
                    l2err[r] = 0;
                    h1err[r] = 0;
                    h2err[r] = 0;

                }
                gsMatrix<real_t> solFull_coarse;
                u_sol.extractFull(solFull_coarse);
                ms_coarse.init(bb2, solFull_coarse);
            }
            else if (method == MethodFlags::NITSCHE)
            {
                if (r != 0) {
                    auto u_coarse = A.getCoeff(sol_nitsche_coarse);
                    l2err[r] = math::sqrt(
                            ev.integral((u_coarse - u_sol).sqNorm() * gg)); // / ev.integral(ff.sqNorm()*meas(G)) );
                    h1err[r] = l2err[r] +
                               math::sqrt(ev.integral((igrad(u_coarse) - igrad(u_sol)).sqNorm() *
                                                      gg)); // /ev.integral( igrad(f).sqNorm()*meas(G) ) );

                    h2err[r] = h1err[r] +
                               math::sqrt(ev.integral((ihess(u_coarse) - ihess(u_sol)).sqNorm() *
                                                      gg)); // /ev.integral( ihess(f).sqNorm()*meas(G) )
                } else {
                    l2err[r] = 0;
                    h1err[r] = 0;
                    h2err[r] = 0;

                }
                u_sol.extract(sol_nitsche_coarse);
            }
        }
        else
        {
            l2err[r]= math::sqrt( ev.integral( (u_ex - u_sol).sqNorm() * gg ) ); // / ev.integral(ff.sqNorm()*meas(G)) );
            h1err[r]= l2err[r] +
                      math::sqrt(ev.integral( ( igrad(u_ex) - (jac(G) * fform(G).inv() * igrad(u_sol).tr()).tr() ).sqNorm() * gg )); // /ev.integral( igrad(ff).sqNorm()*meas(G) ) );

            //h2err[r]= h1err[r] +
            //          math::sqrt(ev.integral( ( ihess(u_ex) - 1.0/gg * (deriv2(u_sol,G0inv)).tr() ).sqNorm() * meas(G) )); // /ev.integral( ihess(ff).sqNorm()*meas(G) )
        }

        // Jump error
        if (residual)
        {
            if (method != MethodFlags::NITSCHE)
            {
                gsMatrix<real_t> solFull;
                u_sol.extractFull(solFull);
                gsMappedSpline<2, real_t> mappedSpline(bb2, solFull);
                auto ms_sol = A.getCoeff(mappedSpline);

                IFaceErr[r] = math::sqrt(ev.integralInterface(
                        (((jac(G.left()) * fform(G.left()).inv() * igrad(ms_sol.left()).tr()).tr() -
                          (jac(G.right()) * fform(G.right()).inv() * igrad(ms_sol.right()).tr()).tr()) *
                         nv(G).normalized()).sqNorm() * gg));
            }
            else if (method == MethodFlags::NITSCHE)
            {
                auto ms_sol = A.getCoeff(sol_nitsche_coarse);
                IFaceErr[r] = math::sqrt(ev.integralInterface(
                        (((jac(G.left()) * fform(G.left()).inv() * igrad(ms_sol.left()).tr()).tr() -
                          (jac(G.right()) * fform(G.right()).inv() * igrad(ms_sol.right()).tr()).tr()) *
                         nv(G).normalized()).sqNorm() * gg));

                // This doesn't work yet. Bug?
                //IFaceErr[r] = math::sqrt(ev.integralInterface((( igrad(u_sol.left(), G.left()) -
                //                                                 igrad(u_sol.right(), G.right())) *
                //                                               nv(G).normalized()).sqNorm() * meas(G)));
            }
        }


        err_time += timer.stop();
        gsInfo << ". " << std::flush; // Error computations done
    } //for loop
    //! [Solver loop]

    timer.stop();
    gsInfo << "\n\nTotal time: "<< setup_time+ma_time+slv_time+err_time <<"\n";
    gsInfo << "     Setup: "<< setup_time <<"\n";
    gsInfo << "  Assembly: "<< ma_time    <<"\n";
    gsInfo << "   Solving: "<< slv_time   <<"\n";
    gsInfo << "     Norms: "<< err_time   <<"\n";
    gsInfo << " Mesh-size: "<< meshsize.transpose() << "\n";
    gsInfo << "      Dofs: "<<dofs.transpose() << "\n";

    //! [Error and convergence rates]
    gsInfo << "\nL2 error: "<<std::scientific<<std::setprecision(3)<<l2err.transpose()<<"\n";
    gsInfo << "H1 error: "<<std::scientific<<h1err.transpose()<<"\n";
    gsInfo << "H2 error: "<<std::scientific<<h2err.transpose()<<"\n";
    gsInfo<< "Deriv Interface error: "<<std::scientific<<IFaceErr.transpose()<<"\n";
    if (method == MethodFlags::NITSCHE)
        gsInfo<< "\nStabilization: " << penalty.transpose() << "\n";

    if (numRefine>0)
    {
        gsInfo<< "\nEoC (L2): " << std::fixed<<std::setprecision(2)
              <<  ( l2err.head(numRefine).array()  /
                    l2err.tail(numRefine).array() ).log().transpose() / std::log(2.0)
              <<"\n";

        gsInfo<<   "EoC (H1): "<< std::fixed<<std::setprecision(2)
              <<( h1err.head(numRefine).array() /
                  h1err.tail(numRefine).array() ).log().transpose() / std::log(2.0) <<"\n";

        gsInfo<<   "EoC (H2): "<< std::fixed<<std::setprecision(2)
              <<( h2err.head(numRefine).array() /
                  h2err.tail(numRefine).array() ).log().transpose() / std::log(2.0) <<"\n";

        gsInfo<<   "EoC (Iface): "<< std::fixed<<std::setprecision(2)
              <<( IFaceErr.head(numRefine).array() /
                  IFaceErr.tail(numRefine).array() ).log().transpose() / std::log(2.0) <<"\n";
    }
    //! [Error and convergence rates]

    //! [Export visualization in ParaView]
    if (plot)
    {
        // Write approximate and exact solution to paraview files
        gsInfo << "Plotting in Paraview...\n";
        ev.options().setSwitch("plot.elements", mesh);
        ev.writeParaview( u_sol   , G, "solution");
        //ev.writeParaview( u_ex - u_sol   , G, "solution_pointwise");
        gsInfo << "Saved with solution.pvd \n";
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";
    //! [Export visualization in ParaView]

    //! [Export data to xml]
    if (!output.empty())
    {
        index_t cols = method == MethodFlags::NITSCHE ? 7+penalty.cols() : 7;
        gsMatrix<real_t> error_collection(l2err.rows(), cols);
        error_collection.col(0) = meshsize;
        error_collection.col(1) = dofs;
        error_collection.col(2) = l2err;
        error_collection.col(3) = h1err;
        error_collection.col(4) = h2err;
        error_collection.col(5) = IFaceErr;
        //error_collection.col(6) = cond_num;
        if (method == MethodFlags::NITSCHE)
            error_collection.block(0,7,penalty.rows(),penalty.cols()) = penalty;

        gsFileData<real_t> xml_out;
        xml_out << error_collection;
        xml_out.addString("Meshsize, dofs, l2err, h1err, h2err, iFaceErr, (penalty)");
        // Add solution
        // [...]
        xml_out.save(output);
        gsInfo << "XML saved to " + output << "\n";
    }
    //! [Export data to xml]

    return  EXIT_SUCCESS;
}
