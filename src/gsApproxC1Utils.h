/** @file gsApproxC1Utils.h

    @brief Provides declaration of Basis abstract interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/
#pragma once


#include <gsCore/gsFunctionSet.h>

namespace gismo
{

void createGluingDataSpace(const gsGeometry<real_t> & patch, const gsBasis<real_t> & basis, index_t dir,
                           gsBSplineBasis<real_t> & result, index_t p_tilde, index_t r_tilde);

void createPlusSpace(const gsGeometry<real_t> & patch, gsBasis<real_t> & basis, index_t dir, gsBSplineBasis<real_t> & res_plus);

void createMinusSpace(const gsGeometry<real_t> & patch, gsBasis<real_t> & basis, index_t dir, gsBSplineBasis<real_t> & res_minus);

void createEdgeSpace(const gsGeometry<real_t> & patch, gsBasis<real_t> & basis, index_t dir, gsBSplineBasis<real_t> & basis_plus,
                     gsBSplineBasis<real_t> & basis_minus, gsBSplineBasis<real_t> & basis_gD, gsTensorBSplineBasis<2, real_t> & result);

void createEdgeSpace(const gsGeometry<real_t> & patch, gsBasis<real_t> & basis, index_t dir, gsBSplineBasis<real_t> & basis_plus,
                     gsBSplineBasis<real_t> & basis_minus, gsTensorBSplineBasis<2, real_t> & result);

void createVertexSpace(const gsGeometry<real_t> & patch, gsBasis<real_t> & basis, bool isInterface_1, bool isInterface_2,
                       gsTensorBSplineBasis<2, real_t> & result, index_t p_tilde, index_t r_tilde);

// Input is parametric coordinates of 1-D \a mp
template <class T>
class gsAlpha : public gismo::gsFunction<T>
{

protected:
    gsGeometry<T> & _geo;
    mutable gsMapData<T> _tmp;
    index_t m_uv;


public:
    /// Shared pointer for gsAlpha
    typedef memory::shared_ptr< gsAlpha > Ptr;

    /// Unique pointer for gsAlpha
    typedef memory::unique_ptr< gsAlpha > uPtr;

    gsAlpha(gsGeometry<T> & geo, index_t uv) :
            _geo(geo), m_uv(uv), _alpha_piece(nullptr)
    {
        _tmp.flags = NEED_JACOBIAN;
    }

    ~gsAlpha() { delete _alpha_piece; }

    GISMO_CLONE_FUNCTION(gsAlpha)

    short_t domainDim() const {return 1;}

    short_t targetDim() const {return 1;}

    mutable gsAlpha<T> * _alpha_piece; // why do we need this?

    const gsFunction<T> & piece(const index_t k) const
    {
        //delete _alpha_piece;
        _alpha_piece = new gsAlpha(*this);
        return *_alpha_piece;
    }

    // Input is parametric coordinates of 1-D \a mp
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize( targetDim() , u.cols() );

        if(_geo.parDim() == _geo.targetDim()) // Planar
        {
            gsMatrix<T> uv, ev, D0;
            uv.setZero(2, u.cols());
            uv.row(m_uv) = u; // u

            T gamma = 1.0;

            for (index_t i = 0; i < uv.cols(); i++) {
                _geo.jacobian_into(uv.col(i), ev);
                uv(0, i) = gamma * ev.determinant();
            }
            result = uv.row(0);
        }
        else if(_geo.parDim() + 1 == _geo.targetDim()) // Surface
        {
            gsMatrix<T> uv, ev, D0, N, Duv;
            uv.setZero(2, u.cols());
            uv.row(m_uv) = u; // u

            for (index_t i = 0; i < uv.cols(); i++) {
                N.setZero(3,1);
                _geo.deriv_into(uv.col(i), ev);
                N.row(0) = ev.row(2) * ev.row(5) - ev.row(4) * ev.row(3);
                N.row(1) = ev.row(4) * ev.row(1) - ev.row(0) * ev.row(5);
                N.row(2) = ev.row(0) * ev.row(3) - ev.row(2) * ev.row(1);

                D0.setZero(3,3);
                Duv.setZero(3,1);
                if (m_uv == 0) // u
                {
                    Duv.row(0) = ev.row(0);
                    Duv.row(1) = ev.row(2);
                    Duv.row(2) = ev.row(4);
                    D0.col(0) = Duv;
                    D0.col(2) = N;
                    D0(0,1) = ev(1,0);
                    D0(1,1) = ev(3,0);
                    D0(2,1) = ev(5,0);
                }
                else if (m_uv == 1) // v
                {
                    Duv.row(0) = ev.row(1);
                    Duv.row(1) = ev.row(3);
                    Duv.row(2) = ev.row(5);
                    D0.col(1) = Duv;
                    D0.col(2) = N;
                    D0(0,0) = ev(0,0);
                    D0(1,0) = ev(2,0);
                    D0(2,0) = ev(4,0);
                }
                uv(0, i) = 1.0 / Duv.norm() * gsEigen::numext::sqrt(D0.determinant());
            }
            result = uv.row(0);
            //gsDebugVar(result);
            //gsDebugVar(m_uv);
        }
    }
};


template <class T>
class gsBeta : public gismo::gsFunction<T>
{

protected:
    gsGeometry<T> & _geo;
    mutable gsMapData<T> _tmp;
    index_t m_uv;


public:
    /// Shared pointer for gsBeta
    typedef memory::shared_ptr< gsBeta > Ptr;

    /// Unique pointer for gsBeta
    typedef memory::unique_ptr< gsBeta > uPtr;

    gsBeta(gsGeometry<T> & geo, index_t uv) :
            _geo(geo), m_uv(uv), _beta_piece(nullptr)
    {
        _tmp.flags = NEED_JACOBIAN;
    }

    ~gsBeta() { delete _beta_piece; }

    GISMO_CLONE_FUNCTION(gsBeta)

    short_t domainDim() const {return 1;}

    short_t targetDim() const {return 1;}

    mutable gsBeta<T> * _beta_piece; // why do we need this?

    const gsFunction<T> & piece(const index_t k) const
    {
        //delete _beta_piece;
        _beta_piece = new gsBeta(*this);
        return *_beta_piece;
    }

    // Input is parametric coordinates of 1-D \a mp
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize( targetDim() , u.cols() );

        if(_geo.parDim() == _geo.targetDim()) // Planar
        {
            gsMatrix<T> uv, ev, D0;

            uv.setZero(2,u.cols());
            uv.row(m_uv) = u; // u

            T gamma = 1.0;

            for(index_t i = 0; i < uv.cols(); i++)
            {
                _geo.jacobian_into(uv.col(i),ev);
                D0  = ev.col(m_uv);
                T D1 = T(1.0)/ D0.norm();
                uv(0,i) = - gamma * D1 * D1 * ev.col(1).transpose() * ev.col(0);
            }
            result = uv.row(0);
        }
        else if(_geo.parDim() + 1 == _geo.targetDim()) // Surface
        {
            gsMatrix<T> uv, ev, D0;

            uv.setZero(2,u.cols());
            uv.row(m_uv) = u; // u

            T gamma = 1.0;

            for(index_t i = 0; i < uv.cols(); i++)
            {
                _geo.jacobian_into(uv.col(i),ev);
                D0  = ev.col(m_uv);
                T D1 = T(1.0)/ D0.norm();
                uv(0,i) = - gamma * D1 * D1 * ev.col(1).transpose() * ev.col(0);
            }
            result = uv.row(0);
            //gsDebugVar(result);
            //gsDebugVar(m_uv);
        }
    }

};

template <class T>
class gsTraceBasis : public gismo::gsFunction<T>
{

protected:
    gsGeometry<T> & _geo;

    gsBSpline<T>  m_basis_beta;
    gsBSplineBasis<T>  m_basis_plus;

    gsBasis<T> &   m_basis;

    mutable gsMapData<T> _tmp;

    bool m_isboundary;
    const index_t m_bfID, m_uv;


public:
    /// Shared pointer for gsTraceBasis
    typedef memory::shared_ptr< gsTraceBasis > Ptr;

    /// Unique pointer for gsTraceBasis
    typedef memory::unique_ptr< gsTraceBasis > uPtr;

    gsTraceBasis(gsGeometry<T> & geo,
                 gsBSpline<T> basis_beta,
                 gsBSplineBasis<T> basis_plus,
                 gsBasis<T> & basis,
                 bool isboundary,
                 const index_t bfID,
                 const index_t uv) :
            _geo(geo), m_basis_beta(basis_beta), m_basis_plus(basis_plus), m_basis(basis),
            m_isboundary(isboundary), m_bfID(bfID), m_uv(uv), _traceBasis_piece(nullptr)
    {
        //_tmp.flags = NEED_JACOBIAN;

        //createPlusSpace(geo, basis, m_uv, m_basis_plus); // Not efficient
    }

    ~gsTraceBasis() { delete _traceBasis_piece; }

    GISMO_CLONE_FUNCTION(gsTraceBasis)

    short_t domainDim() const {return 2;}

    short_t targetDim() const {return 1;}

    mutable gsTraceBasis<T> * _traceBasis_piece; // why do we need this?

    const gsFunction<T> & piece(const index_t k) const
    {
        //delete _traceBasis_piece;
        _traceBasis_piece = new gsTraceBasis(*this);
        return *_traceBasis_piece;
    }

    // Input is parametric coordinates of 2-D \a mp
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize( targetDim() , u.cols() );

        // tau/p
        gsBSplineBasis<T> bsp_temp = dynamic_cast<gsBSplineBasis<T> & >(m_basis.component(1-m_uv));

        T p = bsp_temp.degree();
        T tau_1 = bsp_temp.knots().at(p + 1); // p + 2

        gsMatrix<T> beta, N_0, N_1, N_i_plus, der_N_i_plus;

        if (!m_isboundary)
            m_basis_beta.eval_into(u.row(m_uv),beta); // 1-dir == PatchID
        else
            beta.setZero(1, u.cols());

        m_basis.component(1-m_uv).evalSingle_into(0,u.row(1-m_uv),N_0); // u
        m_basis.component(1-m_uv).evalSingle_into(1,u.row(1-m_uv),N_1); // u

        m_basis_plus.evalSingle_into(m_bfID,u.row(m_uv),N_i_plus); // v
        m_basis_plus.derivSingle_into(m_bfID,u.row(m_uv),der_N_i_plus);

        gsMatrix<T> temp = beta.cwiseProduct(der_N_i_plus);
        result = N_i_plus.cwiseProduct(N_0 + N_1) - temp.cwiseProduct(N_1) * tau_1 / p;
    }

};


template <class T>
class gsNormalDerivBasis : public gismo::gsFunction<T>
{

protected:
    gsGeometry<T> & _geo;

    gsBSpline<T> m_basis_alpha;
    gsBSplineBasis<T> m_basis_minus;

    gsBasis<T> & m_basis;

    mutable gsMapData<T> _tmp;

    bool m_isboundary;
    const index_t m_bfID, m_uv;


public:
    /// Shared pointer for gsNormalDerivBasis
    typedef memory::shared_ptr< gsNormalDerivBasis > Ptr;

    /// Unique pointer for gsNormalDerivBasis
    typedef memory::unique_ptr< gsNormalDerivBasis > uPtr;

    gsNormalDerivBasis(gsGeometry<T> & geo,
                 gsBSpline<T> basis_alpha,
                 gsBSplineBasis<T> basis_minus,
                 gsBasis<T> & basis,
                 bool isboundary,
                 const index_t bfID,
                 const index_t uv) :
            _geo(geo), m_basis_alpha(basis_alpha), m_basis_minus(basis_minus), m_basis(basis),
            m_isboundary(isboundary), m_bfID(bfID), m_uv(uv), _normalDerivBasis_piece(nullptr)
    {
        //_tmp.flags = NEED_JACOBIAN;
        //createMinusSpace(geo, basis, m_uv, m_basis_minus);
    }

    ~gsNormalDerivBasis() { delete _normalDerivBasis_piece; }

    GISMO_CLONE_FUNCTION(gsNormalDerivBasis)

    short_t domainDim() const {return 2;}

    short_t targetDim() const {return 1;}

    mutable gsNormalDerivBasis<T> * _normalDerivBasis_piece; // why do we need this?

    const gsFunction<T> & piece(const index_t k) const
    {
        //delete _normalDerivBasis_piece;
        _normalDerivBasis_piece = new gsNormalDerivBasis(*this);
        return *_normalDerivBasis_piece;
    }

    // Input is parametric coordinates of 2-D \a mp
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize( targetDim() , u.cols() );

        // tau/p
        gsBSplineBasis<T> bsp_temp = dynamic_cast<gsBSplineBasis<T> & >(m_basis.component(1-m_uv));

        T p = bsp_temp.degree();
        T tau_1 = bsp_temp.knots().at(p + 1); // p + 2

        gsMatrix<T> alpha, N_1, N_j_minus;

        if (!m_isboundary)
            m_basis_alpha.eval_into(u.row(m_uv),alpha); // 1-dir == PatchID
        else
            alpha.setOnes(1, u.cols());

        m_basis.component(1-m_uv).evalSingle_into(1,u.row(1-m_uv),N_1); // u

        m_basis_minus.evalSingle_into(m_bfID,u.row(m_uv),N_j_minus); // v

        if (!m_isboundary)
            result = (m_uv == 0 ? T(-1.0) : T(1.0)) * alpha.cwiseProduct(N_j_minus.cwiseProduct(N_1)) * tau_1 / p;
        else
            result = (m_uv == 0 ? T(-1.0) : T(1.0)) * alpha.cwiseProduct(N_j_minus.cwiseProduct(N_1));
    }

};


template <class T>
class gsVertexBasis : public gismo::gsFunction<T>
{

protected:
    const gsGeometry<T> &   m_geo;
    gsBasis<T> & m_basis;

    std::vector<gsBSpline<T>>            m_alpha;
    std::vector<gsBSpline<T>>            m_beta;

    std::vector<gsBSplineBasis<T>>       m_basis_plus;
    std::vector<gsBSplineBasis<T>>       m_basis_minus;

    const gsMatrix<T> m_Phi;
    const std::vector<bool> m_kindOfEdge;

    const index_t m_bfID;

    mutable gsMapData<T> _tmp;

public:
    /// Shared pointer for gsVertexBasis
    typedef memory::shared_ptr< gsVertexBasis > Ptr;

    /// Unique pointer for gsVertexBasis
    typedef memory::unique_ptr< gsVertexBasis > uPtr;

    gsVertexBasis(const gsGeometry<T> &   geo,
                  gsBasis<T> & basis,
                  std::vector<gsBSpline<T>> alpha,
                  std::vector<gsBSpline<T>> beta,
                  std::vector<gsBSplineBasis<T>> basis_plus,
                  std::vector<gsBSplineBasis<T>> basis_minus,
                  const gsMatrix<T> Phi,
                  const std::vector<bool> kindOfEdge,
                  const index_t bfID
            ) : m_geo(geo), m_basis(basis), m_alpha(alpha), m_beta(beta), m_basis_plus(basis_plus), m_basis_minus(basis_minus),
            m_Phi(Phi), m_kindOfEdge(kindOfEdge), m_bfID(bfID),
            _vertexBasis_piece(nullptr)
    {

    }

    ~gsVertexBasis() { delete _vertexBasis_piece; }

    GISMO_CLONE_FUNCTION(gsVertexBasis)

    short_t domainDim() const {return 2;}

    short_t targetDim() const {return 1;}

    mutable gsVertexBasis<T> * _vertexBasis_piece; // why do we need this?

    const gsFunction<T> & piece(const index_t k) const
    {
        //delete _vertexBasis_piece;
        _vertexBasis_piece = new gsVertexBasis(*this);
        return *_vertexBasis_piece;
    }

    // Input is parametric coordinates of 2-D \a mp
    // TODO IMPROVE THE FUNCTION
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize( targetDim() , u.cols() );
        result.setZero();

        // Computing c, c+ and c-
        // Point zero
        gsMatrix<T> zero;
        zero.setZero(2,1);

        std::vector<gsMatrix<T>> c_0, c_1;
        std::vector<gsMatrix <T>> c_0_plus, c_1_plus, c_2_plus;
        std::vector<gsMatrix <T>> c_0_plus_deriv, c_1_plus_deriv, c_2_plus_deriv;
        std::vector<gsMatrix <T>> c_0_minus, c_1_minus;
        for (index_t i = 0; i < 2; i++) // i == 0 == u , i == 1 == v
        {
            gsMatrix<T> b_0, b_1;
            gsMatrix<T> b_0_plus, b_1_plus, b_2_plus;
            gsMatrix<T> b_0_plus_deriv, b_1_plus_deriv, b_2_plus_deriv;
            gsMatrix<T> b_0_minus, b_1_minus;

            //gsBSplineBasis<T> bsp_temp = dynamic_cast<gsBSplineBasis<T> & >(m_basis.component(i));
            //gsBSplineBasis<T> bsp_temp_2 = dynamic_cast<gsBSplineBasis<T> & >(m_basis.component(i));
            //T p = bsp_temp.degree();
            //T h_geo = bsp_temp.knots().at(p + 1);
            //T h_geo_2 = bsp_temp_2.knots().at(p + 1);

            m_basis.component(i).evalSingle_into(0, u.row(i),b_0); // first
            m_basis.component(i).evalSingle_into(1, u.row(i),b_1); // second

            m_basis_plus[i].evalSingle_into(0, u.row(i),b_0_plus);
            m_basis_plus[i].evalSingle_into(1, u.row(i),b_1_plus);
            m_basis_plus[i].evalSingle_into(2, u.row(i),b_2_plus);

            m_basis_plus[i].derivSingle_into(0, u.row(i),b_0_plus_deriv);
            m_basis_plus[i].derivSingle_into(1, u.row(i),b_1_plus_deriv);
            m_basis_plus[i].derivSingle_into(2, u.row(i),b_2_plus_deriv);

            m_basis_minus[i].evalSingle_into(0, u.row(i),b_0_minus);
            m_basis_minus[i].evalSingle_into(1, u.row(i),b_1_minus);

            gsMatrix<T> b_1_0, b_1_minus_0;
            m_basis.component(i).derivSingle_into(1, zero.row(i),b_1_0);
            m_basis_minus[i].derivSingle_into(1, zero.row(i),b_1_minus_0);

            T factor_b_1 = 1.0/b_1_0(0,0);
            c_0.push_back(b_0 + b_1);
            c_1.push_back(factor_b_1 * b_1);

            T factor_b_1_minus = 1.0/b_1_minus_0(0,0);
            c_0_minus.push_back(b_0_minus + b_1_minus);
            c_1_minus.push_back(factor_b_1_minus * b_1_minus);

            gsMatrix<T> der_b_1_plus_0, der2_b_1_plus_0, der2_b_2_plus_0;
            m_basis_plus[i].derivSingle_into(1, zero.row(i), der_b_1_plus_0);
            m_basis_plus[i].deriv2Single_into(1, zero.row(i), der2_b_1_plus_0);
            m_basis_plus[i].deriv2Single_into(2, zero.row(i), der2_b_2_plus_0);

            T factor_c_1_plus = 1.0/der_b_1_plus_0(0,0);
            T factor2_c_1_plus = -der2_b_1_plus_0(0,0)/(der_b_1_plus_0(0,0)*der2_b_2_plus_0(0,0));
            T factor_c_2_plus = 1.0/der2_b_2_plus_0(0,0);

            c_0_plus.push_back(b_0_plus + b_1_plus + b_2_plus);
            c_1_plus.push_back(factor_c_1_plus * b_1_plus + factor2_c_1_plus * b_2_plus);
            c_2_plus.push_back(factor_c_2_plus * b_2_plus );

            c_0_plus_deriv.push_back(b_0_plus_deriv + b_1_plus_deriv + b_2_plus_deriv);
            c_1_plus_deriv.push_back(factor_c_1_plus * b_1_plus_deriv + factor2_c_1_plus * b_2_plus_deriv);
            c_2_plus_deriv.push_back(factor_c_2_plus * b_2_plus_deriv);
        }

        std::vector<gsMatrix<T>> alpha, beta, alpha_0, beta_0, alpha_deriv, beta_deriv;
        gsMatrix < T > temp_mat;
        if (m_kindOfEdge[0])
        {
            m_alpha[0].eval_into(u.row(0),temp_mat); // 1-dir == PatchID
            alpha.push_back(temp_mat); // u

            m_alpha[0].eval_into(zero.row(0),temp_mat); // 1-dir == PatchID
            alpha_0.push_back(temp_mat); // u

            m_alpha[0].deriv_into(zero.row(0),temp_mat); // 1-dir == PatchID
            alpha_deriv.push_back(temp_mat); // u

            m_beta[0].eval_into(u.row(0),temp_mat); // 1-dir == PatchID
            beta.push_back(temp_mat); // u

            m_beta[0].eval_into(zero.row(0),temp_mat); // 1-dir == PatchID
            beta_0.push_back(temp_mat); // u

            m_beta[0].deriv_into(zero.row(0),temp_mat); // 1-dir == PatchID
            beta_deriv.push_back(temp_mat); // u
        }
        else
        {
            temp_mat.setOnes(1, u.cols());
            alpha.push_back(temp_mat); // u

            temp_mat.setOnes(1, zero.cols());
            alpha_0.push_back(temp_mat); // u

            temp_mat.setZero(1, zero.cols());
            alpha_deriv.push_back(temp_mat); // u

            temp_mat.setZero(1, u.cols());
            beta.push_back(temp_mat); // u

            temp_mat.setZero(1, zero.cols());
            beta_0.push_back(temp_mat); // u

            temp_mat.setZero(1, zero.cols());
            beta_deriv.push_back(temp_mat); // u
        }

        if (m_kindOfEdge[1]) {
            m_alpha[1].eval_into(u.row(1), temp_mat); // 1-dir == PatchID
            alpha.push_back(temp_mat); // v

            m_alpha[1].eval_into(zero.row(0), temp_mat); // 1-dir == PatchID
            alpha_0.push_back(temp_mat); // v

            m_alpha[1].deriv_into(zero.row(0), temp_mat); // 1-dir == PatchID
            alpha_deriv.push_back(temp_mat); // v

            m_beta[1].eval_into(u.row(1), temp_mat); // 1-dir == PatchID
            beta.push_back(temp_mat); // v

            m_beta[1].eval_into(zero.row(0), temp_mat); // 1-dir == PatchID
            beta_0.push_back(temp_mat); // v

            m_beta[1].deriv_into(zero.row(0), temp_mat); // 1-dir == PatchID
            beta_deriv.push_back(temp_mat); // v
        }
        else
        {
            temp_mat.setOnes(1, u.cols());
            alpha.push_back(temp_mat); // u

            temp_mat.setOnes(1, zero.cols());
            alpha_0.push_back(temp_mat); // u

            temp_mat.setZero(1, zero.cols());
            alpha_deriv.push_back(temp_mat); // u

            temp_mat.setZero(1, u.cols());
            beta.push_back(temp_mat); // u

            temp_mat.setZero(1, zero.cols());
            beta_0.push_back(temp_mat); // u

            temp_mat.setZero(1, zero.cols());
            beta_deriv.push_back(temp_mat); // u
        }

        // Geo data:
        gsMatrix<T> geo_jac, geo_der2;
        geo_jac = m_geo.jacobian(zero);
        geo_der2 = m_geo.deriv2(zero);



        // Compute dd^^(i_k) and dd^^(i_k-1)
        gsMatrix<T> dd_ik_plus, dd_ik_minus;
        gsMatrix<T> dd_ik_minus_deriv, dd_ik_plus_deriv;
        dd_ik_minus = -1.0/(alpha_0[0](0,0)) * (geo_jac.col(1) +
                                              beta_0[0](0,0) * geo_jac.col(0));

        dd_ik_plus = 1.0/(alpha_0[1](0,0)) * (geo_jac.col(0) +
                                            beta_0[1](0,0) * geo_jac.col(1));

        gsMatrix<T> geo_deriv2_12(m_geo.targetDim(),1), geo_deriv2_11(m_geo.targetDim(),1), geo_deriv2_22(m_geo.targetDim(),1);
        geo_deriv2_12.row(0) = geo_der2.row(2);
        geo_deriv2_12.row(1) = geo_der2.row(5);
        geo_deriv2_11.row(0) = geo_der2.row(0);
        geo_deriv2_11.row(1) = geo_der2.row(3);
        geo_deriv2_22.row(0) = geo_der2.row(1);
        geo_deriv2_22.row(1) = geo_der2.row(4);

        if(m_geo.parDim() + 1 == m_geo.targetDim()) // Surface
        {
            geo_deriv2_12.row(2) = m_geo.deriv2(zero).row(8);
            geo_deriv2_11.row(2) = m_geo.deriv2(zero).row(6);
            geo_deriv2_22.row(2) = m_geo.deriv2(zero).row(7);
        }

        gsMatrix<T> alpha_squared_u = alpha_0[0]*alpha_0[0];
        gsMatrix<T> alpha_squared_v = alpha_0[1]*alpha_0[1];

        dd_ik_minus_deriv = -1.0/(alpha_squared_u(0,0)) * // N^2
                            ((geo_deriv2_12 + (beta_deriv[0](0,0) * geo_jac.col(0) +
                                               beta_0[0](0,0) * geo_deriv2_11))*alpha_0[0](0,0) -
                             (geo_jac.col(1) + beta_0[0](0,0) * geo_jac.col(0)) *
                             alpha_deriv[0](0,0));

        dd_ik_plus_deriv = 1.0/(alpha_squared_v(0,0)) *
                           ((geo_deriv2_12 + (beta_deriv[1](0,0) * geo_jac.col(1) +
                                              beta_0[1](0,0) * geo_deriv2_22))*alpha_0[1](0,0) -
                            (geo_jac.col(0) + beta_0[1](0,0) * geo_jac.col(1)) *
                            alpha_deriv[1](0,0));

        // Comupute d_(0,0)^(i_k), d_(1,0)^(i_k), d_(0,1)^(i_k), d_(1,1)^(i_k) ; i_k == 2
        std::vector<gsMatrix<T>> d_ik;
        d_ik.push_back(m_Phi.col(0));
        d_ik.push_back(m_Phi.block(1, 0, m_geo.targetDim(), 6).transpose() * geo_jac.col(0) ); // deriv into u
        d_ik.push_back(m_Phi.block(1, 0, m_geo.targetDim(), 6).transpose() * geo_jac.col(1) ); // deriv into v

        if(m_geo.parDim() + 1 == m_geo.targetDim()) // In the surface case the dimension of the second derivative vector is 9x1
        {
            d_ik.push_back( (geo_jac(0,0) * m_Phi.row(4).transpose() + geo_jac(1,0) * m_Phi.row(7).transpose() + geo_jac(2,0) * m_Phi.row(10).transpose()) * geo_jac(0,1) +
                            (geo_jac(0,0) * m_Phi.row(5).transpose() + geo_jac(1,0) * m_Phi.row(8).transpose() + geo_jac(2,0) * m_Phi.row(11).transpose()) * geo_jac(1,1) +
                            (geo_jac(0,0) * m_Phi.row(6).transpose() + geo_jac(1,0) * m_Phi.row(9).transpose() + geo_jac(2,0) * m_Phi.row(12).transpose()) * geo_jac(2,1) +
                                    m_Phi.block(1, 0, 1, 6).transpose() * geo_der2.row(2) +
                                    m_Phi.block(2, 0, 1, 6).transpose() * geo_der2.row(5) +
                                    m_Phi.block(3, 0, 1, 6).transpose() * geo_der2.row(8) );

        }
        else
        {
            d_ik.push_back((geo_jac(0, 0) * m_Phi.col(3) + geo_jac(1, 0) * m_Phi.col(4)) * geo_jac(0, 1) +
                           (geo_jac(0, 0) * m_Phi.col(4) + geo_jac(1, 0) * m_Phi.col(5)) * geo_jac(1, 1) +
                           m_Phi.block(0, 1, 6, 1) * geo_der2.row(2) +
                           m_Phi.block(0, 2, 6, 1) * geo_der2.row(5)); // Hessian
        }
        // Compute d_(*,*)^(il,ik)
        std::vector<gsMatrix<T>> d_ilik_minus, d_ilik_plus;
        d_ilik_minus.push_back(m_Phi.col(0));
        d_ilik_minus.push_back(m_Phi.block(1, 0, m_geo.targetDim(), 6).transpose() * geo_jac.col(0));

        if(m_geo.parDim() + 1 == m_geo.targetDim()) // In the surface case the dimension of the second derivative vector is 9x1
        {
            d_ilik_minus.push_back( (geo_jac(0,0) * m_Phi.row(4).transpose() + geo_jac(1,0) * m_Phi.row(7).transpose() + geo_jac(2,0) * m_Phi.row(10).transpose()) * geo_jac(0,0) +
                                    (geo_jac(0,0) * m_Phi.row(5).transpose() + geo_jac(1,0) * m_Phi.row(8).transpose() + geo_jac(2,0) * m_Phi.row(11).transpose()) * geo_jac(1,0) +
                                    (geo_jac(0,0) * m_Phi.row(6).transpose() + geo_jac(1,0) * m_Phi.row(9).transpose() + geo_jac(2,0) * m_Phi.row(12).transpose()) * geo_jac(2,0) +
                                    m_Phi.block(1, 0, 1, 6).transpose() * geo_der2.row(0) +
                                    m_Phi.block(2, 0, 1, 6).transpose() * geo_der2.row(3) +
                                    m_Phi.block(3, 0, 1, 6).transpose() * geo_der2.row(6) );
        }
        else
        {
            d_ilik_minus.push_back((geo_jac(0, 0) * m_Phi.col(3) + geo_jac(1, 0) * m_Phi.col(4)) * geo_jac(0, 0) +
                                   (geo_jac(0, 0) * m_Phi.col(4) + geo_jac(1, 0) * m_Phi.col(5)) * geo_jac(1, 0) +
                                   m_Phi.block(0, 1, 6, 1) * geo_der2.row(0) +
                                   m_Phi.block(0, 2, 6, 1) * geo_der2.row(3));
        }

        d_ilik_minus.push_back(m_Phi.block(1, 0, m_geo.targetDim(), 6).transpose() * dd_ik_minus);

        if(m_geo.parDim() + 1 == m_geo.targetDim()) // In the surface case the dimension of the second derivative vector is 9x1
        {
            d_ilik_minus.push_back( (geo_jac(0,0) * m_Phi.row(4).transpose() + geo_jac(1,0) * m_Phi.row(7).transpose() + geo_jac(2,0) * m_Phi.row(10).transpose()) * dd_ik_minus(0,0) +
                                    (geo_jac(0,0) * m_Phi.row(5).transpose() + geo_jac(1,0) * m_Phi.row(8).transpose() + geo_jac(2,0) * m_Phi.row(11).transpose()) * dd_ik_minus(1,0) +
                                    (geo_jac(0,0) * m_Phi.row(6).transpose() + geo_jac(1,0) * m_Phi.row(9).transpose() + geo_jac(2,0) * m_Phi.row(12).transpose()) * dd_ik_minus(2,0) +
                                    m_Phi.block(1, 0, 1, 6).transpose() * dd_ik_minus_deriv.row(0) +
                                    m_Phi.block(2, 0, 1, 6).transpose() * dd_ik_minus_deriv.row(1) +
                                    m_Phi.block(3, 0, 1, 6).transpose() * dd_ik_minus_deriv.row(2) );
        }
        else
        {
            d_ilik_minus.push_back((geo_jac(0, 0) * m_Phi.col(3) + geo_jac(1, 0) * m_Phi.col(4)) * dd_ik_minus(0, 0) +
                                   (geo_jac(0, 0) * m_Phi.col(4) + geo_jac(1, 0) * m_Phi.col(5)) * dd_ik_minus(1, 0) +
                                   m_Phi.block(0, 1, 6, 1) * dd_ik_minus_deriv.row(0) +
                                   m_Phi.block(0, 2, 6, 1) * dd_ik_minus_deriv.row(1));
        }

        d_ilik_plus.push_back(m_Phi.col(0));
        d_ilik_plus.push_back(m_Phi.block(1, 0, m_geo.targetDim(), 6).transpose() * geo_jac.col(1));

        if(m_geo.parDim() + 1 == m_geo.targetDim()) // In the surface case the dimension of the second derivative vector is 9x1
        {
            d_ilik_plus.push_back( (geo_jac(0,1) * m_Phi.row(4).transpose() + geo_jac(1,1) * m_Phi.row(7).transpose() + geo_jac(2,1) * m_Phi.row(10).transpose()) * geo_jac(0,1) +
                                   (geo_jac(0,1) * m_Phi.row(5).transpose() + geo_jac(1,1) * m_Phi.row(8).transpose() + geo_jac(2,1) * m_Phi.row(11).transpose()) * geo_jac(1,1) +
                                   (geo_jac(0,1) * m_Phi.row(6).transpose() + geo_jac(1,1) * m_Phi.row(9).transpose() + geo_jac(2,1) * m_Phi.row(12).transpose()) * geo_jac(2,1) +
                                   m_Phi.block(1, 0, 1, 6).transpose() * geo_der2.row(1) +
                                   m_Phi.block(2, 0, 1, 6).transpose() * geo_der2.row(4) +
                                   m_Phi.block(3, 0, 1, 6).transpose() * geo_der2.row(7) );
        }
        else
        {
            d_ilik_plus.push_back((geo_jac(0, 1) * m_Phi.col(3) + geo_jac(1, 1) * m_Phi.col(4)) * geo_jac(0, 1) +
                                  (geo_jac(0, 1) * m_Phi.col(4) + geo_jac(1, 1) * m_Phi.col(5)) * geo_jac(1, 1) +
                                  m_Phi.block(0, 1, 6, 1) * geo_der2.row(1) +
                                  m_Phi.block(0, 2, 6, 1) * geo_der2.row(4));
        }

        d_ilik_plus.push_back(m_Phi.block(1, 0, m_geo.targetDim(), 6).transpose()  * dd_ik_plus);

        if(m_geo.parDim() + 1 == m_geo.targetDim()) // In the surface case the dimension of the second derivative vector is 9x1
        {
            d_ilik_plus.push_back( (geo_jac(0,1) * m_Phi.row(4).transpose() + geo_jac(1,1) * m_Phi.row(7).transpose() + geo_jac(2,1) * m_Phi.row(10).transpose()) * dd_ik_plus(0,0) +
                                   (geo_jac(0,1) * m_Phi.row(5).transpose() + geo_jac(1,1) * m_Phi.row(8).transpose() + geo_jac(2,1) * m_Phi.row(11).transpose()) * dd_ik_plus(1,0) +
                                   (geo_jac(0,1) * m_Phi.row(6).transpose() + geo_jac(1,1) * m_Phi.row(9).transpose() + geo_jac(2,1) * m_Phi.row(12).transpose()) * dd_ik_plus(2,0) +
                                   m_Phi.block(1, 0, 1, 6).transpose() * dd_ik_plus_deriv.row(0) +
                                   m_Phi.block(2, 0, 1, 6).transpose() * dd_ik_plus_deriv.row(1) +
                                   m_Phi.block(3, 0, 1, 6).transpose() * dd_ik_plus_deriv.row(2) );
        }
        else
        {
            d_ilik_plus.push_back((geo_jac(0, 1) * m_Phi.col(3) + geo_jac(1, 1) * m_Phi.col(4)) * dd_ik_plus(0, 0) +
                                  (geo_jac(0, 1) * m_Phi.col(4) + geo_jac(1, 1) * m_Phi.col(5)) * dd_ik_plus(1, 0) +
                                  m_Phi.block(0, 1, 6, 1) * dd_ik_plus_deriv.row(0) +
                                  m_Phi.block(0, 2, 6, 1) * dd_ik_plus_deriv.row(1));
        }


        result = d_ilik_minus.at(0)(m_bfID,0) * (c_0_plus.at(0).cwiseProduct(c_0.at(1)) -
                                                   beta[0].cwiseProduct(c_0_plus_deriv.at(0).cwiseProduct(c_1.at(1)))) +
                        d_ilik_minus.at(1)(m_bfID,0) * (c_1_plus.at(0).cwiseProduct(c_0.at(1)) -
                                                   beta[0].cwiseProduct(c_1_plus_deriv.at(0).cwiseProduct(c_1.at(1)))) +
                        d_ilik_minus.at(2)(m_bfID,0) * (c_2_plus.at(0).cwiseProduct(c_0.at(1)) -
                                                   beta[0].cwiseProduct(c_2_plus_deriv.at(0).cwiseProduct(c_1.at(1)))) -
                        d_ilik_minus.at(3)(m_bfID,0) * alpha[0].cwiseProduct(c_0_minus.at(0).cwiseProduct(c_1.at(1))) -
                        d_ilik_minus.at(4)(m_bfID,0) * alpha[0].cwiseProduct(c_1_minus.at(0).cwiseProduct(c_1.at(1))); // f*_(ik-1,ik)

        //if (kindOfEdge[0])
        //rhsVals.at(i).setZero();

        //if (!kindOfEdge[1])
        result += d_ilik_plus.at(0)(m_bfID,0) * (c_0_plus.at(1).cwiseProduct(c_0.at(0)) -
                                                   beta[1].cwiseProduct(c_0_plus_deriv.at(1).cwiseProduct(c_1.at(0)))) +
                         d_ilik_plus.at(1)(m_bfID,0) * (c_1_plus.at(1).cwiseProduct(c_0.at(0)) -
                                                   beta[1].cwiseProduct(c_1_plus_deriv.at(1).cwiseProduct(c_1.at(0)))) +
                         d_ilik_plus.at(2)(m_bfID,0) * (c_2_plus.at(1).cwiseProduct(c_0.at(0)) -
                                                   beta[1].cwiseProduct(c_2_plus_deriv.at(1).cwiseProduct(c_1.at(0)))) +
                         d_ilik_plus.at(3)(m_bfID,0) * alpha[1].cwiseProduct(c_0_minus.at(1).cwiseProduct(c_1.at(0))) +
                         d_ilik_plus.at(4)(m_bfID,0) * alpha[1].cwiseProduct(c_1_minus.at(1).cwiseProduct(c_1.at(0))); // f*_(ik+1,ik)

        result -= d_ik.at(0)(m_bfID,0) * c_0.at(0).cwiseProduct(c_0.at(1)) + d_ik.at(2)(m_bfID,0) * c_0.at(0).cwiseProduct(c_1.at(1)) +
                         d_ik.at(1)(m_bfID,0) * c_1.at(0).cwiseProduct(c_0.at(1)) + d_ik.at(3)(m_bfID,0) * c_1.at(0).cwiseProduct(c_1.at(1)); // f*_(ik)


    }

};

}
