/** @file gsMPBESMapHB2D.h

    @brief Provides declaration of Basis abstract interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

// This class gsMPBESMapHB2D has the sole purpose of creating a mapping of the type gsWeightMapper

#pragma once

#include <gsCore/gsBasis.h>
#include <gsCore/gsBoxTopology.h>
#include <gsNurbs/gsKnotVector.h>

#include <gsUnstructuredSplines/src/gsMPBESBasis.h>
#include <gsUnstructuredSplines/src/gsMPBESMapTensor.h>

#define TO_HTENSOR(x) static_cast<const gsHTensorBasis<d,T> *>(x)
#define TO_BSPLINE(x) static_cast<const gsTensorBSplineBasis<d,T> *>(x)

namespace gismo
{

/** @brief
      A univariate Lagrange basis.

      \tparam T coefficient type

      \ingroup basis
  */

template<short_t d,class T>
class gsMPBESMapHB2D : public gsMPBESMapTensor<d,T>
{
    //static const index_t d = 2;
private:
    typedef gsBasis<T> BasisType;
    typedef gsMPBESMapTensor<d,T> Base;
public:
    gsMPBESMapHB2D(index_t incrSmoothnessDegree, gsBoxTopology * topol, gsMPBESBasis<d,T> * basis) :
        Base(incrSmoothnessDegree,topol,basis)
    { }

    ~gsMPBESMapHB2D() { }

private:
    using Base::m_basis;
    using Base::m_incrSmoothnessDegree;
    using Base::m_topol;
    using Base::m_mapper;
    using Base::m_global;
    using Base::_setTensorMappingOfPatch;
    using Base::_getPatch;
    using Base::_getPatchIndex;
    using Base::_getLocalIndex;

private:
    //////////////////////////////////////////////////
    // general help functions for checking, finalizing and building of the mapping
    //////////////////////////////////////////////////

    bool _checkMapping() const
    {
        return true;
    }

    void _finalize() const
    {
        m_level = _getMaxLevel();
        gsSparseMatrix<T> mat=m_mapper->asMatrix();
        mat.conservativeResize(mat.rows(),m_global);
        delete m_mapper;
        m_mapper=new gsWeightMapper<T>(mat);
        m_mapper->optimize(gsWeightMapper<T>::optTargetToSource);
    }

    void _setMappingOfPatch(index_t const patch) const
    {
        m_level=0;
        for(index_t i = 0;i<=_getMaxLevel();i++)
        {
            if(m_level<=_getMaxLevel(patch))
                _setTensorMappingOfPatch(patch);
            m_level++;
        }
    }

    index_t _getMaxLevel() const
    {
        index_t level = 0;
        for (size_t i = 0; i < m_basis->nPatches(); i++)
        {
            level = std::max(level, _getMaxLevel(i));
        }
        return level;
    }

    index_t _getMaxLevel(index_t patch) const
    {
        return TO_HTENSOR(&(m_basis->getBase(patch)))->maxLevel();
    }


    index_t getDistanceOfVertex(const patchCorner& pc,const patchSide& ps) const
    {
        std::vector<T> endpoints;
        T parametricDistance = m_basis->getParametricDistanceOfVertex(pc,ps);
        if(math::almostEqual<T>(parametricDistance,0.0))
            return 0;
        index_t patch = ps.patch;
        index_t deg = m_basis->degree(patch,1-ps.direction());
        gsTensorBSplineBasis<d,T>* basis = TO_HTENSOR(&(m_basis->getBase(patch)))->getBases()[m_level];
        gsKnotVector<T> knots = basis->knots(1-(ps.direction()));
        gsVector<bool> pars;
        pc.parameters_into(d,pars);
        if(!pars(1-ps.direction()))
            knots.reverse();
        for(size_t i = deg+1;i<knots.size();i++)
            endpoints.push_back(knots.at(i));
        std::sort(endpoints.begin(),endpoints.end());
        index_t nr=0;
        for(;nr<(index_t)(endpoints.size());nr++)
            if(math::almostEqual<T>(endpoints[nr],parametricDistance)||endpoints[nr]>=parametricDistance)
                break;
        return nr+1;
    }

private:
    //////////////////////////////////////////////////
    // functions calculating the weights for the mapping
    //////////////////////////////////////////////////

    gsKnotVector<T> _getKnotVector(index_t const patch,index_t const par) const
    {// todo: remove
        gsKnotVector<T> kvComp = TO_HTENSOR(&(m_basis->getBase(patch)))->getBases()[m_level]->knots(par);
        return gsKnotVector<T>(kvComp);
    }

    index_t _getParMax(index_t patch,bool par) const
    {
        return TO_HTENSOR(&(m_basis->getBase(patch)))->getBases()[m_level]->size(par)-1;
    }

private:
    //////////////////////////////////////////////////
    // functions for working with Indexes
    //////////////////////////////////////////////////
    // localIndex = hierachical index of hierachical splines collected over all patches.

    bool _getLocalIndex_into(index_t const patch,index_t const u,index_t const v,index_t & localIndex) const
    {
        index_t patchIndex = _getPatchIndex(patch,u,v);
        localIndex=_getLocalIndex(patch,patchIndex);
        if(patchIndex==-1)
            return false;
        else
            return true;
    }

    index_t _getLocalIndex(index_t const patch,index_t u, index_t v) const
    {
        return _getLocalIndex(patch,_getPatchIndex(patch,u,v));
    }

    index_t _getPatchIndex(index_t const patch,boxSide const side,bool const flag) const
    {
        index_t u,v;
        index_t level = 0;
        gsVector<index_t,d> vec;
        index_t index, patchindex;
        do
        {
            if(level>(index_t)(TO_HTENSOR(&(m_basis->getBase(patch)))->maxLevel()))
                GISMO_ERROR("did not find the patchindex");

            const index_t u_amount=TO_HTENSOR(&(m_basis->getBase(patch)))->tensorLevel(level).size(0);
            const index_t v_amount=TO_HTENSOR(&(m_basis->getBase(patch)))->tensorLevel(level).size(1);
            if(side.direction())
                if(side.parameter())
                {
                    u=flag ? (u_amount-1) : 0;
                    v=v_amount-1;
                }
                else
                {
                    u=flag ? (u_amount-1) :0;
                    v=0;
                }
            else
                if(side.parameter())
                {
                    u=u_amount-1;
                    v=flag ? (v_amount-1) : 0;
                }
                else
                {
                    u=0;
                    v=flag ? (v_amount-1) : 0;
                }
            vec(0)=u,vec(1)=v;
            index=TO_HTENSOR(&(m_basis->getBase(patch)))->getBases()[level]->index(vec);
            patchindex=TO_HTENSOR(&(m_basis->getBase(patch)))->flatTensorIndexToHierachicalIndex(index,level);
            level++;
        }while(-1==patchindex);
        return patchindex;
    }

    index_t _getPatchIndex(index_t const patch,index_t const u,index_t const v) const
    {
        index_t index=_getTensorIndex(patch,u,v);
        return TO_HTENSOR(&(m_basis->getBase(patch)))->flatTensorIndexToHierachicalIndex(index,m_level);
    }

    // tensorIndex = flat tensor index of one level of hierarchical splines
    index_t _getTensorIndex(index_t const patch,index_t const u, index_t const v) const
    {
        gsVector<index_t,d> vec;
        vec(0)=u,vec(1)=v;
        return TO_HTENSOR(&(m_basis->getBase(patch)))->getBases()[m_level]->index(vec);
    }

    index_t _getTensorIndex(index_t const patch,index_t const patchIndex) const
    {
        index_t combIndex = TO_HTENSOR(&(m_basis->getBase(patch)))->flatTensorIndexOf(patchIndex,m_level);
//        for(index_t i = 0;i<m_level;i++)
//            combIndex-= TO_HTENSOR((*m_bases)[patch])->m_bases[i]->size();
        return combIndex;
    }

    index_t _getPar(index_t localIndex,bool par) const
    {
        index_t patch = _getPatch(localIndex);
        index_t patchIndex = _getPatchIndex(localIndex);
        return _getPar(patch,_getTensorIndex(patch,patchIndex),par);
    }

    index_t _getPar(index_t patch,index_t tensorIndex, bool par) const
    {
        gsVector<index_t,d> vec = TO_HTENSOR(&(m_basis->getBase(patch)))->getBases()[m_level]->tensorIndex(tensorIndex);
        GISMO_ASSERT(static_cast<index_t>(vec(par))<TO_HTENSOR(&(m_basis->getBase(patch)))->getBases()[m_level]->size(par),"wrong tensorIndex");
        GISMO_ASSERT(static_cast<index_t>(vec(!par))<TO_HTENSOR(&(m_basis->getBase(patch)))->getBases()[m_level]->size(!par),"wrong tensorIndex");
        return vec(par);
    }

private:
    mutable index_t m_level; // used in the construction to determine the level, after the construction it is used to say the max_level of the H-Splines
};

}

#undef TO_BSPLINE
#undef TO_HTENSOR
