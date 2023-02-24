/** @file gsMPBESMapB2D.h

    @brief Provides declaration of Basis abstract interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

// This class gsMPBESMapB2D has the sole purpose of creating a mapping of the type gsWeightMapper

#pragma once

#include <gsCore/gsBasis.h>
#include <gsCore/gsBoxTopology.h>
#include <gsNurbs/gsKnotVector.h>

#include <gsUnstructuredSplines/src/gsMPBESBasis.h>
#include <gsUnstructuredSplines/src/gsMPBESMapTensor.h>

#define TO_BSPLINE(x) static_cast<const gsTensorBSplineBasis<d,T> *>(x)

namespace gismo
{

/** @brief
      A univariate Lagrange basis.

      \tparam T coefficient type

      \ingroup basis
  */

template<short_t d,class T>
class gsMPBESMapB2D : public gsMPBESMapTensor<d,T>
{
    //static const index_t d = 2;
private:
    typedef gsBasis<T> BasisType;
    typedef gsMPBESMapTensor<d,T> Base;
public:
    gsMPBESMapB2D(index_t incrSmoothnessDegree, gsBoxTopology * topol, gsMPBESBasis<d,T> * basis) :
        Base(incrSmoothnessDegree,topol,basis)
    { }

    ~gsMPBESMapB2D() { }

private:
    using Base::m_incrSmoothnessDegree;
    using Base::m_topol;
    using Base::m_basis;
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
        bool consistent = true;
//        for(index_t i=0;i<m_local_to_global.rows();i++)
//        {
//            IndexContainer globals = _local2global(i);
//            index_t size=globals.size();
//            if(size!=1&&size!=2&&(size!=degree(0)&&size!=degree(1))&&size>degree(0)&&size>degree(1))
//            {
////                std::cout << "Error in " << i << ".th row:";
////                print(std::cout,i);
////                std::cout << std::endl;
////                print(std::cout,120);
////                GISMO_ERROR("mapping wrong");
//            }
//        }
//        for(index_t i=0;i<m_local_to_global.cols();i++)
//        {
//            IndexContainer locals = _global2local(i);
//            index_t size=locals.size();
//            if(size!=1&&size!=4)
//            {
//                std::cout << "Error in " << i << ".th col:";
//                print(std::cout,-1,i);
//                std::cout << std::endl;
//                //GISMO_ERROR("mapping wrong");
//            }
//        }
        return consistent;
    }

    void _finalize() const
    {
        gsSparseMatrix<T> mat=m_mapper->asMatrix();
        mat.conservativeResize(mat.rows(),m_global);
        delete m_mapper;
        m_mapper=new gsWeightMapper<T>(mat);
        m_mapper->optimize(gsWeightMapper<T>::optTargetToSource);
    }

    void _setMappingOfPatch(index_t const patch) const
    {
        _setTensorMappingOfPatch(patch);
    }

    index_t getDistanceOfVertex(const patchCorner& pc,const patchSide& ps) const
    {
        std::vector<T> endpoints;
        T parametricDistance = m_basis->getParametricDistanceOfVertex(pc,ps);
        index_t patch = ps.patch;

        if(math::almostEqual<T>(parametricDistance,0.0))
            return 0;

        index_t deg = m_basis->degree(patch,1-(ps.direction()));
        gsKnotVector<T> knots = TO_BSPLINE(&(m_basis->getBase(patch)))->knots(1-(ps.direction()));
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
    {
        return TO_BSPLINE(&(m_basis->getBase(patch)))->knots(par);
    }

    index_t _getParMax(index_t patch,bool par) const
    {
        return TO_BSPLINE(&(m_basis->getBase(patch)))->size(par)-1;
    }

private:
    //////////////////////////////////////////////////
    // functions for working with Indexes
    //////////////////////////////////////////////////

    bool _getLocalIndex_into(index_t const patch,index_t const u,index_t const v,index_t & localIndex) const
    {
        localIndex=_getLocalIndex(patch,u,v);
        return true;
    }

    index_t _getLocalIndex(index_t const patch,index_t u, index_t v) const
    {
        return _getLocalIndex(patch,_getPatchIndex(patch,u,v));
    }

    index_t _getPatchIndex(index_t const patch,boxSide const side,bool const flag) const
    {
        index_t u_amount=_getParMax(patch,0);
        index_t v_amount=_getParMax(patch,1);
        if(side.direction())
            if(side.parameter())
                return _getPatchIndex(patch,flag?u_amount:0,v_amount);
            else
                return _getPatchIndex(patch,flag?u_amount:0,0);
        else
            if(side.parameter())
                return _getPatchIndex(patch,u_amount,flag?v_amount:0);
            else
                return _getPatchIndex(patch,0,flag?v_amount:0);
    }

    index_t _getPatchIndex(index_t const patch,index_t const u,index_t const v) const
    {
        gsVector<index_t,d> vec;
        vec(0)=u,vec(1)=v;
        return TO_BSPLINE(&(m_basis->getBase(patch)))->index(vec);
    }

    index_t _getPar(index_t localIndex,bool par) const
    {
        return _getPar(_getPatch(localIndex),_getPatchIndex(localIndex),par);
    }

    index_t _getPar(index_t patch,index_t patchIndex, bool par) const
    {
        gsVector<index_t,d> vec = TO_BSPLINE(&(m_basis->getBase(patch)))->tensorIndex(patchIndex);
        return vec(par);
    }
};

}

#undef TO_BSPLINE
