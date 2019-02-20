/**
* \file bezier_curve.h
* \brief class allowing to create a Bezier curve of dimension 1 <= n <= 3.
* \author Steve T.
* \version 0.1
* \date 06/17/2013
*/


#ifndef _CLASS_QUADRATIC_COST
#define _CLASS_QUADRATIC_COST

#include "hpp/spline/optimization/definitions.h"
#include "hpp/spline/optimization/details.h"

#include <Eigen/Core>

namespace spline
{
namespace  optimization
{

template<typename Point, int Dim, typename Numeric>
quadratic_variable<Numeric> compute_derivative_cost_internal
                    (const problem_data<Point, Dim, Numeric>& pData, const std::size_t num_derivate)
{
    typedef bezier_curve<Numeric, Numeric, Dim, true,linear_variable<Dim, Numeric> > bezier_t;
    typedef typename bezier_t::t_point_t t_point_t;
    typedef typename t_point_t::const_iterator cit_point_t;
    bezier_t acc = pData.bezier->compute_derivate(num_derivate);
    const t_point_t& wps = acc.waypoints();
    return bezier_product<Point, Dim, Numeric, cit_point_t>
            (pData,wps.begin(),wps.end(),wps.begin(),wps.end());
}

template<typename Point, int Dim, typename Numeric>
quadratic_variable<Numeric> compute_derivative_cost
                    (const problem_data<Point, Dim, Numeric>& pData, const derivative_flag flag)
{
    std::size_t size = (std::size_t)(flag) ;
    return compute_derivative_cost_internal<Point,Dim,Numeric>(pData, size);
}

} // namespace optimization
} // namespace spline
#endif //_CLASS_QUADRATIC_COST

