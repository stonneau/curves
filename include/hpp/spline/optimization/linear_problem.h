/**
* \file bezier_curve.h
* \brief class allowing to create a Bezier curve of dimension 1 <= n <= 3.
* \author Steve T.
* \version 0.1
* \date 06/17/2013
*/


#ifndef _CLASS_LINEAR_PROBLEM
#define _CLASS_LINEAR_PROBLEM

#include "hpp/spline/optimization/definitions.h"
#include "hpp/spline/optimization/details.h"
#include "hpp/spline/optimization/quadratic_cost.h"

#include <Eigen/Core>

namespace spline
{
namespace  optimization
{

template<typename Point, int Dim, typename Numeric>
problem<Point, Dim, Numeric> generate_problem(const problem_definition<Point, Dim, Numeric>& pDef)
{
    problem<Point, Dim, Numeric> prob;
    problem_data<Point, Dim, Numeric> pData = setup_control_points<Point, Dim, Numeric>(pDef);
    initInequalityMatrix<Point, Dim, Numeric>(pDef,pData,prob);
    prob.cost = compute_jerk_cost<Point, Dim, Numeric>(pData);
    return prob;
}
} // namespace optimization
} // namespace spline
#endif //_CLASS_LINEAR_PROBLEM

