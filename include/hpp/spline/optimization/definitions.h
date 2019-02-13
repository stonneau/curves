/**
* \file definitions.h
* \brief utils for defining optimization problems
* \author Steve T.
* \version 0.1
* \date 06/17/2013
*/


#ifndef _CLASS_DEFINITIONS_H
#define _CLASS_DEFINITIONS_H


#include <hpp/spline/bezier_curve.h>
#include <hpp/spline/linear_variable.h>
#include <hpp/spline/curve_constraint.h>

#include <Eigen/StdVector>

namespace spline
{
namespace  optimization
{

enum constraint_flag{
    INIT_POS  = 0x001,
    INIT_VEL  = 0x002,
    INIT_ACC  = 0x004,
    END_POS   = 0x008,
    END_VEL   = 0x010,
    END_ACC   = 0x020,
    ALL       = 0x03f,
    NONE      = 0x100
};

template<typename Point, int Dim, typename Numeric>
struct problem_definition
{
    typedef Point  point_t;
    typedef Numeric  num_t;
    typedef curve_constraints<point_t> curve_constraints_t;
    typedef Eigen::Matrix< num_t , Eigen::Dynamic , 1> vectorx_t;
    typedef Eigen::Matrix< num_t , Eigen::Dynamic , Eigen::Dynamic> matrix_x_t;
    typedef Eigen::Matrix< num_t , Eigen::Dynamic , Dim> matrix_dim_t;
    typedef std::vector<matrix_dim_t, Eigen::aligned_allocator<matrix_dim_t> > T_matrix_dim_t;
    typedef std::vector<vectorx_t, Eigen::aligned_allocator<vectorx_t> > T_vectorx_t;
    typedef typename T_matrix_dim_t::const_iterator CIT_matrix_dim_t;
    typedef typename T_vectorx_t::const_iterator CIT_vectorx_t;

    problem_definition()
        : flag(NONE)
        , start(point_t())
        , end(point_t())
        , curveConstraints()
        , degree(5)
        , totalTime(1.)
        , splitTimes_(vectorx_t::Zero(0)) {}


    constraint_flag flag;
    point_t start;
    point_t end;
    curve_constraints_t curveConstraints;
    std::size_t degree;
    num_t totalTime;
    vectorx_t splitTimes_;
    T_matrix_dim_t inequalityMatrices_; // must be of size (splitTimes_ + 1)
    T_vectorx_t    inequalityVectors_;  // must be of size (splitTimes_ + 1)
};

} // namespace optimization
} // namespace spline
#endif //_CLASS_DEFINITIONS_H

