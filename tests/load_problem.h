
#ifndef _CLASS_LOAD_TEST_PROBLEMS
#define _CLASS_LOAD_TEST_PROBLEMS


#include "hpp/spline/exact_cubic.h"
#include "hpp/spline/bezier_curve.h"
#include "hpp/spline/polynom.h"
#include "hpp/spline/spline_deriv_constraint.h"
#include "hpp/spline/helpers/effector_spline.h"
#include "hpp/spline/helpers/effector_spline_rotation.h"
#include "hpp/spline/bezier_polynom_conversion.h"
#include "hpp/spline/optimization/linear_problem.h"
#include "hpp/spline/optimization/details.h"

namespace spline {

typedef Eigen::Vector3d point_t;
typedef std::vector<point_t,Eigen::aligned_allocator<point_t> >  t_point_t;
typedef polynom  <double, double, 3, true, point_t, t_point_t> polynom_t;
typedef exact_cubic <double, double, 3, true, point_t> exact_cubic_t;
typedef spline_deriv_constraint <double, double, 3, true, point_t> spline_deriv_constraint_t;
typedef bezier_curve  <double, double, 3, true, point_t> bezier_curve_t;
typedef spline_deriv_constraint_t::spline_constraints spline_constraints_t;
typedef std::pair<double, point_t> Waypoint;
typedef std::vector<Waypoint> T_Waypoint;


typedef Eigen::Matrix<double,1,1> point_one;
typedef polynom<double, double, 1, true, point_one> polynom_one;
typedef exact_cubic   <double, double, 1, true, point_one> exact_cubic_one;
typedef std::pair<double, point_one> WaypointOne;
typedef std::vector<WaypointOne> T_WaypointOne;

namespace optimization
{
typedef curve_constraints<point_t> constraint_linear;
typedef linear_variable<3, double> linear_variable_t;
typedef variables<linear_variable_t> variables_t;
typedef std::pair<std::size_t, std::size_t >   pair_size_t;
typedef std::pair<variables_t, pair_size_t > var_pair_t;
typedef problem_data<point_t, 3, double> problem_data_t;

}

}

#endif //_CLASS_LOAD_TEST_PROBLEMS
