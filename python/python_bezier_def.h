#include "hpp/spline/bezier_curve.h"
#include "hpp/spline/bezier_curve.h"
#include "hpp/spline/polynom.h"
#include "hpp/spline/exact_cubic.h"
#include "hpp/spline/spline_deriv_constraint.h"
#include "hpp/spline/curve_constraint.h"
#include "hpp/spline/bezier_polynom_conversion.h"
#include "hpp/spline/bernstein.h"
#include "hpp/spline/linear_variable.h"
#include "hpp/spline/optimization/linear_problem.h"


#include "python_definitions.h"

#include <vector>

#ifndef _PYTHON_BEZIER_DEF
#define _PYTHON_BEZIER_DEF


namespace spline
{

typedef spline::bezier_curve  <real, real, 3, true, point_t> bezier_t;
typedef spline::bezier_curve  <real, real, 6, true, point6_t> bezier6_t;
typedef spline::polynom  <real, real, 3, true, point_t, t_point_t> polynom_t;
typedef spline::exact_cubic  <real, real, 3, true, point_t, t_point_t> exact_cubic_t;
typedef polynom_t::coeff_t coeff_t;
typedef std::pair<real, point_t> waypoint_t;
typedef std::vector<waypoint_t, Eigen::aligned_allocator<point_t> > t_waypoint_t;

typedef spline::Bern<double> bernstein_t;


typedef spline::spline_deriv_constraint  <real, real, 3, true, point_t, t_point_t> spline_deriv_constraint_t;
typedef spline::curve_constraints<point_t> curve_constraints_t;
typedef spline::curve_constraints<point6_t> curve_constraints6_t;
/*** TEMPLATE SPECIALIZATION FOR PYTHON ****/


template <typename Bezier, typename PointList, typename T_Point>
Bezier* wrapBezierConstructorTemplate(const PointList& array, const real ub =1.)
{
    T_Point asVector = vectorFromEigenArray<PointList, T_Point>(array);
    return new Bezier(asVector.begin(), asVector.end(), ub);
}

template <typename Bezier, typename PointList, typename T_Point, typename CurveConstraints>
Bezier* wrapBezierConstructorConstraintsTemplate(const PointList& array, const CurveConstraints& constraints, const real ub =1.)
{
    T_Point asVector = vectorFromEigenArray<PointList, T_Point>(array);
    return new Bezier(asVector.begin(), asVector.end(), constraints, ub);
}

/*3D constructors */
bezier_t* wrapBezierConstructor(const point_list_t& array);

bezier_t* wrapBezierConstructorBounds(const point_list_t& array, const real ub);

bezier_t* wrapBezierConstructorConstraints
(const point_list_t& array, const curve_constraints_t& constraints);

bezier_t* wrapBezierConstructorBoundsConstraints
(const point_list_t& array, const curve_constraints_t& constraints, const real ub);

/*END 3D constructors */
/*6D constructors */
bezier6_t* wrapBezierConstructor6(const point_list6_t& array);

bezier6_t* wrapBezierConstructorBounds6(const point_list6_t& array, const real ub);

bezier6_t* wrapBezierConstructor6Constraints
(const point_list6_t& array, const curve_constraints6_t& constraints);

bezier6_t* wrapBezierConstructorBounds6Constraints
(const point_list6_t& array, const curve_constraints6_t& constraints, const real ub);
/*END 6D constructors */

polynom_t* wrapSplineConstructor(const coeff_t& array);


t_waypoint_t getWayPoints(const coeff_t& array, const time_waypoints_t& time_wp);

template <typename BezierType, int dim>
Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> wayPointsToList(const BezierType& self)
{
    typedef typename BezierType::t_point_t t_point;
    typedef typename BezierType::t_point_t::const_iterator cit_point;
    const t_point& wps = self.waypoints();
    Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> res (dim, wps.size());
    int col = 0;
    for(cit_point cit = wps.begin(); cit != wps.end(); ++cit, ++col)
        res.block<dim,1>(0,col) = *cit;
    return res;
}

exact_cubic_t* wrapExactCubicConstructor(const coeff_t& array, const time_waypoints_t& time_wp);

spline_deriv_constraint_t* wrapSplineDerivConstraint(
        const coeff_t& array, const time_waypoints_t& time_wp, const curve_constraints_t& constraints);

spline_deriv_constraint_t* wrapSplineDerivConstraintNoConstraints
(const coeff_t& array, const time_waypoints_t& time_wp);

point_t get_init_vel(const curve_constraints_t& c);

point_t get_init_acc(const curve_constraints_t& c);

point_t get_end_vel(const curve_constraints_t& c);

point_t get_end_acc(const curve_constraints_t& c);

void set_init_vel(curve_constraints_t& c, const point_t& val);

void set_init_acc(curve_constraints_t& c, const point_t& val);

void set_end_vel(curve_constraints_t& c, const point_t& val);

void set_end_acc(curve_constraints_t& c, const point_t& val);

} //namespace spline.


#endif //_PYTHON_BEZIER_DEF
