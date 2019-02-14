#include "python_bezier_def.h"

#include <Eigen/Core>

namespace spline
{
/*3D constructors */
bezier_t* wrapBezierConstructor(const point_list_t& array)
{
    return wrapBezierConstructorTemplate<bezier_t, point_list_t, t_point_t>(array) ;
}
bezier_t* wrapBezierConstructorBounds(const point_list_t& array, const real ub)
{
    return wrapBezierConstructorTemplate<bezier_t, point_list_t, t_point_t>(array, ub) ;
}
bezier_t* wrapBezierConstructorConstraints(const point_list_t& array, const curve_constraints_t& constraints)
{
    return wrapBezierConstructorConstraintsTemplate<bezier_t, point_list_t, t_point_t, curve_constraints_t>(array, constraints) ;
}
bezier_t* wrapBezierConstructorBoundsConstraints(const point_list_t& array, const curve_constraints_t& constraints, const real ub)
{
    return wrapBezierConstructorConstraintsTemplate<bezier_t, point_list_t, t_point_t, curve_constraints_t>(array, constraints, ub) ;
}
/*END 3D constructors */
/*6D constructors */
bezier6_t* wrapBezierConstructor6(const point_list6_t& array)
{
    return wrapBezierConstructorTemplate<bezier6_t, point_list6_t, t_point6_t>(array) ;
}
bezier6_t* wrapBezierConstructorBounds6(const point_list6_t& array, const real ub)
{
    return wrapBezierConstructorTemplate<bezier6_t, point_list6_t, t_point6_t>(array, ub) ;
}
bezier6_t* wrapBezierConstructor6Constraints(const point_list6_t& array, const curve_constraints6_t& constraints)
{
    return wrapBezierConstructorConstraintsTemplate<bezier6_t, point_list6_t, t_point6_t, curve_constraints6_t>(array, constraints) ;
}
bezier6_t* wrapBezierConstructorBounds6Constraints(const point_list6_t& array, const curve_constraints6_t& constraints, const real ub)
{
    return wrapBezierConstructorConstraintsTemplate<bezier6_t, point_list6_t, t_point6_t, curve_constraints6_t>(array, constraints, ub) ;
}
/*END 6D constructors */

polynom_t* wrapSplineConstructor(const coeff_t& array)
{
    return new polynom_t(array, 0., 1.);
}


t_waypoint_t getWayPoints(const coeff_t& array, const time_waypoints_t& time_wp)
{
    t_waypoint_t res;
    for(int i =0;i<array.cols();++i)
        res.push_back(std::make_pair(time_wp(i), array.col(i)));
    return res;
}

exact_cubic_t* wrapExactCubicConstructor(const coeff_t& array, const time_waypoints_t& time_wp)
{
    t_waypoint_t wps = getWayPoints(array, time_wp);
    return new exact_cubic_t(wps.begin(), wps.end());
}


spline_deriv_constraint_t* wrapSplineDerivConstraint(const coeff_t& array, const time_waypoints_t& time_wp, const curve_constraints_t& constraints)
{
    t_waypoint_t wps = getWayPoints(array, time_wp);
    return new spline_deriv_constraint_t(wps.begin(), wps.end(),constraints);
}

spline_deriv_constraint_t* wrapSplineDerivConstraintNoConstraints(const coeff_t& array, const time_waypoints_t& time_wp)
{
    t_waypoint_t wps = getWayPoints(array, time_wp);
    return new spline_deriv_constraint_t(wps.begin(), wps.end());
}

point_t get_init_vel(const curve_constraints_t& c)
{
    return c.init_vel;
}

point_t get_init_acc(const curve_constraints_t& c)
{
    return c.init_acc;
}

point_t get_end_vel(const curve_constraints_t& c)
{
    return c.end_vel;
}

point_t get_end_acc(const curve_constraints_t& c)
{
    return c.end_acc;
}

void set_init_vel(curve_constraints_t& c, const point_t& val)
{
    c.init_vel = val;
}

void set_init_acc(curve_constraints_t& c, const point_t& val)
{
    c.init_acc = val;
}

void set_end_vel(curve_constraints_t& c, const point_t& val)
{
    c.end_vel = val;
}

void set_end_acc(curve_constraints_t& c, const point_t& val)
{
    c.end_acc = val;
}

}
 // namespace spline
