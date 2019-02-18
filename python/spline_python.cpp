#include "hpp/spline/optimization/linear_problem.h"

#include "python_definitions.h"
#include "python_bezier_def.h"
#include "python_variables.h"

#include <vector>

#include <eigenpy/memory.hpp>
#include <eigenpy/eigenpy.hpp>

#include <boost/python.hpp>

using namespace  spline;
using namespace  spline::optimization;
using namespace boost::python;

/*** TEMPLATE SPECIALIZATION FOR PYTHON ****/
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(bernstein_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(bezier_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(bezier6_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(polynom_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(exact_cubic_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(curve_constraints_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(spline_deriv_constraint_t)

BOOST_PYTHON_MODULE(hpp_spline)
{
    /** BEGIN eigenpy init**/
    eigenpy::enableEigenPy();

    eigenpy::enableEigenPySpecific<point_t,point_t>();
    eigenpy::enableEigenPySpecific<ret_point_t,ret_point_t>();
    eigenpy::enableEigenPySpecific<point_list_t,point_list_t>();
    eigenpy::enableEigenPySpecific<point6_t,point6_t>();
    eigenpy::enableEigenPySpecific<ret_point6_t,ret_point6_t>();
    eigenpy::enableEigenPySpecific<point_list6_t,point_list6_t>();
    eigenpy::enableEigenPySpecific<coeff_t,coeff_t>();
    /*eigenpy::exposeAngleAxis();
    eigenpy::exposeQuaternion();*/
    /** END eigenpy init**/

    /** BEGIN bezier curve 6**/
    class_<bezier6_t>
        ("bezier6", no_init)
            .def("__init__", make_constructor(&wrapBezierConstructor6))
            .def("__init__", make_constructor(&wrapBezierConstructorBounds6))
            //.def("__init__", make_constructor(&wrapBezierConstructor6Constraints))
            //.def("__init__", make_constructor(&wrapBezierConstructorBounds6Constraints))
            .def("min", &bezier6_t::min)
            .def("max", &bezier6_t::max)
            .def("__call__", &bezier6_t::operator())
            .def("derivate", &bezier6_t::derivate)
            .def("compute_derivate", &bezier6_t::compute_derivate)
            .def("compute_primitive", &bezier6_t::compute_primitive)
            .def("waypoints", &wayPointsToList<bezier6_t,6>)
            .def_readonly("degree", &bezier6_t::degree_)
            .def_readonly("nbWaypoints", &bezier6_t::size_)
        ;
    /** END bezier curve**/

    /** BEGIN bezier curve**/
    class_<bezier_t>
        ("bezier", no_init)
            .def("__init__", make_constructor(&wrapBezierConstructor))
            .def("__init__", make_constructor(&wrapBezierConstructorBounds))
            .def("__init__", make_constructor(&wrapBezierConstructorConstraints))
            .def("__init__", make_constructor(&wrapBezierConstructorBoundsConstraints))
            .def("min", &bezier_t::min)
            .def("max", &bezier_t::max)
            .def("__call__", &bezier_t::operator())
            .def("derivate", &bezier_t::derivate)
            .def("compute_derivate", &bezier_t::compute_derivate)
            .def("compute_primitive", &bezier_t::compute_primitive)
            .def("waypoints", &wayPointsToList<bezier_t,3>)
            .def_readonly("degree", &bezier_t::degree_)
            .def_readonly("nbWaypoints", &bezier_t::size_)
        ;
    /** END bezier curve**/


    /** BEGIN variable points bezier curve**/
    class_<MatrixVector>
        ("MatrixVector", no_init)
        .def_readonly("A", &MatrixVector::A)
        .def_readonly("b", &MatrixVector::b)
        ;

    class_<LinearBezierVector>
        ("bezierVarVector", no_init)
        .def_readonly("size", &LinearBezierVector::size)
        .def("at", &LinearBezierVector::at, return_value_policy<manage_new_object>())
        ;

    class_<bezier_linear_variable_t>
        ("bezierVar", no_init)
            .def("__init__", make_constructor(&wrapBezierLinearConstructor))
            .def("__init__", make_constructor(&wrapBezierLinearConstructorBounds))
            .def("min", &bezier_linear_variable_t::min)
            .def("max", &bezier_linear_variable_t::max)
            //.def("__call__", &bezier_linear_control_t::operator())
            .def("derivate", &bezier_linear_variable_t::derivate)
            .def("compute_derivate", &bezier_linear_variable_t::compute_derivate)
            .def("compute_primitive", &bezier_linear_variable_t::compute_primitive)
            .def("split", &split_py, return_value_policy<manage_new_object>())
            .def("waypoints", &wayPointsToLists, return_value_policy<manage_new_object>())
            .def_readonly("degree", &bezier_linear_variable_t::degree_)
            .def_readonly("nbWaypoints", &bezier_linear_variable_t::size_)
        ;


    class_<LinearBezierVector>
        ("bezierVarVector", no_init)
        .def_readonly("size", &LinearBezierVector::size)
        .def("at", &LinearBezierVector::at, return_value_policy<manage_new_object>())
        ;

    class_<problem_definition_t>
        ("problemDefinition", init<>())
            .add_property("flag", &get_pd_flag, &set_pd_flag)
            .add_property("start", &get_start, &set_start)
            .add_property("end", &get_end, &set_end)
            .add_property("degree", &get_degree, &set_degree)
            .add_property("totalTime", &get_total_time, &set_total_time)
            .add_property("splits", &get_split_times, &set_split_time)
            .add_property("curveConstraints", &get_constraint, &set_constraint)
            .def("inequality", &get_ineq_at, return_value_policy<manage_new_object>())
            .def("removeInequality", &del_ineq_at)
            .def("addInequality", &add_ineq_at)
        ;


    class_<quadratic_variable_t >
        ("cost", no_init)
        .add_property("A", &cost_t_quad)
        .add_property("b", &cost_t_linear)
        .add_property("c", &cost_t_constant)
        ;

    class_<problem_t>
        ("problem", init<>())
        .add_property("cost", &problem_t_cost)
        .add_property("A", &problem_t_ineqMatrix)
        .add_property("b", &problem_t_ineqVector)
        ;


    def("setupControlPoints", &setup_control_points_3_t);
    def("generate_problem", &generate_problem_3_t);

    class_<problem_data_t>
        ("problemData", no_init)
            .def("bezier", &pDataBezier, return_value_policy<manage_new_object>())
            .def_readonly("numControlPoints", &problem_data_t::numControlPoints)
            .def_readonly("numVariables", &problem_data_t::numVariables)
            .def_readonly("startVariableIndex", &problem_data_t::startVariableIndex)
            .def_readonly("numStateConstraints", &problem_data_t::numStateConstraints)
        ;


    enum_<constraint_flag>("constraint_flag")
            .value("INIT_POS", INIT_POS)
            .value("INIT_VEL", INIT_VEL)
            .value("INIT_ACC", INIT_ACC)
            .value("END_POS", END_POS)
            .value("END_VEL", END_VEL)
            .value("END_ACC", END_ACC)
            .value("ALL", ALL)
            .value("NONE", NONE)
            .export_values();
    /** END variable points bezier curve**/

    /** BEGIN spline curve function**/
    class_<polynom_t>("polynom",  init<const polynom_t::coeff_t, const real, const real >())
            .def("__init__", make_constructor(&wrapSplineConstructor))
            .def("min", &polynom_t::min)
            .def("max", &polynom_t::max)
            .def("__call__", &polynom_t::operator())
            .def("derivate", &polynom_t::derivate)
        ;
    /** END cubic function**/


    /** BEGIN exact_cubic curve**/
    class_<exact_cubic_t>
        ("exact_cubic", no_init)
            .def("__init__", make_constructor(&wrapExactCubicConstructor))
            .def("min", &exact_cubic_t::min)
            .def("max", &exact_cubic_t::max)
            .def("__call__", &exact_cubic_t::operator())
            .def("derivate", &exact_cubic_t::derivate)
        ;
    /** END bezier curve**/


    /** BEGIN curve constraints**/
    class_<curve_constraints_t>
        ("curve_constraints", init<>())
            .add_property("init_vel", &get_init_vel, &set_init_vel)
            .add_property("init_acc", &get_init_acc, &set_init_acc)
            .add_property("end_vel", &get_end_vel, &set_end_vel)
            .add_property("end_acc", &get_end_acc, &set_end_acc)
        ;
    /** END curve constraints**/


    /** BEGIN spline_deriv_constraints**/
    class_<spline_deriv_constraint_t>
        ("spline_deriv_constraint", no_init)
            .def("__init__", make_constructor(&wrapSplineDerivConstraint))
            .def("__init__", make_constructor(&wrapSplineDerivConstraintNoConstraints))
            .def("min", &exact_cubic_t::min)
            .def("max", &exact_cubic_t::max)
            .def("__call__", &exact_cubic_t::operator())
            .def("derivate", &exact_cubic_t::derivate)
        ;
    /** END spline_deriv_constraints**/

    /** BEGIN bernstein polynom**/
    class_<bernstein_t>
        ("bernstein", init<const unsigned int, const unsigned int>())
            .def("__call__", &bernstein_t::operator())
        ;
    /** END bernstein polynom**/

    /** BEGIN Bezier to polynom conversion**/
    def("from_bezier", from_bezier<bezier_t,polynom_t>);
    /** END Bezier to polynom conversion**/


} // namespace spline
