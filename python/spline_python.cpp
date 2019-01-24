#include "spline/bezier_curve.h"
#include "spline/polynom.h"
#include "spline/exact_cubic.h"
#include "spline/spline_deriv_constraint.h"
#include "spline/curve_constraint.h"
#include "spline/bezier_polynom_conversion.h"
#include "spline/bernstein.h"


#include <vector>

#include <eigenpy/memory.hpp>
#include <eigenpy/eigenpy.hpp>

#include <boost/python.hpp>

/*** TEMPLATE SPECIALIZATION FOR PYTHON ****/
typedef double real;
typedef Eigen::Matrix3d matrix3_t;
typedef Eigen::Matrix3d matrixX_t;
typedef Eigen::Vector3d point_t;
typedef Eigen::VectorXd vectorX_t;
typedef Eigen::Matrix<double, 6, 1, 0, 6, 1> point6_t;
typedef Eigen::Matrix<double, 3, 1, 0, 3, 1> ret_point_t;
typedef Eigen::Matrix<double, 6, 1, 0, 6, 1> ret_point6_t;
typedef Eigen::VectorXd time_waypoints_t;
typedef Eigen::Matrix<real, 3, Eigen::Dynamic> point_list_t;
typedef Eigen::Matrix<real, 6, Eigen::Dynamic> point_list6_t;
typedef std::vector<point_t,Eigen::aligned_allocator<matrix3_t> >  t_matrix3_t;
typedef std::vector<point_t,Eigen::aligned_allocator<matrixX_t> >  t_matrixX_t;
typedef std::vector<point_t,Eigen::aligned_allocator<point_t> >  t_point_t;
typedef std::vector<point6_t,Eigen::aligned_allocator<point6_t> >  t_point6_t;
typedef std::pair<real, point_t> Waypoint;
typedef std::vector<Waypoint> T_Waypoint;
typedef std::pair<real, point6_t> Waypoint6;
typedef std::vector<Waypoint6> T_Waypoint6;

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

EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(bernstein_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(bezier_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(bezier6_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(polynom_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(exact_cubic_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(curve_constraints_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(spline_deriv_constraint_t)

namespace spline
{
using namespace boost::python;
template <typename PointList, typename T_Point>
T_Point vectorFromEigenArray(const PointList& array)
{
    T_Point res;
    for(int i =0;i<array.cols();++i)
        res.push_back(array.col(i));
    return res;
}

// y_i = A x_i + 3, dimension 3 for each y_i
struct linear_control_t{
    matrix3_t A_;
    point_t b_;

    linear_control_t(): A_(matrix3_t()), b_(point_t()){}
    linear_control_t(const matrix3_t& A, const point_t& b):A_(A),b_(b) {}


    linear_control_t& operator+=(const linear_control_t& w1)
    {
        this->A_ += w1.A_;
        this->b_ += w1.b_;
        return *this;
    }
    linear_control_t& operator-=(const linear_control_t& w1)
    {
        this->A_ -= w1.A_;
        this->b_ -= w1.b_;
        return *this;
    }

    static linear_control_t Zero(size_t dim=0){
        linear_control_t w;
        w.A_  = matrix3_t::Zero();
        w.b_  = point_t::Zero();
        return w;
    }
};

linear_control_t operator+(const linear_control_t& w1, const linear_control_t& w2)
{
    return linear_control_t(w1.A_ + w2.A_, w1.b_ + w2.b_);
}
linear_control_t operator-(const linear_control_t& w1, const linear_control_t& w2)
{
    return linear_control_t(w1.A_ - w2.A_, w1.b_ - w2.b_);
}

linear_control_t operator*(const double k, const linear_control_t& w){
    return linear_control_t(k*w.A_,k*w.b_);
}

linear_control_t operator*(const linear_control_t& w,const double k){
    return linear_control_t(k*w.A_,k*w.b_);
}

linear_control_t operator/(const linear_control_t& w,const double k){
    return linear_control_t(w.A_/k,w.b_/k);
}

struct t_linear_control_t{
    std::vector<linear_control_t> variables_;

    t_linear_control_t() {}


    t_linear_control_t& operator+=(const t_linear_control_t& w1)
    {
        if(variables_.size() == 0)
            variables_ = w1.variables_;
        else if (w1.variables_.size() !=0)
        {
            assert(variables_.size() == w1.variables_.size());
            std::vector<linear_control_t>::const_iterator cit = w1.variables_.begin();
            for(std::vector<linear_control_t>::iterator it = variables_.begin(); it != variables_.end(); ++it, ++cit)
            {
                (*it)+=(*cit);
            }
        }
        return *this;
    }

    t_linear_control_t& operator-=(const t_linear_control_t& w1)
    {
        if(variables_.size() == 0)
            variables_ = w1.variables_;
        else if (w1.variables_.size() !=0)
        {
            assert(variables_.size() == w1.variables_.size());
            std::vector<linear_control_t>::const_iterator cit = w1.variables_.begin();
            for(std::vector<linear_control_t>::iterator it = variables_.begin(); it != variables_.end(); ++it, ++cit)
            {
                (*it)-=(*cit);
            }
        }
        return *this;
    }

    static t_linear_control_t Zero(size_t /*dim*/){
        t_linear_control_t w;
        return w;
    }
};

t_linear_control_t operator+(const t_linear_control_t& w1, const t_linear_control_t& w2)
{
    if(w2.variables_.size() == 0)
        return w1;
    if(w1.variables_.size() == 0)
        return w2;
    t_linear_control_t res;
    assert(w2.variables_.size() == w1.variables_.size());
    std::vector<linear_control_t>::const_iterator cit = w1.variables_.begin();
    for(std::vector<linear_control_t>::const_iterator cit2 = w2.variables_.begin(); cit2 != w2.variables_.end(); ++cit, ++cit2)
    {
        res.variables_.push_back((*cit)+(*cit2));
    }
    return res;
}

t_linear_control_t operator-(const t_linear_control_t& w1, const t_linear_control_t& w2)
{
    if(w2.variables_.size() == 0)
        return w1;
    if(w1.variables_.size() == 0)
        return w2;
    t_linear_control_t res;
    assert(w2.variables_.size() == w1.variables_.size());
    std::vector<linear_control_t>::const_iterator cit = w1.variables_.begin();
    for(std::vector<linear_control_t>::const_iterator cit2 = w2.variables_.begin(); cit2 != w2.variables_.end(); ++cit, ++cit2)
    {
        res.variables_.push_back((*cit)-(*cit2));
    }
    return res;
}

t_linear_control_t operator*(const double k, const t_linear_control_t& w)
{
    if(w.variables_.size() == 0)
        return w;
    t_linear_control_t res;
    for(std::vector<linear_control_t>::const_iterator cit = w.variables_.begin(); cit != w.variables_.end(); ++cit)
    {
        res.variables_.push_back(k*(*cit));
    }
    return res;
}

t_linear_control_t operator*(const t_linear_control_t& w,const double k)
{
    if(w.variables_.size() == 0)
        return w;
    t_linear_control_t res;
    for(std::vector<linear_control_t>::const_iterator cit = w.variables_.begin(); cit != w.variables_.end(); ++cit)
    {
        res.variables_.push_back((*cit)*k);
    }
    return res;
}

t_linear_control_t operator/(const t_linear_control_t& w,const double k)
{
    if(w.variables_.size() == 0)
        return w;
    t_linear_control_t res;
    for(std::vector<linear_control_t>::const_iterator cit = w.variables_.begin(); cit != w.variables_.end(); ++cit)
    {
        res.variables_.push_back((*cit)/k);
    }
    return res;
}


//typedef spline::curve_constraints<linear_control_t> curve_constraints_linear_t;
typedef spline::bezier_curve  <real, real, 3, true, t_linear_control_t> bezier_linear_control_t;



std::vector<linear_control_t> matrix3DFromEigenArray(const point_list_t& matrices, const point_list_t& vectors)
{
    assert(vectors.cols() * 3  == matrices.cols() ) ;
    std::vector<linear_control_t> res;
    for(int i =0;i<vectors.cols();++i)
        res.push_back(linear_control_t(matrices.block<3,3>(0,i*3), vectors.col(i)));
    return res;
}

t_linear_control_t fillWithZeros(const linear_control_t& var, const std::size_t totalvar, const std::size_t i)
{
    t_linear_control_t res;
    std::vector<linear_control_t>& vars = res.variables_;
    for (std::size_t idx = 0; idx < i; ++idx)
        vars.push_back(linear_control_t::Zero());
    vars.push_back(var);
    for (std::size_t idx = i+1; idx < totalvar; ++idx)
        vars.push_back(linear_control_t::Zero());
    return res;
}

std::vector<t_linear_control_t> computeLinearControlPoints(const point_list_t& matrices, const point_list_t& vectors)
{
    std::vector<t_linear_control_t> res;
    std::vector<linear_control_t> variables = matrix3DFromEigenArray(matrices, vectors);
    // now need to fill all this with zeros...
    std::size_t totalvar = variables.size();
    for (std::size_t i = 0; i < totalvar; ++i)
        res.push_back( fillWithZeros(variables[i],totalvar,i));
    return res;
}


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

/*linear variable control points*/
bezier_linear_control_t* wrapBezierLinearConstructor(const point_list_t& matrices, const point_list_t& vectors)
{
    std::vector<t_linear_control_t> asVector = computeLinearControlPoints(matrices, vectors);
    return new bezier_linear_control_t(asVector.begin(), asVector.end(), 1.) ;
}
bezier_linear_control_t* wrapBezierLinearConstructorBounds(const point_list_t& matrices, const point_list_t& vectors, const real ub)
{
    std::vector<t_linear_control_t> asVector = computeLinearControlPoints(matrices, vectors);
    return new bezier_linear_control_t(asVector.begin(), asVector.end(), ub) ;
}
/*bezier_linear_control_t* wrapBezierLinearConstructorConstraints(const point_list_t& matrices, const point_list_t& vectors, const curve_constraints_linear_t& constraints)
{
    std::vector<linear_control_t> asVector = matrix3DFromEigenArray(matrices, vectors);
    return new bezier_linear_control_t(asVector.begin(), asVector.end(), constraints, 1.) ;
}
bezier_linear_control_t* wrapBezierLinearConstructorBoundsConstraints(const point_list_t& matrices, const point_list_t& vectors, const curve_constraints_linear_t& constraints, const real ub)
{
    std::vector<linear_control_t> asVector = matrix3DFromEigenArray(matrices, vectors);
    return new bezier_linear_control_t(asVector.begin(), asVector.end(), constraints, ub) ;
}*/
/*END 3D constructors */

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

typedef std::pair<Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> > linear_points_t;

struct LinearControlPointsHolder
{
    linear_points_t res;
    Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> A()
    {
        return res.first;
    }
    Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> b()
    {
        return res.second;
    }
};

LinearControlPointsHolder*
        wayPointsToLists(const bezier_linear_control_t& self)
{
    typedef typename bezier_linear_control_t::t_point_t t_point;
    typedef typename bezier_linear_control_t::t_point_t::const_iterator cit_point;
    const t_point& wps = self.waypoints();
    // retrieve num variables.
    std::size_t dim = wps[0].variables_.size()*3;
    Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> matrices (dim,wps.size() * 3);
    Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> vectors  (dim,wps.size());
    int col = 0;
    for(cit_point cit = wps.begin(); cit != wps.end(); ++cit, ++col)
    {
        const std::vector<linear_control_t>& variables = cit->variables_;
        int i = 0;
        for(std::vector<linear_control_t>::const_iterator varit = variables.begin();
            varit != variables.end(); ++varit, i+=3)
        {
            vectors.block<3,1>(i,col)   =  varit->b_;
            matrices.block<3,3>(i,col*3) = varit->A_;
        }
    }
    LinearControlPointsHolder* res (new LinearControlPointsHolder);
    res->res = std::make_pair(matrices, vectors);
    return res;
}

struct LinearBezierVector
{
    std::vector<bezier_linear_control_t> beziers_;
    std::size_t size()
    {
        return beziers_.size();
    }
    bezier_linear_control_t* at(std::size_t i)
    {
        assert (i<size());
        return new bezier_linear_control_t(beziers_[i]);
    }
};

// does not include end time
LinearBezierVector* split(const bezier_linear_control_t& self,  const vectorX_t& times)
{
    LinearBezierVector* res (new LinearBezierVector);
    bezier_linear_control_t current = self;
    real current_time = 0.;
    real tmp;
    for(int i = 0; i < times.rows(); ++i)
    {
        tmp =times[i];
        std::pair<bezier_linear_control_t, bezier_linear_control_t> pairsplit = current.split(tmp-current_time);
        res->beziers_.push_back(pairsplit.first);
        current = pairsplit.second;
        current_time += tmp-current_time;
    }
    res->beziers_.push_back(current);
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



BOOST_PYTHON_MODULE(spline)
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
    eigenpy::enableEigenPySpecific<vectorX_t,vectorX_t>();
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

    /** BEGIN bezier curve**/

    class_<LinearControlPointsHolder>
        ("LinearWaypoint", no_init)
        .def_readonly("A", &LinearControlPointsHolder::A)
        .def_readonly("b", &LinearControlPointsHolder::b)
        ;

    class_<LinearBezierVector>
        ("bezierVarVector", no_init)
        .def_readonly("size", &LinearBezierVector::size)
        .def("at", &LinearBezierVector::at, return_value_policy<manage_new_object>())
        ;

    class_<bezier_linear_control_t>
        ("bezierVar", no_init)
            .def("__init__", make_constructor(&wrapBezierLinearConstructor))
            .def("__init__", make_constructor(&wrapBezierLinearConstructorBounds))
            //.def("__init__", make_constructor(&wrapBezierConstructorConstraints))
            //.def("__init__", make_constructor(&wrapBezierConstructorBoundsConstraints))
            .def("min", &bezier_linear_control_t::min)
            .def("max", &bezier_linear_control_t::max)
            //.def("__call__", &bezier_linear_control_t::operator())
            .def("derivate", &bezier_linear_control_t::derivate)
            .def("compute_derivate", &bezier_linear_control_t::compute_derivate)
            .def("compute_primitive", &bezier_linear_control_t::compute_primitive)
            .def("split", &split, return_value_policy<manage_new_object>())
            .def("waypoints", &wayPointsToLists, return_value_policy<manage_new_object>())
            .def_readonly("degree", &bezier_linear_control_t::degree_)
            .def_readonly("nbWaypoints", &bezier_linear_control_t::size_)
        ;
    /** END bezier curve**/


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


}

} // namespace spline
