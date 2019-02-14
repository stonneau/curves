#include "hpp/spline/bezier_curve.h"
#include "hpp/spline/linear_variable.h"

#include "python_definitions.h"
#include "hpp/spline/optimization/definitions.h"
#include "hpp/spline/optimization/linear_problem.h"

#include <vector>

#ifndef _VARIABLES_PYTHON_BINDINGS
#define _VARIABLES_PYTHON_BINDINGS


namespace spline
{
namespace optimization
{
static const int dim = 3;
typedef linear_variable<dim, real> linear_variable_3_t;
typedef variables<linear_variable_3_t> variables_3_t;
typedef bezier_curve  <real, real, dim, true, variables_3_t> bezier_linear_variable_t;
typedef problem_definition<point_t, dim, real> problem_definition_t;
typedef problem_data<point_t, dim, real>problem_data_t;
typedef problem<point_t, dim, real> problem_t;

problem_data_t setup_control_points_3_t(problem_definition_t &pDef);


/*linear variable control points*/
bezier_linear_variable_t* pDataBezier(const problem_data_t* pData);

bezier_linear_variable_t* wrapBezierLinearConstructor(const point_list_t& matrices, const point_list_t& vectors);

bezier_linear_variable_t* wrapBezierLinearConstructorBounds
    (const point_list_t& matrices, const point_list_t& vectors, const real ub);

typedef std::pair<Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic>,
                  Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> > linear_points_t;

struct MatrixVector
{
    linear_points_t res;
    Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> A() {return res.first;}
    Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> b() {return res.second;}
};

MatrixVector generate_problem_3_t(const problem_definition_t &pDef);

MatrixVector* wayPointsToLists(const bezier_linear_variable_t& self);

struct LinearBezierVector
{
    std::vector<bezier_linear_variable_t> beziers_;
    std::size_t size() {return beziers_.size();}
    bezier_linear_variable_t* at(std::size_t i)
    {
        assert (i<size());
        return new bezier_linear_variable_t(beziers_[i]);
    }
};

// does not include end time
LinearBezierVector* split_py(const bezier_linear_variable_t& self,  const vectorX_t& times);

void set_pd_flag(problem_definition_t* pDef, const int flag);
void set_start(problem_definition_t* pDef, const point_t& val );
void set_end(problem_definition_t* pDef, const point_t& val );
void set_degree(problem_definition_t* pDef, const std::size_t val );
void set_total_time(problem_definition_t* pDef, const std::size_t val );
void set_split_time(problem_definition_t* pDef, const Eigen::VectorXd& val );
Eigen::VectorXd get_split_times(const problem_definition_t* pDef);
constraint_flag get_pd_flag(const problem_definition_t* pDef);
Eigen::Vector3d get_start(const problem_definition_t* pDef);
Eigen::Vector3d get_end(const problem_definition_t* pDef);
std::size_t get_degree(const problem_definition_t* pDef);
double get_total_time(const problem_definition_t* pDef);
Eigen::VectorXd get_split_times(const problem_definition_t* pDef);
MatrixVector* get_ineq_at(const problem_definition_t* pDef, const std::size_t idx);
bool del_ineq_at(problem_definition_t* pDef, const std::size_t idx);
bool add_ineq_at(problem_definition_t* pDef, const Eigen::MatrixXd ineq, const Eigen::VectorXd vec);
} //namespace optimization
} //namespace spline.


#endif //_VARIABLES_PYTHON_BINDINGS
