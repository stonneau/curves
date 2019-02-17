#include "python_variables.h"
#include "python_definitions.h"

#include <Eigen/Core>


namespace spline
{
namespace optimization {




void set_constraint(problem_definition_t &pDef, const curve_constraints_t& constraints)
{
    pDef.curveConstraints = constraints;
}
curve_constraints_t get_constraint(problem_definition_t &pDef)
{
    return pDef.curveConstraints;
}

problem_data_t setup_control_points_3_t(problem_definition_t &pDef)
{
    problem_data_t pData = setup_control_points<point_t,dim,real>(pDef);
    return pData;//return new problem_data_t(pData);
}


cost_function_t problem_t_cost(const problem_t& p)
{
    return p.cost;
}
Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> problem_t_ineqMatrix(const problem_t& p)
{
    return p.ineqMatrix;
}
Eigen::Matrix<real, Eigen::Dynamic, 1> problem_t_ineqVector(const problem_t& p)
{
    return p.ineqVector;
}

Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> cost_t_quad(const cost_function_t& p)
{
    return p.quadratic;
}
Eigen::Matrix<real, Eigen::Dynamic, 1> cost_t_linear(const cost_function_t & p)
{
    return p.linear;
}

problem_t generate_problem_3_t(const problem_definition_t &pDef)
{
    return generate_problem<point_t,dim,real>(pDef);
}

void set_pd_flag(problem_definition_t* pDef, const int flag)
{
    pDef->flag = (constraint_flag)(flag);
}
void set_start(problem_definition_t* pDef, const point_t &val )
{
    pDef->start = val;
}
void set_end(problem_definition_t* pDef, const point_t &val )
{
    pDef->end = val;
}
void set_degree(problem_definition_t* pDef, const std::size_t val )
{
    pDef->degree = val;
}
void set_total_time(problem_definition_t* pDef, const std::size_t val )
{
    pDef->totalTime = val;
}
void set_split_time(problem_definition_t* pDef, const Eigen::VectorXd& val )
{
    pDef->splitTimes_ = val;
}
Eigen::VectorXd get_split_times(const problem_definition_t* pDef)
{
    return pDef->splitTimes_;
}

constraint_flag get_pd_flag(const problem_definition_t* pDef)
{
    return pDef->flag;
}
Eigen::Vector3d get_start(const problem_definition_t* pDef)
{
    return pDef->start;
}
Eigen::Vector3d get_end(const problem_definition_t* pDef)
{
    return pDef->end;
}
std::size_t get_degree(const problem_definition_t* pDef)
{
    return pDef->degree;
}
double get_total_time(const problem_definition_t* pDef)
{
    return pDef->totalTime;
}

MatrixVector* get_ineq_at(const problem_definition_t* pDef, const std::size_t idx)
{
    if (idx > pDef->inequalityMatrices_.size() - 1)
        throw std::runtime_error("required id is beyond number of inequality matrices");
    MatrixVector* res = new MatrixVector();
    res->res = std::make_pair(pDef->inequalityMatrices_[idx], pDef->inequalityVectors_[idx]);
    return res;
}
bool del_ineq_at(problem_definition_t* pDef, const std::size_t idx)
{
    if (idx > pDef->inequalityMatrices_.size() - 1)
        return false;
    pDef->inequalityMatrices_.erase(pDef->inequalityMatrices_.begin() + idx -1);
    pDef->inequalityVectors_.erase(pDef->inequalityVectors_.begin() + idx -1);
    return true;
}
bool add_ineq_at(problem_definition_t* pDef, const Eigen::MatrixXd ineq, const Eigen::VectorXd vec)
{
    if (ineq.rows() != vec.rows())
        throw std::runtime_error("ineq vector and matrix do not have the same number of rows");
    if (!(pDef->inequalityMatrices_.empty()) && ineq.cols() != pDef->inequalityMatrices_.back().cols())
        throw std::runtime_error("inequality matrix does not have the same variable dimension as existing matrices");
    pDef->inequalityMatrices_.push_back(ineq);
    pDef->inequalityVectors_.push_back(vec);
    return true;
}

bezier_linear_variable_t* pDataBezier(const problem_data_t* pData)
{
    const bezier_linear_variable_t& b = *pData->bezier;
    return new bezier_linear_variable_t(b.waypoints().begin(), b.waypoints().end(),b.T_, b.mult_T_);
}

std::vector<linear_variable_3_t> matrix3DFromEigenArray(const point_list_t& matrices, const point_list_t& vectors)
{
    assert(vectors.cols() * 3  == matrices.cols() ) ;
    std::vector<linear_variable_3_t> res;
    for(int i =0;i<vectors.cols();++i)
        res.push_back(linear_variable_3_t(matrices.block<3,3>(0,i*3), vectors.col(i)));
    return res;
}

linear_variable_3_t fillWithZeros(const linear_variable_3_t& var, const std::size_t totalvar, const std::size_t i)
{
    linear_variable_3_t::matrix_dim_x_t B(linear_variable_3_t::matrix_dim_x_t::Zero(dim,totalvar*dim));
    B.block(0,dim*i,dim,dim) = var.B();
    return linear_variable_3_t (B,var.c());
}

std::vector<linear_variable_3_t> computeLinearControlPoints(const point_list_t& matrices, const point_list_t& vectors)
{
    std::vector<linear_variable_3_t> res;
    std::vector<linear_variable_3_t> variables = matrix3DFromEigenArray(matrices, vectors);
    // now need to fill all this with zeros...
    std::size_t totalvar = variables.size();
    for (std::size_t i = 0; i < totalvar; ++i)
        res.push_back( fillWithZeros(variables[i],totalvar,i));
    return res;
}

/*linear variable control points*/
bezier_linear_variable_t* wrapBezierLinearConstructor(const point_list_t& matrices, const point_list_t& vectors)
{
    std::vector<linear_variable_3_t> asVector = computeLinearControlPoints(matrices, vectors);
    return new bezier_linear_variable_t(asVector.begin(), asVector.end(), 1.) ;
}

bezier_linear_variable_t* wrapBezierLinearConstructorBounds(const point_list_t& matrices, const point_list_t& vectors, const real ub)
{
    std::vector<linear_variable_3_t> asVector = computeLinearControlPoints(matrices, vectors);
    return new bezier_linear_variable_t(asVector.begin(), asVector.end(), ub) ;
}


MatrixVector*
        wayPointsToLists(const bezier_linear_variable_t& self)
{
    typedef typename bezier_linear_variable_t::t_point_t t_point;
    typedef typename bezier_linear_variable_t::t_point_t::const_iterator cit_point;
    const t_point& wps = self.waypoints();
    // retrieve num variables.
    std::size_t dim = wps[0].B().cols();
    Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> matrices (dim,wps.size() * 3);
    Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> vectors =
            Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic>::Zero(3,wps.size());
    int i = 0;
    for(cit_point cit = wps.begin(); cit != wps.end(); ++cit, ++i)
    {
        matrices.block(0,i*3,dim,3) = cit->B().transpose();
        vectors.block<3,1>(0,i)   =  cit->c();
    }
    MatrixVector* res (new MatrixVector);
    res->res = std::make_pair(matrices, vectors);
    return res;
}


// does not include end time
LinearBezierVector* split_py(const bezier_linear_variable_t& self,  const vectorX_t& times)
{
    LinearBezierVector* res (new LinearBezierVector);
    bezier_linear_variable_t current = self;
    real current_time = 0.;
    real tmp;
    for(int i = 0; i < times.rows(); ++i)
    {
        tmp =times[i];
        std::pair<bezier_linear_variable_t, bezier_linear_variable_t> pairsplit = current.split(tmp-current_time);
        res->beziers_.push_back(pairsplit.first);
        current = pairsplit.second;
        current_time += tmp-current_time;
    }
    res->beziers_.push_back(current);
    return res;
}

} // namespace optimization
} // namespace spline
