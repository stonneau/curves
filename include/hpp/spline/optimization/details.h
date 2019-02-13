/**
* \file bezier_curve.h
* \brief class allowing to create a Bezier curve of dimension 1 <= n <= 3.
* \author Steve T.
* \version 0.1
* \date 06/17/2013
*/


#ifndef _CLASS_LINEAR_PROBLEM_DETAILS
#define _CLASS_LINEAR_PROBLEM_DETAILS

#include <hpp/spline/bezier_curve.h>
#include <hpp/spline/linear_variable.h>
#include <hpp/spline/curve_constraint.h>
#include <hpp/spline/optimization/definitions.h>

#include <Eigen/StdVector>

namespace spline
{
namespace  optimization
{
template<typename Point, int Dim, typename Numeric>
struct problem_data
{
     problem_data() : bezier(0){}
    ~problem_data() {if (bezier) delete bezier;}

    typedef linear_variable<Dim, Numeric>     var_t;
    typedef Eigen::Matrix<Numeric, Eigen::Dynamic, Dim > matrix_dim_x_t;
    typedef Eigen::Matrix<Numeric, Eigen::Dynamic, 1 >   vector_x_t;
    typedef bezier_curve<Numeric, Numeric, Dim, true,variables<linear_variable<Dim, Numeric> > > bezier_t;

    std::vector<var_t> variables_; // includes constant variables
    std::size_t numVariables; // total number of variable (* DIM for total size)
    std::size_t numControlPoints; // total number of variable (* DIM for total size)
    std::size_t startVariableIndex; //before that index, variables are constant
    std::size_t numStateConstraints;
    matrix_dim_x_t ineqMatrix;
    bezier_t* bezier;
};

inline std::size_t num_active_constraints(const constraint_flag& flag)
{
    long lValue = (long)(flag);
    std::size_t iCount = 0;
    while (lValue != 0)
    {
        lValue = lValue & (lValue - 1);
        iCount++;
    }
    return (flag & NONE) ? iCount-1 : iCount;
}

template < typename LinearVar, typename Variables>
Variables fillWithZeros(const LinearVar& var, const std::size_t totalvar, const std::size_t i)
{
    Variables res;
    std::vector<LinearVar>& vars = res.variables_;
    for (std::size_t idx = 0; idx < i; ++idx)
        vars.push_back(LinearVar::Zero());
    vars.push_back(var);
    for (std::size_t idx = i+1; idx < totalvar; ++idx)
        vars.push_back(LinearVar::Zero());
    return res;
}

template < typename Numeric, typename Bezier, typename LinearVar, typename Variables>
Bezier* computeLinearControlPoints(const std::vector<LinearVar>& linearVars, const Numeric totalTime )
{
    std::vector<Variables> res;
    // now need to fill all this with zeros...
    std::size_t totalvar = linearVars.size();
    for (std::size_t i = 0; i < totalvar; ++i)
        res.push_back( fillWithZeros<LinearVar, Variables>(linearVars[i],totalvar,i));
    return new Bezier(res.begin(),res.end(), totalTime);
}


template<typename Point, int Dim, typename Numeric>
problem_data<Point, Dim, Numeric> setup_control_points(const problem_definition<Point, Dim, Numeric>& pDef)
{
    typedef Numeric num_t;
    typedef Point   point_t;
    typedef linear_variable<Dim, Numeric>     var_t;
    typedef variables<var_t>    vars_t;
    typedef problem_data<Point, Dim, Numeric> problem_data_t;

    const curve_constraints<Point>& constraints = pDef.curveConstraints;
    const std::size_t& degree = pDef.degree;
    const constraint_flag& flag = pDef.flag;

    const std::size_t numControlPoints = pDef.degree +1;
    const std::size_t numActiveConstraints = num_active_constraints(flag);
    if (numActiveConstraints >= numControlPoints)
        throw std::runtime_error("In setup_control_points; too many constraints for the considered degree");


    problem_data_t problemData;
    typename vars_t::T_var_t& variables_ = problemData.variables_;

    std::size_t numConstants = 0;
    std::size_t i =0;
    if(flag & INIT_POS)
    {
        variables_.push_back(var_t(pDef.start));
        ++numConstants;
        ++i;
        if(flag & INIT_VEL)
        {
            point_t vel = pDef.start + constraints.init_vel / (num_t)degree;
            variables_.push_back(var_t(vel));
            ++numConstants;
            ++i;
            if(flag & INIT_ACC)
            {
                point_t acc = constraints.init_acc / (num_t)(degree * (degree-1))
                        + 2* vel- pDef.start;;
                variables_.push_back(var_t(acc));
                ++numConstants;
                ++i;
            }
        }
    }
    const std::size_t first_variable_idx = i;
    // variables
    for(; i + 3< numControlPoints; ++i)
        variables_.push_back(var_t());
    //end constraints
    if(flag & END_POS)
    {
        if(flag & END_VEL)
        {
            point_t vel = pDef.end - constraints.end_vel  / (num_t)degree;
            if(flag & END_ACC)
            {
                point_t acc = constraints.end_acc  / (num_t)(degree * (degree-1))
                        + 2* vel - pDef.end;
                variables_.push_back(var_t(acc));
                ++numConstants; ++i;
            }
            else if(i<numControlPoints-2)
            {
                variables_.push_back(var_t());
                ++i;
            }
            variables_.push_back(var_t(vel));
            ++numConstants; ++i;
        }
        else
        {
            while(i<numControlPoints-1)
            {
                variables_.push_back(var_t());
                ++i;
            }
        }
        variables_.push_back(var_t(pDef.end));
        ++numConstants; ++i;
    }
    // add remaining variables (only if no end_pos constraints)
    for(; i<numControlPoints; ++i)
        variables_.push_back(var_t());

    assert(numControlPoints > numConstants);
    assert(numControlPoints == variables_.size());


    problemData.bezier = computeLinearControlPoints<Numeric,
                                            bezier_curve<Numeric, Numeric, Dim, true,vars_t>,
                                            var_t, vars_t>(variables_,  pDef.totalTime);
    problemData.numControlPoints = numControlPoints;
    problemData.numVariables = numControlPoints-numConstants;
    problemData.startVariableIndex =first_variable_idx;
    problemData.numStateConstraints = numActiveConstraints - problemData.numVariables;
    return problemData;
}


// TODO assumes constant are inside constraints...
template<typename Point, int Dim, typename Numeric>
long compute_num_ineq_control_points
(const problem_definition<Point, Dim, Numeric>& pDef, const problem_data<Point, Dim, Numeric> & pData)
{
    typedef problem_definition<Point, Dim, Numeric> problem_definition_t;
    long rows(0);
    // rows depends on each constraint size, and the number of waypoints
    for (typename problem_definition_t::CIT_vectorx_t cit = pDef.inequalityVectors_.begin();
         cit != pDef.inequalityVectors_.end(); ++cit)
        rows += cit->rows() * pData.numControlPoints;
    return rows;
}

template<typename Point, int Dim, typename Numeric>
long compute_num_ineq_state_constraints
(const problem_definition<Point, Dim, Numeric>& pDef, const problem_data<Point, Dim, Numeric> & pData)
{
    //TODO
    return 0;
}


template< typename Point, int Dim, typename Numeric, typename Bezier, typename T_matrix_t, typename T_vector_t>
void bezierWaypointsToMatrixForm(const std::size_t startVariableIndex, const std::size_t numVariables,
                                 const Bezier& bezier, T_matrix_t& matrices, T_vector_t& vectors)
{
    typedef Eigen::Matrix<Numeric, Dim, Eigen::Dynamic> matrix_t;
    typedef Eigen::Matrix<Numeric, Dim, 1> vector_t;
    typedef linear_variable<Dim, Numeric> linear_variable_t;
    typedef variables<linear_variable_t> variables_t;
    typedef typename Bezier::t_point_t t_point;
    typedef typename Bezier::t_point_t::const_iterator cit_point;
    const t_point& wps = bezier.waypoints();
    // each control has a linear expression depending on all variables
    for(cit_point cit = wps.begin(); cit != wps.end(); ++cit)
    {
        matrix_t matCurrentWp = matrix_t::Zero(Dim, Dim*numVariables);
        vector_t vecCurrentWp = vector_t::Zero();
        const std::vector<linear_variable_t>& variables = cit->variables_;
        for(typename variables_t::CIT_var_t varit = variables.begin();
            varit != variables.end(); ++varit)
                vecCurrentWp +=  varit->b_;
        //assert(variables.begin() + startVariableIndex + numVariables <= variables.end());
        int col = 0;
        // loop only through variables that are not constant
        // TODO??? allow to put constant everywhere ?
        for(typename variables_t::CIT_var_t  varit = variables.begin() + startVariableIndex;
            varit != variables.begin() + startVariableIndex + numVariables; ++varit, col+=Dim)
        {
            matCurrentWp.block(0,col,Dim,Dim) =  varit->A_;
        }
        matrices.push_back(matCurrentWp);
        vectors.push_back(vecCurrentWp);
    }
}

template<typename Point, int Dim, typename Numeric>
std::vector<bezier_curve<Numeric, Numeric, Dim, true,
            variables<linear_variable<Dim, Numeric> > > >
split(const problem_definition<Point, Dim, Numeric>& pDef, problem_data<Point, Dim, Numeric> & pData)
{
    typedef Numeric real;
    typedef linear_variable<Dim, Numeric> linear_variable_t;
    typedef variables<linear_variable_t> variables_t;
    typedef bezier_curve< Numeric, Numeric, Dim, true,variables_t> bezier_t;
    typedef std::vector<bezier_t> T_bezier_t;

    const Eigen::VectorXd& times = pDef.splitTimes_;
    T_bezier_t res;
    bezier_t& current = *pData.bezier;
    real current_time = 0.;
    real tmp;
    for(int i = 0; i < times.rows(); ++i)
    {
        tmp =times[i];
        std::pair<bezier_t, bezier_t> pairsplit = current.split(tmp-current_time);
        res.push_back(pairsplit.first);
        current = pairsplit.second;
        current_time += tmp-current_time;
    }
    res.push_back(current);
    return res;
}


// TODO assumes constant are inside constraints...
template<typename Point, int Dim, typename Numeric>
Eigen::Matrix<Numeric, Eigen::Dynamic, Eigen::Dynamic> initInequalityMatrix
(const problem_definition<Point, Dim, Numeric>& pDef, problem_data<Point, Dim, Numeric> & pData)
{
    typedef problem_definition<Point, Dim, Numeric> problem_definition_t;
    typedef typename problem_definition_t::matrix_x_t matrix_x_t;
    typedef typename problem_definition_t::vectorx_t vectorx_t;
    typedef bezier_curve<Numeric, Numeric, Dim, true,
            variables<linear_variable<Dim, Numeric> > > bezier_t;
    typedef std::vector<bezier_t> T_bezier_t;
    typedef typename T_bezier_t::const_iterator CIT_bezier_t;
    typedef Eigen::Matrix<Numeric, Dim, Eigen::Dynamic> matrix_dimx_t;
    typedef Eigen::Matrix<Numeric, Dim, 1> vector_dim_t;
    typedef std::vector<matrix_dimx_t, Eigen::aligned_allocator<matrix_dimx_t> > T_matrix_dimx_t;
    typedef std::vector<matrix_dimx_t, Eigen::aligned_allocator<vector_dim_t> > T_vector_dim_t;
    typedef typename T_matrix_dimx_t::const_iterator CIT_matrix_dimx_t;
    typedef typename T_vector_dim_t::const_iterator CIT_vector_dim_t;

    long cols =  pData.numVariables * Dim;
    long rows = compute_num_ineq_control_points<Point, Dim, Numeric>(pDef, pData);
    //rows+= compute_num_ineq_state_constraints<Point, Dim, Numeric>(pDef, pData); // TODO
    matrix_x_t ineqMatrix = matrix_x_t::Zero(rows,cols);
    vectorx_t ineqVec = vectorx_t::Zero(rows);

    // compute sub-bezier curves
    T_bezier_t beziers = split<Point, Dim, Numeric>(pDef,pData);

    assert(pDef.inequalityMatrices_.size() == pDef.inequalityVectors_.size());
    assert(pDef.inequalityMatrices_.size() == beziers.size());

    long currentRowIdx = 0;
    typename problem_definition_t::CIT_matrix_dim_t cmit = pDef.inequalityMatrices_.begin();
    typename problem_definition_t::CIT_vectorx_t cvit = pDef.inequalityVectors_.begin();
    // for each bezier split ..
    for (CIT_bezier_t bit = beziers.begin();
         bit != beziers.end(); ++bit, ++cvit, ++cmit)
    {
        //compute vector of linear expressions of each control point
        T_matrix_dimx_t matrices; T_vector_dim_t vectors;
        bezierWaypointsToMatrixForm<Point, Dim, Numeric, bezier_t, T_matrix_dimx_t, T_vector_dim_t>
                (pData.startVariableIndex, pData.numVariables, *bit, matrices, vectors);
        // copy constraints the number of times there are control points
        CIT_vector_dim_t vdim_cit = vectors.begin();
        for(CIT_matrix_dimx_t mdim_cit = matrices.begin();
            mdim_cit != matrices.end(); ++mdim_cit, ++vdim_cit)
        {
            ineqMatrix.block(currentRowIdx, 0,cmit->rows(),cols)
                    = (*cmit)*(*mdim_cit) ; // constraint inequality for current bezier * expression of control point
            ineqVec.segment(currentRowIdx,cmit->rows()) = *cvit - (*cmit)*(*vdim_cit) ;
            currentRowIdx += cmit->rows();
        }
    }
    assert (rows == currentRowIdx); // we filled all the constraints
    return ineqMatrix;
}

inline constraint_flag operator~(constraint_flag a)
{return static_cast<constraint_flag>(~static_cast<const int>(a));}

inline constraint_flag operator|(constraint_flag a, constraint_flag b)
{return static_cast<constraint_flag>(static_cast<const int>(a) | static_cast<const int>(b));}

inline constraint_flag operator&(constraint_flag a, constraint_flag b)
{return static_cast<constraint_flag>(static_cast<const int>(a) & static_cast<const int>(b));}

inline constraint_flag operator^(constraint_flag a, constraint_flag b)
{return static_cast<constraint_flag>(static_cast<const int>(a) ^ static_cast<const int>(b));}

inline constraint_flag& operator|=(constraint_flag& a, constraint_flag b)
{return (constraint_flag&)((int&)(a) |= static_cast<const int>(b));}

inline constraint_flag& operator&=(constraint_flag& a, constraint_flag b)
{return (constraint_flag&)((int&)(a) &= static_cast<const int>(b));}

inline constraint_flag& operator^=(constraint_flag& a, constraint_flag b)
{return (constraint_flag&)((int&)(a) ^= static_cast<const int>(b));}

} // namespace optimization
} // namespace spline
#endif //_CLASS_LINEAR_PROBLEM_DETAILS

