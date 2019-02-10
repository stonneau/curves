/**
* \file bezier_curve.h
* \brief class allowing to create a Bezier curve of dimension 1 <= n <= 3.
* \author Steve T.
* \version 0.1
* \date 06/17/2013
*/


#ifndef _CLASS_LINEAR_PROBLEM
#define _CLASS_LINEAR_PROBLEM


#include <hpp/spline/bezier_curve.h>
#include <hpp/spline/linear_variable.h>
#include <hpp/spline/curve_constraint.h>

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
struct problem_data
{
    typedef linear_variable<Dim, Numeric>     var_t;
    typedef variables<var_t>    vars_t;

    vars_t variables_;
    std::size_t numVariables;
    std::size_t startVariableIndex;
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
    return iCount;
}

template<typename Point, int Dim, typename Numeric>
problem_data<Point, Dim, Numeric> setup_control_points(const std::size_t degree,
                          const constraint_flag flag,
                          const Point& initPos = Point(),
                          const Point& endPos  = Point(),
                          const curve_constraints<Point>& constraints = curve_constraints<Point>())
{
    typedef Numeric num_t;
    typedef Point   point_t;
    typedef linear_variable<Dim, Numeric>     var_t;
    typedef variables<var_t>    vars_t;
    typedef problem_data<Point, Dim, Numeric> problem_data_t;

    const std::size_t numControlPoints = degree +1;
    if (num_active_constraints(flag) >= numControlPoints)
        throw std::runtime_error("In setup_control_points; too many constraints for the considered degree");

    vars_t res;
    typename vars_t::T_var_t& variables_ = res.variables_;

    std::size_t numConstants = 0;
    std::size_t i =0;
    if(flag & INIT_POS)
    {
        variables_.push_back(var_t(initPos));
        ++numConstants;
        ++i;
        if(flag & INIT_VEL)
        {
            point_t vel = initPos + constraints.init_vel / (num_t)degree;
            variables_.push_back(var_t(vel));
            ++numConstants;
            ++i;
            if(flag & INIT_ACC)
            {
                point_t acc = constraints.init_acc / (num_t)(degree * (degree-1))
                        + 2* vel- initPos;;
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
            point_t vel = endPos - constraints.end_vel  / (num_t)degree;
            if(flag & END_ACC)
            {
                point_t acc = constraints.end_acc  / (num_t)(degree * (degree-1))
                        + 2* vel - endPos;
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
        variables_.push_back(var_t(endPos));
        ++numConstants; ++i;
    }
    // add remaining variables (only if no end_pos constraints)
    for(; i<numControlPoints; ++i)
        variables_.push_back(var_t());

    assert(numControlPoints > numConstants);
    assert(numControlPoints == variables_.size());


    problem_data_t problemData;
    problemData.variables_ = res;
    problemData.numVariables = numControlPoints-numConstants;
    problemData.startVariableIndex =first_variable_idx;
    return problemData;
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
#endif //_CLASS_LINEAR_PROBLEM

