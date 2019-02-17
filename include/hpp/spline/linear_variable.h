/**
* \file linear_variable.h
* \brief storage for variable points of the form p_i = a_i x + b_i
* \author Steve T.
* \version 0.1
* \date 07/02/2019
*/


#ifndef _CLASS_LINEAR_VARIABLE
#define _CLASS_LINEAR_VARIABLE

#include "curve_abc.h"

#include "MathDefs.h"

#include <math.h>
#include <vector>
#include <Eigen/Core>
#include <stdexcept>

namespace spline
{
template <int Dim, typename Numeric=double>
struct linear_variable
{
    typedef Eigen::Matrix<Numeric, Dim, Dim> matrix_dim_t;
    typedef Eigen::Matrix<Numeric, Dim, Eigen::Dynamic> matrix_dim_x_t;
    typedef Eigen::Matrix<Numeric, Dim, 1> point_dim_t;
    typedef Eigen::Matrix<Numeric, Eigen::Dynamic, 1> vectord_t;
    typedef linear_variable<Dim, Numeric> linear_variable_t;

    linear_variable(): B_(matrix_dim_t::Identity()), c_(point_dim_t::Zero()){} //variable
    linear_variable(const point_dim_t& c):B_(matrix_dim_t::Zero()),c_(c) {} // constant
    linear_variable(const matrix_dim_x_t& B, const point_dim_t& c):B_(B),c_(c) {} //mixed

    // linear evaluation
    point_dim_t operator()(const Eigen::Ref<const point_dim_t>& val) const
    {
        return B() * val + c();
    }

    linear_variable_t& operator+=(const linear_variable_t& w1)
    {
        // handling zero case
        if(c_.rows() == 0)
        {
            this->B_ = w1.B_;
            this->c_ = w1.c_;
        }
        else
        {
            this->B_ += w1.B_;
            this->c_ += w1.c_;
        }
        return *this;
    }
    linear_variable_t& operator-=(const linear_variable_t& w1)
    {
        // handling zero case
        if(c_.rows() == 0)
        {
            this->B_ = w1.B_;
            this->c_ = w1.c_;
        }
        else
        {
            this->B_ -= w1.B_;
            this->c_ -= w1.c_;
        }
        return *this;
    }

    static linear_variable_t Zero(size_t /*dim=0*/){
        return linear_variable_t(matrix_dim_t::Zero(), vectord_t::Zero(0));
    }

    const matrix_dim_x_t& B() const {return B_;}
    const point_dim_t& c () const {return c_;}

private:
    matrix_dim_x_t B_;
    point_dim_t c_;
};

template <typename Numeric=double>
struct quadratic_variable
{
    typedef Eigen::Matrix<Numeric, Eigen::Dynamic, Eigen::Dynamic> matrix_x_t;
    typedef Eigen::Matrix<Numeric, Eigen::Dynamic, 1> point_t;

    quadratic_variable(const matrix_x_t& A, const point_t& b, const Numeric c = 0):
        c_(c),
        b_(b),
        A_(A){assert(A.cols() == b.rows()) && A.cols() == A.rows();}

    quadratic_variable(const point_t& b, const Numeric c = 0):
        c_(c),
        b_(b),
        A_(matrix_x_t::Identity(b.rows())){}

    quadratic_variable(const matrix_x_t& A, const Numeric c = 0):
        c_(c),
        b_(point_t::Zero(A.cols())),
        A_(matrix_x_t::Identity(b.rows())){assert(A.cols() == A.rows());}

    // linear evaluation
    Numeric operator()(const Eigen::Ref<const point_t>& val) const
    {
        return val.transpose() * A() * val + b().transpose() * val + c();
    }

    quadratic_variable& operator+=(const quadratic_variable& w1)
    {
        this->A_ += w1.A_;
        this->b_ += w1.b_;
        this->c_ += w1.c_;
        return *this;
    }
    quadratic_variable& operator-=(const quadratic_variable& w1)
    {
        this->A_ -= w1.A_;
        this->b_ -= w1.b_;
        this->c_ -= w1.c_;
        return *this;
    }

    matrix_x_t& A() const {return A_;}
    point_t&  b () const {return b_;}
    Numeric&  c () const {return c_;}

private:
    matrix_x_t A_;
    point_t b_;
    Numeric c_;
};

template <int D, typename N>
inline linear_variable<D,N> operator+(const linear_variable<D,N>& w1, const linear_variable<D,N>& w2)
{
    linear_variable<D,N> res(w1.B(), w2.c());
    return res+=w2;
}

template <int D, typename N>
linear_variable<D,N> operator-(const linear_variable<D,N>& w1, const linear_variable<D,N>& w2)
{
    linear_variable<D,N> res(w1.B(), w2.c());
    return res-=w2;
}

template <int D, typename N>
linear_variable<D,N> operator*(const double k, const linear_variable<D,N>& w){
    return linear_variable<D,N>(k*w.B(),k*w.c());
}

template <int D, typename N>
linear_variable<D,N> operator*(const linear_variable<D,N>& w,const double k){
    return linear_variable<D,N>(k*w.B(),k*w.c());
}

template <int D, typename N>
linear_variable<D,N> operator/(const linear_variable<D,N>& w,const double k){
    return linear_variable<D,N>(w.B()/k,w.c()/k);
}

} // namespace spline
#endif //_CLASS_LINEAR_VARIABLE

