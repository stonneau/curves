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
    typedef Eigen::Matrix<Numeric, Eigen::Dynamic, Eigen::Dynamic> matrix_x_t;
    typedef linear_variable<Dim, Numeric> linear_variable_t;

    linear_variable(): B_(matrix_dim_t::Identity()), c_(point_dim_t::Zero()), zero(false){} //variable
    linear_variable(const point_dim_t& c):B_(matrix_dim_t::Zero()),c_(c), zero(false) {} // constant
    linear_variable(const matrix_dim_x_t& B, const point_dim_t& c):B_(B),c_(c), zero(false) {} //mixed

    // linear evaluation
    point_dim_t operator()(const Eigen::Ref<const point_dim_t>& val) const
    {
        assert(!isZero());
        return B() * val + c();
    }

    linear_variable_t& operator+=(const linear_variable_t& w1)
    {
        // handling zero case
        if(isZero())
        {
            this->B_ = w1.B_;
            this->c_ = w1.c_;
            zero = w1.isZero();
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
        if(isZero())
        {
            this->B_ = -w1.B_;
            this->c_ = -w1.c_;
            zero = w1.isZero();
        }
        else
        {
            this->B_ -= w1.B_;
            this->c_ -= w1.c_;
        }
        return *this;
    }
    linear_variable_t& operator/=(const double d)
    {
        // handling zero case
        if(!isZero())
        {
            this->B_ /= d;
            this->c_ /= d;
        }
        return *this;
    }
    linear_variable_t& operator*=(const double d)
    {
        // handling zero case
        if(!isZero())
        {
            this->B_ *= d;
            this->c_ *= d;
        }
        return *this;
    }

    static linear_variable_t Zero(size_t /*dim=0*/){
        return linear_variable_t(matrix_x_t::Identity(Dim,Dim), point_dim_t::Zero());
    }

    const matrix_dim_x_t& B() const {return B_;}
    const point_dim_t& c () const {return c_;}
    const bool isZero () const {return zero;}

private:
    matrix_dim_x_t B_;
    point_dim_t c_;
    bool zero;
};

template <typename Numeric=double>
struct quadratic_variable
{
    typedef Eigen::Matrix<Numeric, Eigen::Dynamic, Eigen::Dynamic> matrix_x_t;
    typedef Eigen::Matrix<Numeric, Eigen::Dynamic, 1> point_t;
    typedef quadratic_variable<Numeric> quadratic_variable_t;

    quadratic_variable()
    {
        c_ = 0.;
        b_ = point_t::Zero(1);
        A_ = matrix_x_t::Zero(1,1);
        zero = true;
    }

    quadratic_variable(const matrix_x_t& A, const point_t& b, const Numeric c = 0):
        c_(c),
        b_(b),
        A_(A),
        zero(false){assert(A.cols() == b.rows() && A.cols() == A.rows());}

    quadratic_variable(const point_t& b, const Numeric c = 0):
        c_(c),
        b_(b),
        A_(matrix_x_t::Zero((int)(b.rows()),(int)(b.rows()))),
        zero(false){}

    static quadratic_variable_t Zero(size_t dim=0){
        return quadratic_variable_t();
    }

    // linear evaluation
    Numeric operator()(const Eigen::Ref<const point_t>& val) const
    {
        assert(!isZero());
        return val.transpose() * A() * val + b().transpose() * val + c();
    }


    quadratic_variable& operator+=(const quadratic_variable& w1)
    {
        if(isZero())
        {
            if(!w1.isZero())
            {
                this->A_ = w1.A_;
                this->b_ = w1.b_;
                this->c_ = w1.c_;
                zero = false;
            }
        }
        else if(!w1.isZero())
        {
            this->A_ += w1.A_;
            this->b_ += w1.b_;
            this->c_ += w1.c_;
        }
        return *this;
    }
    quadratic_variable& operator-=(const quadratic_variable& w1)
    {
        if(isZero())
        {
            if(!w1.isZero())
            {
                this->A_ = -w1.A_;
                this->b_ = -w1.b_;
                this->c_ = -w1.c_;
                zero = false;
            }
        }
        else if(!w1.isZero())
        {
            this->A_ -= w1.A_;
            this->b_ -= w1.b_;
            this->c_ -= w1.c_;
        }
        return *this;
    }

    quadratic_variable& operator/=(const double d)
    {
        // handling zero case
        if(!isZero())
        {
            this->A_ /= d;
            this->b_ /= d;
            this->c_ /= d;
        }
        return *this;
    }
    quadratic_variable& operator*=(const double d)
    {
        // handling zero case
        if(!isZero())
        {
            this->A_ *= d;
            this->b_ *= d;
            this->c_ *= d;
        }
        return *this;
    }

    const matrix_x_t& A() const {return A_;}
    const point_t&  b () const {return b_;}
    const Numeric  c () const {return c_;}
    const bool  isZero() const {return zero;}

private:
    Numeric c_;
    point_t b_;
    matrix_x_t A_;
    bool zero;
};


template <typename N>
Eigen::Matrix<N,Eigen::Dynamic,Eigen::Dynamic>
    to_diagonal(const Eigen::Ref<const Eigen::Matrix<N,Eigen::Dynamic,1> > vec)
{
    typedef typename Eigen::Matrix<N,Eigen::Dynamic,Eigen::Dynamic> matrix_t;
    matrix_t res(matrix_t::Zero(vec.rows(),vec.rows()));
    for(int i =0; i<vec.rows(); ++i)
        res(i,i) = vec(i);
    return res;
}

// only works with diagonal linear variables
template <int D, typename N>
inline quadratic_variable<N> operator*(const linear_variable<D,N>& w1, const linear_variable<D,N>& w2)
{
    typedef quadratic_variable<N> quad_var_t;
    typedef linear_variable<D,N> lin_var_t;
    typedef typename quad_var_t::matrix_x_t matrix_x_t;
    typedef typename quad_var_t::point_t point_t;
    typedef typename lin_var_t::point_dim_t point_dim_t;
    point_dim_t ones = point_dim_t::Ones();
    point_t b1   = w1.B().transpose()*ones, b2 = w2.B().transpose()*ones;
    matrix_x_t B1 = to_diagonal<N>(b1);
    matrix_x_t B2 = to_diagonal<N>(b2); //b1.array().square()
    // omitting all transposes since A matrices are diagonals
    matrix_x_t  A = B1.transpose() * B2;
    //point_t  b = B1.transpose() * b2 + B2.transpose() * b1;
    point_t  b = w1.c().transpose() * w2.B() + w2.c().transpose() * w1.B();
    N c = w1.c().transpose() * w2.c();
    return quad_var_t(A,b,c);
}

template <int D, typename N>
inline linear_variable<D,N> operator+(const linear_variable<D,N>& w1, const linear_variable<D,N>& w2)
{
    linear_variable<D,N> res(w1.B(), w1.c());
    return res+=w2;
}

template <int D, typename N>
linear_variable<D,N> operator-(const linear_variable<D,N>& w1, const linear_variable<D,N>& w2)
{
    linear_variable<D,N> res(w1.B(), w1.c());
    return res-=w2;
}

template <int D, typename N>
linear_variable<D,N> operator*(const double k, const linear_variable<D,N>& w){
    linear_variable<D,N> res(w.B(), w.c());
    return res*=k;
}

template <int D, typename N>
linear_variable<D,N> operator*(const linear_variable<D,N>& w,const double k){
    linear_variable<D,N> res(w.B(), w.c());
    return res*=k;
}

template <int D, typename N>
linear_variable<D,N> operator/(const linear_variable<D,N>& w,const double k){
    linear_variable<D,N> res(w.B(), w.c());
    return res/=k;
}

template <typename N>
inline quadratic_variable<N> operator+(const quadratic_variable<N>& w1, const quadratic_variable<N>& w2)
{
    quadratic_variable<N> res(w1.A(),w1.b(), w1.c());
    return res+=w2;
}

template <typename N>
quadratic_variable<N> operator-(const quadratic_variable<N>& w1, const quadratic_variable<N>& w2)
{
    quadratic_variable<N> res(w1.A(),w1.b(), w1.c());
    return res-=w2;
}

template <typename N>
quadratic_variable<N> operator*(const double k, const quadratic_variable<N>& w){
    quadratic_variable<N> res(w.A(),w.b(), w.c());
    return res*=k;
}

template <typename N>
quadratic_variable<N> operator*(const quadratic_variable<N>& w,const double k){
    quadratic_variable<N> res(w.A(),w.b(), w.c());
    return res*=k;
}

template <typename N>
quadratic_variable<N> operator/(const quadratic_variable<N>& w,const double k){
    quadratic_variable<N> res(w.A(),w.b(), w.c());
    return res/=k;
}
} // namespace spline
#endif //_CLASS_LINEAR_VARIABLE

