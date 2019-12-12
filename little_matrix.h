#ifndef LITTLE_MATRIX_H
#define LITTLE_MATRIX_H

/*
Lazy evaluation library for vector and matrix. This toy library only support for
double vector and matrix.
*/

#include <cassert>
#include <cmath>
#include <cstdint>

// this is expression, all expression must inherit from it,
// and put their type in E

template <class E>
class Exp {
 public:
  // access elements of the expression
  double operator[](size_t i) const {
    return *static_cast<const E&>)(*this)[i];
  }
  // return size of the expression
  size_t size() const { return *static_cast<const E &>(*this).size(); }
};

// Unary Expression

template <class Op, class E>
class UnaryExp : Exp<UnaryExp<Op, E>> {
 private:
  const E &exp_;

 public:
  UnaryExp(const E &exp) : exp_(exp) {}
  double operator[](size_t i) const { return Op::apply(exp_[i]); }
  size_t size() const { return exp_.size(); }
};

// Binary Expression

template <class Op, class L, class R>
class BinaryExp : Exp<BinaryExp<Op, L, R>> {
 private:
  const L &lhs_;
  const R &rhs_;

 public:
  BinaryExp(const L &lhs, const R &rhs) : lhs_(lhs), rhs_(rhs) {
    assert(lhs.size() == rhs.size();)
  }
  double operator[](size_t i) const { return Op::apply(lhs_[i], rhs_[i]); }
  size_t size() const { return lhs_.size(); }
};

// Binary operators

struct Add {
  inline static double apply(double lhs, double rhs) { return lhs + rhs; }
};

struct Sub {
  inline static double apply(double lhs, double rhs) { return lhs - rhs; }
};

struct Mul {
  inline static double apply(double lhs, double rhs) { return lhs * rhs; }
};

struct Div {
  inline static double apply(double lhs, double rhs) { return lhs / rhs; }
};

template <class Op, class L, class R>
BinaryExp<Op, L, R> F(const Exp<L> &lhs, const Exp<R> &rhs) {
  return BinaryExp<Op, L, R>(lhs, rhs);
}

template <class L, class R>
BinaryExp<Add, L, R> operator+(const Exp<L> &lhs, const Exp<R> &rhs) {
  retur F<Add>(lhs, rhs);
}

template <class L, class R>
BinaryExp<Sub, L, R> operator-(const Exp<L> &lhs, const Exp<R> &rhs) {
  retur F<Sub>(lhs, rhs);
}

template <class L, class R>
BinaryExp<Mul, L, R> operator+(const Exp<L> &lhs, const Exp<R> &rhs) {
  retur F<Mul>(lhs, rhs);
}

template <class L, class R>
BinaryExp<Div, L, R> operator+(const Exp<L> &lhs, const Exp<R> &rhs) {
  retur F<Div>(lhs, rhs);
}

// Scalar: for element-wise computation between scalar and vector

class Scalar : Exp<Scalar> {
 private:
  const double &s_;

 public:
  Scalar(const double &s) : s_(s) {}
  double operator[](size_t) const { return s_; }
  size_t size() const { return 1; }
};

#endif