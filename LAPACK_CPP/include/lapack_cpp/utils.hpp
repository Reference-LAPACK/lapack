#ifndef LAPACK_CPP_UTILS_HPP
#define LAPACK_CPP_UTILS_HPP

#include <iostream>
#include <type_traits>

#include "lapack_cpp/base.hpp"

namespace lapack_cpp {

// Code for printing
template <typename T, typename idx_t>
void print(const Vector<T, idx_t>& v)
{
    std::cout << "(" << v.size() << ")[" << std::endl;
    for (idx_t i = 0; i < v.size(); ++i) {
        std::cout << v[i] << std::endl;
    }
    std::cout << "]" << std::endl;
}

// Code for printing
template <typename T, Layout layout, typename idx_t>
void print(const Matrix<T, layout, idx_t>& m)
{
    std::cout << "(" << m.num_rows() << "," << m.num_columns() << ")["
              << std::endl;
    for (idx_t i = 0; i < m.num_rows(); ++i) {
        std::cout << "[";
        for (idx_t j = 0; j < m.num_columns() - 1; ++j) {
            std::cout << m(i, j) << ",";
        }
        std::cout << m(i, m.num_columns() - 1) << "]" << std::endl;
    }
    std::cout << "]" << std::endl;
}

// Code for printing
template <typename T, Layout layout, typename idx_t>
void print(const ConstMatrix<T, layout, idx_t>& m)
{
    std::cout << "(" << m.num_rows() << "," << m.num_columns() << ")["
              << std::endl;
    for (idx_t i = 0; i < m.num_rows(); ++i) {
        std::cout << "[";
        for (idx_t j = 0; j < m.num_columns() - 1; ++j) {
            std::cout << m(i, j) << ",";
        }
        std::cout << m(i, m.num_columns() - 1) << "]" << std::endl;
    }
    std::cout << "]" << std::endl;
}

// Initialize a matrix with random values
template <typename T, Layout layout, typename idx_t>
void randomize(Matrix<T, layout, idx_t>& m)
{
    for (idx_t j = 0; j < m.num_columns(); ++j) {
        for (idx_t i = 0; i < m.num_rows(); ++i) {
            m(i, j) = (T)rand() / (T)RAND_MAX;
        }
    }
}

template <typename T, Layout layout, typename idx_t>
void randomize(Matrix<std::complex<T>, layout, idx_t>& m)
{
    for (idx_t j = 0; j < m.num_columns(); ++j) {
        for (idx_t i = 0; i < m.num_rows(); ++i) {
            m(i, j) = std::complex<T>((T)rand() / (T)RAND_MAX,
                                      (T)rand() / (T)RAND_MAX);
        }
    }
}

// Initialize a vector with random values
template <typename T, typename idx_t>
void randomize(Vector<T, idx_t>& v)
{
    for (idx_t i = 0; i < v.size(); ++i) {
        v[i] = (T)rand() / (T)RAND_MAX;
    }
}

}  // namespace lapack_cpp

#endif