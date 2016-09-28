/*
 * Copyright (c) 2016, Christoph Conrads (https://christoph-conrads.name)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *  * Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *  * Neither the name of the copyright holders nor the
 *    names of its contributors may be used to endorse or promote products
 *    derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */
#ifndef LAPACK_TESTS_GGQRCS_HPP
#define LAPACK_TESTS_GGQRCS_HPP

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <lapack.hpp>

#include <utility>

#include <algorithm>
#include <limits>
#include <cmath>

#include <iostream>
#include <boost/numeric/ublas/io.hpp>



namespace ublas = boost::numeric::ublas;
typedef lapack::integer_t Integer;


template<typename T, class Storage>
ublas::matrix<T, Storage> build_R(
	std::size_t r,
	const ublas::matrix<T, Storage>& X, const ublas::matrix<T, Storage>& Y)
{
	assert( X.size2() == Y.size2() );

	typedef ublas::matrix<T, Storage> Matrix;
	typedef ublas::matrix_range<Matrix> MatrixRange;

	const std::size_t m = X.size1();
	const std::size_t n = X.size2();

	Matrix R(r, n, 0);

	if(r <= m)
	{
		MatrixRange R12 = ublas::subrange(R, 0, r, n-r, n);
		const auto& X1 = ublas::subrange(X, 0, r, 0, r);
		R12 = X1;
	}
	else
	{
		MatrixRange R12 = ublas::subrange(R, 0, m, n-r, n);
		MatrixRange R22 = ublas::subrange(R, m, r, n-r, n);
		const auto& X1 = ublas::subrange(X, 0, m, 0, r);
		const auto& Y1 = ublas::subrange(Y, 0, r-m, 0, r);

		R12 = X1;
		R22 = Y1;
	}

	return R;
}



template<typename T, class Storage>
std::pair< ublas::matrix<T, Storage>, ublas::matrix<T, Storage> >
build_diagonals(
	std::size_t r,
	const ublas::matrix<T, Storage>& A, const ublas::matrix<T, Storage>& B,
	const ublas::vector<T>& theta)
{
	typedef ublas::matrix<T, Storage> Matrix;
	typedef ublas::identity_matrix<T> IdentityMatrix;
	typedef ublas::matrix_range<Matrix> MatrixRange;

	const std::size_t m = A.size1();
	const std::size_t p = B.size1();
	const std::size_t k = std::min( {m, p, r, m + p - r} );
	const std::size_t k1 = (p < r) ? r - p : 0;
	const std::size_t k2 = (m < r) ? r - m : 0;

	Matrix D1(m, r, 0);
	Matrix D2(p, r, 0);

	if(k1 > 0)
	{
		MatrixRange D1_33 = ublas::subrange(D1, m - k1, m, r - k1, r);
		D1_33 = IdentityMatrix(k1);
	}

	if(k2 > 0)
	{
		MatrixRange D2_11 = ublas::subrange(D2, 0, k2, 0, k2);
		D2_11 = IdentityMatrix(k2);
	}

	if(k > 0)
	{
		MatrixRange D1_22 = ublas::subrange(D1, m-k-k1, m-k1, r-k-k1, r-k1);
		MatrixRange D2_22 = ublas::subrange(D2, k2, k2+k, k2, k2+k);

		for(std::size_t i = 0; i < k; ++i)
		{
			D1_22(i, i) = std::sin( theta(i) );
			D2_22(i, i) = std::cos( theta(i) );
		}
	}

	return std::make_pair(D1, D2);
}



template<typename T, class Storage>
ublas::matrix<T, Storage> reconstruct_matrix(
	const ublas::matrix<T, Storage>& U,
	const ublas::matrix<T, Storage>& D,
	const ublas::matrix<T, Storage>& R,
	const ublas::matrix<T, Storage>& Qt)
{
	typedef ublas::matrix<T, Storage> Matrix;

	const Matrix D_R = ublas::prod(D, R);
	const Matrix U_D_R = ublas::prod(U, D_R);
	const Matrix A = ublas::prod(U_D_R, Qt);

	return A;
}



template<typename T, class Storage>
void check_results(
	Integer ret,
	const ublas::matrix<T, Storage>& A, const ublas::matrix<T, Storage>& B,
	T w, Integer l,
	const ublas::vector<T> theta,
	const ublas::matrix<T, Storage>& U1, const ublas::matrix<T, Storage>& U2,
	const ublas::matrix<T, Storage>& Qt,
	const ublas::matrix<T, Storage>& X, const ublas::matrix<T, Storage>& Y)
{
	assert( A.size2() == B.size2() );
	assert( A.size1() == U1.size1() );
	assert( B.size1() == U2.size1() );
	assert( A.size2() == Qt.size1() );

	const std::size_t m = A.size1();
	const std::size_t n = A.size2();
	const std::size_t p = B.size1();
	const T eps = std::numeric_limits<T>::epsilon();


	// check scalars
	BOOST_CHECK_EQUAL( ret, 0 );

	BOOST_CHECK( !std::isnan(w) );
	BOOST_CHECK( std::isfinite(w) );
	BOOST_CHECK_GT( w, 0 );

	BOOST_CHECK_GE( l, 0 );
	BOOST_CHECK_LE( l, std::min(m+p, n) );
	const std::size_t r = l;


	// check that unitary matrices are indeed unitary
	auto measure_unity = [] (const auto& U) -> double
	{
		assert( U.size1() == U.size2() );

		const std::size_t n = U.size1();
		const ublas::identity_matrix<T> id(n);

		double ret =
			ublas::norm_frobenius( ublas::prod( ublas::herm(U), U ) - id );
		return ret;
	};

	BOOST_CHECK_LE( measure_unity(U1), eps );
	BOOST_CHECK_LE( measure_unity(U2), eps );
	BOOST_CHECK_LE( measure_unity(Qt), eps );


	// check the "singular values"
	const std::size_t k = std::min( {m, p, r, m + p - r} );
	BOOST_CHECK_GE( theta.size(), k );

	for(std::size_t i = 0; i < k; ++i)
	{
		BOOST_CHECK_GE( theta[i], 0 );
		BOOST_CHECK_LE( theta[i], M_PI/2 );

		if( i > 0 )
			BOOST_CHECK_LE( theta[i-1], theta[i] );
	}


	// construct R
	typedef ublas::matrix<T, Storage> Matrix;
	Matrix R = build_R(r, X, Y);

	const T frob_R = ublas::norm_frobenius(R);
	const T frob_A = ublas::norm_frobenius(A);
	const T frob_B = ublas::norm_frobenius(B);
	const T frob_G = std::sqrt( std::pow(frob_A, 2) + std::pow(frob_B, 2) );

	BOOST_CHECK_CLOSE( frob_R, frob_G, frob_G * eps );

	// reconstruct A, B from GSVD
	const auto ds = build_diagonals(r, A, B, theta);
	const Matrix& D1 = ds.first;
	const Matrix& D2 = ds.second;

	const Matrix almost_A = reconstruct_matrix(U1, D1, R, Qt);
	const Matrix almost_B = reconstruct_matrix(U2, D2, R, Qt);

	BOOST_CHECK_LE( ublas::norm_frobenius(A - almost_A), frob_A * eps );
	BOOST_CHECK_LE( ublas::norm_frobenius(w*B - almost_B), frob_B * eps );
}



template<typename T>
struct Fixture
{
	typedef ublas::matrix<T, ublas::column_major> Matrix;
	template<typename U> using Vector = ublas::vector<U>;

	static constexpr T NaN = std::numeric_limits<T>::quiet_NaN();


	Fixture(std::size_t m, std::size_t n, std::size_t p) :
		A(m, n, 0),
		B(p, n, 0),
		U1(m, m, NaN),
		U2(p, p, NaN),
		Qt(n, n, NaN),
		theta(n, NaN),
		iwork(m + n + p, -1)
	{
		assert( m > 0 );
		assert( n > 0 );
		assert( p > 0 );

		// query workspace size
		T lwork_opt_f = NaN;
		T w = NaN;
		Integer l = -1;

		Integer ret = lapack::ggqrcs(
			'Y', 'Y', 'Y', m, n, p, &w, &l,
			&A(0, 0), m, &B(0, 0), p,
			&theta(0),
			&U1(0, 0), m, &U2(0, 0), p, &Qt(0, 0), n,
			&lwork_opt_f, -1, &iwork(0) );
		assert( ret == 0 );

		// resize work accordingly
		std::size_t lwork_opt = lwork_opt_f;
		work.resize( lwork_opt );

		std::fill( work.begin(), work.end(), NaN );
	}

	Matrix A, B;
	Matrix U1, U2, Qt;
	Vector<T> theta;
	Vector<T> work;
	Vector<Integer> iwork;
};


template<typename T> constexpr T Fixture<T>::NaN;



BOOST_AUTO_TEST_SUITE(LAPACK_TEST_SUITE_NAME)

BOOST_AUTO_TEST_CASE_TEMPLATE(
	ggqrcs_simple_test, T, test_types)
{
	const std::size_t m = 2;
	const std::size_t n = 2;
	const std::size_t p = 2;

	Fixture<T> fixture(m, n, p);

	auto A = fixture.A;
	auto B = fixture.B;

	A(0,0) = 1;
	B(1,1) = 1;

	fixture.A = A;
	fixture.B = B;

	auto& theta = fixture.theta;
	auto& U1 = fixture.U1;
	auto& U2 = fixture.U2;
	auto& Qt = fixture.Qt;
	auto& work = fixture.work;
	auto& iwork = fixture.iwork;

	const Integer lwork = work.size();
	double w = -1;
	Integer l = -1;

	Integer ret = lapack::ggqrcs(
		'Y', 'Y', 'Y', m, n, p, &w, &l,
		&A(0, 0), m, &B(0, 0), p,
		&theta(0),
		&U1(0, 0), m, &U2(0, 0), p, &Qt(0, 0), n,
		&work(0), lwork, &iwork(0) );

	check_results(ret, fixture.A, fixture.B, w, l, theta, U1, U2, Qt, A, B);

	BOOST_CHECK_EQUAL( w, 1 );
	BOOST_CHECK_EQUAL( l, 2 );
}


BOOST_AUTO_TEST_CASE_TEMPLATE(
	ggqrcs_zero_test, T, test_types)
{
	const std::size_t m = 4;
	const std::size_t n = 3;
	const std::size_t p = 2;

	Fixture<T> fixture(m, n, p);

	auto A = fixture.A;
	auto B = fixture.B;

	auto& theta = fixture.theta;
	auto& U1 = fixture.U1;
	auto& U2 = fixture.U2;
	auto& Qt = fixture.Qt;
	auto& work = fixture.work;
	auto& iwork = fixture.iwork;

	const Integer lwork = work.size();
	double w = -1;
	Integer l = -1;

	Integer ret = lapack::ggqrcs(
		'Y', 'Y', 'Y', m, n, p, &w, &l,
		&A(0, 0), m, &B(0, 0), p,
		&theta(0),
		&U1(0, 0), m, &U2(0, 0), p, &Qt(0, 0), n,
		&work(0), lwork, &iwork(0) );

	check_results(ret, fixture.A, fixture.B, w, l, theta, U1, U2, Qt, A, B);

	BOOST_CHECK_EQUAL( w, 1 );
	BOOST_CHECK_EQUAL( l, 0 );
}


BOOST_AUTO_TEST_CASE_TEMPLATE(
	ggqrcs_rectangular_test, T, test_types)
{
	for(std::size_t m : { 2, 13, 41 })
	{
		for(std::size_t n : {3, 7, 31})
		{
			for(std::size_t p : {5, 11, 17})
			{
				Fixture<T> fixture(m, n, p);

				auto A = fixture.A;
				auto B = fixture.B;

				A(0,0) = 1;
				B(1,1) = 1;

				fixture.A = A;
				fixture.B = B;

				auto& theta = fixture.theta;
				auto& U1 = fixture.U1;
				auto& U2 = fixture.U2;
				auto& Qt = fixture.Qt;
				auto& work = fixture.work;
				auto& iwork = fixture.iwork;

				const Integer lwork = work.size();
				double w = -1;
				Integer l = -1;

				Integer ret = lapack::ggqrcs(
					'Y', 'Y', 'Y', m, n, p, &w, &l,
					&A(0, 0), m, &B(0, 0), p,
					&theta(0),
					&U1(0, 0), m, &U2(0, 0), p, &Qt(0, 0), n,
					&work(0), lwork, &iwork(0) );

				check_results(
					ret, fixture.A, fixture.B, w, l, theta, U1, U2, Qt, A, B);

				BOOST_CHECK_EQUAL( w, 1 );
			}
		}
	}
}

BOOST_AUTO_TEST_SUITE_END()

#endif
