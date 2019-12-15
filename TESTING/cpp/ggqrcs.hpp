/*
 * Copyright (c) 2016, 2019, Christoph Conrads (https://christoph-conrads.name)
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
#include <random>
#include <type_traits>



namespace ublas = boost::numeric::ublas;
typedef lapack::integer_t Integer;


template<
	typename Real,
	std::enable_if<std::is_fundamental<Real>::value>* = nullptr
>
bool nan_p(Real x)
{
	return std::isnan(x);
}

template<typename Real>
bool nan_p(std::complex<Real> z)
{
	return std::isnan(z.real()) or std::isnan(z.imag());
}

template<
	typename Real,
	std::enable_if<std::is_fundamental<Real>::value>* = nullptr
>
bool finite_p(Real x)
{
	return std::isfinite(x);
}

template<typename Real>
bool finite_p(std::complex<Real> z)
{
	return std::isfinite(z.real()) and std::isfinite(z.imag());
}




template<typename T, class Storage>
ublas::matrix<T, Storage> build_R(
	std::size_t r,
	const ublas::matrix<T, Storage>& X, const ublas::matrix<T, Storage>& Y)
{
	BOOST_ASSERT( X.size2() == Y.size2() );

	typedef ublas::matrix<T, Storage> Matrix;
	typedef ublas::matrix_range<Matrix> MatrixRange;
	typedef ublas::matrix_range<const Matrix> ConstMatrixRange;
	typedef ublas::banded_adaptor<ConstMatrixRange> BandedAdaptor;

	const std::size_t m = X.size1();
	const std::size_t n = X.size2();

	Matrix R(r, n, 0);

	if(r <= m)
	{
		MatrixRange R12 = ublas::subrange(R, 0, r, n-r, n);
		ConstMatrixRange X1 = ublas::subrange(X, 0, r, 0, r);
		BandedAdaptor X1U(X1, 0, r);

		R12 = X1U;
	}
	else
	{
		MatrixRange R12 = ublas::subrange(R, 0, m, n-r, n);
		MatrixRange R22 = ublas::subrange(R, m, r, n+m-r, n);

		ConstMatrixRange X1 = ublas::subrange(X, 0, m, 0, r);
		ConstMatrixRange Y1 = ublas::subrange(Y, 0, r-m, 0, r-m);

		BandedAdaptor X1U(X1, 0, r);
		BandedAdaptor Y1U(Y1, 0, r);

		R12 = X1U;
		R22 = Y1U;
	}

	return R;
}



template<typename Real, typename T, class Storage>
std::pair< ublas::matrix<T, Storage>, ublas::matrix<T, Storage> >
build_diagonals(
	std::size_t r,
	const ublas::matrix<T, Storage>& A, const ublas::matrix<T, Storage>& B,
	const ublas::vector<Real>& theta)
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



template<typename Real, typename T, class Storage>
void check_results(
	Integer ret, Real lwkopt,
	const ublas::matrix<T, Storage>& A, const ublas::matrix<T, Storage>& B,
	Real w, Integer l,
	const ublas::vector<Real> theta,
	const ublas::matrix<T, Storage>& U1, const ublas::matrix<T, Storage>& U2,
	const ublas::matrix<T, Storage>& Qt,
	const ublas::matrix<T, Storage>& X, const ublas::matrix<T, Storage>& Y)
{
	BOOST_REQUIRE( A.size2() == B.size2() );
	BOOST_REQUIRE( A.size1() == U1.size1() );
	BOOST_REQUIRE( B.size1() == U2.size1() );
	BOOST_REQUIRE( A.size2() == Qt.size1() );


	typedef ublas::matrix<T, Storage> Matrix;


	const std::size_t m = A.size1();
	const std::size_t n = A.size2();
	const std::size_t p = B.size1();
	const auto eps = std::numeric_limits<Real>::epsilon();


	// check scalars
	BOOST_CHECK_EQUAL( ret, 0 );
	BOOST_REQUIRE( !nan_p(lwkopt) );
	BOOST_REQUIRE( std::isfinite(lwkopt) );
	BOOST_CHECK_GT( lwkopt, (m+p) * n );

	BOOST_REQUIRE( !nan_p(w) );
	BOOST_REQUIRE( std::isfinite(w) );
	BOOST_REQUIRE_GT( w, 0 );

	BOOST_CHECK_GE( l, 0 );
	BOOST_CHECK_LE( l, std::min(m+p, n) );
	const std::size_t r = l;
	const std::size_t k = std::min( {m, p, r, m + p - r} );


	// construct R
	const Matrix R = build_R(r, X, Y);


	// check for NaN
	bool(*isnan)(T) = &nan_p;
	BOOST_REQUIRE( std::none_of( A.data().begin(), A.data().end(), isnan) );
	BOOST_REQUIRE( std::none_of( B.data().begin(), B.data().end(), isnan) );
	BOOST_REQUIRE( std::none_of( R.data().begin(), R.data().end(), isnan) );
	BOOST_REQUIRE( std::none_of( U1.data().begin(), U1.data().end(), isnan) );
	BOOST_REQUIRE( std::none_of( U2.data().begin(), U2.data().end(), isnan) );
	BOOST_REQUIRE( std::none_of( Qt.data().begin(), Qt.data().end(), isnan) );
	if( k > 0 )
	{
		bool(*isnan_r)(Real) = &nan_p;
		BOOST_REQUIRE( std::none_of( &theta(0), &theta(0)+k, isnan_r) );
	}


	// check for infinity
	bool(*is_inf)(T) = &finite_p;
	BOOST_REQUIRE( std::none_of( A.data().begin(), A.data().end(), is_inf) );
	BOOST_REQUIRE( std::none_of( B.data().begin(), B.data().end(), is_inf) );
	BOOST_REQUIRE( std::none_of( R.data().begin(), R.data().end(), is_inf) );
	BOOST_REQUIRE( std::none_of( U1.data().begin(), U1.data().end(), is_inf) );
	BOOST_REQUIRE( std::none_of( U2.data().begin(), U2.data().end(), is_inf) );
	BOOST_REQUIRE( std::none_of( Qt.data().begin(), Qt.data().end(), is_inf) );
	if( k > 0 )
	{
		bool(*is_inf)(Real) = &finite_p;
		BOOST_REQUIRE( std::none_of( &theta(0), &theta(0)+k, is_inf) );
	}


	// check that unitary matrices are indeed unitary
	// The bound is based on Inequality (19.13), Equation (3.8) in
	// Higham: "Accuracy and Stability of Numerical Algorithms".
	// Note that only Qt is the orthogonal factor of a QR decomposition.
	auto measure_unity = [] (const auto& U) -> double
	{
		BOOST_ASSERT( U.size1() == U.size2() );

		const std::size_t n = U.size1();
		const ublas::identity_matrix<T> id(n);

		double ret =
			ublas::norm_frobenius( ublas::prod( ublas::herm(U), U ) - id );
		return ret;
	};

	BOOST_CHECK_LE( measure_unity(U1), 2 * std::sqrt(m) * (m+p) * r * eps );
	BOOST_CHECK_LE( measure_unity(U2), 2 * std::sqrt(p) * (m+p) * r * eps );
	BOOST_CHECK_LE( measure_unity(Qt), 2 * std::sqrt(n) * n * r * eps );


	// check the "singular values"
	BOOST_CHECK_GE( theta.size(), k );

	for(std::size_t i = 0; i < k; ++i)
	{
		BOOST_CHECK_GE( theta[i], 0 );
		BOOST_CHECK_LE( theta[i], Real(M_PI/2) );

		if( i > 0 )
			BOOST_CHECK_LE( theta[i-1], theta[i] );
	}


	// reconstruct A, B from GSVD
	const auto ds = build_diagonals(r, A, B, theta);
	const Matrix& D1 = ds.first;
	const Matrix& D2 = ds.second;

	const Matrix almost_A = reconstruct_matrix(U1, D1, R, Qt);
	const Matrix almost_B = reconstruct_matrix(U2, D2, R, Qt);

	const Real frob_A = ublas::norm_frobenius(A);
	const Real frob_B = ublas::norm_frobenius(B);

	// The tolerance here is based on the backward error bounds for the QR
	// factorization given in Theorem 19.4, Equation (3.8) in
	// Higham: "Accuracy and Stability of Numerical Algorithms".
	BOOST_CHECK_LE(
		ublas::norm_frobenius(A - almost_A), 2 * (m+p) * n * frob_A * eps );
	BOOST_CHECK_LE(
		ublas::norm_frobenius(w*B - almost_B), 3*w * (m+p) * n * frob_B * eps );
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
		BOOST_ASSERT( m > 0 );
		BOOST_ASSERT( n > 0 );
		BOOST_ASSERT( p > 0 );

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
		BOOST_ASSERT( ret == 0 );

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
	T w = -1;
	Integer l = -1;

	Integer ret = lapack::ggqrcs(
		'Y', 'Y', 'Y', m, n, p, &w, &l,
		&A(0, 0), m, &B(0, 0), p,
		&theta(0),
		&U1(0, 0), m, &U2(0, 0), p, &Qt(0, 0), n,
		&work(0), lwork, &iwork(0) );

	check_results(
		ret, work(0), fixture.A, fixture.B, w, l, theta, U1, U2, Qt, A, B);

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
	T w = -1;
	Integer l = -1;

	Integer ret = lapack::ggqrcs(
		'Y', 'Y', 'Y', m, n, p, &w, &l,
		&A(0, 0), m, &B(0, 0), p,
		&theta(0),
		&U1(0, 0), m, &U2(0, 0), p, &Qt(0, 0), n,
		&work(0), lwork, &iwork(0) );

	check_results(
		ret, work(0), fixture.A, fixture.B, w, l, theta, U1, U2, Qt, A, B);

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
				A(1,0) = 1;
				B(0,1) = 1;
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
				T w = -1;
				Integer l = -1;

				Integer ret = lapack::ggqrcs(
					'Y', 'Y', 'Y', m, n, p, &w, &l,
					&A(0, 0), m, &B(0, 0), p,
					&theta(0),
					&U1(0, 0), m, &U2(0, 0), p, &Qt(0, 0), n,
					&work(0), lwork, &iwork(0) );

				check_results(
					ret, work(0), fixture.A, fixture.B,
					w, l, theta, U1, U2, Qt, A, B);

				BOOST_CHECK_EQUAL( w, 1 );
			}
		}
	}
}



BOOST_AUTO_TEST_CASE_TEMPLATE(
	ggqrcs_random_test, T, test_types)
{
	const std::size_t dimensions[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

	for(std::size_t m : dimensions)
	{
		for(std::size_t n : dimensions)
		{
			for(std::size_t p : dimensions)
			{
				std::mt19937_64 gen( 1u );
				std::uniform_real_distribution<T> dist(-1, +1);
				auto make_rand = [&gen, &dist] (T s)
				{
					auto rand = [&gen, &dist, s] () { return s*dist(gen); };
					return rand;
				};


				for(std::size_t iter = 0; iter < 100; ++iter)
				{
					Fixture<T> fixture(m, n, p);

					auto& A = fixture.A;
					auto& B = fixture.B;

					std::generate(
						A.data().begin(), A.data().end(), make_rand(1000) );
					std::generate(
						B.data().begin(), B.data().end(), make_rand(1) );

					auto& theta = fixture.theta;
					auto& U1 = fixture.U1;
					auto& U2 = fixture.U2;
					auto& Qt = fixture.Qt;
					auto& work = fixture.work;
					auto& iwork = fixture.iwork;

					const Integer lwork = work.size();
					auto w = T{-1};
					auto l = Integer{-1};

					Integer ret = lapack::ggqrcs(
						'Y', 'Y', 'Y', m, n, p, &w, &l,
						&A(0, 0), m, &B(0, 0), p,
						&theta(0),
						&U1(0, 0), m, &U2(0, 0), p, &Qt(0, 0), n,
						&work(0), lwork, &iwork(0) );

					check_results(
						ret, work(0), fixture.A, fixture.B,
						w, l, theta, U1, U2, Qt, A, B);
				}
			}
		}
	}
}



template<typename Real>
struct ComplexFixture
{
	using T = std::complex<Real>;

	using Matrix = ublas::matrix<T, ublas::column_major>;
	template<typename U> using Vector = ublas::vector<U>;

	static constexpr Real NaN = std::numeric_limits<Real>::quiet_NaN();

	Matrix A, B;
	Matrix U1, U2, Qt;
	Vector<Real> theta;
	Vector<T> work;
	Vector<Real> rwork;
	Vector<Integer> iwork;


	ComplexFixture(std::size_t m, std::size_t n, std::size_t p) :
		A(m, n, 0),
		B(p, n, 0),
		U1(m, m, NaN),
		U2(p, p, NaN),
		Qt(n, n, NaN),
		theta(n, NaN),
		rwork(2*n, NaN),
		iwork(m + n + p, -1)
	{
		BOOST_ASSERT( m > 0 );
		BOOST_ASSERT( n > 0 );
		BOOST_ASSERT( p > 0 );

		// query workspace size
		auto lwork_opt_f = T(NaN, NaN);
		auto lrwork_opt_f = NaN;
		auto w = NaN;
		auto l = Integer{-1};
		auto lwork = Integer{-1};
		auto lrwork = Integer(2*n);

		Integer ret = lapack::ggqrcs(
			'Y', 'Y', 'Y', m, n, p, &w, &l,
			&A(0, 0), m, &B(0, 0), p,
			&theta(0),
			&U1(0, 0), m, &U2(0, 0), p, &Qt(0, 0), n,
			&lwork_opt_f, lwork, &lrwork_opt_f, lrwork, &iwork(0) );
		BOOST_ASSERT( ret == 0 );
		BOOST_ASSERT( lwork_opt_f.real() >= 0 );
		BOOST_ASSERT( lrwork_opt_f >= 2*n );

		// resize work accordingly
		std::size_t lwork_opt = lwork_opt_f.real();
		work.resize( lwork_opt );

		std::size_t lrwork_opt = lrwork_opt_f;
		rwork.resize( lrwork_opt );

		std::fill( work.begin(), work.end(), NaN );
		std::fill( rwork.begin(), rwork.end(), NaN );
	}
};

template<typename Real> constexpr Real ComplexFixture<Real>::NaN;


BOOST_AUTO_TEST_CASE_TEMPLATE(
	complex_ggqrcs_random_test, Real, test_types)
{
	using T = std::complex<Real>;

	const std::size_t dimensions[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

	for(std::size_t m : dimensions)
	{
		for(std::size_t n : dimensions)
		{
			for(std::size_t p : dimensions)
			{
				using distribution_t = std::uniform_real_distribution<Real>;
				auto gen = std::mt19937_64( 1u );
				auto dist = distribution_t(-M_SQRT1_2, +M_SQRT1_2);
				auto make_rand = [&gen, &dist] (T s)
				{
					auto rand = [&gen, &dist, s] () {
						return s*T(dist(gen), dist(gen));
					};
					return rand;
				};


				for(std::size_t iter = 0; iter < 100; ++iter)
				{
					auto fixture = ComplexFixture<Real>(m, n, p);
					auto A = fixture.A;
					auto B = fixture.B;

					std::generate(
						A.data().begin(), A.data().end(), make_rand(1000) );
					std::generate(
						B.data().begin(), B.data().end(), make_rand(1) );

					fixture.A = A;
					fixture.B = B;

					auto& theta = fixture.theta;
					auto& U1 = fixture.U1;
					auto& U2 = fixture.U2;
					auto& Qt = fixture.Qt;
					auto& work = fixture.work;
					auto& rwork = fixture.rwork;
					auto& iwork = fixture.iwork;

					const Integer lwork = work.size();
					const Integer lrwork = rwork.size();
					auto w = Real{-1};
					auto l = Integer{-1};

					Integer ret = lapack::ggqrcs(
						'Y', 'Y', 'Y', m, n, p, &w, &l,
						&A(0, 0), m, &B(0, 0), p,
						&theta(0),
						&U1(0, 0), m, &U2(0, 0), p, &Qt(0, 0), n,
						&work(0), lwork,
						&rwork(0), lrwork, &iwork(0) );

					check_results(
						ret, work(0).real(), fixture.A, fixture.B,
						w, l, theta, U1, U2, Qt, A, B);
				}
			}
		}
	}
}

BOOST_AUTO_TEST_SUITE_END()

#endif
