/*
 * Copyright (c) 2016, 2019, 2020 Christoph Conrads
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

#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

#include <lapack.hpp>

#include <algorithm>
#include <limits>
#include <cmath>
#include <cstdint>
#include <ctime>
#include <random>
#include <type_traits>
#include <utility>



namespace ublas = boost::numeric::ublas;

using Integer = lapack::integer_t;


template<typename Number>
bool nan_p(Number x)
{
	return std::isnan(std::real(x)) or std::isnan(std::imag(x));
}


template<typename Number>
bool inf_p(Number x)
{
	return std::isinf(std::real(x)) or std::isinf(std::imag(x));
}


template<typename Number>
bool finite_p(Number x)
{
	return std::isfinite(std::real(x)) and std::isfinite(std::imag(x));
}



/**
 * Given a type `T`, `real_from<T>::type` is `T` if `T` is a built-in type,
 * otherwise is returns the numeric type used to implement `T`.
 *
 * Examples:
 * * `real_from<float>::type == float`
 * * `real_from<std::complex<float>>::type == float`
 */
template<typename T> struct real_from {};

template<> struct real_from<float> { using type = float; };
template<> struct real_from<double> { using type = double; };
template<typename Real>
struct real_from<std::complex<Real>> { using type = Real; };

static_assert(std::is_same<typename real_from<float>::type, float>::value, "");
static_assert(std::is_same<typename real_from<double>::type,double>::value, "");
static_assert(
	std::is_same<typename real_from<std::complex<float>>::type,float>::value, ""
);



template<typename Real>
struct not_a_number
{
	static constexpr Real value = std::numeric_limits<Real>::quiet_NaN();
};

template<typename Real>
constexpr Real not_a_number<Real>::value;


template<typename Real>
struct not_a_number<std::complex<Real>>
{
	using type = std::complex<Real>;

	static constexpr Real nan = std::numeric_limits<Real>::quiet_NaN();
	static constexpr type value = type{nan, nan};
};

template<typename Real>
constexpr std::complex<Real> not_a_number<std::complex<Real>>::value;



template<typename Real>
struct UniformDistribution
{
	UniformDistribution() : dist_(-1, +1) {}

	template<typename Engine>
	Real operator() (Engine& gen)
	{
		return dist_(gen);
	}

	std::uniform_real_distribution<Real> dist_;
};

template<typename Real>
struct UniformDistribution<std::complex<Real>>
{
	using result_type = std::complex<Real>;

	UniformDistribution() : dist_(-1, +1) {}

	template<typename Engine>
	result_type operator() (Engine& gen)
	{
		constexpr auto pi = Real{M_PI};
		auto radius = dist_(gen);
		auto angle = 2 * pi * dist_(gen);

		return std::polar(radius, angle);
	}

	std::uniform_real_distribution<Real> dist_;
};



template<typename T, class Storage>
ublas::matrix<T, Storage> build_R(
	std::size_t r,
	const ublas::matrix<T, Storage>& X, const ublas::matrix<T, Storage>& Y)
{
	BOOST_VERIFY( X.size2() == Y.size2() );

	using Matrix = ublas::matrix<T, Storage>;
	using MatrixRange = ublas::matrix_range<Matrix>;
	using ConstMatrixRange = ublas::matrix_range<const Matrix>;
	using BandedAdaptor = ublas::banded_adaptor<ConstMatrixRange>;

	auto m = X.size1();
	auto n = X.size2();
	auto R = Matrix(r, n, 0);

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



template<
	typename Number,
	class Storage = ublas::column_major,
	typename Real = typename real_from<Number>::type
>
std::pair< ublas::matrix<Number, Storage>, ublas::matrix<Number, Storage> >
build_diagonals_like(
	Number,
	std::size_t m, std::size_t p, std::size_t r,
	const ublas::vector<Real>& theta)
{
	using Matrix = ublas::matrix<Number, Storage>;
	using IdentityMatrix = ublas::identity_matrix<Number>;
	using MatrixRange = ublas::matrix_range<Matrix>;

	auto k = std::min( {m, p, r, m + p - r} );
	auto k1 = (p < r) ? r - p : std::size_t{0};
	auto k2 = (m < r) ? r - m : std::size_t{0};

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



template<typename Number, class Storage>
ublas::matrix<Number, Storage> reconstruct_matrix(
	const ublas::matrix<Number, Storage>& U,
	const ublas::matrix<Number, Storage>& D,
	const ublas::matrix<Number, Storage>& R,
	const ublas::matrix<Number, Storage>& Qt)
{
	using Matrix = ublas::matrix<Number, Storage>;

	auto D_R = Matrix(ublas::prod(D, R));
	auto U_D_R = Matrix(ublas::prod(U, D_R));
	auto A = Matrix(ublas::prod(U_D_R, Qt));

	return A;
}



/**
 * @return The Frobenius norm of U* U - I
 */
template<
	typename Number,
	typename Real = typename real_from<Number>::type
>
Real measure_isometry(const ublas::matrix<Number, ublas::column_major>& U)
{
	using Matrix = ublas::matrix<Number, ublas::column_major>;

	BOOST_VERIFY( U.size1() >= U.size2() );

	if(U.size2() == 0)
		return 0;

	auto m = U.size1();
	auto n = U.size2();
	auto J = Matrix(n, n);
	auto I = ublas::identity_matrix<Number>(n);
	auto alpha = Number{1};
	auto beta = Number{0};
	auto ret = lapack::gemm(
		'C', 'N', n, n, m, alpha, &U(0,0), m, &U(0,0), m, beta, &J(0,0), n
	);

	BOOST_VERIFY( ret == 0 );

	return ublas::norm_frobenius(J - I);
}


using real_test_types = boost::mpl::list<float,double>;

BOOST_AUTO_TEST_CASE_TEMPLATE(
	test_measure_isometry_simple_real, Real, real_test_types)
{
	for(auto m = std::size_t{0}; m < 5; ++m)
	{
		for(auto n = std::size_t{0}; n <= m; ++n)
		{
			auto A = ublas::matrix<Real, ublas::column_major>(m, n);

			std::fill( A.data().begin(), A.data().end(), Real{0} );

			for(auto i = std::size_t{0}; i < n; ++i)
			{
				A(i,i) = std::pow(Real{-1}, Real(i));
			}

			BOOST_CHECK_EQUAL( 0, measure_isometry(A) );
		}
	}
}

BOOST_AUTO_TEST_CASE_TEMPLATE(
	test_measure_isometry_real, Real, real_test_types)
{
	auto A = ublas::matrix<Real, ublas::column_major>(4, 2);

	// column-major order is required for this statement to work
	std::iota( A.data().begin(), A.data().end(), 1u );

	auto I = ublas::identity_matrix<Real>(2);
	auto AT_A = ublas::matrix<Real>(2, 2);

	AT_A(0,0) = 30; AT_A(0,1) = 70;
	AT_A(1,0) = 70; AT_A(1,1) =174;

	auto expected_result = ublas::norm_frobenius(AT_A - I);

	BOOST_CHECK_EQUAL( expected_result, measure_isometry(A) );
}


/**
 * This function checks if a matrix A might be considered orthogonal or unitary,
 * respectively, by comparing the Frobenius norm of `A*A - I` to a cut-off
 * value.
 *
 * The cut-off value is based on Inequality (19.13), Equation (3.8) in Higham:
 * "Accuracy and Stability of Numerical Algorithms".
 */
template<
	typename Number,
	class Storage,
	typename Real = typename real_from<Number>::type
>
bool is_almost_isometric(
	const ublas::matrix<Number, Storage>& U, Real multiplier = 4)
{
	BOOST_VERIFY( multiplier >= 1 );

	constexpr auto eps = std::numeric_limits<Real>::epsilon();
	auto m = U.size1();
	auto n = U.size2();
	auto p = std::min(m, n);

	if(p == 0)
		return true;

	auto r = measure_isometry(U);
	auto tol = std::sqrt(p) * m * n * eps;

	return r <= multiplier * tol;
}



template<
	typename Number,
	class Storage,
	typename Real = typename real_from<Number>::type
>
void check_results(
	Integer ret,
	const ublas::matrix<Number, Storage>& A,
	const ublas::matrix<Number, Storage>& B,
	Real w, Integer rank,
	const ublas::vector<Real> theta,
	const ublas::matrix<Number, Storage>& U1,
	const ublas::matrix<Number, Storage>& U2,
	const ublas::matrix<Number, Storage>& Qt,
	const ublas::matrix<Number, Storage>& X,
	const ublas::matrix<Number, Storage>& Y)
{
	BOOST_REQUIRE( A.size2() == B.size2() );
	BOOST_REQUIRE( A.size1() == U1.size1() );
	BOOST_REQUIRE( B.size1() == U2.size1() );
	BOOST_REQUIRE( A.size2() == Qt.size1() );

	using Matrix = ublas::matrix<Number, ublas::column_major>;

	auto m = A.size1();
	auto n = A.size2();
	auto p = B.size1();
	auto eps = std::numeric_limits<Real>::epsilon();

	// check scalars
	BOOST_CHECK_EQUAL( ret, 0 );

	BOOST_REQUIRE( !nan_p(w) );
	BOOST_REQUIRE( std::isfinite(w) );
	BOOST_REQUIRE_GT( w, 0 );

	BOOST_CHECK_GE( rank, 0 );
	BOOST_CHECK_LE( rank, std::min(m+p, n) );

	auto r = static_cast<std::size_t>(rank);
	auto k = std::min( {m, p, r, m + p - r} );

	BOOST_REQUIRE_GE(theta.size(), k);


	// construct R
	auto R = build_R(r, X, Y);

	BOOST_REQUIRE_EQUAL(R.size1(), r);
	BOOST_REQUIRE_EQUAL(R.size2(), n);


	// check for NaN
	bool(*isnan)(Number) = &nan_p;
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
	bool(*is_inf)(Number) = &inf_p;
	BOOST_REQUIRE( std::none_of( A.data().begin(), A.data().end(), is_inf) );
	BOOST_REQUIRE( std::none_of( B.data().begin(), B.data().end(), is_inf) );
	BOOST_REQUIRE( std::none_of( R.data().begin(), R.data().end(), is_inf) );
	BOOST_REQUIRE( std::none_of( U1.data().begin(), U1.data().end(), is_inf) );
	BOOST_REQUIRE( std::none_of( U2.data().begin(), U2.data().end(), is_inf) );
	BOOST_REQUIRE( std::none_of( Qt.data().begin(), Qt.data().end(), is_inf) );
	if( k > 0 )
	{
		bool(*is_inf)(Real) = &inf_p;
		BOOST_REQUIRE( std::none_of( &theta(0), &theta(0)+k, is_inf) );
	}


	// check that unitary matrices are indeed unitary
	// The bound is based on Inequality (19.13), Equation (3.8) in
	// Higham: "Accuracy and Stability of Numerical Algorithms".
	BOOST_CHECK_LE( measure_isometry(U1), 4 * std::sqrt(m) * (m+p) * r * eps );
	BOOST_CHECK_LE( measure_isometry(U2), 4 * std::sqrt(p) * (m+p) * r * eps );
	BOOST_CHECK_LE( measure_isometry(Qt), 4 * std::sqrt(n) * n * r * eps );


	// check the "singular values"
	BOOST_CHECK_GE( theta.size(), k );

	for(auto i = std::size_t{0}; i < k; ++i)
	{
		BOOST_CHECK_GE( theta[i], 0 );
		BOOST_CHECK_LE( theta[i], Real(M_PI/2) );

		if( i > 0 )
			BOOST_CHECK_LE( theta[i-1], theta[i] );
	}


	// reconstruct A, B from GSVD
	auto ds = build_diagonals_like(Number{}, m, p, r, theta);
	auto& D1 = ds.first;
	auto& D2 = ds.second;

	Matrix almost_A = reconstruct_matrix(U1, D1, R, Qt);
	Matrix almost_B = reconstruct_matrix(U2, D2, R, Qt);

	auto frob_A = ublas::norm_frobenius(A);
	auto frob_B = ublas::norm_frobenius(B);

	// The tolerance here is based on the backward error bounds for the QR
	// factorization given in Theorem 19.4, Equation (3.8) in
	// Higham: "Accuracy and Stability of Numerical Algorithms".
	BOOST_CHECK_LE(
		ublas::norm_frobenius(A - almost_A), 10 * (m+p) * n * frob_A * eps );
	BOOST_CHECK_LE(
		ublas::norm_frobenius(w*B - almost_B), 10*w * (m+p) * n * frob_B * eps );
}


/**
 * This structure hides the differences between real and complex implementations
 * of xGGQRCS.
 */
template<typename Real>
struct QrCsCaller
{
	using Number = Real;
	using Matrix = ublas::matrix<Number, ublas::column_major>;
	template<typename U> using Vector = ublas::vector<U>;

	Real w = not_a_number<Real>::value;
	Integer rank = -1;
	std::size_t m, n, p;
	std::size_t ldx, ldy, ldu1, ldu2, ldqt;
	Matrix X, Y;
	Matrix U1, U2, Qt;
	Vector<Number> theta;
	Vector<Number> work;
	Vector<Integer> iwork;


	QrCsCaller(std::size_t m_, std::size_t n_, std::size_t p_)
		: QrCsCaller(m_, n_, p_, m_, p_, m_, p_, n_)
	{}


	QrCsCaller(
		std::size_t m_, std::size_t n_, std::size_t p_,
		std::size_t ldx_, std::size_t ldy_,
		std::size_t ldu1_, std::size_t ldu2_, std::size_t ldqt_
	) :
		m(m_),
		n(n_),
		p(p_),
		ldx(ldx_), ldy(ldy_),
		ldu1(ldu1_), ldu2(ldu2_), ldqt(ldqt_),
		X(ldx, n, 0),
		Y(ldy, n, 0),
		U1(ldu1, m, not_a_number<Number>::value),
		U2(ldu2, p, not_a_number<Number>::value),
		Qt(ldqt, n, not_a_number<Number>::value),
		theta(n, not_a_number<Real>::value),
		iwork(m + n + p, -1)
	{
		BOOST_VERIFY( m > 0 );
		BOOST_VERIFY( n > 0 );
		BOOST_VERIFY( p > 0 );
		BOOST_VERIFY( ldx >= m );
		BOOST_VERIFY( ldy >= p );
		BOOST_VERIFY( ldu1 >= m );
		BOOST_VERIFY( ldu2 >= p );
		BOOST_VERIFY( ldqt >= n );

		auto nan = not_a_number<Number>::value;

		// query workspace size
		auto lwork_opt_f = nan;
		auto w = nan;
		auto rank = Integer{-1};
		auto ret = lapack::ggqrcs(
			'Y', 'Y', 'Y', m, n, p, &w, &rank,
			&X(0, 0), ldx, &Y(0, 0), ldy,
			&theta(0),
			&U1(0, 0), ldu1, &U2(0, 0), ldu2, &Qt(0, 0), ldqt,
			&lwork_opt_f, -1, &iwork(0) );
		BOOST_REQUIRE_EQUAL( ret, 0 );

		// resize workspace accordingly
		auto lwork_opt =
			static_cast<std::size_t>(std::real(lwork_opt_f));

		work.resize( lwork_opt );
		std::fill( work.begin(), work.end(), nan );
	}


	Integer operator() ()
	{
		return lapack::ggqrcs(
			'Y', 'Y', 'Y', m, n, p, &w, &rank,
			&X(0, 0), ldx, &Y(0, 0), ldy,
			&theta(0),
			&U1(0, 0), ldu1, &U2(0, 0), ldu2, &Qt(0, 0), ldqt,
			&work(0), work.size(), &iwork(0)
		);
	}
};


template<typename Real>
struct QrCsCaller<std::complex<Real>>
{
	using Number = std::complex<Real>;
	using Matrix = ublas::matrix<Number, ublas::column_major>;
	template<typename U> using Vector = ublas::vector<U>;

	Real w = not_a_number<Real>::value;
	Integer rank = -1;
	std::size_t m, n, p;
	std::size_t ldx, ldy, ldu1, ldu2, ldqt;
	Matrix X, Y;
	Matrix U1, U2, Qt;
	Vector<Real> theta;
	Vector<Number> work;
	Vector<Real> rwork;
	Vector<Integer> iwork;


	QrCsCaller(std::size_t m_, std::size_t n_, std::size_t p_)
		: QrCsCaller(m_, n_, p_, m_, p_, m_, p_, n_)
	{}

	QrCsCaller(
		std::size_t m_, std::size_t n_, std::size_t p_,
		std::size_t ldx_, std::size_t ldy_,
		std::size_t ldu1_, std::size_t ldu2_, std::size_t ldqt_
	) :
		m(m_),
		n(n_),
		p(p_),
		ldx(ldx_), ldy(ldy_),
		ldu1(ldu1_), ldu2(ldu2_), ldqt(ldqt_),
		X(ldx, n, 0),
		Y(ldy, n, 0),
		U1(ldu1, m, not_a_number<Number>::value),
		U2(ldu2, p, not_a_number<Number>::value),
		Qt(ldqt, n, not_a_number<Number>::value),
		theta(n, not_a_number<Real>::value),
		iwork(m + n + p, -1)
	{
		BOOST_VERIFY( m > 0 );
		BOOST_VERIFY( n > 0 );
		BOOST_VERIFY( p > 0 );
		BOOST_VERIFY( ldx >= m );
		BOOST_VERIFY( ldy >= p );
		BOOST_VERIFY( ldu1 >= m );
		BOOST_VERIFY( ldu2 >= p );
		BOOST_VERIFY( ldqt >= n );

		auto nan = not_a_number<Number>::value;
		auto real_nan = not_a_number<Real>::value;

		// query workspace sizes
		auto lwork_opt_f = nan;
		auto lrwork_opt_f = real_nan;
		auto w = real_nan;
		auto rank = Integer{-1};
		auto ret = lapack::ggqrcs(
			'Y', 'Y', 'Y', m, n, p, &w, &rank,
			&X(0, 0), ldx, &Y(0, 0), ldy,
			&theta(0),
			&U1(0, 0), ldu1, &U2(0, 0), ldu2, &Qt(0, 0), ldqt,
			&lwork_opt_f, -1, &lrwork_opt_f, 1, &iwork(0) );
		BOOST_REQUIRE_EQUAL( ret, 0 );

		auto lwork_opt =
			static_cast<std::size_t>(std::real(lwork_opt_f));

		work.resize( lwork_opt );
		std::fill( work.begin(), work.end(), nan );

		ret = lapack::ggqrcs(
			'Y', 'Y', 'Y', m, n, p, &w, &rank,
			&X(0, 0), ldx, &Y(0, 0), ldy,
			&theta(0),
			&U1(0, 0), ldu1, &U2(0, 0), ldu2, &Qt(0, 0), ldqt,
			&work(0), lwork_opt, &lrwork_opt_f, -1, &iwork(0) );
		BOOST_REQUIRE_EQUAL( ret, 0 );

		auto lrwork_opt =
			static_cast<std::size_t>(std::real(lrwork_opt_f));

		rwork.resize( lrwork_opt );
		std::fill( rwork.begin(), rwork.end(), real_nan );
	}


	Integer operator() ()
	{
		return lapack::ggqrcs(
			'Y', 'Y', 'Y', m, n, p, &w, &rank,
			&X(0, 0), ldx, &Y(0, 0), ldy,
			&theta(0),
			&U1(0, 0), ldu1, &U2(0, 0), ldu2, &Qt(0, 0), ldqt,
			&work(0), work.size(),
			&rwork(0), rwork.size(),
			&iwork(0)
		);
	}
};


template<typename Number, class Matrix>
void check_results(
	Integer ret,
	const Matrix& A, const Matrix& B,
	const QrCsCaller<Number> caller)
{
	auto f = [] (const Matrix& A, std::size_t m, std::size_t n) {
		BOOST_VERIFY( A.size1() >= m );
		BOOST_VERIFY( A.size2() == n );
		return Matrix(ublas::subrange(A, 0, m, 0, n));
	};

	auto m = caller.m;
	auto n = caller.n;
	auto p = caller.p;
	auto X = f(caller.X, m, n);
	auto Y = f(caller.Y, p, n);
	auto U1 = f(caller.U1, m, m);
	auto U2 = f(caller.U2, p, p);
	auto Qt = f(caller.Qt, n, n);

	check_results(
		ret,
		A, B,
		caller.w, caller.rank,
		caller.theta,
		U1, U2, Qt,
		X, Y
	);
}



BOOST_AUTO_TEST_CASE_TEMPLATE(ggqrcs_simple_test, Number, test_types)
{
	auto m = std::size_t{2};
	auto n = std::size_t{2};
	auto p = std::size_t{2};
	auto caller = QrCsCaller<Number>(m, n, p);
	auto A = caller.X;
	auto B = caller.Y;

	A(0,0) = 1;
	B(1,1) = 1;

	caller.X = A;
	caller.Y = B;

	auto ret = caller();
	check_results(ret, A, B, caller);

	BOOST_CHECK_EQUAL( caller.w, 1 );
	BOOST_CHECK_EQUAL( caller.rank, 2 );
}


/**
 * This LAPACK xerbla implementation prints an error message but does not
 * terminate the program thereby allowing the calling LAPACK function to return
 * to its caller.
 *
 * @param[in] caller A string WITHOUT ZERO TERMINATOR
 * @param[in] caller_len The length of the string referenced by f_caller
 */
extern "C" void xerbla_(
	const char* caller, int* p_info, std::size_t caller_len)
{
	BOOST_VERIFY( caller != nullptr );
	BOOST_VERIFY( p_info != nullptr );

	// "sz" prefix taken from hungarian notation (zero-terminated string)
	char szCaller[80];
	auto num_bytes_to_copy = std::min(sizeof(szCaller)-1, caller_len);

	std::memset(szCaller, 0, sizeof(szCaller));
	std::strncpy(szCaller, caller, num_bytes_to_copy);
	std::fprintf(
		stderr, "%s: parameter %d has illegal value\n", szCaller, *p_info
	);
}


template<
	typename Number,
	typename std::enable_if<
		std::is_fundamental<Number>::value, int
	>::type* = nullptr
>
void ggqrcs_zero_dimensions_test_impl(Number)
{
	using Real = typename real_from<Number>::type;

	constexpr auto nan = not_a_number<Number>::value;
	constexpr auto real_nan = not_a_number<Real>::value;

	auto w = real_nan;
	auto rank = Integer{-1};
	auto lda = 1;
	auto A = std::vector<Number>(lda*1, nan);
	auto ldb = 1;
	auto B = std::vector<Number>(ldb*1, nan);
	auto theta = std::vector<Real>(1, real_nan);
	auto lwork = 1;
	auto work = std::vector<Number>(lwork, nan);
	auto iwork = std::vector<Integer>(1, -1);
	auto f = [&] (std::size_t m, std::size_t n, std::size_t p) {
		return lapack::ggqrcs(
			'N', 'N', 'N', m, n, p, &w, &rank,
			&A[0], lda, &B[0], ldb,
			&theta[0],
			nullptr, 1, nullptr, 1, nullptr, 1,
			&work[0], lwork, &iwork[0]
		);
	};

	BOOST_CHECK_EQUAL( f(0, 1, 1), -4 );
	BOOST_CHECK_EQUAL( f(1, 0, 1), -5 );
	BOOST_CHECK_EQUAL( f(1, 1, 0), -6 );
}

template<
	typename Number,
	typename std::enable_if<
		!std::is_fundamental<Number>::value, int
	>::type* = nullptr
>
void ggqrcs_zero_dimensions_test_impl(Number)
{
	using Real = typename real_from<Number>::type;

	constexpr auto nan = not_a_number<Number>::value;
	constexpr auto real_nan = not_a_number<Real>::value;

	auto w = real_nan;
	auto rank = Integer{-1};
	auto lda = 1;
	auto A = std::vector<Number>(lda*1, nan);
	auto ldb = 1;
	auto B = std::vector<Number>(ldb*1, nan);
	auto theta = std::vector<Real>(1, real_nan);
	auto lwork = 1;
	auto work = std::vector<Number>(lwork, nan);
	auto lrwork = 1;
	auto rwork = std::vector<Real>(lrwork, real_nan);
	auto iwork = std::vector<Integer>(1, -1);
	auto f = [&] (std::size_t m, std::size_t n, std::size_t p) {
		return lapack::ggqrcs(
			'N', 'N', 'N', m, n, p, &w, &rank,
			&A[0], lda, &B[0], ldb,
			&theta[0],
			nullptr, 1, nullptr, 1, nullptr, 1,
			&work[0], lwork, &rwork[0], lrwork, &iwork[0]
		);
	};

	BOOST_CHECK_EQUAL( f(0, 1, 1), -4 );
	BOOST_CHECK_EQUAL( f(1, 0, 1), -5 );
	BOOST_CHECK_EQUAL( f(1, 1, 0), -6 );
}

BOOST_AUTO_TEST_CASE_TEMPLATE(ggqrcs_zero_dimensions_test, Number, test_types)
{
	ggqrcs_zero_dimensions_test_impl(Number{0});
}



BOOST_AUTO_TEST_CASE_TEMPLATE(ggqrcs_zero_test, Number, test_types)
{
	auto m = std::size_t{4};
	auto n = std::size_t{3};
	auto p = std::size_t{2};
	auto caller = QrCsCaller<Number>(m, n, p);
	auto A = caller.X;
	auto B = caller.Y;

	auto ret = caller();
	check_results(ret, A, B, caller);

	BOOST_CHECK_EQUAL( caller.w, 1 );
	BOOST_CHECK_EQUAL( caller.rank, 0 );
}


BOOST_AUTO_TEST_CASE_TEMPLATE(ggqrcs_rectangular_test, Number, test_types)
{
	for(std::size_t m : { 2, 13, 41 })
	{
		for(std::size_t n : {3, 7, 31})
		{
			for(std::size_t p : {5, 11, 17})
			{
				auto caller = QrCsCaller<Number>(m, n, p);
				auto A = caller.X;
				auto B = caller.Y;

				A(0,0) = 1;
				A(1,0) = 1;
				B(0,1) = 1;
				B(1,1) = 1;

				caller.X = A;
				caller.Y = B;

				auto ret = caller();
				check_results(ret, A, B, caller);

				BOOST_CHECK_EQUAL( caller.w, 1 );
			}
		}
	}
}



/**
 * @return Matrix A with m rows, n columns, and A^* A = I.
 */
template<
	typename Number,
	class Engine,
	class Storage = ublas::column_major
>
ublas::matrix<Number, Storage> make_isometric_matrix_like(
	Number, std::size_t m, std::size_t n, Engine* gen)
{
	BOOST_VERIFY( m >= n );

	using Matrix = ublas::matrix<Number, Storage>;

	auto p = std::min(m, n);

	if(p <= 1)
		return ublas::identity_matrix<Number, Storage>(m, n);

	auto dist = UniformDistribution<Number>();
	auto rand = [gen, &dist] () { return dist(*gen); };
	auto A = Matrix(m, n);

	std::generate( A.data().begin(), A.data().end(), rand );

	constexpr auto nan = not_a_number<Number>::value;
	auto tau = ublas::vector<Number>(n, nan);
	auto k = std::max(m, n);
	auto lwork = static_cast<Integer>(2 * k * k);
	auto work = ublas::vector<Number>(lwork, nan);
	auto ret = lapack::geqrf(m, n, &A(0,0), m, &tau(0), &work(0), lwork);

	BOOST_VERIFY( ret == 0 );

	ret = lapack::ungqr(m, n, n, &A(0,0), m, &tau(0), &work(0), lwork);

	BOOST_VERIFY( ret == 0 );
	BOOST_VERIFY( is_almost_isometric(A) );

	return A;
}

/**
 * @return A random m x n matrix with spectral condition number cond2.
 */
template<
	typename Number,
	class Engine,
	class Storage = ublas::column_major,
	typename Real = typename real_from<Number>::type
>
ublas::matrix<Number, Storage> make_matrix_like(
	Number dummy, std::size_t m, std::size_t n, Real cond2, Engine* gen)
{
	using Matrix = ublas::matrix<Number, Storage>;
	using BandedMatrix = ublas::banded_matrix<Number>;

	auto p = std::min(m, n);

	if(p <= 1)
		return ublas::identity_matrix<Number, Storage>(m, n);

	auto S = BandedMatrix(p, p);
	// do not sort singular values
	S(0,0) = 1;
	S(1,1) = cond2;
	auto sv_dist = std::uniform_real_distribution<Real>(1, cond2);

	for(auto i = std::size_t{2}; i < p; ++i)
		S(i,i) = sv_dist(*gen);

	auto U = make_isometric_matrix_like(dummy, m, p, gen);
	auto V = make_isometric_matrix_like(dummy, n, p, gen);

	auto US = Matrix(ublas::prod(U, S));
	auto A = Matrix(ublas::prod(US, ublas::trans(V)));

	return A;
}



template<typename Number>
void ggqrcs_random_test_impl(
	Number dummy,
	std::size_t m, std::size_t n, std::size_t p, std::size_t r,
	std::uint64_t seed)
{
	using Real = typename real_from<Number>::type;
	using Matrix = ublas::matrix<Number, ublas::column_major>;

	constexpr auto real_nan = not_a_number<Real>::value;

	BOOST_TEST_CONTEXT("m=" << m) {
	BOOST_TEST_CONTEXT("n=" << n) {
	BOOST_TEST_CONTEXT("p=" << p) {
	BOOST_TEST_CONTEXT("rank=" << r) {
	BOOST_TEST_CONTEXT("seed=" << seed) {

	auto gen = std::mt19937(seed);

	gen.discard(1u << 17);

	auto k = std::min( {m, p, r, m + p - r} );
	auto theta_dist = std::uniform_real_distribution<Real>(0, M_PI/2);
	auto theta = ublas::vector<Real>(k, real_nan);

	std::generate(
		theta.begin(), theta.end(),
		[&gen, &theta_dist](){ return theta_dist(gen); }
	);

	auto min_log_cond_R = Real{0};
	auto max_log_cond_R = static_cast<Real>(std::numeric_limits<Real>::digits);
	auto log_cond_dist =
		std::uniform_real_distribution<Real>(min_log_cond_R, max_log_cond_R);
	auto log_cond_R = log_cond_dist(gen);
	auto cond_R = std::pow(Real{2}, log_cond_R);
	auto R_Qt = make_matrix_like(dummy, r, n, cond_R, &gen);
	auto U1 = make_isometric_matrix_like(dummy, m, m, &gen);
	auto U2 = make_isometric_matrix_like(dummy, p, p, &gen);
	auto ds = build_diagonals_like(dummy, m, p, r, theta);
	auto D1 = ds.first;
	auto D2 = ds.second;
	auto A = Matrix(ublas::prod(Matrix(ublas::prod(U1, D1)), R_Qt));
	auto B = Matrix(ublas::prod(Matrix(ublas::prod(U2, D2)), R_Qt));

	// initialize caller
	auto ldx = m + 11;
	auto ldy = p + 5;
	auto ldu1 = m + 13;
	auto ldu2 = p + 7;
	auto ldqt = n + 17;
	auto caller = QrCsCaller<Number>(m, n, p, ldx, ldy, ldu1, ldu2, ldqt);

	ublas::subrange(caller.X, 0, m, 0, n) = A;
	ublas::subrange(caller.Y, 0, p, 0, n) = B;

	auto ret = caller();

	check_results(ret, A, B, caller);

	BOOST_CHECK_LE( caller.rank, r );
}
}
}
}
}
}


BOOST_AUTO_TEST_CASE_TEMPLATE(ggqrcs_random_test, Number, test_types)
{
	constexpr std::size_t dimensions[] = { 1, 2, 3, 4, 10, 20 };

	auto gen = std::mt19937();
	auto seed_dist = std::uniform_int_distribution<std::uint64_t>(0);

	gen.discard(1u << 17);

	for(auto m : dimensions)
	{
		for(auto n : dimensions)
		{
			for(auto p : dimensions)
			{
				auto max_rank = std::min(m+p, n);
				for(auto rank = std::size_t{0}; rank <= max_rank; ++rank)
				{
					for(auto iteration = 0u; iteration < 10u; ++iteration)
					{
						auto seed = seed_dist(gen);

						ggqrcs_random_test_impl(Number{0}, m, n, p, rank, seed);
					}
				}
			}
		}
	}
}

BOOST_TEST_DECORATOR(* boost::unit_test::disabled())
BOOST_AUTO_TEST_CASE_TEMPLATE(
	infinite_ggqrcs_random_test, Number, test_types)
{
	constexpr auto min_dimension = std::size_t{1};
	constexpr auto max_dimension = std::size_t{1000};

	auto master_seed = std::uintmax_t(std::time(nullptr));

	std::printf("infinite_ggqrcs_random_test master_seed=%ju\n", master_seed);

	auto gen = std::mt19937(master_seed);
	auto dim_dist =
		std::uniform_int_distribution<std::size_t>(min_dimension,max_dimension);
	auto seed_dist = std::uniform_int_distribution<std::uint64_t>(0);

	gen.discard(1u << 17);

	auto start_time_sec = std::time(nullptr);
	auto last_time_sec = std::time_t{0};
	auto iteration = std::uintmax_t{0};

	constexpr char FMT[] = "%7jd %13ju  %3zu %3zu %3zu %4zu  %20zu\n";
	std::printf(
		"%7s %13s  %3s %3s %3s %4s  %20s\n",
		"time(s)", "iteration", "m", "n", "p", "rank", "seed"
	);

	while(true)
	{
		auto m = dim_dist(gen);
		auto n = dim_dist(gen);
		auto p = dim_dist(gen);
		auto max_rank = (m + p <= n) ? m + p : n;

		for(auto rank = std::size_t{0}; rank <= max_rank; ++rank, ++iteration)
		{
			auto seed = seed_dist(gen);
			auto now_sec = std::time(nullptr);
			auto second = std::time_t{1};

			if(last_time_sec + 60*second < now_sec)
			{
				auto time_passed_sec = std::intmax_t{now_sec - start_time_sec};

				std::printf(
					FMT, time_passed_sec, iteration, m, n, p, rank, seed
				);

				last_time_sec = now_sec;
			}

			ggqrcs_random_test_impl(Number{0}, m, n, p, rank, seed);
		}
	}
}



template<
	typename Number,
	typename std::enable_if<
		std::is_fundamental<Number>::value, int
	>::type* = nullptr
>
void uncsd2by1_regression_20200420_impl(Number)
{
	using Real = typename real_from<Number>::type;
	using Matrix = ublas::matrix<Number, ublas::column_major>;

	auto m = std::size_t{130}; // 160
	auto n = std::size_t{131}; // 179
	auto p = std::size_t{130}; // 150
	auto ldq = m + p;
	auto gen = std::minstd_rand();
	auto Q = make_isometric_matrix_like(Real{0}, ldq, n, &gen);

	BOOST_VERIFY( is_almost_isometric(Q) );

	auto nan = not_a_number<Number>::value;
	auto real_nan = not_a_number<Real>::value;
	auto theta = ublas::vector<Real>(n, real_nan);
	auto U1 = Matrix(m, m, nan);
	auto U2 = Matrix(p, p, nan);
	auto Qt = Matrix(n, n, nan);
	auto lwork = 32 * ldq;
	auto work = ublas::vector<Number>(lwork, nan);
	auto iwork = ublas::vector<Integer>(m+p, -1);
	auto ret = lapack::uncsd2by1(
		'Y', 'Y', 'N',
		m + p, m, n,
		&Q(0, 0), ldq, &Q(m, 0), ldq,
		&theta(0),
		&U1(0,0), m, &U2(0,0), p, &Qt(0,0), n,
		&work(0), lwork, &iwork(0)
	);

	BOOST_REQUIRE_EQUAL( ret, 0 );
	BOOST_CHECK( is_almost_isometric(U1) );
	BOOST_CHECK( is_almost_isometric(U2) );
}

template<
	typename Number,
	typename std::enable_if<
		!std::is_fundamental<Number>::value, int
	>::type* = nullptr
>
void uncsd2by1_regression_20200420_impl(Number)
{
	using Real = typename real_from<Number>::type;
	using Matrix = ublas::matrix<Number, ublas::column_major>;

	auto m = std::size_t{130}; // 160
	auto n = std::size_t{131}; // 179
	auto p = std::size_t{130}; // 150
	auto ldq = m + p;
	auto gen = std::minstd_rand();
	auto Q = make_isometric_matrix_like(Number{0}, ldq, n, &gen);

	BOOST_VERIFY( is_almost_isometric(Q) );

	auto nan = not_a_number<Number>::value;
	auto real_nan = not_a_number<Real>::value;
	auto theta = ublas::vector<Real>(n, real_nan);
	auto U1 = Matrix(m, m, nan);
	auto U2 = Matrix(p, p, nan);
	auto Qt = Matrix(n, n, nan);
	auto lwork = 32 * ldq;
	auto work = ublas::vector<Number>(lwork, nan);
	auto lrwork = 32 * ldq;
	auto rwork = ublas::vector<Real>(lrwork, real_nan);
	auto iwork = ublas::vector<Integer>(m+p, -1);
	auto ret = lapack::uncsd2by1(
		'Y', 'Y', 'N',
		m + p, m, n,
		&Q(0, 0), ldq, &Q(m, 0), ldq,
		&theta(0),
		&U1(0,0), m, &U2(0,0), p, &Qt(0,0), n,
		&work(0), lwork, &rwork(0), lrwork, &iwork(0)
	);

	BOOST_REQUIRE_EQUAL( ret, 0 );
	BOOST_CHECK( is_almost_isometric(U1) );
	BOOST_CHECK( is_almost_isometric(U2) );
}

BOOST_AUTO_TEST_CASE_TEMPLATE(uncsd2by1_regression_20200420, Number, test_types)
{
	uncsd2by1_regression_20200420_impl(Number{});
}

#endif
