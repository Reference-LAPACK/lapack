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
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <lapack.hpp>

#include <limits>


namespace ublas = boost::numeric::ublas;
typedef lapack::integer_t Integer;


template<typename T>
struct Fixture
{
	typedef ublas::matrix<T, ublas::column_major> Matrix;
	template<typename U> using Vector = ublas::vector<U>;

	Fixture(std::size_t m, std::size_t p, std::size_t n) :
		A(m, n),
		B(p, n),
		U1(m, m),
		U2(p, p),
		Qt(n, n),
		theta(n),
		iwork(m + n + p)
	{
		const T nan = std::numeric_limits<T>::quiet_NaN();

		// query workspace size
		T lwork_opt_f = nan;
		T w = nan;
		Integer l = -1;

		Integer ret = lapack::ggqrcs(
			'Y', 'Y', 'Y', m, n, p, &w, &l,
			&A(0, 0), m, &B(0, 0), p,
			&theta(0),
			&U1(0, 0), m, &U2(0, 0), p, &Qt(0, 0), n,
			&lwork_opt_f, -1, &iwork(0) );
		assert( ret == 0 );

		std::size_t lwork_opt = lwork_opt_f;
		work.resize( lwork_opt );

		// initialize
		std::fill( A.data().begin(), A.data().end(), 0 );
		std::fill( B.data().begin(), B.data().end(), 0 );
		std::fill( U1.data().begin(), U1.data().end(), nan );
		std::fill( U2.data().begin(), U2.data().end(), nan );
		std::fill( Qt.data().begin(), Qt.data().end(), nan );
		std::fill( work.begin(), work.end(), nan );
		std::fill( iwork.begin(), iwork.end(), -1 );
	}

	Matrix A, B;
	Matrix U1, U2, Qt;
	Vector<T> theta;
	Vector<T> work;
	Vector<Integer> iwork;
};


BOOST_AUTO_TEST_SUITE(LAPACK_TEST_SUITE_NAME)

BOOST_AUTO_TEST_CASE_TEMPLATE(
	ggqrcs_simple_test, T, test_types)
{
	const size_t m = 2;
	const size_t n = 2;
	const size_t p = 2;

	Fixture<T> fixture(m, n, p);

	auto& A = fixture.A;
	auto& B = fixture.B;

	A(0,0) = 1;
	B(1,1) = 1;

	auto& theta = fixture.theta;
	auto& U1 = fixture.U1;
	auto& U2 = fixture.U2;
	auto& Qt = fixture.Qt;
	auto& work = fixture.work;
	auto& iwork = fixture.iwork;

	const Integer lwork = work.size();
	double w = -1;
	Integer rank = -1;

	Integer ret = lapack::ggqrcs(
		'Y', 'Y', 'Y', m, n, p, &w, &rank,
		&A(0, 0), m, &B(0, 0), p,
		&theta(0),
		&U1(0, 0), m, &U2(0, 0), p, &Qt(0, 0), n,
		&work(0), lwork, &iwork(0) );

	BOOST_CHECK_EQUAL( ret, 0 );
	BOOST_CHECK_EQUAL( w, 1 );
	BOOST_CHECK_EQUAL( rank, 2 );
}

BOOST_AUTO_TEST_SUITE_END()

#endif
