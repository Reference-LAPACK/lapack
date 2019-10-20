#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <boost/timer/timer.hpp>

#define HAS_GGSVD3
#include <lapack.hpp>

#include <random>

#include <cstdio>
#include <cstdlib>

#include <cassert>



namespace ublas = boost::numeric::ublas;
typedef lapack::integer_t Integer;


template<typename U> using Vector = ublas::vector<U>;


template<typename T, class Storage>
std::pair< ublas::matrix<T, Storage>, ublas::matrix<T, Storage> >
make_matrices(std::size_t m, std::size_t n, std::size_t p, unsigned seed)
{
	typedef ublas::matrix<T, Storage> Matrix;

	const T nan = std::numeric_limits<T>::quiet_NaN();


	std::mt19937_64 gen( seed );
	std::uniform_real_distribution<T> dist(-1, +1);
	auto make_rand = [&gen, &dist] (T s)
	{
		auto rand = [&gen, &dist, s] () { return s*dist(gen); };
		return rand;
	};


	Matrix A(m, n, nan);
	std::generate( A.data().begin(), A.data().end(), make_rand(1000) );

	Matrix B(p, n, nan);
	std::generate( B.data().begin(), B.data().end(), make_rand(1) );

	return std::make_pair(A, B);
}



enum class Solver { Direct, QR_CS };


template<typename T, class Storage>
void compute_gsvd(
	ublas::matrix<T, Storage> A, ublas::matrix<T, Storage> B,
	Solver solver)
{
	assert( A.size2() == B.size2() );

	typedef ublas::matrix<T, Storage> Matrix;

	const std::size_t m = A.size1();
	const std::size_t n = A.size2();
	const std::size_t p = B.size1();
	const T nan = std::numeric_limits<T>::quiet_NaN();

	Matrix U1(m, m, nan);
	Matrix U2(p, p, nan);
	Vector<Integer> iwork(m + n + p, -1);


	if( solver == Solver::Direct )
	{
		Integer k = -1;
		Integer l = -1;
		T lwork_opt_f = nan;
		Vector<T> alpha(n, nan);
		Vector<T> beta(n, nan);
		Matrix Q(n, n, nan);

		Integer ret = lapack::ggsvd3(
			'U', 'V', 'Q', m, n, p, &k, &l,
			&A(0, 0), m, &B(0, 0), p,
			&alpha(0), &beta(0),
			&U1(0, 0), m, &U2(0, 0), p, &Q(0, 0), n,
			&lwork_opt_f, -1, &iwork(0) );
		assert( ret == 0 ); (void) ret;

		assert( !std::isnan(lwork_opt_f) );
		assert( std::isfinite(lwork_opt_f) );
		assert( lwork_opt_f > 0 );

		Integer lwork = lwork_opt_f;
		Vector<T> work( lwork, nan );

		boost::timer::auto_cpu_timer t( "gsvd %ws real %us user %ss sys\n" );

		ret = lapack::ggsvd3(
			'U', 'V', 'Q', m, n, p, &k, &l,
			&A(0, 0), m, &B(0, 0), p,
			&alpha(0), &beta(0),
			&U1(0, 0), m, &U2(0, 0), p, &Q(0, 0), n,
			&work(0), lwork, &iwork(0) );
		assert( ret == 0 );
	}
	else if( solver == Solver::QR_CS )
	{
		T lwork_opt_f = nan;
		T w = nan;
		Integer l = -1;
		Vector<T> theta(n, nan);
		Matrix Qt(n, n, nan);

		Integer ret = lapack::ggqrcs(
			'Y', 'Y', 'Y', m, n, p, &w, &l,
			&A(0, 0), m, &B(0, 0), p,
			&theta(0),
			&U1(0, 0), m, &U2(0, 0), p, &Qt(0, 0), n,
			&lwork_opt_f, -1, &iwork(0) );
		assert( ret == 0 ); (void) ret;

		assert( !std::isnan(lwork_opt_f) );
		assert( std::isfinite(lwork_opt_f) );
		assert( lwork_opt_f > 0 );

		Integer lwork = lwork_opt_f;
		Vector<T> work( lwork, nan );

		boost::timer::auto_cpu_timer t( "qrcs %ws real %us user %ss sys\n" );

		ret = lapack::ggqrcs(
			'Y', 'Y', 'Y', m, n, p, &w, &l,
			&A(0, 0), m, &B(0, 0), p,
			&theta(0),
			&U1(0, 0), m, &U2(0, 0), p, &Qt(0, 0), n,
			&work(0), lwork, &iwork(0) );
		assert( ret == 0 );
	}
}



int main(int argc, const char** argv)
{
	if( argc != 4 )
	{
		const char* aout = argc > 0 ? argv[0] : "a.out";
		std::fprintf( stderr, "usage: %s <m> <n> <p>\n", aout );
		return 1;
	}

	const std::size_t m = std::atoi( argv[1] );
	const std::size_t n = std::atoi( argv[2] );
	const std::size_t p = std::atoi( argv[3] );

	if( m <= 0 || n <= 0 || p <= 0 )
	{
		std::fprintf(
			stderr,
			"m=%zu, n=%zu, p=%zu must be positive\n",
			m, n, p );
		return 2;
	}

	auto matrices = make_matrices<double, ublas::column_major>(m, n, p, 1u);
	const auto& A = matrices.first;
	const auto& B = matrices.second;

	compute_gsvd( A, B, Solver::QR_CS );
	compute_gsvd( A, B, Solver::Direct );
}
