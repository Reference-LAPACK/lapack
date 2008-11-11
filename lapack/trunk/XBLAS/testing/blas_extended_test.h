#ifndef BLAS_EXTENDED_TEST_H
#define BLAS_EXTENDED_TEST_H

#define MAX_BAD_TESTS 100
#define TOTAL_FAILURE_THRESHOLD 1000

double power(int i1, int i2);
double xrand(int *is);
int FixedBits(double r_true_l, double r_true_t);

void ddmuld(double dda_l, double dda_t, double db, double *ddc_l,
	    double *ddc_t);
void ddadd(double dda_l, double dda_t, double ddb_l, double ddb_t,
	   double *ddc_l, double *ddc_t);
void dddiv(double dda_l, double dda_t, double ddb_l, double ddb_t,
	   double *ddc_l, double *ddc_t);
void z_dddivd(double *dda_l, double *dda_t, double *db, double *ddc_l,
	      double *ddc_t);

void testgen_BLAS_sdot(int n, int n_fix2, int n_mix, int norm,
		       enum blas_conj_type conj, float *alpha, int alpha_flag,
		       float *beta, int beta_flag, float *x, float *y,
		       int *seed, float *r, double *r_true_l,
		       double *r_true_t);

void testgen_BLAS_cdot(int n, int n_fix2, int n_mix, int norm,
		       enum blas_conj_type conj, void *alpha, int alpha_flag,
		       void *beta, int beta_flag, void *x, void *y, int *seed,
		       void *r, double r_true_l[], double r_true_t[]);

void testgen_BLAS_ddot(int n, int n_fix2, int n_mix, int norm,
		       enum blas_conj_type conj, double *alpha,
		       int alpha_flag, double *beta, int beta_flag, double *x,
		       double *y, int *seed, double *r, double *r_true_l,
		       double *r_true_t);

void testgen_BLAS_zdot(int n, int n_fix2, int n_mix, int norm,
		       enum blas_conj_type conj, void *alpha, int alpha_flag,
		       void *beta, int beta_flag, void *x, void *y, int *seed,
		       void *r, double r_true_l[], double r_true_t[]);

void testgen_BLAS_sdot2(int n, int n_fix2, int n_mix, int norm,
			enum blas_conj_type conj, float *alpha,
			int alpha_flag, float *beta, int beta_flag,
			float *head_x, float *tail_x, float *y, int *seed,
			float *r, double *r_true_l, double *r_true_t);

void testgen_BLAS_cdot2(int n, int n_fix2, int n_mix, int norm,
			enum blas_conj_type conj, void *alpha, int alpha_flag,
			void *beta, int beta_flag, void *head_x, void *tail_x,
			void *y, int *seed, void *r, double r_true_l[],
			double r_true_t[]);

void testgen_BLAS_ddot2(int n, int n_fix2, int n_mix, int norm,
			enum blas_conj_type conj, double *alpha,
			int alpha_flag, double *beta, int beta_flag,
			double *head_x, double *tail_x, double *y, int *seed,
			double *r, double *r_true_l, double *r_true_t);

void testgen_BLAS_zdot2(int n, int n_fix2, int n_mix, int norm,
			enum blas_conj_type conj, void *alpha, int alpha_flag,
			void *beta, int beta_flag, void *head_x, void *tail_x,
			void *y, int *seed, void *r, double r_true_l[],
			double r_true_t[]);


void BLAS_sdot_testgen(int n, int n_fix2, int n_mix, int norm,
		       enum blas_conj_type conj, float *alpha, int alpha_flag,
		       float *beta, int beta_flag, float *x, float *y,
		       int *seed, float *r, double *r_true_l,
		       double *r_true_t);
void BLAS_ddot_testgen(int n, int n_fix2, int n_mix, int norm,
		       enum blas_conj_type conj, double *alpha,
		       int alpha_flag, double *beta, int beta_flag, double *x,
		       double *y, int *seed, double *r, double *r_true_l,
		       double *r_true_t);
void BLAS_cdot_testgen(int n, int n_fix2, int n_mix, int norm,
		       enum blas_conj_type conj, void *alpha, int alpha_flag,
		       void *beta, int beta_flag, void *x, void *y, int *seed,
		       void *r, double *r_true_l, double *r_true_t);
void BLAS_zdot_testgen(int n, int n_fix2, int n_mix, int norm,
		       enum blas_conj_type conj, void *alpha, int alpha_flag,
		       void *beta, int beta_flag, void *x, void *y, int *seed,
		       void *r, double *r_true_l, double *r_true_t);
void BLAS_cdot_s_s_testgen(int n, int n_fix2, int n_mix, int norm,
			   enum blas_conj_type conj, void *alpha,
			   int alpha_flag, void *beta, int beta_flag,
			   float *x, float *y, int *seed, void *r,
			   double *r_true_l, double *r_true_t);
void BLAS_cdot_s_c_testgen(int n, int n_fix2, int n_mix, int norm,
			   enum blas_conj_type conj, void *alpha,
			   int alpha_flag, void *beta, int beta_flag,
			   float *x, void *y, int *seed, void *r,
			   double *r_true_l, double *r_true_t);
void BLAS_cdot_c_s_testgen(int n, int n_fix2, int n_mix, int norm,
			   enum blas_conj_type conj, void *alpha,
			   int alpha_flag, void *beta, int beta_flag, void *x,
			   float *y, int *seed, void *r, double *r_true_l,
			   double *r_true_t);
void BLAS_zdot_d_d_testgen(int n, int n_fix2, int n_mix, int norm,
			   enum blas_conj_type conj, void *alpha,
			   int alpha_flag, void *beta, int beta_flag,
			   double *x, double *y, int *seed, void *r,
			   double *r_true_l, double *r_true_t);
void BLAS_zdot_d_z_testgen(int n, int n_fix2, int n_mix, int norm,
			   enum blas_conj_type conj, void *alpha,
			   int alpha_flag, void *beta, int beta_flag,
			   double *x, void *y, int *seed, void *r,
			   double *r_true_l, double *r_true_t);
void BLAS_zdot_z_d_testgen(int n, int n_fix2, int n_mix, int norm,
			   enum blas_conj_type conj, void *alpha,
			   int alpha_flag, void *beta, int beta_flag, void *x,
			   double *y, int *seed, void *r, double *r_true_l,
			   double *r_true_t);
void BLAS_ddot_s_s_testgen(int n, int n_fix2, int n_mix, int norm,
			   enum blas_conj_type conj, double *alpha,
			   int alpha_flag, double *beta, int beta_flag,
			   float *x, float *y, int *seed, double *r,
			   double *r_true_l, double *r_true_t);
void BLAS_ddot_s_d_testgen(int n, int n_fix2, int n_mix, int norm,
			   enum blas_conj_type conj, double *alpha,
			   int alpha_flag, double *beta, int beta_flag,
			   float *x, double *y, int *seed, double *r,
			   double *r_true_l, double *r_true_t);
void BLAS_ddot_d_s_testgen(int n, int n_fix2, int n_mix, int norm,
			   enum blas_conj_type conj, double *alpha,
			   int alpha_flag, double *beta, int beta_flag,
			   double *x, float *y, int *seed, double *r,
			   double *r_true_l, double *r_true_t);
void BLAS_zdot_c_c_testgen(int n, int n_fix2, int n_mix, int norm,
			   enum blas_conj_type conj, void *alpha,
			   int alpha_flag, void *beta, int beta_flag, void *x,
			   void *y, int *seed, void *r, double *r_true_l,
			   double *r_true_t);
void BLAS_zdot_c_z_testgen(int n, int n_fix2, int n_mix, int norm,
			   enum blas_conj_type conj, void *alpha,
			   int alpha_flag, void *beta, int beta_flag, void *x,
			   void *y, int *seed, void *r, double *r_true_l,
			   double *r_true_t);
void BLAS_zdot_z_c_testgen(int n, int n_fix2, int n_mix, int norm,
			   enum blas_conj_type conj, void *alpha,
			   int alpha_flag, void *beta, int beta_flag, void *x,
			   void *y, int *seed, void *r, double *r_true_l,
			   double *r_true_t);
void scopy_vector(const float *x, int n, int incx, float *y, int incy);
void dcopy_vector(const double *x, int n, int incx, double *y, int incy);
void ccopy_vector(const void *x, int n, int incx, void *y, int incy);
void zcopy_vector(const void *x, int n, int incx, void *y, int incy);

void sprint_vector(const float *x, int n, int inc, const char *name);
void dprint_vector(const double *x, int n, int inc, const char *name);
void cprint_vector(const void *x, int n, int inc, const char *name);
void zprint_vector(const void *x, int n, int inc, const char *name);

void test_BLAS_sdot(int n, enum blas_conj_type conj, float alpha, float beta,
		    float rin, float rout, double r_true_l, double r_true_t,
		    float *x, int incx, float *y, int incy, double eps_int,
		    double un_int, double *test_ratio);
void test_BLAS_ddot(int n, enum blas_conj_type conj, double alpha,
		    double beta, double rin, double rout, double r_true_l,
		    double r_true_t, double *x, int incx, double *y, int incy,
		    double eps_int, double un_int, double *test_ratio);
void test_BLAS_cdot(int n, enum blas_conj_type conj, const void *alpha,
		    const void *beta, const void *rin, const void *rout,
		    double *r_true_l, double *r_true_t, void *x, int incx,
		    void *y, int incy, double eps_int, double un_int,
		    double *test_ratio);
void test_BLAS_zdot(int n, enum blas_conj_type conj, const void *alpha,
		    const void *beta, const void *rin, const void *rout,
		    double *r_true_l, double *r_true_t, void *x, int incx,
		    void *y, int incy, double eps_int, double un_int,
		    double *test_ratio);
void test_BLAS_cdot_s_s(int n, enum blas_conj_type conj, const void *alpha,
			const void *beta, const void *rin, const void *rout,
			double *r_true_l, double *r_true_t, float *x,
			int incx, float *y, int incy, double eps_int,
			double un_int, double *test_ratio);
void test_BLAS_cdot_s_c(int n, enum blas_conj_type conj, const void *alpha,
			const void *beta, const void *rin, const void *rout,
			double *r_true_l, double *r_true_t, float *x,
			int incx, void *y, int incy, double eps_int,
			double un_int, double *test_ratio);
void test_BLAS_cdot_c_s(int n, enum blas_conj_type conj, const void *alpha,
			const void *beta, const void *rin, const void *rout,
			double *r_true_l, double *r_true_t, void *x, int incx,
			float *y, int incy, double eps_int, double un_int,
			double *test_ratio);
void test_BLAS_zdot_d_d(int n, enum blas_conj_type conj, const void *alpha,
			const void *beta, const void *rin, const void *rout,
			double *r_true_l, double *r_true_t, double *x,
			int incx, double *y, int incy, double eps_int,
			double un_int, double *test_ratio);
void test_BLAS_zdot_d_z(int n, enum blas_conj_type conj, const void *alpha,
			const void *beta, const void *rin, const void *rout,
			double *r_true_l, double *r_true_t, double *x,
			int incx, void *y, int incy, double eps_int,
			double un_int, double *test_ratio);
void test_BLAS_zdot_z_d(int n, enum blas_conj_type conj, const void *alpha,
			const void *beta, const void *rin, const void *rout,
			double *r_true_l, double *r_true_t, void *x, int incx,
			double *y, int incy, double eps_int, double un_int,
			double *test_ratio);
void test_BLAS_ddot_s_s(int n, enum blas_conj_type conj, double alpha,
			double beta, double rin, double rout, double r_true_l,
			double r_true_t, float *x, int incx, float *y,
			int incy, double eps_int, double un_int,
			double *test_ratio);
void test_BLAS_ddot_s_d(int n, enum blas_conj_type conj, double alpha,
			double beta, double rin, double rout, double r_true_l,
			double r_true_t, float *x, int incx, double *y,
			int incy, double eps_int, double un_int,
			double *test_ratio);
void test_BLAS_ddot_d_s(int n, enum blas_conj_type conj, double alpha,
			double beta, double rin, double rout, double r_true_l,
			double r_true_t, double *x, int incx, float *y,
			int incy, double eps_int, double un_int,
			double *test_ratio);
void test_BLAS_zdot_c_c(int n, enum blas_conj_type conj, const void *alpha,
			const void *beta, const void *rin, const void *rout,
			double *r_true_l, double *r_true_t, void *x, int incx,
			void *y, int incy, double eps_int, double un_int,
			double *test_ratio);
void test_BLAS_zdot_c_z(int n, enum blas_conj_type conj, const void *alpha,
			const void *beta, const void *rin, const void *rout,
			double *r_true_l, double *r_true_t, void *x, int incx,
			void *y, int incy, double eps_int, double un_int,
			double *test_ratio);
void test_BLAS_zdot_z_c(int n, enum blas_conj_type conj, const void *alpha,
			const void *beta, const void *rin, const void *rout,
			double *r_true_l, double *r_true_t, void *x, int incx,
			void *y, int incy, double eps_int, double un_int,
			double *test_ratio);

void BLAS_sdot2_testgen(int n, int n_fix2, int n_mix, int norm,
			enum blas_conj_type conj, float *alpha,
			int alpha_flag, float *beta, int beta_flag,
			float *head_x, float *tail_x, float *y, int *seed,
			float *r, double *r_true_l, double *r_true_t);
void BLAS_ddot2_testgen(int n, int n_fix2, int n_mix, int norm,
			enum blas_conj_type conj, double *alpha,
			int alpha_flag, double *beta, int beta_flag,
			double *head_x, double *tail_x, double *y, int *seed,
			double *r, double *r_true_l, double *r_true_t);
void BLAS_cdot2_testgen(int n, int n_fix2, int n_mix, int norm,
			enum blas_conj_type conj, void *alpha, int alpha_flag,
			void *beta, int beta_flag, void *head_x, void *tail_x,
			void *y, int *seed, void *r, double *r_true_l,
			double *r_true_t);
void BLAS_zdot2_testgen(int n, int n_fix2, int n_mix, int norm,
			enum blas_conj_type conj, void *alpha, int alpha_flag,
			void *beta, int beta_flag, void *head_x, void *tail_x,
			void *y, int *seed, void *r, double *r_true_l,
			double *r_true_t);
void BLAS_cdot2_s_s_testgen(int n, int n_fix2, int n_mix, int norm,
			    enum blas_conj_type conj, void *alpha,
			    int alpha_flag, void *beta, int beta_flag,
			    float *head_x, float *tail_x, float *y, int *seed,
			    void *r, double *r_true_l, double *r_true_t);
void BLAS_cdot2_s_c_testgen(int n, int n_fix2, int n_mix, int norm,
			    enum blas_conj_type conj, void *alpha,
			    int alpha_flag, void *beta, int beta_flag,
			    float *head_x, float *tail_x, void *y, int *seed,
			    void *r, double *r_true_l, double *r_true_t);
void BLAS_cdot2_c_s_testgen(int n, int n_fix2, int n_mix, int norm,
			    enum blas_conj_type conj, void *alpha,
			    int alpha_flag, void *beta, int beta_flag,
			    void *head_x, void *tail_x, float *y, int *seed,
			    void *r, double *r_true_l, double *r_true_t);
void BLAS_zdot2_d_d_testgen(int n, int n_fix2, int n_mix, int norm,
			    enum blas_conj_type conj, void *alpha,
			    int alpha_flag, void *beta, int beta_flag,
			    double *head_x, double *tail_x, double *y,
			    int *seed, void *r, double *r_true_l,
			    double *r_true_t);
void BLAS_zdot2_d_z_testgen(int n, int n_fix2, int n_mix, int norm,
			    enum blas_conj_type conj, void *alpha,
			    int alpha_flag, void *beta, int beta_flag,
			    double *head_x, double *tail_x, void *y,
			    int *seed, void *r, double *r_true_l,
			    double *r_true_t);
void BLAS_zdot2_z_d_testgen(int n, int n_fix2, int n_mix, int norm,
			    enum blas_conj_type conj, void *alpha,
			    int alpha_flag, void *beta, int beta_flag,
			    void *head_x, void *tail_x, double *y, int *seed,
			    void *r, double *r_true_l, double *r_true_t);
void BLAS_ddot2_s_s_testgen(int n, int n_fix2, int n_mix, int norm,
			    enum blas_conj_type conj, double *alpha,
			    int alpha_flag, double *beta, int beta_flag,
			    float *head_x, float *tail_x, float *y, int *seed,
			    double *r, double *r_true_l, double *r_true_t);
void BLAS_ddot2_s_d_testgen(int n, int n_fix2, int n_mix, int norm,
			    enum blas_conj_type conj, double *alpha,
			    int alpha_flag, double *beta, int beta_flag,
			    float *head_x, float *tail_x, double *y,
			    int *seed, double *r, double *r_true_l,
			    double *r_true_t);
void BLAS_ddot2_d_s_testgen(int n, int n_fix2, int n_mix, int norm,
			    enum blas_conj_type conj, double *alpha,
			    int alpha_flag, double *beta, int beta_flag,
			    double *head_x, double *tail_x, float *y,
			    int *seed, double *r, double *r_true_l,
			    double *r_true_t);
void BLAS_zdot2_c_c_testgen(int n, int n_fix2, int n_mix, int norm,
			    enum blas_conj_type conj, void *alpha,
			    int alpha_flag, void *beta, int beta_flag,
			    void *head_x, void *tail_x, void *y, int *seed,
			    void *r, double *r_true_l, double *r_true_t);
void BLAS_zdot2_c_z_testgen(int n, int n_fix2, int n_mix, int norm,
			    enum blas_conj_type conj, void *alpha,
			    int alpha_flag, void *beta, int beta_flag,
			    void *head_x, void *tail_x, void *y, int *seed,
			    void *r, double *r_true_l, double *r_true_t);
void BLAS_zdot2_z_c_testgen(int n, int n_fix2, int n_mix, int norm,
			    enum blas_conj_type conj, void *alpha,
			    int alpha_flag, void *beta, int beta_flag,
			    void *head_x, void *tail_x, void *y, int *seed,
			    void *r, double *r_true_l, double *r_true_t);
void BLAS_sdot2_x(enum blas_conj_type conj, int n, float alpha,
		  const float *x, int incx, float beta, const float *head_y,
		  const float *tail_y, int incy, float *r,
		  enum blas_prec_type prec);
void BLAS_ddot2_x(enum blas_conj_type conj, int n, double alpha,
		  const double *x, int incx, double beta,
		  const double *head_y, const double *tail_y, int incy,
		  double *r, enum blas_prec_type prec);
void BLAS_cdot2_x(enum blas_conj_type conj, int n, const void *alpha,
		  const void *x, int incx, const void *beta,
		  const void *head_y, const void *tail_y, int incy, void *r,
		  enum blas_prec_type prec);
void BLAS_zdot2_x(enum blas_conj_type conj, int n, const void *alpha,
		  const void *x, int incx, const void *beta,
		  const void *head_y, const void *tail_y, int incy, void *r,
		  enum blas_prec_type prec);
void s_r_truth2(enum blas_conj_type conj, int n, float alpha, const float *x,
		int incx, float beta, const float *head_y,
		const float *tail_y, int incy, float *r, double *head_r_true,
		double *tail_r_true);
void d_r_truth2(enum blas_conj_type conj, int n, double alpha,
		const double *x, int incx, double beta, const double *head_y,
		const double *tail_y, int incy, double *r,
		double *head_r_true, double *tail_r_true);
void c_r_truth2(enum blas_conj_type conj, int n, const void *alpha,
		const void *x, int incx, const void *beta, const void *head_y,
		const void *tail_y, int incy, void *r, double *head_r_true,
		double *tail_r_true);
void z_r_truth2(enum blas_conj_type conj, int n, const void *alpha,
		const void *x, int incx, const void *beta, const void *head_y,
		const void *tail_y, int incy, void *r, double *head_r_true,
		double *tail_r_true);

void test_BLAS_sdot2(int n, enum blas_conj_type conj, float alpha, float beta,
		     float rin, float rout, double r_true_l, double r_true_t,
		     float *x, int incx, float *head_y, float *tail_y,
		     int incy, double eps_int, double un_int,
		     double *test_ratio);
void test_BLAS_ddot2(int n, enum blas_conj_type conj, double alpha,
		     double beta, double rin, double rout, double r_true_l,
		     double r_true_t, double *x, int incx, double *head_y,
		     double *tail_y, int incy, double eps_int, double un_int,
		     double *test_ratio);
void test_BLAS_cdot2(int n, enum blas_conj_type conj, const void *alpha,
		     const void *beta, const void *rin, const void *rout,
		     double *r_true_l, double *r_true_t, void *x, int incx,
		     void *head_y, void *tail_y, int incy, double eps_int,
		     double un_int, double *test_ratio);
void test_BLAS_zdot2(int n, enum blas_conj_type conj, const void *alpha,
		     const void *beta, const void *rin, const void *rout,
		     double *r_true_l, double *r_true_t, void *x, int incx,
		     void *head_y, void *tail_y, int incy, double eps_int,
		     double un_int, double *test_ratio);
void test_BLAS_cdot2_s_s(int n, enum blas_conj_type conj, const void *alpha,
			 const void *beta, const void *rin, const void *rout,
			 double *r_true_l, double *r_true_t, float *x,
			 int incx, float *head_y, float *tail_y, int incy,
			 double eps_int, double un_int, double *test_ratio);
void test_BLAS_cdot2_s_c(int n, enum blas_conj_type conj, const void *alpha,
			 const void *beta, const void *rin, const void *rout,
			 double *r_true_l, double *r_true_t, float *x,
			 int incx, void *head_y, void *tail_y, int incy,
			 double eps_int, double un_int, double *test_ratio);
void test_BLAS_cdot2_c_s(int n, enum blas_conj_type conj, const void *alpha,
			 const void *beta, const void *rin, const void *rout,
			 double *r_true_l, double *r_true_t, void *x,
			 int incx, float *head_y, float *tail_y, int incy,
			 double eps_int, double un_int, double *test_ratio);
void test_BLAS_zdot2_d_d(int n, enum blas_conj_type conj, const void *alpha,
			 const void *beta, const void *rin, const void *rout,
			 double *r_true_l, double *r_true_t, double *x,
			 int incx, double *head_y, double *tail_y, int incy,
			 double eps_int, double un_int, double *test_ratio);
void test_BLAS_zdot2_d_z(int n, enum blas_conj_type conj, const void *alpha,
			 const void *beta, const void *rin, const void *rout,
			 double *r_true_l, double *r_true_t, double *x,
			 int incx, void *head_y, void *tail_y, int incy,
			 double eps_int, double un_int, double *test_ratio);
void test_BLAS_zdot2_z_d(int n, enum blas_conj_type conj, const void *alpha,
			 const void *beta, const void *rin, const void *rout,
			 double *r_true_l, double *r_true_t, void *x,
			 int incx, double *head_y, double *tail_y, int incy,
			 double eps_int, double un_int, double *test_ratio);
void test_BLAS_ddot2_s_s(int n, enum blas_conj_type conj, double alpha,
			 double beta, double rin, double rout,
			 double r_true_l, double r_true_t, float *x, int incx,
			 float *head_y, float *tail_y, int incy,
			 double eps_int, double un_int, double *test_ratio);
void test_BLAS_ddot2_s_d(int n, enum blas_conj_type conj, double alpha,
			 double beta, double rin, double rout,
			 double r_true_l, double r_true_t, float *x, int incx,
			 double *head_y, double *tail_y, int incy,
			 double eps_int, double un_int, double *test_ratio);
void test_BLAS_ddot2_d_s(int n, enum blas_conj_type conj, double alpha,
			 double beta, double rin, double rout,
			 double r_true_l, double r_true_t, double *x,
			 int incx, float *head_y, float *tail_y, int incy,
			 double eps_int, double un_int, double *test_ratio);
void test_BLAS_zdot2_c_c(int n, enum blas_conj_type conj, const void *alpha,
			 const void *beta, const void *rin, const void *rout,
			 double *r_true_l, double *r_true_t, void *x,
			 int incx, void *head_y, void *tail_y, int incy,
			 double eps_int, double un_int, double *test_ratio);
void test_BLAS_zdot2_c_z(int n, enum blas_conj_type conj, const void *alpha,
			 const void *beta, const void *rin, const void *rout,
			 double *r_true_l, double *r_true_t, void *x,
			 int incx, void *head_y, void *tail_y, int incy,
			 double eps_int, double un_int, double *test_ratio);
void test_BLAS_zdot2_z_c(int n, enum blas_conj_type conj, const void *alpha,
			 const void *beta, const void *rin, const void *rout,
			 double *r_true_l, double *r_true_t, void *x,
			 int incx, void *head_y, void *tail_y, int incy,
			 double eps_int, double un_int, double *test_ratio);


void BLAS_sgemv_testgen(int norm, enum blas_order_type order,
			enum blas_trans_type trans, int m, int n,
			float *alpha, int alpha_flag, float *A, int lda,
			float *x, float *beta, int beta_flag, float *y,
			int *seed, double *r_true_l, double *r_true_t);
void BLAS_dgemv_testgen(int norm, enum blas_order_type order,
			enum blas_trans_type trans, int m, int n,
			double *alpha, int alpha_flag, double *A, int lda,
			double *x, double *beta, int beta_flag, double *y,
			int *seed, double *r_true_l, double *r_true_t);
void BLAS_cgemv_testgen(int norm, enum blas_order_type order,
			enum blas_trans_type trans, int m, int n, void *alpha,
			int alpha_flag, void *A, int lda, void *x, void *beta,
			int beta_flag, void *y, int *seed, double *r_true_l,
			double *r_true_t);
void BLAS_zgemv_testgen(int norm, enum blas_order_type order,
			enum blas_trans_type trans, int m, int n, void *alpha,
			int alpha_flag, void *A, int lda, void *x, void *beta,
			int beta_flag, void *y, int *seed, double *r_true_l,
			double *r_true_t);
void BLAS_cgemv_s_s_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type trans, int m, int n,
			    void *alpha, int alpha_flag, float *A, int lda,
			    float *x, void *beta, int beta_flag, void *y,
			    int *seed, double *r_true_l, double *r_true_t);
void BLAS_cgemv_s_c_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type trans, int m, int n,
			    void *alpha, int alpha_flag, float *A, int lda,
			    void *x, void *beta, int beta_flag, void *y,
			    int *seed, double *r_true_l, double *r_true_t);
void BLAS_cgemv_c_s_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type trans, int m, int n,
			    void *alpha, int alpha_flag, void *A, int lda,
			    float *x, void *beta, int beta_flag, void *y,
			    int *seed, double *r_true_l, double *r_true_t);
void BLAS_zgemv_d_d_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type trans, int m, int n,
			    void *alpha, int alpha_flag, double *A, int lda,
			    double *x, void *beta, int beta_flag, void *y,
			    int *seed, double *r_true_l, double *r_true_t);
void BLAS_zgemv_d_z_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type trans, int m, int n,
			    void *alpha, int alpha_flag, double *A, int lda,
			    void *x, void *beta, int beta_flag, void *y,
			    int *seed, double *r_true_l, double *r_true_t);
void BLAS_zgemv_z_d_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type trans, int m, int n,
			    void *alpha, int alpha_flag, void *A, int lda,
			    double *x, void *beta, int beta_flag, void *y,
			    int *seed, double *r_true_l, double *r_true_t);
void BLAS_dgemv_s_s_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type trans, int m, int n,
			    double *alpha, int alpha_flag, float *A, int lda,
			    float *x, double *beta, int beta_flag, double *y,
			    int *seed, double *r_true_l, double *r_true_t);
void BLAS_dgemv_s_d_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type trans, int m, int n,
			    double *alpha, int alpha_flag, float *A, int lda,
			    double *x, double *beta, int beta_flag, double *y,
			    int *seed, double *r_true_l, double *r_true_t);
void BLAS_dgemv_d_s_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type trans, int m, int n,
			    double *alpha, int alpha_flag, double *A, int lda,
			    float *x, double *beta, int beta_flag, double *y,
			    int *seed, double *r_true_l, double *r_true_t);
void BLAS_zgemv_c_c_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type trans, int m, int n,
			    void *alpha, int alpha_flag, void *A, int lda,
			    void *x, void *beta, int beta_flag, void *y,
			    int *seed, double *r_true_l, double *r_true_t);
void BLAS_zgemv_c_z_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type trans, int m, int n,
			    void *alpha, int alpha_flag, void *A, int lda,
			    void *x, void *beta, int beta_flag, void *y,
			    int *seed, double *r_true_l, double *r_true_t);
void BLAS_zgemv_z_c_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type trans, int m, int n,
			    void *alpha, int alpha_flag, void *A, int lda,
			    void *x, void *beta, int beta_flag, void *y,
			    int *seed, double *r_true_l, double *r_true_t);



void sge_commit_row(enum blas_order_type order, enum blas_trans_type trans,
		    int m, int n, float *a, int lda, float *a_vec, int row);
void dge_commit_row(enum blas_order_type order, enum blas_trans_type trans,
		    int m, int n, double *a, int lda, double *a_vec, int row);
void cge_commit_row(enum blas_order_type order, enum blas_trans_type trans,
		    int m, int n, void *a, int lda, void *a_vec, int row);
void zge_commit_row(enum blas_order_type order, enum blas_trans_type trans,
		    int m, int n, void *a, int lda, void *a_vec, int row);

void sge_commit_col(enum blas_order_type order, enum blas_trans_type trans,
		    int m, int n, float *a, int lda, float *a_vec, int col);
void dge_commit_col(enum blas_order_type order, enum blas_trans_type trans,
		    int m, int n, double *a, int lda, double *a_vec, int col);
void cge_commit_col(enum blas_order_type order, enum blas_trans_type trans,
		    int m, int n, void *a, int lda, void *a_vec, int col);
void zge_commit_col(enum blas_order_type order, enum blas_trans_type trans,
		    int m, int n, void *a, int lda, void *a_vec, int col);

void sge_copy_row(enum blas_order_type order, enum blas_trans_type trans,
		  int m, int n, float *a, int lda, float *a_vec, int row);
void dge_copy_row(enum blas_order_type order, enum blas_trans_type trans,
		  int m, int n, double *a, int lda, double *a_vec, int row);
void cge_copy_row(enum blas_order_type order, enum blas_trans_type trans,
		  int m, int n, void *a, int lda, void *a_vec, int row);
void zge_copy_row(enum blas_order_type order, enum blas_trans_type trans,
		  int m, int n, void *a, int lda, void *a_vec, int row);

void sge_copy_col(enum blas_order_type order, enum blas_trans_type trans,
		  int m, int n, float *a, int lda, float *a_vec, int col);
void dge_copy_col(enum blas_order_type order, enum blas_trans_type trans,
		  int m, int n, double *a, int lda, double *a_vec, int col);
void cge_copy_col(enum blas_order_type order, enum blas_trans_type trans,
		  int m, int n, void *a, int lda, void *a_vec, int col);
void zge_copy_col(enum blas_order_type order, enum blas_trans_type trans,
		  int m, int n, void *a, int lda, void *a_vec, int col);

void sge_copy_matrix(enum blas_order_type order, int m, int n, float *a,
		     int lda, float *b, int ldb);
void dge_copy_matrix(enum blas_order_type order, int m, int n, double *a,
		     int lda, double *b, int ldb);
void cge_copy_matrix(enum blas_order_type order, int m, int n, void *a,
		     int lda, void *b, int ldb);
void zge_copy_matrix(enum blas_order_type order, int m, int n, void *a,
		     int lda, void *b, int ldb);

void sge_print_matrix(float *a, int m, int n, int lda,
		      enum blas_order_type order, const char *name);
void dge_print_matrix(double *a, int m, int n, int lda,
		      enum blas_order_type order, const char *name);
void cge_print_matrix(void *a, int m, int n, int lda,
		      enum blas_order_type order, const char *name);
void zge_print_matrix(void *a, int m, int n, int lda,
		      enum blas_order_type order, const char *name);


void BLAS_sgbmv_testgen(int norm, enum blas_order_type order,
			enum blas_trans_type trans, int m, int n, int kl,
			int ku, float *alpha, int alpha_flag, float *AB,
			int lda, float *x, float *beta, int beta_flag,
			float *y, int *seed, double *r_true_l,
			double *r_true_t);
void BLAS_dgbmv_testgen(int norm, enum blas_order_type order,
			enum blas_trans_type trans, int m, int n, int kl,
			int ku, double *alpha, int alpha_flag, double *AB,
			int lda, double *x, double *beta, int beta_flag,
			double *y, int *seed, double *r_true_l,
			double *r_true_t);
void BLAS_cgbmv_testgen(int norm, enum blas_order_type order,
			enum blas_trans_type trans, int m, int n, int kl,
			int ku, void *alpha, int alpha_flag, void *AB,
			int lda, void *x, void *beta, int beta_flag, void *y,
			int *seed, double *r_true_l, double *r_true_t);
void BLAS_zgbmv_testgen(int norm, enum blas_order_type order,
			enum blas_trans_type trans, int m, int n, int kl,
			int ku, void *alpha, int alpha_flag, void *AB,
			int lda, void *x, void *beta, int beta_flag, void *y,
			int *seed, double *r_true_l, double *r_true_t);
void BLAS_cgbmv_s_s_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type trans, int m, int n, int kl,
			    int ku, void *alpha, int alpha_flag, float *AB,
			    int lda, float *x, void *beta, int beta_flag,
			    void *y, int *seed, double *r_true_l,
			    double *r_true_t);
void BLAS_cgbmv_s_c_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type trans, int m, int n, int kl,
			    int ku, void *alpha, int alpha_flag, float *AB,
			    int lda, void *x, void *beta, int beta_flag,
			    void *y, int *seed, double *r_true_l,
			    double *r_true_t);
void BLAS_cgbmv_c_s_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type trans, int m, int n, int kl,
			    int ku, void *alpha, int alpha_flag, void *AB,
			    int lda, float *x, void *beta, int beta_flag,
			    void *y, int *seed, double *r_true_l,
			    double *r_true_t);
void BLAS_zgbmv_d_d_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type trans, int m, int n, int kl,
			    int ku, void *alpha, int alpha_flag, double *AB,
			    int lda, double *x, void *beta, int beta_flag,
			    void *y, int *seed, double *r_true_l,
			    double *r_true_t);
void BLAS_zgbmv_d_z_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type trans, int m, int n, int kl,
			    int ku, void *alpha, int alpha_flag, double *AB,
			    int lda, void *x, void *beta, int beta_flag,
			    void *y, int *seed, double *r_true_l,
			    double *r_true_t);
void BLAS_zgbmv_z_d_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type trans, int m, int n, int kl,
			    int ku, void *alpha, int alpha_flag, void *AB,
			    int lda, double *x, void *beta, int beta_flag,
			    void *y, int *seed, double *r_true_l,
			    double *r_true_t);
void BLAS_dgbmv_s_s_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type trans, int m, int n, int kl,
			    int ku, double *alpha, int alpha_flag, float *AB,
			    int lda, float *x, double *beta, int beta_flag,
			    double *y, int *seed, double *r_true_l,
			    double *r_true_t);
void BLAS_dgbmv_s_d_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type trans, int m, int n, int kl,
			    int ku, double *alpha, int alpha_flag, float *AB,
			    int lda, double *x, double *beta, int beta_flag,
			    double *y, int *seed, double *r_true_l,
			    double *r_true_t);
void BLAS_dgbmv_d_s_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type trans, int m, int n, int kl,
			    int ku, double *alpha, int alpha_flag, double *AB,
			    int lda, float *x, double *beta, int beta_flag,
			    double *y, int *seed, double *r_true_l,
			    double *r_true_t);
void BLAS_zgbmv_c_c_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type trans, int m, int n, int kl,
			    int ku, void *alpha, int alpha_flag, void *AB,
			    int lda, void *x, void *beta, int beta_flag,
			    void *y, int *seed, double *r_true_l,
			    double *r_true_t);
void BLAS_zgbmv_c_z_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type trans, int m, int n, int kl,
			    int ku, void *alpha, int alpha_flag, void *AB,
			    int lda, void *x, void *beta, int beta_flag,
			    void *y, int *seed, double *r_true_l,
			    double *r_true_t);
void BLAS_zgbmv_z_c_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type trans, int m, int n, int kl,
			    int ku, void *alpha, int alpha_flag, void *AB,
			    int lda, void *x, void *beta, int beta_flag,
			    void *y, int *seed, double *r_true_l,
			    double *r_true_t);

void sgbmv_prepare(enum blas_order_type order, enum blas_trans_type trans,
		   int m, int n, int kl, int ku, float *AB, int lda, float *y,
		   int row, int *nfix2, int *nmix, int *ysize);
void cgbmv_prepare(enum blas_order_type order, enum blas_trans_type trans,
		   int m, int n, int kl, int ku, void *AB, int lda, void *y,
		   int row, int *nfix2, int *nmix, int *ysize);
void dgbmv_prepare(enum blas_order_type order, enum blas_trans_type trans,
		   int m, int n, int kl, int ku, double *AB, int lda,
		   double *y, int row, int *nfix2, int *nmix, int *ysize);
void zgbmv_prepare(enum blas_order_type order, enum blas_trans_type trans,
		   int m, int n, int kl, int ku, void *AB, int lda, void *y,
		   int row, int *nfix2, int *nmix, int *ysize);


void sgbmv_commit(enum blas_order_type order, enum blas_trans_type trans,
		  int m, int n, int kl, int ku, float *AB, int lda, float *y,
		  int row);
void cgbmv_commit(enum blas_order_type order, enum blas_trans_type trans,
		  int m, int n, int kl, int ku, void *AB, int lda, void *y,
		  int row);
void dgbmv_commit(enum blas_order_type order, enum blas_trans_type trans,
		  int m, int n, int kl, int ku, double *AB, int lda,
		  double *y, int row);
void zgbmv_commit(enum blas_order_type order, enum blas_trans_type trans,
		  int m, int n, int kl, int ku, void *AB, int lda, void *y,
		  int row);


void sgbmv_copy(enum blas_order_type order, enum blas_trans_type trans, int m,
		int n, int kl, int ku, const float *AB, int lda, float *y,
		int row);
void cgbmv_copy(enum blas_order_type order, enum blas_trans_type trans, int m,
		int n, int kl, int ku, const void *AB, int lda, void *y,
		int row);
void dgbmv_copy(enum blas_order_type order, enum blas_trans_type trans, int m,
		int n, int kl, int ku, const double *AB, int lda, double *y,
		int row);
void zgbmv_copy(enum blas_order_type order, enum blas_trans_type trans, int m,
		int n, int kl, int ku, const void *AB, int lda, void *y,
		int row);


void BLAS_sgemm_testgen(int norm, enum blas_order_type order,
			enum blas_trans_type transa,
			enum blas_trans_type transb, int m, int n, int k,
			int randomize, float *alpha, int alpha_flag, float *a,
			int lda, float *beta, int beta_flag, float *b,
			int ldb, float *c, int ldc, int *seed,
			double *head_r_true, double *tail_r_true);
void BLAS_dgemm_testgen(int norm, enum blas_order_type order,
			enum blas_trans_type transa,
			enum blas_trans_type transb, int m, int n, int k,
			int randomize, double *alpha, int alpha_flag,
			double *a, int lda, double *beta, int beta_flag,
			double *b, int ldb, double *c, int ldc, int *seed,
			double *head_r_true, double *tail_r_true);
void BLAS_cgemm_testgen(int norm, enum blas_order_type order,
			enum blas_trans_type transa,
			enum blas_trans_type transb, int m, int n, int k,
			int randomize, void *alpha, int alpha_flag, void *a,
			int lda, void *beta, int beta_flag, void *b, int ldb,
			void *c, int ldc, int *seed, double *head_r_true,
			double *tail_r_true);
void BLAS_zgemm_testgen(int norm, enum blas_order_type order,
			enum blas_trans_type transa,
			enum blas_trans_type transb, int m, int n, int k,
			int randomize, void *alpha, int alpha_flag, void *a,
			int lda, void *beta, int beta_flag, void *b, int ldb,
			void *c, int ldc, int *seed, double *head_r_true,
			double *tail_r_true);
void BLAS_cgemm_s_s_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type transa,
			    enum blas_trans_type transb, int m, int n, int k,
			    int randomize, void *alpha, int alpha_flag,
			    float *a, int lda, void *beta, int beta_flag,
			    float *b, int ldb, void *c, int ldc, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_cgemm_s_c_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type transa,
			    enum blas_trans_type transb, int m, int n, int k,
			    int randomize, void *alpha, int alpha_flag,
			    float *a, int lda, void *beta, int beta_flag,
			    void *b, int ldb, void *c, int ldc, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_cgemm_c_s_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type transa,
			    enum blas_trans_type transb, int m, int n, int k,
			    int randomize, void *alpha, int alpha_flag,
			    void *a, int lda, void *beta, int beta_flag,
			    float *b, int ldb, void *c, int ldc, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_zgemm_d_d_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type transa,
			    enum blas_trans_type transb, int m, int n, int k,
			    int randomize, void *alpha, int alpha_flag,
			    double *a, int lda, void *beta, int beta_flag,
			    double *b, int ldb, void *c, int ldc, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_zgemm_d_z_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type transa,
			    enum blas_trans_type transb, int m, int n, int k,
			    int randomize, void *alpha, int alpha_flag,
			    double *a, int lda, void *beta, int beta_flag,
			    void *b, int ldb, void *c, int ldc, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_zgemm_z_d_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type transa,
			    enum blas_trans_type transb, int m, int n, int k,
			    int randomize, void *alpha, int alpha_flag,
			    void *a, int lda, void *beta, int beta_flag,
			    double *b, int ldb, void *c, int ldc, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_dgemm_s_s_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type transa,
			    enum blas_trans_type transb, int m, int n, int k,
			    int randomize, double *alpha, int alpha_flag,
			    float *a, int lda, double *beta, int beta_flag,
			    float *b, int ldb, double *c, int ldc, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_dgemm_s_d_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type transa,
			    enum blas_trans_type transb, int m, int n, int k,
			    int randomize, double *alpha, int alpha_flag,
			    float *a, int lda, double *beta, int beta_flag,
			    double *b, int ldb, double *c, int ldc, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_dgemm_d_s_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type transa,
			    enum blas_trans_type transb, int m, int n, int k,
			    int randomize, double *alpha, int alpha_flag,
			    double *a, int lda, double *beta, int beta_flag,
			    float *b, int ldb, double *c, int ldc, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_zgemm_c_c_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type transa,
			    enum blas_trans_type transb, int m, int n, int k,
			    int randomize, void *alpha, int alpha_flag,
			    void *a, int lda, void *beta, int beta_flag,
			    void *b, int ldb, void *c, int ldc, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_zgemm_c_z_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type transa,
			    enum blas_trans_type transb, int m, int n, int k,
			    int randomize, void *alpha, int alpha_flag,
			    void *a, int lda, void *beta, int beta_flag,
			    void *b, int ldb, void *c, int ldc, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_zgemm_z_c_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type transa,
			    enum blas_trans_type transb, int m, int n, int k,
			    int randomize, void *alpha, int alpha_flag,
			    void *a, int lda, void *beta, int beta_flag,
			    void *b, int ldb, void *c, int ldc, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_sgemv2_testgen(int norm, enum blas_order_type order,
			 enum blas_trans_type trans, int m, int n,
			 float *alpha, int alpha_flag, float *A, int lda,
			 float *head_x, float *tail_x, float *beta,
			 int beta_flag, float *y, int *seed, double *r_true_l,
			 double *r_true_t);
void BLAS_dgemv2_testgen(int norm, enum blas_order_type order,
			 enum blas_trans_type trans, int m, int n,
			 double *alpha, int alpha_flag, double *A, int lda,
			 double *head_x, double *tail_x, double *beta,
			 int beta_flag, double *y, int *seed,
			 double *r_true_l, double *r_true_t);
void BLAS_cgemv2_testgen(int norm, enum blas_order_type order,
			 enum blas_trans_type trans, int m, int n,
			 void *alpha, int alpha_flag, void *A, int lda,
			 void *head_x, void *tail_x, void *beta,
			 int beta_flag, void *y, int *seed, double *r_true_l,
			 double *r_true_t);
void BLAS_zgemv2_testgen(int norm, enum blas_order_type order,
			 enum blas_trans_type trans, int m, int n,
			 void *alpha, int alpha_flag, void *A, int lda,
			 void *head_x, void *tail_x, void *beta,
			 int beta_flag, void *y, int *seed, double *r_true_l,
			 double *r_true_t);
void BLAS_cgemv2_s_s_testgen(int norm, enum blas_order_type order,
			     enum blas_trans_type trans, int m, int n,
			     void *alpha, int alpha_flag, float *A, int lda,
			     float *head_x, float *tail_x, void *beta,
			     int beta_flag, void *y, int *seed,
			     double *r_true_l, double *r_true_t);
void BLAS_cgemv2_s_c_testgen(int norm, enum blas_order_type order,
			     enum blas_trans_type trans, int m, int n,
			     void *alpha, int alpha_flag, float *A, int lda,
			     void *head_x, void *tail_x, void *beta,
			     int beta_flag, void *y, int *seed,
			     double *r_true_l, double *r_true_t);
void BLAS_cgemv2_c_s_testgen(int norm, enum blas_order_type order,
			     enum blas_trans_type trans, int m, int n,
			     void *alpha, int alpha_flag, void *A, int lda,
			     float *head_x, float *tail_x, void *beta,
			     int beta_flag, void *y, int *seed,
			     double *r_true_l, double *r_true_t);
void BLAS_zgemv2_d_d_testgen(int norm, enum blas_order_type order,
			     enum blas_trans_type trans, int m, int n,
			     void *alpha, int alpha_flag, double *A, int lda,
			     double *head_x, double *tail_x, void *beta,
			     int beta_flag, void *y, int *seed,
			     double *r_true_l, double *r_true_t);
void BLAS_zgemv2_d_z_testgen(int norm, enum blas_order_type order,
			     enum blas_trans_type trans, int m, int n,
			     void *alpha, int alpha_flag, double *A, int lda,
			     void *head_x, void *tail_x, void *beta,
			     int beta_flag, void *y, int *seed,
			     double *r_true_l, double *r_true_t);
void BLAS_zgemv2_z_d_testgen(int norm, enum blas_order_type order,
			     enum blas_trans_type trans, int m, int n,
			     void *alpha, int alpha_flag, void *A, int lda,
			     double *head_x, double *tail_x, void *beta,
			     int beta_flag, void *y, int *seed,
			     double *r_true_l, double *r_true_t);
void BLAS_dgemv2_s_s_testgen(int norm, enum blas_order_type order,
			     enum blas_trans_type trans, int m, int n,
			     double *alpha, int alpha_flag, float *A, int lda,
			     float *head_x, float *tail_x, double *beta,
			     int beta_flag, double *y, int *seed,
			     double *r_true_l, double *r_true_t);
void BLAS_dgemv2_s_d_testgen(int norm, enum blas_order_type order,
			     enum blas_trans_type trans, int m, int n,
			     double *alpha, int alpha_flag, float *A, int lda,
			     double *head_x, double *tail_x, double *beta,
			     int beta_flag, double *y, int *seed,
			     double *r_true_l, double *r_true_t);
void BLAS_dgemv2_d_s_testgen(int norm, enum blas_order_type order,
			     enum blas_trans_type trans, int m, int n,
			     double *alpha, int alpha_flag, double *A,
			     int lda, float *head_x, float *tail_x,
			     double *beta, int beta_flag, double *y,
			     int *seed, double *r_true_l, double *r_true_t);
void BLAS_zgemv2_c_c_testgen(int norm, enum blas_order_type order,
			     enum blas_trans_type trans, int m, int n,
			     void *alpha, int alpha_flag, void *A, int lda,
			     void *head_x, void *tail_x, void *beta,
			     int beta_flag, void *y, int *seed,
			     double *r_true_l, double *r_true_t);
void BLAS_zgemv2_c_z_testgen(int norm, enum blas_order_type order,
			     enum blas_trans_type trans, int m, int n,
			     void *alpha, int alpha_flag, void *A, int lda,
			     void *head_x, void *tail_x, void *beta,
			     int beta_flag, void *y, int *seed,
			     double *r_true_l, double *r_true_t);
void BLAS_zgemv2_z_c_testgen(int norm, enum blas_order_type order,
			     enum blas_trans_type trans, int m, int n,
			     void *alpha, int alpha_flag, void *A, int lda,
			     void *head_x, void *tail_x, void *beta,
			     int beta_flag, void *y, int *seed,
			     double *r_true_l, double *r_true_t);

void BLAS_ssymv2_testgen(int norm, enum blas_order_type order,
			 enum blas_uplo_type uplo, int n, float *alpha,
			 int alpha_flag, float *A, int lda, float *head_x,
			 float *tail_x, float *beta, int beta_flag, float *y,
			 int *seed, double *r_true_l, double *r_true_t);
void BLAS_dsymv2_testgen(int norm, enum blas_order_type order,
			 enum blas_uplo_type uplo, int n, double *alpha,
			 int alpha_flag, double *A, int lda, double *head_x,
			 double *tail_x, double *beta, int beta_flag,
			 double *y, int *seed, double *r_true_l,
			 double *r_true_t);
void BLAS_csymv2_testgen(int norm, enum blas_order_type order,
			 enum blas_uplo_type uplo, int n, void *alpha,
			 int alpha_flag, void *A, int lda, void *head_x,
			 void *tail_x, void *beta, int beta_flag, void *y,
			 int *seed, double *r_true_l, double *r_true_t);
void BLAS_zsymv2_testgen(int norm, enum blas_order_type order,
			 enum blas_uplo_type uplo, int n, void *alpha,
			 int alpha_flag, void *A, int lda, void *head_x,
			 void *tail_x, void *beta, int beta_flag, void *y,
			 int *seed, double *r_true_l, double *r_true_t);
void BLAS_csymv2_s_s_testgen(int norm, enum blas_order_type order,
			     enum blas_uplo_type uplo, int n, void *alpha,
			     int alpha_flag, float *A, int lda, float *head_x,
			     float *tail_x, void *beta, int beta_flag,
			     void *y, int *seed, double *r_true_l,
			     double *r_true_t);
void BLAS_csymv2_s_c_testgen(int norm, enum blas_order_type order,
			     enum blas_uplo_type uplo, int n, void *alpha,
			     int alpha_flag, float *A, int lda, void *head_x,
			     void *tail_x, void *beta, int beta_flag, void *y,
			     int *seed, double *r_true_l, double *r_true_t);
void BLAS_csymv2_c_s_testgen(int norm, enum blas_order_type order,
			     enum blas_uplo_type uplo, int n, void *alpha,
			     int alpha_flag, void *A, int lda, float *head_x,
			     float *tail_x, void *beta, int beta_flag,
			     void *y, int *seed, double *r_true_l,
			     double *r_true_t);
void BLAS_zsymv2_d_d_testgen(int norm, enum blas_order_type order,
			     enum blas_uplo_type uplo, int n, void *alpha,
			     int alpha_flag, double *A, int lda,
			     double *head_x, double *tail_x, void *beta,
			     int beta_flag, void *y, int *seed,
			     double *r_true_l, double *r_true_t);
void BLAS_zsymv2_d_z_testgen(int norm, enum blas_order_type order,
			     enum blas_uplo_type uplo, int n, void *alpha,
			     int alpha_flag, double *A, int lda, void *head_x,
			     void *tail_x, void *beta, int beta_flag, void *y,
			     int *seed, double *r_true_l, double *r_true_t);
void BLAS_zsymv2_z_d_testgen(int norm, enum blas_order_type order,
			     enum blas_uplo_type uplo, int n, void *alpha,
			     int alpha_flag, void *A, int lda, double *head_x,
			     double *tail_x, void *beta, int beta_flag,
			     void *y, int *seed, double *r_true_l,
			     double *r_true_t);
void BLAS_dsymv2_s_s_testgen(int norm, enum blas_order_type order,
			     enum blas_uplo_type uplo, int n, double *alpha,
			     int alpha_flag, float *A, int lda, float *head_x,
			     float *tail_x, double *beta, int beta_flag,
			     double *y, int *seed, double *r_true_l,
			     double *r_true_t);
void BLAS_dsymv2_s_d_testgen(int norm, enum blas_order_type order,
			     enum blas_uplo_type uplo, int n, double *alpha,
			     int alpha_flag, float *A, int lda,
			     double *head_x, double *tail_x, double *beta,
			     int beta_flag, double *y, int *seed,
			     double *r_true_l, double *r_true_t);
void BLAS_dsymv2_d_s_testgen(int norm, enum blas_order_type order,
			     enum blas_uplo_type uplo, int n, double *alpha,
			     int alpha_flag, double *A, int lda,
			     float *head_x, float *tail_x, double *beta,
			     int beta_flag, double *y, int *seed,
			     double *r_true_l, double *r_true_t);
void BLAS_zsymv2_c_c_testgen(int norm, enum blas_order_type order,
			     enum blas_uplo_type uplo, int n, void *alpha,
			     int alpha_flag, void *A, int lda, void *head_x,
			     void *tail_x, void *beta, int beta_flag, void *y,
			     int *seed, double *r_true_l, double *r_true_t);
void BLAS_zsymv2_c_z_testgen(int norm, enum blas_order_type order,
			     enum blas_uplo_type uplo, int n, void *alpha,
			     int alpha_flag, void *A, int lda, void *head_x,
			     void *tail_x, void *beta, int beta_flag, void *y,
			     int *seed, double *r_true_l, double *r_true_t);
void BLAS_zsymv2_z_c_testgen(int norm, enum blas_order_type order,
			     enum blas_uplo_type uplo, int n, void *alpha,
			     int alpha_flag, void *A, int lda, void *head_x,
			     void *tail_x, void *beta, int beta_flag, void *y,
			     int *seed, double *r_true_l, double *r_true_t);


void BLAS_sskmv2_testgen(int norm, enum blas_order_type order,
			 enum blas_uplo_type uplo, int n, float *alpha,
			 float *beta, float *a, int lda, float *x_head,
			 float *x_tail, float *y, int *seed,
			 double *head_r_true, double *tail_r_true);
void BLAS_dskmv2_testgen(int norm, enum blas_order_type order,
			 enum blas_uplo_type uplo, int n, double *alpha,
			 double *beta, double *a, int lda, double *x_head,
			 double *x_tail, double *y, int *seed,
			 double *head_r_true, double *tail_r_true);
void BLAS_dskmv2_testgen_d_s(int norm, enum blas_order_type order,
			     enum blas_uplo_type uplo, int n, double *alpha,
			     double *beta, double *a, int lda, float *x_head,
			     float *x_tail, double *y, int *seed,
			     double *head_r_true, double *tail_r_true);
void BLAS_dskmv2_testgen_s_d(int norm, enum blas_order_type order,
			     enum blas_uplo_type uplo, int n, double *alpha,
			     double *beta, float *a, int lda, double *x_head,
			     double *x_tail, double *y, int *seed,
			     double *head_r_true, double *tail_r_true);
void BLAS_dskmv2_testgen_s_s(int norm, enum blas_order_type order,
			     enum blas_uplo_type uplo, int n, double *alpha,
			     double *beta, float *a, int lda, float *x_head,
			     float *x_tail, double *y, int *seed,
			     double *head_r_true, double *tail_r_true);

void BLAS_chemv2_testgen(int norm, enum blas_order_type order,
			 enum blas_uplo_type uplo, int n, void *alpha,
			 int alpha_flag, void *beta, int beta_flag, void *a,
			 int lda, void *x_head, void *x_tail, void *y,
			 int *seed, double *head_r_true, double *tail_r_true);
void BLAS_zhemv2_testgen(int norm, enum blas_order_type order,
			 enum blas_uplo_type uplo, int n, void *alpha,
			 int alpha_flag, void *beta, int beta_flag, void *a,
			 int lda, void *x_head, void *x_tail, void *y,
			 int *seed, double *head_r_true, double *tail_r_true);
void BLAS_zhemv2_c_z_testgen(int norm, enum blas_order_type order,
			     enum blas_uplo_type uplo, int n, void *alpha,
			     int alpha_flag, void *beta, int beta_flag,
			     void *a, int lda, void *x_head, void *x_tail,
			     void *y, int *seed, double *head_r_true,
			     double *tail_r_true);
void BLAS_zhemv2_z_c_testgen(int norm, enum blas_order_type order,
			     enum blas_uplo_type uplo, int n, void *alpha,
			     int alpha_flag, void *beta, int beta_flag,
			     void *a, int lda, void *x_head, void *x_tail,
			     void *y, int *seed, double *head_r_true,
			     double *tail_r_true);
void BLAS_zhemv2_c_c_testgen(int norm, enum blas_order_type order,
			     enum blas_uplo_type uplo, int n, void *alpha,
			     int alpha_flag, void *beta, int beta_flag,
			     void *a, int lda, void *x_head, void *x_tail,
			     void *y, int *seed, double *head_r_true,
			     double *tail_r_true);
void BLAS_zhemv2_z_d_testgen(int norm, enum blas_order_type order,
			     enum blas_uplo_type uplo, int n, void *alpha,
			     int alpha_flag, void *beta, int beta_flag,
			     void *a, int lda, double *x_head, double *x_tail,
			     void *y, int *seed, double *head_r_true,
			     double *tail_r_true);
void BLAS_chemv2_c_s_testgen(int norm, enum blas_order_type order,
			     enum blas_uplo_type uplo, int n, void *alpha,
			     int alpha_flag, void *beta, int beta_flag,
			     void *a, int lda, float *x_head, float *x_tail,
			     void *y, int *seed, double *head_r_true,
			     double *tail_r_true);


void BLAS_sgbmv2_testgen(int norm, enum blas_order_type order,
			 enum blas_trans_type trans, int m, int n, int kl,
			 int ku, float *alpha, int alpha_flag, float *AB,
			 int lda, float *x_head, float *x_tail, float *beta,
			 int beta_flag, float *y, int *seed, double *r_true_l,
			 double *r_true_t);
void BLAS_dgbmv2_testgen(int norm, enum blas_order_type order,
			 enum blas_trans_type trans, int m, int n, int kl,
			 int ku, double *alpha, int alpha_flag, double *AB,
			 int lda, double *x_head, double *x_tail,
			 double *beta, int beta_flag, double *y, int *seed,
			 double *r_true_l, double *r_true_t);
void BLAS_cgbmv2_testgen(int norm, enum blas_order_type order,
			 enum blas_trans_type trans, int m, int n, int kl,
			 int ku, void *alpha, int alpha_flag, void *AB,
			 int lda, void *x_head, void *x_tail, void *beta,
			 int beta_flag, void *y, int *seed, double *r_true_l,
			 double *r_true_t);
void BLAS_zgbmv2_testgen(int norm, enum blas_order_type order,
			 enum blas_trans_type trans, int m, int n, int kl,
			 int ku, void *alpha, int alpha_flag, void *AB,
			 int lda, void *x_head, void *x_tail, void *beta,
			 int beta_flag, void *y, int *seed, double *r_true_l,
			 double *r_true_t);
void BLAS_cgbmv2_s_s_testgen(int norm, enum blas_order_type order,
			     enum blas_trans_type trans, int m, int n, int kl,
			     int ku, void *alpha, int alpha_flag, float *AB,
			     int lda, float *x_head, float *x_tail,
			     void *beta, int beta_flag, void *y, int *seed,
			     double *r_true_l, double *r_true_t);
void BLAS_cgbmv2_s_c_testgen(int norm, enum blas_order_type order,
			     enum blas_trans_type trans, int m, int n, int kl,
			     int ku, void *alpha, int alpha_flag, float *AB,
			     int lda, void *x_head, void *x_tail, void *beta,
			     int beta_flag, void *y, int *seed,
			     double *r_true_l, double *r_true_t);
void BLAS_cgbmv2_c_s_testgen(int norm, enum blas_order_type order,
			     enum blas_trans_type trans, int m, int n, int kl,
			     int ku, void *alpha, int alpha_flag, void *AB,
			     int lda, float *x_head, float *x_tail,
			     void *beta, int beta_flag, void *y, int *seed,
			     double *r_true_l, double *r_true_t);
void BLAS_zgbmv2_d_d_testgen(int norm, enum blas_order_type order,
			     enum blas_trans_type trans, int m, int n, int kl,
			     int ku, void *alpha, int alpha_flag, double *AB,
			     int lda, double *x_head, double *x_tail,
			     void *beta, int beta_flag, void *y, int *seed,
			     double *r_true_l, double *r_true_t);
void BLAS_zgbmv2_d_z_testgen(int norm, enum blas_order_type order,
			     enum blas_trans_type trans, int m, int n, int kl,
			     int ku, void *alpha, int alpha_flag, double *AB,
			     int lda, void *x_head, void *x_tail, void *beta,
			     int beta_flag, void *y, int *seed,
			     double *r_true_l, double *r_true_t);
void BLAS_zgbmv2_z_d_testgen(int norm, enum blas_order_type order,
			     enum blas_trans_type trans, int m, int n, int kl,
			     int ku, void *alpha, int alpha_flag, void *AB,
			     int lda, double *x_head, double *x_tail,
			     void *beta, int beta_flag, void *y, int *seed,
			     double *r_true_l, double *r_true_t);
void BLAS_dgbmv2_s_s_testgen(int norm, enum blas_order_type order,
			     enum blas_trans_type trans, int m, int n, int kl,
			     int ku, double *alpha, int alpha_flag, float *AB,
			     int lda, float *x_head, float *x_tail,
			     double *beta, int beta_flag, double *y,
			     int *seed, double *r_true_l, double *r_true_t);
void BLAS_dgbmv2_s_d_testgen(int norm, enum blas_order_type order,
			     enum blas_trans_type trans, int m, int n, int kl,
			     int ku, double *alpha, int alpha_flag, float *AB,
			     int lda, double *x_head, double *x_tail,
			     double *beta, int beta_flag, double *y,
			     int *seed, double *r_true_l, double *r_true_t);
void BLAS_dgbmv2_d_s_testgen(int norm, enum blas_order_type order,
			     enum blas_trans_type trans, int m, int n, int kl,
			     int ku, double *alpha, int alpha_flag,
			     double *AB, int lda, float *x_head,
			     float *x_tail, double *beta, int beta_flag,
			     double *y, int *seed, double *r_true_l,
			     double *r_true_t);
void BLAS_zgbmv2_c_c_testgen(int norm, enum blas_order_type order,
			     enum blas_trans_type trans, int m, int n, int kl,
			     int ku, void *alpha, int alpha_flag, void *AB,
			     int lda, void *x_head, void *x_tail, void *beta,
			     int beta_flag, void *y, int *seed,
			     double *r_true_l, double *r_true_t);
void BLAS_zgbmv2_c_z_testgen(int norm, enum blas_order_type order,
			     enum blas_trans_type trans, int m, int n, int kl,
			     int ku, void *alpha, int alpha_flag, void *AB,
			     int lda, void *x_head, void *x_tail, void *beta,
			     int beta_flag, void *y, int *seed,
			     double *r_true_l, double *r_true_t);
void BLAS_zgbmv2_z_c_testgen(int norm, enum blas_order_type order,
			     enum blas_trans_type trans, int m, int n, int kl,
			     int ku, void *alpha, int alpha_flag, void *AB,
			     int lda, void *x_head, void *x_tail, void *beta,
			     int beta_flag, void *y, int *seed,
			     double *r_true_l, double *r_true_t);
void BLAS_sge_sum_mv_testgen(int norm, enum blas_order_type order, int m,
			     int n, int randomize, float *alpha,
			     int alpha_flag, float *beta, int beta_flag,
			     float *a, int lda, float *b, int ldb, float *x,
			     int incx, float *alpha_use_ptr, float *a_use,
			     float *b_use, int *seed, double *head_r_true,
			     double *tail_r_true);
void BLAS_dge_sum_mv_testgen(int norm, enum blas_order_type order, int m,
			     int n, int randomize, double *alpha,
			     int alpha_flag, double *beta, int beta_flag,
			     double *a, int lda, double *b, int ldb,
			     double *x, int incx, double *alpha_use_ptr,
			     double *a_use, double *b_use, int *seed,
			     double *head_r_true, double *tail_r_true);
void BLAS_cge_sum_mv_testgen(int norm, enum blas_order_type order, int m,
			     int n, int randomize, void *alpha,
			     int alpha_flag, void *beta, int beta_flag,
			     void *a, int lda, void *b, int ldb, void *x,
			     int incx, void *alpha_use_ptr, void *a_use,
			     void *b_use, int *seed, double *head_r_true,
			     double *tail_r_true);
void BLAS_zge_sum_mv_testgen(int norm, enum blas_order_type order, int m,
			     int n, int randomize, void *alpha,
			     int alpha_flag, void *beta, int beta_flag,
			     void *a, int lda, void *b, int ldb, void *x,
			     int incx, void *alpha_use_ptr, void *a_use,
			     void *b_use, int *seed, double *head_r_true,
			     double *tail_r_true);
void BLAS_cge_sum_mv_s_s_testgen(int norm, enum blas_order_type order, int m,
				 int n, int randomize, void *alpha,
				 int alpha_flag, void *beta, int beta_flag,
				 float *a, int lda, float *b, int ldb,
				 float *x, int incx, void *alpha_use_ptr,
				 float *a_use, float *b_use, int *seed,
				 double *head_r_true, double *tail_r_true);
void BLAS_cge_sum_mv_s_c_testgen(int norm, enum blas_order_type order, int m,
				 int n, int randomize, void *alpha,
				 int alpha_flag, void *beta, int beta_flag,
				 float *a, int lda, float *b, int ldb,
				 void *x, int incx, void *alpha_use_ptr,
				 float *a_use, float *b_use, int *seed,
				 double *head_r_true, double *tail_r_true);
void BLAS_cge_sum_mv_c_s_testgen(int norm, enum blas_order_type order, int m,
				 int n, int randomize, void *alpha,
				 int alpha_flag, void *beta, int beta_flag,
				 void *a, int lda, void *b, int ldb, float *x,
				 int incx, void *alpha_use_ptr, void *a_use,
				 void *b_use, int *seed, double *head_r_true,
				 double *tail_r_true);
void BLAS_zge_sum_mv_d_d_testgen(int norm, enum blas_order_type order, int m,
				 int n, int randomize, void *alpha,
				 int alpha_flag, void *beta, int beta_flag,
				 double *a, int lda, double *b, int ldb,
				 double *x, int incx, void *alpha_use_ptr,
				 double *a_use, double *b_use, int *seed,
				 double *head_r_true, double *tail_r_true);
void BLAS_zge_sum_mv_d_z_testgen(int norm, enum blas_order_type order, int m,
				 int n, int randomize, void *alpha,
				 int alpha_flag, void *beta, int beta_flag,
				 double *a, int lda, double *b, int ldb,
				 void *x, int incx, void *alpha_use_ptr,
				 double *a_use, double *b_use, int *seed,
				 double *head_r_true, double *tail_r_true);
void BLAS_zge_sum_mv_z_d_testgen(int norm, enum blas_order_type order, int m,
				 int n, int randomize, void *alpha,
				 int alpha_flag, void *beta, int beta_flag,
				 void *a, int lda, void *b, int ldb,
				 double *x, int incx, void *alpha_use_ptr,
				 void *a_use, void *b_use, int *seed,
				 double *head_r_true, double *tail_r_true);
void BLAS_dge_sum_mv_s_s_testgen(int norm, enum blas_order_type order, int m,
				 int n, int randomize, double *alpha,
				 int alpha_flag, double *beta, int beta_flag,
				 float *a, int lda, float *b, int ldb,
				 float *x, int incx, double *alpha_use_ptr,
				 float *a_use, float *b_use, int *seed,
				 double *head_r_true, double *tail_r_true);
void BLAS_dge_sum_mv_s_d_testgen(int norm, enum blas_order_type order, int m,
				 int n, int randomize, double *alpha,
				 int alpha_flag, double *beta, int beta_flag,
				 float *a, int lda, float *b, int ldb,
				 double *x, int incx, double *alpha_use_ptr,
				 float *a_use, float *b_use, int *seed,
				 double *head_r_true, double *tail_r_true);
void BLAS_dge_sum_mv_d_s_testgen(int norm, enum blas_order_type order, int m,
				 int n, int randomize, double *alpha,
				 int alpha_flag, double *beta, int beta_flag,
				 double *a, int lda, double *b, int ldb,
				 float *x, int incx, double *alpha_use_ptr,
				 double *a_use, double *b_use, int *seed,
				 double *head_r_true, double *tail_r_true);
void BLAS_zge_sum_mv_c_c_testgen(int norm, enum blas_order_type order, int m,
				 int n, int randomize, void *alpha,
				 int alpha_flag, void *beta, int beta_flag,
				 void *a, int lda, void *b, int ldb, void *x,
				 int incx, void *alpha_use_ptr, void *a_use,
				 void *b_use, int *seed, double *head_r_true,
				 double *tail_r_true);
void BLAS_zge_sum_mv_c_z_testgen(int norm, enum blas_order_type order, int m,
				 int n, int randomize, void *alpha,
				 int alpha_flag, void *beta, int beta_flag,
				 void *a, int lda, void *b, int ldb, void *x,
				 int incx, void *alpha_use_ptr, void *a_use,
				 void *b_use, int *seed, double *head_r_true,
				 double *tail_r_true);
void BLAS_zge_sum_mv_z_c_testgen(int norm, enum blas_order_type order, int m,
				 int n, int randomize, void *alpha,
				 int alpha_flag, void *beta, int beta_flag,
				 void *a, int lda, void *b, int ldb, void *x,
				 int incx, void *alpha_use_ptr, void *a_use,
				 void *b_use, int *seed, double *head_r_true,
				 double *tail_r_true);

void BLAS_ssymv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo,
			int n, int randomize,
			float *alpha, int alpha_flag, float *beta,
			int beta_flag, float *a, int lda, float *x, int incx,
			float *y, int incy, int *seed, double *head_r_true,
			double *tail_r_true);
void BLAS_dsymv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, int n, int randomize,
			double *alpha, int alpha_flag, double *beta,
			int beta_flag, double *a, int lda, double *x,
			int incx, double *y, int incy, int *seed,
			double *head_r_true, double *tail_r_true);
void BLAS_csymv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, int n, int randomize,
			void *alpha, int alpha_flag, void *beta,
			int beta_flag, void *a, int lda, void *x, int incx,
			void *y, int incy, int *seed, double *head_r_true,
			double *tail_r_true);
void BLAS_zsymv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, int n, int randomize,
			void *alpha, int alpha_flag, void *beta,
			int beta_flag, void *a, int lda, void *x, int incx,
			void *y, int incy, int *seed, double *head_r_true,
			double *tail_r_true);
void BLAS_csymv_s_s_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, float *a, int lda, float *x,
			    int incx, void *y, int incy, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_csymv_s_c_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, float *a, int lda, void *x,
			    int incx, void *y, int incy, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_csymv_c_s_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, void *a, int lda, float *x,
			    int incx, void *y, int incy, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_zsymv_d_d_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, double *a, int lda, double *x,
			    int incx, void *y, int incy, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_zsymv_d_z_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, double *a, int lda, void *x,
			    int incx, void *y, int incy, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_zsymv_z_d_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, void *a, int lda, double *x,
			    int incx, void *y, int incy, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_dsymv_s_s_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    double *alpha, int alpha_flag, double *beta,
			    int beta_flag, float *a, int lda, float *x,
			    int incx, double *y, int incy, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_dsymv_s_d_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    double *alpha, int alpha_flag, double *beta,
			    int beta_flag, float *a, int lda, double *x,
			    int incx, double *y, int incy, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_dsymv_d_s_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    double *alpha, int alpha_flag, double *beta,
			    int beta_flag, double *a, int lda, float *x,
			    int incx, double *y, int incy, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_zsymv_c_c_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, void *a, int lda, void *x,
			    int incx, void *y, int incy, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_zsymv_c_z_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, void *a, int lda, void *x,
			    int incx, void *y, int incy, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_zsymv_z_c_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, void *a, int lda, void *x,
			    int incx, void *y, int incy, int *seed,
			    double *head_r_true, double *tail_r_true);

void ssy_commit_row(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, float *a, int lda, float *a_vec, int row);
void dsy_commit_row(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, double *a, int lda, double *a_vec, int row);
void csy_commit_row(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, void *a, int lda, void *a_vec, int row);
void zsy_commit_row(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, void *a, int lda, void *a_vec, int row);

void ssy_copy_row(enum blas_order_type order, enum blas_uplo_type uplo, int n,
		  float *a, int lda, float *a_vec, int row);
void dsy_copy_row(enum blas_order_type order, enum blas_uplo_type uplo, int n,
		  double *a, int lda, double *a_vec, int row);
void csy_copy_row(enum blas_order_type order, enum blas_uplo_type uplo, int n,
		  void *a, int lda, void *a_vec, int row);
void zsy_copy_row(enum blas_order_type order, enum blas_uplo_type uplo, int n,
		  void *a, int lda, void *a_vec, int row);

void ssy_print_matrix(float *a, int n, int lda,
		      enum blas_order_type order, enum blas_uplo_type uplo);
void dsy_print_matrix(double *a, int n, int lda,
		      enum blas_order_type order, enum blas_uplo_type uplo);
void csy_print_matrix(void *a, int n, int lda,
		      enum blas_order_type order, enum blas_uplo_type uplo);
void zsy_print_matrix(void *a, int n, int lda,
		      enum blas_order_type order, enum blas_uplo_type uplo);


void BLAS_ssbmv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo,
			int n, int randomize,
			float *alpha, int alpha_flag, float *beta,
			int beta_flag, float *a, int k, int lda, float *x,
			int incx, float *y, int incy, int *seed,
			double *head_r_true, double *tail_r_true);
void BLAS_dsbmv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, int n, int randomize,
			double *alpha, int alpha_flag, double *beta,
			int beta_flag, double *a, int k, int lda, double *x,
			int incx, double *y, int incy, int *seed,
			double *head_r_true, double *tail_r_true);
void BLAS_csbmv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, int n, int randomize,
			void *alpha, int alpha_flag, void *beta,
			int beta_flag, void *a, int k, int lda, void *x,
			int incx, void *y, int incy, int *seed,
			double *head_r_true, double *tail_r_true);
void BLAS_zsbmv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, int n, int randomize,
			void *alpha, int alpha_flag, void *beta,
			int beta_flag, void *a, int k, int lda, void *x,
			int incx, void *y, int incy, int *seed,
			double *head_r_true, double *tail_r_true);
void BLAS_csbmv_s_s_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, float *a, int k, int lda, float *x,
			    int incx, void *y, int incy, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_csbmv_s_c_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, float *a, int k, int lda, void *x,
			    int incx, void *y, int incy, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_csbmv_c_s_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, void *a, int k, int lda, float *x,
			    int incx, void *y, int incy, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_zsbmv_d_d_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, double *a, int k, int lda,
			    double *x, int incx, void *y, int incy, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_zsbmv_d_z_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, double *a, int k, int lda, void *x,
			    int incx, void *y, int incy, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_zsbmv_z_d_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, void *a, int k, int lda, double *x,
			    int incx, void *y, int incy, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_dsbmv_s_s_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    double *alpha, int alpha_flag, double *beta,
			    int beta_flag, float *a, int k, int lda, float *x,
			    int incx, double *y, int incy, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_dsbmv_s_d_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    double *alpha, int alpha_flag, double *beta,
			    int beta_flag, float *a, int k, int lda,
			    double *x, int incx, double *y, int incy,
			    int *seed, double *head_r_true,
			    double *tail_r_true);
void BLAS_dsbmv_d_s_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    double *alpha, int alpha_flag, double *beta,
			    int beta_flag, double *a, int k, int lda,
			    float *x, int incx, double *y, int incy,
			    int *seed, double *head_r_true,
			    double *tail_r_true);
void BLAS_zsbmv_c_c_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, void *a, int k, int lda, void *x,
			    int incx, void *y, int incy, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_zsbmv_c_z_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, void *a, int k, int lda, void *x,
			    int incx, void *y, int incy, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_zsbmv_z_c_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, void *a, int k, int lda, void *x,
			    int incx, void *y, int incy, int *seed,
			    double *head_r_true, double *tail_r_true);
void ssbmv_commit_row(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, float *a, int k, int lda, float *a_vec, int row);
void dsbmv_commit_row(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, double *a, int k, int lda, double *a_vec,
		      int row);
void csbmv_commit_row(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, void *a, int k, int lda, void *a_vec, int row);
void zsbmv_commit_row(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, void *a, int k, int lda, void *a_vec, int row);

void ssbmv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, float *a, int k, int lda, float *a_vec, int row);
void dsbmv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, double *a, int k, int lda, double *a_vec, int row);
void csbmv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, void *a, int k, int lda, void *a_vec, int row);
void zsbmv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, void *a, int k, int lda, void *a_vec, int row);

void sprint_sbmv_matrix(float *a, int n, int k, int lda,
			enum blas_order_type order, enum blas_uplo_type uplo);
void dprint_sbmv_matrix(double *a, int n, int k, int lda,
			enum blas_order_type order, enum blas_uplo_type uplo);
void cprint_sbmv_matrix(void *a, int n, int k, int lda,
			enum blas_order_type order, enum blas_uplo_type uplo);
void zprint_sbmv_matrix(void *a, int n, int k, int lda,
			enum blas_order_type order, enum blas_uplo_type uplo);


void BLAS_sskew_testgen_hemv(int norm, enum blas_order_type order,
			     enum blas_uplo_type uplo,
			     int n, int randomize,
			     float *alpha, float *beta,
			     float *a, int lda, float *x, int incx, float *y,
			     int incy, int *seed, double *head_r_true,
			     double *tail_r_true);
void BLAS_dskew_testgen_hemv(int norm, enum blas_order_type order,
			     enum blas_uplo_type uplo, int n, int randomize,
			     double *alpha, double *beta, double *a, int lda,
			     double *x, int incx, double *y, int incy,
			     int *seed, double *head_r_true,
			     double *tail_r_true);
void BLAS_dskew_testgen_hemv_d_s(int norm, enum blas_order_type order,
				 enum blas_uplo_type uplo, int n,
				 int randomize, double *alpha, double *beta,
				 double *a, int lda, float *x, int incx,
				 double *y, int incy, int *seed,
				 double *head_r_true, double *tail_r_true);
void BLAS_dskew_testgen_hemv_s_d(int norm, enum blas_order_type order,
				 enum blas_uplo_type uplo, int n,
				 int randomize, double *alpha, double *beta,
				 float *a, int lda, double *x, int incx,
				 double *y, int incy, int *seed,
				 double *head_r_true, double *tail_r_true);
void BLAS_dskew_testgen_hemv_s_s(int norm, enum blas_order_type order,
				 enum blas_uplo_type uplo, int n,
				 int randomize, double *alpha, double *beta,
				 float *a, int lda, float *x, int incx,
				 double *y, int incy, int *seed,
				 double *head_r_true, double *tail_r_true);

void BLAS_chemv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo,
			int n, int randomize,
			void *alpha, int alpha_flag, void *beta,
			int beta_flag, void *a, int lda, void *x, int incx,
			void *y, int incy, int *seed, double *head_r_true,
			double *tail_r_true);
void BLAS_zhemv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, int n, int randomize,
			void *alpha, int alpha_flag, void *beta,
			int beta_flag, void *a, int lda, void *x, int incx,
			void *y, int incy, int *seed, double *head_r_true,
			double *tail_r_true);
void BLAS_zhemv_c_z_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, void *a, int lda, void *x,
			    int incx, void *y, int incy, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_zhemv_z_c_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, void *a, int lda, void *x,
			    int incx, void *y, int incy, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_zhemv_c_c_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, void *a, int lda, void *x,
			    int incx, void *y, int incy, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_zhemv_z_d_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, void *a, int lda, double *x,
			    int incx, void *y, int incy, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_chemv_c_s_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, void *a, int lda, float *x,
			    int incx, void *y, int incy, int *seed,
			    double *head_r_true, double *tail_r_true);
void che_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_side_type side, int n, void *a, int lda,
		  void *a_vec, int row);
void zhe_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_side_type side, int n, void *a, int lda,
		  void *a_vec, int row);

void che_print_matrix(void *a, int n, int lda,
		      enum blas_order_type order, enum blas_uplo_type uplo);
void zhe_print_matrix(void *a, int n, int lda,
		      enum blas_order_type order, enum blas_uplo_type uplo);

void sskew_commit_row(enum blas_order_type order, enum blas_uplo_type uplo,
		      enum blas_side_type side, int n, float *a, int lda,
		      float *a_vec, int row);
void dskew_commit_row(enum blas_order_type order, enum blas_uplo_type uplo,
		      enum blas_side_type side, int n, double *a, int lda,
		      double *a_vec, int row);

void sskew_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		    enum blas_side_type side, int n, float *a, int lda,
		    float *a_vec, int row);
void dskew_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		    enum blas_side_type side, int n, double *a, int lda,
		    double *a_vec, int row);

void BLAS_sskew_testgen_hbmv(int norm, enum blas_order_type order,
			     enum blas_uplo_type uplo,
			     int n, int randomize,
			     float *alpha, float *beta,
			     float *a, int k, int lda, float *x, int incx,
			     float *y, int incy, int *seed,
			     double *head_r_true, double *tail_r_true);
void BLAS_dskew_testgen_hbmv(int norm, enum blas_order_type order,
			     enum blas_uplo_type uplo, int n, int randomize,
			     double *alpha, double *beta, double *a, int k,
			     int lda, double *x, int incx, double *y,
			     int incy, int *seed, double *head_r_true,
			     double *tail_r_true);
void BLAS_dskew_testgen_hbmv_d_s(int norm, enum blas_order_type order,
				 enum blas_uplo_type uplo, int n,
				 int randomize, double *alpha, double *beta,
				 double *a, int k, int lda, float *x,
				 int incx, double *y, int incy, int *seed,
				 double *head_r_true, double *tail_r_true);
void BLAS_dskew_testgen_hbmv_s_d(int norm, enum blas_order_type order,
				 enum blas_uplo_type uplo, int n,
				 int randomize, double *alpha, double *beta,
				 float *a, int k, int lda, double *x,
				 int incx, double *y, int incy, int *seed,
				 double *head_r_true, double *tail_r_true);
void BLAS_dskew_testgen_hbmv_s_s(int norm, enum blas_order_type order,
				 enum blas_uplo_type uplo, int n,
				 int randomize, double *alpha, double *beta,
				 float *a, int k, int lda, float *x, int incx,
				 double *y, int incy, int *seed,
				 double *head_r_true, double *tail_r_true);

void BLAS_chbmv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, int n, int randomize,
			void *alpha, int alpha_flag, void *beta,
			int beta_flag, void *a, int k, int lda, void *x,
			int incx, void *y, int incy, int *seed,
			double *head_r_true, double *tail_r_true);
void BLAS_zhbmv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, int n, int randomize,
			void *alpha, int alpha_flag, void *beta,
			int beta_flag, void *a, int k, int lda, void *x,
			int incx, void *y, int incy, int *seed,
			double *head_r_true, double *tail_r_true);
void BLAS_zhbmv_c_z_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, void *a, int k, int lda, void *x,
			    int incx, void *y, int incy, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_zhbmv_z_c_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, void *a, int k, int lda, void *x,
			    int incx, void *y, int incy, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_zhbmv_c_c_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, void *a, int k, int lda, void *x,
			    int incx, void *y, int incy, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_zhbmv_z_d_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, void *a, int k, int lda, double *x,
			    int incx, void *y, int incy, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_chbmv_c_s_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, void *a, int k, int lda, float *x,
			    int incx, void *y, int incy, int *seed,
			    double *head_r_true, double *tail_r_true);

void sskew_commit_row_hbmv(enum blas_order_type order,
			   enum blas_uplo_type uplo, int n, float *a, int k,
			   int lda, float *a_vec, int row);
void dskew_commit_row_hbmv(enum blas_order_type order,
			   enum blas_uplo_type uplo, int n, double *a, int k,
			   int lda, double *a_vec, int row);

void chbmv_commit_row(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, void *a, int k, int lda, void *a_vec, int row);
void zhbmv_commit_row(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, void *a, int k, int lda, void *a_vec, int row);

void sskew_copy_row_hbmv(enum blas_order_type order, enum blas_uplo_type uplo,
			 int n, float *a, int k, int lda,
			 float *a_vec, int row);
void dskew_copy_row_hbmv(enum blas_order_type order, enum blas_uplo_type uplo,
			 int n, double *a, int k, int lda,
			 double *a_vec, int row);

void chbmv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, void *a, int k, int lda, void *a_vec, int row);
void zhbmv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, void *a, int k, int lda, void *a_vec, int row);

void cprint_hbmv_matrix(void *a, int n, int k, int lda,
			enum blas_order_type order, enum blas_uplo_type uplo);
void zprint_hbmv_matrix(void *a, int n, int k, int lda,
			enum blas_order_type order, enum blas_uplo_type uplo);

void BLAS_ssymm_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, enum blas_side_type side,
			int m, int n, int randomize, float *alpha,
			int alpha_flag, float *beta, int beta_flag, float *a,
			int lda, float *b, int ldb, float *c, int ldc,
			int *seed, double *head_r_true, double *tail_r_true);
void BLAS_dsymm_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, enum blas_side_type side,
			int m, int n, int randomize, double *alpha,
			int alpha_flag, double *beta, int beta_flag,
			double *a, int lda, double *b, int ldb, double *c,
			int ldc, int *seed, double *head_r_true,
			double *tail_r_true);
void BLAS_csymm_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, enum blas_side_type side,
			int m, int n, int randomize, void *alpha,
			int alpha_flag, void *beta, int beta_flag, void *a,
			int lda, void *b, int ldb, void *c, int ldc,
			int *seed, double *head_r_true, double *tail_r_true);
void BLAS_zsymm_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, enum blas_side_type side,
			int m, int n, int randomize, void *alpha,
			int alpha_flag, void *beta, int beta_flag, void *a,
			int lda, void *b, int ldb, void *c, int ldc,
			int *seed, double *head_r_true, double *tail_r_true);
void BLAS_csymm_s_s_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    enum blas_side_type side, int m, int n,
			    int randomize, void *alpha, int alpha_flag,
			    void *beta, int beta_flag, float *a, int lda,
			    float *b, int ldb, void *c, int ldc, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_csymm_s_c_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    enum blas_side_type side, int m, int n,
			    int randomize, void *alpha, int alpha_flag,
			    void *beta, int beta_flag, float *a, int lda,
			    void *b, int ldb, void *c, int ldc, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_csymm_c_s_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    enum blas_side_type side, int m, int n,
			    int randomize, void *alpha, int alpha_flag,
			    void *beta, int beta_flag, void *a, int lda,
			    float *b, int ldb, void *c, int ldc, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_zsymm_d_d_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    enum blas_side_type side, int m, int n,
			    int randomize, void *alpha, int alpha_flag,
			    void *beta, int beta_flag, double *a, int lda,
			    double *b, int ldb, void *c, int ldc, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_zsymm_d_z_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    enum blas_side_type side, int m, int n,
			    int randomize, void *alpha, int alpha_flag,
			    void *beta, int beta_flag, double *a, int lda,
			    void *b, int ldb, void *c, int ldc, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_zsymm_z_d_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    enum blas_side_type side, int m, int n,
			    int randomize, void *alpha, int alpha_flag,
			    void *beta, int beta_flag, void *a, int lda,
			    double *b, int ldb, void *c, int ldc, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_dsymm_s_s_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    enum blas_side_type side, int m, int n,
			    int randomize, double *alpha, int alpha_flag,
			    double *beta, int beta_flag, float *a, int lda,
			    float *b, int ldb, double *c, int ldc, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_dsymm_s_d_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    enum blas_side_type side, int m, int n,
			    int randomize, double *alpha, int alpha_flag,
			    double *beta, int beta_flag, float *a, int lda,
			    double *b, int ldb, double *c, int ldc, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_dsymm_d_s_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    enum blas_side_type side, int m, int n,
			    int randomize, double *alpha, int alpha_flag,
			    double *beta, int beta_flag, double *a, int lda,
			    float *b, int ldb, double *c, int ldc, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_zsymm_c_c_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    enum blas_side_type side, int m, int n,
			    int randomize, void *alpha, int alpha_flag,
			    void *beta, int beta_flag, void *a, int lda,
			    void *b, int ldb, void *c, int ldc, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_zsymm_c_z_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    enum blas_side_type side, int m, int n,
			    int randomize, void *alpha, int alpha_flag,
			    void *beta, int beta_flag, void *a, int lda,
			    void *b, int ldb, void *c, int ldc, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_zsymm_z_c_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    enum blas_side_type side, int m, int n,
			    int randomize, void *alpha, int alpha_flag,
			    void *beta, int beta_flag, void *a, int lda,
			    void *b, int ldb, void *c, int ldc, int *seed,
			    double *head_r_true, double *tail_r_true);

void BLAS_chemm_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, enum blas_side_type side,
			int m, int n, int randomize, void *alpha,
			int alpha_flag, void *beta, int beta_flag, void *a,
			int lda, void *b, int ldb, void *c, int ldc,
			int *seed, double *head_r_true, double *tail_r_true);
void BLAS_zhemm_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, enum blas_side_type side,
			int m, int n, int randomize, void *alpha,
			int alpha_flag, void *beta, int beta_flag, void *a,
			int lda, void *b, int ldb, void *c, int ldc,
			int *seed, double *head_r_true, double *tail_r_true);
void BLAS_zhemm_c_z_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    enum blas_side_type side, int m, int n,
			    int randomize, void *alpha, int alpha_flag,
			    void *beta, int beta_flag, void *a, int lda,
			    void *b, int ldb, void *c, int ldc, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_zhemm_z_c_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    enum blas_side_type side, int m, int n,
			    int randomize, void *alpha, int alpha_flag,
			    void *beta, int beta_flag, void *a, int lda,
			    void *b, int ldb, void *c, int ldc, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_zhemm_c_c_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    enum blas_side_type side, int m, int n,
			    int randomize, void *alpha, int alpha_flag,
			    void *beta, int beta_flag, void *a, int lda,
			    void *b, int ldb, void *c, int ldc, int *seed,
			    double *head_r_true, double *tail_r_true);

void BLAS_zhemm_z_d_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    enum blas_side_type side, int m, int n,
			    int randomize, void *alpha, int alpha_flag,
			    void *beta, int beta_flag, void *a, int lda,
			    double *b, int ldb, void *c, int ldc, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_chemm_c_s_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    enum blas_side_type side, int m, int n,
			    int randomize, void *alpha, int alpha_flag,
			    void *beta, int beta_flag, void *a, int lda,
			    float *b, int ldb, void *c, int ldc, int *seed,
			    double *head_r_true, double *tail_r_true);

void BLAS_sskew_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, enum blas_side_type side,
			int m, int n, float *alpha, int alpha_flag,
			float *beta, int beta_flag, float *a, int lda,
			float *b, int ldb, float *c, int ldc, int *seed,
			double *head_r_true, double *tail_r_true);
void BLAS_dskew_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, enum blas_side_type side,
			int m, int n, double *alpha, int alpha_flag,
			double *beta, int beta_flag, double *a, int lda,
			double *b, int ldb, double *c, int ldc, int *seed,
			double *head_r_true, double *tail_r_true);
void BLAS_dskew_d_s_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    enum blas_side_type side, int m, int n,
			    double *alpha, int alpha_flag, double *beta,
			    int beta_flag, double *a, int lda, float *b,
			    int ldb, double *c, int ldc, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_dskew_s_d_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    enum blas_side_type side, int m, int n,
			    double *alpha, int alpha_flag, double *beta,
			    int beta_flag, float *a, int lda, double *b,
			    int ldb, double *c, int ldc, int *seed,
			    double *head_r_true, double *tail_r_true);
void BLAS_dskew_s_s_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo,
			    enum blas_side_type side, int m, int n,
			    double *alpha, int alpha_flag, double *beta,
			    int beta_flag, float *a, int lda, float *b,
			    int ldb, double *c, int ldc, int *seed,
			    double *head_r_true, double *tail_r_true);


void BLAS_chpmv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo,
			int n, int randomize,
			void *alpha, int alpha_flag, void *beta,
			int beta_flag, void *a, void *x, int incx, void *y,
			int incy, int *seed, double *r_true_l,
			double *r_true_t);
void BLAS_zhpmv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, int n, int randomize,
			void *alpha, int alpha_flag, void *beta,
			int beta_flag, void *a, void *x, int incx, void *y,
			int incy, int *seed, double *r_true_l,
			double *r_true_t);
void BLAS_zhpmv_c_z_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, void *a, void *x, int incx,
			    void *y, int incy, int *seed, double *r_true_l,
			    double *r_true_t);
void BLAS_zhpmv_z_c_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, void *a, void *x, int incx,
			    void *y, int incy, int *seed, double *r_true_l,
			    double *r_true_t);
void BLAS_zhpmv_c_c_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, void *a, void *x, int incx,
			    void *y, int incy, int *seed, double *r_true_l,
			    double *r_true_t);
void BLAS_zhpmv_z_d_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, void *a, double *x, int incx,
			    void *y, int incy, int *seed, double *r_true_l,
			    double *r_true_t);
void BLAS_chpmv_c_s_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, void *a, float *x, int incx,
			    void *y, int incy, int *seed, double *r_true_l,
			    double *r_true_t);
void chpmv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, void *a, void *a_vec, int row);
void zhpmv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, void *a, void *a_vec, int row);

void chpmv_commit_row(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, void *a, void *a_vec, int row);
void zhpmv_commit_row(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, void *a, void *a_vec, int row);

void chpmv_pack_matrix(enum blas_order_type order, enum blas_uplo_type uplo,
		       int n, void *a_packed, void *a_full, int lda);
void zhpmv_pack_matrix(enum blas_order_type order, enum blas_uplo_type uplo,
		       int n, void *a_packed, void *a_full, int lda);

void cprint_hpmv_matrix(void *a, int n,
			enum blas_order_type order, enum blas_uplo_type uplo);
void zprint_hpmv_matrix(void *a, int n,
			enum blas_order_type order, enum blas_uplo_type uplo);

void BLAS_sspmv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, int n, int randomize,
			float *alpha, int alpha_flag, float *beta,
			int beta_flag, float *a, float *x, int incx, float *y,
			int incy, int *seed, double *r_true_l,
			double *r_true_t);
void BLAS_dspmv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, int n, int randomize,
			double *alpha, int alpha_flag, double *beta,
			int beta_flag, double *a, double *x, int incx,
			double *y, int incy, int *seed, double *r_true_l,
			double *r_true_t);
void BLAS_cspmv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, int n, int randomize,
			void *alpha, int alpha_flag, void *beta,
			int beta_flag, void *a, void *x, int incx, void *y,
			int incy, int *seed, double *r_true_l,
			double *r_true_t);
void BLAS_zspmv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, int n, int randomize,
			void *alpha, int alpha_flag, void *beta,
			int beta_flag, void *a, void *x, int incx, void *y,
			int incy, int *seed, double *r_true_l,
			double *r_true_t);
void BLAS_cspmv_s_s_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, float *a, float *x, int incx,
			    void *y, int incy, int *seed, double *r_true_l,
			    double *r_true_t);
void BLAS_cspmv_s_c_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, float *a, void *x, int incx,
			    void *y, int incy, int *seed, double *r_true_l,
			    double *r_true_t);
void BLAS_cspmv_c_s_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, void *a, float *x, int incx,
			    void *y, int incy, int *seed, double *r_true_l,
			    double *r_true_t);
void BLAS_zspmv_d_d_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, double *a, double *x, int incx,
			    void *y, int incy, int *seed, double *r_true_l,
			    double *r_true_t);
void BLAS_zspmv_d_z_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, double *a, void *x, int incx,
			    void *y, int incy, int *seed, double *r_true_l,
			    double *r_true_t);
void BLAS_zspmv_z_d_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, void *a, double *x, int incx,
			    void *y, int incy, int *seed, double *r_true_l,
			    double *r_true_t);
void BLAS_dspmv_s_s_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    double *alpha, int alpha_flag, double *beta,
			    int beta_flag, float *a, float *x, int incx,
			    double *y, int incy, int *seed, double *r_true_l,
			    double *r_true_t);
void BLAS_dspmv_s_d_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    double *alpha, int alpha_flag, double *beta,
			    int beta_flag, float *a, double *x, int incx,
			    double *y, int incy, int *seed, double *r_true_l,
			    double *r_true_t);
void BLAS_dspmv_d_s_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    double *alpha, int alpha_flag, double *beta,
			    int beta_flag, double *a, float *x, int incx,
			    double *y, int incy, int *seed, double *r_true_l,
			    double *r_true_t);
void BLAS_zspmv_c_c_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, void *a, void *x, int incx,
			    void *y, int incy, int *seed, double *r_true_l,
			    double *r_true_t);
void BLAS_zspmv_c_z_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, void *a, void *x, int incx,
			    void *y, int incy, int *seed, double *r_true_l,
			    double *r_true_t);
void BLAS_zspmv_z_c_testgen(int norm, enum blas_order_type order,
			    enum blas_uplo_type uplo, int n, int randomize,
			    void *alpha, int alpha_flag, void *beta,
			    int beta_flag, void *a, void *x, int incx,
			    void *y, int incy, int *seed, double *r_true_l,
			    double *r_true_t);

void sspmv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, const float *a, float *a_vec, int row);

void dspmv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, const double *a, double *a_vec, int row);

void cspmv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, const void *a, void *a_vec, int row);

void zspmv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, const void *a, void *a_vec, int row);


void sspmv_commit_row(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, float *a, const float *a_vec, int row);

void dspmv_commit_row(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, double *a, const double *a_vec, int row);

void cspmv_commit_row(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, void *a, const void *a_vec, int row);

void zspmv_commit_row(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, void *a, const void *a_vec, int row);

void sspmv_pack_matrix(enum blas_order_type order, enum blas_uplo_type uplo,
		       int n, float *a_packed, float *a_full, int lda);
void dspmv_pack_matrix(enum blas_order_type order, enum blas_uplo_type uplo,
		       int n, double *a_packed, double *a_full, int lda);
void cspmv_pack_matrix(enum blas_order_type order, enum blas_uplo_type uplo,
		       int n, void *a_packed, void *a_full, int lda);
void zspmv_pack_matrix(enum blas_order_type order, enum blas_uplo_type uplo,
		       int n, void *a_packed, void *a_full, int lda);

void sprint_spmv_matrix(float *a, int n,
			enum blas_order_type order, enum blas_uplo_type uplo);
void dprint_spmv_matrix(double *a, int n,
			enum blas_order_type order, enum blas_uplo_type uplo);
void cprint_spmv_matrix(void *a, int n,
			enum blas_order_type order, enum blas_uplo_type uplo);
void zprint_spmv_matrix(void *a, int n,
			enum blas_order_type order, enum blas_uplo_type uplo);

void BLAS_ssum_testgen(int n, int norm, float *x, int *seed,
		       double *sum_true_l, double *sum_true_t);
void BLAS_dsum_testgen(int n, int norm, double *x, int *seed,
		       double *sum_true_l, double *sum_true_t);
void BLAS_csum_testgen(int n, int norm, void *x, int *seed,
		       double *sum_true_l, double *sum_true_t);
void BLAS_zsum_testgen(int n, int norm, void *x, int *seed,
		       double *sum_true_l, double *sum_true_t);
void test_BLAS_ssum(int n, float sum_comp, double sum_true_l,
		    double sum_true_t, float *x, int incx,
		    double eps_int, double un_int, double *test_ratio);
void test_BLAS_dsum(int n, double sum_comp, double sum_true_l,
		    double sum_true_t, double *x, int incx,
		    double eps_int, double un_int, double *test_ratio);
void test_BLAS_csum(int n, const void *sum_comp, double *sum_true_l,
		    double *sum_true_t, void *x, int incx,
		    double eps_int, double un_int, double *test_ratio);
void test_BLAS_zsum(int n, const void *sum_comp, double *sum_true_l,
		    double *sum_true_t, void *x, int incx,
		    double eps_int, double un_int, double *test_ratio);

void BLAS_sdot_x_testgen(int n, int n_fix2, int n_mix, int norm,
			 enum blas_conj_type conj, float *alpha,
			 int alpha_flag, float *beta, int beta_flag,
			 double *x_l, double *x_t, float *y, int *seed,
			 float *r, double *r_true_l, double *r_true_t);

void BLAS_ddot_x_testgen(int n, int n_fix2, int n_mix, int norm,
			 enum blas_conj_type conj, double *alpha,
			 int alpha_flag, double *beta, int beta_flag,
			 double *x_l, double *x_t, double *y, int *seed,
			 double *r, double *r_true_l, double *r_true_t);

void BLAS_ddot_s_x_testgen(int n, int n_fix2, int n_mix, int norm,
			   enum blas_conj_type conj, double *alpha,
			   int alpha_flag, double *beta, int beta_flag,
			   double *x_l, double *x_t, float *y, int *seed,
			   double *r, double *r_true_l, double *r_true_t);

void testgen_BLAS_sdot_x(int n, int n_fix2, int n_mix, int norm,
			 enum blas_conj_type conj, float *alpha,
			 int alpha_flag, float *beta, int beta_flag,
			 double *x_l, double *x_t, float *y, int *seed,
			 float *r, double *r_true_l, double *r_true_t);

void BLAS_strsv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, enum blas_trans_type trans,
			enum blas_diag_type diag, int n, float *alpha,
			int alpha_flag, float *T, int lda, float *x,
			int *seed, double *head_r_true, double *tail_r_true,
			int row, enum blas_prec_type prec);
void BLAS_dtrsv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, enum blas_trans_type trans,
			enum blas_diag_type diag, int n, double *alpha,
			int alpha_flag, double *T, int lda, double *x,
			int *seed, double *head_r_true, double *tail_r_true,
			int row, enum blas_prec_type prec);
void BLAS_dtrsv_s_testgen(int norm, enum blas_order_type order,
			  enum blas_uplo_type uplo,
			  enum blas_trans_type trans,
			  enum blas_diag_type diag, int n, double *alpha,
			  int alpha_flag, float *T, int lda, double *x,
			  int *seed, double *head_r_true, double *tail_r_true,
			  int row, enum blas_prec_type prec);
void BLAS_ctrsv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, enum blas_trans_type trans,
			enum blas_diag_type diag, int n, void *alpha,
			int alpha_flag, void *T, int lda, void *x, int *seed,
			double *head_r_true, double *tail_r_true, int row,
			enum blas_prec_type prec);
void BLAS_ztrsv_c_testgen(int norm, enum blas_order_type order,
			  enum blas_uplo_type uplo,
			  enum blas_trans_type trans,
			  enum blas_diag_type diag, int n, void *alpha,
			  int alpha_flag, void *T, int lda, void *x,
			  int *seed, double *head_r_true, double *tail_r_true,
			  int row, enum blas_prec_type prec);
void BLAS_ztrsv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, enum blas_trans_type trans,
			enum blas_diag_type diag, int n, void *alpha,
			int alpha_flag, void *T, int lda, void *x, int *seed,
			double *head_r_true, double *tail_r_true, int row,
			enum blas_prec_type prec);
void BLAS_ctrsv_s_testgen(int norm, enum blas_order_type order,
			  enum blas_uplo_type uplo,
			  enum blas_trans_type trans,
			  enum blas_diag_type diag, int n, void *alpha,
			  int alpha_flag, float *T, int lda, void *x,
			  int *seed, double *head_r_true, double *tail_r_true,
			  int row, enum blas_prec_type prec);
void BLAS_ztrsv_d_testgen(int norm, enum blas_order_type order,
			  enum blas_uplo_type uplo,
			  enum blas_trans_type trans,
			  enum blas_diag_type diag, int n, void *alpha,
			  int alpha_flag, double *T, int lda, void *x,
			  int *seed, double *head_r_true, double *tail_r_true,
			  int row, enum blas_prec_type prec);
void strsv_commit(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, int length, float *T, int lda,
		  const float *y, int row);
void dtrsv_commit(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, int length, double *T, int lda,
		  const double *y, int row);

void BLAS_stbsv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, enum blas_trans_type trans,
			enum blas_diag_type diag, int n, int k, int randomize,
			float *alpha, int alpha_flag, float *T, int ldt,
			float *x, int *seed, double *head_r_true,
			double *tail_r_true, int row,
			enum blas_prec_type prec);
void BLAS_dtbsv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, enum blas_trans_type trans,
			enum blas_diag_type diag, int n, int k, int randomize,
			double *alpha, int alpha_flag, double *T, int ldt,
			double *x, int *seed, double *head_r_true,
			double *tail_r_true, int row,
			enum blas_prec_type prec);
void BLAS_dtbsv_s_testgen(int norm, enum blas_order_type order,
			  enum blas_uplo_type uplo,
			  enum blas_trans_type trans,
			  enum blas_diag_type diag, int n, int k,
			  int randomize, double *alpha, int alpha_flag,
			  float *T, int ldt, double *x, int *seed,
			  double *head_r_true, double *tail_r_true, int row,
			  enum blas_prec_type prec);
void BLAS_ctbsv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, enum blas_trans_type trans,
			enum blas_diag_type diag, int n, int k, int randomize,
			void *alpha, int alpha_flag, void *T, int ldt,
			void *x, int *seed, double *head_r_true,
			double *tail_r_true, int row,
			enum blas_prec_type prec);
void BLAS_ztbsv_c_testgen(int norm, enum blas_order_type order,
			  enum blas_uplo_type uplo,
			  enum blas_trans_type trans,
			  enum blas_diag_type diag, int n, int k,
			  int randomize, void *alpha, int alpha_flag, void *T,
			  int ldt, void *x, int *seed, double *head_r_true,
			  double *tail_r_true, int row,
			  enum blas_prec_type prec);
void BLAS_ztbsv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, enum blas_trans_type trans,
			enum blas_diag_type diag, int n, int k, int randomize,
			void *alpha, int alpha_flag, void *T, int ldt,
			void *x, int *seed, double *head_r_true,
			double *tail_r_true, int row,
			enum blas_prec_type prec);
void BLAS_ctbsv_s_testgen(int norm, enum blas_order_type order,
			  enum blas_uplo_type uplo,
			  enum blas_trans_type trans,
			  enum blas_diag_type diag, int n, int k,
			  int randomize, void *alpha, int alpha_flag,
			  float *T, int ldt, void *x, int *seed,
			  double *head_r_true, double *tail_r_true, int row,
			  enum blas_prec_type prec);
void BLAS_ztbsv_d_testgen(int norm, enum blas_order_type order,
			  enum blas_uplo_type uplo,
			  enum blas_trans_type trans,
			  enum blas_diag_type diag, int n, int k,
			  int randomize, void *alpha, int alpha_flag,
			  double *T, int ldt, void *x, int *seed,
			  double *head_r_true, double *tail_r_true, int row,
			  enum blas_prec_type prec);
void stbsv_copy(enum blas_order_type order, enum blas_uplo_type uplo,
		enum blas_trans_type trans, int n, int k, const float *T,
		int ldt, float *y, int row);
void dtbsv_copy(enum blas_order_type order, enum blas_uplo_type uplo,
		enum blas_trans_type trans, int n, int k, const double *T,
		int ldt, double *y, int row);
void ctbsv_copy(enum blas_order_type order, enum blas_uplo_type uplo,
		enum blas_trans_type trans, int n, int k, const void *T,
		int ldt, void *y, int row);
void ztbsv_copy(enum blas_order_type order, enum blas_uplo_type uplo,
		enum blas_trans_type trans, int n, int k, const void *T,
		int ldt, void *y, int row);

void stbsv_commit(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, int n, int k, float *T, int ldt,
		  float *y, int row);
void dtbsv_commit(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, int n, int k, double *T,
		  int ldt, double *y, int row);
void ctbsv_commit(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, int n, int k, void *T, int ldt,
		  void *y, int row);
void ztbsv_commit(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, int n, int k, void *T, int ldt,
		  void *y, int row);

void sprint_tbsv_matrix(float *T, int n, int k, int ldt,
			enum blas_order_type order, enum blas_uplo_type uplo,
			enum blas_trans_type trans);
void dprint_tbsv_matrix(double *T, int n, int k, int ldt,
			enum blas_order_type order, enum blas_uplo_type uplo,
			enum blas_trans_type trans);
void cprint_tbsv_matrix(void *T, int n, int k, int ldt,
			enum blas_order_type order, enum blas_uplo_type uplo,
			enum blas_trans_type trans);
void zprint_tbsv_matrix(void *T, int n, int k, int ldt,
			enum blas_order_type order, enum blas_uplo_type uplo,
			enum blas_trans_type trans);


void BLAS_stpmv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, enum blas_trans_type trans,
			enum blas_diag_type diag, int n, float *alpha,
			int alpha_flag, float *tp, float *x, int *seed,
			double *head_r_true, double *tail_r_true);
void BLAS_dtpmv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, enum blas_trans_type trans,
			enum blas_diag_type diag, int n, double *alpha,
			int alpha_flag, double *tp, double *x, int *seed,
			double *head_r_true, double *tail_r_true);
void BLAS_ctpmv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, enum blas_trans_type trans,
			enum blas_diag_type diag, int n, void *alpha,
			int alpha_flag, void *tp, void *x, int *seed,
			double *head_r_true, double *tail_r_true);
void BLAS_ztpmv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, enum blas_trans_type trans,
			enum blas_diag_type diag, int n, void *alpha,
			int alpha_flag, void *tp, void *x, int *seed,
			double *head_r_true, double *tail_r_true);
void BLAS_dtpmv_s_testgen(int norm, enum blas_order_type order,
			  enum blas_uplo_type uplo,
			  enum blas_trans_type trans,
			  enum blas_diag_type diag, int n, double *alpha,
			  int alpha_flag, float *tp, double *x, int *seed,
			  double *head_r_true, double *tail_r_true);
void BLAS_ztpmv_c_testgen(int norm, enum blas_order_type order,
			  enum blas_uplo_type uplo,
			  enum blas_trans_type trans,
			  enum blas_diag_type diag, int n, void *alpha,
			  int alpha_flag, void *tp, void *x, int *seed,
			  double *head_r_true, double *tail_r_true);
void BLAS_ctpmv_s_testgen(int norm, enum blas_order_type order,
			  enum blas_uplo_type uplo,
			  enum blas_trans_type trans,
			  enum blas_diag_type diag, int n, void *alpha,
			  int alpha_flag, float *tp, void *x, int *seed,
			  double *head_r_true, double *tail_r_true);
void BLAS_ztpmv_d_testgen(int norm, enum blas_order_type order,
			  enum blas_uplo_type uplo,
			  enum blas_trans_type trans,
			  enum blas_diag_type diag, int n, void *alpha,
			  int alpha_flag, double *tp, void *x, int *seed,
			  double *head_r_true, double *tail_r_true);
void stpmv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		    enum blas_trans_type trans, int n, const float *a,
		    float *a_vec, int row);
void dtpmv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		    enum blas_trans_type trans, int n, const double *a,
		    double *a_vec, int row);
void ctpmv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		    enum blas_trans_type trans, int n, const void *a,
		    void *a_vec, int row);
void ztpmv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		    enum blas_trans_type trans, int n, const void *a,
		    void *a_vec, int row);

void stpmv_commit_row(enum blas_order_type order, enum blas_uplo_type uplo,
		      enum blas_trans_type trans, int n, float *a,
		      const float *a_vec, int row);
void dtpmv_commit_row(enum blas_order_type order, enum blas_uplo_type uplo,
		      enum blas_trans_type trans, int n, double *a,
		      const double *a_vec, int row);
void ctpmv_commit_row(enum blas_order_type order, enum blas_uplo_type uplo,
		      enum blas_trans_type trans, int n, void *a,
		      const void *a_vec, int row);
void ztpmv_commit_row(enum blas_order_type order, enum blas_uplo_type uplo,
		      enum blas_trans_type trans, int n, void *a,
		      const void *a_vec, int row);

void BLAS_strmv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, enum blas_trans_type trans,
			enum blas_diag_type diag, int n, float *alpha,
			int alpha_flag, float *T, int ldt, float *x,
			int *seed, double *head_r_true, double *tail_r_true);
void BLAS_dtrmv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, enum blas_trans_type trans,
			enum blas_diag_type diag, int n, double *alpha,
			int alpha_flag, double *T, int ldt, double *x,
			int *seed, double *head_r_true, double *tail_r_true);
void BLAS_dtrmv_s_testgen(int norm, enum blas_order_type order,
			  enum blas_uplo_type uplo,
			  enum blas_trans_type trans,
			  enum blas_diag_type diag, int n, double *alpha,
			  int alpha_flag, float *T, int ldt, double *x,
			  int *seed, double *head_r_true,
			  double *tail_r_true);
void BLAS_ctrmv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, enum blas_trans_type trans,
			enum blas_diag_type diag, int n, void *alpha,
			int alpha_flag, void *T, int ldt, void *x, int *seed,
			double *head_r_true, double *tail_r_true);
void BLAS_ztrmv_c_testgen(int norm, enum blas_order_type order,
			  enum blas_uplo_type uplo,
			  enum blas_trans_type trans,
			  enum blas_diag_type diag, int n, void *alpha,
			  int alpha_flag, void *T, int ldt, void *x,
			  int *seed, double *head_r_true,
			  double *tail_r_true);
void BLAS_ztrmv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, enum blas_trans_type trans,
			enum blas_diag_type diag, int n, void *alpha,
			int alpha_flag, void *T, int ldt, void *x, int *seed,
			double *head_r_true, double *tail_r_true);
void BLAS_ctrmv_s_testgen(int norm, enum blas_order_type order,
			  enum blas_uplo_type uplo,
			  enum blas_trans_type trans,
			  enum blas_diag_type diag, int n, void *alpha,
			  int alpha_flag, float *T, int ldt, void *x,
			  int *seed, double *head_r_true,
			  double *tail_r_true);
void BLAS_ztrmv_d_testgen(int norm, enum blas_order_type order,
			  enum blas_uplo_type uplo,
			  enum blas_trans_type trans,
			  enum blas_diag_type diag, int n, void *alpha,
			  int alpha_flag, double *T, int ldt, void *x,
			  int *seed, double *head_r_true,
			  double *tail_r_true);
void str_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, int n, const float *T, int ldt,
		  float *y, int row);
void dtr_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, int n, const double *T, int ldt,
		  double *y, int row);
void ctr_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, int n, const void *T, int ldt,
		  void *y, int row);
void ztr_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, int n, const void *T, int ldt,
		  void *y, int row);

#endif /* BLAS_EXTENDED_TEST_H */
