#ifndef BLAS_EXTENDED_PROTO_H
#define BLAS_EXTENDED_PROTO_H


void BLAS_ddot_d_s(enum blas_conj_type conj, int n, double alpha,
		   const double *x, int incx, double beta,
		   const float *y, int incy, double *r);
void BLAS_ddot_s_d(enum blas_conj_type conj, int n, double alpha,
		   const float *x, int incx, double beta,
		   const double *y, int incy, double *r);
void BLAS_ddot_s_s(enum blas_conj_type conj, int n, double alpha,
		   const float *x, int incx, double beta,
		   const float *y, int incy, double *r);
void BLAS_zdot_z_c(enum blas_conj_type conj, int n, const void *alpha,
		   const void *x, int incx, const void *beta,
		   const void *y, int incy, void *r);
void BLAS_zdot_c_z(enum blas_conj_type conj, int n, const void *alpha,
		   const void *x, int incx, const void *beta,
		   const void *y, int incy, void *r);
void BLAS_zdot_c_c(enum blas_conj_type conj, int n, const void *alpha,
		   const void *x, int incx, const void *beta,
		   const void *y, int incy, void *r);
void BLAS_cdot_c_s(enum blas_conj_type conj, int n, const void *alpha,
		   const void *x, int incx, const void *beta,
		   const float *y, int incy, void *r);
void BLAS_cdot_s_c(enum blas_conj_type conj, int n, const void *alpha,
		   const float *x, int incx, const void *beta,
		   const void *y, int incy, void *r);
void BLAS_cdot_s_s(enum blas_conj_type conj, int n, const void *alpha,
		   const float *x, int incx, const void *beta,
		   const float *y, int incy, void *r);
void BLAS_zdot_z_d(enum blas_conj_type conj, int n, const void *alpha,
		   const void *x, int incx, const void *beta,
		   const double *y, int incy, void *r);
void BLAS_zdot_d_z(enum blas_conj_type conj, int n, const void *alpha,
		   const double *x, int incx, const void *beta,
		   const void *y, int incy, void *r);
void BLAS_zdot_d_d(enum blas_conj_type conj, int n, const void *alpha,
		   const double *x, int incx, const void *beta,
		   const double *y, int incy, void *r);
void BLAS_sdot_x(enum blas_conj_type conj, int n, float alpha,
		 const float *x, int incx, float beta,
		 const float *y, int incy,
		 float *r, enum blas_prec_type prec);
void BLAS_ddot_x(enum blas_conj_type conj, int n, double alpha,
		 const double *x, int incx, double beta,
		 const double *y, int incy,
		 double *r, enum blas_prec_type prec);
void BLAS_cdot_x(enum blas_conj_type conj, int n, const void *alpha,
		 const void *x, int incx, const void *beta,
		 const void *y, int incy, void *r, enum blas_prec_type prec);
void BLAS_zdot_x(enum blas_conj_type conj, int n, const void *alpha,
		 const void *x, int incx, const void *beta,
		 const void *y, int incy, void *r, enum blas_prec_type prec);
void BLAS_ddot_d_s_x(enum blas_conj_type conj, int n, double alpha,
		     const double *x, int incx, double beta,
		     const float *y, int incy,
		     double *r, enum blas_prec_type prec);
void BLAS_ddot_s_d_x(enum blas_conj_type conj, int n, double alpha,
		     const float *x, int incx, double beta,
		     const double *y, int incy,
		     double *r, enum blas_prec_type prec);
void BLAS_ddot_s_s_x(enum blas_conj_type conj, int n, double alpha,
		     const float *x, int incx, double beta,
		     const float *y, int incy,
		     double *r, enum blas_prec_type prec);
void BLAS_zdot_z_c_x(enum blas_conj_type conj, int n, const void *alpha,
		     const void *x, int incx, const void *beta,
		     const void *y, int incy,
		     void *r, enum blas_prec_type prec);
void BLAS_zdot_c_z_x(enum blas_conj_type conj, int n, const void *alpha,
		     const void *x, int incx, const void *beta,
		     const void *y, int incy,
		     void *r, enum blas_prec_type prec);
void BLAS_zdot_c_c_x(enum blas_conj_type conj, int n, const void *alpha,
		     const void *x, int incx, const void *beta,
		     const void *y, int incy,
		     void *r, enum blas_prec_type prec);
void BLAS_cdot_c_s_x(enum blas_conj_type conj, int n, const void *alpha,
		     const void *x, int incx, const void *beta,
		     const float *y, int incy,
		     void *r, enum blas_prec_type prec);
void BLAS_cdot_s_c_x(enum blas_conj_type conj, int n, const void *alpha,
		     const float *x, int incx, const void *beta,
		     const void *y, int incy,
		     void *r, enum blas_prec_type prec);
void BLAS_cdot_s_s_x(enum blas_conj_type conj, int n, const void *alpha,
		     const float *x, int incx, const void *beta,
		     const float *y, int incy,
		     void *r, enum blas_prec_type prec);
void BLAS_zdot_z_d_x(enum blas_conj_type conj, int n, const void *alpha,
		     const void *x, int incx, const void *beta,
		     const double *y, int incy,
		     void *r, enum blas_prec_type prec);
void BLAS_zdot_d_z_x(enum blas_conj_type conj, int n, const void *alpha,
		     const double *x, int incx, const void *beta,
		     const void *y, int incy,
		     void *r, enum blas_prec_type prec);
void BLAS_zdot_d_d_x(enum blas_conj_type conj, int n, const void *alpha,
		     const double *x, int incx, const void *beta,
		     const double *y, int incy,
		     void *r, enum blas_prec_type prec);


void BLAS_ssum_x(int n, const float *x, int incx,
		 float *sum, enum blas_prec_type prec);
void BLAS_dsum_x(int n, const double *x, int incx,
		 double *sum, enum blas_prec_type prec);
void BLAS_csum_x(int n, const void *x, int incx,
		 void *sum, enum blas_prec_type prec);
void BLAS_zsum_x(int n, const void *x, int incx,
		 void *sum, enum blas_prec_type prec);


void BLAS_daxpby_s(int n, double alpha, const float *x, int incx,
		   double beta, double *y, int incy);
void BLAS_caxpby_s(int n, const void *alpha, const float *x, int incx,
		   const void *beta, void *y, int incy);
void BLAS_zaxpby_c(int n, const void *alpha, const void *x, int incx,
		   const void *beta, void *y, int incy);
void BLAS_zaxpby_d(int n, const void *alpha, const double *x, int incx,
		   const void *beta, void *y, int incy);
void BLAS_saxpby_x(int n, float alpha, const float *x, int incx,
		   float beta, float *y, int incy, enum blas_prec_type prec);
void BLAS_daxpby_x(int n, double alpha, const double *x, int incx,
		   double beta, double *y,
		   int incy, enum blas_prec_type prec);
void BLAS_caxpby_x(int n, const void *alpha, const void *x, int incx,
		   const void *beta, void *y,
		   int incy, enum blas_prec_type prec);
void BLAS_zaxpby_x(int n, const void *alpha, const void *x, int incx,
		   const void *beta, void *y,
		   int incy, enum blas_prec_type prec);
void BLAS_daxpby_s_x(int n, double alpha, const float *x, int incx,
		     double beta, double *y,
		     int incy, enum blas_prec_type prec);
void BLAS_zaxpby_c_x(int n, const void *alpha, const void *x, int incx,
		     const void *beta, void *y,
		     int incy, enum blas_prec_type prec);
void BLAS_caxpby_s_x(int n, const void *alpha, const float *x, int incx,
		     const void *beta, void *y,
		     int incy, enum blas_prec_type prec);
void BLAS_zaxpby_d_x(int n, const void *alpha, const double *x, int incx,
		     const void *beta, void *y,
		     int incy, enum blas_prec_type prec);


void BLAS_dwaxpby_d_s(int n, double alpha, const double *x, int incx,
		      double beta, const float *y, int incy, double *w,
		      int incw);
void BLAS_dwaxpby_s_d(int n, double alpha, const float *x, int incx,
		      double beta, const double *y, int incy, double *w,
		      int incw);
void BLAS_dwaxpby_s_s(int n, double alpha, const float *x, int incx,
		      double beta, const float *y, int incy, double *w,
		      int incw);
void BLAS_zwaxpby_z_c(int n, const void *alpha, const void *x, int incx,
		      const void *beta, const void *y, int incy, void *w,
		      int incw);
void BLAS_zwaxpby_c_z(int n, const void *alpha, const void *x, int incx,
		      const void *beta, const void *y, int incy, void *w,
		      int incw);
void BLAS_zwaxpby_c_c(int n, const void *alpha, const void *x, int incx,
		      const void *beta, const void *y, int incy, void *w,
		      int incw);
void BLAS_cwaxpby_c_s(int n, const void *alpha, const void *x, int incx,
		      const void *beta, const float *y, int incy, void *w,
		      int incw);
void BLAS_cwaxpby_s_c(int n, const void *alpha, const float *x, int incx,
		      const void *beta, const void *y, int incy, void *w,
		      int incw);
void BLAS_cwaxpby_s_s(int n, const void *alpha, const float *x, int incx,
		      const void *beta, const float *y, int incy, void *w,
		      int incw);
void BLAS_zwaxpby_z_d(int n, const void *alpha, const void *x, int incx,
		      const void *beta, const double *y, int incy, void *w,
		      int incw);
void BLAS_zwaxpby_d_z(int n, const void *alpha, const double *x, int incx,
		      const void *beta, const void *y, int incy, void *w,
		      int incw);
void BLAS_zwaxpby_d_d(int n, const void *alpha, const double *x, int incx,
		      const void *beta, const double *y, int incy, void *w,
		      int incw);
void BLAS_swaxpby_x(int n, float alpha, const float *x, int incx,
		    float beta, const float *y, int incy, float *w,
		    int incw, enum blas_prec_type prec);
void BLAS_dwaxpby_x(int n, double alpha, const double *x, int incx,
		    double beta, const double *y, int incy, double *w,
		    int incw, enum blas_prec_type prec);
void BLAS_cwaxpby_x(int n, const void *alpha, const void *x, int incx,
		    const void *beta, const void *y, int incy, void *w,
		    int incw, enum blas_prec_type prec);
void BLAS_zwaxpby_x(int n, const void *alpha, const void *x, int incx,
		    const void *beta, const void *y, int incy, void *w,
		    int incw, enum blas_prec_type prec);
void BLAS_dwaxpby_d_s_x(int n, double alpha, const double *x, int incx,
			double beta, const float *y, int incy, double *w,
			int incw, enum blas_prec_type prec);
void BLAS_dwaxpby_s_d_x(int n, double alpha, const float *x, int incx,
			double beta, const double *y, int incy, double *w,
			int incw, enum blas_prec_type prec);
void BLAS_dwaxpby_s_s_x(int n, double alpha, const float *x, int incx,
			double beta, const float *y, int incy, double *w,
			int incw, enum blas_prec_type prec);
void BLAS_zwaxpby_z_c_x(int n, const void *alpha, const void *x, int incx,
			const void *beta, const void *y, int incy, void *w,
			int incw, enum blas_prec_type prec);
void BLAS_zwaxpby_c_z_x(int n, const void *alpha, const void *x, int incx,
			const void *beta, const void *y, int incy, void *w,
			int incw, enum blas_prec_type prec);
void BLAS_zwaxpby_c_c_x(int n, const void *alpha, const void *x, int incx,
			const void *beta, const void *y, int incy, void *w,
			int incw, enum blas_prec_type prec);
void BLAS_cwaxpby_c_s_x(int n, const void *alpha, const void *x, int incx,
			const void *beta, const float *y, int incy, void *w,
			int incw, enum blas_prec_type prec);
void BLAS_cwaxpby_s_c_x(int n, const void *alpha, const float *x, int incx,
			const void *beta, const void *y, int incy, void *w,
			int incw, enum blas_prec_type prec);
void BLAS_cwaxpby_s_s_x(int n, const void *alpha, const float *x, int incx,
			const void *beta, const float *y, int incy, void *w,
			int incw, enum blas_prec_type prec);
void BLAS_zwaxpby_z_d_x(int n, const void *alpha, const void *x, int incx,
			const void *beta, const double *y, int incy, void *w,
			int incw, enum blas_prec_type prec);
void BLAS_zwaxpby_d_z_x(int n, const void *alpha, const double *x, int incx,
			const void *beta, const void *y, int incy, void *w,
			int incw, enum blas_prec_type prec);
void BLAS_zwaxpby_d_d_x(int n, const void *alpha, const double *x, int incx,
			const void *beta, const double *y, int incy, void *w,
			int incw, enum blas_prec_type prec);


void BLAS_dgemv_d_s(enum blas_order_type order, enum blas_trans_type trans,
		    int m, int n, double alpha, const double *a, int lda,
		    const float *x, int incx, double beta, double *y,
		    int incy);
void BLAS_dgemv_s_d(enum blas_order_type order, enum blas_trans_type trans,
		    int m, int n, double alpha, const float *a, int lda,
		    const double *x, int incx, double beta, double *y,
		    int incy);
void BLAS_dgemv_s_s(enum blas_order_type order, enum blas_trans_type trans,
		    int m, int n, double alpha, const float *a, int lda,
		    const float *x, int incx, double beta, double *y,
		    int incy);
void BLAS_zgemv_z_c(enum blas_order_type order, enum blas_trans_type trans,
		    int m, int n, const void *alpha, const void *a, int lda,
		    const void *x, int incx, const void *beta, void *y,
		    int incy);
void BLAS_zgemv_c_z(enum blas_order_type order, enum blas_trans_type trans,
		    int m, int n, const void *alpha, const void *a, int lda,
		    const void *x, int incx, const void *beta, void *y,
		    int incy);
void BLAS_zgemv_c_c(enum blas_order_type order, enum blas_trans_type trans,
		    int m, int n, const void *alpha, const void *a, int lda,
		    const void *x, int incx, const void *beta, void *y,
		    int incy);
void BLAS_cgemv_c_s(enum blas_order_type order, enum blas_trans_type trans,
		    int m, int n, const void *alpha, const void *a, int lda,
		    const float *x, int incx, const void *beta, void *y,
		    int incy);
void BLAS_cgemv_s_c(enum blas_order_type order, enum blas_trans_type trans,
		    int m, int n, const void *alpha, const float *a, int lda,
		    const void *x, int incx, const void *beta, void *y,
		    int incy);
void BLAS_cgemv_s_s(enum blas_order_type order, enum blas_trans_type trans,
		    int m, int n, const void *alpha, const float *a, int lda,
		    const float *x, int incx, const void *beta, void *y,
		    int incy);
void BLAS_zgemv_z_d(enum blas_order_type order, enum blas_trans_type trans,
		    int m, int n, const void *alpha, const void *a, int lda,
		    const double *x, int incx, const void *beta, void *y,
		    int incy);
void BLAS_zgemv_d_z(enum blas_order_type order, enum blas_trans_type trans,
		    int m, int n, const void *alpha, const double *a, int lda,
		    const void *x, int incx, const void *beta, void *y,
		    int incy);
void BLAS_zgemv_d_d(enum blas_order_type order, enum blas_trans_type trans,
		    int m, int n, const void *alpha, const double *a, int lda,
		    const double *x, int incx, const void *beta, void *y,
		    int incy);
void BLAS_sgemv_x(enum blas_order_type order, enum blas_trans_type trans,
		  int m, int n, float alpha, const float *a, int lda,
		  const float *x, int incx, float beta, float *y,
		  int incy, enum blas_prec_type prec);
void BLAS_dgemv_x(enum blas_order_type order, enum blas_trans_type trans,
		  int m, int n, double alpha, const double *a, int lda,
		  const double *x, int incx, double beta, double *y,
		  int incy, enum blas_prec_type prec);
void BLAS_cgemv_x(enum blas_order_type order, enum blas_trans_type trans,
		  int m, int n, const void *alpha, const void *a, int lda,
		  const void *x, int incx, const void *beta, void *y,
		  int incy, enum blas_prec_type prec);
void BLAS_zgemv_x(enum blas_order_type order, enum blas_trans_type trans,
		  int m, int n, const void *alpha, const void *a, int lda,
		  const void *x, int incx, const void *beta, void *y,
		  int incy, enum blas_prec_type prec);
void BLAS_dgemv_d_s_x(enum blas_order_type order, enum blas_trans_type trans,
		      int m, int n, double alpha, const double *a, int lda,
		      const float *x, int incx, double beta, double *y,
		      int incy, enum blas_prec_type prec);
void BLAS_dgemv_s_d_x(enum blas_order_type order, enum blas_trans_type trans,
		      int m, int n, double alpha, const float *a, int lda,
		      const double *x, int incx, double beta, double *y,
		      int incy, enum blas_prec_type prec);
void BLAS_dgemv_s_s_x(enum blas_order_type order, enum blas_trans_type trans,
		      int m, int n, double alpha, const float *a, int lda,
		      const float *x, int incx, double beta, double *y,
		      int incy, enum blas_prec_type prec);
void BLAS_zgemv_z_c_x(enum blas_order_type order, enum blas_trans_type trans,
		      int m, int n, const void *alpha, const void *a, int lda,
		      const void *x, int incx, const void *beta, void *y,
		      int incy, enum blas_prec_type prec);
void BLAS_zgemv_c_z_x(enum blas_order_type order, enum blas_trans_type trans,
		      int m, int n, const void *alpha, const void *a, int lda,
		      const void *x, int incx, const void *beta, void *y,
		      int incy, enum blas_prec_type prec);
void BLAS_zgemv_c_c_x(enum blas_order_type order, enum blas_trans_type trans,
		      int m, int n, const void *alpha, const void *a, int lda,
		      const void *x, int incx, const void *beta, void *y,
		      int incy, enum blas_prec_type prec);
void BLAS_cgemv_c_s_x(enum blas_order_type order, enum blas_trans_type trans,
		      int m, int n, const void *alpha, const void *a, int lda,
		      const float *x, int incx, const void *beta, void *y,
		      int incy, enum blas_prec_type prec);
void BLAS_cgemv_s_c_x(enum blas_order_type order, enum blas_trans_type trans,
		      int m, int n, const void *alpha, const float *a,
		      int lda, const void *x, int incx, const void *beta,
		      void *y, int incy, enum blas_prec_type prec);
void BLAS_cgemv_s_s_x(enum blas_order_type order, enum blas_trans_type trans,
		      int m, int n, const void *alpha, const float *a,
		      int lda, const float *x, int incx, const void *beta,
		      void *y, int incy, enum blas_prec_type prec);
void BLAS_zgemv_z_d_x(enum blas_order_type order, enum blas_trans_type trans,
		      int m, int n, const void *alpha, const void *a, int lda,
		      const double *x, int incx, const void *beta, void *y,
		      int incy, enum blas_prec_type prec);
void BLAS_zgemv_d_z_x(enum blas_order_type order, enum blas_trans_type trans,
		      int m, int n, const void *alpha, const double *a,
		      int lda, const void *x, int incx, const void *beta,
		      void *y, int incy, enum blas_prec_type prec);
void BLAS_zgemv_d_d_x(enum blas_order_type order, enum blas_trans_type trans,
		      int m, int n, const void *alpha, const double *a,
		      int lda, const double *x, int incx, const void *beta,
		      void *y, int incy, enum blas_prec_type prec);


void BLAS_dge_sum_mv_d_s(enum blas_order_type order, int m, int n,
			 double alpha, const double *a, int lda,
			 const float *x, int incx,
			 double beta, const double *b, int ldb,
			 double *y, int incy);
void BLAS_dge_sum_mv_s_d(enum blas_order_type order, int m, int n,
			 double alpha, const float *a, int lda,
			 const double *x, int incx,
			 double beta, const float *b, int ldb,
			 double *y, int incy);
void BLAS_dge_sum_mv_s_s(enum blas_order_type order, int m, int n,
			 double alpha, const float *a, int lda,
			 const float *x, int incx,
			 double beta, const float *b, int ldb,
			 double *y, int incy);
void BLAS_zge_sum_mv_z_c(enum blas_order_type order, int m, int n,
			 const void *alpha, const void *a, int lda,
			 const void *x, int incx,
			 const void *beta, const void *b, int ldb,
			 void *y, int incy);
void BLAS_zge_sum_mv_c_z(enum blas_order_type order, int m, int n,
			 const void *alpha, const void *a, int lda,
			 const void *x, int incx,
			 const void *beta, const void *b, int ldb,
			 void *y, int incy);
void BLAS_zge_sum_mv_c_c(enum blas_order_type order, int m, int n,
			 const void *alpha, const void *a, int lda,
			 const void *x, int incx,
			 const void *beta, const void *b, int ldb,
			 void *y, int incy);
void BLAS_cge_sum_mv_c_s(enum blas_order_type order, int m, int n,
			 const void *alpha, const void *a, int lda,
			 const float *x, int incx,
			 const void *beta, const void *b, int ldb,
			 void *y, int incy);
void BLAS_cge_sum_mv_s_c(enum blas_order_type order, int m, int n,
			 const void *alpha, const float *a, int lda,
			 const void *x, int incx,
			 const void *beta, const float *b, int ldb,
			 void *y, int incy);
void BLAS_cge_sum_mv_s_s(enum blas_order_type order, int m, int n,
			 const void *alpha, const float *a, int lda,
			 const float *x, int incx,
			 const void *beta, const float *b, int ldb,
			 void *y, int incy);
void BLAS_zge_sum_mv_z_d(enum blas_order_type order, int m, int n,
			 const void *alpha, const void *a, int lda,
			 const double *x, int incx,
			 const void *beta, const void *b, int ldb,
			 void *y, int incy);
void BLAS_zge_sum_mv_d_z(enum blas_order_type order, int m, int n,
			 const void *alpha, const double *a, int lda,
			 const void *x, int incx,
			 const void *beta, const double *b, int ldb,
			 void *y, int incy);
void BLAS_zge_sum_mv_d_d(enum blas_order_type order, int m, int n,
			 const void *alpha, const double *a, int lda,
			 const double *x, int incx,
			 const void *beta, const double *b, int ldb,
			 void *y, int incy);
void BLAS_sge_sum_mv_x(enum blas_order_type order, int m, int n,
		       float alpha, const float *a, int lda,
		       const float *x, int incx,
		       float beta, const float *b, int ldb,
		       float *y, int incy, enum blas_prec_type prec);
void BLAS_dge_sum_mv_x(enum blas_order_type order, int m, int n,
		       double alpha, const double *a, int lda,
		       const double *x, int incx,
		       double beta, const double *b, int ldb,
		       double *y, int incy, enum blas_prec_type prec);
void BLAS_cge_sum_mv_x(enum blas_order_type order, int m, int n,
		       const void *alpha, const void *a, int lda,
		       const void *x, int incx,
		       const void *beta, const void *b, int ldb,
		       void *y, int incy, enum blas_prec_type prec);
void BLAS_zge_sum_mv_x(enum blas_order_type order, int m, int n,
		       const void *alpha, const void *a, int lda,
		       const void *x, int incx,
		       const void *beta, const void *b, int ldb,
		       void *y, int incy, enum blas_prec_type prec);
void BLAS_dge_sum_mv_d_s_x(enum blas_order_type order, int m, int n,
			   double alpha, const double *a, int lda,
			   const float *x, int incx,
			   double beta, const double *b, int ldb,
			   double *y, int incy, enum blas_prec_type prec);
void BLAS_dge_sum_mv_s_d_x(enum blas_order_type order, int m, int n,
			   double alpha, const float *a, int lda,
			   const double *x, int incx,
			   double beta, const float *b, int ldb,
			   double *y, int incy, enum blas_prec_type prec);
void BLAS_dge_sum_mv_s_s_x(enum blas_order_type order, int m, int n,
			   double alpha, const float *a, int lda,
			   const float *x, int incx,
			   double beta, const float *b, int ldb,
			   double *y, int incy, enum blas_prec_type prec);
void BLAS_zge_sum_mv_z_c_x(enum blas_order_type order, int m, int n,
			   const void *alpha, const void *a, int lda,
			   const void *x, int incx,
			   const void *beta, const void *b, int ldb,
			   void *y, int incy, enum blas_prec_type prec);
void BLAS_zge_sum_mv_c_z_x(enum blas_order_type order, int m, int n,
			   const void *alpha, const void *a, int lda,
			   const void *x, int incx,
			   const void *beta, const void *b, int ldb,
			   void *y, int incy, enum blas_prec_type prec);
void BLAS_zge_sum_mv_c_c_x(enum blas_order_type order, int m, int n,
			   const void *alpha, const void *a, int lda,
			   const void *x, int incx,
			   const void *beta, const void *b, int ldb,
			   void *y, int incy, enum blas_prec_type prec);
void BLAS_cge_sum_mv_c_s_x(enum blas_order_type order, int m, int n,
			   const void *alpha, const void *a, int lda,
			   const float *x, int incx,
			   const void *beta, const void *b, int ldb,
			   void *y, int incy, enum blas_prec_type prec);
void BLAS_cge_sum_mv_s_c_x(enum blas_order_type order, int m, int n,
			   const void *alpha, const float *a, int lda,
			   const void *x, int incx,
			   const void *beta, const float *b, int ldb,
			   void *y, int incy, enum blas_prec_type prec);
void BLAS_cge_sum_mv_s_s_x(enum blas_order_type order, int m, int n,
			   const void *alpha, const float *a, int lda,
			   const float *x, int incx,
			   const void *beta, const float *b, int ldb,
			   void *y, int incy, enum blas_prec_type prec);
void BLAS_zge_sum_mv_z_d_x(enum blas_order_type order, int m, int n,
			   const void *alpha, const void *a, int lda,
			   const double *x, int incx,
			   const void *beta, const void *b, int ldb,
			   void *y, int incy, enum blas_prec_type prec);
void BLAS_zge_sum_mv_d_z_x(enum blas_order_type order, int m, int n,
			   const void *alpha, const double *a, int lda,
			   const void *x, int incx,
			   const void *beta, const double *b, int ldb,
			   void *y, int incy, enum blas_prec_type prec);
void BLAS_zge_sum_mv_d_d_x(enum blas_order_type order, int m, int n,
			   const void *alpha, const double *a, int lda,
			   const double *x, int incx,
			   const void *beta, const double *b, int ldb,
			   void *y, int incy, enum blas_prec_type prec);


void BLAS_dgbmv_d_s(enum blas_order_type order, enum blas_trans_type trans,
		    int m, int n, int kl, int ku, double alpha,
		    const double *a, int lda, const float *x, int incx,
		    double beta, double *y, int incy);
void BLAS_dgbmv_s_d(enum blas_order_type order, enum blas_trans_type trans,
		    int m, int n, int kl, int ku, double alpha,
		    const float *a, int lda, const double *x, int incx,
		    double beta, double *y, int incy);
void BLAS_dgbmv_s_s(enum blas_order_type order, enum blas_trans_type trans,
		    int m, int n, int kl, int ku, double alpha,
		    const float *a, int lda, const float *x, int incx,
		    double beta, double *y, int incy);
void BLAS_zgbmv_z_c(enum blas_order_type order, enum blas_trans_type trans,
		    int m, int n, int kl, int ku, const void *alpha,
		    const void *a, int lda, const void *x, int incx,
		    const void *beta, void *y, int incy);
void BLAS_zgbmv_c_z(enum blas_order_type order, enum blas_trans_type trans,
		    int m, int n, int kl, int ku, const void *alpha,
		    const void *a, int lda, const void *x, int incx,
		    const void *beta, void *y, int incy);
void BLAS_zgbmv_c_c(enum blas_order_type order, enum blas_trans_type trans,
		    int m, int n, int kl, int ku, const void *alpha,
		    const void *a, int lda, const void *x, int incx,
		    const void *beta, void *y, int incy);
void BLAS_cgbmv_c_s(enum blas_order_type order, enum blas_trans_type trans,
		    int m, int n, int kl, int ku, const void *alpha,
		    const void *a, int lda, const float *x, int incx,
		    const void *beta, void *y, int incy);
void BLAS_cgbmv_s_c(enum blas_order_type order, enum blas_trans_type trans,
		    int m, int n, int kl, int ku, const void *alpha,
		    const float *a, int lda, const void *x, int incx,
		    const void *beta, void *y, int incy);
void BLAS_cgbmv_s_s(enum blas_order_type order, enum blas_trans_type trans,
		    int m, int n, int kl, int ku, const void *alpha,
		    const float *a, int lda, const float *x, int incx,
		    const void *beta, void *y, int incy);
void BLAS_zgbmv_z_d(enum blas_order_type order, enum blas_trans_type trans,
		    int m, int n, int kl, int ku, const void *alpha,
		    const void *a, int lda, const double *x, int incx,
		    const void *beta, void *y, int incy);
void BLAS_zgbmv_d_z(enum blas_order_type order, enum blas_trans_type trans,
		    int m, int n, int kl, int ku, const void *alpha,
		    const double *a, int lda, const void *x, int incx,
		    const void *beta, void *y, int incy);
void BLAS_zgbmv_d_d(enum blas_order_type order, enum blas_trans_type trans,
		    int m, int n, int kl, int ku, const void *alpha,
		    const double *a, int lda, const double *x, int incx,
		    const void *beta, void *y, int incy);
void BLAS_sgbmv_x(enum blas_order_type order, enum blas_trans_type trans,
		  int m, int n, int kl, int ku, float alpha,
		  const float *a, int lda, const float *x, int incx,
		  float beta, float *y, int incy, enum blas_prec_type prec);
void BLAS_dgbmv_x(enum blas_order_type order, enum blas_trans_type trans,
		  int m, int n, int kl, int ku, double alpha,
		  const double *a, int lda, const double *x, int incx,
		  double beta, double *y, int incy, enum blas_prec_type prec);
void BLAS_cgbmv_x(enum blas_order_type order, enum blas_trans_type trans,
		  int m, int n, int kl, int ku, const void *alpha,
		  const void *a, int lda, const void *x, int incx,
		  const void *beta, void *y, int incy,
		  enum blas_prec_type prec);
void BLAS_zgbmv_x(enum blas_order_type order, enum blas_trans_type trans,
		  int m, int n, int kl, int ku, const void *alpha,
		  const void *a, int lda, const void *x, int incx,
		  const void *beta, void *y, int incy,
		  enum blas_prec_type prec);
void BLAS_dgbmv_d_s_x(enum blas_order_type order, enum blas_trans_type trans,
		      int m, int n, int kl, int ku, double alpha,
		      const double *a, int lda, const float *x, int incx,
		      double beta, double *y, int incy,
		      enum blas_prec_type prec);
void BLAS_dgbmv_s_d_x(enum blas_order_type order, enum blas_trans_type trans,
		      int m, int n, int kl, int ku, double alpha,
		      const float *a, int lda, const double *x, int incx,
		      double beta, double *y, int incy,
		      enum blas_prec_type prec);
void BLAS_dgbmv_s_s_x(enum blas_order_type order, enum blas_trans_type trans,
		      int m, int n, int kl, int ku, double alpha,
		      const float *a, int lda, const float *x, int incx,
		      double beta, double *y, int incy,
		      enum blas_prec_type prec);
void BLAS_zgbmv_z_c_x(enum blas_order_type order, enum blas_trans_type trans,
		      int m, int n, int kl, int ku, const void *alpha,
		      const void *a, int lda, const void *x, int incx,
		      const void *beta, void *y, int incy,
		      enum blas_prec_type prec);
void BLAS_zgbmv_c_z_x(enum blas_order_type order, enum blas_trans_type trans,
		      int m, int n, int kl, int ku, const void *alpha,
		      const void *a, int lda, const void *x, int incx,
		      const void *beta, void *y, int incy,
		      enum blas_prec_type prec);
void BLAS_zgbmv_c_c_x(enum blas_order_type order, enum blas_trans_type trans,
		      int m, int n, int kl, int ku, const void *alpha,
		      const void *a, int lda, const void *x, int incx,
		      const void *beta, void *y, int incy,
		      enum blas_prec_type prec);
void BLAS_cgbmv_c_s_x(enum blas_order_type order, enum blas_trans_type trans,
		      int m, int n, int kl, int ku, const void *alpha,
		      const void *a, int lda, const float *x, int incx,
		      const void *beta, void *y, int incy,
		      enum blas_prec_type prec);
void BLAS_cgbmv_s_c_x(enum blas_order_type order, enum blas_trans_type trans,
		      int m, int n, int kl, int ku, const void *alpha,
		      const float *a, int lda, const void *x, int incx,
		      const void *beta, void *y, int incy,
		      enum blas_prec_type prec);
void BLAS_cgbmv_s_s_x(enum blas_order_type order, enum blas_trans_type trans,
		      int m, int n, int kl, int ku, const void *alpha,
		      const float *a, int lda, const float *x, int incx,
		      const void *beta, void *y, int incy,
		      enum blas_prec_type prec);
void BLAS_zgbmv_z_d_x(enum blas_order_type order, enum blas_trans_type trans,
		      int m, int n, int kl, int ku, const void *alpha,
		      const void *a, int lda, const double *x, int incx,
		      const void *beta, void *y, int incy,
		      enum blas_prec_type prec);
void BLAS_zgbmv_d_z_x(enum blas_order_type order, enum blas_trans_type trans,
		      int m, int n, int kl, int ku, const void *alpha,
		      const double *a, int lda, const void *x, int incx,
		      const void *beta, void *y, int incy,
		      enum blas_prec_type prec);
void BLAS_zgbmv_d_d_x(enum blas_order_type order, enum blas_trans_type trans,
		      int m, int n, int kl, int ku, const void *alpha,
		      const double *a, int lda, const double *x, int incx,
		      const void *beta, void *y, int incy,
		      enum blas_prec_type prec);


void BLAS_dsymv_d_s(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, double alpha, const double *a, int lda,
		    const float *x, int incx, double beta,
		    double *y, int incy);
void BLAS_dsymv_s_d(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, double alpha, const float *a, int lda,
		    const double *x, int incx, double beta,
		    double *y, int incy);
void BLAS_dsymv_s_s(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, double alpha, const float *a, int lda,
		    const float *x, int incx, double beta,
		    double *y, int incy);
void BLAS_zsymv_z_c(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, const void *alpha, const void *a, int lda,
		    const void *x, int incx, const void *beta,
		    void *y, int incy);
void BLAS_zsymv_c_z(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, const void *alpha, const void *a, int lda,
		    const void *x, int incx, const void *beta,
		    void *y, int incy);
void BLAS_zsymv_c_c(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, const void *alpha, const void *a, int lda,
		    const void *x, int incx, const void *beta,
		    void *y, int incy);
void BLAS_csymv_c_s(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, const void *alpha, const void *a, int lda,
		    const float *x, int incx, const void *beta,
		    void *y, int incy);
void BLAS_csymv_s_c(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, const void *alpha, const float *a, int lda,
		    const void *x, int incx, const void *beta,
		    void *y, int incy);
void BLAS_csymv_s_s(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, const void *alpha, const float *a, int lda,
		    const float *x, int incx, const void *beta,
		    void *y, int incy);
void BLAS_zsymv_z_d(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, const void *alpha, const void *a, int lda,
		    const double *x, int incx, const void *beta,
		    void *y, int incy);
void BLAS_zsymv_d_z(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, const void *alpha, const double *a, int lda,
		    const void *x, int incx, const void *beta,
		    void *y, int incy);
void BLAS_zsymv_d_d(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, const void *alpha, const double *a, int lda,
		    const double *x, int incx, const void *beta,
		    void *y, int incy);
void BLAS_ssymv_x(enum blas_order_type order, enum blas_uplo_type uplo,
		  int n, float alpha, const float *a, int lda,
		  const float *x, int incx, float beta,
		  float *y, int incy, enum blas_prec_type prec);
void BLAS_dsymv_x(enum blas_order_type order, enum blas_uplo_type uplo,
		  int n, double alpha, const double *a, int lda,
		  const double *x, int incx, double beta,
		  double *y, int incy, enum blas_prec_type prec);
void BLAS_csymv_x(enum blas_order_type order, enum blas_uplo_type uplo,
		  int n, const void *alpha, const void *a, int lda,
		  const void *x, int incx, const void *beta,
		  void *y, int incy, enum blas_prec_type prec);
void BLAS_zsymv_x(enum blas_order_type order, enum blas_uplo_type uplo,
		  int n, const void *alpha, const void *a, int lda,
		  const void *x, int incx, const void *beta,
		  void *y, int incy, enum blas_prec_type prec);
void BLAS_dsymv_d_s_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, double alpha, const double *a, int lda,
		      const float *x, int incx, double beta,
		      double *y, int incy, enum blas_prec_type prec);
void BLAS_dsymv_s_d_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, double alpha, const float *a, int lda,
		      const double *x, int incx, double beta,
		      double *y, int incy, enum blas_prec_type prec);
void BLAS_dsymv_s_s_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, double alpha, const float *a, int lda,
		      const float *x, int incx, double beta,
		      double *y, int incy, enum blas_prec_type prec);
void BLAS_zsymv_z_c_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, const void *alpha, const void *a, int lda,
		      const void *x, int incx, const void *beta,
		      void *y, int incy, enum blas_prec_type prec);
void BLAS_zsymv_c_z_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, const void *alpha, const void *a, int lda,
		      const void *x, int incx, const void *beta,
		      void *y, int incy, enum blas_prec_type prec);
void BLAS_zsymv_c_c_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, const void *alpha, const void *a, int lda,
		      const void *x, int incx, const void *beta,
		      void *y, int incy, enum blas_prec_type prec);
void BLAS_csymv_c_s_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, const void *alpha, const void *a, int lda,
		      const float *x, int incx, const void *beta,
		      void *y, int incy, enum blas_prec_type prec);
void BLAS_csymv_s_c_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, const void *alpha, const float *a, int lda,
		      const void *x, int incx, const void *beta,
		      void *y, int incy, enum blas_prec_type prec);
void BLAS_csymv_s_s_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, const void *alpha, const float *a, int lda,
		      const float *x, int incx, const void *beta,
		      void *y, int incy, enum blas_prec_type prec);
void BLAS_zsymv_z_d_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, const void *alpha, const void *a, int lda,
		      const double *x, int incx, const void *beta,
		      void *y, int incy, enum blas_prec_type prec);
void BLAS_zsymv_d_z_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, const void *alpha, const double *a, int lda,
		      const void *x, int incx, const void *beta,
		      void *y, int incy, enum blas_prec_type prec);
void BLAS_zsymv_d_d_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, const void *alpha, const double *a, int lda,
		      const double *x, int incx, const void *beta,
		      void *y, int incy, enum blas_prec_type prec);


void BLAS_dspmv_d_s(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, double alpha, const double *ap,
		    const float *x, int incx, double beta,
		    double *y, int incy);
void BLAS_dspmv_s_d(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, double alpha, const float *ap,
		    const double *x, int incx, double beta,
		    double *y, int incy);
void BLAS_dspmv_s_s(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, double alpha, const float *ap,
		    const float *x, int incx, double beta,
		    double *y, int incy);
void BLAS_zspmv_z_c(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, const void *alpha, const void *ap,
		    const void *x, int incx, const void *beta,
		    void *y, int incy);
void BLAS_zspmv_c_z(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, const void *alpha, const void *ap,
		    const void *x, int incx, const void *beta,
		    void *y, int incy);
void BLAS_zspmv_c_c(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, const void *alpha, const void *ap,
		    const void *x, int incx, const void *beta,
		    void *y, int incy);
void BLAS_cspmv_c_s(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, const void *alpha, const void *ap,
		    const float *x, int incx, const void *beta,
		    void *y, int incy);
void BLAS_cspmv_s_c(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, const void *alpha, const float *ap,
		    const void *x, int incx, const void *beta,
		    void *y, int incy);
void BLAS_cspmv_s_s(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, const void *alpha, const float *ap,
		    const float *x, int incx, const void *beta,
		    void *y, int incy);
void BLAS_zspmv_z_d(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, const void *alpha, const void *ap,
		    const double *x, int incx, const void *beta,
		    void *y, int incy);
void BLAS_zspmv_d_z(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, const void *alpha, const double *ap,
		    const void *x, int incx, const void *beta,
		    void *y, int incy);
void BLAS_zspmv_d_d(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, const void *alpha, const double *ap,
		    const double *x, int incx, const void *beta,
		    void *y, int incy);
void BLAS_sspmv_x(enum blas_order_type order, enum blas_uplo_type uplo,
		  int n, float alpha, const float *ap,
		  const float *x, int incx, float beta,
		  float *y, int incy, enum blas_prec_type prec);
void BLAS_dspmv_x(enum blas_order_type order, enum blas_uplo_type uplo,
		  int n, double alpha, const double *ap,
		  const double *x, int incx, double beta,
		  double *y, int incy, enum blas_prec_type prec);
void BLAS_cspmv_x(enum blas_order_type order, enum blas_uplo_type uplo,
		  int n, const void *alpha, const void *ap,
		  const void *x, int incx, const void *beta,
		  void *y, int incy, enum blas_prec_type prec);
void BLAS_zspmv_x(enum blas_order_type order, enum blas_uplo_type uplo,
		  int n, const void *alpha, const void *ap,
		  const void *x, int incx, const void *beta,
		  void *y, int incy, enum blas_prec_type prec);
void BLAS_dspmv_d_s_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, double alpha, const double *ap,
		      const float *x, int incx, double beta,
		      double *y, int incy, enum blas_prec_type prec);
void BLAS_dspmv_s_d_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, double alpha, const float *ap,
		      const double *x, int incx, double beta,
		      double *y, int incy, enum blas_prec_type prec);
void BLAS_dspmv_s_s_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, double alpha, const float *ap,
		      const float *x, int incx, double beta,
		      double *y, int incy, enum blas_prec_type prec);
void BLAS_zspmv_z_c_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, const void *alpha, const void *ap,
		      const void *x, int incx, const void *beta,
		      void *y, int incy, enum blas_prec_type prec);
void BLAS_zspmv_c_z_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, const void *alpha, const void *ap,
		      const void *x, int incx, const void *beta,
		      void *y, int incy, enum blas_prec_type prec);
void BLAS_zspmv_c_c_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, const void *alpha, const void *ap,
		      const void *x, int incx, const void *beta,
		      void *y, int incy, enum blas_prec_type prec);
void BLAS_cspmv_c_s_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, const void *alpha, const void *ap,
		      const float *x, int incx, const void *beta,
		      void *y, int incy, enum blas_prec_type prec);
void BLAS_cspmv_s_c_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, const void *alpha, const float *ap,
		      const void *x, int incx, const void *beta,
		      void *y, int incy, enum blas_prec_type prec);
void BLAS_cspmv_s_s_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, const void *alpha, const float *ap,
		      const float *x, int incx, const void *beta,
		      void *y, int incy, enum blas_prec_type prec);
void BLAS_zspmv_z_d_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, const void *alpha, const void *ap,
		      const double *x, int incx, const void *beta,
		      void *y, int incy, enum blas_prec_type prec);
void BLAS_zspmv_d_z_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, const void *alpha, const double *ap,
		      const void *x, int incx, const void *beta,
		      void *y, int incy, enum blas_prec_type prec);
void BLAS_zspmv_d_d_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, const void *alpha, const double *ap,
		      const double *x, int incx, const void *beta,
		      void *y, int incy, enum blas_prec_type prec);


void BLAS_dsbmv_d_s(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, int k, double alpha, const double *a, int lda,
		    const float *x, int incx, double beta,
		    double *y, int incy);
void BLAS_dsbmv_s_d(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, int k, double alpha, const float *a, int lda,
		    const double *x, int incx, double beta,
		    double *y, int incy);
void BLAS_dsbmv_s_s(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, int k, double alpha, const float *a, int lda,
		    const float *x, int incx, double beta,
		    double *y, int incy);
void BLAS_zsbmv_z_c(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, int k, const void *alpha, const void *a, int lda,
		    const void *x, int incx, const void *beta,
		    void *y, int incy);
void BLAS_zsbmv_c_z(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, int k, const void *alpha, const void *a, int lda,
		    const void *x, int incx, const void *beta,
		    void *y, int incy);
void BLAS_zsbmv_c_c(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, int k, const void *alpha, const void *a, int lda,
		    const void *x, int incx, const void *beta,
		    void *y, int incy);
void BLAS_csbmv_c_s(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, int k, const void *alpha, const void *a, int lda,
		    const float *x, int incx, const void *beta,
		    void *y, int incy);
void BLAS_csbmv_s_c(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, int k, const void *alpha, const float *a, int lda,
		    const void *x, int incx, const void *beta,
		    void *y, int incy);
void BLAS_csbmv_s_s(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, int k, const void *alpha, const float *a, int lda,
		    const float *x, int incx, const void *beta,
		    void *y, int incy);
void BLAS_zsbmv_z_d(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, int k, const void *alpha, const void *a, int lda,
		    const double *x, int incx, const void *beta,
		    void *y, int incy);
void BLAS_zsbmv_d_z(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, int k, const void *alpha, const double *a, int lda,
		    const void *x, int incx, const void *beta,
		    void *y, int incy);
void BLAS_zsbmv_d_d(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, int k, const void *alpha, const double *a, int lda,
		    const double *x, int incx, const void *beta,
		    void *y, int incy);
void BLAS_ssbmv_x(enum blas_order_type order, enum blas_uplo_type uplo,
		  int n, int k, float alpha, const float *a, int lda,
		  const float *x, int incx, float beta,
		  float *y, int incy, enum blas_prec_type prec);
void BLAS_dsbmv_x(enum blas_order_type order, enum blas_uplo_type uplo,
		  int n, int k, double alpha, const double *a, int lda,
		  const double *x, int incx, double beta,
		  double *y, int incy, enum blas_prec_type prec);
void BLAS_csbmv_x(enum blas_order_type order, enum blas_uplo_type uplo,
		  int n, int k, const void *alpha, const void *a, int lda,
		  const void *x, int incx, const void *beta,
		  void *y, int incy, enum blas_prec_type prec);
void BLAS_zsbmv_x(enum blas_order_type order, enum blas_uplo_type uplo,
		  int n, int k, const void *alpha, const void *a, int lda,
		  const void *x, int incx, const void *beta,
		  void *y, int incy, enum blas_prec_type prec);
void BLAS_dsbmv_d_s_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, int k, double alpha, const double *a, int lda,
		      const float *x, int incx, double beta,
		      double *y, int incy, enum blas_prec_type prec);
void BLAS_dsbmv_s_d_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, int k, double alpha, const float *a, int lda,
		      const double *x, int incx, double beta,
		      double *y, int incy, enum blas_prec_type prec);
void BLAS_dsbmv_s_s_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, int k, double alpha, const float *a, int lda,
		      const float *x, int incx, double beta,
		      double *y, int incy, enum blas_prec_type prec);
void BLAS_zsbmv_z_c_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, int k, const void *alpha, const void *a, int lda,
		      const void *x, int incx, const void *beta,
		      void *y, int incy, enum blas_prec_type prec);
void BLAS_zsbmv_c_z_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, int k, const void *alpha, const void *a, int lda,
		      const void *x, int incx, const void *beta,
		      void *y, int incy, enum blas_prec_type prec);
void BLAS_zsbmv_c_c_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, int k, const void *alpha, const void *a, int lda,
		      const void *x, int incx, const void *beta,
		      void *y, int incy, enum blas_prec_type prec);
void BLAS_csbmv_c_s_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, int k, const void *alpha, const void *a, int lda,
		      const float *x, int incx, const void *beta,
		      void *y, int incy, enum blas_prec_type prec);
void BLAS_csbmv_s_c_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, int k, const void *alpha, const float *a,
		      int lda, const void *x, int incx, const void *beta,
		      void *y, int incy, enum blas_prec_type prec);
void BLAS_csbmv_s_s_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, int k, const void *alpha, const float *a,
		      int lda, const float *x, int incx, const void *beta,
		      void *y, int incy, enum blas_prec_type prec);
void BLAS_zsbmv_z_d_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, int k, const void *alpha, const void *a, int lda,
		      const double *x, int incx, const void *beta, void *y,
		      int incy, enum blas_prec_type prec);
void BLAS_zsbmv_d_z_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, int k, const void *alpha, const double *a,
		      int lda, const void *x, int incx, const void *beta,
		      void *y, int incy, enum blas_prec_type prec);
void BLAS_zsbmv_d_d_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, int k, const void *alpha, const double *a,
		      int lda, const double *x, int incx, const void *beta,
		      void *y, int incy, enum blas_prec_type prec);


void BLAS_zhemv_z_c(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, const void *alpha, const void *a, int lda,
		    const void *x, int incx, const void *beta,
		    void *y, int incy);
void BLAS_zhemv_c_z(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, const void *alpha, const void *a, int lda,
		    const void *x, int incx, const void *beta,
		    void *y, int incy);
void BLAS_zhemv_c_c(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, const void *alpha, const void *a, int lda,
		    const void *x, int incx, const void *beta,
		    void *y, int incy);
void BLAS_chemv_c_s(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, const void *alpha, const void *a, int lda,
		    const float *x, int incx, const void *beta,
		    void *y, int incy);
void BLAS_zhemv_z_d(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, const void *alpha, const void *a, int lda,
		    const double *x, int incx, const void *beta,
		    void *y, int incy);
void BLAS_chemv_x(enum blas_order_type order, enum blas_uplo_type uplo,
		  int n, const void *alpha, const void *a, int lda,
		  const void *x, int incx, const void *beta,
		  void *y, int incy, enum blas_prec_type prec);
void BLAS_zhemv_x(enum blas_order_type order, enum blas_uplo_type uplo,
		  int n, const void *alpha, const void *a, int lda,
		  const void *x, int incx, const void *beta,
		  void *y, int incy, enum blas_prec_type prec);
void BLAS_zhemv_z_c_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, const void *alpha, const void *a, int lda,
		      const void *x, int incx, const void *beta,
		      void *y, int incy, enum blas_prec_type prec);
void BLAS_zhemv_c_z_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, const void *alpha, const void *a, int lda,
		      const void *x, int incx, const void *beta,
		      void *y, int incy, enum blas_prec_type prec);
void BLAS_zhemv_c_c_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, const void *alpha, const void *a, int lda,
		      const void *x, int incx, const void *beta,
		      void *y, int incy, enum blas_prec_type prec);
void BLAS_chemv_c_s_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, const void *alpha, const void *a, int lda,
		      const float *x, int incx, const void *beta,
		      void *y, int incy, enum blas_prec_type prec);
void BLAS_zhemv_z_d_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, const void *alpha, const void *a, int lda,
		      const double *x, int incx, const void *beta,
		      void *y, int incy, enum blas_prec_type prec);


void BLAS_zhpmv_z_c(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, const void *alpha, const void *ap,
		    const void *x, int incx, const void *beta, void *y,
		    int incy);
void BLAS_zhpmv_c_z(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, const void *alpha, const void *ap,
		    const void *x, int incx, const void *beta, void *y,
		    int incy);
void BLAS_zhpmv_c_c(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, const void *alpha, const void *ap,
		    const void *x, int incx, const void *beta, void *y,
		    int incy);
void BLAS_chpmv_c_s(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, const void *alpha, const void *ap,
		    const float *x, int incx, const void *beta, void *y,
		    int incy);
void BLAS_zhpmv_z_d(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, const void *alpha, const void *ap,
		    const double *x, int incx, const void *beta, void *y,
		    int incy);
void BLAS_chpmv_x(enum blas_order_type order, enum blas_uplo_type uplo,
		  int n, const void *alpha, const void *ap,
		  const void *x, int incx, const void *beta, void *y,
		  int incy, enum blas_prec_type prec);
void BLAS_zhpmv_x(enum blas_order_type order, enum blas_uplo_type uplo,
		  int n, const void *alpha, const void *ap,
		  const void *x, int incx, const void *beta, void *y,
		  int incy, enum blas_prec_type prec);
void BLAS_zhpmv_z_c_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, const void *alpha, const void *ap,
		      const void *x, int incx, const void *beta, void *y,
		      int incy, enum blas_prec_type prec);
void BLAS_zhpmv_c_z_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, const void *alpha, const void *ap,
		      const void *x, int incx, const void *beta, void *y,
		      int incy, enum blas_prec_type prec);
void BLAS_zhpmv_c_c_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, const void *alpha, const void *ap,
		      const void *x, int incx, const void *beta, void *y,
		      int incy, enum blas_prec_type prec);
void BLAS_chpmv_c_s_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, const void *alpha, const void *ap,
		      const float *x, int incx, const void *beta, void *y,
		      int incy, enum blas_prec_type prec);
void BLAS_zhpmv_z_d_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, const void *alpha, const void *ap,
		      const double *x, int incx, const void *beta, void *y,
		      int incy, enum blas_prec_type prec);


void BLAS_zhbmv_z_c(enum blas_order_type order,
		    enum blas_uplo_type uplo, int n, int k,
		    const void *alpha, const void *a, int lda,
		    const void *x, int incx, const void *beta,
		    void *y, int incy);
void BLAS_zhbmv_c_z(enum blas_order_type order,
		    enum blas_uplo_type uplo, int n, int k,
		    const void *alpha, const void *a, int lda,
		    const void *x, int incx, const void *beta,
		    void *y, int incy);
void BLAS_zhbmv_c_c(enum blas_order_type order,
		    enum blas_uplo_type uplo, int n, int k,
		    const void *alpha, const void *a, int lda,
		    const void *x, int incx, const void *beta,
		    void *y, int incy);
void BLAS_chbmv_c_s(enum blas_order_type order,
		    enum blas_uplo_type uplo, int n, int k,
		    const void *alpha, const void *a, int lda,
		    const float *x, int incx, const void *beta,
		    void *y, int incy);
void BLAS_zhbmv_z_d(enum blas_order_type order,
		    enum blas_uplo_type uplo, int n, int k,
		    const void *alpha, const void *a, int lda,
		    const double *x, int incx, const void *beta,
		    void *y, int incy);
void BLAS_chbmv_x(enum blas_order_type order,
		  enum blas_uplo_type uplo, int n, int k,
		  const void *alpha, const void *a, int lda,
		  const void *x, int incx, const void *beta,
		  void *y, int incy, enum blas_prec_type prec);
void BLAS_zhbmv_x(enum blas_order_type order,
		  enum blas_uplo_type uplo, int n, int k,
		  const void *alpha, const void *a, int lda,
		  const void *x, int incx, const void *beta,
		  void *y, int incy, enum blas_prec_type prec);
void BLAS_zhbmv_z_c_x(enum blas_order_type order,
		      enum blas_uplo_type uplo, int n, int k,
		      const void *alpha, const void *a, int lda,
		      const void *x, int incx, const void *beta,
		      void *y, int incy, enum blas_prec_type prec);
void BLAS_zhbmv_c_z_x(enum blas_order_type order,
		      enum blas_uplo_type uplo, int n, int k,
		      const void *alpha, const void *a, int lda,
		      const void *x, int incx, const void *beta,
		      void *y, int incy, enum blas_prec_type prec);
void BLAS_zhbmv_c_c_x(enum blas_order_type order,
		      enum blas_uplo_type uplo, int n, int k,
		      const void *alpha, const void *a, int lda,
		      const void *x, int incx, const void *beta,
		      void *y, int incy, enum blas_prec_type prec);
void BLAS_chbmv_c_s_x(enum blas_order_type order,
		      enum blas_uplo_type uplo, int n, int k,
		      const void *alpha, const void *a, int lda,
		      const float *x, int incx, const void *beta,
		      void *y, int incy, enum blas_prec_type prec);
void BLAS_zhbmv_z_d_x(enum blas_order_type order,
		      enum blas_uplo_type uplo, int n, int k,
		      const void *alpha, const void *a, int lda,
		      const double *x, int incx, const void *beta,
		      void *y, int incy, enum blas_prec_type prec);


void BLAS_dtrmv_s(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag, int n,
		  double alpha, const float *T, int ldt, double *x, int incx);
void BLAS_ztrmv_c(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag, int n,
		  const void *alpha, const void *T, int ldt,
		  void *x, int incx);
void BLAS_ctrmv_s(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag, int n,
		  const void *alpha, const float *T, int ldt,
		  void *x, int incx);
void BLAS_ztrmv_d(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag, int n,
		  const void *alpha, const double *T, int ldt,
		  void *x, int incx);
void BLAS_strmv_x(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag, int n,
		  float alpha, const float *T, int ldt,
		  float *x, int incx, enum blas_prec_type prec);
void BLAS_dtrmv_x(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag, int n,
		  double alpha, const double *T, int ldt,
		  double *x, int incx, enum blas_prec_type prec);
void BLAS_ctrmv_x(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag, int n,
		  const void *alpha, const void *T, int ldt,
		  void *x, int incx, enum blas_prec_type prec);
void BLAS_ztrmv_x(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag, int n,
		  const void *alpha, const void *T, int ldt,
		  void *x, int incx, enum blas_prec_type prec);
void BLAS_dtrmv_s_x(enum blas_order_type order, enum blas_uplo_type uplo,
		    enum blas_trans_type trans, enum blas_diag_type diag,
		    int n, double alpha, const float *T, int ldt, double *x,
		    int incx, enum blas_prec_type prec);
void BLAS_ztrmv_c_x(enum blas_order_type order, enum blas_uplo_type uplo,
		    enum blas_trans_type trans, enum blas_diag_type diag,
		    int n, const void *alpha, const void *T, int ldt, void *x,
		    int incx, enum blas_prec_type prec);
void BLAS_ctrmv_s_x(enum blas_order_type order, enum blas_uplo_type uplo,
		    enum blas_trans_type trans, enum blas_diag_type diag,
		    int n, const void *alpha, const float *T, int ldt,
		    void *x, int incx, enum blas_prec_type prec);
void BLAS_ztrmv_d_x(enum blas_order_type order, enum blas_uplo_type uplo,
		    enum blas_trans_type trans, enum blas_diag_type diag,
		    int n, const void *alpha, const double *T, int ldt,
		    void *x, int incx, enum blas_prec_type prec);


void BLAS_dtpmv_s(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag,
		  int n, double alpha, const float *tp, double *x, int incx);
void BLAS_ztpmv_c(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag,
		  int n, const void *alpha, const void *tp,
		  void *x, int incx);
void BLAS_ctpmv_s(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag,
		  int n, const void *alpha, const float *tp,
		  void *x, int incx);
void BLAS_ztpmv_d(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag,
		  int n, const void *alpha, const double *tp,
		  void *x, int incx);
void BLAS_stpmv_x(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag,
		  int n, float alpha, const float *tp,
		  float *x, int incx, enum blas_prec_type prec);
void BLAS_dtpmv_x(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag,
		  int n, double alpha, const double *tp,
		  double *x, int incx, enum blas_prec_type prec);
void BLAS_ctpmv_x(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag,
		  int n, const void *alpha, const void *tp,
		  void *x, int incx, enum blas_prec_type prec);
void BLAS_ztpmv_x(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag,
		  int n, const void *alpha, const void *tp,
		  void *x, int incx, enum blas_prec_type prec);
void BLAS_dtpmv_s_x(enum blas_order_type order, enum blas_uplo_type uplo,
		    enum blas_trans_type trans, enum blas_diag_type diag,
		    int n, double alpha, const float *tp,
		    double *x, int incx, enum blas_prec_type prec);
void BLAS_ztpmv_c_x(enum blas_order_type order, enum blas_uplo_type uplo,
		    enum blas_trans_type trans, enum blas_diag_type diag,
		    int n, const void *alpha, const void *tp,
		    void *x, int incx, enum blas_prec_type prec);
void BLAS_ctpmv_s_x(enum blas_order_type order, enum blas_uplo_type uplo,
		    enum blas_trans_type trans, enum blas_diag_type diag,
		    int n, const void *alpha, const float *tp,
		    void *x, int incx, enum blas_prec_type prec);
void BLAS_ztpmv_d_x(enum blas_order_type order, enum blas_uplo_type uplo,
		    enum blas_trans_type trans, enum blas_diag_type diag,
		    int n, const void *alpha, const double *tp,
		    void *x, int incx, enum blas_prec_type prec);


void BLAS_dtrsv_s(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag,
		  int n, double alpha, const float *T, int ldt,
		  double *x, int incx);
void BLAS_ztrsv_c(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag,
		  int n, const void *alpha, const void *T, int ldt,
		  void *x, int incx);
void BLAS_ctrsv_s(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag,
		  int n, const void *alpha, const float *T, int ldt,
		  void *x, int incx);
void BLAS_ztrsv_d(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag,
		  int n, const void *alpha, const double *T, int ldt,
		  void *x, int incx);
void BLAS_strsv_x(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag,
		  int n, float alpha, const float *T, int ldt,
		  float *x, int incx, enum blas_prec_type prec);
void BLAS_dtrsv_x(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag,
		  int n, double alpha, const double *T, int ldt,
		  double *x, int incx, enum blas_prec_type prec);
void BLAS_dtrsv_s_x(enum blas_order_type order, enum blas_uplo_type uplo,
		    enum blas_trans_type trans, enum blas_diag_type diag,
		    int n, double alpha, const float *T, int ldt,
		    double *x, int incx, enum blas_prec_type prec);
void BLAS_ctrsv_x(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag,
		  int n, const void *alpha, const void *T, int ldt,
		  void *x, int incx, enum blas_prec_type prec);
void BLAS_ztrsv_x(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag,
		  int n, const void *alpha, const void *T, int ldt,
		  void *x, int incx, enum blas_prec_type prec);
void BLAS_ztrsv_c_x(enum blas_order_type order, enum blas_uplo_type uplo,
		    enum blas_trans_type trans, enum blas_diag_type diag,
		    int n, const void *alpha, const void *T, int ldt,
		    void *x, int incx, enum blas_prec_type prec);
void BLAS_ctrsv_s_x(enum blas_order_type order, enum blas_uplo_type uplo,
		    enum blas_trans_type trans, enum blas_diag_type diag,
		    int n, const void *alpha, const float *T, int ldt,
		    void *x, int incx, enum blas_prec_type prec);
void BLAS_ztrsv_d_x(enum blas_order_type order, enum blas_uplo_type uplo,
		    enum blas_trans_type trans, enum blas_diag_type diag,
		    int n, const void *alpha, const double *T, int ldt,
		    void *x, int incx, enum blas_prec_type prec);


void BLAS_dtbsv_s(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag,
		  int n, int k, double alpha, const float *t, int ldt,
		  double *x, int incx);
void BLAS_ztbsv_c(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag,
		  int n, int k, const void *alpha, const void *t, int ldt,
		  void *x, int incx);
void BLAS_ctbsv_s(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag,
		  int n, int k, const void *alpha, const float *t, int ldt,
		  void *x, int incx);
void BLAS_ztbsv_d(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag,
		  int n, int k, const void *alpha, const double *t, int ldt,
		  void *x, int incx);
void BLAS_stbsv_x(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag,
		  int n, int k, float alpha, const float *t, int ldt,
		  float *x, int incx, enum blas_prec_type prec);
void BLAS_dtbsv_x(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag,
		  int n, int k, double alpha, const double *t, int ldt,
		  double *x, int incx, enum blas_prec_type prec);
void BLAS_dtbsv_s_x(enum blas_order_type order, enum blas_uplo_type uplo,
		    enum blas_trans_type trans, enum blas_diag_type diag,
		    int n, int k, double alpha, const float *t, int ldt,
		    double *x, int incx, enum blas_prec_type prec);
void BLAS_ctbsv_x(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag,
		  int n, int k, const void *alpha, const void *t, int ldt,
		  void *x, int incx, enum blas_prec_type prec);
void BLAS_ztbsv_x(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag,
		  int n, int k, const void *alpha, const void *t, int ldt,
		  void *x, int incx, enum blas_prec_type prec);
void BLAS_ztbsv_c_x(enum blas_order_type order, enum blas_uplo_type uplo,
		    enum blas_trans_type trans, enum blas_diag_type diag,
		    int n, int k, const void *alpha, const void *t, int ldt,
		    void *x, int incx, enum blas_prec_type prec);
void BLAS_ctbsv_s_x(enum blas_order_type order, enum blas_uplo_type uplo,
		    enum blas_trans_type trans, enum blas_diag_type diag,
		    int n, int k, const void *alpha, const float *t, int ldt,
		    void *x, int incx, enum blas_prec_type prec);
void BLAS_ztbsv_d_x(enum blas_order_type order, enum blas_uplo_type uplo,
		    enum blas_trans_type trans, enum blas_diag_type diag,
		    int n, int k, const void *alpha, const double *t, int ldt,
		    void *x, int incx, enum blas_prec_type prec);


void BLAS_dgemm_d_s(enum blas_order_type order, enum blas_trans_type transa,
		    enum blas_trans_type transb, int m, int n, int k,
		    double alpha, const double *a, int lda, const float *b,
		    int ldb, double beta, double *c, int ldc);
void BLAS_dgemm_s_d(enum blas_order_type order, enum blas_trans_type transa,
		    enum blas_trans_type transb, int m, int n, int k,
		    double alpha, const float *a, int lda, const double *b,
		    int ldb, double beta, double *c, int ldc);
void BLAS_dgemm_s_s(enum blas_order_type order, enum blas_trans_type transa,
		    enum blas_trans_type transb, int m, int n, int k,
		    double alpha, const float *a, int lda, const float *b,
		    int ldb, double beta, double *c, int ldc);
void BLAS_zgemm_z_c(enum blas_order_type order, enum blas_trans_type transa,
		    enum blas_trans_type transb, int m, int n, int k,
		    const void *alpha, const void *a, int lda, const void *b,
		    int ldb, const void *beta, void *c, int ldc);
void BLAS_zgemm_c_z(enum blas_order_type order, enum blas_trans_type transa,
		    enum blas_trans_type transb, int m, int n, int k,
		    const void *alpha, const void *a, int lda, const void *b,
		    int ldb, const void *beta, void *c, int ldc);
void BLAS_zgemm_c_c(enum blas_order_type order, enum blas_trans_type transa,
		    enum blas_trans_type transb, int m, int n, int k,
		    const void *alpha, const void *a, int lda, const void *b,
		    int ldb, const void *beta, void *c, int ldc);
void BLAS_cgemm_c_s(enum blas_order_type order, enum blas_trans_type transa,
		    enum blas_trans_type transb, int m, int n, int k,
		    const void *alpha, const void *a, int lda, const float *b,
		    int ldb, const void *beta, void *c, int ldc);
void BLAS_cgemm_s_c(enum blas_order_type order, enum blas_trans_type transa,
		    enum blas_trans_type transb, int m, int n, int k,
		    const void *alpha, const float *a, int lda, const void *b,
		    int ldb, const void *beta, void *c, int ldc);
void BLAS_cgemm_s_s(enum blas_order_type order, enum blas_trans_type transa,
		    enum blas_trans_type transb, int m, int n, int k,
		    const void *alpha, const float *a, int lda,
		    const float *b, int ldb, const void *beta, void *c,
		    int ldc);
void BLAS_zgemm_z_d(enum blas_order_type order, enum blas_trans_type transa,
		    enum blas_trans_type transb, int m, int n, int k,
		    const void *alpha, const void *a, int lda,
		    const double *b, int ldb, const void *beta, void *c,
		    int ldc);
void BLAS_zgemm_d_z(enum blas_order_type order, enum blas_trans_type transa,
		    enum blas_trans_type transb, int m, int n, int k,
		    const void *alpha, const double *a, int lda,
		    const void *b, int ldb, const void *beta, void *c,
		    int ldc);
void BLAS_zgemm_d_d(enum blas_order_type order, enum blas_trans_type transa,
		    enum blas_trans_type transb, int m, int n, int k,
		    const void *alpha, const double *a, int lda,
		    const double *b, int ldb, const void *beta, void *c,
		    int ldc);
void BLAS_sgemm_x(enum blas_order_type order, enum blas_trans_type transa,
		  enum blas_trans_type transb, int m, int n, int k,
		  float alpha, const float *a, int lda, const float *b,
		  int ldb, float beta, float *c, int ldc,
		  enum blas_prec_type prec);
void BLAS_dgemm_x(enum blas_order_type order, enum blas_trans_type transa,
		  enum blas_trans_type transb, int m, int n, int k,
		  double alpha, const double *a, int lda, const double *b,
		  int ldb, double beta, double *c, int ldc,
		  enum blas_prec_type prec);
void BLAS_cgemm_x(enum blas_order_type order, enum blas_trans_type transa,
		  enum blas_trans_type transb, int m, int n, int k,
		  const void *alpha, const void *a, int lda, const void *b,
		  int ldb, const void *beta, void *c, int ldc,
		  enum blas_prec_type prec);
void BLAS_zgemm_x(enum blas_order_type order, enum blas_trans_type transa,
		  enum blas_trans_type transb, int m, int n, int k,
		  const void *alpha, const void *a, int lda, const void *b,
		  int ldb, const void *beta, void *c, int ldc,
		  enum blas_prec_type prec);
void BLAS_dgemm_d_s_x(enum blas_order_type order, enum blas_trans_type transa,
		      enum blas_trans_type transb, int m, int n, int k,
		      double alpha, const double *a, int lda, const float *b,
		      int ldb, double beta, double *c, int ldc,
		      enum blas_prec_type prec);
void BLAS_dgemm_s_d_x(enum blas_order_type order, enum blas_trans_type transa,
		      enum blas_trans_type transb, int m, int n, int k,
		      double alpha, const float *a, int lda, const double *b,
		      int ldb, double beta, double *c, int ldc,
		      enum blas_prec_type prec);
void BLAS_dgemm_s_s_x(enum blas_order_type order, enum blas_trans_type transa,
		      enum blas_trans_type transb, int m, int n, int k,
		      double alpha, const float *a, int lda, const float *b,
		      int ldb, double beta, double *c, int ldc,
		      enum blas_prec_type prec);
void BLAS_zgemm_z_c_x(enum blas_order_type order, enum blas_trans_type transa,
		      enum blas_trans_type transb, int m, int n, int k,
		      const void *alpha, const void *a, int lda,
		      const void *b, int ldb, const void *beta, void *c,
		      int ldc, enum blas_prec_type prec);
void BLAS_zgemm_c_z_x(enum blas_order_type order, enum blas_trans_type transa,
		      enum blas_trans_type transb, int m, int n, int k,
		      const void *alpha, const void *a, int lda,
		      const void *b, int ldb, const void *beta, void *c,
		      int ldc, enum blas_prec_type prec);
void BLAS_zgemm_c_c_x(enum blas_order_type order, enum blas_trans_type transa,
		      enum blas_trans_type transb, int m, int n, int k,
		      const void *alpha, const void *a, int lda,
		      const void *b, int ldb, const void *beta, void *c,
		      int ldc, enum blas_prec_type prec);
void BLAS_cgemm_c_s_x(enum blas_order_type order, enum blas_trans_type transa,
		      enum blas_trans_type transb, int m, int n, int k,
		      const void *alpha, const void *a, int lda,
		      const float *b, int ldb, const void *beta, void *c,
		      int ldc, enum blas_prec_type prec);
void BLAS_cgemm_s_c_x(enum blas_order_type order, enum blas_trans_type transa,
		      enum blas_trans_type transb, int m, int n, int k,
		      const void *alpha, const float *a, int lda,
		      const void *b, int ldb, const void *beta, void *c,
		      int ldc, enum blas_prec_type prec);
void BLAS_cgemm_s_s_x(enum blas_order_type order, enum blas_trans_type transa,
		      enum blas_trans_type transb, int m, int n, int k,
		      const void *alpha, const float *a, int lda,
		      const float *b, int ldb, const void *beta, void *c,
		      int ldc, enum blas_prec_type prec);
void BLAS_zgemm_z_d_x(enum blas_order_type order, enum blas_trans_type transa,
		      enum blas_trans_type transb, int m, int n, int k,
		      const void *alpha, const void *a, int lda,
		      const double *b, int ldb, const void *beta, void *c,
		      int ldc, enum blas_prec_type prec);
void BLAS_zgemm_d_z_x(enum blas_order_type order, enum blas_trans_type transa,
		      enum blas_trans_type transb, int m, int n, int k,
		      const void *alpha, const double *a, int lda,
		      const void *b, int ldb, const void *beta, void *c,
		      int ldc, enum blas_prec_type prec);
void BLAS_zgemm_d_d_x(enum blas_order_type order, enum blas_trans_type transa,
		      enum blas_trans_type transb, int m, int n, int k,
		      const void *alpha, const double *a, int lda,
		      const double *b, int ldb, const void *beta, void *c,
		      int ldc, enum blas_prec_type prec);


void BLAS_dsymm_d_s(enum blas_order_type order, enum blas_side_type side,
		    enum blas_uplo_type uplo, int m, int n,
		    double alpha, const double *a, int lda,
		    const float *b, int ldb, double beta, double *c, int ldc);
void BLAS_dsymm_s_d(enum blas_order_type order, enum blas_side_type side,
		    enum blas_uplo_type uplo, int m, int n,
		    double alpha, const float *a, int lda,
		    const double *b, int ldb, double beta,
		    double *c, int ldc);
void BLAS_dsymm_s_s(enum blas_order_type order, enum blas_side_type side,
		    enum blas_uplo_type uplo, int m, int n,
		    double alpha, const float *a, int lda,
		    const float *b, int ldb, double beta, double *c, int ldc);
void BLAS_zsymm_z_c(enum blas_order_type order, enum blas_side_type side,
		    enum blas_uplo_type uplo, int m, int n,
		    const void *alpha, const void *a, int lda,
		    const void *b, int ldb, const void *beta,
		    void *c, int ldc);
void BLAS_zsymm_c_z(enum blas_order_type order, enum blas_side_type side,
		    enum blas_uplo_type uplo, int m, int n,
		    const void *alpha, const void *a, int lda,
		    const void *b, int ldb, const void *beta,
		    void *c, int ldc);
void BLAS_zsymm_c_c(enum blas_order_type order, enum blas_side_type side,
		    enum blas_uplo_type uplo, int m, int n,
		    const void *alpha, const void *a, int lda,
		    const void *b, int ldb, const void *beta,
		    void *c, int ldc);
void BLAS_csymm_c_s(enum blas_order_type order, enum blas_side_type side,
		    enum blas_uplo_type uplo, int m, int n,
		    const void *alpha, const void *a, int lda,
		    const float *b, int ldb, const void *beta,
		    void *c, int ldc);
void BLAS_csymm_s_c(enum blas_order_type order, enum blas_side_type side,
		    enum blas_uplo_type uplo, int m, int n,
		    const void *alpha, const float *a, int lda,
		    const void *b, int ldb, const void *beta,
		    void *c, int ldc);
void BLAS_csymm_s_s(enum blas_order_type order, enum blas_side_type side,
		    enum blas_uplo_type uplo, int m, int n,
		    const void *alpha, const float *a, int lda,
		    const float *b, int ldb, const void *beta,
		    void *c, int ldc);
void BLAS_zsymm_z_d(enum blas_order_type order, enum blas_side_type side,
		    enum blas_uplo_type uplo, int m, int n,
		    const void *alpha, const void *a, int lda,
		    const double *b, int ldb, const void *beta,
		    void *c, int ldc);
void BLAS_zsymm_d_z(enum blas_order_type order, enum blas_side_type side,
		    enum blas_uplo_type uplo, int m, int n,
		    const void *alpha, const double *a, int lda,
		    const void *b, int ldb, const void *beta,
		    void *c, int ldc);
void BLAS_zsymm_d_d(enum blas_order_type order, enum blas_side_type side,
		    enum blas_uplo_type uplo, int m, int n,
		    const void *alpha, const double *a, int lda,
		    const double *b, int ldb, const void *beta,
		    void *c, int ldc);
void BLAS_ssymm_x(enum blas_order_type order, enum blas_side_type side,
		  enum blas_uplo_type uplo, int m, int n,
		  float alpha, const float *a, int lda,
		  const float *b, int ldb, float beta,
		  float *c, int ldc, enum blas_prec_type prec);
void BLAS_dsymm_x(enum blas_order_type order, enum blas_side_type side,
		  enum blas_uplo_type uplo, int m, int n,
		  double alpha, const double *a, int lda,
		  const double *b, int ldb, double beta,
		  double *c, int ldc, enum blas_prec_type prec);
void BLAS_csymm_x(enum blas_order_type order, enum blas_side_type side,
		  enum blas_uplo_type uplo, int m, int n,
		  const void *alpha, const void *a, int lda,
		  const void *b, int ldb, const void *beta,
		  void *c, int ldc, enum blas_prec_type prec);
void BLAS_zsymm_x(enum blas_order_type order, enum blas_side_type side,
		  enum blas_uplo_type uplo, int m, int n,
		  const void *alpha, const void *a, int lda,
		  const void *b, int ldb, const void *beta,
		  void *c, int ldc, enum blas_prec_type prec);
void BLAS_dsymm_d_s_x(enum blas_order_type order, enum blas_side_type side,
		      enum blas_uplo_type uplo, int m, int n,
		      double alpha, const double *a, int lda,
		      const float *b, int ldb, double beta,
		      double *c, int ldc, enum blas_prec_type prec);
void BLAS_dsymm_s_d_x(enum blas_order_type order, enum blas_side_type side,
		      enum blas_uplo_type uplo, int m, int n,
		      double alpha, const float *a, int lda,
		      const double *b, int ldb, double beta,
		      double *c, int ldc, enum blas_prec_type prec);
void BLAS_dsymm_s_s_x(enum blas_order_type order, enum blas_side_type side,
		      enum blas_uplo_type uplo, int m, int n,
		      double alpha, const float *a, int lda,
		      const float *b, int ldb, double beta,
		      double *c, int ldc, enum blas_prec_type prec);
void BLAS_zsymm_z_c_x(enum blas_order_type order, enum blas_side_type side,
		      enum blas_uplo_type uplo, int m, int n,
		      const void *alpha, const void *a, int lda,
		      const void *b, int ldb, const void *beta,
		      void *c, int ldc, enum blas_prec_type prec);
void BLAS_zsymm_c_z_x(enum blas_order_type order, enum blas_side_type side,
		      enum blas_uplo_type uplo, int m, int n,
		      const void *alpha, const void *a, int lda,
		      const void *b, int ldb, const void *beta,
		      void *c, int ldc, enum blas_prec_type prec);
void BLAS_zsymm_c_c_x(enum blas_order_type order, enum blas_side_type side,
		      enum blas_uplo_type uplo, int m, int n,
		      const void *alpha, const void *a, int lda,
		      const void *b, int ldb, const void *beta,
		      void *c, int ldc, enum blas_prec_type prec);
void BLAS_csymm_c_s_x(enum blas_order_type order, enum blas_side_type side,
		      enum blas_uplo_type uplo, int m, int n,
		      const void *alpha, const void *a, int lda,
		      const float *b, int ldb, const void *beta,
		      void *c, int ldc, enum blas_prec_type prec);
void BLAS_csymm_s_c_x(enum blas_order_type order, enum blas_side_type side,
		      enum blas_uplo_type uplo, int m, int n,
		      const void *alpha, const float *a, int lda,
		      const void *b, int ldb, const void *beta,
		      void *c, int ldc, enum blas_prec_type prec);
void BLAS_csymm_s_s_x(enum blas_order_type order, enum blas_side_type side,
		      enum blas_uplo_type uplo, int m, int n,
		      const void *alpha, const float *a, int lda,
		      const float *b, int ldb, const void *beta,
		      void *c, int ldc, enum blas_prec_type prec);
void BLAS_zsymm_z_d_x(enum blas_order_type order, enum blas_side_type side,
		      enum blas_uplo_type uplo, int m, int n,
		      const void *alpha, const void *a, int lda,
		      const double *b, int ldb, const void *beta,
		      void *c, int ldc, enum blas_prec_type prec);
void BLAS_zsymm_d_z_x(enum blas_order_type order, enum blas_side_type side,
		      enum blas_uplo_type uplo, int m, int n,
		      const void *alpha, const double *a, int lda,
		      const void *b, int ldb, const void *beta,
		      void *c, int ldc, enum blas_prec_type prec);
void BLAS_zsymm_d_d_x(enum blas_order_type order, enum blas_side_type side,
		      enum blas_uplo_type uplo, int m, int n,
		      const void *alpha, const double *a, int lda,
		      const double *b, int ldb, const void *beta,
		      void *c, int ldc, enum blas_prec_type prec);


void BLAS_zhemm_z_c(enum blas_order_type order, enum blas_side_type side,
		    enum blas_uplo_type uplo, int m, int n,
		    const void *alpha, const void *a, int lda,
		    const void *b, int ldb, const void *beta,
		    void *c, int ldc);
void BLAS_zhemm_c_z(enum blas_order_type order, enum blas_side_type side,
		    enum blas_uplo_type uplo, int m, int n,
		    const void *alpha, const void *a, int lda,
		    const void *b, int ldb, const void *beta,
		    void *c, int ldc);
void BLAS_zhemm_c_c(enum blas_order_type order, enum blas_side_type side,
		    enum blas_uplo_type uplo, int m, int n,
		    const void *alpha, const void *a, int lda,
		    const void *b, int ldb, const void *beta,
		    void *c, int ldc);
void BLAS_chemm_c_s(enum blas_order_type order, enum blas_side_type side,
		    enum blas_uplo_type uplo, int m, int n,
		    const void *alpha, const void *a, int lda,
		    const float *b, int ldb, const void *beta,
		    void *c, int ldc);
void BLAS_zhemm_z_d(enum blas_order_type order, enum blas_side_type side,
		    enum blas_uplo_type uplo, int m, int n,
		    const void *alpha, const void *a, int lda,
		    const double *b, int ldb, const void *beta,
		    void *c, int ldc);
void BLAS_chemm_x(enum blas_order_type order, enum blas_side_type side,
		  enum blas_uplo_type uplo, int m, int n,
		  const void *alpha, const void *a, int lda,
		  const void *b, int ldb, const void *beta,
		  void *c, int ldc, enum blas_prec_type prec);
void BLAS_zhemm_x(enum blas_order_type order, enum blas_side_type side,
		  enum blas_uplo_type uplo, int m, int n,
		  const void *alpha, const void *a, int lda,
		  const void *b, int ldb, const void *beta,
		  void *c, int ldc, enum blas_prec_type prec);
void BLAS_zhemm_z_c_x(enum blas_order_type order, enum blas_side_type side,
		      enum blas_uplo_type uplo, int m, int n,
		      const void *alpha, const void *a, int lda,
		      const void *b, int ldb, const void *beta,
		      void *c, int ldc, enum blas_prec_type prec);
void BLAS_zhemm_c_z_x(enum blas_order_type order, enum blas_side_type side,
		      enum blas_uplo_type uplo, int m, int n,
		      const void *alpha, const void *a, int lda,
		      const void *b, int ldb, const void *beta,
		      void *c, int ldc, enum blas_prec_type prec);
void BLAS_zhemm_c_c_x(enum blas_order_type order, enum blas_side_type side,
		      enum blas_uplo_type uplo, int m, int n,
		      const void *alpha, const void *a, int lda,
		      const void *b, int ldb, const void *beta,
		      void *c, int ldc, enum blas_prec_type prec);
void BLAS_chemm_c_s_x(enum blas_order_type order, enum blas_side_type side,
		      enum blas_uplo_type uplo, int m, int n,
		      const void *alpha, const void *a, int lda,
		      const float *b, int ldb, const void *beta,
		      void *c, int ldc, enum blas_prec_type prec);
void BLAS_zhemm_z_d_x(enum blas_order_type order, enum blas_side_type side,
		      enum blas_uplo_type uplo, int m, int n,
		      const void *alpha, const void *a, int lda,
		      const double *b, int ldb, const void *beta,
		      void *c, int ldc, enum blas_prec_type prec);


void BLAS_dgemv2_d_s(enum blas_order_type order, enum blas_trans_type trans,
		     int m, int n, double alpha, const double *a, int lda,
		     const float *head_x, const float *tail_x, int incx,
		     double beta, double *y, int incy);
void BLAS_dgemv2_s_d(enum blas_order_type order, enum blas_trans_type trans,
		     int m, int n, double alpha, const float *a, int lda,
		     const double *head_x, const double *tail_x, int incx,
		     double beta, double *y, int incy);
void BLAS_dgemv2_s_s(enum blas_order_type order, enum blas_trans_type trans,
		     int m, int n, double alpha, const float *a, int lda,
		     const float *head_x, const float *tail_x, int incx,
		     double beta, double *y, int incy);
void BLAS_zgemv2_z_c(enum blas_order_type order, enum blas_trans_type trans,
		     int m, int n, const void *alpha, const void *a, int lda,
		     const void *head_x, const void *tail_x, int incx,
		     const void *beta, void *y, int incy);
void BLAS_zgemv2_c_z(enum blas_order_type order, enum blas_trans_type trans,
		     int m, int n, const void *alpha, const void *a, int lda,
		     const void *head_x, const void *tail_x, int incx,
		     const void *beta, void *y, int incy);
void BLAS_zgemv2_c_c(enum blas_order_type order, enum blas_trans_type trans,
		     int m, int n, const void *alpha, const void *a, int lda,
		     const void *head_x, const void *tail_x, int incx,
		     const void *beta, void *y, int incy);
void BLAS_cgemv2_c_s(enum blas_order_type order, enum blas_trans_type trans,
		     int m, int n, const void *alpha, const void *a, int lda,
		     const float *head_x, const float *tail_x, int incx,
		     const void *beta, void *y, int incy);
void BLAS_cgemv2_s_c(enum blas_order_type order, enum blas_trans_type trans,
		     int m, int n, const void *alpha, const float *a, int lda,
		     const void *head_x, const void *tail_x, int incx,
		     const void *beta, void *y, int incy);
void BLAS_cgemv2_s_s(enum blas_order_type order, enum blas_trans_type trans,
		     int m, int n, const void *alpha, const float *a, int lda,
		     const float *head_x, const float *tail_x, int incx,
		     const void *beta, void *y, int incy);
void BLAS_zgemv2_z_d(enum blas_order_type order, enum blas_trans_type trans,
		     int m, int n, const void *alpha, const void *a, int lda,
		     const double *head_x, const double *tail_x, int incx,
		     const void *beta, void *y, int incy);
void BLAS_zgemv2_d_z(enum blas_order_type order, enum blas_trans_type trans,
		     int m, int n, const void *alpha, const double *a,
		     int lda, const void *head_x, const void *tail_x,
		     int incx, const void *beta, void *y, int incy);
void BLAS_zgemv2_d_d(enum blas_order_type order, enum blas_trans_type trans,
		     int m, int n, const void *alpha, const double *a,
		     int lda, const double *head_x, const double *tail_x,
		     int incx, const void *beta, void *y, int incy);
void BLAS_sgemv2_x(enum blas_order_type order, enum blas_trans_type trans,
		   int m, int n, float alpha, const float *a, int lda,
		   const float *head_x, const float *tail_x, int incx,
		   float beta, float *y, int incy, enum blas_prec_type prec);
void BLAS_dgemv2_x(enum blas_order_type order, enum blas_trans_type trans,
		   int m, int n, double alpha, const double *a, int lda,
		   const double *head_x, const double *tail_x, int incx,
		   double beta, double *y, int incy,
		   enum blas_prec_type prec);
void BLAS_cgemv2_x(enum blas_order_type order, enum blas_trans_type trans,
		   int m, int n, const void *alpha, const void *a, int lda,
		   const void *head_x, const void *tail_x, int incx,
		   const void *beta, void *y, int incy,
		   enum blas_prec_type prec);
void BLAS_zgemv2_x(enum blas_order_type order, enum blas_trans_type trans,
		   int m, int n, const void *alpha, const void *a, int lda,
		   const void *head_x, const void *tail_x, int incx,
		   const void *beta, void *y, int incy,
		   enum blas_prec_type prec);
void BLAS_dgemv2_d_s_x(enum blas_order_type order, enum blas_trans_type trans,
		       int m, int n, double alpha, const double *a, int lda,
		       const float *head_x, const float *tail_x, int incx,
		       double beta, double *y, int incy,
		       enum blas_prec_type prec);
void BLAS_dgemv2_s_d_x(enum blas_order_type order, enum blas_trans_type trans,
		       int m, int n, double alpha, const float *a, int lda,
		       const double *head_x, const double *tail_x, int incx,
		       double beta, double *y, int incy,
		       enum blas_prec_type prec);
void BLAS_dgemv2_s_s_x(enum blas_order_type order, enum blas_trans_type trans,
		       int m, int n, double alpha, const float *a, int lda,
		       const float *head_x, const float *tail_x, int incx,
		       double beta, double *y, int incy,
		       enum blas_prec_type prec);
void BLAS_zgemv2_z_c_x(enum blas_order_type order, enum blas_trans_type trans,
		       int m, int n, const void *alpha, const void *a,
		       int lda, const void *head_x, const void *tail_x,
		       int incx, const void *beta, void *y, int incy,
		       enum blas_prec_type prec);
void BLAS_zgemv2_c_z_x(enum blas_order_type order, enum blas_trans_type trans,
		       int m, int n, const void *alpha, const void *a,
		       int lda, const void *head_x, const void *tail_x,
		       int incx, const void *beta, void *y, int incy,
		       enum blas_prec_type prec);
void BLAS_zgemv2_c_c_x(enum blas_order_type order, enum blas_trans_type trans,
		       int m, int n, const void *alpha, const void *a,
		       int lda, const void *head_x, const void *tail_x,
		       int incx, const void *beta, void *y, int incy,
		       enum blas_prec_type prec);
void BLAS_cgemv2_c_s_x(enum blas_order_type order, enum blas_trans_type trans,
		       int m, int n, const void *alpha, const void *a,
		       int lda, const float *head_x, const float *tail_x,
		       int incx, const void *beta, void *y, int incy,
		       enum blas_prec_type prec);
void BLAS_cgemv2_s_c_x(enum blas_order_type order, enum blas_trans_type trans,
		       int m, int n, const void *alpha, const float *a,
		       int lda, const void *head_x, const void *tail_x,
		       int incx, const void *beta, void *y, int incy,
		       enum blas_prec_type prec);
void BLAS_cgemv2_s_s_x(enum blas_order_type order, enum blas_trans_type trans,
		       int m, int n, const void *alpha, const float *a,
		       int lda, const float *head_x, const float *tail_x,
		       int incx, const void *beta, void *y, int incy,
		       enum blas_prec_type prec);
void BLAS_zgemv2_z_d_x(enum blas_order_type order, enum blas_trans_type trans,
		       int m, int n, const void *alpha, const void *a,
		       int lda, const double *head_x, const double *tail_x,
		       int incx, const void *beta, void *y, int incy,
		       enum blas_prec_type prec);
void BLAS_zgemv2_d_z_x(enum blas_order_type order, enum blas_trans_type trans,
		       int m, int n, const void *alpha, const double *a,
		       int lda, const void *head_x, const void *tail_x,
		       int incx, const void *beta, void *y, int incy,
		       enum blas_prec_type prec);
void BLAS_zgemv2_d_d_x(enum blas_order_type order, enum blas_trans_type trans,
		       int m, int n, const void *alpha, const double *a,
		       int lda, const double *head_x, const double *tail_x,
		       int incx, const void *beta, void *y, int incy,
		       enum blas_prec_type prec);


void BLAS_dsymv2_d_s(enum blas_order_type order, enum blas_uplo_type uplo,
		     int n, double alpha, const double *a, int lda,
		     const float *x_head, const float *x_tail, int incx,
		     double beta, double *y, int incy);
void BLAS_dsymv2_s_d(enum blas_order_type order, enum blas_uplo_type uplo,
		     int n, double alpha, const float *a, int lda,
		     const double *x_head, const double *x_tail, int incx,
		     double beta, double *y, int incy);
void BLAS_dsymv2_s_s(enum blas_order_type order, enum blas_uplo_type uplo,
		     int n, double alpha, const float *a, int lda,
		     const float *x_head, const float *x_tail, int incx,
		     double beta, double *y, int incy);
void BLAS_zsymv2_z_c(enum blas_order_type order, enum blas_uplo_type uplo,
		     int n, const void *alpha, const void *a, int lda,
		     const void *x_head, const void *x_tail, int incx,
		     const void *beta, void *y, int incy);
void BLAS_zsymv2_c_z(enum blas_order_type order, enum blas_uplo_type uplo,
		     int n, const void *alpha, const void *a, int lda,
		     const void *x_head, const void *x_tail, int incx,
		     const void *beta, void *y, int incy);
void BLAS_zsymv2_c_c(enum blas_order_type order, enum blas_uplo_type uplo,
		     int n, const void *alpha, const void *a, int lda,
		     const void *x_head, const void *x_tail, int incx,
		     const void *beta, void *y, int incy);
void BLAS_csymv2_c_s(enum blas_order_type order, enum blas_uplo_type uplo,
		     int n, const void *alpha, const void *a, int lda,
		     const float *x_head, const float *x_tail, int incx,
		     const void *beta, void *y, int incy);
void BLAS_csymv2_s_c(enum blas_order_type order, enum blas_uplo_type uplo,
		     int n, const void *alpha, const float *a, int lda,
		     const void *x_head, const void *x_tail, int incx,
		     const void *beta, void *y, int incy);
void BLAS_csymv2_s_s(enum blas_order_type order, enum blas_uplo_type uplo,
		     int n, const void *alpha, const float *a, int lda,
		     const float *x_head, const float *x_tail, int incx,
		     const void *beta, void *y, int incy);
void BLAS_zsymv2_z_d(enum blas_order_type order, enum blas_uplo_type uplo,
		     int n, const void *alpha, const void *a, int lda,
		     const double *x_head, const double *x_tail, int incx,
		     const void *beta, void *y, int incy);
void BLAS_zsymv2_d_z(enum blas_order_type order, enum blas_uplo_type uplo,
		     int n, const void *alpha, const double *a, int lda,
		     const void *x_head, const void *x_tail, int incx,
		     const void *beta, void *y, int incy);
void BLAS_zsymv2_d_d(enum blas_order_type order, enum blas_uplo_type uplo,
		     int n, const void *alpha, const double *a, int lda,
		     const double *x_head, const double *x_tail, int incx,
		     const void *beta, void *y, int incy);
void BLAS_ssymv2_x(enum blas_order_type order, enum blas_uplo_type uplo,
		   int n, float alpha, const float *a, int lda,
		   const float *x_head, const float *x_tail, int incx,
		   float beta, float *y, int incy, enum blas_prec_type prec);
void BLAS_dsymv2_x(enum blas_order_type order, enum blas_uplo_type uplo,
		   int n, double alpha, const double *a, int lda,
		   const double *x_head, const double *x_tail, int incx,
		   double beta, double *y, int incy,
		   enum blas_prec_type prec);
void BLAS_csymv2_x(enum blas_order_type order, enum blas_uplo_type uplo,
		   int n, const void *alpha, const void *a, int lda,
		   const void *x_head, const void *x_tail, int incx,
		   const void *beta, void *y, int incy,
		   enum blas_prec_type prec);
void BLAS_zsymv2_x(enum blas_order_type order, enum blas_uplo_type uplo,
		   int n, const void *alpha, const void *a, int lda,
		   const void *x_head, const void *x_tail, int incx,
		   const void *beta, void *y, int incy,
		   enum blas_prec_type prec);
void BLAS_dsymv2_d_s_x(enum blas_order_type order, enum blas_uplo_type uplo,
		       int n, double alpha, const double *a, int lda,
		       const float *x_head, const float *x_tail, int incx,
		       double beta, double *y, int incy,
		       enum blas_prec_type prec);
void BLAS_dsymv2_s_d_x(enum blas_order_type order, enum blas_uplo_type uplo,
		       int n, double alpha, const float *a, int lda,
		       const double *x_head, const double *x_tail, int incx,
		       double beta, double *y, int incy,
		       enum blas_prec_type prec);
void BLAS_dsymv2_s_s_x(enum blas_order_type order, enum blas_uplo_type uplo,
		       int n, double alpha, const float *a, int lda,
		       const float *x_head, const float *x_tail, int incx,
		       double beta, double *y, int incy,
		       enum blas_prec_type prec);
void BLAS_zsymv2_z_c_x(enum blas_order_type order, enum blas_uplo_type uplo,
		       int n, const void *alpha, const void *a, int lda,
		       const void *x_head, const void *x_tail, int incx,
		       const void *beta, void *y, int incy,
		       enum blas_prec_type prec);
void BLAS_zsymv2_c_z_x(enum blas_order_type order, enum blas_uplo_type uplo,
		       int n, const void *alpha, const void *a, int lda,
		       const void *x_head, const void *x_tail, int incx,
		       const void *beta, void *y, int incy,
		       enum blas_prec_type prec);
void BLAS_zsymv2_c_c_x(enum blas_order_type order, enum blas_uplo_type uplo,
		       int n, const void *alpha, const void *a, int lda,
		       const void *x_head, const void *x_tail, int incx,
		       const void *beta, void *y, int incy,
		       enum blas_prec_type prec);
void BLAS_csymv2_c_s_x(enum blas_order_type order, enum blas_uplo_type uplo,
		       int n, const void *alpha, const void *a, int lda,
		       const float *x_head, const float *x_tail, int incx,
		       const void *beta, void *y, int incy,
		       enum blas_prec_type prec);
void BLAS_csymv2_s_c_x(enum blas_order_type order, enum blas_uplo_type uplo,
		       int n, const void *alpha, const float *a, int lda,
		       const void *x_head, const void *x_tail, int incx,
		       const void *beta, void *y, int incy,
		       enum blas_prec_type prec);
void BLAS_csymv2_s_s_x(enum blas_order_type order, enum blas_uplo_type uplo,
		       int n, const void *alpha, const float *a, int lda,
		       const float *x_head, const float *x_tail, int incx,
		       const void *beta, void *y, int incy,
		       enum blas_prec_type prec);
void BLAS_zsymv2_z_d_x(enum blas_order_type order, enum blas_uplo_type uplo,
		       int n, const void *alpha, const void *a, int lda,
		       const double *x_head, const double *x_tail, int incx,
		       const void *beta, void *y, int incy,
		       enum blas_prec_type prec);
void BLAS_zsymv2_d_z_x(enum blas_order_type order, enum blas_uplo_type uplo,
		       int n, const void *alpha, const double *a, int lda,
		       const void *x_head, const void *x_tail, int incx,
		       const void *beta, void *y, int incy,
		       enum blas_prec_type prec);
void BLAS_zsymv2_d_d_x(enum blas_order_type order, enum blas_uplo_type uplo,
		       int n, const void *alpha, const double *a, int lda,
		       const double *x_head, const double *x_tail, int incx,
		       const void *beta, void *y, int incy,
		       enum blas_prec_type prec);


void BLAS_zhemv2_z_c(enum blas_order_type order, enum blas_uplo_type uplo,
		     int n, const void *alpha, const void *a, int lda,
		     const void *x_head, const void *x_tail, int incx,
		     const void *beta, const void *y, int incy);
void BLAS_zhemv2_c_z(enum blas_order_type order, enum blas_uplo_type uplo,
		     int n, const void *alpha, const void *a, int lda,
		     const void *x_head, const void *x_tail, int incx,
		     const void *beta, const void *y, int incy);
void BLAS_zhemv2_c_c(enum blas_order_type order, enum blas_uplo_type uplo,
		     int n, const void *alpha, const void *a, int lda,
		     const void *x_head, const void *x_tail, int incx,
		     const void *beta, const void *y, int incy);
void BLAS_chemv2_c_s(enum blas_order_type order, enum blas_uplo_type uplo,
		     int n, const void *alpha, const void *a, int lda,
		     const float *x_head, const float *x_tail, int incx,
		     const void *beta, const float *y, int incy);
void BLAS_zhemv2_z_d(enum blas_order_type order, enum blas_uplo_type uplo,
		     int n, const void *alpha, const void *a, int lda,
		     const double *x_head, const double *x_tail, int incx,
		     const void *beta, const double *y, int incy);
void BLAS_chemv2_x(enum blas_order_type order, enum blas_uplo_type uplo,
		   int n, const void *alpha, const void *a, int lda,
		   const void *x_head, const void *x_tail, int incx,
		   const void *beta, const void *y, int incy,
		   enum blas_prec_type prec);
void BLAS_zhemv2_x(enum blas_order_type order, enum blas_uplo_type uplo,
		   int n, const void *alpha, const void *a, int lda,
		   const void *x_head, const void *x_tail, int incx,
		   const void *beta, const void *y, int incy,
		   enum blas_prec_type prec);
void BLAS_zhemv2_z_c_x(enum blas_order_type order, enum blas_uplo_type uplo,
		       int n, const void *alpha, const void *a, int lda,
		       const void *x_head, const void *x_tail, int incx,
		       const void *beta, const void *y, int incy,
		       enum blas_prec_type prec);
void BLAS_zhemv2_c_z_x(enum blas_order_type order, enum blas_uplo_type uplo,
		       int n, const void *alpha, const void *a, int lda,
		       const void *x_head, const void *x_tail, int incx,
		       const void *beta, const void *y, int incy,
		       enum blas_prec_type prec);
void BLAS_zhemv2_c_c_x(enum blas_order_type order, enum blas_uplo_type uplo,
		       int n, const void *alpha, const void *a, int lda,
		       const void *x_head, const void *x_tail, int incx,
		       const void *beta, const void *y, int incy,
		       enum blas_prec_type prec);
void BLAS_chemv2_c_s_x(enum blas_order_type order, enum blas_uplo_type uplo,
		       int n, const void *alpha, const void *a, int lda,
		       const float *x_head, const float *x_tail, int incx,
		       const void *beta, const float *y, int incy,
		       enum blas_prec_type prec);
void BLAS_zhemv2_z_d_x(enum blas_order_type order, enum blas_uplo_type uplo,
		       int n, const void *alpha, const void *a, int lda,
		       const double *x_head, const double *x_tail, int incx,
		       const void *beta, const double *y, int incy,
		       enum blas_prec_type prec);


void BLAS_dgbmv2_d_s(enum blas_order_type order, enum blas_trans_type trans,
		     int m, int n, int kl, int ku, double alpha,
		     const double *a, int lda, const float *head_x,
		     const float *tail_x, int incx, double beta,
		     double *y, int incy);
void BLAS_dgbmv2_s_d(enum blas_order_type order, enum blas_trans_type trans,
		     int m, int n, int kl, int ku, double alpha,
		     const float *a, int lda, const double *head_x,
		     const double *tail_x, int incx, double beta,
		     double *y, int incy);
void BLAS_dgbmv2_s_s(enum blas_order_type order, enum blas_trans_type trans,
		     int m, int n, int kl, int ku, double alpha,
		     const float *a, int lda, const float *head_x,
		     const float *tail_x, int incx, double beta,
		     double *y, int incy);
void BLAS_zgbmv2_z_c(enum blas_order_type order, enum blas_trans_type trans,
		     int m, int n, int kl, int ku, const void *alpha,
		     const void *a, int lda, const void *head_x,
		     const void *tail_x, int incx, const void *beta,
		     void *y, int incy);
void BLAS_zgbmv2_c_z(enum blas_order_type order, enum blas_trans_type trans,
		     int m, int n, int kl, int ku, const void *alpha,
		     const void *a, int lda, const void *head_x,
		     const void *tail_x, int incx, const void *beta,
		     void *y, int incy);
void BLAS_zgbmv2_c_c(enum blas_order_type order, enum blas_trans_type trans,
		     int m, int n, int kl, int ku, const void *alpha,
		     const void *a, int lda, const void *head_x,
		     const void *tail_x, int incx, const void *beta,
		     void *y, int incy);
void BLAS_cgbmv2_c_s(enum blas_order_type order, enum blas_trans_type trans,
		     int m, int n, int kl, int ku, const void *alpha,
		     const void *a, int lda, const float *head_x,
		     const float *tail_x, int incx, const void *beta,
		     void *y, int incy);
void BLAS_cgbmv2_s_c(enum blas_order_type order, enum blas_trans_type trans,
		     int m, int n, int kl, int ku, const void *alpha,
		     const float *a, int lda, const void *head_x,
		     const void *tail_x, int incx, const void *beta,
		     void *y, int incy);
void BLAS_cgbmv2_s_s(enum blas_order_type order, enum blas_trans_type trans,
		     int m, int n, int kl, int ku, const void *alpha,
		     const float *a, int lda, const float *head_x,
		     const float *tail_x, int incx, const void *beta,
		     void *y, int incy);
void BLAS_zgbmv2_z_d(enum blas_order_type order, enum blas_trans_type trans,
		     int m, int n, int kl, int ku, const void *alpha,
		     const void *a, int lda, const double *head_x,
		     const double *tail_x, int incx, const void *beta,
		     void *y, int incy);
void BLAS_zgbmv2_d_z(enum blas_order_type order, enum blas_trans_type trans,
		     int m, int n, int kl, int ku, const void *alpha,
		     const double *a, int lda, const void *head_x,
		     const void *tail_x, int incx, const void *beta,
		     void *y, int incy);
void BLAS_zgbmv2_d_d(enum blas_order_type order, enum blas_trans_type trans,
		     int m, int n, int kl, int ku, const void *alpha,
		     const double *a, int lda, const double *head_x,
		     const double *tail_x, int incx, const void *beta,
		     void *y, int incy);
void BLAS_sgbmv2_x(enum blas_order_type order, enum blas_trans_type trans,
		   int m, int n, int kl, int ku, float alpha,
		   const float *a, int lda, const float *head_x,
		   const float *tail_x, int incx, float beta,
		   float *y, int incy, enum blas_prec_type prec);
void BLAS_dgbmv2_x(enum blas_order_type order, enum blas_trans_type trans,
		   int m, int n, int kl, int ku, double alpha,
		   const double *a, int lda, const double *head_x,
		   const double *tail_x, int incx, double beta,
		   double *y, int incy, enum blas_prec_type prec);
void BLAS_cgbmv2_x(enum blas_order_type order, enum blas_trans_type trans,
		   int m, int n, int kl, int ku, const void *alpha,
		   const void *a, int lda, const void *head_x,
		   const void *tail_x, int incx, const void *beta,
		   void *y, int incy, enum blas_prec_type prec);
void BLAS_zgbmv2_x(enum blas_order_type order, enum blas_trans_type trans,
		   int m, int n, int kl, int ku, const void *alpha,
		   const void *a, int lda, const void *head_x,
		   const void *tail_x, int incx, const void *beta,
		   void *y, int incy, enum blas_prec_type prec);
void BLAS_dgbmv2_d_s_x(enum blas_order_type order, enum blas_trans_type trans,
		       int m, int n, int kl, int ku, double alpha,
		       const double *a, int lda, const float *head_x,
		       const float *tail_x, int incx, double beta,
		       double *y, int incy, enum blas_prec_type prec);
void BLAS_dgbmv2_s_d_x(enum blas_order_type order, enum blas_trans_type trans,
		       int m, int n, int kl, int ku, double alpha,
		       const float *a, int lda, const double *head_x,
		       const double *tail_x, int incx, double beta,
		       double *y, int incy, enum blas_prec_type prec);
void BLAS_dgbmv2_s_s_x(enum blas_order_type order, enum blas_trans_type trans,
		       int m, int n, int kl, int ku, double alpha,
		       const float *a, int lda, const float *head_x,
		       const float *tail_x, int incx, double beta,
		       double *y, int incy, enum blas_prec_type prec);
void BLAS_zgbmv2_z_c_x(enum blas_order_type order, enum blas_trans_type trans,
		       int m, int n, int kl, int ku, const void *alpha,
		       const void *a, int lda, const void *head_x,
		       const void *tail_x, int incx, const void *beta,
		       void *y, int incy, enum blas_prec_type prec);
void BLAS_zgbmv2_c_z_x(enum blas_order_type order, enum blas_trans_type trans,
		       int m, int n, int kl, int ku, const void *alpha,
		       const void *a, int lda, const void *head_x,
		       const void *tail_x, int incx, const void *beta,
		       void *y, int incy, enum blas_prec_type prec);
void BLAS_zgbmv2_c_c_x(enum blas_order_type order, enum blas_trans_type trans,
		       int m, int n, int kl, int ku, const void *alpha,
		       const void *a, int lda, const void *head_x,
		       const void *tail_x, int incx, const void *beta,
		       void *y, int incy, enum blas_prec_type prec);
void BLAS_cgbmv2_c_s_x(enum blas_order_type order, enum blas_trans_type trans,
		       int m, int n, int kl, int ku, const void *alpha,
		       const void *a, int lda, const float *head_x,
		       const float *tail_x, int incx, const void *beta,
		       void *y, int incy, enum blas_prec_type prec);
void BLAS_cgbmv2_s_c_x(enum blas_order_type order, enum blas_trans_type trans,
		       int m, int n, int kl, int ku, const void *alpha,
		       const float *a, int lda, const void *head_x,
		       const void *tail_x, int incx, const void *beta,
		       void *y, int incy, enum blas_prec_type prec);
void BLAS_cgbmv2_s_s_x(enum blas_order_type order, enum blas_trans_type trans,
		       int m, int n, int kl, int ku, const void *alpha,
		       const float *a, int lda, const float *head_x,
		       const float *tail_x, int incx, const void *beta,
		       void *y, int incy, enum blas_prec_type prec);
void BLAS_zgbmv2_z_d_x(enum blas_order_type order, enum blas_trans_type trans,
		       int m, int n, int kl, int ku, const void *alpha,
		       const void *a, int lda, const double *head_x,
		       const double *tail_x, int incx, const void *beta,
		       void *y, int incy, enum blas_prec_type prec);
void BLAS_zgbmv2_d_z_x(enum blas_order_type order, enum blas_trans_type trans,
		       int m, int n, int kl, int ku, const void *alpha,
		       const double *a, int lda, const void *head_x,
		       const void *tail_x, int incx, const void *beta,
		       void *y, int incy, enum blas_prec_type prec);
void BLAS_zgbmv2_d_d_x(enum blas_order_type order, enum blas_trans_type trans,
		       int m, int n, int kl, int ku, const void *alpha,
		       const double *a, int lda, const double *head_x,
		       const double *tail_x, int incx, const void *beta,
		       void *y, int incy, enum blas_prec_type prec);



int BLAS_fpinfo_x(enum blas_cmach_type cmach, enum blas_prec_type prec);
void BLAS_error(const char *rname, int iflag, int ival, char *form, ...);

#endif /* BLAS_EXTENDED_PROTO_H */
