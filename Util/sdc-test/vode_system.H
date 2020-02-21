#ifndef VODE_SYSTEM_H
#define VODE_SYSTEM_H

#include "RealVector.H"
#include "IntegerVector.H"
#include "RealSparseMatrix.H"

class VodeSystem {
public:
    static const size_t neqs = 3;
    static const size_t nnz = 8;

    IntegerVector<nnz> csr_col_index;
    IntegerVector<neqs+1> csr_row_count;

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    VodeSystem() {
      
        csr_col_index.data[0] = 0;
        csr_col_index.data[1] = 1;
        csr_col_index.data[2] = 2;
        csr_col_index.data[3] = 0;
        csr_col_index.data[4] = 1;
        csr_col_index.data[5] = 2;
        csr_col_index.data[6] = 1;
        csr_col_index.data[7] = 2;

        csr_row_count.data[0] = 0;
        csr_row_count.data[1] = 3;
        csr_row_count.data[2] = 6;
        csr_row_count.data[3] = 8;
	/*
        csr_col_index.data[0] = 0;

        csr_row_count.data[0] = 0;
      */
    }

    AMREX_GPU_DEVICE AMREX_INLINE
    static void evaluate(Real time, RealVector<neqs>& y, RealVector<neqs>& rhs) {

      Real rpar[4];
      Real* rpar_ptr=&(rpar[0]);
      Real a = 0.0688707121;
      rpar[0]=11026.08482;
      rpar[1]=0.008728890037;
      rpar[2]=2.587879236e+10;
      rpar[3]=1/a-1;
      
      RhsFnReal(time,y.data,rhs.data, rpar_ptr, 1);
    }    

    AMREX_GPU_DEVICE AMREX_INLINE
    static void evaluate(Real time, RealVector<neqs>& y, RealVector<neqs>& rhs,
                         RealSparseMatrix<neqs,nnz>& jac) {

      //      printf("before evaluate y=%g f=%g\n",y.data[0],rhs.data[0]);
      
      Real rpar[4];
      Real* rpar_ptr= &(rpar[0]);
      Real a = 0.0688707121;
      for(int i=0;i<neqs;i++)
	{
      rpar[4*i+0]=11026.08482;
      rpar[4*i+1]=0.008728890037;
      rpar[4*i+2]=2.587879236e+10;
      rpar[4*i+3]=1/a-1;
	}
      int one_in=3;
      RhsFnReal(time,y.data,rhs.data, rpar_ptr, one_in);

      
	Real fract = 1e-2*y.data[0];
	//	printf("mid evaluate fract=%g y=%g f=%g\n",fract,y.data[0],rhs.data[0]);
	RealVector<neqs> y_pred;
	RealVector<neqs> rhs_pred;

	for (size_t i = 0; i < neqs; i++) {
	  y_pred.data[i] = y.data[i]+fract;
	  rhs_pred.data[i] = 0.0;
	}

      Real rpar2[4];
      Real* rpar2_ptr= &(rpar2[0]);
      for(int i=0;i<neqs;i++)
	{
      rpar2[4*i+0]=11026.08482;
      rpar2[4*i+1]=0.008728890037;
      rpar2[4*i+2]=2.587879236e+10;
      rpar2[4*i+3]=1/a-1;
	}

      RhsFnReal(time,y_pred.data,rhs_pred.data, rpar2_ptr, one_in);

      jac.data[0] =  0.0e0;
      jac.data[1] =  0.0e0;
      jac.data[2] =  0.0e0;
      
      jac.data[3] =  0.0e0;
      jac.data[4] =  0.0e0;
      jac.data[5] =  0.0e0;

      jac.data[6] = 0.0e0;
      jac.data[7] = 0.0e0;

      jac.data[0] = (rhs_pred.data[0]-rhs.data[0])/fract;
      jac.data[3] = (rhs_pred.data[0]-rhs.data[0])/fract;
      jac.data[7] = (rhs_pred.data[0]-rhs.data[0])/fract;

      printf("after evaluate y=%g f=%g jac=%g\n",y.data[0],rhs.data[0],jac.data[0]);
      
    }
};
#endif
