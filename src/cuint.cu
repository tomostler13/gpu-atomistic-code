// File: cuint.cu
// Author:Tom Ostler
// Last-modified: 09 Jan 2013 10:11:08
#include "../inc/cufields.h"
#include "../inc/cuda.h"
#include "../inc/config.h"
#include "../inc/spins.h"
#include "../inc/mat.h"
#include "../inc/geom.h"
#include "../inc/config.h"
#include "../inc/random.h"
#include "../inc/sim.h"
#include "../inc/intmat.h"
#include "../inc/util.h"
#include "../inc/mat.h"
#include "../inc/cuint.h"
#include "../inc/mf.h"
#include "../inc/tdp.h"
#include "../inc/fields.h"
//Cuda headers
#include <cuda.h>
#include <cuda_runtime.h>
#include <curand.h>
#include <curand_kernel.h>
#include <cufft.h>
//Library headers
#include <cmath>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <iostream>
#ifdef DEBUG
#define CUDA_CALL(x) do { if((x) != cudaSuccess) {\
    printf("Error at %s:%d\n",__FILE__,__LINE__);\
        exit(EXIT_FAILURE);}} while(0)
#else
#define CUDA_CALL(x) (x)
#endif /*DEBUG*/

#define FIXOUT(a,b) a.width(75);a << std::left << b;
#define NR_END 1
#define FREE_ARG char*
#define FREERETURN {free_dmatrix(fjac,1,NLAT,1,NLAT);free_dvector(fvec,1,NLAT);\
	free_dvector(p,1,NLAT);free_ivector(indx,1,NLAT);return;}
#define NSLAT 2
namespace cuint
{
    //constant memory declarations
    __constant__ int NLAT=2;
    __constant__ double J0[NSLAT][NSLAT] = { {0,0},{0,0} };
    __constant__ double D[NSLAT]={0,0};
    __constant__ double conc[NSLAT]={0,0};
    __constant__ double Ctssf=0.0;
    __constant__ double Cchippf=0.0;
    __constant__ double mu[NSLAT];
    __constant__ double lambda=0.0;
	__constant__ double Ms=0.0;
	__constant__ double llbpf=0.0;
	__constant__ double rdt=0.0;
    void copyConstData()
    {
        //copy the contents of the Mean Field arrays to the device
        assert(mf::mfi);
        assert(mat::mi);
        assert(tdp::tdpi);
        assert(sim::simi);
        cudaMemcpyToSymbol(J0, mf::J0.ptr(), mf::J0.size()*sizeof(double));
        cudaMemcpyToSymbol(D,mf::D.ptr(),mf::D.size()*sizeof(double));
        cudaMemcpyToSymbol(conc,mat::conc.ptr(),mat::conc.size()*sizeof(double));
        cudaMemcpyToSymbol(*(&Ctssf),&tdp::tssf,sizeof(double));
        cudaMemcpyToSymbol(*(&Cchippf),&mat::chippf,sizeof(double));
        cudaMemcpyToSymbol(mu,mat::mu.ptr(),mat::mu.size()*sizeof(double));
        cudaMemcpyToSymbol(*(&NLAT),&mat::nlat,sizeof(int));
        cudaMemcpyToSymbol(*(&lambda),&mat::lambda,sizeof(double));
		cudaMemcpyToSymbol(*(&Ms),&mat::Ms,sizeof(double));
		cudaMemcpyToSymbol(*(&llbpf),&sim::llbpf,sizeof(double));
		cudaMemcpyToSymbol(*(&rdt),&sim::rdt,sizeof(double));
    }
    //Numerical recipes arrays that we need
    __device__ float *vector(long nl,long nh)
    {
        float *v;
        v=(float *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(float)));
        return(v-nl+NR_END);
    }
    __device__ double *dvector(long nl,long nh)
    /* allocate a double vector with subscript range v[nl..nh] */
    {
        double *v;

        v=(double *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(double)));
        return v-nl+NR_END;
    }
    __device__ int *ivector(long nl,long nh)
    /* allocate an int vector with subscript range v[nl..nh] */
    {
        int *v;

        v=(int *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(int)));
        return v-nl+NR_END;
    }
    __device__ double **dmatrix(long nrl,long nrh,long ncl,long nch)
    /* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
    {
        long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        double **m;

        /* allocate pointers to rows */
        m=(double **) malloc((unsigned int)((nrow+NR_END)*sizeof(double*)));
        m += NR_END;
        m -= nrl;

        /* allocate rows and set pointers to them */
        m[nrl]=(double *) malloc((unsigned int)((nrow*ncol+NR_END)*sizeof(double)));
        m[nrl] += NR_END;
        m[nrl] -= ncl;

        for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

        /* return pointer to array of pointers to rows */
        return m;
    }
    __device__ void free_ivector(int *v,long nl,long nh)
    /* free an int vector allocated with ivector() */
    {
        free((FREE_ARG) (v+nl-NR_END));
    }
    __device__ void free_dvector(double *v,long nl,long nh)
    /* free a double vector allocated with dvector() */
    {
        free((FREE_ARG) (v+nl-NR_END));
    }
    __device__ void free_dmatrix(double **m,long nrl,long nrh,long ncl,long nch)
    /* free a double matrix allocated by dmatrix() */
    {
        free((FREE_ARG) (m[nrl]+ncl-NR_END));
        free((FREE_ARG) (m+nrl-NR_END));
    }
    __device__ double langevin(double x)
    {
        return((1./tanh(x)-(1./x)));
    }
    //derivative of the langevin function
    __device__ double dlangevin(double x)
    {
        return((1./(x*x))-(1./(sinh(x)*sinh(x))));
    }

    __device__ double (*chiparfp)(double,float,double *);
    __device__ double chipar0(double Tc,float T,double *hmf)
    {
        //Fit Parameter
        const double a0 = 0.8;
        const double a1 =-2.2e-07;
        const double a2 = 1.95e-13;
        const double a3 =-1.3e-17;
        const double a4 =-4e-23;
        const double a5 =-6.5076312364e-32;

        double chi_CGS = 0.0;
        //double chi_SI  = 0.0;
        double chi = 0.0;
        const double PI= 3.1415926535897932384626433832795028841971693993;;
        //double factor =  0.75947907;

        const double factor=Tc-T;

        if(T<Tc) chi_CGS =(a0/660.*Tc)/(4.*PI)/(factor)+a1*factor+ a2*factor*factor*factor+a3*factor*factor*factor*factor+ a4*factor*factor*factor*factor*factor*factor+ a5*factor*factor*factor*factor*factor*factor*factor*factor*factor;
        else chi_CGS = (1.1*1.4/660.*Tc)/(4*PI)/(T-Tc);
        //chi_SI = 4*PI*chi_CGS;     // CHI_SI = 4*PI Chi_CGS
        //chi = chi_SI*4*A_FePt*A_FePt*A_FePt/MU_S_FePt/MU_0;
        chi = chi_CGS*9.54393845712027+0.308e-14; // (Tesla)
        return(chi);
    }

    __device__ double chiparmf(double Tc,float T,double *hmf)
    {
        const double lb=1./(1.38e-16*double(T));
        const double Bp=dlangevin(hmf[0]*lb);
        const double J=J0[0][0];
        //There is a factor of 10^4 here because
        //1T=10^4 Oe
        return((mu[0]*9.27e-21*1e4/J)*Bp*lb*J/(1.0-Bp*lb*J));

    }



    __device__ double (*chiperpfp)(double,float,const double,const double);
    __device__ double chiperp0(double Tc,float T,const double me,const double chipar)
    {
        //Fit Parameter
        const double a0 = 0.00211549427182711;
        const double a1 = 0.110224660661792;
        const double a2 = -0.855153260915204;
        const double a3 = 3.42088365387997;
        const double a4 = -7.84821585896818;
        const double a5 = 10.3247035469514;
        const double a6 = -6.85608273303224;
        const double a7 = 0.797453198330591;
        const double a8 = 1.53787854178089;
        const double a9 = -0.627128148404525;

        double chi_CGS = 0.0;
        //double chi_SI  = 0.0;
        double chi     = 0.0;
        const double PI= 3.1415926535897932384626433832795028841971693993;;
        if(T<(1.065*Tc)) chi_CGS = a0+ a1*pow(pow(((1.068*Tc-T)/(1.068*Tc)),0.5),2.)+ a2*pow((((1.068*Tc)-T)/(1.068*Tc)),2.)+ a3*pow((((1.068*Tc)-T)/(1.068*Tc)),3.)+ a4*pow((((1.068*Tc)-T)/(1.068*Tc)),4.) + a5*pow((((1.068*Tc)-T)/(1.068*Tc)),5.) + a6*pow((((1.068*Tc)-T)/(1.068*Tc)),6.) + a7*pow((((1.068*Tc)-T)/(1.068*Tc)),7.)+ a8*pow((((1.068*Tc)-T)/(1.068*Tc)),8.) + a9*pow((((1.068*Tc)-T)/(1.068*Tc)),9.);
        else chi_CGS = (0.8*1.4/660.*Tc)/(4*PI)/(T-Tc);

        //chi_SI = 4*PI*chi_CGS;     // CHI_SI = 4*PI Chi_CGS
        //chi = chi_SI*4*A_FePt*A_FePt*A_FePt/MU_S_FePt/MU_0;
        chi = chi_CGS*9.54393845712027; // (Tesla)

        return(chi*Ctssf);

    }

	__device__ double (*exchstifffp)(double,float);
    __device__ double exchstiffNULL(double Tc,float T)
    {
        return(0);
    }
	__device__ double exchstiff0(double Tc,float T)
	{
		const double ex0 = 3.90858143659231e-13;
        const double ex1 = 5.65571902911896e-11;
        const double ex2 = -1.11221431025254e-10;
        const double ex3 = 1.67761522644194e-10;
        const double ex4 = -1.38437771856782e-10;
        const double ex5 = 4.6483423884759e-11;

		const double factor=(Tc-double(T))/Tc;

		if(T<Tc)
		{
			return(ex0+ex1*factor+ex2*factor*factor+ex3*factor*factor*factor+ex4*factor*factor*factor*factor+ex5*factor*factor*factor*factor*factor);
		}
		else
		{
			return(0.0);
		}

	}
	//function pointer for anisotropy
	__device__ void (*anisfp)(float *,double *,double);
	__device__ void anis0(float *h,double *s,double chiperp)
	{
		return;
	}
	__device__ void anis1(float *h,double *s,double chiperp)
	{
		h[0]+=-s[0]/chiperp;
		h[1]+=-s[1]/chiperp;
	}
	//function pointer for applied field
	__device__ void (*appliedfp)(float *,float,float,float);
	__device__ void applied0(float *h,float Bx,float By,float Bz)
	{
        h[0]+=Bx;
        h[1]+=By;
        h[2]+=Bz;
	}

    __device__ double chiperpmf(double Tc,float T,const double me,const double chipar)
    {
        if(T<Tc)
        {
            const double oome=1./me;
            return(Cchippf*oome*oome);
        }
        else
        {
            return(chipar);
        }

    }

    __device__ double (*mefp)(double,float,double*,double);
    __device__ double me0(double Tc,float T,double *hmf,double meprev)
    {
        const double Temp=double(T);
        const double TCmT=Tc-Temp;
        const double Tc_m_T_o_Tc=TCmT/Tc;
        const double factor=(1.068*Tc-Temp)/(1.068*Tc);

        if(Temp<Tc)
        {
            return(1.3*sqrt(Tc_m_T_o_Tc)-0.12*Tc_m_T_o_Tc-0.51*Tc_m_T_o_Tc*Tc_m_T_o_Tc + 0.37*Tc_m_T_o_Tc*Tc_m_T_o_Tc*Tc_m_T_o_Tc - 0.01*Tc_m_T_o_Tc*Tc_m_T_o_Tc*Tc_m_T_o_Tc*Tc_m_T_o_Tc-0.03*Tc_m_T_o_Tc*Tc_m_T_o_Tc*Tc_m_T_o_Tc*Tc_m_T_o_Tc*Tc_m_T_o_Tc*Tc_m_T_o_Tc*Tc_m_T_o_Tc*Tc_m_T_o_Tc);
        }
        else
        {
            return(0);
        }
    }
    //NR functions
    __device__ void usrfun(double *hmf,double *msub,int n,double *fvec,double **fjac,double& lbeta)
    {
        //calculate the sublattice magnetization Langevin function
        for(unsigned int i = 0 ; i < n ; i++)
        {
            hmf[i]=D[i];
            for(unsigned int j = 0 ; j < n ; j++)
            {
                hmf[i] += conc[j]*J0[i][j]*msub[j+1];
            }
//            hmfg(atom,i)=hmf[i];
            fvec[i+1]=msub[i+1]-langevin(hmf[i]*lbeta);
        }
        double I[NSLAT][NSLAT];
        for(unsigned int  i = 0 ; i < NLAT ; i++)
        {
            for(unsigned int j = 0 ; j < NLAT ; j++)
            {
                I[i][j]=0.0;
                if(i==j)
                {
                    I[i][j]=1.0;
                }
            }
        }
        //calculate the Jacobian matrix
        for(unsigned int i = 0 ; i < n ; i++)
        {
            for(unsigned int j = 0 ; j < n ; j++)
            {
                fjac[i+1][j+1] = I[i][j];
                fjac[i+1][j+1] -= (  (2.0*D[i]*I[i][j]+conc[j]*J0[i][j])*lbeta*dlangevin(hmf[i]*lbeta));
            }
        }

    }
    __device__ void ludcmp(double **a, int n, int *indx, double& d)
    {
        int imax;
        double big=0.0,dum=0.0,sum=0.0,temp=0.0;
        double *vv;
        //vv stores the implicit scaling of each row
        vv=dvector(1,n);//this declaration requires the NR utilities
        d=1.0;//No row interchanges yet.
        for (int i=1;i<=n;i++)
        {
            //Loop over rows to get the implicit scaling information
            big=0.0;
            for(int j=1;j<=n;j++)
            {
                if ((temp=fabs(a[i][j])) > big) big=temp;
            }
            //No nonzero largest element
            if (big == 0.0)
            {
//                std::cerr << "Error: singular matrix in routine ludcmp" << std::endl;
            }
            vv[i]=1.0/big;//save the scaling
        }//end of i loop
        for (int j=1;j<=n;j++)
        {
            //This is the loop over the columns of Crout's method.
            for (int i=1;i<j;i++)
            {
                //This is equation (2.3.12) except for i=j from NR in C
        //printf("\ni=%d",i);
                sum=a[i][j];
                for (int k=1;k<i;k++)
                {     /*printf("\nk=%d",k);*/
                    sum -= a[i][k]*a[k][j];
                }
                    a[i][j]=sum;
            }
            big=0.0;//Initialize for the search for largest pivot element.
            for (int i=j;i<=n;i++)
            {
                //This is the i=j part of equation (2.3.12) and i = k +1 . . . N, of equation (2.3.12).
                sum=a[i][j];
                for (int k=1;k<j;k++)
                {
                    sum -= a[i][k]*a[k][j];
                }
                a[i][j]=sum;
                //Is the figure of merit for the pivot better than the best so far?
                if ( (dum=vv[i]*fabs(sum)) >= big)
                {
                    big=dum;
                    imax=i;
                }
            }
            if (j != imax)//Do we need to interchange rows?
            {
                //printf("%d\t%d\n",j,imax);
                for (int k=1;k<=n;k++)
                {
                    //Yes, do so...
                    dum=a[imax][k];
                    a[imax][k]=a[j][k];
                    a[j][k]=dum;
                }
                d = -(d);//...and change the parity of d.
                vv[imax]=vv[j];//Also interchange the scale factor.
            }
            indx[j]=imax;
            if (a[j][j] == 0.0)
            {
                a[j][j]=1e-40;
            }
            //If the pivot element is zero the matrix is singular (at least to the precision of the
            //algorithm). For some applications on singular matrices, it is desirable to substitute
            //TINY for zero.

            if (j != n)
            {
            //printf("Strangely satisfied...%d\t%d\n",j,n);
                dum=1.0/(a[j][j]);
                for (int i=j+1;i<=n;i++)
                {
                    a[i][j] *= dum;
                }
            }
        }//go back for the next column in the reduction
        free_dvector(vv,1,n);
    }
    __device__ void lubksb(double **a, int n, int *indx, double b[])
    {
        int i,ii=0,ip,j;
        double sum=0.0;

        for (i=1;i<=n;i++) {
            ip=indx[i];
            sum=b[ip];
            b[ip]=b[i];
            if (ii)
                for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
            else if (sum) ii=i;
            b[i]=sum;
        }
        for (i=n;i>=1;i--) {
            sum=b[i];
            for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
            b[i]=sum/a[i][i];
        }
    }
    /* (C) Copr. 1986-92 Numerical Recipes Software ?421.1-9. */

    //the arguements are as follows:
    //ntrail - number of trial steps
    //x - the sublattice magnetizations
    //n - the number of sublattices
    //tolx - convergence parameter
    //tolf - function sum convergence parameter
    //T - temperature
    __device__ void mnewt(int& ntrial, double *x, double& tolx, double& tolf, double T , double *hmf)
    {
        double beta = 1./(1.38e-16*T);
        int k,i,*indx;
        double errx,errf,d,*fvec,**fjac,*p;

        indx=ivector(1,NLAT);
        for(unsigned int i = 1 ; i <= NLAT ; i++)
        {
            indx[i]=0;
        }
        p=dvector(1,NLAT);
        fvec=dvector(1,NLAT);
        fjac=dmatrix(1,NLAT,1,NLAT);

        for (k=1;k<=ntrial;k++) {
            usrfun(hmf,x,NLAT,fvec,fjac,beta);
            errf=0.0;
            for (i=1;i<=NLAT;i++) errf += fabs(fvec[i]);
            if (errf <= tolf) FREERETURN
            for (i=1;i<=NLAT;i++) p[i] = -fvec[i];
            ludcmp(fjac,NLAT,indx,d);
            lubksb(fjac,NLAT,indx,p);
            errx=0.0;
            for (i=1;i<=NLAT;i++) {
                errx += fabs(p[i]);
                x[i] += p[i];
            }
            if (errx <= tolx) FREERETURN
        }
        FREERETURN

        #undef NRANSI
        #undef FREERETURN
    }
    //mean field calculation of equilibrium magnetization
    __device__ double memf(double Tc,float T,double *hmf,double meprev)
    {
        double tolx=1e-5,tolf=1e-5;
        double *x;
        x=dvector(1,NSLAT);
        //initial guess of magnetization
        for(unsigned int j = 0 ; j <= NLAT ; j++)
        {
            x[j]=meprev;
        }
        //Newton-Rhapson method (from NR)
        int ntrial=1000;
        mnewt(ntrial,x,tolx,tolf,double(T),hmf);
        return(x[1]);
    }
    __global__ void CSetFunctionPointer(int metype,int mefunc,int chipartype,int chiparfunc,int chiperptype,int chiperpfunc,int exchstifftype,int exchstifffunc,unsigned int anist,unsigned int ff)
    {
        const int i = blockDim.x*blockIdx.x + threadIdx.x;
        if(i==0)
        {
			if(anist==0)
			{
				//no anisotropy
				anisfp=anis0;
			}
			else if(anist==1)
			{
				//hard axis -(m_x e_x - m_y e_y)/chiperp
				anisfp=anis1;
			}
			//applied field
			if(ff==0)
			{
				appliedfp=applied0;
			}

            if(metype>=0)
            {
                if(mefunc==0)
                {
                    mefp=me0;
                }
            }
            else
            {
                mefp=memf;
            }
            if(chipartype>=0)
            {
                if(chiparfunc==0)
                {
                    chiparfp=chipar0;
                }
            }
            else
            {
                chiparfp=chiparmf;
            }
            if(chiperptype>=0)
            {
                if(chiperpfunc==0)
                {
                    chiperpfp=chiperp0;
                }
            }
            else
            {
                chiperpfp=chiperpmf;
            }
			if(exchstifftype>=0)
			{
				if(exchstifffunc==0)
				{
					exchstifffp=exchstiff0;
				}
			}
			else
			{
                if(exchstifftype==-1)
                {
                    exchstifffp=exchstiffNULL;
                }
			}
        }
    }

    __global__ void CHeun1(int N,float Bx,float By,float Bz,double Tc,float *CHDemag,float *Cfspin,double *Cspin,double *Cespin,float *CTemp,int *Cxadj,int *Cadjncy,double *CsurfArea,double *Csigma,float *Crand,double *Cfn)
    {
        register const int i = blockDim.x*blockIdx.x + threadIdx.x;
        if(i<N)
        {
			const float lrn[6]={Crand[6*i],Crand[6*i+1],Crand[6*i+2],Crand[6*i+3],Crand[6*i+4],Crand[6*i+5]};
            float h[3]={CHDemag[3*i],CHDemag[3*i+1],CHDemag[3*i+2]};
            double s[3]={Cspin[3*i],Cspin[3*i+1],Cspin[3*i+2]};
            float T=CTemp[i];
            double hmf[NSLAT];
            const double m_squared=s[0]*s[0]+s[1]*s[1]+s[2]*s[2];
			const double oomi2=1./(m_squared);
            for(unsigned int j = 0 ; j < NSLAT ; j++){hmf[j]=0.0;}
            const double me=mefp(Tc,T,hmf,sqrt(m_squared));
            const double m_e_squared=me*me;
            const double chipar=chiparfp(Tc,T,hmf);
            const double chiperp=chiperpfp(Tc,T,me,chipar);
			const double exchstiffness=exchstifffp(Tc,T);
            const double AlpPar=lambda*0.6666666666666666666666666*double(T)/Tc;
            double AlpPerp=0;
            const double CW1pf=sqrt(double(T)*AlpPar);
            double pf=0.0;
            if(T<Tc)
            {
                pf=(0.5/(chipar))*(1.0-m_squared/m_e_squared);
                AlpPerp=lambda*(1.0-T/(3.0*Tc));
            }
            else
            {
                pf=-(1.0/(chipar))*(1.0 + Tc*3.0*m_squared/(5.0*(T-Tc)));
                AlpPerp=AlpPar;
            }
            const double CW2pf = sqrt(double(T)*(AlpPerp-AlpPar));
			anisfp(h,s,chiperp);
			appliedfp(h,Bx,By,Bz);
			//longitudinal field component
			for(unsigned int j = 0 ; j < 3 ; j++)
			{
				h[j]+=pf*s[j];
			}


			//exchange field calculation
			for(unsigned int n = Cxadj[i] ; n < Cxadj[i+1] ; n++)
			{
				const unsigned int neighbour=Cadjncy[n];
				//lookup spin information
				float sj[3]={Cfspin[3*neighbour],Cfspin[3*neighbour+1],Cfspin[3*neighbour+2]};
				double area=CsurfArea[n];
				double exchpf=2.0*exchstiffness/(m_e_squared*Ms);
				for(unsigned int j = 0 ; j < 3 ; j++)
				{
					h[j]+=exchpf*((sj[j]-s[j])/area);
				}
			}

			const double sdoth = double(s[0]*h[0]+s[1]*h[1]+s[2]*h[2]);
			const double tpf=CW2pf*Csigma[i];
			const float hpt[3]={h[0]+tpf*lrn[0],h[1]+tpf*lrn[1],h[2]+lrn[2]};
			const float sxh[3]={s[1]*h[2] - s[2]*h[1],s[2]*h[0]-s[0]*h[2],s[0]*h[1]-s[1]*h[0]};
            const float sxhpt[3]={s[1]*hpt[2] - s[2]*hpt[1],s[2]*hpt[0]-s[0]*hpt[2],s[0]*hpt[1]-s[1]*hpt[0]};
            const float sxsxhpt[3]={s[1]*sxhpt[2]-s[2]*sxhpt[1],s[2]*sxhpt[0]-s[0]*sxhpt[2],s[0]*sxhpt[1]-s[1]*sxhpt[0]};
			const double tpf1=CW1pf*Csigma[i];
			double lfn[3]={0,0,0};
			for(unsigned int j = 0 ; j < 3 ; j++)
			{
				//use lrn[3+j] because we want second half of lrn array for second Weiner process
				lfn[j] = llbpf*(double(sxh[j]) + (AlpPerp*oomi2)*double(sxsxhpt[j]) - (AlpPar*oomi2)*sdoth*double(s[j])) + lrn[3+j]*tpf1;
				Cfn[3*i+j]=lfn[j];
				Cespin[3*i+j]=s[j]+lfn[j]*rdt;
			}

        }
    }

    __global__ void CHeun2(int N,float Bx,float By,float Bz,double Tc,float *Cfspin,double *Cspin,float *CHDemag,float *CTemp,float *Crand,double *Cespin,double *Cfn,int *Cxadj,int *Cadjncy,double *CsurfArea,double *Csigma)
    {
        register const int i = blockDim.x*blockIdx.x + threadIdx.x;
        if(i<N)
        {
            float lrn[6]={Crand[6*i],Crand[6*i+1],Crand[6*i+2],Crand[6*i+3],Crand[6*i+4],Crand[6*i+5]};
            float h[3]={CHDemag[3*i],CHDemag[3*i+1],CHDemag[3*i+2]};
            double s[3]={Cespin[3*i],Cespin[3*i+1],Cespin[3*i+2]};
            float T=CTemp[i];
            double hmf[NSLAT];
            const double m_squared=s[0]*s[0]+s[1]*s[1]+s[2]*s[2];
			const double oomi2=1./(m_squared);
            for(unsigned int j = 0 ; j < NSLAT ; j++){hmf[j]=0.0;}
            const double me=mefp(Tc,T,hmf,sqrt(m_squared));
            const double m_e_squared=me*me;
            const double chipar=chiparfp(Tc,T,hmf);
            const double chiperp=chiperpfp(Tc,T,me,chipar);
			const double exchstiffness=exchstifffp(Tc,T);
            const double AlpPar=lambda*0.6666666666666666666666666*double(T)/Tc;
            double AlpPerp=0;
            const double CW1pf=sqrt(double(T)*AlpPar);
            double pf=0.0;
            const double fn[3]={Cfn[3*i],Cfn[3*i+1],Cfn[3*i+2]};
            double fnp1[3]={0,0,0};

            if(T<Tc)
            {
                pf=(0.5/(chipar))*(1.0-m_squared/m_e_squared);
                AlpPerp=lambda*(1.0-T/(3.0*Tc));
            }
            else
            {
                pf=-(1.0/(chipar))*(1.0 + Tc*3.0*m_squared/(5.0*(T-Tc)));
                AlpPerp=AlpPar;
            }
            const double CW2pf = sqrt(double(T)*(AlpPerp-AlpPar));
			anisfp(h,s,chiperp);
			appliedfp(h,Bx,By,Bz);
			//longitudinal field component
			for(unsigned int j = 0 ; j < 3 ; j++)
			{
				h[j]+=pf*s[j];
			}


			//exchange field calculation
			for(unsigned int n = Cxadj[i] ; n < Cxadj[i+1] ; n++)
			{
				const unsigned int neighbour=Cadjncy[n];
				//lookup spin information
                double sj[3]={Cespin[3*neighbour],Cespin[3*neighbour+1],Cespin[3*neighbour+2]};
                double area=CsurfArea[n];
                double exchpf=2.0*exchstiffness/(m_e_squared*Ms);
				for(unsigned int j = 0 ; j < 3 ; j++)
				{
					h[j]+=exchpf*((sj[j]-s[j])/area);
				}
			}
			const double sdoth = double(s[0]*h[0]+s[1]*h[1]+s[2]*h[2]);
			const double tpf=CW2pf*Csigma[i];
			const float hpt[3]={h[0]+tpf*lrn[0],h[1]+tpf*lrn[1],h[2]+lrn[2]};
			const float sxh[3]={s[1]*h[2] - s[2]*h[1],s[2]*h[0]-s[0]*h[2],s[0]*h[1]-s[1]*h[0]};
            const float sxhpt[3]={s[1]*hpt[2] - s[2]*hpt[1],s[2]*hpt[0]-s[0]*hpt[2],s[0]*hpt[1]-s[1]*hpt[0]};
            const float sxsxhpt[3]={s[1]*sxhpt[2]-s[2]*sxhpt[1],s[2]*sxhpt[0]-s[0]*sxhpt[2],s[0]*sxhpt[1]-s[1]*sxhpt[0]};
			const double tpf1=CW1pf*Csigma[i];

            double os[3];
            for(unsigned int j = 0 ; j < 3 ; j++)
			{
                os[j]=Cspin[3*i+j];
				//use lrn[3+j] because we want second half of lrn array for second Weiner process
				fnp1[j] = llbpf*(double(sxh[j]) + (AlpPerp*oomi2)*double(sxsxhpt[j]) - (AlpPar*oomi2)*sdoth*double(s[j])) + lrn[3+j]*tpf1;
                os[j]+=(fn[j]+fnp1[j])*rdt*0.5;
                Cspin[3*i+j]=os[j];
                //store a floating point copy of the spin array as well
                Cfspin[3*i+j]=float(os[j]);
            }

        }
    }
}
