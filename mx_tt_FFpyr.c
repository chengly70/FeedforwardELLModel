#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "mex.h"

/* Box-Muller; faster than using sin,cos */
 void z_randn(double *normal_rv, double *nrm_rv2)
 {
	 double u, v, s=2.0, z, z2;
	 while(s>=1 || s==0){
		 u=(double)(rand())/RAND_MAX*2-1;
		 v=(double)(rand())/RAND_MAX*2-1;
		 s=u*u+v*v;
	 }
	 z=sqrt(-2.0*log(s)/s)*u;
	 z2=z/u*v;
	 
	 *normal_rv=z;
	 *nrm_rv2=z2;
 }

void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]) 
{
/* 
[nu_e,f_rate,f_isi,t_vc]=mx_tt_FFpyr(I_aff,I_ie,I_ee,tV,W_ie,g_ie,W_ee,g_ee,Thres,qPreF);
called in Matlab script
must compile first: >> mex mx_tt_FFpyr.c
ALL TIMES IN SECONDS
The LIF model
 nu_e = Ne x Lts (psth for each neuron); f_rate = Ne x 1 (f-rate); (t_vc,f_isi) is POPULATION isi
*/  

/*  Check for proper number of arguments */
if (nrhs !=  10) {
	mexErrMsgTxt("10 input arguments required.");
} else if (nlhs > 4) {
	mexErrMsgTxt("Too many output arguments.");
}

int Ne, Nei, i, j, k, ind1, N=1000;
double tau_m=0.01, Isyn=-0.5, Esyn=6.5, v_min=-0.5, ref_s=.001, sig_BP=0.75, tau_eta=0.005, sqtn; /* intrins parms  */
double rndUnif, isi_max=7.0;

int l_wie, l_wee, len_isi, numSpk=0;
double *I_aff, *I_ie, *I_ee, *tV, *Wie, *g_iep, *Wee, *g_eep, *Thres, *qPreF; /* passed-in */
double g_ie, Gie, g_ee, Gee; /* (some) passed-in */
	
double randc1, randc2=0.0, dt=0.0001, dts=0.001; /* used in for-loop; */
int Lt, Lts, spct, nrnSpace; /* 1 s for each realiz */
  
/* outputs */
double *nu_e, *f_rate, *f_isi, *t_vc;
	
/* Using mxGetScalar to retreive input arguments */
I_aff=mxGetPr(prhs[0]);     /* Lt x 1 vector */
I_ie=mxGetPr(prhs[1]);      /* Nei x Lt matrix */
I_ee=mxGetPr(prhs[2]);      /* Nei x Lt matrix */
tV=mxGetPr(prhs[3]);        /* Lt x 1 vector */
Wie=mxGetPr(prhs[4]);       /* Ne x l_wie matrix */
g_iep=mxGetPr(prhs[5]);     /* scalar */
Wee=mxGetPr(prhs[6]);       /* Ne x l_wee matrix */
g_eep=mxGetPr(prhs[7]);     /* scalar */
Thres=mxGetPr(prhs[8]);    /* Ne x 1 vector */
qPreF=mxGetPr(prhs[9]);    /* Ne x 1 vector */
    
g_ie=g_iep[0];
g_ee=g_eep[0];

Lt=mxGetNumberOfElements(prhs[0]);    /* length of vector */
Nei=mxGetM(prhs[1]);
Ne=mxGetM(prhs[4]);    /* # rows in connect matrix */
l_wie=mxGetN(prhs[4]); /* # cols in conn matrix; smaller than Ne and Nei to save time */
l_wee=mxGetN(prhs[6]); /* # cols in conn matrix; smaller than Ne and Nei to save time */

    if (Lt!=mxGetN(prhs[1]) || Lt!=mxGetN(prhs[2]) || Lt!=mxGetNumberOfElements(prhs[3])) {
        mexErrMsgTxt("Size of I_ie,I_ee or tV does not match I_aff (time, Lt).");
    }
    if (Ne!=mxGetNumberOfElements(prhs[8]) || Ne!=mxGetM(prhs[6])) {
        mexErrMsgTxt("Length of Thres (Ne) does not match conn. or Row(Wee)~=Row(Wie).");
    }
    
//Calc the len_isi AND Lts before allocating output
len_isi=(int)(isi_max/dts+1);
    
spct=(int)(dts/dt); /* save spikes every spct-time steps */
Lts=Lt/spct;   /* only sample ms; !Lts is total time (ms) */
    

/* OUTPUTs ... */
plhs[0]=mxCreateDoubleMatrix( Ne, Lts, mxREAL); /* only mxREAL part */
nu_e = mxGetPr(plhs[0]);
plhs[1]=mxCreateDoubleMatrix( Ne, 1, mxREAL); /* only mxREAL part */
f_rate = mxGetPr(plhs[1]);
plhs[2]=mxCreateDoubleMatrix( len_isi, 1, mxREAL);
f_isi = mxGetPr(plhs[2]);
plhs[3]=mxCreateDoubleMatrix( len_isi, 1, mxREAL);
t_vc = mxGetPr(plhs[3]);

//Setting (t_vc) ISI state variable now
for (j=0; j<len_isi; j++) {
    t_vc[j]=j*dts; //t_vc=(0:dts:isi_max)';
}
    
/* number of time steps before out of refractory */
nrnSpace=(int)(ref_s/dt);


/* used in for-loop; ODEs */
double vlt_E[Ne];
double etaV[Ne]; /* indiv noise */
int t_lstSpk[Ne]; /* indx time of last spike */

int W_ie[Ne*l_wie];
int W_ee[Ne*l_wee];
    
/*  initialization */
srand ( time(NULL) ); /* seeing random # gen. */
sqtn=1/sqrt(tau_eta*dt);
for (i=0; i<Lts; i++) {
    for (j=0; j<Ne; j++) {
        nu_e[i*Ne+j]=0.0;
    }
    if (i==0) { //only do it once
        for (j=0; j<Ne; j++) {
            f_rate[j]=0.0;
        }
    }
}
    
for (j=0; j<l_wie; j++) {
	for (k=0; k<Ne; k++)
		W_ie[j*Ne+k]=(int)Wie[j*Ne+k];
}
for (j=0; j<l_wee; j++) {
    for (k=0; k<Ne; k++)
        W_ee[j*Ne+k]=(int)Wee[j*Ne+k];
    }
    
srand ( time(NULL) ); /* seeing random # gen. */

/* set for FIRST run only; then Initi Cond. determined from t(end) of previous run */
for (j=0; j<Ne; j++) {
    vlt_E[j] = (double)(rand())/RAND_MAX; /* rand unif i.c. */
	etaV[j] = ((double)(rand())/RAND_MAX-0.5)*2; /* rand unif i.c. [-1,1] */
    t_lstSpk[j]=-nrnSpace;
}
		

/*  MAIN....N realizations. */
for (i=0; i<N; i++){
    if (i>0) {
        for (j=0; j<Ne; j++) { /* SET INITIAL CONDITIONS */
            t_lstSpk[j]=t_lstSpk[j]-Lt;
        }
    }
	
/* start of time-loop */
for (j=0; j<Lt; j++){
    
	for (k=0; k<Ne; k++){		
		if((j-t_lstSpk[k]) >= nrnSpace ){ /* v change only if not in refrac */
			Gie=0.0; /* recalc */
            ind1=0;
            while (ind1<l_wie && W_ie[ind1*Ne+k]!=-1) {
                Gie+=(g_ie*I_ie[j*Nei+ W_ie[ind1*Ne+k] ]*qPreF[k] );
                ind1++;
            }
            Gee=0.0; /* recalc */
            ind1=0;
            while (ind1<l_wee && W_ee[ind1*Ne+k]!=-1) {
                Gee+=(g_ee*I_ee[j*Nei+ W_ee[ind1*Ne+k] ]*qPreF[k] );
                ind1++;
            }
            
            vlt_E[k]=vlt_E[k]+dt/tau_m*(-vlt_E[k]+I_aff[j]-Gie*(vlt_E[k]-Isyn)-Gee*(vlt_E[k]-Esyn)+sig_BP*etaV[k]);
			if (vlt_E[k]<v_min) { /* lower barrier */
				vlt_E[k]=2*v_min-vlt_E[k];
			}
		}
        
		if(k%2 == 0){
			z_randn(&randc1,&randc2); /*k starts at 0, so randc2 will be assigned */
			etaV[k]=etaV[k]+dt*(-etaV[k]/tau_eta+sqtn*randc1);
		}
		else
			etaV[k]=etaV[k]+dt*(-etaV[k]/tau_eta+sqtn*randc2);
		
		/* spiked */
		if(vlt_E[k]>=Thres[k]){
			vlt_E[k]=0.0;		 // reset voltage
			nu_e[j/spct*Ne+k]+=1; // record spike; rely on integer division
			
            ind1=(int)round(dt*(j-t_lstSpk[k])/dts);
            if (ind1>len_isi) {
                ind1=len_isi;
            }
            f_isi[ind1] += 1;        // update f_isi; rely on integer division
            numSpk += 1;            //count of how many spikes total
            
			t_lstSpk[k]=j;
		}

	}//end of k-loop (over Ne)
    
}/* ending j-loop (time) */

	
/* UPDATE running sum so unbiased-estim of var,cov, etc */
	
} /* ending i-loop (realizations) */

for (j=0; j<(Lts); j++) { /* normalize now so don't have to keep running sum [not often spiking] */
	for (k=0; k<Ne; k++) {
		nu_e[j*Ne+k]=nu_e[j*Ne+k]/(N*dts);
	}
}
    
for (k=0; k<Ne; k++) {
    for (j=0; j<Lts; j++) {
        f_rate[k] += nu_e[j*Ne+k]/Lts;
    }
}
    
for (j=0; j<len_isi; j++) {
    f_isi[j] /= (dts*numSpk); //normalize so integrates to 1
}
	
return;
                
}
