/* ctomo.c - to invert travel times , this version also updates the dtimes in phases.b */
/* this version for Pn, it accounts for the ray offset */
/* this does gauss siedel iteration */
/* & this version does not do the scaleing */
/* & this version bins the stations */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
/*next to lines are for mmap()*/
#include <sys/types.h>
#include <sys/mman.h>
/*next two for lseek()*/
#include <sys/types.h>
#include <unistd.h>


/* inversion parameters*/
#define NLSQR	25			/*number of lsqr iterations*/
#define NBOOTV	0			/*number of bootstrap iterations for variance*/
#define WEIGHT	500		/*amplitude of laplacian damping*/
#define ANIWEI	500		/*include anisotropy, weight for anisotropy*/

#define PIO180	0.01745329252


float xmin=0.0,xmax=360.0,ymin=-90.0,ymax=90.0;
float xmesh=1.0,ymesh=1.0;
float weight=WEIGHT, aniwei=ANIWEI;
int nlsqr=NLSQR, nbootv=NBOOTV;
int nx, ny;
int npha;
int nstat, nevn;
struct rec *phase;
double *utemp,*xtemp,*vtemp;


#include "data14.h"

double sqrt(), fabs(), cos(), sin(), exp2(), pow(), amin(), amax(), exp();
int intt();
/*Compiler error*/
void *malloc();
//caddr_t mmap();*/


struct locations *station, *event;
struct solution *sol;

main(argc,argv)
int argc; char **argv;
{
	int             i, j, k,l, ii, jj;
	int             iin = 0;

	double          wei;
	double          res1, res2, res3, res4;
	double			res2slo,res2cos,res2sin;
	double          time;
	int             fdrec, fdsol;
	int             fdevn, fdstn;
	FILE           *fopen();
	double **u, **x, **v;
	int nstat2,nevn2,nsol2;
	double ran1();
	double temp, temp2;

	int ncell;
	int xcell[999],ycell[999];
	double lencell[999];
	double length;
	int ntimes,nparams;
	double anorm;

	struct rec *phase2;

	if(getarg(argc,argv,"-region","%f/%f/%f/%f",&xmin,&xmax,&ymin,&ymax)<4) errprint("ctomo: cant getarg region");
	if(getarg(argc,argv,"-mesh","%f/%f",&xmesh,&ymesh)<2) errprint("ctomo: cant getarg xmesh");
	getarg(argc,argv,"-weight","%f",&weight);
	getarg(argc,argv,"-aniwei","%f",&aniwei);
	getarg(argc,argv,"-nlsqr","%d",&nlsqr);
	getarg(argc,argv,"-nbootv","%d",&nbootv);
	if(nbootv!=NBOOTV&&NBOOTV*nbootv==0) errprint ("ctomo: must recompile with NBOOTV reset");
	if(aniwei!=ANIWEI&&ANIWEI*aniwei==0) errprint ("ctomo: must recompile with ANIWEI reset");
	fprintf(stderr,"xmin,xmax,ymin,ymax,xmesh,ymesh= %f %f %f %f %f %f\n",xmin,xmax,ymin,ymax,xmesh,ymesh);
	fprintf(stderr,"weight,aniwei,nlsqr,nbootv= %f %f %d %d\n",weight,aniwei,nlsqr,nbootv);

	nx= (xmax-xmin)*xmesh;
	ny= (ymax-ymin)*ymesh;
	fprintf(stderr,"nx,ny= %d %d\n",nx,ny);

	/* open all files*/
	fdrec = open("../data/phases.b", 2);
	fdevn = open("../data/events.b", 2);
	fdstn = open("../data/stations.b", 2);
	fdsol = creat("../data/sol.b", 0644);
	fprintf(stderr,"fdrec,fdevn,fdstn,fdsol= %d %d %d %d\n", fdrec, fdevn, fdstn, fdsol);

	/* mmap memory spaces */
	phase = (struct rec *) mmap(0,lseek(fdrec,0L,SEEK_END),PROT_READ|PROT_WRITE, MAP_SHARED,fdrec,0);
	npha = lseek(fdrec,0L,SEEK_END)/sizeof(struct rec);
	station = (struct locations *) mmap(0,lseek(fdstn,0L,SEEK_END),PROT_READ|PROT_WRITE, MAP_SHARED,fdstn,0);
	nstat = lseek(fdstn,0L,SEEK_END)/sizeof(struct locations);
	event = (struct locations *) mmap(0,lseek(fdevn,0L,SEEK_END),PROT_READ|PROT_WRITE, MAP_SHARED,fdevn,0);
	nevn = lseek(fdevn,0L,SEEK_END)/sizeof(struct locations);
	fprintf(stderr,"station,event,phase= %d %d %d\n",station,event,phase);
	fprintf(stderr,"nstat,nevn,npha = %d %d %d\n", nstat, nevn, npha);

	/* allocate storage spaces */
	sol = (struct solution *) malloc(nx*ny * sizeof(struct solution));
	phase2 = (struct rec *) malloc(lseek(fdrec,0L,SEEK_END));
	fprintf(stderr,"sol= %d\n", sol);
	x = (double **) malloc((3*nx*ny+nstat+nevn)*8);
	v = (double **) malloc((3*nx*ny+nstat+nevn)*8);
	u = (double **) malloc((npha+3*nx*ny+1)*8);
	fprintf(stderr,"x,v,u= %d %d %d\n",x,v,u);
	xtemp = (double *) malloc((3*nx*ny+nstat+nevn)*8);
	vtemp = (double *) malloc((3*nx*ny+nstat+nevn)*8);
	utemp = (double *) malloc((npha+3*nx*ny+1)*8);
	fprintf(stderr,"xtemp,vtemp,utemp= %d %d %d\n",xtemp,vtemp,utemp);

	/* zero everything here */
	for(i=0;i<nstat;i++) station[i].delay=station[i].wei=station[i].sum=station[i].sum2=0.0;
	for(i=0;i<nevn;i++) event[i].delay=event[i].wei=event[i].sum=event[i].sum2=0.0;
	for(i=0;i<nx*ny;i++) {
		sol[i].ds=sol[i].cos=sol[i].sin=0.0;
		sol[i].wei=sol[i].weicos=sol[i].weisin=0.0;
		sol[i].sum=sol[i].sumcos=sol[i].sumsin=0.0;
		sol[i].sum2=sol[i].sum2cos=sol[i].sum2sin=sol[i].sum2cossin=0.0;
		sol[i].nu=0;
	}
	for(i=0;i<nstat+nevn+3*nx*ny;i++) xtemp[i]=vtemp[i]=0.0;
	for(i=0;i<npha+3*nx*ny;i++) utemp[i]=0.0;

	/*get preconditioners for the conjugate gradient*/
	/*assume you already know station & event counts*/
	for(i=0;i<npha;i++) {
		if(phase[i].quality<=0.0) continue;
		iin++;
		traceit(phase[i].slat-ymin, phase[i].slon-xmin, phase[i].rlat-ymin, phase[i].rlon-xmin,
			phase[i].scale, xcell,ycell,lencell,&ncell,&length);
		for(j=0;j<ncell;j++) {
			jj=xcell[j]+ycell[j]*nx;
			sol[jj].wei += lencell[j]*lencell[j];
			sol[jj].weicos += lencell[j]*cos(2.0*PIO180*phase[i].azse)*lencell[j]*cos(2.0*PIO180*phase[i].azse);
			sol[jj].weisin += lencell[j]*sin(2.0*PIO180*phase[i].azse)*lencell[j]*sin(2.0*PIO180*phase[i].azse);
			sol[jj].nu++;
		}
	}

	/*get estimate of rank via trace(ata/diag(ata+dtd) = anorm*/
	/*the 20 is from 4**2+1+1+1+1 for the lapacian damping*/
	anorm=0.0;
	for(i=0;i<nstat;i++) if(station[i].count!=0) anorm+=1;
	for(i=0;i<nevn;i++) if(event[i].count!=0) anorm+=1;
	for(i=0;i<nx*ny;i++) if(sol[i].nu!=0) anorm += sol[i].wei/(sol[i].wei+weight*weight*20);
#if ANIWEI
	for(i=0;i<nx*ny;i++) if(sol[i].nu!=0) {
		anorm += sol[i].weicos/(sol[i].weicos+aniwei*aniwei*20);
		anorm += sol[i].weisin/(sol[i].weisin+aniwei*aniwei*20);
	}
#endif
	fprintf(stderr,"\nctomo: anorm= %f\n\n",anorm);


#if ANIWEI
	for(i=0;i<nx*ny;i++) sol[i].weicos = sqrt(sol[i].weicos+20*aniwei*aniwei);
	for(i=0;i<nx*ny;i++) sol[i].weisin = sqrt(sol[i].weisin+20*aniwei*aniwei);
#endif
	for(i=0;i<nx*ny;i++) sol[i].wei = sqrt(sol[i].wei+20*weight*weight);
	for(i=0;i<nstat;i++) station[i].wei = sqrt((float)station[i].count);
	for(i=0;i<nevn;i++) event[i].wei = sqrt((float)event[i].count);
	fprintf(stderr,"weights found\n");

	/*figure out how many parameters*/
	for(i=0,ii=0;i<nstat;i++) if(station[i].count!=0) ii++;
	nstat2=ii;
	for(i=0,ii=0;i<nevn;i++) if(event[i].count!=0) ii++;
	nevn2=ii;
	for(i=0,ii=0;i<nx*ny;i++) if(sol[i].nu!=0) ii++;
	nsol2=ii;
	fprintf(stderr,"total stations,events,sols = %d %d %d\n",nstat2,nevn2,nsol2);

	/*assign spots to x pointer*/

	/*stations*/
	for(i=0,ii=0;i<nstat;i++) {
		if(station[i].count==0) continue;
		x[ii]= &xtemp[i];
		v[ii]= &vtemp[i];
		ii++;
	}
	/*event*/
	for(i=0,ii=ii;i<nevn;i++) {
		if(event[i].count==0) continue;
		x[ii]= &xtemp[i+nstat];
		v[ii]= &vtemp[i+nstat];

		ii++;
	}
	/*slowness*/
	for(i=0,ii=ii;i<nx*ny;i++) {
		if(sol[i].nu==0) continue;
		x[ii]= &xtemp[i+nstat+nevn];
		v[ii]= &vtemp[i+nstat+nevn];

		ii++;
	}
#if ANIWEI

	/*cos*/
	for(i=0,ii=ii;i<nx*ny;i++) {
		if(sol[i].nu==0) continue;

		x[ii]=&xtemp[i+nstat+nevn+nx*ny];
		v[ii]=&vtemp[i+nstat+nevn+nx*ny];
		ii++;
	}
	/*sin*/
	for(i=0,ii=ii;i<nx*ny;i++) {
		if(sol[i].nu==0) continue;
		x[ii]=&xtemp[i+nstat+nevn+2*nx*ny];
		v[ii]=&vtemp[i+nstat+nevn+2*nx*ny];
		ii++;
	}
#endif
	nparams=ii;

	/*put dtimes in u pointer*/
	for(i=0,ii=0;i<npha;i++) {
		if(phase[i].quality<=0.0) continue;
		u[ii] = &utemp[i];
		utemp[i] = phase[i].dtime;

		ii++;
	}
	ntimes=ii;
	fprintf(stderr,"number of good times = %d\n",ii);
	/*put in slowness damp eqns*/
	for(ii=ii,i=0;i<nx*ny;i++) {
		if(sol[i].nu==0) continue;
		u[ii] = &utemp[i+npha];
		utemp[i+npha] = 0.0;
		ii++;
	}
#if ANIWEI

	/*put in cos damp eqns*/
	for(ii=ii,i=0;i<nx*ny;i++) {
		if(sol[i].nu==0) continue;
		u[ii] = &utemp[i+npha+nx*ny];
		utemp[i+npha+nx*ny] = 0.0;
		ii++;
	}
	/*put in sin damp eqns*/
	for(ii=ii,i=0;i<nx*ny;i++) {
		if(sol[i].nu==0) continue;
		u[ii] = &utemp[i+npha+2*nx*ny];

		utemp[i+npha+2*nx*ny] = 0.0;
		ii++;
	}
#endif

	/* now add an eqn to zero station delays*/
#if ANIWEI
	u[ii] = &utemp[1+npha+3*nx*ny];
	utemp[1+npha+3*nx*ny] = 0.0;
#else
	u[ii] = &utemp[1+npha+nx*ny];
	utemp[1+npha+nx*ny] = 0.0;
#endif
	ii++;

	ntimes=ii;

	/*now do LSQR*/
	fprintf(stderr,"starting lsqr,ntimes,nparams= %d %d\n",ntimes,nparams);
	/******************************************************/
	lsqr(ntimes, nparams, x, u, v, nlsqr);

	/******************************************************/

	/*undo the preconditioning*/
	for(i=0;i<nstat;i++) if(station[i].count!=0) station[i].delay = xtemp[i]/station[i].wei;
	for(i=0;i<nevn;i++) if(event[i].count!=0) event[i].delay = xtemp[i+nstat]/event[i].wei;
	for(i=0;i<nx*ny;i++) if(sol[i].nu!=0) sol[i].ds = xtemp[i+nstat+nevn]/sol[i].wei;
#if ANIWEI

	for(i=0;i<nx*ny;i++) if(sol[i].nu!=0) sol[i].cos = xtemp[i+nstat+nevn+nx*ny]/sol[i].weicos;
	for(i=0;i<nx*ny;i++) if(sol[i].nu!=0) sol[i].sin = xtemp[i+nstat+nevn+2*nx*ny]/sol[i].weisin;
#endif


	/* update travel time for slowness & compute residuals */
	fprintf(stderr,"start residual computations\n");
	iin=res1=res2=res3=res4 = 0.0;
	for(i=0,iin=0; i<npha; i++) {
		if(phase[i].quality<=0.0) continue;
		iin++;

		traceit(phase[i].slat-ymin, phase[i].slon-xmin, phase[i].rlat-ymin, phase[i].rlon-xmin,
			phase[i].scale, xcell,ycell,lencell,&ncell,&length);
		time=phase[i].dtime;

		for(j=0;j<ncell;j++) {
			jj=xcell[j]+ycell[j]*nx;
			time -= lencell[j]*sol[jj].ds;
#if ANIWEI
			time -= lencell[j]*sol[jj].cos*cos(2.0*PIO180*phase[i].azse);
			time -= lencell[j]*sol[jj].sin*sin(2.0*PIO180*phase[i].azse);
#endif

		}

		time -= station[phase[i].stnno].delay;
		time -= event[phase[i].evnno].delay;
		phase[i].dtime=time;
		wei = phase[i].quality;
		res1 += fabs(time) * wei;
		res2 += time * time * wei;

		res3 += time * wei;
		res4 += wei;
	}

	fprintf(stderr,"\nres2 w/o damping residuals= %f\n",res2);
	fprintf(stderr,"anorm,log(res2/res4),penalty = %f, %f, %f\n",anorm,log(res2/res4),2*(anorm+1)/(res4-anorm-2));
	fprintf(stderr,"weight,aniwei= %f, %f\n",weight,aniwei);
	fprintf(stderr,"AIC= %f\n", log(res2/res4)+2*(anorm+1)/res4 );
	fprintf(stderr,"AICc= %f\n", log(res2/res4)+2*(anorm+1)/(res4-anorm-2) );
	fprintf(stderr,"FPE= %f\n", (res2/res4)*(res4+anorm)/(res4-anorm) );
	fprintf(stderr,"SIC= %f\n", log(res2/res4)+log(res4)*anorm/res4 );
	fprintf(stderr,"HQ= %f\n", log(res2/res4)+log(log(res4))*anorm/res4 );
	fprintf(stderr,"Cp= %f\n", res2/(.7065)+res4-2*anorm );
	fprintf(stderr,"AICth= %f\n", log(res2/res4)+6*(anorm+1)/res4 );
	fprintf(stderr,"AICs= %f\n\n", log(res2/res4)+(ntimes/nparams)*2*(anorm+1)/res4 );

	/*update res2 to account for the damping equations*/
	res2slo=res2cos=res2sin=0.0;
	for(j=0;j<ny;j++) for(i=0;i<nx;i++) {
		if(sol[i+nx*j].nu==0) continue;
		/*slo terms*/
		temp = -4.0*weight*sol[i+nx*j].ds;
		if(i-1>=0 && sol[i+nx*j-1].nu!=0)
			temp += weight*sol[i+nx*j-1].ds;
		if(i+1<nx && sol[i+nx*j+1].nu!=0)
			temp += weight*sol[i+nx*j+1].ds;
		if(j-1>=0 && sol[i+nx*(j-1)].nu!=0)
			temp += weight*sol[i+nx*(j-1)].ds;

		if(j+1<ny && sol[i+nx*(j+1)].nu!=0)
			temp += weight*sol[i+nx*(j+1)].ds;
		res2slo += temp*temp/(weight*weight);
		res2 += temp*temp;
#if ANIWEI
		/*cos terms*/
		temp = -4.0*aniwei*sol[i+nx*j].cos;
		if(i-1>=0 && sol[i+nx*j-1].nu!=0)
			temp += aniwei*sol[i+nx*j-1].cos;
		if(i+1<nx && sol[i+nx*j+1].nu!=0)
			temp += aniwei*sol[i+nx*j+1].cos;
		if(j-1>=0 && sol[i+nx*(j-1)].nu!=0)
			temp += aniwei*sol[i+nx*(j-1)].cos;
		if(j+1<ny && sol[i+nx*(j+1)].nu!=0)
			temp += aniwei*sol[i+nx*(j+1)].cos;
		res2cos += temp*temp/(aniwei*aniwei);
		res2 += temp*temp;
		/*sin terms*/

		temp = -4.0*aniwei*sol[i+nx*j].sin;
		if(i-1>=0 && sol[i+nx*j-1].nu!=0)
			temp += aniwei*sol[i+nx*j-1].sin;
		if(i+1<nx && sol[i+nx*j+1].nu!=0)
			temp += aniwei*sol[i+nx*j+1].sin;
		if(j-1>=0 && sol[i+nx*(j-1)].nu!=0)
			temp += aniwei*sol[i+nx*(j-1)].sin;
		if(j+1<ny && sol[i+nx*(j+1)].nu!=0)
			temp += aniwei*sol[i+nx*(j+1)].sin;
		res2sin += temp*temp/(aniwei*aniwei);
		res2 += temp*temp;
#endif
	}
	fprintf(stderr,"res2 with damping residuals= %f\n",res2);

	fprintf(stderr,"weight,aniwei= %f %f\n",weight,aniwei);

	/*print it out*/
	/*res1 & res2 are weighted residuals*/
	fprintf(stderr,"iin,res1,res2, %10d %10.6f %10.6f\n", iin, res1, res2);
	fprintf(stderr,"res3,res4 = %10.6f %10.5f\n", res3, res4);
	fprintf(stderr,"ntimes,nparams= %d, %d\n",ntimes,nparams);
	fprintf(stderr,"rms = %10.5f\n",sqrt(res2/(ntimes-nparams)));
	fprintf(stderr,"average delay= %20.15f\n", (res3 / res4));


#if !NBOOTV
	/*this uses the lrsq estimates of variance*/
	/*the sin-cos covariances are not estimated*/

	/*undo the preconditioning, and assign*/
	for(i=0;i<nstat;i++) if(station[i].count!=0) station[i].sum2 = vtemp[i]/station[i].wei/station[i].wei;
	for(i=0;i<nevn;i++) if(event[i].count!=0) event[i].sum2 = vtemp[i+nstat]/event[i].wei/event[i].wei;
	for(i=0;i<nx*ny;i++) if(sol[i].nu!=0) sol[i].sum2 = vtemp[i+nstat+nevn]/sol[i].wei/sol[i].wei;
#if ANIWEI
	for(i=0;i<nx*ny;i++) if(sol[i].nu!=0) sol[i].sum2cos = vtemp[i+nstat+nevn+nx*ny]/sol[i].weicos/sol[i].weicos;
	for(i=0;i<nx*ny;i++) if(sol[i].nu!=0) sol[i].sum2sin = vtemp[i+nstat+nevn+2*nx*ny]/sol[i].weisin/sol[i].weisin;
#endif
#endif

#if NBOOTV
	/*do bootstrap errors, put results into sum & sum2*/
	fprintf(stderr,"\nstart bootstrap resolution\n");
	/*put into sum2 and sin*/
	for(i=0;i<nstat;i++) station[i].sum=station[i].sum2=0.0;

	for(i=0;i<nevn;i++) event[i].sum=event[i].sum2=0.0;
	for(i=0;i<nx*ny;i++) sol[i].sum=sol[i].sum2=0.0;

	/*store a copy of phase in phase2*/
	for(i=0;i<npha;i++) phase2[i]=phase[i];

	/*do the bootstrap iterations*/
	for(k=0;k<nbootv;k++) {
		fprintf(stderr,"bootstrap errors k= %d\n",k);


		/*fill up the traveltimes in utemp*/
		for(i=0;i<npha;i++) {
			if(phase2[i].quality<=0) continue;
			while(phase2[l=npha*ran1(4321)].quality<=0);
			phase[i] = phase2[l];
			utemp[i]=phase[i].ttime;
		}
		/*put in damp eqns*/
		for(i=0;i<nx*ny;i++) {
			if(sol[i].nu==0) continue;
			utemp[i+npha]=0.0;
#if ANIWEI
			utemp[i+npha+nx*ny]=utemp[i+npha+2*nx*ny]=0.0;
#endif
		}

		/*now do LSQR*/
		fprintf(stderr,"starting lsqr,ntimes,nparams,u[0]= %d %d %f\n",ntimes,nparams,*u[0]);
		lsqr(ntimes, nparams, x, u, v, nlsqr);

		/*undo the preconditioning*/
		for(j=0;j<nstat;j++) if(station[j].count!=0) xtemp[j] = xtemp[j]/station[j].wei;
		for(j=0;j<nevn;j++) if(event[j].count!=0) xtemp[j+nstat] = xtemp[j+nstat]/event[j].wei;
		for(j=0;j<nx*ny;j++) if(sol[j].nu!=0) xtemp[j+nstat+nevn] = xtemp[j+nstat+nevn]/sol[j].wei;
#if ANIWEI
		for(i=0;i<nx*ny;i++) if(sol[i].nu!=0) xtemp[i+nstat+nevn+nx*ny] /= sol[i].weicos;
		for(i=0;i<nx*ny;i++) if(sol[i].nu!=0) xtemp[i+nstat+nevn+2*nx*ny] /= sol[i].weisin;
#endif

		/*update the variance in sum*/
		/*stations*/
		for(j=0;j<nstat;j++) if(station[j].count!=0) {

			station[j].sum2 += (xtemp[j]*xtemp[j]);
			station[j].sum += xtemp[j];
		}
		/*events*/
		for(j=0;j<nevn;j++) if(event[j].count!=0) {
			event[j].sum2 += (xtemp[j+nstat]*xtemp[j+nstat]);
			event[j].sum += xtemp[j+nstat];
		}
		/*slownesses*/
		for(j=0;j<nx*ny;j++) if(sol[j].nu!=0) {
			sol[j].sum2 += (xtemp[j+nstat+nevn]*xtemp[j+nstat+nevn]);
			sol[j].sum += xtemp[j+nstat+nevn];
		}
#if ANIWEI
		/*cosines*/

		for(j=0;j<nx*ny;j++) if(sol[j].nu!=0) {
			sol[j].sum2cos += xtemp[j+nstat+nevn+nx*ny]*xtemp[j+nstat+nevn+nx*ny];
			sol[j].sumcos += xtemp[j+nstat+nevn+nx*ny];
		}
		/*sines*/
		for(j=0;j<nx*ny;j++) if(sol[j].nu!=0) {
			sol[j].sum2sin += xtemp[j+nstat+nevn+2*nx*ny]*xtemp[j+nstat+nevn+2*nx*ny];
			sol[j].sumsin += xtemp[j+nstat+nevn+2*nx*ny];
		}
		/*cos-sin*/
		for(j=0;j<nx*ny;j++) if(sol[j].nu!=0) {
			sol[j].sum2cossin += xtemp[j+nstat+nevn+nx*ny]*xtemp[j+nstat+nevn+2*nx*ny];
		}
#endif
	}

	/*compute mean in sum & variance in sum2*/
	/*stations*/
	for(j=0;j<nstat;j++) {
		if(station[j].count==0) continue;
		station[j].sum /= nbootv; /*bootstrap estimate ave*/
		station[j].sum2 = (station[j].sum2 - nbootv*station[j].sum*station[j].sum)/(nbootv-1);
	}
	/*events*/
	for(j=0;j<nevn;j++) {
		if(event[j].count==0) continue;
		event[j].sum /= nbootv;
		event[j].sum2 = (event[j].sum2 - nbootv*event[j].sum*event[j].sum)/(nbootv-1);
	}
	/*slowness*/
	for(j=0;j<nx*ny;j++) {
		if(sol[j].nu==0) continue;
		sol[j].sum /= nbootv;
		sol[j].sum2 = (sol[j].sum2 - nbootv*sol[j].sum*sol[j].sum)/(nbootv-1);;
	}
#if ANIWEI
	/*cosines*/
	for(j=0;j<nx*ny;j++) {
		if(sol[j].nu==0) continue;
		sol[j].sumcos /= nbootv;
		sol[j].sum2cos = (sol[j].sum2cos - nbootv*sol[j].sumcos*sol[j].sumcos)/(nbootv-1);;
	}
	/*sines*/

	for(j=0;j<nx*ny;j++) {
		if(sol[j].nu==0) continue;
		sol[j].sumsin /= nbootv;
		sol[j].sum2sin = (sol[j].sum2sin - nbootv*sol[j].sumsin*sol[j].sumsin)/(nbootv-1);;
	}
	/*cos-sin*/
	for(j=0;j<nx*ny;j++) {
		if(sol[j].nu==0) continue;
		sol[j].sum2cossin = (sol[j].sum2cossin - nbootv*sol[j].sumcos*sol[j].sumsin)/(nbootv-1);;
	}
#endif
	fprintf(stderr,"sum of variances = %g\n",time);
#endif
/*end of NBOOTV if*/

	/* write it all out now */
	fprintf(stderr,"start writing it out\n");
	i = write(fdsol, sol, nx * ny * sizeof(struct solution));
	fprintf(stderr,"sol write= %d\n", i);

	/* all done */
	fprintf(stderr,"ctomo all done now\n\n");
}

/*******************************************************************/
/*****************************************************************/
/*these are the traceing subroutines*/


traceit(slat, slon, rlat, rlon, scale, xcell,ycell,lencell,ncell,length)
/* computes the travel time from the locations */
/* assumes straight lines */
/*xcell,ycell,are matrixes that contain the cell number*/
/*ncell is the total number of cells the ray crosses*/
float	slat, slon, rlat, rlon, scale;
int *ncell;
int *xcell,*ycell;
double *lencell, *length;
{
	double	x1, y1, x2, y2, dx, dy;
	double	len, x, y;
	int		ix, iy;
	double	xnew, ynew;
	int		ixnew, iynew;
	double	alpha, sgnalp;
	double	amax();

	/*zero some stuff*/
	*ncell=0;

	/* set rlon to be minimum, if not switch */
	if(slon<rlon) {
		x1 = x = slon;
		x2 = rlon;
		y1 = y = slat;
		y2 = rlat;
	} else {
		x1 = x = rlon;
		x2 = slon;
		y1 = y = rlat;
		y2 = slat;
	}

	/* set some parameters */
	dx = x2 - x1;
	dy = y2 - y1;
	if(dy == 0.0) dy = 0.0001;
	*length = sqrt(dx*dx + dy*dy) * scale;
	ix = intt(x1 * xmesh);
	iy = intt(y1 * ymesh);
	alpha = dx / dy;
	sgnalp = sign(alpha);

	/* now trace */
	while (x<x2) {
		ynew = (iy + amax(sgnalp, 0.0)) / ymesh;
		xnew = x + alpha * (ynew - y);
		if(xnew>x2) {
			xnew = x2;
			ynew = y2;
		}
		ixnew = intt(xnew * xmesh);
		iynew = iy + sgnalp;
		if(ixnew>ix) {
			ixnew = ix + 1;
			xnew = ((float) ixnew) / xmesh;
			ynew = y + (xnew - x) / alpha;
			iynew = iy;
		}
		dx = xnew - x;
		dy = ynew - y;
		len = sqrt(dx*dx + dy*dy) * scale;

		/* this is the important loop */
		if(ix >= 0 && ix<nx && iy >= 0 && iy<ny) {
			xcell[*ncell]=ix;
			ycell[*ncell]=iy;
			lencell[*ncell]=len;
			(*ncell)++;
		}

		ix = ixnew;
		iy = iynew;
		x = xnew;
		y = ynew;
	}
/* urbi: compiler error
	/* return; */
}

/*********************lsqr routines *********************************/
lsqr(m, n, x, u, v,itmax)
/* subroutine to solve the linear tomographic problem Ax=u using */
/* the lsqr algorithm from G.Nolet, Seismic Tomography, p18 */
/*this version has more error control stuff in it*/
/* m is number of data */
/* n is number of unknown */
/* x(n) is the solution */
/* u(m) is the data (overwritten), passed in as phase*/
/* scratch array v(n) is zeroed and then overwritten by the variances*/
/* the array w(n) & v(n) is dynamically allocated */
/* avpu(m,n,u,v) computes u=u+A*v, ie the forward problem */
/* atupv(m,n,u,v) computes v=v+At*u, ie the backprojection*/
int			m, n;
int			itmax;
double		**x,**u,**v;
/*note that x and v are pointers to pointers*/
{
	double damp,anorm,acond,rnorm,arnorm,xnorm;
	double	*w;
	double *se;

	int		i, j;

	double sqrt(), fabs(), normlz();

	double alfa,bbnorm,beta;
	double cs,cs1,cs2,dampsq,ddnorm,delta;
	double gamma,gambar,phi,phibar,psi;
	double res1,res2,rho,rhobar,rhbar1,rhbar2,rhs;
	double sn,sn1,sn2,t,tau,test2;
	double theta,t1,t2,t3,xxnorm,z,zbar;

	fprintf(stderr,"lsqr: m,n= %d %d\n",m,n);

	/* define data, soln, & scratch areas*/
	w = (double *) malloc(8*n);
	se = (double *) malloc(8*n);

	/*initialize everything that needs it*/
	damp=anorm=acond=bbnorm=ddnorm=res2=xnorm=xxnorm=sn2=z=0.0;
	dampsq=damp*damp;
	cs2= -1.0;
	for(i=0; i<n; i++) *v[i]=*x[i]=se[i]=0.0;

	/*start of lsqr*/
	/*set up first bidiagonalizition vectors*/
	/*these satisfy beta*u=b, alfa*v=atrans*u */
	beta = normlz(m, u);
	atupv(m,n,u,v);
	alfa = normlz(n, v);
	for(i=0; i<n; i++) w[i] = *v[i];
	rhobar = alfa;
	phibar = beta;
	rnorm = beta;
	arnorm = alfa*beta;
	test2=alfa/beta;

	fprintf(stderr,"iter,arnorm,rnorm,anorm,acond,xnorm,test2 S %f %f %f %f %f %f\n",
		arnorm,rnorm,anorm,acond,xnorm,test2);

	/*the main iteration loop*/
	for (i = 0; i < itmax; i++) {

		/*do the next bidiagonalization for beta,u,alfa,v updates*/
		/*then beta*u=a*v-alfa*u and alfa*v=At*u-beta*v*/
		for (j=0;j<m;j++) *u[j] *= -alfa;	/* bidiagonalization */
		avpu(m,n,u,v);
		beta = normlz(m, u);
		bbnorm += alfa*alfa + beta*beta + dampsq;
		for(j=0;j<n;j++) *v[j] *= -beta;
		atupv(m,n,u,v);
		alfa = normlz(n, v);

		/*use plate rotation to eliminate damping parameter*/
		rhbar2 = rhobar*rhobar + dampsq;
		rhbar1 = sqrt(rhbar2);
		cs1 = rhobar/rhbar1;
		sn1 = damp/rhbar1;
		psi= sn1*phibar;
		phibar *= cs1;

		/*use plane rotation to elimniate the subdiagonal element (beta) */
		rho = sqrt(rhbar2 + beta * beta);	/* modified QR */
		cs = rhbar1 / rho;
		sn = beta / rho;
		theta = sn * alfa;
		rhobar = -cs * alfa;
		phi = cs * phibar;
		phibar = sn * phibar;	/* phibar is the sqrt of the sum-of-squares */
		tau = sn*phi;

		/*update x,w and errors*/
		t1 = phi / rho;
		t2 = -theta / rho;
		t3 = 1.0/rho;
		for (j = 0; j < n; j++) {
			t= w[j];
			*x[j] = t1 * t + (*x[j]);
			w[j] = t2 * t + (*v[j]);
			t = t3*t*t3*t;
			se[j] += t;
			ddnorm += t;
		}

		/*plane rotate on right to eleminate theta*/
		/*and use to estimate xnorm */
		delta = sn2*rho;
		gambar = -cs2*rho;
		rhs = phi-delta*z;
		zbar = rhs/gambar;
		xnorm = sqrt(xxnorm+zbar*zbar);
		gamma = sqrt(gambar*gambar+theta*theta);
		cs2 = gambar/gamma;
		sn2 = theta/gamma;
		z = rhs/gamma;
		xxnorm += z*z;

		/*estimates norms */
		anorm = sqrt(bbnorm);	/*anorm is estimate of frobenius norm*/
		acond = anorm*sqrt(ddnorm);
		res1 = phibar*phibar;
		res2 += psi*psi;
		rnorm = sqrt(res1+res2);	/*estimate of sqrt(norm)*/
		arnorm = alfa*fabs(tau); /*arnorm is Atranspose X, should be 0 if converged*/

		/*get test2 from norms (machine limit at convergence)*/
		test2=arnorm/(anorm*rnorm);

		fprintf(stderr,"iter,arnorm,rnorm,anorm,acond,xnorm,test2 %d %f %f %f %f %f %f\n",
			i,arnorm,rnorm,anorm,acond,xnorm,test2);

		/* arnorm is size of  Atrans*r */
		/* rnorm is sqrt of sos of residual vector r */
		/* anorm is norm of A */
		/* acond is condition number */
		/* xnorm is norm of sln vector x */
		/* test2 < machine limit at convergence */
	}

	for(i=0;i<n;i++) *v[i] = se[i]*rnorm*rnorm/(m-n);
/* nothing to return : compiler error:urbi
	/* return; */
}

/*************************************************************************/
double normlz(n, x)
int		n;
double	**x;

/* normalizes vector x */
{
	int i,j;
	double ss,s ;

	/*fprintf(stderr,"normlz: start, n= %d\n",n);*/

	for(i=0,ss=0.0; i<n; i++) ss += (*x[i])*(*x[i]);
	s= sqrt(ss);
	ss=1.0/s;
	for(i=0; i<n; i++) *x[i] *= ss;

	/*fprintf(stderr,"normlz: s, s**2= %f %f\n",s,s*s);*/
	return(s);
}

/*********************************************************************/

avpu(m,n,u,v)
int m,n; double **u,**v;
/*the forward projection*/
{
	int		i,j,ii,jj,k,kk,l,iii;
	int ncell;
	int xcell[999],ycell[999];
	double lencell[999];
	double length;
	double exp();

	/*fprintf(stderr,"avpu: start\n");*/

	/*for each ray call trace*/
	for(i=0;i<npha;i++) {
		if(phase[i].quality<=0.0) continue;
		traceit(phase[i].slat-ymin, phase[i].slon-xmin, phase[i].rlat-ymin, phase[i].rlon-xmin,
				phase[i].scale, xcell,ycell,lencell,&ncell,&length);
		for(j=0;j<ncell;j++) {
			jj=xcell[j]+ycell[j]*nx;
#if ANIWEI
			utemp[i] += lencell[j]*vtemp[jj+nstat+nevn+nx*ny]*cos(2.0*PIO180*phase[i].azse)/sol[jj].weicos
				+ lencell[j]*vtemp[jj+nstat+nevn+2*nx*ny]*sin(2.0*PIO180*phase[i].azse)/sol[jj].weisin
				+ lencell[j]*vtemp[jj+nstat+nevn]/sol[jj].wei;
#else
			utemp[i] += lencell[j]*vtemp[jj+nstat+nevn]/sol[jj].wei;
#endif
		}
		utemp[i] += vtemp[phase[i].stnno]/station[phase[i].stnno].wei
			+ vtemp[phase[i].evnno+nstat]/event[phase[i].evnno].wei;
	}

	/*zero station delays*/
	for(i=0;i<nstat;i++) {
		if(station[i].count==0) continue;
#if ANIWEI
		utemp[npha+1+3*nx*ny] += vtemp[i]/station[i].wei;
#else
		utemp[npha+1+nx*ny] += vtemp[i]/station[i].wei;
#endif
	}


	/*damp slowness here*/
	for(j=0;j<ny;j++) for(i=0;i<nx;i++) {
		if(sol[i+nx*j].nu==0) continue;
		utemp[npha+i+nx*j] += -4.0*weight*vtemp[i+nx*j+nstat+nevn]/sol[i+nx*j].wei;
		if(i-1>=0 && sol[i+nx*j-1].nu!=0) utemp[npha+i+nx*j] += weight*vtemp[i+nx*j-1+nstat+nevn]/sol[i+nx*j-1].wei;
		if(i+1<nx && sol[i+nx*j+1].nu!=0) utemp[npha+i+nx*j] += weight*vtemp[i+nx*j+1+nstat+nevn]/sol[i+nx*j+1].wei;
		if(j-1>=0 && sol[i+nx*(j-1)].nu!=0) utemp[npha+i+nx*j] += weight*vtemp[i+nx*(j-1)+nstat+nevn]/sol[i+nx*(j-1)].wei;
		if(j+1<ny && sol[i+nx*(j+1)].nu!=0) utemp[npha+i+nx*j] += weight*vtemp[i+nx*(j+1)+nstat+nevn]/sol[i+nx*(j+1)].wei;
	}
#if ANIWEI
	/*damp cos here*/
	for(j=0;j<ny;j++) for(i=0;i<nx;i++) {
		if(sol[i+nx*j].nu==0) continue;
		utemp[npha+i+nx*j+nx*ny] += -4.0*aniwei*vtemp[i+nx*j+nstat+nevn+nx*ny]/sol[i+nx*j].weicos;
		if(i-1>=0 && sol[i+nx*j-1].nu!=0)
			utemp[npha+i+nx*j+nx*ny] += aniwei*vtemp[i+nx*j-1+nstat+nevn+nx*ny]/sol[i+nx*j-1].weicos;
		if(i+1<nx && sol[i+nx*j+1].nu!=0)
			utemp[npha+i+nx*j+nx*ny] += aniwei*vtemp[i+nx*j+1+nstat+nevn+nx*ny]/sol[i+nx*j+1].weicos;
		if(j-1>=0 && sol[i+nx*(j-1)].nu!=0)
			utemp[npha+i+nx*j+nx*ny] += aniwei*vtemp[i+nx*(j-1)+nstat+nevn+nx*ny]/sol[i+nx*(j-1)].weicos;
		if(j+1<ny && sol[i+nx*(j+1)].nu!=0)
			utemp[npha+i+nx*j+nx*ny] += aniwei*vtemp[i+nx*(j+1)+nstat+nevn+nx*ny]/sol[i+nx*(j+1)].weicos;
	}
	/*damp sin here*/
	for(j=0;j<ny;j++) for(i=0;i<nx;i++) {
		if(sol[i+nx*j].nu==0) continue;
		utemp[npha+i+nx*j+2*nx*ny] += -4.0*aniwei*vtemp[i+nx*j+nstat+nevn+2*nx*ny]/sol[i+nx*j].weisin;
		if(i-1>=0 && sol[i+nx*j-1].nu!=0)
			utemp[npha+i+nx*j+2*nx*ny] += aniwei*vtemp[i+nx*j-1+nstat+nevn+2*nx*ny]/sol[i+nx*j-1].weisin;
		if(i+1<nx && sol[i+nx*j+1].nu!=0)
			utemp[npha+i+nx*j+2*nx*ny] += aniwei*vtemp[i+nx*j+1+nstat+nevn+2*nx*ny]/sol[i+nx*j+1].weisin;
		if(j-1>=0 && sol[i+nx*(j-1)].nu!=0)
			utemp[npha+i+nx*j+2*nx*ny] += aniwei*vtemp[i+nx*(j-1)+nstat+nevn+2*nx*ny]/sol[i+nx*(j-1)].weisin;
		if(j+1<ny && sol[i+nx*(j+1)].nu!=0)
			utemp[npha+i+nx*j+2*nx*ny] += aniwei*vtemp[i+nx*(j+1)+nstat+nevn+2*nx*ny]/sol[i+nx*(j+1)].weisin;
	}
#endif

	/*fprintf(stderr,"avpu done\n");*/

}

atupv(m,n,u,v)
int m,n; double **u,**v;
/*backproject celldistance * time*/
/*u is data vector in phase[].empty*/
/*v is temp solution vector is ds2 & delay2*/
/*this is the gradient operator*/
{

	int	i,j,ii,jj,k,kk,l,iii;
	int ncell;
	int xcell[999],ycell[999];
	double lencell[999];
	double length;
	double exp();

	/*fprintf(stderr,"atupv: start\n");*/

	/*backproject the ray*/
	for(i=0,ii=0;i<npha;i++) {
		if(phase[i].quality<=0.0) continue;
		/*fprintf(stderr,"atupv: start traceit, i,ii= %d %d",i,ii);*/
		traceit(phase[i].slat-ymin, phase[i].slon-xmin, phase[i].rlat-ymin, phase[i].rlon-xmin,
			phase[i].scale, xcell,ycell,lencell,&ncell,&length);
		for(j=0;j<ncell;j++) {
			jj=xcell[j]+ycell[j]*nx;
			vtemp[jj+nstat+nevn] += lencell[j]*utemp[i]/sol[jj].wei;
#if ANIWEI
			vtemp[jj+nstat+nevn+nx*ny] += lencell[j]*utemp[i]*cos(2.0*PIO180*phase[i].azse)/sol[jj].weicos;
			vtemp[jj+nstat+nevn+2*nx*ny] += lencell[j]*utemp[i]*sin(2.0*PIO180*phase[i].azse)/sol[jj].weisin;
#endif
		}
		vtemp[phase[i].stnno] += utemp[i]/station[phase[i].stnno].wei;
		vtemp[phase[i].evnno+nstat] += utemp[i]/event[phase[i].evnno].wei;
		ii++;
	}

	/*zero station delays*/
	for(i=0;i<nstat;i++) {
		if(station[i].count==0) continue;
#if ANIWEI
		vtemp[i] += utemp[npha+3*nx*ny+1]/station[i].wei;
#else
		vtemp[i] += utemp[npha+nx*ny+1]/station[i].wei;
#endif
	}


	iii=ii;

 	/*damp slowness here*/
 	/*divide by weight*/
	for(ii=0,j=0;j<ny;j++) for(i=0;i<nx;i++) {
		if(sol[i+nx*j].nu==0) continue;
		vtemp[i+nx*j+nstat+nevn] += -4.0*weight*utemp[npha+i+nx*j]/sol[i+nx*j].wei;
		if(i-1>=0 && sol[i+nx*j-1].nu!=0) vtemp[i+nx*j-1+nstat+nevn] += weight*utemp[npha+i+nx*j]/sol[i+nx*j-1].wei;
		if(i+1<nx && sol[i+nx*j+1].nu!=0) vtemp[i+nx*j+1+nstat+nevn] += weight*utemp[npha+i+nx*j]/sol[i+nx*j+1].wei;
		if(j-1>=0 && sol[i+nx*(j-1)].nu!=0) vtemp[i+nx*(j-1)+nstat+nevn] += weight*utemp[npha+i+nx*j]/sol[i+nx*(j-1)].wei;
		if(j+1<ny && sol[i+nx*(j+1)].nu!=0) vtemp[i+nx*(j+1)+nstat+nevn] += weight*utemp[npha+i+nx*j]/sol[i+nx*(j+1)].wei;
		ii++;
	}
#if ANIWEI
	/*damp cos here*/
	for(ii=0,j=0;j<ny;j++) for(i=0;i<nx;i++) {
		if(sol[i+nx*j].nu==0) continue;
		vtemp[i+nx*j+nstat+nevn+nx*ny] += -4.0*aniwei*utemp[npha+i+nx*j+nx*ny]/sol[i+nx*j].weicos;
		if(i-1>=0 && sol[i+nx*j-1].nu!=0)
			vtemp[i+nx*j-1+nstat+nevn+nx*ny] += aniwei*utemp[npha+i+nx*j+nx*ny]/sol[i+nx*j-1].weicos;
		if(i+1<nx && sol[i+nx*j+1].nu!=0)
			vtemp[i+nx*j+1+nstat+nevn+nx*ny] += aniwei*utemp[npha+i+nx*j+nx*ny]/sol[i+nx*j+1].weicos;
		if(j-1>=0 && sol[i+nx*(j-1)].nu!=0)
			vtemp[i+nx*(j-1)+nstat+nevn+nx*ny] += aniwei*utemp[npha+i+nx*j+nx*ny]/sol[i+nx*(j-1)].weicos;
		if(j+1<ny && sol[i+nx*(j+1)].nu!=0)
			vtemp[i+nx*(j+1)+nstat+nevn+nx*ny] += aniwei*utemp[npha+i+nx*j+nx*ny]/sol[i+nx*(j+1)].weicos;
		ii++;
	}
	/*damp sin here*/
	for(ii=0,j=0;j<ny;j++) for(i=0;i<nx;i++) {
		if(sol[i+nx*j].nu==0) continue;
		vtemp[i+nx*j+nstat+nevn+2*nx*ny] += -4.0*aniwei*utemp[npha+i+nx*j+2*nx*ny]/sol[i+nx*j].weisin;
		if(i-1>=0 && sol[i+nx*j-1].nu!=0)
			vtemp[i+nx*j-1+nstat+nevn+2*nx*ny] += aniwei*utemp[npha+i+nx*j+2*nx*ny]/sol[i+nx*j-1].weisin;
		if(i+1<nx && sol[i+nx*j+1].nu!=0)
			vtemp[i+nx*j+1+nstat+nevn+2*nx*ny] += aniwei*utemp[npha+i+nx*j+2*nx*ny]/sol[i+nx*j+1].weisin;
		if(j-1>=0 && sol[i+nx*(j-1)].nu!=0)
			vtemp[i+nx*(j-1)+nstat+nevn+2*nx*ny] += aniwei*utemp[npha+i+nx*j+2*nx*ny]/sol[i+nx*(j-1)].weisin;
		if(j+1<ny && sol[i+nx*(j+1)].nu!=0)
			vtemp[i+nx*(j+1)+nstat+nevn+2*nx*ny] += aniwei*utemp[npha+i+nx*j+2*nx*ny]/sol[i+nx*(j+1)].weisin;
		ii++;
	}
#endif


	/*fprintf(stderr,"atupv done\n");*/
}

/*************************************************/
double ran1(idum)
int idum;
/*returns random number 0.0 inclusive to 1.0 exclusive*/
/*thearn 6/90*/
{
	static int seed=0;
	float x;

	if(seed==0) {
		srand(idum);
		seed=1;
	}
	x = (float)rand()/(2147483647);
	/*max on sun is 2**31-1*/
	/*max under gcc is 32767 */
	/*max under linux is 2147483647*/
	/*fprintf(stderr,"ran1:x,seed,idum %f %d %d\n",x,seed,idum);*/
	return x;
}

/***************************************************************/
