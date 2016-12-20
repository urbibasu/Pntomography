/* syphase.c - to use for Pn data */
/* this creates test data */
/* this version does impulse response */
/* cellmesh is cells/square */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fcntl.h>
/*next to lines are for mmap()*/
#include <sys/types.h>
#include <sys/mman.h>
/*next two for lseek()*/
#include <sys/types.h>
#include <unistd.h>

#define PIO180	0.0175329252
#define PI		3.141592654
/*slowness levels: .00390625=>0.25k/s */
#define ANI		0.00390625
#define SLO		0.00390625
/*#define SLO		0.008*/
#define STA		0.0


float xmin=0.0,xmax=360.0,ymin=-90.0,ymax=90.0;
float xmesh=5.0,ymesh=5.0;
int nx, ny;

#include "data14.h"

double sqrt(),fabs(),cos(),sin(),atan2(),log(),exp();
double gasdev();
int npha, nevn, nstat, nsol;

struct rec *phase;
struct locations *event,*station;
struct solution *sol;

main(argc,argv)
int argc; char **argv;
{
	int i,j,k,jj;
	int knt=0;
	double res1=0.0, res2=0.0, res3=0.0, res4=0.0;
	float cellmesh=0.0, animesh=0.0, stamesh=0.0;
	int fdrec, fdevn, fdstn, fdsol;
	FILE *fopen();
	double	inter,slope;
	int xcell[999],ycell[999];
	double lencell[999],length;
	int ncell;
	
	fprintf(stderr,"ANI,SLO,STA= %f %f %f\n",ANI,SLO,STA);

	if(getarg(argc,argv,"-region","%f/%f/%f/%f",&xmin,&xmax,&ymin,&ymax)<4) errprint("readsol: cant getarg region");
	if(getarg(argc,argv,"-mesh","%f/%f",&xmesh,&ymesh)<2) errprint("readsol: cant getarg xmesh");
	getarg(argc,argv,"-cellmesh","%f",&cellmesh);
	getarg(argc,argv,"-animesh","%f",&animesh);
	getarg(argc,argv,"-stamesh","%f",&stamesh);
	fprintf(stderr,"xmin,xmax,ymin,ymax,xmesh,ymesh %f %f %f %f %f %f\n",xmin,xmax,ymin,ymax,xmesh,ymesh);
	fprintf(stderr,"cellmesh,animesh,stamesh= %f %f %f\n",cellmesh,animesh,stamesh);
	nx= (xmax-xmin)*xmesh;
	ny= (ymax-ymin)*ymesh;

	printf("starting syphase\n");
	
	/* open data files */
	fdrec=open("../data/phases.b",O_RDWR);
	fdevn=open("../data/events.b",O_RDWR);
	fdstn=open("../data/stations.b",O_RDWR);
	fdsol=open("../data/sol.b",O_RDWR);
	printf("fdrec= %d %d %d %d\n",fdrec,fdevn,fdsol,fdstn);

	/* set the data pointers */
	sol = (struct solution *) mmap(0,lseek(fdsol,0L,SEEK_END),PROT_READ|PROT_WRITE,MAP_SHARED,fdsol,0);
	nsol = lseek(fdsol,0L,SEEK_END)/sizeof(struct solution);
	phase = (struct rec *) mmap(0,lseek(fdrec,0L,SEEK_END),PROT_READ|PROT_WRITE, MAP_SHARED,fdrec,0);
	npha = lseek(fdrec,0L,SEEK_END)/sizeof(struct rec);
	station = (struct locations *) mmap(0,lseek(fdstn,0L,SEEK_END),PROT_READ|PROT_WRITE, MAP_SHARED,fdstn,0);
	nstat = lseek(fdstn,0L,SEEK_END)/sizeof(struct locations);
	event = (struct locations *) mmap(0,lseek(fdevn,0L,SEEK_END),PROT_READ|PROT_WRITE, MAP_SHARED,fdevn,0);
	nevn = lseek(fdevn,0L,SEEK_END)/sizeof(struct locations);
	printf("station,event,phase,sol= %d %d %d %d\n",station,event,phase,sol);
	printf("nstat,nevn,nphan,sol = %d %d %d %d\n", nstat, nevn, npha, nsol);
	printf("nx,ny= %d %d\n",nx,ny);

	/* fix sol by zeroing*/

	for(i=0;i<nx*ny;i++) sol[i].sum=sol[i].cos=sol[i].sin=sol[i].ds=sol[i].sumcos=sol[i].sumsin=0.0;

	/* fix stations by zeroing*/
	for(i=0; i<nstat; i++) station[i].delay=station[i].sum=0.0;

	/*fix events with zeroing*/
	for(i=0; i<nevn; i++) event[i].delay=event[i].sum=0.0;

	/***********get the generating model*********/
	/*comment out to get previous model*/
	/*getmod  - puts spike in center */
	/*getmod2 - chekerboard */
	/*getmod3 - sinusoids */
	/*getmod4 - stripes */
	/*getmod5 - Dinaride velocity anomaly test */
	/*getmod6 - dinaride anisotropy anomoly test */
	/*getmod7 - constant anisotropy*/
	/*getmod8 - constant velocity*/
	/*getmod9 - anisotropic rectangle*/
	printf("starting syndat for synthetic data\n");
	/*getmod(cellmesh,stamesh,animesh); */
	getmod2(cellmesh,stamesh,animesh);
	/********************************************/

	/*now fix records*/
	for(i=0; i<npha; i++) {
		if(i%1000==0) printf("%d records\r",i);
		if(phase[i].quality <= 0.0) continue;
		knt++;
		phase[i].quality=1.0;
		phase[i].delay=0.0;
		phase[i].empty=0.0;
		phase[i].ttime =phase[i].dtime=0.0;
		/*add NOISE here*/
		phase[i].ttime=gasdev()*1.4;
		/*do the forward trace*/
		traceit(phase[i].slat-ymin,phase[i].slon-xmin, phase[i].rlat-ymin, phase[i].rlon-xmin,
			phase[i].scale, xcell,ycell,lencell,&ncell,&length);
		for(j=0;j<ncell;j++) {
			jj=xcell[j]+ycell[j]*nx;
			phase[i].ttime += lencell[j]*sol[jj].ds;
			phase[i].ttime += lencell[j]*sol[jj].cos*cos(2.0*PIO180*phase[i].azse);
			phase[i].ttime += lencell[j]*sol[jj].sin*sin(2.0*PIO180*phase[i].azse);
		}
		phase[i].ttime += station[phase[i].stnno].delay+event[phase[i].evnno].delay;

		/*calculate norms*/
		res1 += fabs(phase[i].ttime)*phase[i].quality;
		res2 += phase[i].ttime*phase[i].ttime*phase[i].quality;
		res3 += phase[i].ttime*phase[i].quality;
		res4 += phase[i].quality;
	}
	/*reset dtime to be ttime*/
	for(i=0;i<npha;i++) phase[i].dtime = phase[i].ttime;
	
	printf("phases fixed\n");
	printf("knt,res1,res2,res3,res4= %d %e %e %e %f\n",
		knt,res1,res2,res3,res4);

	/* write it all out now*/
	printf("done- knt= %8d\n",knt);
	printf("res1, res2 = %10.5f %14.9f\n",res1,res2);
	printf("ave, res4 = %10.5f %10.5f\n\n",(res3/res4),res4);
}

/*************************************************************************/
/*************************************************************************/

getmod(cellmesh,stamesh,animesh)
double cellmesh,stamesh,animesh;
/*this one just puts a spike in the middle*/
{
	int i,j;
	
	if(cellmesh!=0.0)
		for(i=-xmesh/cellmesh/2; i<=xmesh/cellmesh/2; i++)
			for(j=-ymesh/cellmesh/2; j<=ymesh/cellmesh/2; j++)
				if(i*i+j*j <= xmesh*ymesh/cellmesh/cellmesh/4.0) {
					if(sol[i+nx/3+nx*(j+ny/3)].nu>0) sol[i+nx/3+nx*(j+ny/3)].ds = -SLO;
					if(sol[i+nx*2/3+nx*(j+ny/3)].nu>0) sol[i+nx*2/3+nx*(j+ny/3)].ds = SLO;
					if(sol[i+nx/3+nx*(j+ny*2/3)].nu>0) sol[i+nx/3+nx*(j+ny*2/3)].ds = SLO;
					if(sol[i+nx*2/3+nx*(j+ny*2/3)].nu>0) sol[i+nx*2/3+nx*(j+ny*2/3)].ds = -SLO;
				}

	if(animesh!=0.0)
		for(i=-xmesh/animesh/2; i<=xmesh/animesh/2; i++)
			for(j=-ymesh/animesh/2; j<=ymesh/animesh/2; j++)
				if(i*i+j*j <= xmesh*ymesh/animesh/animesh/4.0) {
					if(sol[i+nx/3+nx*(j+ny/3)].nu>0) sol[i+nx/3+nx*(j+ny/3)].sin = -ANI;
					if(sol[i+nx*2/3+nx*(j+ny/3)].nu>0) sol[i+nx*2/3+nx*(j+ny/3)].sin = ANI;
					if(sol[i+nx/3+nx*(j+ny*2/3)].nu>0) sol[i+nx/3+nx*(j+ny*2/3)].sin = ANI;
					if(sol[i+nx*2/3+nx*(j+ny*2/3)].nu>0) sol[i+nx*2/3+nx*(j+ny*2/3)].sin = -ANI;
				}

}

/********************************************************************/

getmod2(cellmesh,stamesh,animesh)
double cellmesh,stamesh,animesh;
/*this one does the checkerboard stuff*/
{
	int i;

	printf("getmod: cellmesh,animesh= %f %f\n",cellmesh,animesh);
	if(cellmesh!=0.0) slobox(SLO,cellmesh,cellmesh);
	if(stamesh!=0.0) stabox(STA,stamesh);
	if(animesh!=0.0) anibox(ANI,animesh,animesh);

}

/*************************************************************************/

getmod3(cellmesh,stamesh,animesh)
double cellmesh,stamesh,animesh;
/*this one does the checkerboard sinusoids*/
{
	int i;

	printf("getmod: cellmesh,animesh= %f %f\n",cellmesh,animesh);
	/*mult slo by 2 so that slo is the rms sinusoid amplitude*/
	if(cellmesh!=0.0) slosin(SLO*2.0,cellmesh,cellmesh);
	if(stamesh!=0.0) stasin(STA,stamesh);
	if(animesh!=0.0) anisin(ANI*2.0,animesh,animesh);

}

/***********************************************************************/

getmod4(cellmesh,stamesh,animesh)
double cellmesh,stamesh,animesh;
/*this one does the stripes stuff*/
{
	int i;

	/*make stripes of opposite directions*/
	printf("getmod: cellmesh,animesh= %f %f\n",cellmesh,animesh);
	if(cellmesh!=0.0) slobox(SLO,cellmesh,0);
	if(stamesh!=0.0) stabox(STA,stamesh);
	if(animesh!=0.0) anibox(ANI,0,animesh);

}

/***********************************************************************/

getmod5(cellmesh,stamesh,animesh)
double cellmesh,stamesh,animesh;
/*this one does a test of the velocity diaride anomaly*/
{
	int i,j,k;

	/*zero it all first*/
	for(i=0;i<nx*ny;i++) sol[i].ds = 0.0;
	/*a box from lat=37 to 43, lon=20 to 21*/
	/*for(i=80; i<84; i++) for(j=28; j<52; j++)
		sol[i+nx*j].ds = SLO;*/
	
	/*chinese box, 105-107, 30-40*/
	/*at 1/2 deg cells*/
	for (i=80; i<84; i++) for(j=40; j<80; j++) {
		sol[i+nx*j].ds = SLO;
		sol[i+nx*j].cos = -(36./44.)*SLO;
	}
	
}

/*******************************************************************/
getmod6(cellmesh,stamesh,animesh)
double cellmesh,stamesh,animesh;
/*this one does a test of the diaride anisotropy anomaly*/
{
	int i,j,k;

	/*zero it all first*/
	for(i=0;i<nx*ny;i++) sol[i].ds = 0.0;
	/*a box from lat=37 to 43, lon=20 to 21*/
	for(i=80; i<84; i++) for(j=28; j<52; j++)
		sol[i+nx*j].cos = -ANI;

}
/****************************************************************/
getmod7(cellmesh,stamesh,animesh)
double cellmesh,stamesh,animesh;
/*constant anisotropy, no velocity*/
{
	int i,j,k;

	/*zero*/
	for(i=0;i<nx*ny;i++) sol[i].ds=0.0;
	for(i=0;i<nx*ny;i++) sol[i].cos= ANI;
}
/************************************************************/
getmod8(cellmesh,stamesh,animesh)
double cellmesh,stamesh,animesh;
/*constant velocity, no anisotropy*/
{
	int i,j,k;

	for(i=0;i<nx*ny;i++) sol[i].ds=SLO;
}

/***********************************************************************/

getmod9(cellmesh,stamesh,animesh)
double cellmesh,stamesh,animesh;
/*this one does a test of the anisotropic rectangle*/
{
	int i,j,k;

	/*zero it all first*/
	for(i=0;i<nx*ny;i++) sol[i].ds = 0.0;
	
	/*chinese box, 105-107, 30-40*/
	/*at 1/2 deg cells*/
	/*for (i=80; i<84; i++) for(j=40; j<80; j++) {
		sol[i+nx*j].ds = SLO;
		sol[i+nx*j].cos = -(36./44.)*SLO;
	}
	*/
	
	/*chinese box, 105-109, 30-42*/
	/*at 1/4 deg cells*/
	for (i=160; i<176; i++) for(j=80; j<128; j++) {
		sol[i+nx*j].ds = SLO;
		sol[i+nx*j].cos = -SLO;
	}

	
}





/**********************support routines ***************************/

anibox(slo,xxmesh,yymesh)
/*checkerboard boxe*/
/*the boxes dont add*/
/*the mesh sizes are in mesh/degree */
float slo;
float xxmesh,yymesh;
{
	int i, j, k;
	int nxmin,nxmax,nymin,nymax;
	float x,y;
		 
	/*do anisotropy */
	for(j=0; j<ny; j++)
		for(i=0; i<nx; i++)
			if(sol[i+nx*j].nu>0) {
				/*sol[i+nx*j].sin =
							slo * ( (( ((int)(i*xxmesh/xmesh)) + ((int)(j*yymesh/ymesh)) )%2)*2-1);
				sol[i+nx*j].cos = slo;*/
				/*k is either -1 or 1*/
				k= ( (( ((int)(i*xxmesh/xmesh)) + ((int)(j*yymesh/ymesh)) )%2)*2-1);
				if(k==1) {
					sol[i+nx*j].sin = 0.0;
					sol[i+nx*j].cos = slo;
				}
				if(k==-1) {
					sol[i+nx*j].sin = 0.0;
					sol[i+nx*j].cos = -slo;
				}
	
			}
}
		
slobox(slo,xxmesh,yymesh)
/*checkerboard boxe*/
/*the boxes dont add*/
/*the mesh sizes are in mesh/degree */
float slo;
float xxmesh,yymesh;
{
	int i, j;
	int nxmin,nxmax,nymin,nymax;
	float x,y;

	/*do the slownesses */
	printf("xxmesh,yymesh= %f %f\n",xxmesh,yymesh);
	for(j=0; j<ny; j++)
		for(i=0; i<nx; i++)
			if(sol[i+nx*j].nu>0) 
				sol[i+nx*j].ds =
					slo* ( (( ((int)(i*xxmesh/xmesh)) + ((int)(j*yymesh/ymesh)) )%2)*2-1);

	
}

stabox(stnamp,mesh)
/*checkerboard stations*/
double stnamp;
double mesh;
{
	int i,j,k;
	fprintf(stderr,"stabox: stnamp,mesh %f %f\n",stnamp,mesh);
	for(i=0;i<nstat;i++) {
		if(station[i].lat<ymin) continue;
		if(station[i].lat>ymax) continue;
		if(station[i].lon<xmin) continue;
		if(station[i].lon>xmax) continue;
		station[i].delay = stnamp*(int)(((unsigned)(station[i].lat*mesh)%2)*2-1)*(int)(((unsigned)(station[i].lon*mesh)%2)*2-1);
	}
}
		
anisin(slo,xxmesh,yymesh)
/*checkerboard boxe*/
/*the boxes dont add*/
/*the mesh sizes are in mesh/degree */
float slo;
float xxmesh,yymesh;
{
	int i, j, k;
	int nxmin,nxmax,nymin,nymax;
	float x,y;
		 
	/*do anisotropy */
	for(j=0; j<ny; j++)
		for(i=0; i<nx; i++)
			if(sol[i+nx*j].nu>0) {				
				sol[i+nx*j].sin= 0.0;
				sol[i+nx*j].cos= slo*sin(PI*i*xxmesh/xmesh)*sin(PI*j*yymesh/ymesh);
			}
}
		
slosin(slo,xxmesh,yymesh)
/*checkerboard boxe*/
/*the boxes dont add*/
/*the mesh sizes are in mesh/degree */
float slo;
float xxmesh,yymesh;
{
	int i, j;
	int nxmin,nxmax,nymin,nymax;
	float x,y;

	/*do the slownesses */
	for(j=0; j<ny; j++)
		for(i=0; i<nx; i++)
			if(sol[i+nx*j].nu>0) 
				sol[i+nx*j].ds = slo*sin(PI*i*xxmesh/xmesh)*sin(PI*j*yymesh/ymesh);

	
}

stasin(stnamp,mesh)
/*checkerboard stations*/
double stnamp;
double mesh;
{
	int i,j,k;
	fprintf(stderr,"stabox: stnamp,mesh %f %f\n",stnamp,mesh);
	for(i=0;i<nstat;i++) {
		if(station[i].lat<ymin) continue;
		if(station[i].lat>ymax) continue;
		if(station[i].lon<xmin) continue;
		if(station[i].lon>xmax) continue;
		station[i].delay = stnamp*sin(PI*station[i].lat*mesh)*sin(PI*station[i].lon*mesh);
	}
}


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
	*length = sqrt(dx * dx + dy * dy) * scale;
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
		len = sqrt(dx * dx + dy * dy) * scale;
		
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

	return;
}
