#include <stdio.h>
#include <math.h>

double sqrt(),fabs(),atan(),atan2();

#include "data14.h"

double sqrt();

main(argc, argv)
int argc; char **argv;
{
	struct solution soln;
	double vel,velbias,velerr;
	double ani,anibias,anierr;
	double dir,dirbias,direrr;
	int i, j;
	int count=0;
	int kntzer=0;
	int rmsnm=0;
	double rms=0.0;
	double rav=0.0;
	double ave=0.0;
	int fdrec;
	double xloc,yloc;
	int nx,ny;
	
	char file[150]="../data/sol.b";

	float xmin=0.0,xmax=360.0,ymin=-90.0,ymax=90.0;
	float xmesh=1.0,ymesh=1.0;
	float vn=8.0;
	int nxskip=1,nyskip=1;

	getarg(argc,argv,"-file","%s",&file);
	fprintf(stderr,"file= %s\n",&file);
	getarg(argc,argv,"-region","%f/%f/%f/%f",&xmin,&xmax,&ymin,&ymax);
	getarg(argc,argv,"-mesh","%f/%f",&xmesh,&ymesh);
	getarg(argc,argv,"-vn","%f",&vn);
	getarg(argc,argv,"-nxskip","%d",&nxskip);
	getarg(argc,argv,"-nyskip","%d",&nyskip);

	nx= (xmax-xmin)*xmesh;
	ny= (ymin-ymax)*ymesh;
	fprintf(stderr,"file,xmin,xmax,ymin,ymax,xmesh,ymesh,vn,nxskip,nyskip,nx,ny=\n %s %f %f %f %f  %f %f %f %d %d %d %d\n",
		&file,xmin,xmax,ymin,ymax,xmesh,ymesh,vn,nxskip,nyskip,nx,ny);

	fprintf(stderr,"file= %s\n",&file);
	fdrec=open(&file,2);
	fprintf(stderr,"fdrec= %d\n",fdrec);
	if(fdrec<=0) errprint("bad file; exiting\n");
	fprintf(stderr,"starting\n");
	
	/*loop thru sol.b*/
	for(i=0;read(fdrec,&soln,sizeof(struct solution));i++)
	{
		/*get rid of undetermined cells*/
		if(soln.nu==0) continue;
		
		/*skip every nskip cells*/
		if((i%nx)%nxskip!=0) continue;
		if((i/nx)%nyskip!=0) continue;

		/*fprintf(stdout,"ds,cos,sin= %f %f %f %f\n",soln.ds,soln.cos,soln.sin);*/
		
		/*update statistics of slowness perturbations*/
		count++;
		vel=soln.ds;
		if(vel==0.0) kntzer++;
		else
		{
			ave += vel;
			rms += vel*vel;
			rav += fabs(vel);
			rmsnm++;
		}
		
		/*xloc & yloc are now the cell centers*/
		xloc=xmin + (float)(i%nx)/xmesh + 0.5/xmesh;
		yloc=ymin + (float)(i/nx)/ymesh + 0.5/ymesh;
		
		/*vel is the velocity perturbation, velbias is velocity bias, velerr is standard error*/
		vel *= -vn*vn;
		velbias = -vn*vn*(soln.sum-soln.ds);
		velerr = vn*vn*sqrt(soln.sum2);

		/*fprintf(stdout,"sum2,sum2cos,sum2sin= %g %g %g\n",soln.sum2,soln.sum2cos,soln.sum2sin);*/
		
		/*ani is the anisotropy velocity perturbation*/
		ani= vn*vn*sqrt(soln.cos*soln.cos+soln.sin*soln.sin);
		anibias = vn*vn*sqrt(soln.sumcos*soln.sumcos + soln.sumsin*soln.sumsin);
		anibias -= ani;
		anierr = soln.sum2cos*soln.cos*soln.cos + soln.sum2sin*soln.sin*soln.sin + 2*soln.sum2cossin*soln.cos*soln.sin;
		anierr = vn*vn*sqrt(anierr)/sqrt(soln.cos*soln.cos+soln.sin*soln.sin);
		
		/*dir is the atan(sin/cos) w/ coefficient 180/pi/2*/
		dir=(28.6479)*atan2(soln.sin,(soln.cos+.000000000001));
		dir = 90.0 - dir; /*convert to gmt coords (pos counterclockwise from E*/
		dir -= 90.0; /*plot max velocity not max slowness*/
		dirbias = (28.6479)*atan2(soln.sumsin,soln.sumcos+.000000000001);
		dirbias = 90.0 - dirbias; /*convert to gmt coords (pos counterclockwise from E*/
		dirbias -= 90.0; /*plot max velocity not max slowness*/
		dirbias -= dir;
		direrr = 28.6419*sqrt(soln.sum2sin*soln.cos*soln.cos + soln.sum2cos*soln.sin*soln.sin + 2*soln.sum2cossin*soln.cos*soln.sin);
		direrr /= (soln.cos*soln.cos + soln.sin*soln.sin);

		fprintf(stdout,"%9.4f %9.4f %5d %7.7f %7.7f  %7.7f %7.7f %7.7f %7.7f %7.7f  %7.7f %7.7f\n",
			xloc,yloc,soln.nu,vel,velbias,
			velerr,ani,anibias,anierr,dir,
			dirbias,direrr);		
	}
	fprintf(stderr,"\n %d zeros of %d total\n",kntzer,count);
	fprintf(stderr,"i= %d\n",i);
	fprintf(stderr,"rms= %f , rmsnm= %d \n",sqrt(rms/((float)rmsnm+.0001)),rmsnm);
	fprintf(stderr,"aveabs= %f , rmsnm= %d \n",(rav/((float)rmsnm+.0001)),rmsnm);
	fprintf(stderr,"ave = %f\n\n",ave/(rmsnm+.0001));
}

