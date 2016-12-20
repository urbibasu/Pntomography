#include <stdio.h>
#include <math.h>
#include <string.h>
/*next to lines are for mmap()*/
#include <sys/types.h>
#include <sys/mman.h>
/*next two for lseek()*/
#include <sys/types.h>
#include <unistd.h>

#define COMMAND	"usage readpha [-file filename] [-rays]"

double sqrt(), fabs();

#include "data14.h"

struct rec *phase;

main(argc, argv)
int argc; char **argv;
{
	int i, j;
	int count=0;
	int fdrec;
	float temp;
	int npha;
	int rays=0;
	char file[80]="../data/phases.b";
	
	getarg(argc,argv,"-file","%s",&file);
	rays=getarg(argc,argv,"-rays","") +1;
	fprintf(stderr,"file= %s\n",&file);
	fprintf(stderr,"rays= %d\n",rays);

	fprintf(stderr,"starting\n");
	fdrec=open(&file,2);
	phase = (struct rec *) mmap(0,lseek(fdrec,0L,SEEK_END),PROT_READ|PROT_WRITE, MAP_SHARED,fdrec,0);
	npha = lseek(fdrec,0L,SEEK_END)/sizeof(struct rec);
	fprintf(stderr,"fdrec,phase,npha= %d %d %d\n",fdrec,phase,npha);

	for(i=j=0;i<npha;i++) {
		if(phase[i].quality>0) {
			j++;
			if(!rays)
				printf("%d %d %f %f %f %f %f %f\n", phase[i].stnno,phase[i].evnno,phase[i].dtime,phase[i].quality,phase[i].ttime,phase[i].offset,phase[i].azse,phase[i].quality);
			else {
				printf("%f %f\n", phase[i].slon,phase[i].slat);
				printf("%f %f\n", phase[i].rlon,phase[i].rlat);
				printf(">\n");
			}
		}
	}
	fprintf(stderr,"total rays %d\n",i);
	fprintf(stderr,"rays used %d\n",j);
	fprintf(stderr,"rdpha done now \n\n");
}


