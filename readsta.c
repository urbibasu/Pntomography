#include <stdio.h>
#include <math.h>
#include <string.h>
/*next to lines are for mmap()*/
#include <sys/types.h>
#include <sys/mman.h>
/*next two for lseek()*/
#include <sys/types.h>
#include <unistd.h>

#define COMMAND	"usage readsta [-file filename]"

double sqrt(),fabs(),lgamma(),log(),exp();

#include "data14.h"

struct locations *station;

main(argc, argv)
int argc; char **argv;
{
	double value;
	int i, j;
	int count=0;
	int kntzer=0;
	int rmsnm=0;
	float rms=0.0;
	float rav=0.0;
	float sum2=0.0,sum=0.0;
	int sumcount=0;
	int fdsta;
	float fval,chival;
	double fprob(),chiprob();
	int nstat;
	char file[80]="../data/stations.b";

	getarg(argc, argv,"-file","%s",&file);
	fprintf(stderr,"%d, done with getarg\n",argc);
	fprintf(stderr,"file= %s\n",&file);

	fprintf(stderr,"starting\n");
	fdsta=open(&file,2);
	station = (struct locations *) mmap(0,lseek(fdsta,0L,SEEK_END),PROT_READ|PROT_WRITE, MAP_SHARED,fdsta,0);
	nstat = lseek(fdsta,0L,SEEK_END)/sizeof(struct locations);
	fprintf(stderr,"station,nstat,fdsta= %d %d %d\n",station,nstat,fdsta);

	/*first get the total count and the total sum of squares*/
	/*then the F calculation can be done*/
	count=0; sum2=0.0; sumcount=0;
	for(i=0;i<nstat;i++) {
		if(station[i].count!=0) {
			count++;
			sum += station[i].delay;
			sum2 += station[i].delay*station[i].delay;
			sumcount += station[i].count;
		}
	}
	
	/*print out 1)name 2)number 3)lat 4)lon 5)elev 6)delay 7)count 8)bias 9)weight 10)sd */
	for(i=0;i<nstat;i++) {
		if(station[i].count!=0) {
			printf("%7s %5d %5.2f %5.2f %7.2f %5.2f %4d %8.3f %8.3f %6.4f %g %g\n",
				station[i].name,station[i].number,station[i].lat,station[i].lon,station[i].elev,
				station[i].delay,station[i].count,(station[i].sum-station[i].delay),station[i].wei,sqrt(fabs(station[i].sum2)),
				station[i].sum,station[i].sum2);
		}
	}
	fclose(stdout);
	/*note: F-stat must be compared to F-table w/ nbootv and nbootv*count dof*/
	fprintf(stderr, "count,sum(delay),sum(sum2),sum(count) = %d %f %f %d\n",count,sum,sum2,sumcount);
	fprintf(stderr,"average= %f\n",sum/count);
	fprintf(stderr,"rms station delay = %f\n",sqrt(sum2/count));
	fprintf(stderr,"rdstn done now \n\n");
}


