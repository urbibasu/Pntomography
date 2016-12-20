#include <stdarg.h>
#include <stdio.h>
#include <string.h>

//getarg(va_alist)
getarg(int argc, ...)
//va_dcl
/* thearn, 8/90 */
/* gets parameters out of command line, and interprets followin commands via format */
/* example calls: */
/*	getarg(argc,argv,"-count","%d",&count); */
/*	getarg(argc,argv,"-flag",""); */
/*	getarg(argc,argv,"-par2","%d%d",&i,&j); */
/* first sting is search string */
/* second string is format for everything else */
/* remaining strings are the parmaters to assign */
/* 10 arguments max per format*/
/* 10 arguments maximum per format */

/* returns:	option found:		number of parameters found (0 to 10) */
/*								11 if more than 10 arguments in format*/
/*			option not found:	-1*/

{
	va_list         ap;
	int             i, j, k;
//	int             argc;
	char          **argv;
	char           *format, *string;
	int             rtnval = 0;	/* return value */
	int 			nargs=0;
	void			*args[10];

	/* get starting parameters */
	va_start(ap, argc);
	//argc = va_arg(ap, int);
	argv = va_arg(ap, char **);
	string = va_arg(ap, char *);
	format = va_arg(ap, char *);
	fprintf(stderr,"string,format= %s %s\n",string,format);
    	fprintf(stderr,"So what are we looking actually ???????\n",string,format);
    	fprintf(stderr,"R u sure ?????\n",string,format);
    	fprintf(stderr,"heyyyyyy\n",string,format);
    

	/*count the parameters needed*/
	for(i=0;format[i]!='\0';i++) if(format[i]=='%') nargs++;
	/*fprintf(stderr,"nargs=%d\n",nargs);*/
	if(nargs>10) errprint("getarg: too many arguments for this parameter");

	/*copy the argument pointers over, they aren't filled yet*/
	for(i=0;i<nargs;i++) args[i] = va_arg(ap,int*);
	fprintf(stderr,"i,nargs=%d,%d\n",i,nargs);

	/* now look at argv parameters */
	rtnval=-1;
	for(i=1; i<argc; i++) {
		if (strcmp(argv[i],string) == 0) {	/* see if argument exists */
			if(nargs==0)
				rtnval=0;	/*error, too few args*/
			else if(nargs==1)
				rtnval=sscanf(argv[i+1],format,
					args[0]);
			else if(nargs==2)
				rtnval=sscanf(argv[i+1],format,
					args[0],args[1]);
			else if(nargs==3)
				rtnval=sscanf(argv[i+1],format,
					args[0],args[1],args[2]);
			else if(nargs==4)
				rtnval=sscanf(argv[i+1],format,
					args[0],args[1],args[2],args[3]);
			else if(nargs==5)
				rtnval=sscanf(argv[i+1],format,
					args[0],args[1],args[2],args[3],args[4]);
			else if(nargs==6)
				rtnval=sscanf(argv[i+1],format,
					args[0],args[1],args[2],args[3],args[4],
					args[5]);
			else if(nargs==7)
				rtnval=sscanf(argv[i+1],format,
					args[0],args[1],args[2],args[3],args[4],
					args[5],args[6]);
			else if(nargs==8)
				rtnval=sscanf(argv[i+1],format,
					args[0],args[1],args[2],args[3],args[4],
					args[5],args[6],args[7]);
			else if(nargs==9)
				rtnval=sscanf(argv[i+1],format,
					args[0],args[1],args[2],args[3],args[4],
					args[5],args[6],args[7],args[8]);
			else if(nargs==10)
				rtnval=sscanf(argv[i+1],format,
					args[0],args[1],args[2],args[3],args[4],
					args[5],args[6],args[7],args[8],args[9]);
			else
				rtnval=sscanf(argv[i+1],format,
					args[0],args[1],args[2],args[3],args[4],
					args[5],args[6],args[7],args[8],args[9])+1;
		}
	}
	va_end(ap);
	return rtnval;
}
