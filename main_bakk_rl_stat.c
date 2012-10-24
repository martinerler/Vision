/*
 * rl.c
 *
 *  Created on: 14.05.2011
 *      Author: martin
 */

#include <stdlib.h>
#include <stdio.h>
#include "pgmio.h"
//#include "functions.h"
#include "RLopti.h"
#include <sys/time.h>

int main (int argc, char *argv[]) {

	float *pixel;
	float *pixelkern1;
	float *out;
	int nx,ny, d;
	int hx,hy;
	struct timeval start, end;

	printf("usage: rl [input image] [kernel image1] [output image] [k] [iteration number]   \n\n");

	// if no 6. command line parameter: set iterations default 5
	if(argv[6]==0)
			argv[6] = "5";
	printf("Input Image: %s \n", argv[1]);
	printf("Kernel: %s \n", argv[2]);
	printf("Output Image: %s \n", argv[3]);
	printf("k: %s \n", argv[4]);
	printf("Iterations: %s \n", argv[5]);

	float k  = (float)atof(argv[4]);
	readpgm(argv[1], &nx, &ny, &pixel);
	readpgm(argv[2], &hx, &hy, &pixelkern1);



	out = (float *) malloc((nx*ny) * sizeof(float));

	gettimeofday(&start, NULL);

	performRL_stat(	out, nx, ny, pixel,
					pixelkern1,hx,hy,
					k, atoi(argv[5]));


	//normalize(hx,hy,pixelkern1);
	//performRL_convolve(out,nx,ny,pixel,pixelkern1,hx,hy,atoi(argv[6]));
	// normaler RL
	/*float* pixelkernbig = (float *) malloc((nx*ny) * sizeof(float));
	expandkernel(d, d, pixelkern1, nx, ny, pixelkernbig);
	pixelSwitch(nx,ny,pixelkernbig);
	performRL(out,nx,ny,pixel,pixelkernbig,nx,ny,10);*/

	// benchmark code
	gettimeofday(&end, NULL);
	long long timegone =   (end.tv_sec * (unsigned int)1e6 +   end.tv_usec) -
					 (start.tv_sec * (unsigned int)1e6 + start.tv_usec);
	printf("\ntime = %lld ms \n", timegone/1000);

	//writepgm16(argv[5],nx,ny,65535,&out);
	writepgm(argv[3],nx,ny,255,&out);
	printf("Iterations: %s \n", argv[5]);
	fflush(stdout); /* force it to go out */
	return 0;
}
