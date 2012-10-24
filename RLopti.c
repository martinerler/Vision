/*
 * RLopti.c
 *
 *  Created on: 24.10.2012
 *      Author: martin
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "RLopti.h"


// zum debuggen:
int ratio1 = 0 , ratio2 = 0;


/*
 * summiert die Grauwerte des gegebenen Bilder anhand der Pixelliste und eines Referenzpunktes
 *
 * pix...Bild indem die Summation erfolgen soll
 * nx,ny..Bilddimension von pix
 * listx,listy...jeweils ein Array mit den Pixelpositionen
 * refx,refy...der gegebene Referenzpunkt: 	die Pixelposition auf die sich die Summation bezieht
 * 											oder: da wo die Kernposition=0,0 ist
 * asize...die listengröße
 * return...die Summe
 */
float inline sumfromArray(float* pix, int nx, int ny, int* listx, int* listy, int refx,int refy, int asize)
{
	int i;
	int index;
	float g=0;
	for(i=0;i<asize;i++) {
		index = getPixelIndexBend(refx + listx[i], refy + listy[i], nx,ny);
		g += pix[index];

	}
	return g;
}


/*
 * calculates the number of Pixels with g > 0 of a given Image
 *
 * float* pix  		input image
 * nx, ny			image dimensions
 * return value		number of pixels > 0
 */
int getKernelPixelsNotNull(float* pix, int nx, int ny) {
	int i,n=0;
	for(i=0;i<nx*ny;++i) {
		if(pix[i]> 0.001)
			n++;
	}
	return n;
}


/*
 * erzeugt eine Punkteliste aus einem Bild mit allen Pixeln mit Grauwert größer 0.
 *
 * pix...bildpointer
 * nx, ny..bildgröße
 * listx, listy...integer array der Punkte, aufgeteilt in x und y koordinaten
 * g....ein float array mit den enthaltenen grauwerten der pixel > 0.1
 * hier wird keine struktur verwendet, sondern eben die zwei integer arrays
 *
 */
void calcArrayfromKernel(float* pix, int nx, int ny, int* listx, int* listy, float* g, int asize) {

	int i,j,n=0;
	int hnx, hny;
	if(nx%2==0 || ny%2==0) {
		printf("calcArrayfromKernel: Kerngröße nicht ungerade\n");
		return;
	}
	hnx=nx/2;
	hny=ny/2;
	for(j=0;j<ny;++j) {
		for(i=0;i<nx;++i) {
			if( pix[getPixelIndex(i,j,nx,ny)] > 0.001) {
				listx[n] = i - hnx;
				listy[n] = j - hny;
				g[n] = pix[getPixelIndex(i,j,nx,ny)];
				n++;
			}


		}
	}
}


/**************************************************
normalizes the kernel image to the sum of pixels = 1
parameter:
			int   : width of the kernel
			int   : height of the kernel
			float*: the matrix to be normalized
Author: Thomas Schwarzbauer
Date: 01.04.2011
***************************************************/
void normalize(int nx, int ny, float* pixel){
	int i;
	int performance;
	float sum = 0;
	performance = nx*ny;
	for(i = 0; i < performance; i++){
		sum+=pixel[i];
	}
	for(i = 0; i < performance; i++){
		pixel[i] /= sum;
	}
}




void listfilter_stat(	float* upixel, float* vpixel, float* konto,float k, int unx, int uny,
						int* listx1, int* listy1, int asize1) {
	int x,y;
	int index;
	float hsize_inv = 1 / (float)(asize1);	//der quotient zur normalisierung entspricht der kerngröße, also asize1
											//hier gleich den kehrwert bilden, weil die division dann billiger ist
	float sum=0.0;

	int debug;

	for(y=0; y<uny;y++) {
		for(x=0;x<unx;x++) {
			index = getPixelIndex(x,y,unx,uny);

			if( konto[index]>=k ) {
				sum = sumfromArray(upixel, unx,uny,listx1,listy1,x,y,asize1);
				vpixel[index] = sum * hsize_inv;
				ratio1++;
			}
			else {
				vpixel[index] = upixel[index];
				ratio2++;
			}
		}
	}
}


void convolve_stat(	int unx, int uny, float* upixel,
					int hnx, int hny, float* hpixel,
					int vnx, int vny, float* vpixel,
					float* konto, float k)
{
	int x,y;	  				// Picture iterator
	int index1,index2,index3;	// array index optimization
    int kerny,kernx;	    	// kernel iterator x and y
    float temp;
	for(y=0;y<vny;y++) {
       for(x=0;x<vnx;x++)
       {
    	   index1 = getPixelIndexBend(x,y,vnx,vny);

    	   // nur falten, wenn der Schwellwert erreicht
    	   if(konto[index1] >= k) {

    		   for(kerny=0;kerny<hny;kerny++) {
				   for(kernx=0;kernx<hnx;kernx++) {
					   index2 =	getPixelIndexBend(x -hnx/2 + kernx, y - hny/2 + kerny,unx,uny);
					   index3 = getPixelIndexBend(kernx,kerny,hnx,hny);
					   temp += upixel[index2] * hpixel[index3];
				   }
			   }
			   vpixel[index1] = temp;
			   temp=0.0;
			   ratio1++;	// debugwert
    	   }
    	   else {
    		   // sonst nur überspringen
    		   vpixel[index1] = upixel[index1];
    		   ratio2++;
    	   }
       }
	}
}

/**************************************************************************
perform Richardson-Lucy deconvolution with statistic improvement
parameters:
			float *u	deblurred image (result)
			int ux		columns of the image
			int uy		rows of the image
			float *f	blurred image (initial image of iteration)

			int* listx1 coordinate array of generic mainkernel
			int* listy1 coordinate array of generic mainkernel
			int asize1  size of coordinate array of generic mainkernel
			int* listx1 coordinate array of generic overlap kernel left
			int* listy1 coordinate array of generic overlap kernel left
			int asize1  size of coordinate array of generic overlap kernel left
			int* listx1 coordinate array of generic overlap kernel right
			int* listy1 coordinate array of generic overlap kernel right
			int asize1  size of coordinate array of generic overlap kernel left

			int iter	number of iterations of RL algorithm

return:		1 if successful
			0 else

**************************************************************************/
int performRL_stat		(	float *u, int unx, int uny, float *f,
							float *h, int hnx, int hny,
							float k, int iter)
{

	int i, j;		// counter variables
	char filename[80];
	//float *hAdj = (float*) malloc((hnx*hny)*sizeof(float));		// storage for adjoint of matrix h
	float *denom = (float*) malloc((unx*uny)*sizeof(float));		// storage for denominator
	float *fraction = (float*) malloc((unx*uny)*sizeof(float));		// storage for fraction
	float *fractionp = (float*) malloc((unx*uny)*sizeof(float));	// storage for fraction
	float *product = (float*) malloc((unx*uny)*sizeof(float));		// storage for product
	float *tempResult = (float*) malloc((unx*uny)*sizeof(float));	// storage for temporary result
	//transpose(hnx, hny, h, hAdj);									// calculate adjoint of matrix h

	// Liste aus Kern generieren
	int asize1= getKernelPixelsNotNull(h,hnx,hny);
	float* g = (float*) malloc(asize1*sizeof(float)); 	// kernel-grauwerte werden ignoriert, alle als konstant angenommen
	int* listx1 = (int*) malloc(asize1*sizeof(int));
	int* listy1 = (int*) malloc(asize1*sizeof(int));
	calcArrayfromKernel(h,hnx,hny,listx1,listy1,g,asize1);
	// wenn der Kern weiterverwendet würde, normalisieren
	normalize(hnx,hny,h);

	float *diff = (float*) malloc((unx*uny)*sizeof(float));				//storage for difference - used for statistical performance improvement
	float *tempResultp = (float*) malloc((unx*uny)*sizeof(float));	// storage for temporary result
	int overallratio1=0,overallratio2=0;

	for(i=0;i<unx*uny;i++)
		diff[i]=100;

	//float k=0.1;

	// initialize RL: u0 = f
	for(i=0; i<(unx*uny); i++){
		u[i] = f[i];
	}

	// core iterative RL algorithm
	for(i=0; i<iter; i++){

		// convolution: h * uk
		//listfilter_stat(u,denom,diff,q,unx,uny,listx1,listy1,asize1);
		convolve_stat(unx,uny,u,hnx,hny,h,unx,uny,denom,diff,k);


		// calculate fraction: f / (h * uk)
		for(j=0; j<(unx*uny); j++){

			if( diff[j] >= k) {
				fraction[j] = (float)(f[j] / denom[j] );	// pixelwise division
				fractionp[j] = fraction[j];
			}
			else {
				fraction[j] = fractionp[j];
			}
		}

		// convolution: hAdj * (f / (h* uk))
		//listfilter_stat(fraction,product,diff,q,unx,uny,listx1,listy1,asize1);
		convolve_stat(unx,uny,fraction, hnx,hny, h, unx,uny,product,diff,k);

		// calculate product: (hAdj * (f / (h* uk))) * uk
		for(j=0; j<(unx*uny); j++){

			// update new uk
			if(diff[j] >= k) {
				tempResult[j] = (product[j] * u[j]);	// pixelwise multiplication
				tempResultp[j] = tempResult[j];
			}
			else
				tempResult[j] = tempResultp[j];
		}

		for(j=0; j<(unx*uny); j++) {

			// calc the greyvalue diff from uk-1 to uk
			diff[j] = fabs( u[j] - tempResult[j]);
			// jede 10te iteration voll falten
			if(i%10==0) diff[j] = 100;

			// update the result
			u[j] = tempResult[j];
		}

		overallratio1+=ratio1;
		overallratio2+=ratio2;
		printf("ratio1 = %i ratio2 = %i ratio = %f overallratio= %f diff0 = %f\n", ratio1,ratio2, ratio1/(float)ratio2, overallratio1/(float)overallratio2, diff[0]);
		ratio1 = 0;
		ratio2 = 0;

		//sprintf(filename, "debugbilder/rl_%i.pgm", i);
		//writepgm(filename,unx,uny,255,&diff);
	}

	free(denom);
	free(fraction);
	free(product);
	free(tempResult);

	printf("Richardson-Lucy Stat Deconvolution finished");

	return 1;
}


