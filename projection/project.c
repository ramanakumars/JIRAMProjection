#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "cspice/include/SpiceUsr.h"

#define JUPITER 599
#define FRAME_HEIGHT 128
#define FRAME_WIDTH  432

double FOCAL_LENGTH=159.82033e-3; // in m
double PIXEL_SIZE=38.0e-6; // 38 microns
double CENTER_X=216.5;
double CENTER_Y=64.5;
double PIXEL_SCALE=0.00023778; // rad/pixel

double f1;

#define M_BAND 0
#define L_BAND 1

void vec2pix(double *v, double *pix){
	register double alpha, x, y;

	alpha  = v[2]/f1;
	x      = CENTER_X - v[1]/alpha;
	y      = CENTER_Y - v[0]/alpha;

	pix[0] = x; pix[1] = y;
}

// subtract two 3D vectors: out = x - y
void subtract3D(double *x, double *y, double *out) {
	for(int i=0; i<3; i++) {
		out[i] = x[i] - y[i];
	}
}

// 3D matrix dot product: out = A dot b
void matmul3D(double *A, double *b, double *out) {
	for(int i=0; i<3; i++) {
		out[i] = 0.;
		for(int j=0; j<3; j++) {
			out[i] += A[i*3+j]*b[j];
		}
	}
}


int* get_image_mask(double *lat, double *lon, int nlat, int nlon, 
		double *scloc, double *jup2cam, double *et, int nframes) {

	f1 = (FOCAL_LENGTH)/(PIXEL_SIZE);

	int *mask;
	double surface_point[3], vec[3], junocam_vec[3], pix[2];
	double vnorm, sdotvec;
	int xx, yy;

	mask = malloc((nlat*nlon)*sizeof(int));
	
	for(int jj=0; jj<nlat; jj++) {
		double latj, loni; 
		if((jj%100)==0) printf("\r %4d/%4d ", jj, nlat);
		latj = lat[jj];
		for(int ii=0; ii<nlon; ii++) {
			mask[jj*nlon+ii] = 0;
			loni = lon[ii];

			// find the vector to the point on the surface
			// in the Jupiter frame
			srfrec_c(JUPITER, loni, latj, surface_point);
			for(int i=0; i<nframes; i++) {
				// the vector from the spacecraft to the surface
				// in the JUPITER frame
				subtract3D(surface_point, &scloc[i*3], junocam_vec);

				// make sure we are not looking through the planet
				sdotvec = surface_point[0]*junocam_vec[0] + 
					surface_point[1]*junocam_vec[1] + 
					surface_point[2]*junocam_vec[2];

				// project this vector into the JUNOCAM frame
				matmul3D(&jup2cam[i*9], junocam_vec, vec);

				// normalize
				vnorm = sqrtf(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
				vec[0] = vec[0]/vnorm;
				vec[1] = vec[1]/vnorm;
				vec[2] = vec[2]/vnorm;

				// find the pixel coordinate corresponding
				// to this vector
				vec2pix(vec, pix);
			
				// check if the pixel falls inside the photoactive area
				if((pix[0]>=0.)&(pix[0]<432.)&(pix[1]>=0.)&(pix[1]<128.0)&(sdotvec<0.)) {
					mask[jj*nlon+ii] = 1;
					break;
				}
			}
		}
		fflush(stdout);
	}
	printf("\n");

	return mask;
}

void process(double et, int band, double *scloc, double *lat, double *lon, double *phase, double *inc, double *emission) {
	double lt, trgepoch, disti, loni, lati, plate_scale;
	int found;
	SpiceDouble pixvec[3], pos_jup[3], srfvec[3], spoint[3];
	SpiceDouble cam2jup[3][3];
	
	f1 = (FOCAL_LENGTH)/(PIXEL_SIZE);

	plate_scale = (1./FOCAL_LENGTH)/(PIXEL_SIZE);

	// get the spacecraft position first
	spkpos_c("JUNO", et, "IAU_JUPITER", "CN+S", "JUPITER", scloc, &lt);
	
	if (band==M_BAND) {
		pxform_c("JUNO_JIRAM_I_MBAND", "IAU_JUPITER", et, cam2jup);
	} else {
		pxform_c("JUNO_JIRAM_I_LBAND", "IAU_JUPITER", et, cam2jup);
	}

	for(int j=0; j<FRAME_HEIGHT; j++) {
		for(int i=0; i<FRAME_WIDTH; i++) {
			double xx, yy, rr;
			xx = (double) i;
			xx = CENTER_X - xx;
			
			yy = (double) j;
			yy = CENTER_Y - yy;
			
			rr = sqrt(xx*xx + yy*yy + f1*f1);


			pixvec[0] = yy/rr; pixvec[1] = xx/rr; pixvec[2] = f1/rr;

			mxv_c(cam2jup, pixvec, pos_jup);

			sincpt_c("Ellipsoid", "JUPITER", et, "IAU_JUPITER", "CN+S", \
					"JUNO", "IAU_JUPITER", pos_jup, spoint, &trgepoch, srfvec, &found);

			if(found) {
				reclat_c(spoint, &disti, &loni, &lati);

				lati = lati*180./M_PI;
				loni = loni*180./M_PI;

				ilumin_c("Ellipsoid", "JUPITER", et, "IAU_JUPITER", "CN+S", "JUNO", spoint, \
							&trgepoch, srfvec, 
							&phase[j*FRAME_WIDTH+i], 
							&inc[j*FRAME_WIDTH+i], 
							&emission[j*FRAME_WIDTH+i]);
			} else {
				loni     = -1000.;
				lati     = -1000.;
			}

			lat[j*FRAME_WIDTH + i] = lati;
			lon[j*FRAME_WIDTH + i] = loni;
		}
	}

}

void furnish(char *file) {
	furnsh_c(file);
}

