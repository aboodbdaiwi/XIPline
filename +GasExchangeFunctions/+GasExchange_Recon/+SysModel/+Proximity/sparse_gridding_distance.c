#include <math.h>
#include <stdio.h>
#include "mex.h"


/**************************************************************************
 *  SPARSE_GRIDDING_DISTANCE.C
 *
 *  Author: Scott Haile Robertson
 *  Website: www.ScottHaileRobertson.com
 *  Date: August 31, 2014
 *
 *  An N-Dimmensional convolution based gridding algorithm. Motivated by
 *  code written by Gary glover
 *  (http://www-mrsrl.stanford.edu/~brian/gridding/)
 *  and also from code by Nick Zwart.
 *
 *  For background reading, I suggest:
 *      1. A fast Sinc Function Gridding Algorithm for Fourier Inversion in
 *         Computer Tomography. O'Sullivan. 1985.
 *      2. Selection of a Convolution Function for Fourier Inversion using
 *         Gridding. Jackson et al. 1991.
 *      3. Rapid Gridding Reconstruction With a Minimal Oversampling Ratio.
 *         Beatty et al. 2005.
 *
 *  NOTE #1: You must compile with the -largeArrayDims option to enable
 *           the creation of sparse matrices
 *
 *  NOTE #2: compiling in debug mode (mex -largeArrayDims -g sparse_gridding_distance.c) turns on error
 *           checking of inputs. I recommend that you compile in debug mode
 *           to ensure that you are passing input arguments correctly, then
 *           recompile without debug mode once you have it wired correctly and
 *           enjoy the speed.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  This code
 *  is for research and academic purposes and is not intended for
 *  clinical use.
 *
 **************************************************************************/

/* In case c math.h libraries dont have min, max, and round defined, we'll define them. */
#ifndef min
double min(double a, double b) {
	return ((a < b) ? a : b);
}
#endif

#ifndef max
double max(double a, double b) {
	return ((a > b) ? a : b);
}
#endif

// #define DEBUG 1;

/*#ifndef round
 * int round(double num) {
 * return (num < 0.0) ? ceil(num - 0.5) : floor(num + 0.5);
 * }
 * #endif*/

/* Recursive function that loops through a bounded section of the output *
 * grid, convolving the ungridded point's data according to the provided *
 * convolution kernel, density compensation value, and the ungridded     *
 * point's value. The recursion allows for n-dimmensional data           *
 * reconstruction.                                                       */
void grid_point(double *sample_loc, unsigned int *idx_convert,
		double kernel_halfwidth_sqr,
		unsigned int ndims,
		int cur_dim, unsigned int *bounds,
		unsigned int *seed_pt, double kern_dist_sq,
		unsigned int *output_dims, unsigned int sample_index,
		unsigned int *n_nonsparse_entries,
		double *sparse_sample_indices,
		double *sparse_voxel_indices,
		double *sparse_distances,
		int force_dim){
	
	/* DEFINITIONS */
	unsigned int i;
	unsigned int j;
	unsigned int lower;
	unsigned int upper;
	unsigned int idx_;
	unsigned int kern_idx;
	double new_kern_dist_sq;
	
	lower = bounds[2 * cur_dim];
	upper = bounds[2 * cur_dim + 1];
	
	for(i = lower; i<=upper; i++){
		/* UPDATE SEED PT FOR RECURSIVE LOOPS */
		seed_pt[cur_dim]=i;
		
		if((cur_dim == force_dim) | (force_dim == -1)){
			new_kern_dist_sq = ((double)i-sample_loc[cur_dim]);
#ifdef DEBUG
			printf("\t\tSaving distance for dim %u = %f (adding %f to %f)\n",cur_dim+1, new_kern_dist_sq, new_kern_dist_sq*new_kern_dist_sq, kern_dist_sq);
#endif
			new_kern_dist_sq *= new_kern_dist_sq; // square it
			new_kern_dist_sq += kern_dist_sq;     // add to growing sum
		}else{
#ifdef DEBUG
			printf("\t\tWrong dimension dim %u, maintaining distance %f\n",cur_dim+1, kern_dist_sq);
#endif
			new_kern_dist_sq = kern_dist_sq;     // add to growing sum
		}
		
		if(cur_dim > 0){
#ifdef DEBUG
			printf("\tRecursing dim %u - gridding [",cur_dim+1, seed_pt[0],seed_pt[1],seed_pt[2]);
			for(j=0; j<ndims; j++){
				if(j < cur_dim){
					printf("%u:%u",bounds[2*j],bounds[2*j+1]);
				} else {
					printf("%u",seed_pt[j]);
				}
				if(j < (ndims-1)){
					printf(",");
				}
			}
			printf("]\n");
#endif
			
			// RECURSE THROUGH OTHER DIMMENSIONS
			grid_point(sample_loc, idx_convert, kernel_halfwidth_sqr, ndims, cur_dim-1,
					bounds, seed_pt, new_kern_dist_sq, output_dims, sample_index,
					n_nonsparse_entries,sparse_sample_indices,
					sparse_voxel_indices, sparse_distances, force_dim);
		} else {
			if(new_kern_dist_sq <= kernel_halfwidth_sqr){
// 				printf("FOUND SAMPLE %i, nkd=%f\n",*n_nonsparse_entries, new_kern_dist_sq);
				// CALCULATE INDEX FROM X,Y,Z COORDINATES
				idx_ = 0;//seed_pt[0];
				for(j=0; j<ndims; j++){
					idx_ = idx_ + (seed_pt[j])*idx_convert[j];
				}
				
				// Only write indices once per sample location
				sparse_sample_indices[*n_nonsparse_entries] = sample_index+1;
				sparse_voxel_indices[*n_nonsparse_entries] = (double)idx_+1;
				sparse_distances[*n_nonsparse_entries] = sqrt(new_kern_dist_sq);
				
				// Increment number of non spares entries
				(*n_nonsparse_entries)++;
			}
 #ifdef DEBUG
			else{
				printf("\tVoxel [%u, %u, %u] is too far (%f > %f) from sample point [%u, %u, %u] (index %u)!\n",
						seed_pt[0],seed_pt[1],seed_pt[2],new_kern_dist_sq,kernel_halfwidth_sqr,sample_loc[0],sample_loc[1],sample_loc[2], sample_index);
			}
#endif
		}

		
		/* RESET THIS DIMENSION */
		seed_pt[cur_dim]=lower;
	}
}


/* Performs convolution based gridding. Loops through a set of
 * n-dimensionalsample points and convolves them onto a grid. */
void sparse_gridding_distance(double *coords, double kernel_width,
		mwSize npts, unsigned int ndims,
		unsigned int *output_dims,
		double *nonsparse_sample_indices,
		double  *nonsparse_voxel_indices,
		double *nonsparse_distances,
		unsigned int *n_nonsparse_entries,
		unsigned int max_size, 
		int force_dim){
	
	/* DEFINITIONS */
	unsigned int *seed_pt;	// seed indices within subarray recursion loops
	unsigned int dim;       // dimmension loop index variable
	unsigned int p;         // point loop index variable
	unsigned int *bounds;   // array defining mininum and maximum boundaries of subarray
// 	unsigned int *sparse_voxel_indices;
	double *sample_loc;
	unsigned int *output_halfwidth;
	double index_multiplier;
	unsigned int *idx_convert;
	double kernel_halfwidth;
	double kernel_halfwidth_sqr;
	
	/* CALCULATE KERNEL HALFWIDTH AND OUTPUT_HALFWIDTH */
	kernel_halfwidth = kernel_width*0.5;
	kernel_halfwidth_sqr = kernel_halfwidth * kernel_halfwidth;
	output_halfwidth = calloc(ndims, sizeof(unsigned int));
	
	if(output_halfwidth == NULL){printf("Error allocating memory for output_halfwidth. Crashing... :)\n");}
	for(dim=0; dim<ndims; dim++){
		output_halfwidth[dim] = (unsigned int) ceil(((double)output_dims[dim])*0.5);
	}
	
	/* CALCULATE X, Y, Z TO INDEX CONVERSIONS */
	idx_convert = calloc(ndims, sizeof(unsigned int));
	if(idx_convert == NULL){printf("Error allocating memory for idx_convert. Crashing... :)\n");}
	for(dim=0; dim<ndims; dim++){
		idx_convert[dim] = 1;
		if(dim>0){
			for(p=0; p<dim; p++){
				idx_convert[dim] = idx_convert[dim] * output_dims[p];
			}
		}
	}
	
	/* ALLOCATE MEMORY */
	bounds = calloc(2 * ndims, sizeof(unsigned int));
	seed_pt = calloc(ndims, sizeof(unsigned int));
	sample_loc = calloc(ndims, sizeof(double));
	if(bounds                == NULL){printf("Error allocating memory for bounds. Crashing... :)\n");}
	if(seed_pt               == NULL){printf("Error allocating memory for seed_pt. Crashing... :)\n");}
	if(sample_loc            == NULL){printf("Error allocating memory for sample_loc. Crashing... :)\n");}

	/* LOOP THROUGH SAMPLE POINTS */
	for (p=0; p<npts; p++){
		
		/* CALCULATE SUBARRAY BOUNDARIES */
		for(dim=0; dim<ndims; dim++){
			sample_loc[dim] = (coords[ndims*p+dim]*(double)output_dims[dim]+(double)output_halfwidth[dim]); // Convert to voxel space indices - zero is upper left hand corner
			
			/* CALCULATE BOUNDS */
			bounds[2*dim] = (unsigned int) max(ceil(sample_loc[dim]-kernel_halfwidth),0);                   // Lower boundary
			bounds[2*dim+1] = (unsigned int) min(floor(sample_loc[dim]+kernel_halfwidth),output_dims[dim]-1);  // Upper boundary
			/* INITIALIZE RECURSIVE SEED POINT AS UPPER LEFT CORNER OF SUBARRAY */
			seed_pt[dim] = bounds[2*dim];
		}
		
#ifdef DEBUG
		printf("GRIDDING ");
		printf("Sample loc = [%f,%f,%f], ",sample_loc[0],sample_loc[1],sample_loc[2]);
		printf("Bounds     = [%u:%u,%u:%u,%u:%u]\n",bounds[0],bounds[1],bounds[2],bounds[3],bounds[4],bounds[5]);
#endif
		
		grid_point(sample_loc, idx_convert, kernel_halfwidth_sqr, ndims, 
				ndims-1, bounds,
				seed_pt, 0, output_dims, p, n_nonsparse_entries,
				nonsparse_sample_indices, nonsparse_voxel_indices, 
				nonsparse_distances, force_dim);
	}
	
	/* */
#ifdef DEBUG
	for (p=0; p<*n_nonsparse_entries; p++){
		printf("nsp=%u  vox=%f, sample %f, distance=[%f]\n",p,nonsparse_voxel_indices[p],
				nonsparse_sample_indices[p], nonsparse_distances[p*ndims]);
		
		//Need to sort by Vox Idx, then (sub sort) by Sample idx
		//ir = sample_idx - 1;
		//jc = first nonzero sample index(row) in vodel index (col)
	}
#endif
			// Return sparse_sample_indices, sparse_voxel_indices, sparse_distances
	
	/* Free up some memory */
	free(bounds);
	free(seed_pt);
	free(output_halfwidth);
}

