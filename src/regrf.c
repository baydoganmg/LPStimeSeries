/*******************************************************************
   Copyright (C) 2001-2012 Leo Breiman, Adele Cutler and Merck & Co., Inc.

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
*******************************************************************/

#include <R.h>
#include "rf.h"
#include "rng.h"

#ifdef _OPENMP
#include <omp.h>
#endif

void regRF_time_series(double *x, double *seglength, int *isRand, int *tardiff, int *segdiff,
		int *xdim, int *sampsize, int *nthsize, int *nrnodes, int *nTree, int *mtry,
		int *cat, int *jprint, int *oobpred, int *treeErrors, int *target, int *targetType, int *treeSize,
		int *nodedepth, int *nodestatus, int *splitType, int *lDaughter, int *rDaughter,
		double *avnode, int *mbest, double *upper, int *keepf, int *replace, double *oobpredictions,
		double *ooberrors, int *inbag, double *errorTree, int *nthreads) {

    double  *xb, *predictions, xrand, temp;
    int *in, *targetcount, k, m, n, j, nsample, idx, mdim;
    int segmentlength, keepF, keepInbag, maxdepth, oobcount=0,oobcount2=0;
    int nthreads_use;

    nsample = xdim[0];    		//number of series
    mdim = xdim[1];	      		//length series
    keepF = keepf[0];     		//keep forest
    keepInbag = keepf[1]; 		//keep inbag data
    maxdepth=log2(*nrnodes+1);	//maximum depth

    nthreads_use = *nthreads;
#ifdef _OPENMP
    if (nthreads_use < 1) nthreads_use = omp_get_max_threads();
    {
        int max_threads = omp_get_max_threads();
        if (nthreads_use > max_threads) nthreads_use = max_threads;
    }
#else
    nthreads_use = 1;
#endif

    GetRNGstate();

    /* Generate per-tree RNG seeds from R's RNG (sequential) */
    rng_state_t *tree_rngs = (rng_state_t *) R_Calloc(*nTree, rng_state_t);
    for (j = 0; j < *nTree; j++) {
        uint64_t seed = (uint64_t)(unif_rand() * 4294967296.0);
        seed = (seed << 32) | (uint64_t)(unif_rand() * 4294967296.0);
        rng_seed(&tree_rngs[j], seed);
    }

	if (*replace) { /* sample */
		xb = (double *) R_Calloc(mdim*(*sampsize), double); //inbag x info
		in = (int *) R_Calloc(nsample, int);
		targetcount = (int *) R_Calloc(mdim*nsample, int);
		zeroInt(targetcount, mdim*nsample);
		predictions = (double *) R_Calloc(mdim*nsample, double);
	} else if (*treeErrors) {
		targetcount = (int *) R_Calloc(mdim, int);
		predictions = (double *) R_Calloc(mdim*nsample, double);
	}

    /*************************************
     * Start the loop over trees.
     *************************************/

    if (*replace) {
        /* OOB path: keep sequential due to shared buffer updates */
        for (j = 0; j < *nTree; ++j) {
            idx = keepF ? j * *nrnodes : 0;
            zeroInt(in, nsample);

            for (n = 0; n < *sampsize; ++n) {
                xrand = rng_uniform(&tree_rngs[j]);
                k = xrand * nsample;
                in[k] = 1;
                for(m = 0; m < mdim; ++m) {
                    xb[m + n * mdim] = x[m + k * mdim];
                }
            }
            if (keepInbag) {
                for (n = 0; n < nsample; ++n) inbag[n + j * nsample] = in[n];
            }
            /* grow the regression tree */
            regTree_time_series(xb, seglength + j, *tardiff, *segdiff, maxdepth, mdim, *sampsize, lDaughter + idx, rDaughter + idx,
                upper + idx, avnode + idx, nodedepth + idx, nodestatus + idx, splitType + idx, *nrnodes,
                treeSize + j, *nthsize, *mtry, mbest + idx, target + j, targetType + j, cat, *isRand,
                &tree_rngs[j]);

            if (*oobpred) {
                segmentlength=(int) (mdim*seglength[j]);
                oobcount=0;
                errorTree[j]=0;
                ooberrors[j]=0;
                zeroDouble(predictions, mdim*nsample);
                oobcount2=0;
                for (n = 0; n < nsample; ++n) {
                    if(in[n]==0){
                        predict_time_series(x + n*mdim, segmentlength, 1, mdim, lDaughter + idx, rDaughter + idx,
                            nodedepth + idx, nodestatus + idx, upper + idx, mbest + idx, splitType + idx,
                            avnode + idx, maxdepth, target[j], predictions + n*mdim, targetcount + n*mdim, 0);

                        for (m = 0; m < mdim; ++m) {
                            if(m>=(target[j]-1)&&m<(target[j]+segmentlength-1)){
                                errorTree[j]=errorTree[j]+pow(x[m+n*mdim]-predictions[m+n*mdim],2);
                                oobpredictions[m+n*mdim]=oobpredictions[m+n*mdim]+predictions[m+n*mdim];
                                oobcount2++;
                            }
                        }
                    }
                    for (m = 0; m < mdim; ++m) {
                        if(targetcount[m+n*mdim]>0){
                            temp=oobpredictions[m+n*mdim]/targetcount[m+n*mdim];
                            ooberrors[j]=ooberrors[j]+pow(x[m+n*mdim]-temp,2);
                            oobcount++;
                        }
                    }
                }
                ooberrors[j]=ooberrors[j]/oobcount;
                errorTree[j]=errorTree[j]/oobcount2;
            }

            if(*jprint>0&&(j+1)%(*jprint)==0)
                Rprintf("Tree %d over\n",j+1);
        }
    } else {
        /* No replacement: parallelize tree building */
#ifdef _OPENMP
        #pragma omp parallel for num_threads(nthreads_use) schedule(dynamic) private(idx)
#endif
        for (j = 0; j < *nTree; ++j) {
            idx = keepF ? j * *nrnodes : 0;

            regTree_time_series(x, seglength + j, *tardiff, *segdiff, maxdepth, mdim, nsample, lDaughter + idx, rDaughter + idx,
                upper + idx, avnode + idx, nodedepth + idx, nodestatus + idx, splitType + idx, *nrnodes,
                treeSize + j, *nthsize, *mtry, mbest + idx, target + j, targetType + j, cat, *isRand,
                &tree_rngs[j]);
        }

        /* Tree errors: sequential since it uses shared buffers */
        if (*treeErrors) {
            for (j = 0; j < *nTree; ++j) {
                idx = keepF ? j * *nrnodes : 0;
                if(*isRand==0) {
                    zeroInt(targetcount, mdim);
                    zeroDouble(predictions, mdim*nsample);
                    segmentlength=(int) (mdim*seglength[j]);
                    predict_time_series(x, segmentlength, nsample, mdim, lDaughter + idx, rDaughter + idx,
                            nodedepth + idx, nodestatus + idx, upper + idx, mbest + idx, splitType + idx,
                            avnode + idx, maxdepth, target[j], predictions, targetcount, 0);

                    oobcount=0;
                    for (m = 0; m < mdim; ++m) {
                        if(targetcount[m]>0){
                            for (n = 0; n < nsample; ++n) {
                                errorTree[j]=errorTree[j]+pow(x[m+n*mdim]-predictions[m+n*mdim],2);
                                oobcount++;
                            }
                        }
                    }
                    errorTree[j]=errorTree[j]/oobcount;
                }
            }
        }
    }

    PutRNGstate();
    /* end of tree iterations=======================================*/

	if (*replace) { /* free memory */
	   if (*oobpred) {
		   for (n = 0; n < nsample; ++n) {
			   for (m = 0; m < mdim; ++m) {
					if(targetcount[m+n*mdim]>0){
						oobpredictions[m+n*mdim]=oobpredictions[m+n*mdim]/targetcount[m+n*mdim];
					} else {
						oobpredictions[m+n*mdim]=-999;
					}
			   }
		   }
	   }
       R_Free(xb);
       R_Free(in);
       R_Free(targetcount);
       R_Free(predictions);
	}
	else if (*treeErrors) {
	   R_Free(predictions);
       R_Free(targetcount);
	}

    R_Free(tree_rngs);
}

void regForest_similarity(double *x, double *y, int *n, int *ny,
			double *seglength, int *mdim, int *ntree, int *usedtrees,
			int *lDaughter, int *rDaughter, int *nodestatus, int *nodedepth,
			int *nrnodes, double *xsplit, int *mbest, int *splitType,
			int *treeSize, int *maxdepth, int *simType, int *similarity,
			int *nthreads) {

    int i, j, k, m, idx1, segmentlength;
    int nthreads_use;

    nthreads_use = *nthreads;
#ifdef _OPENMP
    if (nthreads_use < 1) nthreads_use = omp_get_max_threads();
    {
        int max_threads = omp_get_max_threads();
        if (nthreads_use > max_threads) nthreads_use = max_threads;
    }
#else
    nthreads_use = 1;
#endif

    zeroInt(similarity, (*ny) * (*n));

    /* Allocate per-thread working buffers */
    int totx = (*n) * (*nrnodes);
    int totxtst = (*ny) * (*nrnodes);

    int **t_noderef = (int **) R_Calloc(nthreads_use, int *);
    int **t_nodetest = (int **) R_Calloc(nthreads_use, int *);
    int **t_tempnodestatus = (int **) R_Calloc(nthreads_use, int *);
    int **t_similarity = (int **) R_Calloc(nthreads_use, int *);

    for (i = 0; i < nthreads_use; i++) {
        t_noderef[i] = (int *) R_Calloc(totx, int);
        t_nodetest[i] = (int *) R_Calloc(totxtst, int);
        t_tempnodestatus[i] = (int *) R_Calloc(*nrnodes, int);
        t_similarity[i] = (int *) R_Calloc((*ny) * (*n), int);
        zeroInt(t_similarity[i], (*ny) * (*n));
    }

    /* Count used trees and build index */
    int nused = 0;
    for (i = 0; i < *ntree; i++) {
        if (usedtrees[i] == 1) nused++;
    }
    int *used_idx = (int *) R_Calloc(nused, int);
    int *used_offset = (int *) R_Calloc(nused, int);
    j = 0;
    idx1 = 0;
    for (i = 0; i < *ntree; i++) {
        if (usedtrees[i] == 1) {
            used_idx[j] = i;
            used_offset[j] = idx1;
            j++;
        }
        idx1 += *nrnodes;
    }

#ifdef _OPENMP
    #pragma omp parallel for num_threads(nthreads_use) schedule(dynamic) private(i, j, k, m, idx1, segmentlength)
#endif
    for (i = 0; i < nused; i++) {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        idx1 = used_offset[i];
        segmentlength = (int)(*mdim * seglength[used_idx[i]]);

        zeroInt(t_noderef[tid], totx);
        zeroInt(t_nodetest[tid], totxtst);
        zeroInt(t_tempnodestatus[tid], *nrnodes);

        // based on the maxdepth setting identify terminal nodes
        for (k = 0; k < *nrnodes; k++) {
            if(nodedepth[idx1+k]==*maxdepth) {
                t_tempnodestatus[tid][k]=NODE_TERMINAL;
            }
            else if(nodestatus[idx1+k]==NODE_TERMINAL){
                t_tempnodestatus[tid][k]=NODE_TERMINAL;
            }
        }

        predictRepresentation_time_series(x, segmentlength, *n, *mdim, lDaughter + idx1, rDaughter + idx1,
                nodedepth + idx1, t_tempnodestatus[tid], xsplit + idx1, mbest + idx1, splitType + idx1,
                t_noderef[tid], *maxdepth);
        predictRepresentation_time_series(y, segmentlength, *ny, *mdim, lDaughter + idx1, rDaughter + idx1,
                nodedepth + idx1, t_tempnodestatus[tid], xsplit + idx1, mbest + idx1, splitType + idx1,
                t_nodetest[tid], *maxdepth);

        for (k = 0; k < *nrnodes; k++) {
            if(t_tempnodestatus[tid][k]==NODE_TERMINAL){
                for(j = 0; j < (*ny); j++){
                    for(m = 0; m < (*n); m++){
                        if(*simType==0){
                            t_similarity[tid][j+(*ny)*m] += abs(t_noderef[tid][(*n)*k+m]-t_nodetest[tid][(*ny)*k+j]);
                        } else {
                            if(t_noderef[tid][(*n)*k+m]>t_nodetest[tid][(*ny)*k+j]) {
                                t_similarity[tid][j+(*ny)*m] += t_nodetest[tid][(*ny)*k+j];
                            } else {
                                t_similarity[tid][j+(*ny)*m] += t_noderef[tid][(*n)*k+m];
                            }
                        }
                    }
                }
            }
        }
    }

    /* Reduce per-thread similarity into final result */
    for (i = 0; i < nthreads_use; i++) {
        for (j = 0; j < (*ny) * (*n); j++) {
            similarity[j] += t_similarity[i][j];
        }
    }

    /* Free per-thread buffers */
    for (i = 0; i < nthreads_use; i++) {
        R_Free(t_noderef[i]);
        R_Free(t_nodetest[i]);
        R_Free(t_tempnodestatus[i]);
        R_Free(t_similarity[i]);
    }
    R_Free(t_noderef);
    R_Free(t_nodetest);
    R_Free(t_tempnodestatus);
    R_Free(t_similarity);
    R_Free(used_idx);
    R_Free(used_offset);
}

void regForest_represent(double *x, int *n, int *whichtree,
			double *seglength, int *mdim, int *ntree, int *usedtrees,
			int *lDaughter, int *rDaughter, int *nodestatus, int *nodedepth,
			int *nrnodes, double *xsplit, int *mbest, int *splitType,
			int *treeSize, int *maxdepth, int *representation, int *repLength) {

	int i, k, m, idx1,nodecount, termnodecount, tempidx;
    int *noderef, *tempnodestatus, segmentlength, totx;

	totx=(*n)*(*nrnodes);
	idx1 = 0;
	// compute the size of the representation  by computing
	// the total number of terminal nodes (termnodecount) over trees
	tempidx = idx1;
	termnodecount = 0;
	for (i = 0; i < *ntree; i++) {
		if(usedtrees[i]==1) {
			for (k = 0; k < *nrnodes; k++) {
				if(nodedepth[tempidx+k]==*maxdepth) {
					termnodecount++;
				}
				else if(nodedepth[tempidx+k]<=*maxdepth&&nodestatus[tempidx+k]==NODE_TERMINAL){
					termnodecount++;
				}
			}
		}
		tempidx += *nrnodes;
	}
	*repLength = termnodecount;


    noderef = (int *) R_Calloc(totx, int);
    tempnodestatus = (int *) R_Calloc(*nrnodes, int);

	nodecount=0;
	for (i = 0; i < *ntree; i++) {
		if(usedtrees[i]==1) {
			segmentlength=(int) (*mdim*seglength[i]);
			zeroInt(noderef, totx);
			zeroInt(tempnodestatus, *nrnodes);

			// based on the maxdepth setting identify terminal nodes
			for (k = 0; k < *nrnodes; k++) {
				if(nodedepth[idx1+k]==*maxdepth) {
					tempnodestatus[k]=NODE_TERMINAL;
				}
				else if(nodedepth[idx1+k]<=*maxdepth&&nodestatus[idx1+k]==NODE_TERMINAL){
					tempnodestatus[k]=NODE_TERMINAL;
				}
			}
			predictRepresentation_time_series(x, segmentlength, *n, *mdim, lDaughter + idx1, rDaughter + idx1,
					nodedepth + idx1, tempnodestatus, xsplit + idx1, mbest + idx1, splitType + idx1,
					noderef, *maxdepth);


			for (k = 0; k < *nrnodes; k++) {
				if(tempnodestatus[k]==NODE_TERMINAL){
					for(m = 0; m < (*n); m++){
						representation[termnodecount*m+nodecount]=noderef[(*n)*k+m];
					}
					nodecount++;
				}
		   }
	   }
	   idx1 += *nrnodes;
	}

	R_Free(noderef);
	R_Free(tempnodestatus);
}

void regForest_predict(double *x, int *n, int *whichtree,
			double *seglength, int *mdim, int *ntree, int *usedtrees,
			int *lDaughter, int *rDaughter, int *nodestatus, int *nodedepth,
			int *nrnodes, double *xsplit, int *mbest, int *splitType,
			double *nodepred, int *treeSize, int *target, int *maxdepth,
			double *prediction, int *targetcount) {

	int i, j, idx1;
    int segmentlength;

	idx1 = 0;
	zeroDouble(prediction,(*n)*(*mdim));
	zeroInt(targetcount,(*mdim));

	for (i = 0; i < *ntree; i++) {
		if(usedtrees[i]==1) {
			segmentlength=(int) (*mdim*seglength[i]);
			predict_time_series(x, segmentlength, *n, *mdim, lDaughter + idx1, rDaughter + idx1,
					nodedepth + idx1, nodestatus + idx1, xsplit + idx1, mbest + idx1, splitType + idx1,
					nodepred + idx1, *maxdepth, target[i], prediction, targetcount, 1);
	   }
	   idx1 += *nrnodes;
	}

	for (i = 0; i < *n; i++) {
		for (j = 0; j < *mdim; j++) {
			if(targetcount[j]>0){
				prediction[j+i*(*mdim)]=prediction[j+i*(*mdim)]/targetcount[j];
			} else {
				prediction[j+i*(*mdim)]=NA_indicator;
			}
		}
	}

}


void regForest_pattern(double *x, int *n, int *whichtree,
			int *whichterminal, double *seglength, int *mdim, int *ntree,
			int *lDaughter, int *rDaughter, int *nodestatus, int *nodedepth,
			int *nrnodes, double *xsplit, int *mbest, int *splitType,
			int *treeSize, int *maxdepth, int *target, int *targetType,
			double *predictpattern, double *targetpattern) {

	int i, j, k, m, idx1, termnodecount;
    int segmentlength, termid, curtree;

	curtree=*whichtree-1;
	idx1=*nrnodes*curtree;

	for (i = 0; i < *n; i++) {
		for (j = 0; j < *mdim; j++) {
			predictpattern[j+i*(*mdim)]=NA_indicator;
			targetpattern[j+i*(*mdim)]=NA_indicator;
		}
	}

	termnodecount = 0;
	for (k = 0; k < *nrnodes; k++) {
		if(nodedepth[idx1+k]==*maxdepth || nodestatus[idx1+k]==NODE_TERMINAL)
			termnodecount++;

		if(termnodecount==*whichterminal)
			break;
	}
	termid=k;
	segmentlength=(int) (*mdim*seglength[curtree]);
	for (j = 0; j < segmentlength ; j++) {
		for (i = 0; i < *n; i++) {
			k = 0;
			while (nodestatus[idx1+k] != NODE_TERMINAL && nodedepth[idx1+k] < *maxdepth) {
				m = mbest[idx1+k] - 1;
				if(splitType[idx1+k]==OBS_SERIES){
					if(m+j>(*mdim)-1){
						k = (x[m+j-(*mdim)+i*(*mdim)] <= xsplit[idx1+k]) ?
						lDaughter[idx1+k] - 1 : rDaughter[idx1+k] - 1;
					} else {
						k = (x[m+j+i*(*mdim)] <= xsplit[idx1+k]) ?
						lDaughter[idx1+k] - 1 : rDaughter[idx1+k] - 1;
					}
				}  else if(splitType[idx1+k]==DIFF_SERIES){
					if(m+j>(*mdim)-2){
						k = ((x[m+(j+2)-(*mdim)+i*(*mdim)]-x[m+j-(*mdim)+1+i*(*mdim)]) <= xsplit[idx1+k]) ?
						lDaughter[idx1+k] - 1 : rDaughter[idx1+k] - 1;
					} else {
						k = ((x[m+(j+1)+i*(*mdim)]-x[m+j+i*(*mdim)]) <= xsplit[idx1+k]) ?
						lDaughter[idx1+k] - 1 : rDaughter[idx1+k] - 1;
					}
				}
			}

			if(termid==k){
				m = target[curtree] - 1;
				if(m+j>(*mdim)-1){
					targetpattern[m+j-(*mdim)+i*(*mdim)] = x[m+j-(*mdim)+i*(*mdim)];
				} else {
					targetpattern[m+j+i*(*mdim)] = x[m+j+i*(*mdim)];
				}

				m = mbest[idx1] - 1;
				if(m+j>(*mdim)-1){
					predictpattern[m+j-(*mdim)+i*(*mdim)] = x[m+j-(*mdim)+i*(*mdim)];
				} else {
					predictpattern[m+j+i*(*mdim)] = x[m+j+i*(*mdim)];
				}

			}
		}
	}

}
