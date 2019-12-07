/*
  This is a Version 2.0 MPI + OpenMP implementation of LULESH

                 Copyright (c) 2010-2013.
      Lawrence Livermore National Security, LLC.
Produced at the Lawrence Livermore National Laboratory.
                  LLNL-CODE-461231
                All rights reserved.

This file is part of LULESH, Version 2.0.
Please also read this link -- http://www.opensource.org/licenses/index.php

//////////////
DIFFERENCES BETWEEN THIS VERSION (2.x) AND EARLIER VERSIONS:
* Addition of regions to make work more representative of multi-material codes
* Default size of each domain is 30^3 (27000 elem) instead of 45^3. This is
  more representative of our actual working set sizes
* Single source distribution supports pure serial, pure OpenMP, MPI-only, 
  and MPI+OpenMP
* Addition of ability to visualize the mesh using VisIt 
  https://wci.llnl.gov/codes/visit/download.html
* Various command line options (see ./lulesh2.0 -h)
 -q              : quiet mode - suppress stdout
 -i <iterations> : number of cycles to run
 -s <size>       : length of cube mesh along side
 -r <numregions> : Number of distinct regions (def: 11)
 -b <balance>    : Load balance between regions of a domain (def: 1)
 -c <cost>       : Extra cost of more expensive regions (def: 1)
 -f <filepieces> : Number of file parts for viz output (def: np/9)
 -p              : Print out progress
 -v              : Output viz file (requires compiling with -DVIZ_MESH
 -h              : This message

 printf("Usage: %s [opts]\n", execname);
      printf(" where [opts] is one or more of:\n");
      printf(" -q              : quiet mode - suppress all stdout\n");
      printf(" -i <iterations> : number of cycles to run\n");
      printf(" -s <size>       : length of cube mesh along side\n");
      printf(" -r <numregions> : Number of distinct regions (def: 11)\n");
      printf(" -b <balance>    : Load balance between regions of a domain (def: 1)\n");
      printf(" -c <cost>       : Extra cost of more expensive regions (def: 1)\n");
      printf(" -f <numfiles>   : Number of files to split viz dump into (def: (np+10)/9)\n");
      printf(" -p              : Print out progress\n");
      printf(" -v              : Output viz file (requires compiling with -DVIZ_MESH\n");
      printf(" -h              : This message\n");
      printf("\n\n");

*Notable changes in LULESH 2.0

* Split functionality into different files
lulesh.cc - where most (all?) of the timed functionality lies
lulesh-comm.cc - MPI functionality
lulesh-init.cc - Setup code
lulesh-viz.cc  - Support for visualization option
lulesh-util.cc - Non-timed functions
*
* The concept of "regions" was added, although every region is the same ideal
*    gas material, and the same sedov blast wave problem is still the only
*    problem its hardcoded to solve.
* Regions allow two things important to making this proxy app more representative:
*   Four of the LULESH routines are now performed on a region-by-region basis,
*     making the memory access patterns non-unit stride
*   Artificial load imbalances can be easily introduced that could impact
*     parallelization strategies.  
* The load balance flag changes region assignment.  Region number is raised to
*   the power entered for assignment probability.  Most likely regions changes
*   with MPI process id.
* The cost flag raises the cost of ~45% of the regions to evaluate EOS by the
*   entered multiple. The cost of 5% is 10x the entered multiple.
* MPI and OpenMP were added, and coalesced into a single version of the source
*   that can support serial builds, MPI-only, OpenMP-only, and MPI+OpenMP
* Added support to write plot files using "poor mans parallel I/O" when linked
*   with the silo library, which in turn can be read by VisIt.
* Enabled variable timestep calculation by default (courant condition), which
*   results in an additional reduction.
* Default domain (mesh) size reduced from 45^3 to 30^3
* Command line options to allow numerous test cases without needing to recompile
* Performance optimizations and code cleanup beyond LULESH 1.0
* Added a "Figure of Merit" calculation (elements solved per microsecond) and
*   output in support of using LULESH 2.0 for the 2017 CORAL procurement
*
* Possible Differences in Final Release (other changes possible)
*
* High Level mesh structure to allow data structure transformations
* Different default parameters
* Minor code performance changes and cleanup

TODO in future versions
* Add reader for (truly) unstructured meshes, probably serial only
* CMake based build system

//////////////

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

   * Redistributions of source code must retain the above copyright
     notice, this list of conditions and the disclaimer below.

   * Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the disclaimer (as noted below)
     in the documentation and/or other materials provided with the
     distribution.

   * Neither the name of the LLNS/LLNL nor the names of its contributors
     may be used to endorse or promote products derived from this software
     without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


Additional BSD Notice

1. This notice is required to be provided under our contract with the U.S.
   Department of Energy (DOE). This work was produced at Lawrence Livermore
   National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.

2. Neither the United States Government nor Lawrence Livermore National
   Security, LLC nor any of their employees, makes any warranty, express
   or implied, or assumes any liability or responsibility for the accuracy,
   completeness, or usefulness of any information, apparatus, product, or
   process disclosed, or represents that its use would not infringe
   privately-owned rights.

3. Also, reference herein to any specific commercial products, process, or
   services by trade name, trademark, manufacturer or otherwise does not
   necessarily constitute or imply its endorsement, recommendation, or
   favoring by the United States Government or Lawrence Livermore National
   Security, LLC. The views and opinions of authors expressed herein do not
   necessarily state or reflect those of the United States Government or
   Lawrence Livermore National Security, LLC, and shall not be used for
   advertising or product endorsement purposes.

*/

#include <climits>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <sys/time.h>
#include <iostream>
#include <unistd.h>
#include <assert.h>

#if _OPENMP
# include <omp.h>
#endif

/*laik headers*/
extern "C" {
#include "laik.h"
#include "laik-backend-mpi.h"
#include "laik-backend-tcp.h"
#include "../util/fault-tolerance-options.h"
}

//porting to laik
//#include "laik_vector.h" // includes lulesh.h
#include "lulesh.h"
#include "laik_partitioners.h"
#include "laik_vector.h"

/*********************************/
/* Data structure implementation */
/*********************************/

/* might want to add access methods so that memory can be */

/* better managed, as in luleshFT */

template<typename T>
T *Allocate(size_t size) {
    return static_cast<T *>(malloc(sizeof(T) * size));
}

template<typename T>
void Release(T **ptr) {
    if (*ptr != NULL) {
        free(*ptr);
        *ptr = NULL;
    }
}



/******************************************/

/* Work Routines */

double laik_reduce(double input, Int_t myRank, Laik_Data *laikTimer, Laik_Partitioning *all, Laik_Partitioning *master);

double getMPI_Wtime(cmdLineOpts* opts);

static inline
void TimeIncrement(Domain &domain, Laik_Data *laikDt, Laik_Partitioning *allPartitioning, Real_t &gnewdt) {
    Real_t targetdt = domain.stoptime() - domain.time();

    if ((domain.dtfixed() <= Real_t(0.0)) && (domain.cycle() != Int_t(0))) {
        Real_t ratio;
        Real_t olddt = domain.deltatime();

        /* This will require a reduction in parallel */
        gnewdt = Real_t(1.0e+20);
        Real_t newdt;
        if (domain.dtcourant() < gnewdt) {
            gnewdt = domain.dtcourant() / Real_t(2.0);
        }
        if (domain.dthydro() < gnewdt) {
            gnewdt = domain.dthydro() * Real_t(2.0) / Real_t(3.0);
        }

#if USE_MPI
        /*
      MPI_Allreduce(&gnewdt, &newdt, 1,
                    ((sizeof(Real_t) == 4) ? MPI_FLOAT : MPI_DOUBLE),
                    MPI_MIN, MPI_COMM_WORLD) ;

      */
        laik_switchto_partitioning(laikDt, allPartitioning, LAIK_DF_Preserve, LAIK_RO_Min);
        newdt = gnewdt;

#else
        newdt = gnewdt;
#endif

        ratio = newdt / olddt;
        if (ratio >= Real_t(1.0)) {
            if (ratio < domain.deltatimemultlb()) {
                newdt = olddt;
            } else if (ratio > domain.deltatimemultub()) {
                newdt = olddt * domain.deltatimemultub();
            }
        }

        if (newdt > domain.dtmax()) {
            newdt = domain.dtmax();
        }
        domain.deltatime() = newdt;
    }

    /* TRY TO PREVENT VERY SMALL SCALING ON THE NEXT CYCLE */
    if ((targetdt > domain.deltatime()) &&
        (targetdt < (Real_t(4.0) * domain.deltatime() / Real_t(3.0)))) {
        targetdt = Real_t(2.0) * domain.deltatime() / Real_t(3.0);
    }

    if (targetdt < domain.deltatime()) {
        domain.deltatime() = targetdt;
    }

    domain.time() += domain.deltatime();

    ++domain.cycle();
}

/******************************************/

static inline
void CollectDomainNodesToElemNodes(Domain &domain,
                                   const Index_t *elemToNode,
                                   Real_t elemX[8],
                                   Real_t elemY[8],
                                   Real_t elemZ[8]) {
    Index_t nd0i = elemToNode[0];
    Index_t nd1i = elemToNode[1];
    Index_t nd2i = elemToNode[2];
    Index_t nd3i = elemToNode[3];
    Index_t nd4i = elemToNode[4];
    Index_t nd5i = elemToNode[5];
    Index_t nd6i = elemToNode[6];
    Index_t nd7i = elemToNode[7];

    elemX[0] = domain.x(nd0i);
    elemX[1] = domain.x(nd1i);
    elemX[2] = domain.x(nd2i);
    elemX[3] = domain.x(nd3i);
    elemX[4] = domain.x(nd4i);
    elemX[5] = domain.x(nd5i);
    elemX[6] = domain.x(nd6i);
    elemX[7] = domain.x(nd7i);

    elemY[0] = domain.y(nd0i);
    elemY[1] = domain.y(nd1i);
    elemY[2] = domain.y(nd2i);
    elemY[3] = domain.y(nd3i);
    elemY[4] = domain.y(nd4i);
    elemY[5] = domain.y(nd5i);
    elemY[6] = domain.y(nd6i);
    elemY[7] = domain.y(nd7i);

    elemZ[0] = domain.z(nd0i);
    elemZ[1] = domain.z(nd1i);
    elemZ[2] = domain.z(nd2i);
    elemZ[3] = domain.z(nd3i);
    elemZ[4] = domain.z(nd4i);
    elemZ[5] = domain.z(nd5i);
    elemZ[6] = domain.z(nd6i);
    elemZ[7] = domain.z(nd7i);

}

/******************************************/

static inline
void InitStressTermsForElems(Domain &domain,
                             Real_t *sigxx, Real_t *sigyy, Real_t *sigzz,
                             Index_t numElem) {
    //
    // pull in the stresses appropriate to the hydro integration
    //

#pragma omp parallel for firstprivate(numElem)
    for (Index_t i = 0; i < numElem; ++i) {
        sigxx[i] = sigyy[i] = sigzz[i] = -domain.p(i) - domain.q(i);
    }
}

/******************************************/

static inline
void CalcElemShapeFunctionDerivatives(Real_t const x[],
                                      Real_t const y[],
                                      Real_t const z[],
                                      Real_t b[][8],
                                      Real_t *const volume) {
    const Real_t x0 = x[0];
    const Real_t x1 = x[1];
    const Real_t x2 = x[2];
    const Real_t x3 = x[3];
    const Real_t x4 = x[4];
    const Real_t x5 = x[5];
    const Real_t x6 = x[6];
    const Real_t x7 = x[7];

    const Real_t y0 = y[0];
    const Real_t y1 = y[1];
    const Real_t y2 = y[2];
    const Real_t y3 = y[3];
    const Real_t y4 = y[4];
    const Real_t y5 = y[5];
    const Real_t y6 = y[6];
    const Real_t y7 = y[7];

    const Real_t z0 = z[0];
    const Real_t z1 = z[1];
    const Real_t z2 = z[2];
    const Real_t z3 = z[3];
    const Real_t z4 = z[4];
    const Real_t z5 = z[5];
    const Real_t z6 = z[6];
    const Real_t z7 = z[7];

    Real_t fjxxi, fjxet, fjxze;
    Real_t fjyxi, fjyet, fjyze;
    Real_t fjzxi, fjzet, fjzze;
    Real_t cjxxi, cjxet, cjxze;
    Real_t cjyxi, cjyet, cjyze;
    Real_t cjzxi, cjzet, cjzze;

    fjxxi = Real_t(.125) * ((x6 - x0) + (x5 - x3) - (x7 - x1) - (x4 - x2));
    fjxet = Real_t(.125) * ((x6 - x0) - (x5 - x3) + (x7 - x1) - (x4 - x2));
    fjxze = Real_t(.125) * ((x6 - x0) + (x5 - x3) + (x7 - x1) + (x4 - x2));

    fjyxi = Real_t(.125) * ((y6 - y0) + (y5 - y3) - (y7 - y1) - (y4 - y2));
    fjyet = Real_t(.125) * ((y6 - y0) - (y5 - y3) + (y7 - y1) - (y4 - y2));
    fjyze = Real_t(.125) * ((y6 - y0) + (y5 - y3) + (y7 - y1) + (y4 - y2));

    fjzxi = Real_t(.125) * ((z6 - z0) + (z5 - z3) - (z7 - z1) - (z4 - z2));
    fjzet = Real_t(.125) * ((z6 - z0) - (z5 - z3) + (z7 - z1) - (z4 - z2));
    fjzze = Real_t(.125) * ((z6 - z0) + (z5 - z3) + (z7 - z1) + (z4 - z2));

    /* compute cofactors */
    cjxxi = (fjyet * fjzze) - (fjzet * fjyze);
    cjxet = -(fjyxi * fjzze) + (fjzxi * fjyze);
    cjxze = (fjyxi * fjzet) - (fjzxi * fjyet);

    cjyxi = -(fjxet * fjzze) + (fjzet * fjxze);
    cjyet = (fjxxi * fjzze) - (fjzxi * fjxze);
    cjyze = -(fjxxi * fjzet) + (fjzxi * fjxet);

    cjzxi = (fjxet * fjyze) - (fjyet * fjxze);
    cjzet = -(fjxxi * fjyze) + (fjyxi * fjxze);
    cjzze = (fjxxi * fjyet) - (fjyxi * fjxet);

    /* calculate partials :
     this need only be done for l = 0,1,2,3   since , by symmetry ,
     (6,7,4,5) = - (0,1,2,3) .
  */
    b[0][0] = -cjxxi - cjxet - cjxze;
    b[0][1] = cjxxi - cjxet - cjxze;
    b[0][2] = cjxxi + cjxet - cjxze;
    b[0][3] = -cjxxi + cjxet - cjxze;
    b[0][4] = -b[0][2];
    b[0][5] = -b[0][3];
    b[0][6] = -b[0][0];
    b[0][7] = -b[0][1];

    b[1][0] = -cjyxi - cjyet - cjyze;
    b[1][1] = cjyxi - cjyet - cjyze;
    b[1][2] = cjyxi + cjyet - cjyze;
    b[1][3] = -cjyxi + cjyet - cjyze;
    b[1][4] = -b[1][2];
    b[1][5] = -b[1][3];
    b[1][6] = -b[1][0];
    b[1][7] = -b[1][1];

    b[2][0] = -cjzxi - cjzet - cjzze;
    b[2][1] = cjzxi - cjzet - cjzze;
    b[2][2] = cjzxi + cjzet - cjzze;
    b[2][3] = -cjzxi + cjzet - cjzze;
    b[2][4] = -b[2][2];
    b[2][5] = -b[2][3];
    b[2][6] = -b[2][0];
    b[2][7] = -b[2][1];

    /* calculate jacobian determinant (volume) */
    //laik_log((Laik_LogLevel)2,"Debug %f\n", *volume);

    *volume = Real_t(8.) * (fjxet * cjxet + fjyet * cjyet + fjzet * cjzet);
}

/******************************************/

static inline
void SumElemFaceNormal(Real_t *normalX0, Real_t *normalY0, Real_t *normalZ0,
                       Real_t *normalX1, Real_t *normalY1, Real_t *normalZ1,
                       Real_t *normalX2, Real_t *normalY2, Real_t *normalZ2,
                       Real_t *normalX3, Real_t *normalY3, Real_t *normalZ3,
                       const Real_t x0, const Real_t y0, const Real_t z0,
                       const Real_t x1, const Real_t y1, const Real_t z1,
                       const Real_t x2, const Real_t y2, const Real_t z2,
                       const Real_t x3, const Real_t y3, const Real_t z3) {
    Real_t bisectX0 = Real_t(0.5) * (x3 + x2 - x1 - x0);
    Real_t bisectY0 = Real_t(0.5) * (y3 + y2 - y1 - y0);
    Real_t bisectZ0 = Real_t(0.5) * (z3 + z2 - z1 - z0);
    Real_t bisectX1 = Real_t(0.5) * (x2 + x1 - x3 - x0);
    Real_t bisectY1 = Real_t(0.5) * (y2 + y1 - y3 - y0);
    Real_t bisectZ1 = Real_t(0.5) * (z2 + z1 - z3 - z0);
    Real_t areaX = Real_t(0.25) * (bisectY0 * bisectZ1 - bisectZ0 * bisectY1);
    Real_t areaY = Real_t(0.25) * (bisectZ0 * bisectX1 - bisectX0 * bisectZ1);
    Real_t areaZ = Real_t(0.25) * (bisectX0 * bisectY1 - bisectY0 * bisectX1);

    *normalX0 += areaX;
    *normalX1 += areaX;
    *normalX2 += areaX;
    *normalX3 += areaX;

    *normalY0 += areaY;
    *normalY1 += areaY;
    *normalY2 += areaY;
    *normalY3 += areaY;

    *normalZ0 += areaZ;
    *normalZ1 += areaZ;
    *normalZ2 += areaZ;
    *normalZ3 += areaZ;
}

/******************************************/

static inline
void CalcElemNodeNormals(Real_t pfx[8],
                         Real_t pfy[8],
                         Real_t pfz[8],
                         const Real_t x[8],
                         const Real_t y[8],
                         const Real_t z[8]) {
    for (Index_t i = 0; i < 8; ++i) {
        pfx[i] = Real_t(0.0);
        pfy[i] = Real_t(0.0);
        pfz[i] = Real_t(0.0);
    }
    /* evaluate face one: nodes 0, 1, 2, 3 */
    SumElemFaceNormal(&pfx[0], &pfy[0], &pfz[0],
                      &pfx[1], &pfy[1], &pfz[1],
                      &pfx[2], &pfy[2], &pfz[2],
                      &pfx[3], &pfy[3], &pfz[3],
                      x[0], y[0], z[0], x[1], y[1], z[1],
                      x[2], y[2], z[2], x[3], y[3], z[3]);
    /* evaluate face two: nodes 0, 4, 5, 1 */
    SumElemFaceNormal(&pfx[0], &pfy[0], &pfz[0],
                      &pfx[4], &pfy[4], &pfz[4],
                      &pfx[5], &pfy[5], &pfz[5],
                      &pfx[1], &pfy[1], &pfz[1],
                      x[0], y[0], z[0], x[4], y[4], z[4],
                      x[5], y[5], z[5], x[1], y[1], z[1]);
    /* evaluate face three: nodes 1, 5, 6, 2 */
    SumElemFaceNormal(&pfx[1], &pfy[1], &pfz[1],
                      &pfx[5], &pfy[5], &pfz[5],
                      &pfx[6], &pfy[6], &pfz[6],
                      &pfx[2], &pfy[2], &pfz[2],
                      x[1], y[1], z[1], x[5], y[5], z[5],
                      x[6], y[6], z[6], x[2], y[2], z[2]);
    /* evaluate face four: nodes 2, 6, 7, 3 */
    SumElemFaceNormal(&pfx[2], &pfy[2], &pfz[2],
                      &pfx[6], &pfy[6], &pfz[6],
                      &pfx[7], &pfy[7], &pfz[7],
                      &pfx[3], &pfy[3], &pfz[3],
                      x[2], y[2], z[2], x[6], y[6], z[6],
                      x[7], y[7], z[7], x[3], y[3], z[3]);
    /* evaluate face five: nodes 3, 7, 4, 0 */
    SumElemFaceNormal(&pfx[3], &pfy[3], &pfz[3],
                      &pfx[7], &pfy[7], &pfz[7],
                      &pfx[4], &pfy[4], &pfz[4],
                      &pfx[0], &pfy[0], &pfz[0],
                      x[3], y[3], z[3], x[7], y[7], z[7],
                      x[4], y[4], z[4], x[0], y[0], z[0]);
    /* evaluate face six: nodes 4, 7, 6, 5 */
    SumElemFaceNormal(&pfx[4], &pfy[4], &pfz[4],
                      &pfx[7], &pfy[7], &pfz[7],
                      &pfx[6], &pfy[6], &pfz[6],
                      &pfx[5], &pfy[5], &pfz[5],
                      x[4], y[4], z[4], x[7], y[7], z[7],
                      x[6], y[6], z[6], x[5], y[5], z[5]);
}

/******************************************/

static inline
void SumElemStressesToNodeForces(const Real_t B[][8],
                                 const Real_t stress_xx,
                                 const Real_t stress_yy,
                                 const Real_t stress_zz,
                                 Real_t fx[], Real_t fy[], Real_t fz[]) {
    for (Index_t i = 0; i < 8; i++) {
        fx[i] = -(stress_xx * B[0][i]);
        fy[i] = -(stress_yy * B[1][i]);
        fz[i] = -(stress_zz * B[2][i]);
    }
}

/******************************************/

static inline
void IntegrateStressForElems(Domain &domain,
                             Real_t *sigxx, Real_t *sigyy, Real_t *sigzz,
                             Real_t *determ, Index_t numElem, Index_t numNode) {
#if _OPENMP
    Index_t numthreads = omp_get_max_threads();
#else
    Index_t numthreads = 1;
#endif

    Index_t numElem8 = numElem * 8;
    Real_t *fx_elem;
    Real_t *fy_elem;
    Real_t *fz_elem;
    Real_t fx_local[8];
    Real_t fy_local[8];
    Real_t fz_local[8];


    if (numthreads > 1) {
        fx_elem = Allocate<Real_t>(numElem8);
        fy_elem = Allocate<Real_t>(numElem8);
        fz_elem = Allocate<Real_t>(numElem8);
    }
    // loop over all elements

#pragma omp parallel for firstprivate(numElem)
    for (Index_t k = 0; k < numElem; ++k) {
        const Index_t *const elemToNode = domain.nodelist(k);
        Real_t B[3][8];// shape function derivatives
        Real_t x_local[8];
        Real_t y_local[8];
        Real_t z_local[8];

        // get nodal coordinates from global arrays and copy into local arrays.
        CollectDomainNodesToElemNodes(domain, elemToNode, x_local, y_local, z_local);

        // Volume calculation involves extra work for numerical consistency
        CalcElemShapeFunctionDerivatives(x_local, y_local, z_local,
                                         B, &determ[k]);

        CalcElemNodeNormals(B[0], B[1], B[2],
                            x_local, y_local, z_local);

        if (numthreads > 1) {
            // Eliminate thread writing conflicts at the nodes by giving
            // each element its own copy to write to
            SumElemStressesToNodeForces(B, sigxx[k], sigyy[k], sigzz[k],
                                        &fx_elem[k * 8],
                                        &fy_elem[k * 8],
                                        &fz_elem[k * 8]);
        } else {
            SumElemStressesToNodeForces(B, sigxx[k], sigyy[k], sigzz[k],
                                        fx_local, fy_local, fz_local);

            // copy nodal force contributions to global force arrray.
            for (Index_t lnode = 0; lnode < 8; ++lnode) {
                Index_t gnode = elemToNode[lnode];
                domain.fx(gnode) += fx_local[lnode];
                domain.fy(gnode) += fy_local[lnode];
                domain.fz(gnode) += fz_local[lnode];
            }
        }
    }

    if (numthreads > 1) {
        // If threaded, then we need to copy the data out of the temporary
        // arrays used above into the final forces field
#pragma omp parallel for firstprivate(numNode)
        for (Index_t gnode = 0; gnode < numNode; ++gnode) {
            Index_t count = domain.nodeElemCount(gnode);
            Index_t *cornerList = domain.nodeElemCornerList(gnode);
            Real_t fx_tmp = Real_t(0.0);
            Real_t fy_tmp = Real_t(0.0);
            Real_t fz_tmp = Real_t(0.0);
            for (Index_t i = 0; i < count; ++i) {
                Index_t elem = cornerList[i];
                fx_tmp += fx_elem[elem];
                fy_tmp += fy_elem[elem];
                fz_tmp += fz_elem[elem];
            }
            domain.fx(gnode) = fx_tmp;
            domain.fy(gnode) = fy_tmp;
            domain.fz(gnode) = fz_tmp;
        }
        Release(&fz_elem);
        Release(&fy_elem);
        Release(&fx_elem);
    }
}

/******************************************/

static inline
void VoluDer(const Real_t x0, const Real_t x1, const Real_t x2,
             const Real_t x3, const Real_t x4, const Real_t x5,
             const Real_t y0, const Real_t y1, const Real_t y2,
             const Real_t y3, const Real_t y4, const Real_t y5,
             const Real_t z0, const Real_t z1, const Real_t z2,
             const Real_t z3, const Real_t z4, const Real_t z5,
             Real_t *dvdx, Real_t *dvdy, Real_t *dvdz) {
    const Real_t twelfth = Real_t(1.0) / Real_t(12.0);

    *dvdx =
            (y1 + y2) * (z0 + z1) - (y0 + y1) * (z1 + z2) +
            (y0 + y4) * (z3 + z4) - (y3 + y4) * (z0 + z4) -
            (y2 + y5) * (z3 + z5) + (y3 + y5) * (z2 + z5);
    *dvdy =
            -(x1 + x2) * (z0 + z1) + (x0 + x1) * (z1 + z2) -
            (x0 + x4) * (z3 + z4) + (x3 + x4) * (z0 + z4) +
            (x2 + x5) * (z3 + z5) - (x3 + x5) * (z2 + z5);

    *dvdz =
            -(y1 + y2) * (x0 + x1) + (y0 + y1) * (x1 + x2) -
            (y0 + y4) * (x3 + x4) + (y3 + y4) * (x0 + x4) +
            (y2 + y5) * (x3 + x5) - (y3 + y5) * (x2 + x5);

    *dvdx *= twelfth;
    *dvdy *= twelfth;
    *dvdz *= twelfth;
}

/******************************************/

static inline
void CalcElemVolumeDerivative(Real_t dvdx[8],
                              Real_t dvdy[8],
                              Real_t dvdz[8],
                              const Real_t x[8],
                              const Real_t y[8],
                              const Real_t z[8]) {
    VoluDer(x[1], x[2], x[3], x[4], x[5], x[7],
            y[1], y[2], y[3], y[4], y[5], y[7],
            z[1], z[2], z[3], z[4], z[5], z[7],
            &dvdx[0], &dvdy[0], &dvdz[0]);
    VoluDer(x[0], x[1], x[2], x[7], x[4], x[6],
            y[0], y[1], y[2], y[7], y[4], y[6],
            z[0], z[1], z[2], z[7], z[4], z[6],
            &dvdx[3], &dvdy[3], &dvdz[3]);
    VoluDer(x[3], x[0], x[1], x[6], x[7], x[5],
            y[3], y[0], y[1], y[6], y[7], y[5],
            z[3], z[0], z[1], z[6], z[7], z[5],
            &dvdx[2], &dvdy[2], &dvdz[2]);
    VoluDer(x[2], x[3], x[0], x[5], x[6], x[4],
            y[2], y[3], y[0], y[5], y[6], y[4],
            z[2], z[3], z[0], z[5], z[6], z[4],
            &dvdx[1], &dvdy[1], &dvdz[1]);
    VoluDer(x[7], x[6], x[5], x[0], x[3], x[1],
            y[7], y[6], y[5], y[0], y[3], y[1],
            z[7], z[6], z[5], z[0], z[3], z[1],
            &dvdx[4], &dvdy[4], &dvdz[4]);
    VoluDer(x[4], x[7], x[6], x[1], x[0], x[2],
            y[4], y[7], y[6], y[1], y[0], y[2],
            z[4], z[7], z[6], z[1], z[0], z[2],
            &dvdx[5], &dvdy[5], &dvdz[5]);
    VoluDer(x[5], x[4], x[7], x[2], x[1], x[3],
            y[5], y[4], y[7], y[2], y[1], y[3],
            z[5], z[4], z[7], z[2], z[1], z[3],
            &dvdx[6], &dvdy[6], &dvdz[6]);
    VoluDer(x[6], x[5], x[4], x[3], x[2], x[0],
            y[6], y[5], y[4], y[3], y[2], y[0],
            z[6], z[5], z[4], z[3], z[2], z[0],
            &dvdx[7], &dvdy[7], &dvdz[7]);
}

/******************************************/

static inline
void CalcElemFBHourglassForce(Real_t *xd, Real_t *yd, Real_t *zd, Real_t hourgam[][4],
                              Real_t coefficient,
                              Real_t *hgfx, Real_t *hgfy, Real_t *hgfz) {
    Real_t hxx[4];
    for (Index_t i = 0; i < 4; i++) {
        hxx[i] = hourgam[0][i] * xd[0] + hourgam[1][i] * xd[1] +
                 hourgam[2][i] * xd[2] + hourgam[3][i] * xd[3] +
                 hourgam[4][i] * xd[4] + hourgam[5][i] * xd[5] +
                 hourgam[6][i] * xd[6] + hourgam[7][i] * xd[7];
    }
    for (Index_t i = 0; i < 8; i++) {
        hgfx[i] = coefficient *
                  (hourgam[i][0] * hxx[0] + hourgam[i][1] * hxx[1] +
                   hourgam[i][2] * hxx[2] + hourgam[i][3] * hxx[3]);
    }
    for (Index_t i = 0; i < 4; i++) {
        hxx[i] = hourgam[0][i] * yd[0] + hourgam[1][i] * yd[1] +
                 hourgam[2][i] * yd[2] + hourgam[3][i] * yd[3] +
                 hourgam[4][i] * yd[4] + hourgam[5][i] * yd[5] +
                 hourgam[6][i] * yd[6] + hourgam[7][i] * yd[7];
    }
    for (Index_t i = 0; i < 8; i++) {
        hgfy[i] = coefficient *
                  (hourgam[i][0] * hxx[0] + hourgam[i][1] * hxx[1] +
                   hourgam[i][2] * hxx[2] + hourgam[i][3] * hxx[3]);
    }
    for (Index_t i = 0; i < 4; i++) {
        hxx[i] = hourgam[0][i] * zd[0] + hourgam[1][i] * zd[1] +
                 hourgam[2][i] * zd[2] + hourgam[3][i] * zd[3] +
                 hourgam[4][i] * zd[4] + hourgam[5][i] * zd[5] +
                 hourgam[6][i] * zd[6] + hourgam[7][i] * zd[7];
    }
    for (Index_t i = 0; i < 8; i++) {
        hgfz[i] = coefficient *
                  (hourgam[i][0] * hxx[0] + hourgam[i][1] * hxx[1] +
                   hourgam[i][2] * hxx[2] + hourgam[i][3] * hxx[3]);
    }
}

/******************************************/

static inline
void CalcFBHourglassForceForElems(Domain &domain,
                                  Real_t *determ,
                                  Real_t *x8n, Real_t *y8n, Real_t *z8n,
                                  Real_t *dvdx, Real_t *dvdy, Real_t *dvdz,
                                  Real_t hourg, Index_t numElem,
                                  Index_t numNode) {

#if _OPENMP
    Index_t numthreads = omp_get_max_threads();
#else
    Index_t numthreads = 1;
#endif
    /*************************************************
    *
    *     FUNCTION: Calculates the Flanagan-Belytschko anti-hourglass
    *               force.
    *
    *************************************************/

    Index_t numElem8 = numElem * 8;

    Real_t *fx_elem;
    Real_t *fy_elem;
    Real_t *fz_elem;

    if (numthreads > 1) {
        fx_elem = Allocate<Real_t>(numElem8);
        fy_elem = Allocate<Real_t>(numElem8);
        fz_elem = Allocate<Real_t>(numElem8);
    }

    Real_t gamma[4][8];

    gamma[0][0] = Real_t(1.);
    gamma[0][1] = Real_t(1.);
    gamma[0][2] = Real_t(-1.);
    gamma[0][3] = Real_t(-1.);
    gamma[0][4] = Real_t(-1.);
    gamma[0][5] = Real_t(-1.);
    gamma[0][6] = Real_t(1.);
    gamma[0][7] = Real_t(1.);
    gamma[1][0] = Real_t(1.);
    gamma[1][1] = Real_t(-1.);
    gamma[1][2] = Real_t(-1.);
    gamma[1][3] = Real_t(1.);
    gamma[1][4] = Real_t(-1.);
    gamma[1][5] = Real_t(1.);
    gamma[1][6] = Real_t(1.);
    gamma[1][7] = Real_t(-1.);
    gamma[2][0] = Real_t(1.);
    gamma[2][1] = Real_t(-1.);
    gamma[2][2] = Real_t(1.);
    gamma[2][3] = Real_t(-1.);
    gamma[2][4] = Real_t(1.);
    gamma[2][5] = Real_t(-1.);
    gamma[2][6] = Real_t(1.);
    gamma[2][7] = Real_t(-1.);
    gamma[3][0] = Real_t(-1.);
    gamma[3][1] = Real_t(1.);
    gamma[3][2] = Real_t(-1.);
    gamma[3][3] = Real_t(1.);
    gamma[3][4] = Real_t(1.);
    gamma[3][5] = Real_t(-1.);
    gamma[3][6] = Real_t(1.);
    gamma[3][7] = Real_t(-1.);

/*************************************************/
/*    compute the hourglass modes */


#pragma omp parallel for firstprivate(numElem, hourg)
    for (Index_t i2 = 0; i2 < numElem; ++i2) {
        Real_t *fx_local, *fy_local, *fz_local;
        Real_t hgfx[8], hgfy[8], hgfz[8];

        Real_t coefficient;

        Real_t hourgam[8][4];
        Real_t xd1[8], yd1[8], zd1[8];

        const Index_t *elemToNode = domain.nodelist(i2);
        Index_t i3 = 8 * i2;
        Real_t volinv = Real_t(1.0) / determ[i2];
        Real_t ss1, mass1, volume13;
        for (Index_t i1 = 0; i1 < 4; ++i1) {

            Real_t hourmodx =
                    x8n[i3] * gamma[i1][0] + x8n[i3 + 1] * gamma[i1][1] +
                    x8n[i3 + 2] * gamma[i1][2] + x8n[i3 + 3] * gamma[i1][3] +
                    x8n[i3 + 4] * gamma[i1][4] + x8n[i3 + 5] * gamma[i1][5] +
                    x8n[i3 + 6] * gamma[i1][6] + x8n[i3 + 7] * gamma[i1][7];

            Real_t hourmody =
                    y8n[i3] * gamma[i1][0] + y8n[i3 + 1] * gamma[i1][1] +
                    y8n[i3 + 2] * gamma[i1][2] + y8n[i3 + 3] * gamma[i1][3] +
                    y8n[i3 + 4] * gamma[i1][4] + y8n[i3 + 5] * gamma[i1][5] +
                    y8n[i3 + 6] * gamma[i1][6] + y8n[i3 + 7] * gamma[i1][7];

            Real_t hourmodz =
                    z8n[i3] * gamma[i1][0] + z8n[i3 + 1] * gamma[i1][1] +
                    z8n[i3 + 2] * gamma[i1][2] + z8n[i3 + 3] * gamma[i1][3] +
                    z8n[i3 + 4] * gamma[i1][4] + z8n[i3 + 5] * gamma[i1][5] +
                    z8n[i3 + 6] * gamma[i1][6] + z8n[i3 + 7] * gamma[i1][7];

            hourgam[0][i1] = gamma[i1][0] - volinv * (dvdx[i3] * hourmodx +
                                                      dvdy[i3] * hourmody +
                                                      dvdz[i3] * hourmodz);

            hourgam[1][i1] = gamma[i1][1] - volinv * (dvdx[i3 + 1] * hourmodx +
                                                      dvdy[i3 + 1] * hourmody +
                                                      dvdz[i3 + 1] * hourmodz);

            hourgam[2][i1] = gamma[i1][2] - volinv * (dvdx[i3 + 2] * hourmodx +
                                                      dvdy[i3 + 2] * hourmody +
                                                      dvdz[i3 + 2] * hourmodz);

            hourgam[3][i1] = gamma[i1][3] - volinv * (dvdx[i3 + 3] * hourmodx +
                                                      dvdy[i3 + 3] * hourmody +
                                                      dvdz[i3 + 3] * hourmodz);

            hourgam[4][i1] = gamma[i1][4] - volinv * (dvdx[i3 + 4] * hourmodx +
                                                      dvdy[i3 + 4] * hourmody +
                                                      dvdz[i3 + 4] * hourmodz);

            hourgam[5][i1] = gamma[i1][5] - volinv * (dvdx[i3 + 5] * hourmodx +
                                                      dvdy[i3 + 5] * hourmody +
                                                      dvdz[i3 + 5] * hourmodz);

            hourgam[6][i1] = gamma[i1][6] - volinv * (dvdx[i3 + 6] * hourmodx +
                                                      dvdy[i3 + 6] * hourmody +
                                                      dvdz[i3 + 6] * hourmodz);

            hourgam[7][i1] = gamma[i1][7] - volinv * (dvdx[i3 + 7] * hourmodx +
                                                      dvdy[i3 + 7] * hourmody +
                                                      dvdz[i3 + 7] * hourmodz);

        }

        /* compute forces */
        /* store forces into h arrays (force arrays) */

        ss1 = domain.ss(i2);
        mass1 = domain.elemMass(i2);
        volume13 = CBRT(determ[i2]);

        Index_t n0si2 = elemToNode[0];
        Index_t n1si2 = elemToNode[1];
        Index_t n2si2 = elemToNode[2];
        Index_t n3si2 = elemToNode[3];
        Index_t n4si2 = elemToNode[4];
        Index_t n5si2 = elemToNode[5];
        Index_t n6si2 = elemToNode[6];
        Index_t n7si2 = elemToNode[7];

        xd1[0] = domain.xd(n0si2);
        xd1[1] = domain.xd(n1si2);
        xd1[2] = domain.xd(n2si2);
        xd1[3] = domain.xd(n3si2);
        xd1[4] = domain.xd(n4si2);
        xd1[5] = domain.xd(n5si2);
        xd1[6] = domain.xd(n6si2);
        xd1[7] = domain.xd(n7si2);

        yd1[0] = domain.yd(n0si2);
        yd1[1] = domain.yd(n1si2);
        yd1[2] = domain.yd(n2si2);
        yd1[3] = domain.yd(n3si2);
        yd1[4] = domain.yd(n4si2);
        yd1[5] = domain.yd(n5si2);
        yd1[6] = domain.yd(n6si2);
        yd1[7] = domain.yd(n7si2);

        zd1[0] = domain.zd(n0si2);
        zd1[1] = domain.zd(n1si2);
        zd1[2] = domain.zd(n2si2);
        zd1[3] = domain.zd(n3si2);
        zd1[4] = domain.zd(n4si2);
        zd1[5] = domain.zd(n5si2);
        zd1[6] = domain.zd(n6si2);
        zd1[7] = domain.zd(n7si2);

        coefficient = -hourg * Real_t(0.01) * ss1 * mass1 / volume13;

        CalcElemFBHourglassForce(xd1, yd1, zd1,
                                 hourgam,
                                 coefficient, hgfx, hgfy, hgfz);

        // With the threaded version, we write into local arrays per elem
        // so we don't have to worry about race conditions
        if (numthreads > 1) {
            fx_local = &fx_elem[i3];
            fx_local[0] = hgfx[0];
            fx_local[1] = hgfx[1];
            fx_local[2] = hgfx[2];
            fx_local[3] = hgfx[3];
            fx_local[4] = hgfx[4];
            fx_local[5] = hgfx[5];
            fx_local[6] = hgfx[6];
            fx_local[7] = hgfx[7];

            fy_local = &fy_elem[i3];
            fy_local[0] = hgfy[0];
            fy_local[1] = hgfy[1];
            fy_local[2] = hgfy[2];
            fy_local[3] = hgfy[3];
            fy_local[4] = hgfy[4];
            fy_local[5] = hgfy[5];
            fy_local[6] = hgfy[6];
            fy_local[7] = hgfy[7];

            fz_local = &fz_elem[i3];
            fz_local[0] = hgfz[0];
            fz_local[1] = hgfz[1];
            fz_local[2] = hgfz[2];
            fz_local[3] = hgfz[3];
            fz_local[4] = hgfz[4];
            fz_local[5] = hgfz[5];
            fz_local[6] = hgfz[6];
            fz_local[7] = hgfz[7];
        } else {
            domain.fx(n0si2) += hgfx[0];
            domain.fy(n0si2) += hgfy[0];
            domain.fz(n0si2) += hgfz[0];

            domain.fx(n1si2) += hgfx[1];
            domain.fy(n1si2) += hgfy[1];
            domain.fz(n1si2) += hgfz[1];

            domain.fx(n2si2) += hgfx[2];
            domain.fy(n2si2) += hgfy[2];
            domain.fz(n2si2) += hgfz[2];

            domain.fx(n3si2) += hgfx[3];
            domain.fy(n3si2) += hgfy[3];
            domain.fz(n3si2) += hgfz[3];

            domain.fx(n4si2) += hgfx[4];
            domain.fy(n4si2) += hgfy[4];
            domain.fz(n4si2) += hgfz[4];

            domain.fx(n5si2) += hgfx[5];
            domain.fy(n5si2) += hgfy[5];
            domain.fz(n5si2) += hgfz[5];

            domain.fx(n6si2) += hgfx[6];
            domain.fy(n6si2) += hgfy[6];
            domain.fz(n6si2) += hgfz[6];

            domain.fx(n7si2) += hgfx[7];
            domain.fy(n7si2) += hgfy[7];
            domain.fz(n7si2) += hgfz[7];
        }
    }

    if (numthreads > 1) {
        // Collect the data from the local arrays into the final force arrays
#pragma omp parallel for firstprivate(numNode)
        for (Index_t gnode = 0; gnode < numNode; ++gnode) {
            Index_t count = domain.nodeElemCount(gnode);
            Index_t *cornerList = domain.nodeElemCornerList(gnode);
            Real_t fx_tmp = Real_t(0.0);
            Real_t fy_tmp = Real_t(0.0);
            Real_t fz_tmp = Real_t(0.0);
            for (Index_t i = 0; i < count; ++i) {
                Index_t elem = cornerList[i];
                fx_tmp += fx_elem[elem];
                fy_tmp += fy_elem[elem];
                fz_tmp += fz_elem[elem];
            }
            domain.fx(gnode) += fx_tmp;
            domain.fy(gnode) += fy_tmp;
            domain.fz(gnode) += fz_tmp;
        }
        Release(&fz_elem);
        Release(&fy_elem);
        Release(&fx_elem);
    }
}

/******************************************/

static inline
void CalcHourglassControlForElems(Domain &domain,
                                  Real_t determ[], Real_t hgcoef) {
    Index_t numElem = domain.numElem();
    Index_t numElem8 = numElem * 8;
    Real_t *dvdx = Allocate<Real_t>(numElem8);
    Real_t *dvdy = Allocate<Real_t>(numElem8);
    Real_t *dvdz = Allocate<Real_t>(numElem8);
    Real_t *x8n = Allocate<Real_t>(numElem8);
    Real_t *y8n = Allocate<Real_t>(numElem8);
    Real_t *z8n = Allocate<Real_t>(numElem8);

    /* start loop over elements */
#pragma omp parallel for firstprivate(numElem)
    for (Index_t i = 0; i < numElem; ++i) {
        Real_t x1[8], y1[8], z1[8];
        Real_t pfx[8], pfy[8], pfz[8];

        Index_t *elemToNode = domain.nodelist(i);
        CollectDomainNodesToElemNodes(domain, elemToNode, x1, y1, z1);

        CalcElemVolumeDerivative(pfx, pfy, pfz, x1, y1, z1);

        /* load into temporary storage for FB Hour Glass control */
        for (Index_t ii = 0; ii < 8; ++ii) {
            Index_t jj = 8 * i + ii;

            dvdx[jj] = pfx[ii];
            dvdy[jj] = pfy[ii];
            dvdz[jj] = pfz[ii];
            x8n[jj] = x1[ii];
            y8n[jj] = y1[ii];
            z8n[jj] = z1[ii];
        }

        determ[i] = domain.volo(i) * domain.v(i);

        /* Do a check for negative volumes */
        if (domain.v(i) <= Real_t(0.0)) {
#if USE_MPI
            laik_log(LAIK_LL_Warning, "Debug 1\n");
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Abort(MPI_COMM_WORLD, VolumeError);
#else
            exit(VolumeError);
#endif
        }
    }

    if (hgcoef > Real_t(0.)) {
        CalcFBHourglassForceForElems(domain,
                                     determ, x8n, y8n, z8n, dvdx, dvdy, dvdz,
                                     hgcoef, numElem, domain.numNode());
    }

    Release(&z8n);
    Release(&y8n);
    Release(&x8n);
    Release(&dvdz);
    Release(&dvdy);
    Release(&dvdx);

    return;
}

/******************************************/

static inline
void CalcVolumeForceForElems(Domain &domain) {
    Index_t numElem = domain.numElem();
    if (numElem != 0) {
        Real_t hgcoef = domain.hgcoef();
        Real_t *sigxx = Allocate<Real_t>(numElem);
        Real_t *sigyy = Allocate<Real_t>(numElem);
        Real_t *sigzz = Allocate<Real_t>(numElem);
        Real_t *determ = Allocate<Real_t>(numElem);

        /* Sum contributions to total stress tensor */
        InitStressTermsForElems(domain, sigxx, sigyy, sigzz, numElem);

        // call elemlib stress integration loop to produce nodal forces from
        // material stresses.
        IntegrateStressForElems(domain,
                                sigxx, sigyy, sigzz, determ, numElem,
                                domain.numNode());

        /*
      if (domain.cycle() == 6){
          laik_log((Laik_LogLevel)2,"Debug: numNodes: %d\n",domain.numNode());
          for (int k = 0; k < numElem; ++k) {
              laik_log((Laik_LogLevel)2,"Debug: k: %d %f\n",k, determ[k]);
          }
      }
      */

        // check for negative element volume
#pragma omp parallel for firstprivate(numElem)
        for (Index_t k = 0; k < numElem; ++k) {
            //laik_log((Laik_LogLevel)2,"Debug: %f\n", determ[k]);
            if (determ[k] <= Real_t(0.0)) {
#if USE_MPI
                laik_log(LAIK_LL_Warning, "Debug 2\n");
                MPI_Barrier(MPI_COMM_WORLD);
                MPI_Abort(MPI_COMM_WORLD, VolumeError);
#else
                exit(VolumeError);
#endif
            }
        }

        CalcHourglassControlForElems(domain, determ, hgcoef);

        Release(&determ);
        Release(&sigzz);
        Release(&sigyy);
        Release(&sigxx);
    }
}

/******************************************/

static inline void CalcForceForNodes(Domain &domain) {
    Index_t numNode = domain.numNode();

#if USE_MPI
    //CommRecv(domain, MSG_COMM_SBN, 3,
    //         domain.sizeX() + 1, domain.sizeY() + 1, domain.sizeZ() + 1,
    //         true, false) ;

    domain.get_fx().switch_to_p1();
    domain.get_fy().switch_to_p1();
    domain.get_fz().switch_to_p1();
#endif

#pragma omp parallel for firstprivate(numNode)
    for (Index_t i = 0; i < numNode; ++i) {
        domain.fx(i) = Real_t(0.0);
        domain.fy(i) = Real_t(0.0);
        domain.fz(i) = Real_t(0.0);
    }

    /* Calcforce calls partial, force, hourq */
    CalcVolumeForceForElems(domain);


#if USE_MPI
    domain.get_fx().switch_to_p2();
    //domain.get_fx().switch_to_p1();
    domain.get_fy().switch_to_p2();
    //domain.get_fy().switch_to_p1();
    domain.get_fz().switch_to_p2();
    //domain.get_fz().switch_to_p1();
/*
  Domain_member fieldData[3] ;
  fieldData[0] = &Domain::fx ;
  fieldData[1] = &Domain::fy ;
  fieldData[2] = &Domain::fz ;
  
  CommSend(domain, MSG_COMM_SBN, 3, fieldData,
           domain.sizeX() + 1, domain.sizeY() + 1, domain.sizeZ() +  1,
           true, false) ;
  CommSBN(domain, 3, fieldData) ;
*/
#endif
}

/******************************************/

static inline
void CalcAccelerationForNodes(Domain &domain, Index_t numNode) {
    //laik_log((Laik_LogLevel)1, "in calculation of the acceleration:");

#pragma omp parallel for firstprivate(numNode)
    for (Index_t i = 0; i < numNode; ++i) {
        domain.xdd(i) = domain.fx(i) / domain.nodalMass(i);
        domain.ydd(i) = domain.fy(i) / domain.nodalMass(i);
        domain.zdd(i) = domain.fz(i) / domain.nodalMass(i);
        //laik_log((Laik_LogLevel)1, "%d i, %f", i, domain.nodalMass(i));
    }
}

/******************************************/

static inline
void ApplyAccelerationBoundaryConditionsForNodes(Domain &domain) {
    Index_t size = domain.sizeX();
    Index_t numNodeBC = (size + 1) * (size + 1);

#pragma omp parallel
    {
        if (!domain.symmXempty() != 0) {
#pragma omp for nowait firstprivate(numNodeBC)
            for (Index_t i = 0; i < numNodeBC; ++i)
                domain.xdd(domain.symmX(i)) = Real_t(0.0);
        }

        if (!domain.symmYempty() != 0) {
#pragma omp for nowait firstprivate(numNodeBC)
            for (Index_t i = 0; i < numNodeBC; ++i)
                domain.ydd(domain.symmY(i)) = Real_t(0.0);
        }

        if (!domain.symmZempty() != 0) {
#pragma omp for nowait firstprivate(numNodeBC)
            for (Index_t i = 0; i < numNodeBC; ++i)
                domain.zdd(domain.symmZ(i)) = Real_t(0.0);
        }
    }
}

/******************************************/

static inline
void CalcVelocityForNodes(Domain &domain, const Real_t dt, const Real_t u_cut,
                          Index_t numNode) {

#pragma omp parallel for firstprivate(numNode)
    for (Index_t i = 0; i < numNode; ++i) {
        Real_t xdtmp, ydtmp, zdtmp;

        xdtmp = domain.xd(i) + domain.xdd(i) * dt;
        if (FABS(xdtmp) < u_cut) xdtmp = Real_t(0.0);
        domain.xd(i) = xdtmp;

        ydtmp = domain.yd(i) + domain.ydd(i) * dt;
        if (FABS(ydtmp) < u_cut) ydtmp = Real_t(0.0);
        domain.yd(i) = ydtmp;

        zdtmp = domain.zd(i) + domain.zdd(i) * dt;
        if (FABS(zdtmp) < u_cut) zdtmp = Real_t(0.0);
        domain.zd(i) = zdtmp;
    }
}

/******************************************/

static inline
void CalcPositionForNodes(Domain &domain, const Real_t dt, Index_t numNode) {
#pragma omp parallel for firstprivate(numNode)
    for (Index_t i = 0; i < numNode; ++i) {
        domain.x(i) += domain.xd(i) * dt;
        domain.y(i) += domain.yd(i) * dt;
        domain.z(i) += domain.zd(i) * dt;
    }
}

/******************************************/

static inline
void LagrangeNodal(Domain &domain) {
#ifdef SEDOV_SYNC_POS_VEL_EARLY
    Domain_member fieldData[6] ;
#endif

    const Real_t delt = domain.deltatime();
    Real_t u_cut = domain.u_cut();
    /* time of boundary condition evaluation is beginning of step for force and
   * acceleration boundary conditions. */
    CalcForceForNodes(domain);

#if USE_MPI
#ifdef SEDOV_SYNC_POS_VEL_EARLY
                                                                                                                            //CommRecv(domain, MSG_SYNC_POS_VEL, 6,
   //         domain.sizeX() + 1, domain.sizeY() + 1, domain.sizeZ() + 1,
   //         false, false) ;
#endif
#endif

    CalcAccelerationForNodes(domain, domain.numNode());

    ApplyAccelerationBoundaryConditionsForNodes(domain);

    CalcVelocityForNodes(domain, delt, u_cut, domain.numNode());

    CalcPositionForNodes(domain, delt, domain.numNode());

#if USE_MPI
#ifdef SEDOV_SYNC_POS_VEL_EARLY
                                                                                                                            fieldData[0] = &Domain::x ;
  fieldData[1] = &Domain::y ;
  fieldData[2] = &Domain::z ;
  fieldData[3] = &Domain::xd ;
  fieldData[4] = &Domain::yd ;
  fieldData[5] = &Domain::zd ;

   //CommSend(domain, MSG_SYNC_POS_VEL, 6, fieldData,
   //         domain.sizeX() + 1, domain.sizeY() + 1, domain.sizeZ() + 1,
   //         false, false) ;
   //CommSyncPosVel(domain) ;
#endif
#endif

    //laik_log((Laik_LogLevel)2, "new iteration");
    //for (Index_t i=0; i<domain.numNode(); ++i) {
    //    laik_log((Laik_LogLevel)2, "%f", domain.x(i));
    //}

    return;
}

/******************************************/

static inline
Real_t CalcElemVolume(const Real_t x0, const Real_t x1,
                      const Real_t x2, const Real_t x3,
                      const Real_t x4, const Real_t x5,
                      const Real_t x6, const Real_t x7,
                      const Real_t y0, const Real_t y1,
                      const Real_t y2, const Real_t y3,
                      const Real_t y4, const Real_t y5,
                      const Real_t y6, const Real_t y7,
                      const Real_t z0, const Real_t z1,
                      const Real_t z2, const Real_t z3,
                      const Real_t z4, const Real_t z5,
                      const Real_t z6, const Real_t z7) {
    Real_t twelveth = Real_t(1.0) / Real_t(12.0);

    Real_t dx61 = x6 - x1;
    Real_t dy61 = y6 - y1;
    Real_t dz61 = z6 - z1;

    Real_t dx70 = x7 - x0;
    Real_t dy70 = y7 - y0;
    Real_t dz70 = z7 - z0;

    Real_t dx63 = x6 - x3;
    Real_t dy63 = y6 - y3;
    Real_t dz63 = z6 - z3;

    Real_t dx20 = x2 - x0;
    Real_t dy20 = y2 - y0;
    Real_t dz20 = z2 - z0;

    Real_t dx50 = x5 - x0;
    Real_t dy50 = y5 - y0;
    Real_t dz50 = z5 - z0;

    Real_t dx64 = x6 - x4;
    Real_t dy64 = y6 - y4;
    Real_t dz64 = z6 - z4;

    Real_t dx31 = x3 - x1;
    Real_t dy31 = y3 - y1;
    Real_t dz31 = z3 - z1;

    Real_t dx72 = x7 - x2;
    Real_t dy72 = y7 - y2;
    Real_t dz72 = z7 - z2;

    Real_t dx43 = x4 - x3;
    Real_t dy43 = y4 - y3;
    Real_t dz43 = z4 - z3;

    Real_t dx57 = x5 - x7;
    Real_t dy57 = y5 - y7;
    Real_t dz57 = z5 - z7;

    Real_t dx14 = x1 - x4;
    Real_t dy14 = y1 - y4;
    Real_t dz14 = z1 - z4;

    Real_t dx25 = x2 - x5;
    Real_t dy25 = y2 - y5;
    Real_t dz25 = z2 - z5;

#define TRIPLE_PRODUCT(x1, y1, z1, x2, y2, z2, x3, y3, z3) \
   ((x1)*((y2)*(z3) - (z2)*(y3)) + (x2)*((z1)*(y3) - (y1)*(z3)) + (x3)*((y1)*(z2) - (z1)*(y2)))

    Real_t volume =
            TRIPLE_PRODUCT(dx31 + dx72, dx63, dx20,
                           dy31 + dy72, dy63, dy20,
                           dz31 + dz72, dz63, dz20) +
            TRIPLE_PRODUCT(dx43 + dx57, dx64, dx70,
                           dy43 + dy57, dy64, dy70,
                           dz43 + dz57, dz64, dz70) +
            TRIPLE_PRODUCT(dx14 + dx25, dx61, dx50,
                           dy14 + dy25, dy61, dy50,
                           dz14 + dz25, dz61, dz50);

#undef TRIPLE_PRODUCT

    volume *= twelveth;

    return volume;
}

/******************************************/

//inline
Real_t CalcElemVolume(const Real_t x[8], const Real_t y[8], const Real_t z[8]) {
    return CalcElemVolume(x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7],
                          y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7],
                          z[0], z[1], z[2], z[3], z[4], z[5], z[6], z[7]);
}

/******************************************/

static inline
Real_t AreaFace(const Real_t x0, const Real_t x1,
                const Real_t x2, const Real_t x3,
                const Real_t y0, const Real_t y1,
                const Real_t y2, const Real_t y3,
                const Real_t z0, const Real_t z1,
                const Real_t z2, const Real_t z3) {
    Real_t fx = (x2 - x0) - (x3 - x1);
    Real_t fy = (y2 - y0) - (y3 - y1);
    Real_t fz = (z2 - z0) - (z3 - z1);
    Real_t gx = (x2 - x0) + (x3 - x1);
    Real_t gy = (y2 - y0) + (y3 - y1);
    Real_t gz = (z2 - z0) + (z3 - z1);
    Real_t area =
            (fx * fx + fy * fy + fz * fz) *
            (gx * gx + gy * gy + gz * gz) -
            (fx * gx + fy * gy + fz * gz) *
            (fx * gx + fy * gy + fz * gz);
    return area;
}

/******************************************/

static inline
Real_t CalcElemCharacteristicLength(const Real_t x[8],
                                    const Real_t y[8],
                                    const Real_t z[8],
                                    const Real_t volume) {
    Real_t a, charLength = Real_t(0.0);

    a = AreaFace(x[0], x[1], x[2], x[3],
                 y[0], y[1], y[2], y[3],
                 z[0], z[1], z[2], z[3]);
    charLength = std::max(a, charLength);

    a = AreaFace(x[4], x[5], x[6], x[7],
                 y[4], y[5], y[6], y[7],
                 z[4], z[5], z[6], z[7]);
    charLength = std::max(a, charLength);

    a = AreaFace(x[0], x[1], x[5], x[4],
                 y[0], y[1], y[5], y[4],
                 z[0], z[1], z[5], z[4]);
    charLength = std::max(a, charLength);

    a = AreaFace(x[1], x[2], x[6], x[5],
                 y[1], y[2], y[6], y[5],
                 z[1], z[2], z[6], z[5]);
    charLength = std::max(a, charLength);

    a = AreaFace(x[2], x[3], x[7], x[6],
                 y[2], y[3], y[7], y[6],
                 z[2], z[3], z[7], z[6]);
    charLength = std::max(a, charLength);

    a = AreaFace(x[3], x[0], x[4], x[7],
                 y[3], y[0], y[4], y[7],
                 z[3], z[0], z[4], z[7]);
    charLength = std::max(a, charLength);

    charLength = Real_t(4.0) * volume / SQRT(charLength);

    return charLength;
}

/******************************************/

static inline
void CalcElemVelocityGradient(const Real_t *const xvel,
                              const Real_t *const yvel,
                              const Real_t *const zvel,
                              const Real_t b[][8],
                              const Real_t detJ,
                              Real_t *const d) {
    const Real_t inv_detJ = Real_t(1.0) / detJ;
    Real_t dyddx, dxddy, dzddx, dxddz, dzddy, dyddz;
    const Real_t *const pfx = b[0];
    const Real_t *const pfy = b[1];
    const Real_t *const pfz = b[2];

    d[0] = inv_detJ * (pfx[0] * (xvel[0] - xvel[6])
                       + pfx[1] * (xvel[1] - xvel[7])
                       + pfx[2] * (xvel[2] - xvel[4])
                       + pfx[3] * (xvel[3] - xvel[5]));

    d[1] = inv_detJ * (pfy[0] * (yvel[0] - yvel[6])
                       + pfy[1] * (yvel[1] - yvel[7])
                       + pfy[2] * (yvel[2] - yvel[4])
                       + pfy[3] * (yvel[3] - yvel[5]));

    d[2] = inv_detJ * (pfz[0] * (zvel[0] - zvel[6])
                       + pfz[1] * (zvel[1] - zvel[7])
                       + pfz[2] * (zvel[2] - zvel[4])
                       + pfz[3] * (zvel[3] - zvel[5]));

    dyddx = inv_detJ * (pfx[0] * (yvel[0] - yvel[6])
                        + pfx[1] * (yvel[1] - yvel[7])
                        + pfx[2] * (yvel[2] - yvel[4])
                        + pfx[3] * (yvel[3] - yvel[5]));

    dxddy = inv_detJ * (pfy[0] * (xvel[0] - xvel[6])
                        + pfy[1] * (xvel[1] - xvel[7])
                        + pfy[2] * (xvel[2] - xvel[4])
                        + pfy[3] * (xvel[3] - xvel[5]));

    dzddx = inv_detJ * (pfx[0] * (zvel[0] - zvel[6])
                        + pfx[1] * (zvel[1] - zvel[7])
                        + pfx[2] * (zvel[2] - zvel[4])
                        + pfx[3] * (zvel[3] - zvel[5]));

    dxddz = inv_detJ * (pfz[0] * (xvel[0] - xvel[6])
                        + pfz[1] * (xvel[1] - xvel[7])
                        + pfz[2] * (xvel[2] - xvel[4])
                        + pfz[3] * (xvel[3] - xvel[5]));

    dzddy = inv_detJ * (pfy[0] * (zvel[0] - zvel[6])
                        + pfy[1] * (zvel[1] - zvel[7])
                        + pfy[2] * (zvel[2] - zvel[4])
                        + pfy[3] * (zvel[3] - zvel[5]));

    dyddz = inv_detJ * (pfz[0] * (yvel[0] - yvel[6])
                        + pfz[1] * (yvel[1] - yvel[7])
                        + pfz[2] * (yvel[2] - yvel[4])
                        + pfz[3] * (yvel[3] - yvel[5]));
    d[5] = Real_t(.5) * (dxddy + dyddx);
    d[4] = Real_t(.5) * (dxddz + dzddx);
    d[3] = Real_t(.5) * (dzddy + dyddz);
}

/******************************************/

//static inline
void CalcKinematicsForElems(Domain &domain, Real_t *vnew,
                            Real_t deltaTime, Index_t numElem) {

    // loop over all elements
#pragma omp parallel for firstprivate(numElem, deltaTime)
    for (Index_t k = 0; k < numElem; ++k) {
        Real_t B[3][8]; /** shape function derivatives */
        Real_t D[6];
        Real_t x_local[8];
        Real_t y_local[8];
        Real_t z_local[8];
        Real_t xd_local[8];
        Real_t yd_local[8];
        Real_t zd_local[8];
        Real_t detJ = Real_t(0.0);

        Real_t volume;
        Real_t relativeVolume;
        const Index_t *const elemToNode = domain.nodelist(k);

        // get nodal coordinates from global arrays and copy into local arrays.
        CollectDomainNodesToElemNodes(domain, elemToNode, x_local, y_local, z_local);

        // volume calculations
        volume = CalcElemVolume(x_local, y_local, z_local);
        relativeVolume = volume / domain.volo(k);
        vnew[k] = relativeVolume;
        domain.delv(k) = relativeVolume - domain.v(k);

        // set characteristic length
        domain.arealg(k) = CalcElemCharacteristicLength(x_local, y_local, z_local,
                                                        volume);

        // get nodal velocities from global array and copy into local arrays.
        for (Index_t lnode = 0; lnode < 8; ++lnode) {
            Index_t gnode = elemToNode[lnode];
            xd_local[lnode] = domain.xd(gnode);
            yd_local[lnode] = domain.yd(gnode);
            zd_local[lnode] = domain.zd(gnode);
        }

        Real_t dt2 = Real_t(0.5) * deltaTime;
        for (Index_t j = 0; j < 8; ++j) {
            x_local[j] -= dt2 * xd_local[j];
            y_local[j] -= dt2 * yd_local[j];
            z_local[j] -= dt2 * zd_local[j];
        }

        CalcElemShapeFunctionDerivatives(x_local, y_local, z_local,
                                         B, &detJ);

        CalcElemVelocityGradient(xd_local, yd_local, zd_local,
                                 B, detJ, D);

        // put velocity gradient quantities into their global arrays.
        domain.dxx(k) = D[0];
        domain.dyy(k) = D[1];
        domain.dzz(k) = D[2];
    }
}

/******************************************/

static inline
void CalcLagrangeElements(Domain &domain, Real_t *vnew) {

    Index_t numElem = domain.numElem();
    if (numElem > 0) {
        const Real_t deltatime = domain.deltatime();

        //domain.AllocateStrains(numElem);

        CalcKinematicsForElems(domain, vnew, deltatime, numElem);

        // element loop to do some stuff not included in the elemlib function.
#pragma omp parallel for firstprivate(numElem)
        for (Index_t k = 0; k < numElem; ++k) {
            // calc strain rate and apply as constraint (only done in FB element)
            Real_t vdov = domain.dxx(k) + domain.dyy(k) + domain.dzz(k);
            Real_t vdovthird = vdov / Real_t(3.0);

            // make the rate of deformation tensor deviatoric
            domain.vdov(k) = vdov;
            domain.dxx(k) -= vdovthird;
            domain.dyy(k) -= vdovthird;
            domain.dzz(k) -= vdovthird;

            // See if any volumes are negative, and take appropriate action.
            if (vnew[k] <= Real_t(0.0)) {
#if USE_MPI
                laik_log(LAIK_LL_Warning, "Debug 3\n");
                MPI_Barrier(MPI_COMM_WORLD);
                MPI_Abort(MPI_COMM_WORLD, VolumeError);
#else
                exit(VolumeError);
#endif
            }
        }
        //domain.DeallocateStrains();
    }

}

/******************************************/

static inline
void CalcMonotonicQGradientsForElems(Domain &domain, Real_t vnew[]) {
    Index_t numElem = domain.numElem();

#pragma omp parallel for firstprivate(numElem)
    for (Index_t i = 0; i < numElem; ++i) {
        const Real_t ptiny = Real_t(1.e-36);
        Real_t ax, ay, az;
        Real_t dxv, dyv, dzv;

        const Index_t *elemToNode = domain.nodelist(i);

        Index_t n0 = elemToNode[0];
        Index_t n1 = elemToNode[1];
        Index_t n2 = elemToNode[2];
        Index_t n3 = elemToNode[3];
        Index_t n4 = elemToNode[4];
        Index_t n5 = elemToNode[5];
        Index_t n6 = elemToNode[6];
        Index_t n7 = elemToNode[7];

        Real_t x0 = domain.x(n0);
        Real_t x1 = domain.x(n1);
        Real_t x2 = domain.x(n2);
        Real_t x3 = domain.x(n3);
        Real_t x4 = domain.x(n4);
        Real_t x5 = domain.x(n5);
        Real_t x6 = domain.x(n6);
        Real_t x7 = domain.x(n7);

        Real_t y0 = domain.y(n0);
        Real_t y1 = domain.y(n1);
        Real_t y2 = domain.y(n2);
        Real_t y3 = domain.y(n3);
        Real_t y4 = domain.y(n4);
        Real_t y5 = domain.y(n5);
        Real_t y6 = domain.y(n6);
        Real_t y7 = domain.y(n7);

        Real_t z0 = domain.z(n0);
        Real_t z1 = domain.z(n1);
        Real_t z2 = domain.z(n2);
        Real_t z3 = domain.z(n3);
        Real_t z4 = domain.z(n4);
        Real_t z5 = domain.z(n5);
        Real_t z6 = domain.z(n6);
        Real_t z7 = domain.z(n7);

        Real_t xv0 = domain.xd(n0);
        Real_t xv1 = domain.xd(n1);
        Real_t xv2 = domain.xd(n2);
        Real_t xv3 = domain.xd(n3);
        Real_t xv4 = domain.xd(n4);
        Real_t xv5 = domain.xd(n5);
        Real_t xv6 = domain.xd(n6);
        Real_t xv7 = domain.xd(n7);

        Real_t yv0 = domain.yd(n0);
        Real_t yv1 = domain.yd(n1);
        Real_t yv2 = domain.yd(n2);
        Real_t yv3 = domain.yd(n3);
        Real_t yv4 = domain.yd(n4);
        Real_t yv5 = domain.yd(n5);
        Real_t yv6 = domain.yd(n6);
        Real_t yv7 = domain.yd(n7);

        Real_t zv0 = domain.zd(n0);
        Real_t zv1 = domain.zd(n1);
        Real_t zv2 = domain.zd(n2);
        Real_t zv3 = domain.zd(n3);
        Real_t zv4 = domain.zd(n4);
        Real_t zv5 = domain.zd(n5);
        Real_t zv6 = domain.zd(n6);
        Real_t zv7 = domain.zd(n7);

        Real_t vol = domain.volo(i) * vnew[i];
        Real_t norm = Real_t(1.0) / (vol + ptiny);

        Real_t dxj = Real_t(-0.25) * ((x0 + x1 + x5 + x4) - (x3 + x2 + x6 + x7));
        Real_t dyj = Real_t(-0.25) * ((y0 + y1 + y5 + y4) - (y3 + y2 + y6 + y7));
        Real_t dzj = Real_t(-0.25) * ((z0 + z1 + z5 + z4) - (z3 + z2 + z6 + z7));

        Real_t dxi = Real_t(0.25) * ((x1 + x2 + x6 + x5) - (x0 + x3 + x7 + x4));
        Real_t dyi = Real_t(0.25) * ((y1 + y2 + y6 + y5) - (y0 + y3 + y7 + y4));
        Real_t dzi = Real_t(0.25) * ((z1 + z2 + z6 + z5) - (z0 + z3 + z7 + z4));

        Real_t dxk = Real_t(0.25) * ((x4 + x5 + x6 + x7) - (x0 + x1 + x2 + x3));
        Real_t dyk = Real_t(0.25) * ((y4 + y5 + y6 + y7) - (y0 + y1 + y2 + y3));
        Real_t dzk = Real_t(0.25) * ((z4 + z5 + z6 + z7) - (z0 + z1 + z2 + z3));

        /* find delvk and delxk ( i cross j ) */

        ax = dyi * dzj - dzi * dyj;
        ay = dzi * dxj - dxi * dzj;
        az = dxi * dyj - dyi * dxj;

        domain.delx_zeta(i) = vol / SQRT(ax * ax + ay * ay + az * az + ptiny);

        ax *= norm;
        ay *= norm;
        az *= norm;

        dxv = Real_t(0.25) * ((xv4 + xv5 + xv6 + xv7) - (xv0 + xv1 + xv2 + xv3));
        dyv = Real_t(0.25) * ((yv4 + yv5 + yv6 + yv7) - (yv0 + yv1 + yv2 + yv3));
        dzv = Real_t(0.25) * ((zv4 + zv5 + zv6 + zv7) - (zv0 + zv1 + zv2 + zv3));

        domain.delv_zeta(i) = ax * dxv + ay * dyv + az * dzv;

        /* find delxi and delvi ( j cross k ) */

        ax = dyj * dzk - dzj * dyk;
        ay = dzj * dxk - dxj * dzk;
        az = dxj * dyk - dyj * dxk;

        domain.delx_xi(i) = vol / SQRT(ax * ax + ay * ay + az * az + ptiny);

        ax *= norm;
        ay *= norm;
        az *= norm;

        dxv = Real_t(0.25) * ((xv1 + xv2 + xv6 + xv5) - (xv0 + xv3 + xv7 + xv4));
        dyv = Real_t(0.25) * ((yv1 + yv2 + yv6 + yv5) - (yv0 + yv3 + yv7 + yv4));
        dzv = Real_t(0.25) * ((zv1 + zv2 + zv6 + zv5) - (zv0 + zv3 + zv7 + zv4));

        domain.delv_xi(i) = ax * dxv + ay * dyv + az * dzv;

        /* find delxj and delvj ( k cross i ) */

        ax = dyk * dzi - dzk * dyi;
        ay = dzk * dxi - dxk * dzi;
        az = dxk * dyi - dyk * dxi;

        domain.delx_eta(i) = vol / SQRT(ax * ax + ay * ay + az * az + ptiny);

        ax *= norm;
        ay *= norm;
        az *= norm;

        dxv = Real_t(-0.25) * ((xv0 + xv1 + xv5 + xv4) - (xv3 + xv2 + xv6 + xv7));
        dyv = Real_t(-0.25) * ((yv0 + yv1 + yv5 + yv4) - (yv3 + yv2 + yv6 + yv7));
        dzv = Real_t(-0.25) * ((zv0 + zv1 + zv5 + zv4) - (zv3 + zv2 + zv6 + zv7));

        domain.delv_eta(i) = ax * dxv + ay * dyv + az * dzv;
    }

}

/******************************************/

static inline
void CalcMonotonicQRegionForElems(Domain &domain, Int_t r,
                                  Real_t vnew[], Real_t ptiny) {
    Real_t monoq_limiter_mult = domain.monoq_limiter_mult();
    Real_t monoq_max_slope = domain.monoq_max_slope();
    Real_t qlc_monoq = domain.qlc_monoq();
    Real_t qqc_monoq = domain.qqc_monoq();

#pragma omp parallel for firstprivate(qlc_monoq, qqc_monoq, monoq_limiter_mult, monoq_max_slope, ptiny)
    for (Index_t ielem = 0; ielem < domain.regElemSize(r); ++ielem) {
        Index_t i = domain.regElemlist(r, ielem);
        Real_t qlin, qquad;
        Real_t phixi, phieta, phizeta;
        Int_t bcMask = domain.elemBC(i);
        Real_t delvm = 0.0, delvp = 0.0;

        /*  phixi     */
        Real_t norm = Real_t(1.) / (domain.delv_xi(i) + ptiny);

        switch (bcMask & XI_M) {
            case XI_M_COMM: /* needs comm data */
            case 0:
                delvm = domain.delv_xi(domain.lxim(i));
                break;
            case XI_M_SYMM:
                delvm = domain.delv_xi(i);
                break;
            case XI_M_FREE:
                delvm = Real_t(0.0);
                break;
            default:
                fprintf(stderr, "Error in switch at %s line %d\n",
                        __FILE__, __LINE__);
                delvm = 0; /* ERROR - but quiets the compiler */
                break;
        }
        switch (bcMask & XI_P) {
            case XI_P_COMM: /* needs comm data */
            case 0:
                delvp = domain.delv_xi(domain.lxip(i));
                break;
            case XI_P_SYMM:
                delvp = domain.delv_xi(i);
                break;
            case XI_P_FREE:
                delvp = Real_t(0.0);
                break;
            default:
                fprintf(stderr, "Error in switch at %s line %d\n",
                        __FILE__, __LINE__);
                delvp = 0; /* ERROR - but quiets the compiler */
                break;
        }

        /*
      laik_log((Laik_LogLevel)2, "r:%d, size:%d, ielem:%d, i:%d, lxim:%d, lxip:%d", r, domain.regElemSize(r), ielem, i, domain.lxim(i), domain.lxip(i));
      */

        delvm = delvm * norm;
        delvp = delvp * norm;

        phixi = Real_t(.5) * (delvm + delvp);

        delvm *= monoq_limiter_mult;
        delvp *= monoq_limiter_mult;

        if (delvm < phixi) phixi = delvm;
        if (delvp < phixi) phixi = delvp;
        if (phixi < Real_t(0.)) phixi = Real_t(0.);
        if (phixi > monoq_max_slope) phixi = monoq_max_slope;


        /*  phieta     */
        norm = Real_t(1.) / (domain.delv_eta(i) + ptiny);

        switch (bcMask & ETA_M) {
            case ETA_M_COMM: /* needs comm data */
            case 0:
                delvm = domain.delv_eta(domain.letam(i));
                break;
            case ETA_M_SYMM:
                delvm = domain.delv_eta(i);
                break;
            case ETA_M_FREE:
                delvm = Real_t(0.0);
                break;
            default:
                fprintf(stderr, "Error in switch at %s line %d\n",
                        __FILE__, __LINE__);
                delvm = 0; /* ERROR - but quiets the compiler */
                break;
        }
        switch (bcMask & ETA_P) {
            case ETA_P_COMM: /* needs comm data */
            case 0:
                delvp = domain.delv_eta(domain.letap(i));
                break;
            case ETA_P_SYMM:
                delvp = domain.delv_eta(i);
                break;
            case ETA_P_FREE:
                delvp = Real_t(0.0);
                break;
            default:
                fprintf(stderr, "Error in switch at %s line %d\n",
                        __FILE__, __LINE__);
                delvp = 0; /* ERROR - but quiets the compiler */
                break;
        }

        delvm = delvm * norm;
        delvp = delvp * norm;

        phieta = Real_t(.5) * (delvm + delvp);

        delvm *= monoq_limiter_mult;
        delvp *= monoq_limiter_mult;

        if (delvm < phieta) phieta = delvm;
        if (delvp < phieta) phieta = delvp;
        if (phieta < Real_t(0.)) phieta = Real_t(0.);
        if (phieta > monoq_max_slope) phieta = monoq_max_slope;

        /*  phizeta     */
        norm = Real_t(1.) / (domain.delv_zeta(i) + ptiny);

        switch (bcMask & ZETA_M) {
            case ZETA_M_COMM: /* needs comm data */
            case 0:
                delvm = domain.delv_zeta(domain.lzetam(i));
                break;
            case ZETA_M_SYMM:
                delvm = domain.delv_zeta(i);
                break;
            case ZETA_M_FREE:
                delvm = Real_t(0.0);
                break;
            default:
                fprintf(stderr, "Error in switch at %s line %d\n",
                        __FILE__, __LINE__);
                delvm = 0; /* ERROR - but quiets the compiler */
                break;
        }
        switch (bcMask & ZETA_P) {
            case ZETA_P_COMM: /* needs comm data */
            case 0:
                delvp = domain.delv_zeta(domain.lzetap(i));
                break;
            case ZETA_P_SYMM:
                delvp = domain.delv_zeta(i);
                break;
            case ZETA_P_FREE:
                delvp = Real_t(0.0);
                break;
            default:
                fprintf(stderr, "Error in switch at %s line %d\n",
                        __FILE__, __LINE__);
                delvp = 0; /* ERROR - but quiets the compiler */
                break;
        }

        delvm = delvm * norm;
        delvp = delvp * norm;

        phizeta = Real_t(.5) * (delvm + delvp);

        delvm *= monoq_limiter_mult;
        delvp *= monoq_limiter_mult;

        if (delvm < phizeta) phizeta = delvm;
        if (delvp < phizeta) phizeta = delvp;
        if (phizeta < Real_t(0.)) phizeta = Real_t(0.);
        if (phizeta > monoq_max_slope) phizeta = monoq_max_slope;

        /* Remove length scale */

        if (domain.vdov(i) > Real_t(0.)) {
            qlin = Real_t(0.);
            qquad = Real_t(0.);
        } else {
            Real_t delvxxi = domain.delv_xi(i) * domain.delx_xi(i);
            Real_t delvxeta = domain.delv_eta(i) * domain.delx_eta(i);
            Real_t delvxzeta = domain.delv_zeta(i) * domain.delx_zeta(i);

            if (delvxxi > Real_t(0.)) delvxxi = Real_t(0.);
            if (delvxeta > Real_t(0.)) delvxeta = Real_t(0.);
            if (delvxzeta > Real_t(0.)) delvxzeta = Real_t(0.);

            Real_t rho = domain.elemMass(i) / (domain.volo(i) * vnew[i]);

            qlin = -qlc_monoq * rho *
                   (delvxxi * (Real_t(1.) - phixi) +
                    delvxeta * (Real_t(1.) - phieta) +
                    delvxzeta * (Real_t(1.) - phizeta));

            qquad = qqc_monoq * rho *
                    (delvxxi * delvxxi * (Real_t(1.) - phixi * phixi) +
                     delvxeta * delvxeta * (Real_t(1.) - phieta * phieta) +
                     delvxzeta * delvxzeta * (Real_t(1.) - phizeta * phizeta));
        }

        domain.qq(i) = qquad;
        domain.ql(i) = qlin;
    }
}

/******************************************/

static inline
void CalcMonotonicQForElems(Domain &domain, Real_t vnew[]) {
    //
    // initialize parameters
    // 
    const Real_t ptiny = Real_t(1.e-36);

    //
    // calculate the monotonic q for all regions
    //
    for (Index_t r = 0; r < domain.numReg(); ++r) {

        if (domain.regElemSize(r) > 0) {
            CalcMonotonicQRegionForElems(domain, r, vnew, ptiny);
        }
    }
}

/******************************************/

static inline
void CalcQForElems(Domain &domain, Real_t vnew[]) {
    //
    // MONOTONIC Q option
    //

    Index_t numElem = domain.numElem();

    if (numElem != 0) {

        //Int_t allElem = numElem +  /* local elem */
        //     2*domain.sizeX()*domain.sizeY() + /* plane ghosts */
        //     2*domain.sizeX()*domain.sizeZ() + /* row ghosts */
        //     2*domain.sizeY()*domain.sizeZ() ; /* col ghosts */
        //domain.AllocateGradients(numElem, numRanks, allElem);

#if USE_MPI
        //CommRecv(domain, MSG_MONOQ, 3,
        //         domain.sizeX(), domain.sizeY(), domain.sizeZ(),
        //         true, true) ;
        domain.get_delv_xi().switch_to_p1();
        domain.get_delv_eta().switch_to_p1();
        domain.get_delv_zeta().switch_to_p1();
#endif

        /* Calculate velocity gradients */
        CalcMonotonicQGradientsForElems(domain, vnew);

#if USE_MPI

        domain.get_delv_xi().switch_to_p2();
        domain.get_delv_eta().switch_to_p2();
        domain.get_delv_zeta().switch_to_p2();

        //Domain_member fieldData[3] ;

        /* Transfer veloctiy gradients in the first order elements */
        /* problem->commElements->Transfer(CommElements::monoQ) ; */

        //fieldData[0] = &Domain::delv_xi ;
        //fieldData[1] = &Domain::delv_eta ;
        //fieldData[2] = &Domain::delv_zeta ;

        //CommSend(domain, MSG_MONOQ, 3, fieldData,
        //         domain.sizeX(), domain.sizeY(), domain.sizeZ(),
        //         true, true) ;

        //CommMonoQ(domain) ;
#endif


        CalcMonotonicQForElems(domain, vnew);

        // Free up memory
        //domain.DeallocateGradients();

        /*
      for (Index_t i=0; i<20; ++i) {
          laik_log((Laik_LogLevel)2, "i: %d, value: %f", i, domain.delv_xi(i));
      }
      */


        /* Don't allow excessive artificial viscosity */
        Index_t idx = -1;
        for (Index_t i = 0; i < numElem; ++i) {
            if (domain.q(i) > domain.qstop()) {
                idx = i;
                break;
            }
        }

        if (idx >= 0) {
#if USE_MPI
            laik_log(LAIK_LL_Warning, "Debug 4\n");
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Abort(MPI_COMM_WORLD, QStopError);
#else
            exit(QStopError);
#endif
        }
    }
}

/******************************************/

static inline
void CalcPressureForElems(Real_t *p_new, Real_t *bvc,
                          Real_t *pbvc, Real_t *e_old,
                          Real_t *compression, Real_t *vnewc,
                          Real_t pmin,
                          Real_t p_cut, Real_t eosvmax,
                          Index_t length, std::vector<Index_t> regElemList) {
#pragma omp parallel for firstprivate(length)
    for (Index_t i = 0; i < length; ++i) {
        Real_t c1s = Real_t(2.0) / Real_t(3.0);
        bvc[i] = c1s * (compression[i] + Real_t(1.));
        pbvc[i] = c1s;
    }

#pragma omp parallel for firstprivate(length, pmin, p_cut, eosvmax)
    for (Index_t i = 0; i < length; ++i) {
        Index_t elem = regElemList[i];

        p_new[i] = bvc[i] * e_old[i];

        if (FABS(p_new[i]) < p_cut)
            p_new[i] = Real_t(0.0);

        if (vnewc[elem] >= eosvmax) /* impossible condition here? */
            p_new[i] = Real_t(0.0);

        if (p_new[i] < pmin)
            p_new[i] = pmin;
    }
}

/******************************************/

static inline
void CalcEnergyForElems(Real_t *p_new, Real_t *e_new, Real_t *q_new,
                        Real_t *bvc, Real_t *pbvc,
                        Real_t *p_old, Real_t *e_old, Real_t *q_old,
                        Real_t *compression, Real_t *compHalfStep,
                        Real_t *vnewc, Real_t *work, Real_t *delvc, Real_t pmin,
                        Real_t p_cut, Real_t e_cut, Real_t q_cut, Real_t emin,
                        Real_t *qq_old, Real_t *ql_old,
                        Real_t rho0,
                        Real_t eosvmax,
                        Index_t length, std::vector<Index_t> regElemList) {
    Real_t *pHalfStep = Allocate<Real_t>(length);

#pragma omp parallel for firstprivate(length, emin)
    for (Index_t i = 0; i < length; ++i) {
        e_new[i] = e_old[i] - Real_t(0.5) * delvc[i] * (p_old[i] + q_old[i])
                   + Real_t(0.5) * work[i];

        if (e_new[i] < emin) {
            e_new[i] = emin;
        }
    }

    CalcPressureForElems(pHalfStep, bvc, pbvc, e_new, compHalfStep, vnewc,
                         pmin, p_cut, eosvmax, length, regElemList);

#pragma omp parallel for firstprivate(length, rho0)
    for (Index_t i = 0; i < length; ++i) {
        Real_t vhalf = Real_t(1.) / (Real_t(1.) + compHalfStep[i]);

        if (delvc[i] > Real_t(0.)) {
            q_new[i] /* = qq_old[i] = ql_old[i] */ = Real_t(0.);
        } else {
            Real_t ssc = (pbvc[i] * e_new[i]
                          + vhalf * vhalf * bvc[i] * pHalfStep[i]) / rho0;

            if (ssc <= Real_t(.1111111e-36)) {
                ssc = Real_t(.3333333e-18);
            } else {
                ssc = SQRT(ssc);
            }

            q_new[i] = (ssc * ql_old[i] + qq_old[i]);
        }

        e_new[i] = e_new[i] + Real_t(0.5) * delvc[i]
                              * (Real_t(3.0) * (p_old[i] + q_old[i])
                                 - Real_t(4.0) * (pHalfStep[i] + q_new[i]));
    }

#pragma omp parallel for firstprivate(length, emin, e_cut)
    for (Index_t i = 0; i < length; ++i) {

        e_new[i] += Real_t(0.5) * work[i];

        if (FABS(e_new[i]) < e_cut) {
            e_new[i] = Real_t(0.);
        }
        if (e_new[i] < emin) {
            e_new[i] = emin;
        }
    }

    CalcPressureForElems(p_new, bvc, pbvc, e_new, compression, vnewc,
                         pmin, p_cut, eosvmax, length, regElemList);

#pragma omp parallel for firstprivate(length, rho0, emin, e_cut)
    for (Index_t i = 0; i < length; ++i) {
        const Real_t sixth = Real_t(1.0) / Real_t(6.0);
        Index_t elem = regElemList[i];
        Real_t q_tilde;

        if (delvc[i] > Real_t(0.)) {
            q_tilde = Real_t(0.);
        } else {
            Real_t ssc = (pbvc[i] * e_new[i]
                          + vnewc[elem] * vnewc[elem] * bvc[i] * p_new[i]) / rho0;

            if (ssc <= Real_t(.1111111e-36)) {
                ssc = Real_t(.3333333e-18);
            } else {
                ssc = SQRT(ssc);
            }

            q_tilde = (ssc * ql_old[i] + qq_old[i]);
        }

        e_new[i] = e_new[i] - (Real_t(7.0) * (p_old[i] + q_old[i])
                               - Real_t(8.0) * (pHalfStep[i] + q_new[i])
                               + (p_new[i] + q_tilde)) * delvc[i] * sixth;

        if (FABS(e_new[i]) < e_cut) {
            e_new[i] = Real_t(0.);
        }
        if (e_new[i] < emin) {
            e_new[i] = emin;
        }
    }

    CalcPressureForElems(p_new, bvc, pbvc, e_new, compression, vnewc,
                         pmin, p_cut, eosvmax, length, regElemList);

#pragma omp parallel for firstprivate(length, rho0, q_cut)
    for (Index_t i = 0; i < length; ++i) {
        Index_t elem = regElemList[i];

        if (delvc[i] <= Real_t(0.)) {
            Real_t ssc = (pbvc[i] * e_new[i]
                          + vnewc[elem] * vnewc[elem] * bvc[i] * p_new[i]) / rho0;

            if (ssc <= Real_t(.1111111e-36)) {
                ssc = Real_t(.3333333e-18);
            } else {
                ssc = SQRT(ssc);
            }

            q_new[i] = (ssc * ql_old[i] + qq_old[i]);

            if (FABS(q_new[i]) < q_cut) q_new[i] = Real_t(0.);
        }
    }

    Release(&pHalfStep);

    return;
}

/******************************************/

static inline
void CalcSoundSpeedForElems(Domain &domain,
                            Real_t *vnewc, Real_t rho0, Real_t *enewc,
                            Real_t *pnewc, Real_t *pbvc,
                            Real_t *bvc, Real_t ss4o3,
                            Index_t len, std::vector<Index_t> regElemList) {
#pragma omp parallel for firstprivate(rho0, ss4o3)
    for (Index_t i = 0; i < len; ++i) {
        Index_t elem = regElemList[i];
        Real_t ssTmp = (pbvc[i] * enewc[i] + vnewc[elem] * vnewc[elem] *
                                             bvc[i] * pnewc[i]) / rho0;
        if (ssTmp <= Real_t(.1111111e-36)) {
            ssTmp = Real_t(.3333333e-18);
        } else {
            ssTmp = SQRT(ssTmp);
        }
        domain.ss(elem) = ssTmp;
    }
}

/******************************************/

static inline
void EvalEOSForElems(Domain &domain, Real_t *vnewc,
                     Int_t numElemReg, std::vector<Index_t> regElemList, Int_t rep) {
    Real_t e_cut = domain.e_cut();
    Real_t p_cut = domain.p_cut();
    Real_t ss4o3 = domain.ss4o3();
    Real_t q_cut = domain.q_cut();

    Real_t eosvmax = domain.eosvmax();
    Real_t eosvmin = domain.eosvmin();
    Real_t pmin = domain.pmin();
    Real_t emin = domain.emin();
    Real_t rho0 = domain.refdens();

    // These temporaries will be of different size for 
    // each call (due to different sized region element
    // lists)
    Real_t *e_old = Allocate<Real_t>(numElemReg);
    Real_t *delvc = Allocate<Real_t>(numElemReg);
    Real_t *p_old = Allocate<Real_t>(numElemReg);
    Real_t *q_old = Allocate<Real_t>(numElemReg);
    Real_t *compression = Allocate<Real_t>(numElemReg);
    Real_t *compHalfStep = Allocate<Real_t>(numElemReg);
    Real_t *qq_old = Allocate<Real_t>(numElemReg);
    Real_t *ql_old = Allocate<Real_t>(numElemReg);
    Real_t *work = Allocate<Real_t>(numElemReg);
    Real_t *p_new = Allocate<Real_t>(numElemReg);
    Real_t *e_new = Allocate<Real_t>(numElemReg);
    Real_t *q_new = Allocate<Real_t>(numElemReg);
    Real_t *bvc = Allocate<Real_t>(numElemReg);
    Real_t *pbvc = Allocate<Real_t>(numElemReg);

    //loop to add load imbalance based on region number 
    for (Int_t j = 0; j < rep; j++) {
        /* compress data, minimal set */
#pragma omp parallel
        {
#pragma omp for nowait firstprivate(numElemReg)
            for (Index_t i = 0; i < numElemReg; ++i) {
                Index_t elem = regElemList[i];
                e_old[i] = domain.e(elem);
                delvc[i] = domain.delv(elem);
                p_old[i] = domain.p(elem);
                q_old[i] = domain.q(elem);
                qq_old[i] = domain.qq(elem);
                ql_old[i] = domain.ql(elem);
            }

#pragma omp for firstprivate(numElemReg)
            for (Index_t i = 0; i < numElemReg; ++i) {
                Index_t elem = regElemList[i];
                Real_t vchalf;
                compression[i] = Real_t(1.) / vnewc[elem] - Real_t(1.);
                vchalf = vnewc[elem] - delvc[i] * Real_t(.5);
                compHalfStep[i] = Real_t(1.) / vchalf - Real_t(1.);
            }

            /* Check for v > eosvmax or v < eosvmin */
            if (eosvmin != Real_t(0.)) {
#pragma omp for nowait firstprivate(numElemReg, eosvmin)
                for (Index_t i = 0; i < numElemReg; ++i) {
                    Index_t elem = regElemList[i];
                    if (vnewc[elem] <= eosvmin) { /* impossible due to calling func? */
                        compHalfStep[i] = compression[i];
                    }
                }
            }
            if (eosvmax != Real_t(0.)) {
#pragma omp for nowait firstprivate(numElemReg, eosvmax)
                for (Index_t i = 0; i < numElemReg; ++i) {
                    Index_t elem = regElemList[i];
                    if (vnewc[elem] >= eosvmax) { /* impossible due to calling func? */
                        p_old[i] = Real_t(0.);
                        compression[i] = Real_t(0.);
                        compHalfStep[i] = Real_t(0.);
                    }
                }
            }

#pragma omp for nowait firstprivate(numElemReg)
            for (Index_t i = 0; i < numElemReg; ++i) {
                work[i] = Real_t(0.);
            }
        }
        CalcEnergyForElems(p_new, e_new, q_new, bvc, pbvc,
                           p_old, e_old, q_old, compression, compHalfStep,
                           vnewc, work, delvc, pmin,
                           p_cut, e_cut, q_cut, emin,
                           qq_old, ql_old, rho0, eosvmax,
                           numElemReg, regElemList);
    }

#pragma omp parallel for firstprivate(numElemReg)
    for (Index_t i = 0; i < numElemReg; ++i) {
        Index_t elem = regElemList[i];
        domain.p(elem) = p_new[i];
        domain.e(elem) = e_new[i];
        domain.q(elem) = q_new[i];
    }

    CalcSoundSpeedForElems(domain,
                           vnewc, rho0, e_new, p_new,
                           pbvc, bvc, ss4o3,
                           numElemReg, regElemList);

    Release(&pbvc);
    Release(&bvc);
    Release(&q_new);
    Release(&e_new);
    Release(&p_new);
    Release(&work);
    Release(&ql_old);
    Release(&qq_old);
    Release(&compHalfStep);
    Release(&compression);
    Release(&q_old);
    Release(&p_old);
    Release(&delvc);
    Release(&e_old);
}

/******************************************/

static inline
void ApplyMaterialPropertiesForElems(Domain &domain, Real_t vnew[]) {
    Index_t numElem = domain.numElem();

    if (numElem != 0) {
        /* Expose all of the variables needed for material evaluation */
        Real_t eosvmin = domain.eosvmin();
        Real_t eosvmax = domain.eosvmax();

#pragma omp parallel
        {
            // Bound the updated relative volumes with eosvmin/max
            if (eosvmin != Real_t(0.)) {
#pragma omp for firstprivate(numElem)
                for (Index_t i = 0; i < numElem; ++i) {
                    if (vnew[i] < eosvmin)
                        vnew[i] = eosvmin;
                }
            }

            if (eosvmax != Real_t(0.)) {
#pragma omp for nowait firstprivate(numElem)
                for (Index_t i = 0; i < numElem; ++i) {
                    if (vnew[i] > eosvmax)
                        vnew[i] = eosvmax;
                }
            }

            // This check may not make perfect sense in LULESH, but
            // it's representative of something in the full code -
            // just leave it in, please
#pragma omp for nowait firstprivate(numElem)
            for (Index_t i = 0; i < numElem; ++i) {
                Real_t vc = domain.v(i);
                if (eosvmin != Real_t(0.)) {
                    if (vc < eosvmin)
                        vc = eosvmin;
                }
                if (eosvmax != Real_t(0.)) {
                    if (vc > eosvmax)
                        vc = eosvmax;
                }
                if (vc <= 0.) {
#if USE_MPI
                    laik_log(LAIK_LL_Warning, "Debug 5\n");
                    MPI_Barrier(MPI_COMM_WORLD);
                    MPI_Abort(MPI_COMM_WORLD, VolumeError);
#else
                    exit(VolumeError);
#endif
                }
            }
        }

        for (Int_t r = 0; r < domain.numReg(); r++) {
            Index_t numElemReg = domain.regElemSize(r);
            std::vector<Index_t> regElemList = domain.regElemlist(r);
            Int_t rep;
            //Determine load imbalance for this region
            //round down the number with lowest cost
            if (r < domain.numReg() / 2)
                rep = 1;
                //you don't get an expensive region unless you at least have 5 regions
            else if (r < (domain.numReg() - (domain.numReg() + 15) / 20))
                rep = 1 + domain.cost();
                //very expensive regions
            else
                rep = 10 * (1 + domain.cost());
            EvalEOSForElems(domain, vnew, numElemReg, regElemList, rep);
        }

    }
}

/******************************************/

static inline
void UpdateVolumesForElems(Domain &domain, Real_t *vnew,
                           Real_t v_cut, Index_t length) {
    if (length != 0) {
#pragma omp parallel for firstprivate(length, v_cut)
        for (Index_t i = 0; i < length; ++i) {
            Real_t tmpV = vnew[i];

            if (FABS(tmpV - Real_t(1.0)) < v_cut)
                tmpV = Real_t(1.0);

            domain.v(i) = tmpV;
        }
    }

    return;
}

/******************************************/

static inline
void LagrangeElements(Domain &domain, Index_t numElem) {
    Real_t *vnew = Allocate<Real_t>(numElem);  /* new relative vol -- temp */

    CalcLagrangeElements(domain, vnew);

    /* Calculate Q.  (Monotonic q option requires communication) */
    CalcQForElems(domain, vnew);

    ApplyMaterialPropertiesForElems(domain, vnew);

    UpdateVolumesForElems(domain, vnew,
                          domain.v_cut(), numElem);

    Release(&vnew);
}

/******************************************/

static inline
void CalcCourantConstraintForElems(Domain &domain, Index_t length,
                                   std::vector<Index_t> regElemlist,
                                   Real_t qqc, Real_t &dtcourant) {
#if _OPENMP
                                                                                                                            Index_t threads = omp_get_max_threads();
   static Index_t *courant_elem_per_thread;
   static Real_t *dtcourant_per_thread;
   static bool first = true;
   if (first) {
     courant_elem_per_thread = new Index_t[threads];
     dtcourant_per_thread = new Real_t[threads];
     first = false;
   }
#else
    Index_t threads = 1;
    Index_t courant_elem_per_thread[1];
    Real_t dtcourant_per_thread[1];
#endif


#pragma omp parallel firstprivate(length, qqc)
    {
        Real_t qqc2 = Real_t(64.0) * qqc * qqc;
        Real_t dtcourant_tmp = dtcourant;
        Index_t courant_elem = -1;

#if _OPENMP
        Index_t thread_num = omp_get_thread_num();
#else
        Index_t thread_num = 0;
#endif

#pragma omp for
        for (Index_t i = 0; i < length; ++i) {
            Index_t indx = regElemlist[i];
            Real_t dtf = domain.ss(indx) * domain.ss(indx);

            if (domain.vdov(indx) < Real_t(0.)) {
                dtf = dtf
                      + qqc2 * domain.arealg(indx) * domain.arealg(indx)
                        * domain.vdov(indx) * domain.vdov(indx);
            }

            dtf = SQRT(dtf);
            dtf = domain.arealg(indx) / dtf;

            if (domain.vdov(indx) != Real_t(0.)) {
                if (dtf < dtcourant_tmp) {
                    dtcourant_tmp = dtf;
                    courant_elem = indx;
                }
            }
        }

        dtcourant_per_thread[thread_num] = dtcourant_tmp;
        courant_elem_per_thread[thread_num] = courant_elem;
    }

    for (Index_t i = 1; i < threads; ++i) {
        if (dtcourant_per_thread[i] < dtcourant_per_thread[0]) {
            dtcourant_per_thread[0] = dtcourant_per_thread[i];
            courant_elem_per_thread[0] = courant_elem_per_thread[i];
        }
    }

    if (courant_elem_per_thread[0] != -1) {
        dtcourant = dtcourant_per_thread[0];
    }

    return;

}

/******************************************/

static inline
void CalcHydroConstraintForElems(Domain &domain, Index_t length,
                                 std::vector<Index_t> regElemlist, Real_t dvovmax, Real_t &dthydro) {
#if _OPENMP
                                                                                                                            Index_t threads = omp_get_max_threads();
   static Index_t *hydro_elem_per_thread;
   static Real_t *dthydro_per_thread;
   static bool first = true;
   if (first) {
     hydro_elem_per_thread = new Index_t[threads];
     dthydro_per_thread = new Real_t[threads];
     first = false;
   }
#else
    Index_t threads = 1;
    Index_t hydro_elem_per_thread[1];
    Real_t dthydro_per_thread[1];
#endif

#pragma omp parallel firstprivate(length, dvovmax)
    {
        Real_t dthydro_tmp = dthydro;
        Index_t hydro_elem = -1;

#if _OPENMP
        Index_t thread_num = omp_get_thread_num();
#else
        Index_t thread_num = 0;
#endif

#pragma omp for
        for (Index_t i = 0; i < length; ++i) {
            Index_t indx = regElemlist[i];

            if (domain.vdov(indx) != Real_t(0.)) {
                Real_t dtdvov = dvovmax / (FABS(domain.vdov(indx)) + Real_t(1.e-20));

                if (dthydro_tmp > dtdvov) {
                    dthydro_tmp = dtdvov;
                    hydro_elem = indx;
                }
            }
        }

        dthydro_per_thread[thread_num] = dthydro_tmp;
        hydro_elem_per_thread[thread_num] = hydro_elem;
    }

    for (Index_t i = 1; i < threads; ++i) {
        if (dthydro_per_thread[i] < dthydro_per_thread[0]) {
            dthydro_per_thread[0] = dthydro_per_thread[i];
            hydro_elem_per_thread[0] = hydro_elem_per_thread[i];
        }
    }

    if (hydro_elem_per_thread[0] != -1) {
        dthydro = dthydro_per_thread[0];
    }

    return;
}

/******************************************/

static inline
void CalcTimeConstraintsForElems(Domain &domain) {

    // Initialize conditions to a very large value
    domain.dtcourant() = 1.0e+20;
    domain.dthydro() = 1.0e+20;

    for (Index_t r = 0; r < domain.numReg(); ++r) {
        /* evaluate time constraint */
        CalcCourantConstraintForElems(domain, domain.regElemSize(r),
                                      domain.regElemlist(r),
                                      domain.qqc(),
                                      domain.dtcourant());

        /* check hydro constraint */
        CalcHydroConstraintForElems(domain, domain.regElemSize(r),
                                    domain.regElemlist(r),
                                    domain.dvovmax(),
                                    domain.dthydro());
    }
}

/******************************************/

static inline
void LagrangeLeapFrog(Domain &domain) {
#ifdef SEDOV_SYNC_POS_VEL_LATE
    Domain_member fieldData[6] ;
#endif

    /* calculate nodal forces, accelerations, velocities, positions, with
    * applied boundary conditions and slide surface considerations */
    LagrangeNodal(domain);

#ifdef SEDOV_SYNC_POS_VEL_LATE
#endif
    /* calculate element quantities (i.e. velocity gradient & q), and update
    * material states */
    LagrangeElements(domain, domain.numElem());

#if USE_MPI
#ifdef SEDOV_SYNC_POS_VEL_LATE
                                                                                                                            CommRecv(domain, MSG_SYNC_POS_VEL, 6,
            domain.sizeX() + 1, domain.sizeY() + 1, domain.sizeZ() + 1,
            false, false) ;

   fieldData[0] = &Domain::x ;
   fieldData[1] = &Domain::y ;
   fieldData[2] = &Domain::z ;
   fieldData[3] = &Domain::xd ;
   fieldData[4] = &Domain::yd ;
   fieldData[5] = &Domain::zd ;
   
   CommSend(domain, MSG_SYNC_POS_VEL, 6, fieldData,
            domain.sizeX() + 1, domain.sizeY() + 1, domain.sizeZ() + 1,
            false, false) ;
#endif
#endif

    CalcTimeConstraintsForElems(domain);

#if USE_MPI
#ifdef SEDOV_SYNC_POS_VEL_LATE
    CommSyncPosVel(domain) ;
#endif
#endif
}

/******************************************/

int main(int argc, char *argv[]) {
    Domain *locDom;
    Int_t numRanks;
    Int_t myRank;
    struct cmdLineOpts opts;


#if USE_MPI
    //Domain_member fieldData ;

    Laik_Instance *inst = laik_init_mpi(&argc, &argv);
    Laik_Group *world = laik_world(inst);

    TRACE_INIT(laik_myid(world));
    TRACE_EVENT_START("INIT", "");

    numRanks = laik_size(world);
    myRank = laik_myid(world);

    //MPI_Comm_size(MPI_COMM_WORLD, &numRanks) ;
    //MPI_Comm_rank(MPI_COMM_WORLD, &myRank) ;
#else
                                                                                                                            numRanks = 1;
   myRank = 0;
#endif

#ifdef FAULT_TOLERANCE
    std::vector<Laik_Checkpoint*> checkpoints;
#endif

    /* Set defaults that can be overridden by command line opts */
    opts.its = 9999999;
    opts.nx = 30;
    opts.numReg = 11;
    opts.numFiles = (int) (numRanks + 10) / 9;
    opts.showProg = 0;
    opts.quiet = 0;
    opts.viz = 0;
    opts.balance = 1;
    opts.cost = 1;
    opts.repart = 0;
    opts.cycle = 5;
    opts.noTimer = false;

    FaultToleranceOptions ftOptions = FaultToleranceOptionsDefault;

    ParseCommandLineOptions(argc, argv, myRank, &opts, &ftOptions);

    if ((myRank == 0) && (opts.quiet == 0)) {
        printf("Running problem size %d^3 per domain until completion\n", opts.nx);
        printf("Num processors: %d\n", numRanks);
#if _OPENMP
        printf("Num threads: %d\n", omp_get_max_threads());
#endif
        printf("Total number of elements: %lld\n\n", (long long int) (numRanks * opts.nx * opts.nx * opts.nx));
        printf("To run other sizes, use -s <integer>.\n");
        printf("To run a fixed number of iterations, use -i <integer>.\n");
        printf("To run a more or less balanced region set, use -b <integer>.\n");
        printf("To change the relative costs of regions, use -c <integer>.\n");
        printf("To print out progress, use -p\n");
        printf("To write an output file for VisIt, use -v\n");
        printf("See help (-h) for more options\n\n");

#ifdef PERFORMANCE
        printf(" \n=== TESTS ON PERFORMANCE MODE === \n");
#endif
#ifdef REPARTITIONING
        printf(" \n=== TESTS ON REPARTITIONING MODE === \n");
#endif
    }

    // Set up the mesh and decompose. Assumes regular cubes for now
    Int_t col, row, plane, side;
    InitMeshDecomp(numRanks, myRank, &col, &row, &plane, &side);

    // laik related code
    int numElem = opts.nx * opts.nx * opts.nx;                           // local number of elements
    int NumNodes = (opts.nx * side + 1) * (opts.nx * side + 1) * (opts.nx * side + 1);     // global number of nodes

    // index spaces used for partitioning of Element and Node data structures
    Laik_Space *indexSpaceElements = nullptr;                            // Elements
    Laik_Space *indexSpaceNodes = nullptr;                               // Nodes
    Laik_Space *indexSapceDt = nullptr;                                  // for reduction of dt

    // partitioning of Element and Node data structures
    Laik_Partitioning *exclusivePartitioning = nullptr;
    Laik_Partitioning *haloPartitioning = nullptr;
    Laik_Partitioning *overlapingPartitioning = nullptr;
    Laik_Partitioning *allPartitioning = nullptr;
    // partitioning objects used after re-partitioning
    Laik_Partitioning *exclusivePartitioning2 = nullptr;
    Laik_Partitioning *haloPartitioning2 = nullptr;
    Laik_Partitioning *overlapingPartitioning2 = nullptr;
    Laik_Partitioning *allPartitioning2 = nullptr;

    // transition objects among partitionings.
    Laik_Transition *transitionToExclusive = nullptr;
    Laik_Transition *transitionToHalo = nullptr;
    Laik_Transition *transitionToOverlappingInit = nullptr;
    Laik_Transition *transitionToOverlappingReduce = nullptr;
    // transition objects used after repartitioning
    Laik_Transition *transitionToExclusive2 = nullptr;
    Laik_Transition *transitionToHalo2 = nullptr;
    Laik_Transition *transitionToOverlappingInit2 = nullptr;
    Laik_Transition *transitionToOverlappingReduce2 = nullptr;

    // Laik data container for dt
    Laik_Data *laikDt = nullptr;

    Laik_Data *laikReductionAuxiliaryContainer = nullptr;

    // create index spaces for elements, nodes and dt
    indexSpaceElements = laik_new_space_1d(inst, numRanks * numElem);
    indexSpaceNodes = laik_new_space_1d(inst, NumNodes);
    indexSapceDt = laik_new_space_1d(inst, 1);

    // create partitionings and transition objects
    create_partitionings_and_transitions(world,
                                         indexSpaceElements, indexSpaceNodes, indexSapceDt,
                                         exclusivePartitioning, haloPartitioning, overlapingPartitioning,
                                         allPartitioning,
                                         transitionToExclusive, transitionToHalo, transitionToOverlappingInit,
                                         transitionToOverlappingReduce);

    // create laik_data for dt
    // Note the rest of containers are
    // created inside inside Domain by
    // initializing laik_vectors
    laikDt = laik_new_data(indexSapceDt, laik_Double);

    laikReductionAuxiliaryContainer = laik_new_data_1d(inst, laik_Double, 1);
    Laik_Partitioning *all, *master;
    all = laik_new_partitioning(laik_All, world, laik_data_get_space(laikReductionAuxiliaryContainer), nullptr);
    master = laik_new_partitioning(laik_Master, world, laik_data_get_space(laikReductionAuxiliaryContainer), nullptr);
    laik_switchto_partitioning(laikReductionAuxiliaryContainer, all, LAIK_DF_None, LAIK_RO_None);

    // for dt we need to get the base pointer to and
    // send it to TimeIncrement function to
    // perform the local updates and reduction
    laik_switchto_partitioning(laikDt, allPartitioning, LAIK_DF_Init, LAIK_RO_Min);
    uint64_t dt_count;
    double *dt_base;
    assert(laik_my_slicecount(laik_data_get_partitioning(laikDt)) == 1);
    laik_get_map_1d(laikDt, 0, (void **) &dt_base, &dt_count);

    // Build the main data structure and initialize it
    // pass the laik_inst and laik_world and partitionings
    // and transitions to Domain since Domain contains
    // all the data structures and inititlizes them.
    locDom = new Domain(numRanks, col, row, plane, opts.nx,
                        side, opts.numReg, opts.balance, opts.cost,
                        inst, world,
                        indexSpaceElements, indexSpaceNodes,
                        exclusivePartitioning, haloPartitioning, overlapingPartitioning,
                        transitionToExclusive, transitionToHalo, transitionToOverlappingInit,
                        transitionToOverlappingReduce);
#if USE_MPI

    /*
   fieldData = &Domain::nodalMass ;


   // Initial domain boundary communication 
   CommRecv(*locDom, MSG_COMM_SBN, 1,
            locDom->sizeX() + 1, locDom->sizeY() + 1, locDom->sizeZ() + 1,
            true, false) ;
   CommSend(*locDom, MSG_COMM_SBN, 1, &fieldData,
            locDom->sizeX() + 1, locDom->sizeY() + 1, locDom->sizeZ() +  1,
            true, false) ;
   CommSBN(*locDom, 1, &fieldData) ;


   // End initialization
   MPI_Barrier(MPI_COMM_WORLD);
    */

#endif

    // BEGIN timestep to solution */
#if USE_MPI
    double start = getMPI_Wtime(&opts);
    double start2 = getMPI_Wtime(&opts);
    int repart_iter = 0;
#else
                                                                                                                            timeval start;
   gettimeofday(&start, NULL) ;
#endif

    TRACE_EVENT_END("INIT", "");

    //TODO do not terminate early in case of fault
    while ( /*(locDom->time() < locDom->stoptime()) && */ (locDom->cycle() < opts.its)) {

// TODO remove the following def
#define USE_MPI 1
#if USE_MPI

        laik_set_iteration(inst, locDom->cycle());

        // Execute any pre planned failures
        exitIfFailureIteration(locDom->cycle(), &ftOptions, inst);

        // check if repartitioning has to be done in this iteration
        // if so, do repartitioning before the actual iteration
        // after repartitioning continue from the current iteration
        if ((opts.repart > 0 && locDom->cycle() == opts.cycle) || (isFaultToleranceActive(&ftOptions) && locDom->cycle() % ftOptions.failureCheckFrequency == 0)) {
            std::vector<int> nodeStatuses;
            int failedCount = -1;
            // Check if a node has failed, then do restore. Else, do a checkpoint.
            if(isFaultToleranceActive(&ftOptions)){
                if(myRank == 0) {
                    laik_log(LAIK_LL_Info, "Checking for failed nodes.");
                }
                TRACE_EVENT_S("FAILURE-DETECT", "");
                nodeStatuses.reserve(laik_size(world));
                failedCount = laik_failure_check_nodes(inst, world, &nodeStatuses[0]);
                laik_log(LAIK_LL_Info, "Rank %i detected %i node failures.", myRank, failedCount);
            }
            if (!isFaultToleranceActive(&ftOptions) || failedCount > 0) {
                TRACE_EVENT_START("RESTORE", "");
                double intermediate_timer = getMPI_Wtime(&opts) - start2;
                double itG = -1;

//                MPI_Reduce(&intermediate_timer, &itG, 1, MPI_DOUBLE,
//                           MPI_MAX, 0, MPI_COMM_WORLD);
                itG = laik_reduce(intermediate_timer, myRank, laikReductionAuxiliaryContainer, all, master);

                if ((myRank == 0) && (opts.quiet == 0)) {
                    printf("Starting Repartitioning, current runtime = %f s\n", itG);
                    printf("Previous Time per Iteration: %f s \n", itG / locDom->cycle());
                }

                laik_log((Laik_LogLevel) 2, "Before repart\n");
                laik_log((Laik_LogLevel) 2, "My ID in main world: %d\n", laik_myid(world));

                double repart_start = getMPI_Wtime(&opts);

                // calculate number of removing tasks and a list of them
                double newside;
                int diffsize;
                int *removeList = nullptr;
                Laik_Group *shrinked_group;

                //Todo this was switched logic. Based on previous working version, this was changed back.
                if (isFaultToleranceActive(&ftOptions)) {
                    laik_log(LAIK_LL_Info, "Fault tolerance recovery repartitioning pre-step initiated");
//                    calculate_removing_list_ft(world, opts, side, newside, diffsize, removeList, &nodeStatuses[0]);

                    // Only remove the actually failed before restoring
                    laik_failure_eliminate_nodes(inst, failedCount, &nodeStatuses[0]);
                    shrinked_group = laik_world_fault_tolerant(inst);


                    for(auto& checkpoint : checkpoints) {
                        laik_checkpoint_remove_failed_slices(checkpoint, world, &nodeStatuses[0]);
                    }

                    laik_log(LAIK_LL_Info, "Restoring checkpoints.");
                    int restored = locDom->restore(checkpoints, shrinked_group);
                    laik_log(LAIK_LL_Info, "Restored %i checkpoints.", restored);

                    // Make sure the dt container is also switched to the temporary world before doing the second switch.
                    laik_log(LAIK_LL_Debug, "Switching the dt container to temporary group.");
                    laik_switchto_new_partitioning(laikDt, shrinked_group, laik_All, LAIK_DF_None, LAIK_RO_None);
                    laik_log(LAIK_LL_Debug, "Switched the dt container to temporary group.");

                    laik_log(LAIK_LL_Debug, "Switching the auxiliary container to temporary group.");
                    laik_switchto_new_partitioning(laikReductionAuxiliaryContainer, shrinked_group, laik_All, LAIK_DF_None, LAIK_RO_None);
                    laik_log(LAIK_LL_Debug, "Switched the auxiliary container to temporary group.");

                    // Fake new temporary world
                    world = shrinked_group;
                    laik_log(LAIK_LL_Info, "Installed fake world.");
                }

                calculate_removing_list(world, opts, side, newside, diffsize, removeList);

                // create the shrinked group by removing tasks in removeList
                shrinked_group = laik_new_shrinked_group(world, diffsize, removeList);

                if(laik_size(world) - diffsize == 1) {
                    laik_log(LAIK_LL_Info, "Only 1 process left in group, disabling redundancy.");
                    ftOptions.redundancyCount = 0;
                }

                // update number of local number of elements per edge in each task
                opts.nx = opts.nx * side / (int) newside;


                // create partitionings and transition objects
                // on the target (shrinked) process group
                // we switch to these partitionings and
                // the old objects are not used anymore.
                laik_log(LAIK_LL_Info, "Creating partitionings and transitions.");
                create_partitionings_and_transitions(shrinked_group,
                                                     indexSpaceElements, indexSpaceNodes, indexSapceDt,
                                                     exclusivePartitioning2, haloPartitioning2, overlapingPartitioning2,
                                                     allPartitioning2,
                                                     transitionToExclusive2, transitionToHalo2,
                                                     transitionToOverlappingInit2, transitionToOverlappingReduce2);
                master = laik_new_partitioning(laik_Master, world, laik_data_get_space(laikReductionAuxiliaryContainer), nullptr);
                all = laik_new_partitioning(laik_All, world, laik_data_get_space(laikReductionAuxiliaryContainer), nullptr);

                // migrate data for all the data structures
                laik_log(LAIK_LL_Info, "Redistributing data structures.");
                //Todo: This was also switched sign, see other TODO
                locDom->re_distribute_data_structures(shrinked_group, exclusivePartitioning2, haloPartitioning2,
                                                      overlapingPartitioning2, transitionToExclusive2,
                                                      transitionToHalo2, transitionToOverlappingInit2,
                                                      transitionToOverlappingReduce2,
                                                      isFaultToleranceActive(&ftOptions));
                laik_switchto_partitioning(laikReductionAuxiliaryContainer, all, LAIK_DF_None, LAIK_RO_None);

                // processes that are not part of the new (shrinked)
                // process group have to exit the main loop
                if (laik_myid(shrinked_group) == -1) {
                    laik_log(LAIK_LL_Info, "Rank %i eliminated. Good bye.", myRank);
                    break;
                }

                // repartition and geting the base pointer for only for dt
                laik_switchto_partitioning(laikDt, allPartitioning2, LAIK_DF_None, LAIK_RO_Min);
                laik_switchto_partitioning(laikDt, allPartitioning2, LAIK_DF_Preserve, LAIK_RO_Min);
                assert(laik_my_slicecount(laik_data_get_partitioning(laikDt)) == 1);
                laik_get_map_1d(laikDt, 0, (void **) &dt_base, &dt_count);

                // update the working process group in codes
                world = shrinked_group;
                std::cout << "New world size: " << laik_size(world) << std::endl;

                // update required values according to the new process group
                InitMeshDecomp(laik_size(world), laik_myid(world), &col, &row, &plane, &side);
                // update Domain according to the new process group
                locDom->re_init_domain(laik_size(world), col, row, plane, opts.nx,
                                       side, opts.numReg, opts.balance, opts.cost);

                double duration = getMPI_Wtime(&opts) - repart_start;

                laik_log((Laik_LogLevel) 2, "After repart\n");
                if (opts.quiet == 0)
                    printf("Repartition Done in %f s on Rank %d\n", duration, laik_myid(world));

                start2 = getMPI_Wtime(&opts);
                repart_iter = locDom->cycle();

                // remove the old partitioning and transition objects
                remove_partitionings_and_transitions(exclusivePartitioning, haloPartitioning, overlapingPartitioning,
                                                     allPartitioning,
                                                     transitionToExclusive, transitionToHalo,
                                                     transitionToOverlappingInit, transitionToOverlappingReduce);

                // update the working allPartitioning used later in
                // reduction for dt
                allPartitioning = allPartitioning2;

                free(removeList);

                TRACE_EVENT_END("RESTORE", "");
            } else if(failedCount == 0 && locDom->cycle() % ftOptions.checkpointFrequency == 0) {
                if(!ftOptions.delayCheckpointRelease) {
                    if(myRank == 0) {
                        laik_log(LAIK_LL_Info, "Freeing %zu checkpoints.\n", checkpoints.size());
                    }

                    for (auto & checkpoint : checkpoints) {
                        laik_checkpoint_free(checkpoint);
                    }
                    checkpoints.clear();
                }

                if(myRank == 0) {
                    laik_log(LAIK_LL_Info, "Creating checkpoints.");
                }
                TRACE_EVENT_START("CHECKPOINT", "");
                std::vector<Laik_Checkpoint*> newCheckpoints;
                locDom->createCheckpoints(newCheckpoints, ftOptions.redundancyCount, ftOptions.rotationDistance);
                TRACE_EVENT_END("CHECKPOINT", "");
                if(myRank == 0) {
                    laik_log(LAIK_LL_Info, "Finished creating checkpoints.");
                }

                if(ftOptions.delayCheckpointRelease) {
                    if(myRank == 0) {
                        laik_log(LAIK_LL_Info, "Freeing %zu checkpoints.\n", checkpoints.size());
                    }

                    for (auto & checkpoint : checkpoints) {
                        laik_checkpoint_free(checkpoint);
                    }
                    checkpoints.clear();
                }

                checkpoints = newCheckpoints;

                // Make sure error handler is installed so we can make use of the checkpoints
                laik_error_handler_set(inst, laik_failure_default_error_handler);
            }
        }
#endif

/*

       if (locDom->cycle() == 5){
            //locDom -> get_vdov().test_print();
           int dbg=1;
           if (laik_myid(world) == 0)
              {
                dbg = 0;
                while (!dbg) {}
              }
               
       }
*/


        TimeIncrement(*locDom, laikDt, allPartitioning, *dt_base);
        LagrangeLeapFrog(*locDom);

        if ((opts.showProg != 0) && (opts.quiet == 0) && (myRank == 0)) {
            printf("cycle = %d, time = %e, dt=%e\n",
                   locDom->cycle(), double(locDom->time()), double(locDom->deltatime()));
        }
    }

    TRACE_EVENT_START("FINALIZE", "");

    // Use reduced max elapsed time
    double elapsed_time;
#if USE_MPI
    elapsed_time = getMPI_Wtime(&opts) - start;
    double alternate_time = getMPI_Wtime(&opts) - start2;
#else
                                                                                                                            timeval end;
   gettimeofday(&end, NULL) ;
   elapsed_time = (double)(end.tv_sec - start.tv_sec) + ((double)(end.tv_usec - start.tv_usec))/1000000 ;
#endif
    double elapsed_timeG;
#if USE_MPI
    double elapsed_timeG2;
//    MPI_Reduce(&elapsed_time, &elapsed_timeG, 1, MPI_DOUBLE,
//               MPI_MAX, 0, MPI_COMM_WORLD);
//    MPI_Reduce(&alternate_time, &elapsed_timeG2, 1, MPI_DOUBLE,
//               MPI_MAX, 0, MPI_COMM_WORLD);
    elapsed_timeG = laik_reduce(elapsed_time, myRank, laikReductionAuxiliaryContainer, all, master);
    elapsed_timeG2 = laik_reduce(alternate_time, myRank, laikReductionAuxiliaryContainer, all, master);
#else
    elapsed_timeG = elapsed_time;
#endif

    // Write out final viz file */
    if (opts.viz) {
        DumpToVisit(*locDom, opts.numFiles, myRank, numRanks);
    }

    if ((myRank == 0) && (opts.quiet == 0)) {
        VerifyAndWriteFinalOutput(elapsed_timeG, *locDom, opts.nx, numRanks);
#ifdef USE_MPI
        if (opts.repart > 0) {
            printf("Time After Repartitioning: %f, per Iteration %f\n", elapsed_timeG2,
                   elapsed_timeG2 / (locDom->cycle() - repart_iter));
        }
#endif
    }

#if USE_MPI
    // free laik_objects
    if (exclusivePartitioning2 != nullptr)
        remove_partitionings_and_transitions(exclusivePartitioning2, haloPartitioning2, overlapingPartitioning2,
                                             allPartitioning2,
                                             transitionToExclusive2, transitionToHalo2, transitionToOverlappingInit2,
                                             transitionToOverlappingReduce2);
    else
        remove_partitionings_and_transitions(exclusivePartitioning, haloPartitioning, overlapingPartitioning,
                                             allPartitioning,
                                             transitionToExclusive, transitionToHalo, transitionToOverlappingInit,
                                             transitionToOverlappingReduce);

    laik_finalize(inst);
#endif

    TRACE_EVENT_END("FINALIZE", "");
    return 0;
}

double getMPI_Wtime(cmdLineOpts* opts) {
    if(opts->noTimer) {
        return 0;
    }
    return MPI_Wtime();
}

double laik_reduce(double input, Int_t myRank, Laik_Data *laikTimer, Laik_Partitioning *all, Laik_Partitioning *master) {
    double output;
    double* buffer;
    uint64_t count;
    laik_switchto_partitioning(laikTimer, all, LAIK_DF_None, LAIK_RO_None);
    assert(laik_my_slicecount(laik_data_get_partitioning(laikTimer)) == 1);
    laik_get_map_1d(laikTimer, 0, reinterpret_cast<void **>(&buffer), &count);
    if(count != 1) {
        std::cerr << "Laik reduce received incorrect partitioning size\n" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    buffer[0] = input;
    laik_switchto_partitioning(laikTimer, master, LAIK_DF_Preserve, LAIK_RO_Max);
    if(myRank == 0) {
        assert(laik_my_slicecount(laik_data_get_partitioning(laikTimer)) == 1);
        laik_get_map_1d(laikTimer, 0, reinterpret_cast<void **>(buffer), &count);
        if(count != 1) {
            std::cerr << "Laik reduce received incorrect partitioning size\n" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
        output = buffer[0];
    } else {
        output = -1;
    }
    return output;
}

bool parseFaultToleranceOptionsProxy(int argc, char **argv, int *arg, int rank, FaultToleranceOptions *ftOptions) {
    return parseFaultToleranceOptions(argc, argv, arg, rank, ftOptions);
}