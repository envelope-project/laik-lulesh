#include <laik_port.h>
#include <lulesh.h>

// directly copied from the lulesh_init.cc -> BuildMesh
// signature is modified to get rid of the "domain"
int * build_element_corner_list(int edgeElems, int edgeNodes, int m_rowLoc, int m_colLoc, int m_planeLoc)
{
  // NOTE: this function does not initialize the mesh 
  // it only creates local node indexes for local
  // nodes
  int numElem = edgeElems*edgeElems*edgeElems;
  int *localNodeList = (int*) malloc (sizeof(int)*numElem*8);
  int *localNode = 0;
  int nidx = 0;
  int zidx = 0;
  for (int plane = 0; plane < edgeElems; ++plane)
  {
      for (int row = 0; row < edgeElems; ++row)
      {
          for (int col = 0; col < edgeElems; ++col)
          {
              localNode = localNodeList + 8*zidx ;
              localNode[0] = nidx;
              localNode[1] = nidx + 1;
              localNode[2] = nidx + edgeNodes + 1;
              localNode[3] = nidx + edgeNodes;
              localNode[4] = nidx + edgeNodes * edgeNodes;
              localNode[5] = nidx + edgeNodes * edgeNodes + 1;
              localNode[6] = nidx + edgeNodes * edgeNodes + edgeNodes + 1;
              localNode[7] = nidx + edgeNodes * edgeNodes + edgeNodes;
              ++zidx;
              ++nidx;
                laik_log((Laik_LogLevel)1,"elem:%d, corners:%d, %d, %d, %d, %d, %d, %d, %d", zidx, localNode[0], localNode[1], localNode[2], localNode[3], localNode[4], localNode[5] 
                , localNode[6], localNode[7]);
          }
          ++nidx;
      }
      nidx += edgeNodes;
  }
  return localNodeList;
}

//TODO free the list in the program
void free_local_corner_list(int* list){
    free(list);
}

void runExclusivePartitioner(Laik_Partitioner* pr,
                                   Laik_BorderArray* ba, Laik_BorderArray* otherBA)
{
    Laik_Group* world = laik_borderarray_getgroup(ba);
    int numRanks = laik_size(world);
    int myRank = laik_myid(world);
    int col, row, plane, side;
    InitMeshDecomp(numRanks, myRank, &col, &row, &plane, &side);
    
    // get the size of the
    Laik_Space* space = laik_borderarray_getspace(ba);
    const Laik_Slice* slice = laik_space_getslice(space);
    int edgeElems= cbrt( (slice->to.i[0]+1) / numRanks);

    int Nx=edgeElems;
    int Ny=edgeElems;
    int Nz=edgeElems;
    int Rx= side;
    int Ry= side;
    int Rz= side;
    int Lx = Rx*Nx;
    int Ly = Ry*Ny;
    int Pxy= Lx*Ly;

    // sine all the tasks run the same partitioning algorithm
    // we should loop over all the tasks and not just this 
    // task
    Laik_Slice slc;
    int r=0;
    int nx=0;
    for (int rz = 0; rz < Rz; rz++)
    {
        for (int ry = 0; ry < Ry; ry++)
        {
            for (int rx = 0; rx < Rx; rx++)
            {
                r = ry + rx*Ry + rz*Rx*Ry; // task number
                // loop over y and z  to create the slices in the
                // partitioning
                for (int ny = 0; ny < Ny; ny++)
                {
                    for (int nz = 0; nz < Nz; nz++)
                    {
                        nx=0;
                        slc.from.i[0]=nx + Lx*ny + Pxy*nz +  rx*Nx + ry*Lx*Ny + Pxy*Nz*rz;
                        nx=Nx;
                        slc.to.i[0]=nx + Lx*ny + Pxy*nz  +  rx*Nx + ry*Lx*Ny + Pxy*Nz*rz;
                        laik_append_slice(ba,r,&slc,0,0);
                    }
                }
            }
        }    
    }


    //int edgeNodes=edgeElems+1;
    //int numElem = edgeElems*edgeElems*edgeElems;

    // later on will be used for element partitioner
    //int* nodeList = build_element_corner_list(edgeElems, edgeNodes, row, col, plane);
    
    
    //slc.from.i[0]=0;
    //slc.to.i[0]=numElem*numRanks;
    //laik_append_slice(ba, 0, &slc, 0, 0);

    //free_local_corner_list(nodeList);

}

Laik_Partitioner* exclusive_partitioner()
{
    return laik_new_partitioner("exclusive", runExclusivePartitioner, 0, LAIK_PF_Merge);
}

void runOverlapingPartitioner(Laik_Partitioner* pr,
                                   Laik_BorderArray* ba, Laik_BorderArray* otherBA)
{
    Laik_Group* world = laik_borderarray_getgroup(ba);
    int numRanks = laik_size(world);
    int myRank = laik_myid(world);
    int col, row, plane, side;
    InitMeshDecomp(numRanks, myRank, &col, &row, &plane, &side); // could be useful for distributed partitioner
    
    // get the size of the
    Laik_Space* space = laik_borderarray_getspace(ba);
    const Laik_Slice* slice = laik_space_getslice(space);
    int edgeElems= cbrt( (slice->to.i[0]+1) / numRanks);

    // the number of halos in each boundary
    int d = *(int*) laik_partitioner_data(pr);

    if (d>edgeElems)
    {
        laik_log (LAIK_LL_Error, "number of halo is too large! fix your application");
        exit(0);
    }

    int Nx=edgeElems;
    int Ny=edgeElems;
    int Nz=edgeElems;
    int Rx= side;
    int Ry= side;
    int Rz= side;
    int Lx = Rx*Nx;
    int Ly = Ry*Ny;
    int Pxy= Lx*Ly;

    // sine all the tasks run the same partitioning algorithm
    // we should loop over all the tasks and not just this 
    // task
    Laik_Slice slc;
    int r=0;
    int nx=0;
    for (int rz = 0; rz < Rz; rz++)
    {
        for (int ry = 0; ry < Ry; ry++)
        {
            for (int rx = 0; rx < Rx; rx++)
            {
                r = ry + rx*Ry + rz*Rx*Ry; // task number
                // loop over y and z  to create the slices in the
                // partitioning
                for (int ny = ((ry==0) ?0:-d); ny < ((ry==Ry-1)?Ny:Ny+d) ; ny++)
                //for (int ny = 0 ; ny < Ny; ny++)
                {
                    for (int nz = ( (rz==0)?0:-d ) ; nz < ( (rz==Rz-1)?Nz:Nz+d ) ; nz++)
                    //for (int nz = 0 ; nz < Nz; nz++)
                    {
                        nx= (rx==0)?0:-d;
                        slc.from.i[0]=nx + Lx*ny + Pxy*nz +  rx*Nx + ry*Lx*Ny + Pxy*Nz*rz;
                        //laik_log((Laik_LogLevel)1,"rank:%d, from %d, ", r, slc.from.i[0]);
                        nx= (rx==Rx-1)?Nx:Nx+d;
                        slc.to.i[0]=nx + Lx*ny + Pxy*nz  +  rx*Nx + ry*Lx*Ny + Pxy*Nz*rz;
                        laik_append_slice(ba,r,&slc,0,0);
                    }
                }
            }
        }    
    }


    //int edgeNodes=edgeElems+1;
    //int numElem = edgeElems*edgeElems*edgeElems;

    // later on will be used for element partitioner
    //int* nodeList = build_element_corner_list(edgeElems, edgeNodes, row, col, plane);
    
    
    //slc.from.i[0]=0;
    //slc.to.i[0]=numElem*numRanks;
    //laik_append_slice(ba, 0, &slc, 0, 0);

    //free_local_corner_list(nodeList);

}

Laik_Partitioner* overlaping_partitioner(int &depth)
{
    //void* data = (void*) &depth;
    return laik_new_partitioner("overlaping", runOverlapingPartitioner, (void*) &depth, LAIK_PF_Merge);
}

void runOverlapingReductionPartitioner(Laik_Partitioner* pr, Laik_BorderArray* ba, Laik_BorderArray* otherBA){

    Laik_Group* world = laik_borderarray_getgroup(ba);
    int numRanks = laik_size(world);
    int myRank = laik_myid(world);
    int col, row, plane, side;
    InitMeshDecomp(numRanks, myRank, &col, &row, &plane, &side);

    // get the size of the
    Laik_Space* space = laik_borderarray_getspace(ba);
    const Laik_Slice* slice = laik_space_getslice(space);
    int edgeNodes= (cbrt( (slice->to.i[0]+1)) -1)/cbrt(numRanks) + 1 ;

    //laik_log ((Laik_LogLevel)2, "elems: %d", edgeNodes);

    // the number of halos in each boundary
    int d = *(int*) laik_partitioner_data(pr);

    if (d>edgeNodes)
    {
        laik_log (LAIK_LL_Error, "number of halo is too large! fix your application");
        exit(0);
    }

    int Nx=edgeNodes;
    int Ny=edgeNodes;
    int Nz=edgeNodes;
    int Rx= side;
    int Ry= side;
    int Rz= side;
    int Lx = Rx*(Nx-1)+1;
    int Ly = Ry*(Ny-1)+1;
    int Pxy= Lx*Ly;


    // sine all the tasks run the same partitioning algorithm
    // we should loop over all the tasks and not just this
    // task
    Laik_Slice slc;
    int r=0;
    int nx=0;
    for (int rz = 0; rz < Rz; rz++)
    {
        for (int ry = 0; ry < Ry; ry++)
        {
            for (int rx = 0; rx < Rx; rx++)
            {
                r = ry + rx*Ry + rz*Rx*Ry; // task number
                // loop over y and z  to create the slices in the
                // partitioning
                //for (int ny = ((ry==0) ?0:-d); ny < ((ry==Ry-1)?Ny:Ny+d) ; ny++)
                for (int ny = 0 ; ny < Ny; ny++)
                {
                    //for (int nz = ( (rz==0)?0:-d ) ; nz < ( (rz==Rz-1)?Nz:Nz+d ) ; nz++)
                    for (int nz = 0 ; nz < Nz; nz++)
                    {
                        nx=0;
                        slc.from.i[0]=nx + Lx*ny + Pxy*nz +  ry*(Nx-1) + rx*Lx*(Ny-1) + rz*Pxy*(Nz-1);
                        //laik_log((Laik_LogLevel)2,"rank:%d, from %d, ", r, slc.from.i[0]);
                        nx=Nx;
                        slc.to.i[0]=nx + Lx*ny + Pxy*nz +  ry*(Nx-1) + rx*Lx*(Ny-1) + rz*Pxy*(Nz-1);
                        laik_append_slice(ba,r,&slc,0,0);
                    }
                }
            }
        }
    }

}

Laik_Partitioner* overlaping_reduction_partitioner(int &depth)
{
    //void* data = (void*) &depth;
    return laik_new_partitioner("overlapingReduction", runOverlapingReductionPartitioner, (void*) &depth, LAIK_PF_Merge);
}
