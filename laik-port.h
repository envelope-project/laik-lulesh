#include "lulesh.h"
#include "laik.h"

// directly copied from the lulesh_init.cc -> BuildMesh
// signature is modified to get rid of the "domain"
Index_t * build_element_corner_list(Int_t edgeElems, Int_t edgeNodes, Int_t m_rowLoc, Int_t m_colLoc, Int_t m_planeLoc)
{
  // NOTE: this function does not initialize the mesh 
  // it only creates local node indexes for local
  // nodes
  Index_t numElem = edgeElems*edgeElems*edgeElems;
  Index_t *localNodeList = (Index_t*) malloc (sizeof(Index_t)*numElem*8);
  Index_t *localNode = 0;
  Index_t nidx = 0;
  Index_t zidx = 0;
  for (Index_t plane = 0; plane < edgeElems; ++plane)
  {
      for (Index_t row = 0; row < edgeElems; ++row)
      {
          for (Index_t col = 0; col < edgeElems; ++col)
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

void runElementPartitionerExclusive(Laik_Partitioner* pr,
                                   Laik_BorderArray* ba, Laik_BorderArray* otherBA)
{
    Laik_Group* world = laik_borderarray_getgroup(ba);
    Int_t numRanks = laik_size(world);
    Int_t myRank = laik_myid(world);
    Int_t col, row, plane, side;
    InitMeshDecomp(numRanks, myRank, &col, &row, &plane, &side);
    Index_t edgeElems= 4;

    Index_t Nx=edgeElems;
    Index_t Ny=edgeElems;
    Index_t Nz=edgeElems;
    Index_t Rx= side;
    Index_t Ry= side;
    Index_t Rz= side;
    Index_t Lx = Rx*Nx;
    Index_t Ly = Ry*Ny;
    Index_t Pxy= Lx*Ly;

    // sine all the tasks run the same partitioning algorithm
    // we should loop over all the tasks and not just this 
    // task
    Laik_Slice slc;
    Index_t r=0;
    Index_t nx=0;
    for (Index_t rz = 0; rz < Rz; rz++)
    {
        for (Index_t ry = 0; ry < Ry; ry++)
        {
            for (Index_t rx = 0; rx < Rx; rx++)
            {
                r = ry + rx*Ry + rz*Rx*Ry; // task number
                // loop over y and z  to create the slices in the
                // partitioning
                for (Index_t ny = 0; ny < Ny; ny++)
                {
                    for (Index_t nz = 0; nz < Nz; nz++)
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


    //Index_t edgeNodes=edgeElems+1;
    //Index_t numElem = edgeElems*edgeElems*edgeElems;

    // later on will be used for element partitioner
    //Int_t* nodeList = build_element_corner_list(edgeElems, edgeNodes, row, col, plane);
    
    
    //slc.from.i[0]=0;
    //slc.to.i[0]=numElem*numRanks;
    //laik_append_slice(ba, 0, &slc, 0, 0);

    //free_local_corner_list(nodeList);

}

Laik_Partitioner* element_partitioner_exclusive()
{
    return laik_new_partitioner("exclusive-element", runElementPartitionerExclusive, 0, LAIK_PF_Merge);
}

void runElementPartitionerOverlaping(Laik_Partitioner* pr,
                                   Laik_BorderArray* ba, Laik_BorderArray* otherBA)
{
    Laik_Group* world = laik_borderarray_getgroup(ba);
    Int_t numRanks = laik_size(world);
    Int_t myRank = laik_myid(world);
    Int_t col, row, plane, side;
    InitMeshDecomp(numRanks, myRank, &col, &row, &plane, &side);
    Index_t edgeElems= 4;

    Index_t Nx=edgeElems;
    Index_t Ny=edgeElems;
    Index_t Nz=edgeElems;
    Index_t Rx= side;
    Index_t Ry= side;
    Index_t Rz= side;
    Index_t Lx = Rx*Nx;
    Index_t Ly = Ry*Ny;
    Index_t Pxy= Lx*Ly;

    // sine all the tasks run the same partitioning algorithm
    // we should loop over all the tasks and not just this 
    // task
    Laik_Slice slc;
    Index_t r=0;
    Index_t nx=0;
    for (Index_t rz = 0; rz < Rz; rz++)
    {
        for (Index_t ry = 0; ry < Ry; ry++)
        {
            for (Index_t rx = 0; rx < Rx; rx++)
            {
                r = ry + rx*Ry + rz*Rx*Ry; // task number
                // loop over y and z  to create the slices in the
                // partitioning
                for (Index_t ny = ((ry==0) ?0:-1); ny < ((ry==Ry-1)?Ny:Ny+1) ; ny++)
                //for (Index_t ny = 0 ; ny < Ny; ny++)
                {
                    for (Index_t nz = ( (rz==0)?0:-1 ) ; nz < ( (rz==Rz-1)?Nz:Nz+1 ) ; nz++)
                    //for (Index_t nz = 0 ; nz < Nz; nz++)
                    {
                        nx= (rx==0)?0:-1;
                        slc.from.i[0]=nx + Lx*ny + Pxy*nz +  rx*Nx + ry*Lx*Ny + Pxy*Nz*rz;
                        laik_log((Laik_LogLevel)2,"rank:%d, from %d, ", r, slc.from.i[0]);
                        nx= (rx==Rx-1)?Nx:Nx+1;
                        slc.to.i[0]=nx + Lx*ny + Pxy*nz  +  rx*Nx + ry*Lx*Ny + Pxy*Nz*rz;
                        laik_append_slice(ba,r,&slc,0,0);
                    }
                }
            }
        }    
    }


    //Index_t edgeNodes=edgeElems+1;
    //Index_t numElem = edgeElems*edgeElems*edgeElems;

    // later on will be used for element partitioner
    //Int_t* nodeList = build_element_corner_list(edgeElems, edgeNodes, row, col, plane);
    
    
    //slc.from.i[0]=0;
    //slc.to.i[0]=numElem*numRanks;
    //laik_append_slice(ba, 0, &slc, 0, 0);

    //free_local_corner_list(nodeList);

}

Laik_Partitioner* element_partitioner_overlaping()
{
    return laik_new_partitioner("overlaping-element", runElementPartitionerOverlaping, 0, LAIK_PF_Merge);
}