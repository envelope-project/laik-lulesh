#include <laik_partitioners.h>
#include <lulesh.h>

void runExclusivePartitioner(Laik_Partitioner* pr,
                                   Laik_Partitioning* p,
                             Laik_Partitioning* oldP)
{
    Laik_Group* world = laik_partitioning_get_group(p);
    int numRanks = laik_size(world);
    int myRank = laik_myid(world);
    int col, row, plane, side;
    InitMeshDecomp(numRanks, myRank, &col, &row, &plane, &side);

    // get the size of the
    Laik_Space* space = laik_partitioning_get_space(p);
    const Laik_Slice* slice = laik_space_asslice(space);
    int edgeElems= (int) (cbrt( (slice->to.i[0]+1) / numRanks ) + 0.1 );

    int Nx=edgeElems;
    int Ny=edgeElems;
    int Nz=edgeElems;
    int Rx= side;
    int Ry= side;
    int Rz= side;
    int Lx = Rx*Nx;
    int Ly = Ry*Ny;
    int Lz = Rz*Nz;
    int Pxy= Lx*Ly;
    int Pxz= Lx*Lz;
    int Pyz= Ly*Lz;

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
                r = rx + ry*Rx + rz*Rx*Ry; // task number
                // loop over z and x  to create the slices in the
                // partitioning
                for (int ny = 0; ny < Ny; ny++)
                {
                    for (int nz = 0; nz < Nz; nz++)
                    {
                        nx=0;
                        slc.from.i[0]=nx + Lx*ny + Pxy*nz +
                                rx*Nx + ry*Lx*Ny + Pxy*Nz*rz;
                        nx=Nx;
                        slc.to.i[0]=nx + Lx*ny + Pxy*nz +
                                rx*Nx + ry*Lx*Ny + Pxy*Nz*rz;
                        laik_append_slice(p,r,&slc,0,0);
                    }
                }
            }
        }
    }

}

Laik_Partitioner* exclusive_partitioner()
{
    return laik_new_partitioner("exclusive", runExclusivePartitioner, 0,
                                LAIK_PF_Merge);
}

void runOverlapingPartitioner(Laik_Partitioner* pr,
                                   Laik_Partitioning* p,
                              Laik_Partitioning* oldP)
{
    Laik_Group* world = laik_partitioning_get_group(p);
    int numRanks = laik_size(world);
    int myRank = laik_myid(world);
    int col, row, plane, side;
    InitMeshDecomp(numRanks, myRank, &col, &row, &plane, &side);

    // get the size of the
    Laik_Space* space = laik_partitioning_get_space(p);
    const Laik_Slice* slice = laik_space_asslice(space);
    int edgeElems= (int) ( cbrt( (slice->to.i[0]+1) / numRanks) + 0.1 );

    // the number of halos in each boundary
    int d = *(int*) laik_partitioner_data(pr);

    if (d>edgeElems)
    {
        laik_log (LAIK_LL_Error, "number of halo is too large!"
                                 " fix your application");
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
    int Lz = Rz*Nz;
    int Pxy= Lx*Ly;
    int Pxz= Lx*Lz;
    int Pyz= Ly*Lz;

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
                r = rx + ry*Rx + rz*Rx*Ry; // task number
                // loop over y and z  to create the slices in the
                // partitioning
                for (int ny = ((ry==0) ?0:-d); ny < ((ry==Ry-1)?Ny:Ny+d) ; ny++)
                {
                    for (int nz = ( (rz==0)?0:-d ) ;
                         nz < ( (rz==Rz-1)?Nz:Nz+d ) ; nz++)
                    {
                        nx= (rx==0)?0:-d;
                        slc.from.i[0]=nx + Lx*ny + Pxy*nz +
                                rx*Nx + ry*Lx*Ny + Pxy*Nz*rz;
                        nx= (rx==Rx-1)?Nx:Nx+d;
                        slc.to.i[0]=nx + Lx*ny + Pxy*nz +
                                rx*Nx + ry*Lx*Ny + Pxy*Nz*rz;
                        laik_append_slice(p,r,&slc,0,0);
                    }
                }
            }
        }
    }

}

Laik_Partitioner* overlaping_partitioner(int &depth)
{
    //void* data = (void*) &depth;
    return laik_new_partitioner("overlaping", runOverlapingPartitioner,
                                (void*) &depth, LAIK_PF_Merge);
}

void runOverlapingReductionPartitioner(Laik_Partitioner* pr,
                                       Laik_Partitioning* p,
                                       Laik_Partitioning* oldP)
{

    Laik_Group* world = laik_partitioning_get_group(p);
    int numRanks = laik_size(world);
    int myRank = laik_myid(world);
    int col, row, plane, side;
    InitMeshDecomp(numRanks, myRank, &col, &row, &plane, &side);

    // get the size of the
    Laik_Space* space = laik_partitioning_get_space(p);
    const Laik_Slice* slice = laik_space_asslice(space);
    int edgeNodes= (int) ( (cbrt( (slice->to.i[0]+1)) -1)/
            cbrt(numRanks) + 1   + 0.1);

    //laik_log ((Laik_LogLevel)2, "elems: %d", edgeNodes);

    // the number of halos in each boundary
    int d = *(int*) laik_partitioner_data(pr);

    if (d>edgeNodes)
    {
        laik_log (LAIK_LL_Error, "number of halo is too large!"
                                 " fix your application");
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
                //r = ry + rx*Ry + rz*Rx*Ry; // (yxz)
                r = rx + ry*Rx + rz*Rx*Ry; // (xyz)
                //r = rx + rz*Rx + ry*Rx*Rz; // (xzy)
                //r = ry + rz*Ry + rx*Rz*Ry; // (yzx)
                //r = rz + ry*Rz + rx*Rz*Ry; // (zyx)
                //r = rz + rx*Rz + ry*Rx*Rz; // (zxy)

                // loop over y and z to create the slices in the
                // partitioning
                for (int ny = 0 ; ny < Ny; ny++)
                {
                    for (int nz = 0 ; nz < Nz; nz++)
                    {
                        nx=0;
                        slc.from.i[0]=
                                nx + Lx*ny + Pxy*nz +  rx*(Nx-1) + ry*Lx*(Ny-1) + rz*Pxy*(Nz-1);
                        nx=Nx;
                        slc.to.i[0]= nx + Lx*ny + Pxy*nz +  rx*(Nx-1) + ry*Lx*(Ny-1) + rz*Pxy*(Nz-1);
                        laik_append_slice(p,r,&slc,0,0);
                    }
                }
            }
        }
    }

}

Laik_Partitioner* overlaping_reduction_partitioner(int &depth)
{
    //void* data = (void*) &depth;
    return laik_new_partitioner("overlapingReduction",
                                runOverlapingReductionPartitioner,
                                (void*) &depth, LAIK_PF_Merge);
}
