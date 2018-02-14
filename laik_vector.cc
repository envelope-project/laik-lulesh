#include <laik_vector.h>
#include <laik_partitioners.h>
#include <lulesh.h>
#include<limits.h>

laik_vector::laik_vector(Laik_Instance* inst, Laik_Group* world){
    this -> inst = inst;
    this -> world = world;


    int numRanks = laik_size(world);
    int myRank = laik_myid(world);
    int col, row, plane, side;
    InitMeshDecomp(numRanks, myRank, &col, &row, &plane, &side);

    b=1;
    f=1;
    d=1;
    u=1;
    l=1;
    r=1;

    if (col==0) {
        d=0;
    }

    if (col==side-1) {
        u=0;
    }

    if (row==0) {
        l=0;
    }

    if (row==side-1) {
        r=0;
    }

    if (plane==0) {
        b=0;
    }

    if (plane==side-1) {
        f=0;
    }

    state=0;
    zero=0;
}

void laik_vector::test_print(){
    double *base;
    uint64_t count;
    int nSlices = laik_phase_my_slicecount(haloPartitioning);
    for (size_t s = 0; s < nSlices; s++)
    {
        laik_map_def(data, s, (void**) &base, &count);
        for (uint64_t i = 0; i < count; i++)
        {
            laik_log(Laik_LogLevel(2),"%f\n", base[i]);
        }
        laik_log(Laik_LogLevel(2),"\n");
    }
}

// ////////////////////////////////////////////////////////////////////////
// implementation of laik_vector with halo partitioning (node partitioning)
// ////////////////////////////////////////////////////////////////////////
laik_vector_halo::laik_vector_halo(Laik_Instance *inst,
                                   Laik_Group *world):laik_vector(inst,world){}

void laik_vector_halo::resize(int count){

    size = count;
    int halo_depth = 1;
    indexSpace = laik_new_space_1d(inst, size);
    exclusivePartitioning = laik_new_accessphase(world,
                                                 indexSpace,
                                                 exclusive_partitioner(), 0);
    haloPartitioning = laik_new_accessphase(world, indexSpace,
                                            overlaping_partitioner(halo_depth), 0);
    overlapingPartitioning = 0;
    data = laik_new_data(indexSpace, laik_Double);

    // go through the slices to just allocate the memory
    uint64_t cnt;
    double* base;

    laik_switchto_phase(data, exclusivePartitioning, LAIK_DF_CopyOut);
    int nSlices = laik_phase_my_slicecount(exclusivePartitioning);
    for (int n = 0; n < nSlices; ++n)
    {
       laik_map_def(data, n, (void **)&base, &cnt);
    }
    laik_switchto_phase(data, haloPartitioning, LAIK_DF_CopyIn);

    this -> count = cnt;
}

double& laik_vector_halo::operator [](int idx){
    uint64_t cnt;
    double* base;

    int i=0;
    int j=0;
    int k=0;
    int index=0;
    int slice=-1;
    int side =-1;
    int numElem=count*count*count;

    Index_t ghostIdx[6] ;  // offsets to ghost locations

    for (Index_t i=0; i<6; ++i) {
      ghostIdx[i] = INT_MAX;
    }

    Int_t pidx = numElem ;
    if (b) {
      ghostIdx[0] = pidx ;
      pidx += count*count ;
    }
    if (f) {
      ghostIdx[1] = pidx ;
      pidx += count*count ;
    }
    if (d) {
      ghostIdx[2] = pidx ;
      pidx += count*count ;
    }
    if (u) {
      ghostIdx[3] = pidx ;
      pidx += count*count ;
    }
    if (l) {
      ghostIdx[4] = pidx ;
      pidx += count*count ;
    }
    if (r) {
      ghostIdx[5] = pidx ;
    }

    /*
    for (int i = 0; i < 6; ++i) {
            laik_log(Laik_LogLevel(2), "%d", ghostIdx[i]);
    }
    */

    for (int i = 0; i < 6; ++i) {
        if  (idx>=ghostIdx[i]){
            side =i;
        }
    }
    // tested OK
    if (state) {
        i = idx % count;
        slice = idx/count;
    }
    else{
        // tested OK
        if (idx < numElem) {
            i = idx%count;
            j = (idx/count)%count;
            k = idx/(count*count);
            slice = (count+d+u)*(k+b)+(j+d);
            i += l;
        }

        // back
        // tested OK
        else if (side==0) {
            index =idx-ghostIdx[side];
            i = index%count;
            j = index/count;
            k = -1;
            slice = (count+d+u)*(k+b)+(j+d);
            i += l;
        }
        // front
        // tested OK
        else if (side==1) {
            index =idx-ghostIdx[side];
            i = index%count;
            j = index/count;
            k = count;
            slice = (count+d+u)*(k+b)+(j+d);
            i += l;
        }
        //down
        //tested OK
        else if (side==2) {
            index =idx-ghostIdx[side];
            i = index%count;
            j = -1;
            k = index/count;
            slice = (count+d+u)*(k+b)+(j+d);
            i += l;
        }
        //up
        //tested OK
        else if (side==3) {
            index =idx-ghostIdx[side];
            i = index%count;
            j = count;
            k = index/count;
            slice = (count+d+u)*(k+b)+(j+d);
            i += l;
        }
        //left
        //tested OK
        else if (side==4) {
            index =idx-ghostIdx[side];
            i = -1;
            j = index%count;
            k = index/count;
            slice = (count+d+u)*(k+b)+(j+d);
            i += l;
        }
        //right
        //tested OK
        else if (side==5) {
            index =idx-ghostIdx[side];
            i = count;
            j = index%count;
            k = index/count;
            slice = (count+d+u)*(k+b)+(j+d);
            i += l;
        }

    }

    //laik_log(Laik_LogLevel(2),"state: %d, idx: %d, side: %d,
    //slice: %d, i:%d", state, idx, side, slice, i);

    if (slice>-1) {
        laik_map_def(data, slice, (void **)&base, &cnt);
        return base[i];
    }

    return zero;

}

void laik_vector_halo::switch_to_write_phase(){
    laik_switchto_phase(data, exclusivePartitioning, LAIK_DF_CopyOut);
    state=1;
}

void laik_vector_halo::switch_to_read_phase(){
    laik_switchto_phase(data, haloPartitioning, LAIK_DF_CopyIn);
    state=0;
}

// ///////////////////////////////////////////////////////////
// implementation of laik_vector with overlapping partitioning
// (node partitioning)
// ///////////////////////////////////////////////////////////

laik_vector_overlapping::laik_vector_overlapping(Laik_Instance *inst,
                                                 Laik_Group *world):laik_vector(inst,world){}

void laik_vector_overlapping::resize(int count){
    size = count;
    int halo_depth = 1;
    indexSpace = laik_new_space_1d(inst, size);
    exclusivePartitioning = 0;
    haloPartitioning = 0;
    overlapingPartitioning =laik_new_accessphase(world,
                                 indexSpace,
                                 overlaping_reduction_partitioner(halo_depth), 0);
    data = laik_new_data(indexSpace, laik_Double);

    // go through the slices to just allocate the memory
    laik_switchto_phase(data, overlapingPartitioning, Laik_DataFlow
                        ( LAIK_DF_Init | LAIK_DF_ReduceOut | LAIK_DF_Sum ) );
    uint64_t cnt;
    double* base;
    int nSlices = laik_phase_my_slicecount(overlapingPartitioning);
    for (int n = 0; n < nSlices; ++n)
    {
       laik_map_def(data, n, (void **)&base, &cnt);
    }
    laik_switchto_phase(data, overlapingPartitioning, LAIK_DF_CopyIn);

    this -> count = cnt;

}

double& laik_vector_overlapping::operator [](int idx){
    uint64_t cnt;
    double* base;

    int i = (idx % (count*count) ) % count;
    int l_idx = i;
    int slice = idx/count;
    laik_map_def(data, slice, (void **)&base, &cnt);

    return base[l_idx];
}

void laik_vector_overlapping::switch_to_write_phase(){
    laik_switchto_phase(data, overlapingPartitioning, Laik_DataFlow
                        ( LAIK_DF_ReduceOut | LAIK_DF_Sum ) );
}

void laik_vector_overlapping::switch_to_read_phase(){
    laik_switchto_phase(data, overlapingPartitioning, LAIK_DF_CopyIn);
}
