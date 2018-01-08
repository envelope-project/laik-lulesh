#include <laik_vector.h>
#include <laik_port.h>
#include <lulesh.h>

//laik_vector::laik_vector(){}

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
        f=0;
    }

    if (plane==side-1) {
        b=0;
    }

    state=0;
}

void laik_vector::resize(int count){

    size = count;
    int halo_depth = 1;
    indexSpace = laik_new_space_1d(inst, size);
    exclusivePartitioning = laik_new_partitioning(world, indexSpace, exclusive_partitioner(), 0);
    haloPartitioning = laik_new_partitioning(world, indexSpace, overlaping_partitioner(halo_depth), 0);
    overlapingPartitioning = laik_new_partitioning(world, indexSpace, overlaping_reduction_partitioner(halo_depth), 0);
    data = laik_new_data(world, indexSpace, laik_Double);

    // go through the slices to just allocate the memory
    uint64_t cnt;
    double* base;

    laik_switchto(data, exclusivePartitioning, LAIK_DF_CopyOut);
    int nSlices = laik_my_slicecount(exclusivePartitioning);
    for (int n = 0; n < nSlices; ++n)
    {
       laik_map_def(data, n, (void **)&base, &cnt);
    }
    laik_switchto(data, haloPartitioning, LAIK_DF_CopyIn);

    this -> count = cnt;
    //laik_switchto(data, exclusivePartitioning, LAIK_DF_CopyOut);

    /*
    laik_switchto(data, haloPartitioning, LAIK_DF_CopyIn);
    nSlices = laik_my_slicecount(haloPartitioning);
    for (int n = 0; n < nSlices; ++n)
    {
       laik_map_def(data, n, (void **)&base, &cnt);
    }
    laik_switchto(data, haloPartitioning, LAIK_DF_CopyIn);
    */

    //laik_log(Laik_LogLevel(2), "here");

}

double& laik_vector::operator [](int idx){
    uint64_t cnt;
    double* base;

    int i=0;
    int index=0;
    int slice=0;

    if (state) {
        i = (idx % (count*count) ) % count;
        slice = idx/count;
    }
    else{

        if (idx < count*count*count) {
            i = idx%count;
            int s = (idx/count)%count;
            int k = idx/(count*count);
            i += l;
            slice = (count+d+u)*(d+k)+(s+d);
        }

        if ( ( idx >= count*count*count && idx < (count*count*count + count*count) ) && (f) ) {
            index =idx-count*count*count;
            i = index%count;
            slice = index/count;
            i += l;
            slice+=d;
        }

        if ( (idx >= (count*count*count + count*count) && idx < (count*count*count + 2*count*count) ) && (b)) {
            index=idx-(count*count*count + count*count);
            i = index%count;
            slice = index/count;
            i += l;
            slice += (count+u+d)*(count+f)+d;
        }


        if ( (idx >= (count*count*count + 2*count*count) && idx < (count*count*count + 3*count*count) ) && (d) ) {
            index=idx-(count*count*count + 2*count*count);
            i = index%count;
            slice = index/count;
            i += l;
            slice = (slice+f)*(count+d+u);
        }


        if ( (idx >= (count*count*count + 3*count*count) && idx < (count*count*count + 4*count*count) ) && (u) ) {
            index=idx-(count*count*count + 3*count*count);
            i = index%count;
            slice = index/count;
            i += l;
            slice = (slice+f+1)*(count+d+u)-1;
        }

        if ( (idx >= (count*count*count + 4*count*count) && idx < (count*count*count + 5*count*count) ) && (l) ){
            index=idx-(count*count*count + 4*count*count);
            int s = index%count;
            int k = index/count;
            slice = (count+d+u)*(d+k)+(s+d);
            i=0;
        }

        if ( (idx >= (count*count*count + 5*count*count) && idx < (count*count*count + 6*count*count) ) && (r) ){
            index=idx-(count*count*count + 5*count*count);
            int s = index%count;
            int k = index/count;
            slice = (count+d+u)*(d+k)+(s+d);
            i=count+l+r-1;
        }

    }

    //laik_log(Laik_LogLevel(2),"%d %d",slice, i);

    laik_map_def(data, slice, (void **)&base, &cnt);

    return base[i];
}

void laik_vector::switch_to_exclusive_partitioning(){
    laik_switchto(data, exclusivePartitioning, LAIK_DF_CopyOut);
    state=1;
}

void laik_vector::switch_to_halo_partitioning(){
    laik_switchto(data, haloPartitioning, LAIK_DF_CopyIn);
    state=0;
}

void laik_vector::test_print(){
    double *base;
    uint64_t count;
    int nSlices = laik_my_slicecount(haloPartitioning);
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

laik_vector_halo::laik_vector_halo(Laik_Instance *inst, Laik_Group *world):
    laik_vector(inst,world){}

laik_vector_overlapping::laik_vector_overlapping(Laik_Instance *inst, Laik_Group *world):
    laik_vector(inst,world){}


void laik_vector_overlapping::resize(int count){

    size = count;
    int halo_depth = 1;
    indexSpace = laik_new_space_1d(inst, size);
    exclusivePartitioning = 0;
    haloPartitioning = 0;
    overlapingPartitioning = laik_new_partitioning(world, indexSpace, overlaping_reduction_partitioner(halo_depth), 0);
    data = laik_new_data(world, indexSpace, laik_Double);

    // go through the slices to just allocate the memory
    laik_switchto(data, overlapingPartitioning, Laik_DataFlow ( LAIK_DF_Init | LAIK_DF_ReduceOut | LAIK_DF_Sum ) );
    uint64_t cnt;
    double* base;
    int nSlices = laik_my_slicecount(overlapingPartitioning);
    for (int n = 0; n < nSlices; ++n)
    {
       laik_map_def(data, n, (void **)&base, &cnt);
    }
    laik_switchto(data, overlapingPartitioning, LAIK_DF_CopyIn);

    this -> count = cnt;

}

double& laik_vector_overlapping::operator [](int idx){
    uint64_t cnt;
    double* base;

    int i = (idx % (count*count) ) % count;
    int l_idx = i;
    int slice = idx/count;

    //laik_log(Laik_LogLevel(2),"%d",idx);

    laik_map_def(data, slice, (void **)&base, &cnt);

    return base[l_idx];
}

void laik_vector_overlapping::switch_to_write_phase(){
    laik_switchto(data, overlapingPartitioning, Laik_DataFlow ( LAIK_DF_ReduceOut | LAIK_DF_Sum ) );
}

void laik_vector_overlapping::switch_to_reduction(){
    laik_switchto(data, overlapingPartitioning, LAIK_DF_CopyIn);
}
