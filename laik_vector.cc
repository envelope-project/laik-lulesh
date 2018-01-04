#include <laik_vector.h>
#include <laik_port.h>
#include <lulesh.h>

//laik_vector::laik_vector(){}

laik_vector::laik_vector(Laik_Instance* inst, Laik_Group* world){
    this -> inst = inst;
    this -> world = world;
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
    this -> switch_to_exclusive_partitioning();
    uint64_t cnt;
    double* base;
    int nSlices = laik_my_slicecount(exclusivePartitioning);
    for (int n = 0; n < nSlices; ++n)
    {
       laik_map_def(data, n, (void **)&base, &cnt);
    }
    this -> switch_to_halo_partitioning();

    this -> count = cnt;

}

double& laik_vector::operator [](int idx){
    uint64_t cnt;
    double* base;

    int i = (idx % (count*count) ) % count;
    int l_idx = i;
    int slice = idx/count;

    laik_map_def(data, slice, (void **)&base, &cnt);

    return base[l_idx];
}

void laik_vector::switch_to_exclusive_partitioning(){
    laik_switchto(data, exclusivePartitioning, LAIK_DF_CopyOut);
}

void laik_vector::switch_to_halo_partitioning(){
    laik_switchto(data, haloPartitioning, LAIK_DF_CopyIn);
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

void laik_vector_overlapping::switch_to_write_phase(){
    laik_switchto(data, overlapingPartitioning, Laik_DataFlow ( LAIK_DF_ReduceOut | LAIK_DF_Sum ) );
}

void laik_vector_overlapping::switch_to_reduction(){
    laik_switchto(data, overlapingPartitioning, LAIK_DF_CopyIn);
}
