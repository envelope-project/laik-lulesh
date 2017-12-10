#include <laik_vector.h>
#include <laik_port.h>

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
    data = laik_new_data(world, indexSpace, laik_Double);

    this -> switch_to_exclusive_partitioning();
    // go through the slices to just allocate the memory
    uint64_t cnt;
    double* base;
    int nSlices = laik_my_slicecount(exclusivePartitioning);
    for (int n = 0; n < nSlices; ++n)
    {
       laik_map_def(data, n, (void **)&base, &cnt);
    }
    this->switch_to_halo_partitioning();
}

double& laik_vector::operator [](int idx){
    uint64_t cnt;
    double* base;
    int slice = 0;
    int l_idx = 0;
    laik_map_def(data, slice, (void **)&base, &cnt);
    slice = idx / cnt;
    laik_map_def(data, slice, (void **)&base, &cnt);
    l_idx = idx % cnt;

    return base[l_idx];
}

void laik_vector::switch_to_exclusive_partitioning(){
    laik_switchto(data, exclusivePartitioning, LAIK_DF_CopyOut);
}

void laik_vector::switch_to_halo_partitioning(){
    laik_switchto(data, haloPartitioning, LAIK_DF_CopyIn);
}

