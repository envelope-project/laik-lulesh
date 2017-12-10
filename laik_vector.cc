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
}

