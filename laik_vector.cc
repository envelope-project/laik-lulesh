#include <laik_vector.h>
#include "laik_vector_comm_exclusive_halo.h"
#include "laik_vector_comm_overlapping_overlapping.h"
#include "laik_vector_repart_exclusive.h"
#include "laik_vector_repart_overlapping.h"

#include <laik_partitioners.h>
#include <lulesh.h>
#include <limits.h>
#include <type_traits>
#include <string.h>

template <typename T>
laik_vector<T>::laik_vector(Laik_Instance* inst, Laik_Group* world, Laik_Space* indexSpace, Laik_Partitioning *p1, Laik_Partitioning *p2, Laik_Transition* t1, Laik_Transition* t2, Laik_ReductionOperation operation):reduction_operation(operation){
    this -> inst = inst;
    this -> world = world;
    this-> indexSpace = indexSpace;
    this -> p1 = p1;
    this -> p2 = p2;
    this -> t1 = t1;
    this -> t2 = t2;
    this -> pointer_cache = nullptr;

    this -> init_config_params(world);
}

template <typename T>
void laik_vector<T>::init_config_params(Laik_Group* group){
    int col, row, plane, side;
    InitMeshDecomp(laik_size(group), laik_myid(group), &col, &row, &plane, &side);

    b=1;
    f=1;
    d=1;
    u=1;
    l=1;
    r=1;

    if (col==0) {
        l=0;
    }

    if (col==side-1) {
        r=0;
    }

    if (row==0) {
        d=0;
    }

    if (row==side-1) {
        u=0;
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

template <typename T>
void laik_vector<T>::test_print(){
    T *base;
    uint64_t count;
    int nSlices = laik_my_slicecount(p1);
    for (int s = 0; s < nSlices; s++)
    {
        laik_map_def(data, s, (void**) &base, &count);
        for (uint64_t i = 0; i < count; i++)
        {
            laik_log(Laik_LogLevel(2),"%f\n", base[i]);
        }
        laik_log(Laik_LogLevel(2),"\n");
    }
}



template class laik_vector<double>;

