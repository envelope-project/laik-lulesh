#ifndef LAIK_VECTOR
#define LAIK_VECTOR

extern "C"{
#include "laik.h"
#include "laik-backend-mpi.h"
}
#include<vector>

class laik_vector
{
public:
    laik_vector(Laik_Instance* inst, Laik_Group* world);
    void resize(int count);
    double& operator [](int idx);
    void calculate_pointers();
    void switch_to_write_phase();
    void switch_to_read_phase();
    void test_print();

protected:
    Laik_Instance* inst;
    Laik_Group* world;
    int size;
    Laik_Space* indexSpace;
    Laik_AccessPhase *exclusivePartitioning;
    Laik_AccessPhase *haloPartitioning;
    Laik_AccessPhase *overlapingPartitioning;
    Laik_Transition *toW;
    Laik_Transition *toR;
    Laik_Data* data;
    int count;
    int f,b,u,d,l,r;
    int state;
    double zero;
    double **exclusive_pointers;
    double **halo_pointers;
    double **overlapping_pointers;
};

class laik_vector_halo:public laik_vector
{
public:
    laik_vector_halo(Laik_Instance* inst, Laik_Group* world);
    double& operator [](int idx);
    double* calc_pointer(int idx, int state);
    void calculate_pointers();
    void resize(int count);
    void switch_to_write_phase();
    void switch_to_read_phase();
};

class laik_vector_overlapping:public laik_vector
{
public:
    laik_vector_overlapping(Laik_Instance* inst, Laik_Group* world);
    double& operator [](int idx);
    double* calc_pointer(int idx);
    void calculate_pointers();
    void resize(int count);
    void switch_to_write_phase();
    void switch_to_read_phase();
};

#endif // LAIK_VECTOR
