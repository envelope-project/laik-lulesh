#ifndef LAIK_VECTOR
#define LAIK_VECTOR

extern "C"{
#include "laik.h"
#include "laik-backend-mpi.h"
}
#include<vector>

class laik_vector:public std::vector<double>
{
public:
    //laik_vector();
    laik_vector(Laik_Instance* inst, Laik_Group* world);
    void resize(int count);
    // TODO accessors
    double& operator [](int idx);
    void switch_to_exclusive_partitioning();
    void switch_to_halo_partitioning();
    void test_print();

protected:
    Laik_Instance* inst;
    Laik_Group* world;
    int size;
    Laik_Space* indexSpace;
    Laik_Partitioning *exclusivePartitioning;
    Laik_Partitioning *haloPartitioning;
    Laik_Partitioning *overlapingPartitioning;
    Laik_Data* data;
    int count;
    int f,b,u,d,l,r;
    int state;
    double zero;

};

class laik_vector_halo:public laik_vector
{
public:
    laik_vector_halo(Laik_Instance* inst, Laik_Group* world);
};

class laik_vector_overlapping:public laik_vector
{
public:
    laik_vector_overlapping(Laik_Instance* inst, Laik_Group* world);
    double& operator [](int idx);
    void resize(int count);
    void switch_to_write_phase();
    void switch_to_reduction();
};

#endif // LAIK_VECTOR
