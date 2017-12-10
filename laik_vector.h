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

private:
    Laik_Instance* inst;
    Laik_Group* world;
    size_t size;
    Laik_Space* indexSpace;
    Laik_Partitioning *exclusivePartitioning;
    Laik_Partitioning *haloPartitioning;
    Laik_Data* data;
};

#endif // LAIK_VECTOR
