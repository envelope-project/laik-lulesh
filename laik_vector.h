#ifndef LAIK_VECTOR
#define LAIK_VECTOR

extern "C"{
#include "laik.h"
#include "laik-backend-mpi.h"
}
#include<vector>

template <typename T>
class laik_vector
{
public:
    laik_vector(Laik_Instance* inst, Laik_Group* world, Laik_Space* indexSpace,
 Laik_Partitioning *p1, Laik_Partitioning *p2, Laik_ReductionOperation operation = LAIK_RO_None);
    void resize(int count);
    inline T& operator [](int idx);
    void calculate_pointers();
    void switch_to_write_phase();
    void switch_to_read_phase();
    void test_print();

protected:
    Laik_Instance* inst;
    Laik_Group* world;
    int size;
    Laik_Space* indexSpace;
    Laik_Partitioning *p1;
    Laik_Partitioning *p2;
    Laik_Transition *toW;
    Laik_Transition *toR;
    Laik_ActionSeq* asW;
    Laik_ActionSeq* asR;
    Laik_Data* data;
    int count;
    int f,b,u,d,l,r;
    int state;
    T zero;
    T **pointer_cache;
    Laik_ReductionOperation reduction_operation;
};

template <typename T>
class laik_vector_halo:public laik_vector<T>
{
public:
    laik_vector_halo(Laik_Instance* inst, Laik_Group* world, Laik_Space* indexSpace, Laik_Partitioning *p1, Laik_Partitioning *p2, Laik_ReductionOperation operation = LAIK_RO_None);
    inline T& operator [](int idx);
    T* calc_pointer(int idx, int state);
    void calculate_pointers();
    void resize(int count);
    void switch_to_write_phase();
    void switch_to_read_phase();

protected:
    using laik_vector<T>::inst;
    using laik_vector<T>::world;
    using laik_vector<T>::size;
    using laik_vector<T>::indexSpace;
    using laik_vector<T>::p1;
    using laik_vector<T>::p2;
    using laik_vector<T>::toW;
    using laik_vector<T>::toR;
    using laik_vector<T>::asW;
    using laik_vector<T>::asR;
    using laik_vector<T>::data;
    using laik_vector<T>::count;
    using laik_vector<T>::f;
    using laik_vector<T>::b;
    using laik_vector<T>::u;
    using laik_vector<T>::d;
    using laik_vector<T>::l;
    using laik_vector<T>::r;
    using laik_vector<T>::state;
    using laik_vector<T>::zero;
    using laik_vector<T>::pointer_cache;
    using laik_vector<T>::reduction_operation;
};

template <typename T>
class laik_vector_overlapping:public laik_vector<T>
{
public:
    laik_vector_overlapping(Laik_Instance* inst, Laik_Group* world, Laik_Space* indexSpace, Laik_Partitioning *p1, Laik_Partitioning *p2, Laik_ReductionOperation operation = LAIK_RO_Sum);
    inline T& operator [](int idx);
    T* calc_pointer(int idx);
    void calculate_pointers();
    void resize(int count);
    void switch_to_write_phase();
    void switch_to_read_phase();

protected:
    using laik_vector<T>::inst;
    using laik_vector<T>::world;
    using laik_vector<T>::size;
    using laik_vector<T>::indexSpace;
    using laik_vector<T>::p1;
    using laik_vector<T>::p2;
    using laik_vector<T>::toW;
    using laik_vector<T>::toR;
    using laik_vector<T>::asW;
    using laik_vector<T>::asR;
    using laik_vector<T>::data;
    using laik_vector<T>::count;
    using laik_vector<T>::f;
    using laik_vector<T>::b;
    using laik_vector<T>::u;
    using laik_vector<T>::d;
    using laik_vector<T>::l;
    using laik_vector<T>::r;
    using laik_vector<T>::state;
    using laik_vector<T>::zero;
    using laik_vector<T>::pointer_cache;
    using laik_vector<T>::reduction_operation;
};


template <typename T>
inline
T& laik_vector_halo<T>::operator [](int idx){
    return *(this -> pointer_cache[idx]);
}

template <typename T>
inline
T& laik_vector_overlapping<T>::operator [](int idx){
    return *(this -> pointer_cache[idx]);
}

#endif // LAIK_VECTOR
