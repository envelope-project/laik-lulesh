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
 Laik_Partitioning *p1, Laik_Partitioning *p2, Laik_Transition* t1, Laik_Transition* t2, Laik_ReductionOperation operation = LAIK_RO_None);
    void resize(int count);
    inline T& operator [](int idx);
    virtual void calculate_pointers() = 0;
    void switch_to_write_phase();
    void switch_to_read_phase();
    virtual void migrate(Laik_Group* new_group, Laik_Partitioning* p_new_1, Laik_Partitioning* p_new_2, Laik_Transition* t_new_1, Laik_Transition* t_new_2) = 0;
    void test_print();
    void init_config_params(Laik_Group* group);

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
    laik_vector_halo(Laik_Instance* inst, Laik_Group* world, Laik_Space* indexSpace, Laik_Partitioning *p1, Laik_Partitioning *p2, Laik_Transition* t1, Laik_Transition* t2, Laik_ReductionOperation operation = LAIK_RO_None);
    inline T& operator [](int idx);
    T* calc_pointer(int idx, int state);
    void calculate_pointers();
    void resize(int count);
    void switch_to_write_phase();
    void switch_to_read_phase();
    void migrate(Laik_Group* new_group, Laik_Partitioning* p_new_1, Laik_Partitioning* p_new_2, Laik_Transition* t_new_1, Laik_Transition* t_new_2);
};

template <typename T>
class laik_vector_overlapping:public laik_vector<T>
{
public:
    laik_vector_overlapping(Laik_Instance* inst, Laik_Group* world, Laik_Space* indexSpace, Laik_Partitioning *p1, Laik_Partitioning *p2, Laik_Transition* t1, Laik_Transition* t2, Laik_ReductionOperation operation = LAIK_RO_Sum);
    inline T& operator [](int idx);
    T* calc_pointer(int idx);
    void calculate_pointers();
    void resize(int count);
    void switch_to_write_phase();
    void switch_to_read_phase();
    void migrate(Laik_Group* new_group, Laik_Partitioning* p_new_1, Laik_Partitioning* p_new_2, Laik_Transition* t_new_1, Laik_Transition* t_new_2);
};

template <typename T>
class laik_vector_ex_repart:public laik_vector<T>
{
public:
    laik_vector_ex_repart(Laik_Instance* inst, Laik_Group* world, Laik_Space* indexSpace, Laik_Partitioning *p1, Laik_Partitioning *p2, Laik_Transition* t1, Laik_Transition* t2, Laik_ReductionOperation operation = LAIK_RO_None);
    inline T& operator [](int idx);
    T* calc_pointer(int idx, int state);
    void calculate_pointers();
    void resize(int count);
    void switch_to_write_phase();
    void switch_to_read_phase();
    void migrate(Laik_Group* new_group, Laik_Partitioning* p_new_1, Laik_Partitioning* p_new_2, Laik_Transition* t_new_1, Laik_Transition* t_new_2);

private:
    std::vector<T> data_vector;
};

template <typename T>
class laik_vector_overlapping_repart:public laik_vector<T>
{
public:
    laik_vector_overlapping_repart(Laik_Instance* inst, Laik_Group* world, Laik_Space* indexSpace, Laik_Partitioning *p1, Laik_Partitioning *p2, Laik_Transition* t1, Laik_Transition* t2, Laik_ReductionOperation operation = LAIK_RO_None);
    inline T& operator [](int idx);
    T* calc_pointer(int idx, int state);
    void calculate_pointers();
    void resize(int count);
    void switch_to_write_phase();
    void switch_to_read_phase();
    void migrate(Laik_Group* new_group, Laik_Partitioning* p_new_1, Laik_Partitioning* p_new_2, Laik_Transition* t_new_1, Laik_Transition* t_new_2);

private:
    std::vector<T> data_vector;
};

template <typename T>
inline
T& laik_vector_halo<T>::operator [](int idx){
    return *(this -> pointer_cache[idx]);
}

template <typename T>
inline
T& laik_vector_ex_repart<T>::operator [](int idx){
    return this -> data_vector[idx];
}

template <typename T>
inline
T& laik_vector_overlapping<T>::operator [](int idx){
    return *(this -> pointer_cache[idx]);
}

template <typename T>
inline
T& laik_vector_overlapping_repart<T>::operator [](int idx){
    return this -> data_vector[idx];
}

#endif // LAIK_VECTOR
