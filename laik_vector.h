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
    
    // virtual member functions to be implemented in conrete laik_vectors
    virtual void resize(int count) = 0;
    inline virtual T& operator [](int idx) = 0;
    virtual void calculate_pointers() = 0;
    virtual void switch_to_write_phase() = 0;
    virtual void switch_to_read_phase() = 0;
    virtual void migrate(Laik_Group* new_group, Laik_Partitioning* p_new_1, Laik_Partitioning* p_new_2, Laik_Transition* t_new_1, Laik_Transition* t_new_2) = 0;
    void test_print();
    void init_config_params(Laik_Group* group);

protected:
    // members from laik
    Laik_Instance* inst;
    Laik_Group* world;
    Laik_ReductionOperation reduction_operation;
    Laik_Space* indexSpace;
    Laik_Partitioning *p1;
    Laik_Partitioning *p2;
    Laik_Transition *toW;
    Laik_Transition *toR;
    Laik_ActionSeq* asW;
    Laik_ActionSeq* asR;
    Laik_Data* data;

    T **pointer_cache;
    
    int size;

    int count;
    int f,b,u,d,l,r;
    int state;
    T zero;
};

template <typename T>
class laik_vector_halo:public laik_vector<T>
{
public:
    laik_vector_halo(Laik_Instance* inst, Laik_Group* world, Laik_Space* indexSpace, Laik_Partitioning *p1, Laik_Partitioning *p2, Laik_Transition* t1, Laik_Transition* t2, Laik_ReductionOperation operation = LAIK_RO_None);
    inline T& operator [](int idx) override;
    T* calc_pointer(int idx, int state);
    void calculate_pointers() override;
    void resize(int count) override;
    void switch_to_write_phase() override;
    void switch_to_read_phase() override;
    void migrate(Laik_Group* new_group, Laik_Partitioning* p_new_1, Laik_Partitioning* p_new_2, Laik_Transition* t_new_1, Laik_Transition* t_new_2) override;
};

template <typename T>
class laik_vector_overlapping:public laik_vector<T>
{
public:
    laik_vector_overlapping(Laik_Instance* inst, Laik_Group* world, Laik_Space* indexSpace, Laik_Partitioning *p1, Laik_Partitioning *p2, Laik_Transition* t1, Laik_Transition* t2, Laik_ReductionOperation operation = LAIK_RO_Sum);
    inline T& operator [](int idx) override;
    T* calc_pointer(int idx);
    void calculate_pointers() override;
    void resize(int count) override;
    void switch_to_write_phase() override;
    void switch_to_read_phase() override;
    void migrate(Laik_Group* new_group, Laik_Partitioning* p_new_1, Laik_Partitioning* p_new_2, Laik_Transition* t_new_1, Laik_Transition* t_new_2) override;
};

template <typename T>
class laik_vector_ex_repart:public laik_vector<T>
{
public:
    laik_vector_ex_repart(Laik_Instance* inst, Laik_Group* world, Laik_Space* indexSpace, Laik_Partitioning *p1, Laik_Partitioning *p2, Laik_Transition* t1, Laik_Transition* t2, Laik_ReductionOperation operation = LAIK_RO_None);
    inline T& operator [](int idx) override;
    T* calc_pointer(int idx, int state);
    void calculate_pointers() override;
    void resize(int count) override;
    void switch_to_write_phase() override;
    void switch_to_read_phase() override;
    void migrate(Laik_Group* new_group, Laik_Partitioning* p_new_1, Laik_Partitioning* p_new_2, Laik_Transition* t_new_1, Laik_Transition* t_new_2) override;

private:
    std::vector<T> data_vector;
};

template <typename T>
class laik_vector_overlapping_repart:public laik_vector<T>
{
public:
    laik_vector_overlapping_repart(Laik_Instance* inst, Laik_Group* world, Laik_Space* indexSpace, Laik_Partitioning *p1, Laik_Partitioning *p2, Laik_Transition* t1, Laik_Transition* t2, Laik_ReductionOperation operation = LAIK_RO_None);
    inline T& operator [](int idx) override;
    T* calc_pointer(int idx, int state);
    void calculate_pointers() override;
    void resize(int count) override;
    void switch_to_write_phase() override;
    void switch_to_read_phase() override;
    void migrate(Laik_Group* new_group, Laik_Partitioning* p_new_1, Laik_Partitioning* p_new_2, Laik_Transition* t_new_1, Laik_Transition* t_new_2) override;

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
