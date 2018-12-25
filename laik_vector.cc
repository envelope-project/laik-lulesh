#include <laik_vector.h>
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
// ////////////////////////////////////////////////////////////////////////
// implementation of laik_vector with halo partitioning (node partitioning)
// ////////////////////////////////////////////////////////////////////////
template <typename T>
void laik_vector_halo<T>::migrate(Laik_Group* new_group, Laik_Partitioning* p_new_1, Laik_Partitioning* p_new_2, Laik_Transition* t_new_1, Laik_Transition* t_new_2){
    uint64_t cnt;
    int* base;
    //int slice = 0;

    this->init_config_params(new_group);

    laik_switchto_partitioning(this->data, this->p1, LAIK_DF_None, LAIK_RO_None);

    // use the reservation API to precalculate the pointers
    Laik_Reservation* reservation = laik_reservation_new(this->data);
    laik_reservation_add(reservation, p_new_2);
    laik_reservation_add(reservation, p_new_1);
    laik_reservation_alloc(reservation);
    laik_data_use_reservation(this->data, reservation);

    laik_switchto_partitioning(this->data, p_new_1, LAIK_DF_Preserve, LAIK_RO_None);

    if (laik_myid(new_group)< 0) {
        return;
    }

    this->as1 = laik_calc_actions(this->data, t_new_1, reservation, reservation);
    this->as2 = laik_calc_actions(this->data, t_new_2, reservation, reservation);

    this -> p1=p_new_1;
    this -> p2=p_new_2;
    this -> t1=t_new_1;
    this -> t2=t_new_2;
    this -> world = new_group;
    if (laik_myid(this->world)<0)
        return ;

    laik_switchto_partitioning(this->data, this->p1, LAIK_DF_None, LAIK_RO_None);
    int nSlices = laik_my_slicecount(this->p1);
    for (int n = 0; n < nSlices; n++)
    {
        laik_map_def(this->data, n, (void **)&base, &cnt);
    }
    laik_switchto_partitioning(this->data, this->p2, LAIK_DF_Preserve, LAIK_RO_None);

    this -> state = 0;
    this -> count = cnt;

    this -> calculate_pointers();
}


template <typename T>
laik_vector_halo<T>::laik_vector_halo(Laik_Instance *inst,
                                   Laik_Group *world,
                                      Laik_Space* indexSpace,
                                      Laik_Partitioning *p1, Laik_Partitioning *p2, Laik_Transition* t1, Laik_Transition* t2, Laik_ReductionOperation operation):laik_vector<T>(inst,world, indexSpace, p1, p2, t1, t2, operation){}
template <typename T>
void laik_vector_halo<T>::resize(int count){

    this->size = count;

    if (std::is_same <T, double>::value) {
        this->data = laik_new_data(this->indexSpace, laik_Double );

    }
    else if (std::is_same <T, int>::value){
        this->data = laik_new_data(this->indexSpace, laik_Int64 );
    }

    // use the reservation API to precalculate the pointers
    Laik_Reservation* reservation = laik_reservation_new(this->data);
    laik_reservation_add(reservation, this->p2);
    laik_reservation_add(reservation, this->p1);

    laik_reservation_alloc(reservation);
    laik_data_use_reservation(this->data, reservation);

    this->as1 = laik_calc_actions(this->data, this->t1, reservation, reservation);
    this->as2 = laik_calc_actions(this->data, this->t2, reservation, reservation);

    // go through the slices to just allocate the memory
    uint64_t cnt;
    T* base;

    //laik_exec_transition(data, toExclusive);
    laik_switchto_partitioning(this->data, this->p1, LAIK_DF_None, this->reduction_operation);
    int nSlices = laik_my_slicecount(this->p1);
    for (int n = 0; n < nSlices; ++n)
    {
       laik_map_def(this->data, n, (void **)&base, &cnt);
    }
    //laik_exec_transition(data, toHalo);
    laik_switchto_partitioning(this->data, this->p2, LAIK_DF_Preserve, this->reduction_operation);

    this -> count = cnt;
    this -> calculate_pointers();
}

template <typename T>
T* laik_vector_halo<T>::calc_pointer(int idx, int state){
    uint64_t cnt;
    T* base;

    int i=0;
    int j=0;
    int k=0;
    int index=0;
    int slice=-1;
    int side =-1;
    int numElem=this->count*this->count*this->count;

    Index_t ghostIdx[6] ;  // offsets to ghost locations

    for (Index_t i=0; i<6; ++i) {
      ghostIdx[i] = INT_MAX;
    }

    Int_t pidx = numElem ;
    if (this->b) {
      ghostIdx[0] = pidx ;
      pidx += this->count*this->count ;
    }
    if (this->f) {
      ghostIdx[1] = pidx ;
      pidx += this->count*this->count ;
    }
    if (this->d) {
      ghostIdx[2] = pidx ;
      pidx += this->count*this->count ;
    }
    if (this->u) {
      ghostIdx[3] = pidx ;
      pidx += this->count*this->count ;
    }
    if (this->l) {
      ghostIdx[4] = pidx ;
      pidx += this->count*this->count ;
    }
    if (this->r) {
      ghostIdx[5] = pidx ;
    }

    for (int i = 0; i < 6; ++i) {
        if  (idx>=ghostIdx[i]){
            side =i;
        }
    }
    // tested OK
    if (state) {
        i = idx % this->count;
        slice = idx/this->count;
    }
    else{
        // tested OK
        if (idx < numElem) {
            i = idx%this->count;
            j = (idx/this->count)%this->count;
            k = idx/(this->count*this->count);
            slice = (this->count+this->d+this->u)*(k+this->b)+(j+this->d);
            i += this->l;
        }

        // back
        // tested OK
        else if (side==0) {
            index =idx-ghostIdx[side];
            i = index%this->count;
            j = index/this->count;
            k = -1;
            slice = (this->count+this->d+this->u)*(k+this->b)+(j+this->d);
            i += this->l;
        }
        // front
        // tested OK
        else if (side==1) {
            index =idx-ghostIdx[side];
            i = index%this->count;
            j = index/this->count;
            k = this->count;
            slice = (this->count+this->d+this->u)*(k+this->b)+(j+this->d);
            i += this->l;
        }
        //down
        //tested OK
        else if (side==2) {
            index =idx-ghostIdx[side];
            i = index%this->count;
            j = -1;
            k = index/this->count;
            slice = (this->count+this->d+this->u)*(k+this->b)+(j+this->d);
            i += this->l;
        }
        //up
        //tested OK
        else if (side==3) {
            index =idx-ghostIdx[side];
            i = index%this->count;
            j = this->count;
            k = index/this->count;
            slice = (this->count+this->d+this->u)*(k+this->b)+(j+this->d);
            i += this->l;
        }
        //left
        //tested OK
        else if (side==4) {
            index =idx-ghostIdx[side];
            i = -1;
            j = index%this->count;
            k = index/this->count;
            slice = (this->count+this->d+this->u)*(k+this->b)+(j+this->d);
            i += this->l;
        }
        //right
        //tested OK
        else if (side==5) {
            index =idx-ghostIdx[side];
            i = this->count;
            j = index%this->count;
            k = index/this->count;
            slice = (this->count+this->d+this->u)*(k+this->b)+(j+this->d);
            i += this->l;
        }
    }

    if (slice>-1) {
        laik_map_def(this->data, slice, (void **)&base, &cnt);
        //laik_log(Laik_LogLevel(2),"state: %d, idx: %d, pointer: %x", state, idx, base+i);
        return base+i;
    }

    return &(this->zero);
}

template <typename T>
void laik_vector_halo<T>::calculate_pointers(){
    int numElems = this->count*this->count*this->count;

    if (this->pointer_cache != nullptr)   free (this->pointer_cache);

    this->pointer_cache= (T**) malloc (numElems * sizeof(T*));
    laik_switchto_partitioning(this->data, this->p1, LAIK_DF_None, this->reduction_operation);

    for (int i = 0; i < numElems; ++i) {
        this->pointer_cache [i] = calc_pointer(i,1);
    }
    // test
    /*
    for (int i = 0; i < numElems; ++i) {
        laik_log((Laik_LogLevel)2,"exclusive pointers: %d %x\n", i, exclusive_pointers [i]);
    }
    */

    int numElemsTotal = numElems + (this->b+this->f+this->d+this->u+this->l+this->r)*this->count*this->count;
    this->pointer_cache= (T**) malloc (numElemsTotal * sizeof(T*));
    laik_switchto_partitioning(this->data, this->p2, LAIK_DF_Preserve, this->reduction_operation);

    for (int i = 0; i < numElemsTotal; ++i) {
        this->pointer_cache [i] = calc_pointer(i,0);
    }
    // test
    /*
    for (int i = 0; i < numElemsTotal; ++i) {
        laik_log((Laik_LogLevel)2,"halo pointers: %d %x\n", i, halo_pointers [i]);
    }
    */

    //laik_log((Laik_LogLevel)2,"debug for the reservation api");
}

template <typename T>
void laik_vector_halo<T>::switch_to_write_phase(){
    laik_exec_actions(this->as1);
    //laik_exec_transition(data,t1);
    //laik_switchto_phase(data, p1, LAIK_DF_CopyOut);
    this->state=1;
}

template <typename T>
void laik_vector_halo<T>::switch_to_read_phase(){
    laik_exec_actions(this->as2);
    //laik_exec_transition(data,t2);
    //laik_switchto_phase(data, p2, LAIK_DF_CopyIn);
    this->state=0;
}


// ////////////////////////////////////////////////////////////////////////
// implementation of laik_vector with exclusive partitioning (elem partitioning)
// for repartitioning of exclusive data structs
// ////////////////////////////////////////////////////////////////////////
template <typename T>
void laik_vector_ex_repart<T>::migrate(Laik_Group* new_group, Laik_Partitioning* p_new_1, Laik_Partitioning* p_new_2, Laik_Transition* t_new_1, Laik_Transition* t_new_2){
    uint64_t cnt;
    T* base;
    int nSlices;

    this->init_config_params(new_group);

    laik_switchto_partitioning(this->data, this->p1, LAIK_DF_Preserve, LAIK_RO_None);
    // copy the data from stl vector into the laik container
    nSlices = laik_my_slicecount(this->p1);
    for (int n = 0; n < nSlices; n++)
    {
        laik_map_def(this->data, n, (void **)&base, &cnt);
        memcpy(base, &data_vector[0] + n*cnt, cnt*sizeof(T));
        //std::copy( base, base + cnt, data_vector.begin() + n*count );
    }

    // perform switches for communication

    laik_switchto_partitioning(this->data, p_new_1, LAIK_DF_Preserve, LAIK_RO_None);

    this -> world = new_group;
    if (laik_myid(this->world)<0)
        return ;

    this -> p1=p_new_1;
    this -> p2=p_new_2;
    this -> t1=t_new_1;
    this -> t2=t_new_2;

    // resize vector
    laik_map_def(this->data, 0, (void **)&base, &cnt);
    int s = cnt*cnt*cnt;
    data_vector.resize(s);

    // copy the data back into the stl vecotrs
    nSlices = laik_my_slicecount(this->p1);
    for (int n = 0; n < nSlices; n++)
    {
        laik_map_def(this->data, n, (void **)&base, &cnt);
        memcpy(&data_vector[0] + n*cnt, base, cnt*sizeof(T));
        //std::copy(data_vector.begin() + n*count ,data_vector.begin() + (n+1)*count-1 , base);
    }
}

template <typename T>
laik_vector_ex_repart<T>::laik_vector_ex_repart(Laik_Instance *inst,
                                   Laik_Group *world,
                                      Laik_Space* indexSpace,
                                      Laik_Partitioning *p1, Laik_Partitioning *p2, Laik_Transition* t1, Laik_Transition* t2, Laik_ReductionOperation operation):laik_vector<T>(inst,world, indexSpace, p1, p2, t1, t2, operation){}
template <typename T>
void laik_vector_ex_repart<T>::resize(int count){

    int s =  count / laik_size(this->world);
    data_vector.resize(s);

    this -> size = count;

    if (std::is_same <T, double>::value) {
        this->data = laik_new_data(this->indexSpace, laik_Double );

    }
    else if (std::is_same <T, int>::value){
        this->data = laik_new_data(this->indexSpace, laik_Int64 );
    }

    laik_switchto_partitioning(this->data, this->p1, LAIK_DF_None, this->reduction_operation);
    uint64_t cnt;
    T* base;
    int nSlices = laik_my_slicecount(this->p1);
    for (int n = 0; n < nSlices; ++n)
    {
       laik_map_def(this->data, n, (void **)&base, &cnt);
    }
    this -> count = cnt;
}

template <typename T>
T* laik_vector_ex_repart<T>::calc_pointer(int idx, int state){
     return &(this->zero);
}

template <typename T>
void laik_vector_ex_repart<T>::calculate_pointers(){

}

template <typename T>
void laik_vector_ex_repart<T>::switch_to_write_phase(){

}

template <typename T>
void laik_vector_ex_repart<T>::switch_to_read_phase(){

}


// ///////////////////////////////////////////////////////////
// implementation of laik_vector with overlapping partitioning
// (node partitioning)
// ///////////////////////////////////////////////////////////

template <typename T>
laik_vector_overlapping<T>::laik_vector_overlapping(Laik_Instance *inst,
                                                 Laik_Group*world,
                                                    Laik_Space* indexSpace,
                                                   Laik_Partitioning *p1, Laik_Partitioning *p2,
                                                    Laik_Transition* t1, Laik_Transition* t2,
                                                    Laik_ReductionOperation operation):laik_vector<T>(inst,world, indexSpace, p1, p2, t1,t2, operation){}
template <typename T>
void laik_vector_overlapping<T>::resize(int count){

    this->size = count;

    if (std::is_same <T, double>::value) {
        this->data = laik_new_data(this->indexSpace, laik_Double );

    }
    else if (std::is_same <T, int>::value){
        this->data = laik_new_data(this->indexSpace, laik_Int64 );
    }
    // use the reservation API to precalculate the pointers

    Laik_Reservation* reservation = laik_reservation_new(this->data);
    laik_reservation_add(reservation, this->p1);

    laik_reservation_alloc(reservation);
    laik_data_use_reservation(this->data, reservation);

    // precalculate the transition object

    this->as1 = laik_calc_actions(this->data, this->t1, reservation, reservation);
    this->as2 = laik_calc_actions(this->data, this->t2, reservation, reservation);

    // go through the slices to just allocate the memory
    laik_switchto_partitioning(this->data, this->p1, LAIK_DF_None, this->reduction_operation );
    Laik_TaskSlice* ts = laik_my_slice(this->p1, 0);
    const Laik_Slice* s = laik_taskslice_get_slice(ts);
    this -> count = laik_slice_size(s);
    this -> calculate_pointers();
}

/*
double& laik_vector_overlapping::operator [](int idx){
    return *(this -> overlapping_pointers[idx]);
}
*/

template <typename T>
T* laik_vector_overlapping<T>::calc_pointer(int idx){
    uint64_t cnt;
    T* base;

    int i = (idx % (this->count*this->count) ) % this->count;
    int l_idx = i;
    int slice = idx/this->count;
    laik_map_def(this->data, slice, (void **)&base, &cnt);

    return base+l_idx;
}

template <typename T>
void laik_vector_overlapping<T>::calculate_pointers(){

    if (this->pointer_cache != nullptr)   free (this->pointer_cache);

    this->pointer_cache = (T**) malloc (this->size * sizeof(T*));

    for (int i = 0; i < this->size; ++i) {
        this->pointer_cache [i] = calc_pointer(i);
    }
}

template <typename T>
void laik_vector_overlapping<T>::switch_to_write_phase(){
    laik_exec_actions(this->as1);
    //laik_exec_transition(data,t1);
    //laik_switchto_phase(data, p1, Laik_DataFlow
    //                    ( LAIK_DF_ReduceOut | LAIK_DF_Sum ) );
}

template <typename T>
void laik_vector_overlapping<T>::switch_to_read_phase(){
    laik_exec_actions(this->as2);
    //laik_exec_transition(data,t2);
    //laik_switchto_phase(data, p1, LAIK_DF_CopyIn);
}

template <typename T>
void laik_vector_overlapping<T>::migrate(Laik_Group* new_group, Laik_Partitioning* p_new_1, Laik_Partitioning* p_new_2, Laik_Transition* t_new_1, Laik_Transition* t_new_2){
    uint64_t cnt;
    int* base;
    //int slice = 0;

    this->init_config_params(new_group);

    laik_switchto_partitioning(this->data, this->p1, LAIK_DF_None, LAIK_RO_Min);

    Laik_Reservation* reservation = laik_reservation_new(this->data);
    laik_reservation_add(reservation, p_new_1);
    laik_reservation_alloc(reservation);
    laik_data_use_reservation(this->data, reservation);

    laik_switchto_partitioning(this->data, p_new_1, LAIK_DF_Preserve, LAIK_RO_Min);

    if (laik_myid(new_group)< 0) {
        return;
    }

    this->as1 = laik_calc_actions(this->data, t_new_1, reservation, reservation);
    this->as2 = laik_calc_actions(this->data, t_new_2, reservation, reservation);

    this -> p1=p_new_1;
    this -> p2=p_new_2;
    this -> t1=t_new_1;
    this -> t2=t_new_2;
    this -> world = new_group;
    if (laik_myid(this->world)<0)
        return ;

    laik_switchto_partitioning(this->data, this->p1, LAIK_DF_None, LAIK_RO_Min );
    int nSlices = laik_my_slicecount(this->p1);
    for (int n = 0; n < nSlices; n++)
    {
        laik_map_def(this->data, n, (void **)&base, &cnt);
    }
    laik_switchto_partitioning(this->data, this->p2, LAIK_DF_Preserve, LAIK_RO_Min);

    this -> count = cnt;
    this -> calculate_pointers();
}

// ////////////////////////////////////////////////////////////////////////
// implementation of laik_vector with overlapping partitioning (node partitioning)
// for repartitioning of overlapping data structs
// ////////////////////////////////////////////////////////////////////////
template <typename T>
void laik_vector_overlapping_repart<T>::migrate(Laik_Group* new_group, Laik_Partitioning* p_new_1, Laik_Partitioning* p_new_2, Laik_Transition* t_new_1, Laik_Transition* t_new_2){
    uint64_t cnt;
    T* base;
    int nSlices;

    this->init_config_params(new_group);

    laik_switchto_partitioning(this->data, this->p1, LAIK_DF_None, LAIK_RO_Min);
    // copy the data from stl vector into the laik container
    nSlices = laik_my_slicecount(this->p1);
    for (int n = 0; n < nSlices; n++)
    {
        laik_map_def(this->data, n, (void **)&base, &cnt);
        memcpy(base, &data_vector[0] + n*cnt, cnt*sizeof(T));
        //std::copy( base, base + cnt, data_vector.begin() + n*count );
    }

    // perform switches for communication
    laik_switchto_partitioning(this->data, p_new_1, LAIK_DF_Preserve, LAIK_RO_Min);

    this -> world = new_group;
    if (laik_myid(this->world)<0)
        return ;

    this -> p1=p_new_1;
    this -> p2=p_new_2;
    this -> t1=t_new_1;
    this -> t2=t_new_2;

    // resize vector
    laik_map_def(this->data, 0, (void **)&base, &cnt);
    int s = cnt*cnt*cnt;
    data_vector.resize(s);

    // copy the data back into the stl vecotrs
    nSlices = laik_my_slicecount(this->p1);
    for (int n = 0; n < nSlices; n++)
    {
        laik_map_def(this->data, n, (void **)&base, &cnt);
        memcpy(&data_vector[0] + n*cnt, base, cnt*sizeof(T));
        //std::copy(data_vector.begin() + n*count ,data_vector.begin() + (n+1)*count-1 , base);
    }

}

template <typename T>
laik_vector_overlapping_repart<T>::laik_vector_overlapping_repart(Laik_Instance *inst,
                                   Laik_Group *world,
                                      Laik_Space* indexSpace,
                                      Laik_Partitioning *p1, Laik_Partitioning *p2, Laik_Transition* t1, Laik_Transition* t2, Laik_ReductionOperation operation):laik_vector<T>(inst,world, indexSpace, p1, p2, t1, t2, operation){}
template <typename T>
void laik_vector_overlapping_repart<T>::resize(int count){

    int side = cbrt (laik_size(this->world));
    int s = (int) ((cbrt (count)  -  1 ) / side + 1 + 0.1 );
    s = s*s*s;
    data_vector.resize(s);

    this -> size = count;

    if (std::is_same <T, double>::value) {
        this->data = laik_new_data(this->indexSpace, laik_Double );

    }
    else if (std::is_same <T, int>::value){
        this->data = laik_new_data(this->indexSpace, laik_Int64 );
    }

    laik_switchto_partitioning(this->data, this->p1, LAIK_DF_None, LAIK_RO_Min);
    Laik_TaskSlice* ts = laik_my_slice(this->p1, 0);
    const Laik_Slice* sl = laik_taskslice_get_slice(ts);
    this -> count = laik_slice_size(sl);
}

template <typename T>
T* laik_vector_overlapping_repart<T>::calc_pointer(int idx, int state){
     return &(this->zero);
}

template <typename T>
void laik_vector_overlapping_repart<T>::calculate_pointers(){

}

template <typename T>
void laik_vector_overlapping_repart<T>::switch_to_write_phase(){

}

template <typename T>
void laik_vector_overlapping_repart<T>::switch_to_read_phase(){

}


template class laik_vector<double>;
template class laik_vector_halo<int>;
template class laik_vector_halo<double>;
//template class laik_vector_ex_repart<int>;
template class laik_vector_ex_repart<double>;
template class laik_vector_overlapping<double>;
template class laik_vector_overlapping_repart<double>;

