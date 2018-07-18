#include <laik_vector.h>
#include <laik_partitioners.h>
#include <lulesh.h>
#include<limits.h>
#include <type_traits>

template <typename T>
laik_vector<T>::laik_vector(Laik_Instance* inst, Laik_Group* world, Laik_Space* indexSpace, Laik_Partitioning *p1, Laik_Partitioning *p2, Laik_Transition* t1, Laik_Transition* t2, Laik_ReductionOperation operation):reduction_operation(operation){
    this -> inst = inst;
    this -> world = world;
    this-> indexSpace = indexSpace;
    this -> p1 = p1;
    this -> p2 = p2;
    this -> toW = t1;
    this -> toR = t2;

    int numRanks = laik_size(world);
    int myRank = laik_myid(world);
    int col, row, plane, side;
    InitMeshDecomp(numRanks, myRank, &col, &row, &plane, &side);

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
    for (size_t s = 0; s < nSlices; s++)
    {
        laik_map_def(data, s, (void**) &base, &count);
        for (uint64_t i = 0; i < count; i++)
        {
            laik_log(Laik_LogLevel(2),"%f\n", base[i]);
        }
        laik_log(Laik_LogLevel(2),"\n");
    }
}

template <typename T>
void laik_vector_halo<T>::migrate(Laik_Group* new_group, Laik_Partitioning* p_new_1, Laik_Partitioning* p_new_2, Laik_Transition* t_new_1, Laik_Transition* t_new_2){
    uint64_t cnt;
    int* base;
    int slice = 0;



    laik_switchto_partitioning(data, p1, LAIK_DF_None, LAIK_RO_None);

    // use the reservation API to precalculate the pointers
    Laik_Reservation* reservation = laik_reservation_new(data);
    laik_reservation_add(reservation, p_new_2);
    laik_reservation_add(reservation, p_new_1);
    laik_reservation_alloc(reservation);
    laik_data_use_reservation(data, reservation);

    laik_switchto_partitioning(data, p_new_1, LAIK_DF_Preserve, LAIK_RO_None);

    if (laik_myid(new_group)< 0) {
        return;
    }

    asW = laik_calc_actions(data, t_new_1, reservation, reservation);
    asR = laik_calc_actions(data, t_new_2, reservation, reservation);

    this -> state = 1;

    this -> p1=p_new_1;
    this -> p2=p_new_2;
    this -> toW=t_new_1;
    this -> toR=t_new_2;
    this -> world = new_group;
    if (laik_myid(world)<0)
        return ;

    laik_switchto_partitioning(data, p1, LAIK_DF_None, LAIK_RO_None);
    int nSlices = laik_my_slicecount(p1);
    for (int n = 0; n < nSlices; n++)
    {
        laik_map_def(data, n, (void **)&base, &cnt);
    }
    laik_switchto_partitioning(data, p2, LAIK_DF_Preserve, LAIK_RO_None);

    this -> state = 0;
    this -> count = cnt;

    this -> calculate_pointers();
}

// ////////////////////////////////////////////////////////////////////////
// implementation of laik_vector with halo partitioning (node partitioning)
// ////////////////////////////////////////////////////////////////////////
template <typename T>
laik_vector_halo<T>::laik_vector_halo(Laik_Instance *inst,
                                   Laik_Group *world,
                                      Laik_Space* indexSpace,
                                      Laik_Partitioning *p1, Laik_Partitioning *p2, Laik_Transition* t1, Laik_Transition* t2, Laik_ReductionOperation operation):laik_vector<T>(inst,world, indexSpace, p1, p2, t1, t2, operation){}
template <typename T>
void laik_vector_halo<T>::resize(int count){

    size = count;

    if (std::is_same <T, double>::value) {
        data = laik_new_data(indexSpace, laik_Double );

    }
    else if (std::is_same <T, int>::value){
        data = laik_new_data(indexSpace, laik_Int64 );
    }

    // use the reservation API to precalculate the pointers
    Laik_Reservation* reservation = laik_reservation_new(data);
    laik_reservation_add(reservation, p2);
    laik_reservation_add(reservation, p1);

    laik_reservation_alloc(reservation);
    laik_data_use_reservation(data, reservation);

    asW = laik_calc_actions(data, toW, reservation, reservation);
    asR = laik_calc_actions(data, toR, reservation, reservation);

    // go through the slices to just allocate the memory
    uint64_t cnt;
    T* base;

    //laik_exec_transition(data, toExclusive);
    laik_switchto_partitioning(data, p1, LAIK_DF_None, reduction_operation);
    int nSlices = laik_my_slicecount(p1);
    for (int n = 0; n < nSlices; ++n)
    {
       laik_map_def(data, n, (void **)&base, &cnt);
    }
    //laik_exec_transition(data, toHalo);
    laik_switchto_partitioning(data, p2, LAIK_DF_Preserve, reduction_operation);

    this -> count = cnt;
    this -> calculate_pointers();
}

/*
double& laik_vector_halo::operator [](int idx){
    return *(this -> halo_pointers[idx]);
}
*/
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
    int numElem=count*count*count;

    Index_t ghostIdx[6] ;  // offsets to ghost locations

    for (Index_t i=0; i<6; ++i) {
      ghostIdx[i] = INT_MAX;
    }

    Int_t pidx = numElem ;
    if (b) {
      ghostIdx[0] = pidx ;
      pidx += count*count ;
    }
    if (f) {
      ghostIdx[1] = pidx ;
      pidx += count*count ;
    }
    if (d) {
      ghostIdx[2] = pidx ;
      pidx += count*count ;
    }
    if (u) {
      ghostIdx[3] = pidx ;
      pidx += count*count ;
    }
    if (l) {
      ghostIdx[4] = pidx ;
      pidx += count*count ;
    }
    if (r) {
      ghostIdx[5] = pidx ;
    }

    /*
    for (int i = 0; i < 6; ++i) {
            laik_log(Laik_LogLevel(2), "%d", ghostIdx[i]);
    }
    */

    for (int i = 0; i < 6; ++i) {
        if  (idx>=ghostIdx[i]){
            side =i;
        }
    }
    // tested OK
    if (state) {
        i = idx % count;
        slice = idx/count;
    }
    else{
        // tested OK
        if (idx < numElem) {
            i = idx%count;
            j = (idx/count)%count;
            k = idx/(count*count);
            slice = (count+d+u)*(k+b)+(j+d);
            i += l;
        }

        // back
        // tested OK
        else if (side==0) {
            index =idx-ghostIdx[side];
            i = index%count;
            j = index/count;
            k = -1;
            slice = (count+d+u)*(k+b)+(j+d);
            i += l;
        }
        // front
        // tested OK
        else if (side==1) {
            index =idx-ghostIdx[side];
            i = index%count;
            j = index/count;
            k = count;
            slice = (count+d+u)*(k+b)+(j+d);
            i += l;
        }
        //down
        //tested OK
        else if (side==2) {
            index =idx-ghostIdx[side];
            i = index%count;
            j = -1;
            k = index/count;
            slice = (count+d+u)*(k+b)+(j+d);
            i += l;
        }
        //up
        //tested OK
        else if (side==3) {
            index =idx-ghostIdx[side];
            i = index%count;
            j = count;
            k = index/count;
            slice = (count+d+u)*(k+b)+(j+d);
            i += l;
        }
        //left
        //tested OK
        else if (side==4) {
            index =idx-ghostIdx[side];
            i = -1;
            j = index%count;
            k = index/count;
            slice = (count+d+u)*(k+b)+(j+d);
            i += l;
        }
        //right
        //tested OK
        else if (side==5) {
            index =idx-ghostIdx[side];
            i = count;
            j = index%count;
            k = index/count;
            slice = (count+d+u)*(k+b)+(j+d);
            i += l;
        }
    }

    if (slice>-1) {
        laik_map_def(data, slice, (void **)&base, &cnt);
        //laik_log(Laik_LogLevel(2),"state: %d, idx: %d, pointer: %x", state, idx, base+i);
        return base+i;
    }

    return &zero;
}

template <typename T>
void laik_vector_halo<T>::calculate_pointers(){
    int numElems = count*count*count;
    pointer_cache= (T**) malloc (numElems * sizeof(T*));
    laik_switchto_partitioning(data, p1, LAIK_DF_None, reduction_operation);

    for (int i = 0; i < numElems; ++i) {
        pointer_cache [i] = calc_pointer(i,1);
    }
    // test
    /*
    for (int i = 0; i < numElems; ++i) {
        laik_log((Laik_LogLevel)2,"exclusive pointers: %d %x\n", i, exclusive_pointers [i]);
    }
    */

    int numElemsTotal = numElems + (b+f+d+u+l+r)*count*count;
    pointer_cache= (T**) malloc (numElemsTotal * sizeof(T*));
    laik_switchto_partitioning(data, p2, LAIK_DF_Preserve, reduction_operation);

    for (int i = 0; i < numElemsTotal; ++i) {
        pointer_cache [i] = calc_pointer(i,0);
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
    laik_exec_actions(asW);    
    //laik_exec_transition(data,toW);
    //laik_switchto_phase(data, p1, LAIK_DF_CopyOut);
    state=1;
}

template <typename T>
void laik_vector_halo<T>::switch_to_read_phase(){
    laik_exec_actions(asR);    
    //laik_exec_transition(data,toR);
    //laik_switchto_phase(data, p2, LAIK_DF_CopyIn);
    state=0;
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

    size = count;

    if (std::is_same <T, double>::value) {
        data = laik_new_data(indexSpace, laik_Double );

    }
    else if (std::is_same <T, int>::value){
        data = laik_new_data(indexSpace, laik_Int64 );
    }
    // use the reservation API to precalculate the pointers

    Laik_Reservation* reservation = laik_reservation_new(data);
    laik_reservation_add(reservation, p1);

    laik_reservation_alloc(reservation);
    laik_data_use_reservation(data, reservation);

    // precalculate the transition object

    asW = laik_calc_actions(data, toW, reservation, reservation);
    asR = laik_calc_actions(data, toR, reservation, reservation);

    // go through the slices to just allocate the memory
    laik_switchto_partitioning(data, p1, LAIK_DF_None, reduction_operation );
    Laik_TaskSlice* ts = laik_my_slice(p1, 0);
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

    int i = (idx % (count*count) ) % count;
    int l_idx = i;
    int slice = idx/count;
    laik_map_def(data, slice, (void **)&base, &cnt);

    return base+l_idx;
}

template <typename T>
void laik_vector_overlapping<T>::calculate_pointers(){
    pointer_cache = (T**) malloc (size * sizeof(T*));

    for (int i = 0; i < size; ++i) {
        pointer_cache [i] = calc_pointer(i);
    }
}

template <typename T>
void laik_vector_overlapping<T>::switch_to_write_phase(){
    laik_exec_actions(asW);
    //laik_exec_transition(data,toW);
    //laik_switchto_phase(data, p1, Laik_DataFlow
    //                    ( LAIK_DF_ReduceOut | LAIK_DF_Sum ) );
}

template <typename T>
void laik_vector_overlapping<T>::switch_to_read_phase(){
    laik_exec_actions(asR);
    //laik_exec_transition(data,toR);
    //laik_switchto_phase(data, p1, LAIK_DF_CopyIn);
}

template <typename T>
void laik_vector_overlapping<T>::migrate(Laik_Group* new_group, Laik_Partitioning* p_new_1, Laik_Partitioning* p_new_2, Laik_Transition* t_new_1, Laik_Transition* t_new_2){
    uint64_t cnt;
    int* base;
    int slice = 0;

    laik_switchto_partitioning(data, p1, LAIK_DF_None, LAIK_RO_Min);

    Laik_Reservation* reservation = laik_reservation_new(data);
    laik_reservation_add(reservation, p_new_1);
    laik_reservation_alloc(reservation);
    laik_data_use_reservation(data, reservation);

    laik_switchto_partitioning(data, p_new_1, LAIK_DF_Preserve, LAIK_RO_Min);

    if (laik_myid(new_group)< 0) {
        return;
    }

    asW = laik_calc_actions(data, t_new_1, reservation, reservation);
    asR = laik_calc_actions(data, t_new_2, reservation, reservation);

    this -> p1=p_new_1;
    this -> p2=p_new_2;
    this -> toW=t_new_1;
    this -> toR=t_new_2;
    this -> world = new_group;
    if (laik_myid(world)<0)
        return ;

    laik_switchto_partitioning(data, p1, LAIK_DF_None, reduction_operation );
    int nSlices = laik_my_slicecount(p1);
    for (int n = 0; n < nSlices; n++)
    {
        laik_map_def(data, n, (void **)&base, &cnt);
    }
    laik_switchto_partitioning(data, p2, LAIK_DF_Preserve, LAIK_RO_Min);

    this -> count = cnt;
    this -> calculate_pointers();
}

template class laik_vector<double>;
template class laik_vector_halo<int>;
template class laik_vector_halo<double>;
template class laik_vector_overlapping<double>;
