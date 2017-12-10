#ifndef LAIK_PORT
#define LAIK_PORT

/*laik headers*/
extern "C"{
#include "laik.h"
#include "laik-backend-mpi.h"
}

int * build_element_corner_list(int edgeElems, int edgeNodes, int m_rowLoc, int m_colLoc, int m_planeLoc);
void free_local_corner_list(int* list);

void runExclusivePartitioner(Laik_Partitioner* pr, Laik_BorderArray* ba, Laik_BorderArray* otherBA);
Laik_Partitioner* exclusive_partitioner();
void runOverlapingPartitioner(Laik_Partitioner* pr, Laik_BorderArray* ba, Laik_BorderArray* otherBA);
Laik_Partitioner* overlaping_partitioner(int &depth);

#endif // LAIK_PORT

