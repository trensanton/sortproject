#include <cmath>
#include <algorithm>
#include <iostream>
#include <utility>
#include <queue>
#include <vector>

#include "basic_defs.h"
#include "databasics.h"
#include "solution.h"

using namespace std;

void rebalance(const dist_sort_t *data, const dist_sort_size_t myDataCount, dist_sort_t **rebalancedData, dist_sort_size_t *rCount) {
	dist_sort_t myStartGlobal=0;
	int nProcs;
	dist_sort_t totalCount;
	int rank;

	MPI_Comm_size(MPI_COMM_WORLD, &nProcs); //get number of processes
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);   //MPI rank

	MPI_Exscan(&myDataCount,&myStartGlobal,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);  // finding total count  
	MPI_Allreduce(&myDataCount, &totalCount, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD); //finding local count
	
	dist_sort_t sizeDataPerProcess = ceil(double(totalCount)/double(nProcs));
	dist_sort_t* sharedData = (dist_sort_t *)(malloc(sizeDataPerProcess * sizeof(dist_sort_t)));

	MPI_Win win;
	MPI_Win_create(sharedData,sizeDataPerProcess * sizeof(dist_sort_t), sizeof(dist_sort_t),MPI_INFO_NULL, MPI_COMM_WORLD, &win);
	MPI_Win_fence(MPI_MODE_NOPRECEDE,win); 

	dist_sort_t mystart = myStartGlobal;
	dist_sort_t myend = myStartGlobal + myDataCount - 1;

	while(mystart<=myend)
	{
		dist_sort_t dest = mystart/sizeDataPerProcess;
		dist_sort_t displacement = mystart%sizeDataPerProcess;
		dist_sort_t size = std::min(sizeDataPerProcess*(dest+1),myend+1)-mystart;
		MPI_Put(&data[mystart-myStartGlobal], size, MPI_UNSIGNED_LONG_LONG, dest, displacement, size, MPI_UNSIGNED_LONG_LONG, win);
		mystart+=size;
	}

	MPI_Win_fence(0,win);
	MPI_Win_fence(MPI_MODE_NOSUCCEED,win); 
    
	MPI_Barrier(MPI_COMM_WORLD);

	int Count = 0;
	*rebalancedData = sharedData;
	if(rank==nProcs-1 && totalCount%sizeDataPerProcess!=0)
	{
		Count = totalCount%sizeDataPerProcess;
	}
	else
	{
		Count = sizeDataPerProcess;
	}

	*rCount = Count;

	MPI_Win_free(&win);
}

void findSplitters(const dist_sort_t *data, const dist_sort_size_t data_size, dist_sort_t *splitters, dist_sort_size_t *counts, int numSplitters) {

}

void moveData(const dist_sort_t *const sendData, const dist_sort_size_t sDataCount,
		dist_sort_t **recvData, dist_sort_size_t *rDataCount,
		const dist_sort_t *const splitters, const dist_sort_t *const counts, int numSplitters) {

}

void sort(dist_sort_t *data, dist_sort_size_t size) {
	// You are welcome to use this sort function.
	// If you wish to implement your own, you can do that too.
	// Don't use bubblesort.
	std::sort(data,data+size);
}
