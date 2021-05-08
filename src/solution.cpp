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

bool tolerance(dist_sort_t a,dist_sort_t b){
	dist_sort_t lowerbound = b*0.99999999;
	dist_sort_t upperbound = b*1.00000001;

	if(a>=lowerbound && a<=upperbound)
	{
		return true;
	}
	else
	{
		return false;
	}   
}

void findSplitters(const dist_sort_t *data, const dist_sort_size_t data_size, dist_sort_t *splitters, dist_sort_size_t *counts, int numSplitters) {
	int rank;
	int nProcs;
	int totalcount;
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs); //get number of processes
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);   //my rank


	for(int i=0;i<numSplitters;i++)
	{
		counts[i]=rank;
	}

	//dist_sort_t globalcount[numSplitters];
	MPI_Allreduce(&data_size,&totalcount, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
     
	

	dist_sort_t localmin,localmax,globalmin,globalmax;
	localmin = DIST_SORT_MAX;
	localmax = 0;

	for(dist_sort_size_t i=0;i<data_size;i++)
	{
		if(data[i]>localmax){localmax=data[i];}
		if(data[i]<localmin){localmin = data[i];}
	}
	 
	MPI_Allreduce(&localmin,&globalmin,1,MPI_UNSIGNED_LONG_LONG,MPI_MIN,MPI_COMM_WORLD);
	MPI_Allreduce(&localmax,&globalmax,1,MPI_UNSIGNED_LONG_LONG,MPI_MAX,MPI_COMM_WORLD);

	dist_sort_t probs[numSplitters];
	dist_sort_t L[numSplitters];
	dist_sort_t R[numSplitters];
	dist_sort_t globalcount[numSplitters];
	dist_sort_t localcount[numSplitters];
	dist_sort_t globalprefixcount[numSplitters];
	dist_sort_size_t dataperregion = totalcount/(numSplitters);
	probs[numSplitters-1]=globalmax;
	L[numSplitters]=globalmax;
	R[numSplitters]=globalmax;
	dist_sort_t isexpected = 0;
    
	if(rank==0)
	{
		for(int i = 0;i<numSplitters-1;i++) //init
		{
			probs[i]= (i+1)*(globalmax/nProcs);
			L[i]=0;
			R[i]=globalmax;
		}

	}

	while(isexpected ==0)
	{
        for(int i = 0;i<numSplitters;i++) //init
		{
			localcount[i]= 0;
			globalcount[i]=0;
		}

		MPI_Bcast(probs,numSplitters,MPI_UNSIGNED_LONG_LONG,0,MPI_COMM_WORLD);
		dist_sort_size_t probcounter = 0;
		
		for(dist_sort_size_t i=0;i<data_size;i++)
		{
			if(data[i]<=probs[probcounter])
			{
				localcount[probcounter]++;
			}
			else
			{
				probcounter++;
				localcount[probcounter]++;
			}     
		}

		MPI_Allreduce(&localcount,&globalcount,numSplitters,MPI_UNSIGNED_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);

		if(rank == 0)
		{ 
        	for(dist_sort_size_t i=0;i<numSplitters;i++)
			{
            	globalprefixcount[i]=(i>0)?globalcount[i]+globalprefixcount[i-1]:globalcount[i];
			}

        	isexpected = 1;
			for(dist_sort_size_t i=0;i<numSplitters-1;i++)
			{
				if(!tolerance(globalprefixcount[i],(i+1)*dataperregion ))
				{
                	isexpected=0;
					dist_sort_size_t leftprob = (i>0)?probs[i-1]:0;
					dist_sort_size_t rightprob = probs[i+1];
					if(globalprefixcount[i]>(i+1)*dataperregion) //move to left
					{
	   					R[i] = currentprob;
       					probs[i] = max(L,leftprob) / 2  + probs[i]/2;
   					}
   					else  //move to right
   					{
						L[i] = currentprob;
						probs[i] = min(R,rightprob) /2 + probs[i]/2;
					}
				}
			}
		}
		MPI_Bcast(&isexpected,1,MPI_UNSIGNED_LONG_LONG,0,MPI_COMM_WORLD);
	}

	MPI_Bcast(probs,numSplitters,MPI_UNSIGNED_LONG_LONG,0,MPI_COMM_WORLD);
	MPI_Bcast(globalcount,numSplitters,MPI_UNSIGNED_LONG_LONG,0,MPI_COMM_WORLD);

	for(dist_sort_size_t i = 0; i<numSplitters;i++)
	{
		counts[i]=globalcount[i];
		splitters[i]=probs[i];
	}
}

void moveData(const dist_sort_t *const sendData, const dist_sort_size_t sDataCount,
		dist_sort_t **recvData, dist_sort_size_t *rDataCount,
		const dist_sort_t *const splitters, const dist_sort_t *const counts, int numSplitters) {

	int nProcs;
	int rank;
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs); //get number of processes
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);   //my rank

	dist_sort_t localcounts[numSplitters];
	dist_sort_t myoffset[numSplitters];

	vector <vector<dist_sort_t>> sendtoprocess(numSplitters);


	for(dist_sort_size_t i = 0;i<numSplitters;i++) //init
		{
			localcounts[i]= 0;
			myoffset[i]=0;
		}
	dist_sort_size_t splittercounter = 0;
	for(dist_sort_size_t i=0;i<sDataCount;i++)
		{
			if(sendData[i]<=splitters[splittercounter])
			{
				localcounts[splittercounter] ++;
				sendtoprocess[splittercounter].push_back(sendData[i]);
			}
			else
			{
				splittercounter++;
				localcounts[splittercounter]++;
				sendtoprocess[splittercounter].push_back(sendData[i]);
				
			}     
		}

	dist_sort_t* receivedbuffer = (dist_sort_t *)(malloc(counts[rank] * sizeof(dist_sort_t)));

	MPI_Exscan(&localcounts,&myoffset,numSplitters,MPI_UNSIGNED_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);

	MPI_Win win;
	MPI_Win_create(receivedbuffer,counts[rank] * sizeof(dist_sort_t), sizeof(dist_sort_t),MPI_INFO_NULL, MPI_COMM_WORLD, &win);
	MPI_Win_fence(MPI_MODE_NOPRECEDE,win); 

	//init data sending

	for(dist_sort_size_t i = 0;i<numSplitters;i++)
	{
		MPI_Put(&sendtoprocess[i][0],sendtoprocess[i].size(), MPI_UNSIGNED_LONG_LONG,i,myoffset[i], sendtoprocess[i].size(), MPI_UNSIGNED_LONG_LONG, win);
	}
    
	MPI_Win_fence(0,win);
	MPI_Win_fence(MPI_MODE_NOSUCCEED,win); 

	*recvData = receivedbuffer;
	*rDataCount = counts[rank];

    MPI_Win_free(&win);
}

void sort(dist_sort_t *data, dist_sort_size_t size) {
	// You are welcome to use this sort function.
	// If you wish to implement your own, you can do that too.
	// Don't use bubblesort.
	std::sort(data,data+size);
}
