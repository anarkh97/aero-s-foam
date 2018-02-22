//
// Created by Michel Lesoinne on 2/13/18.
//
#include <mpi.h>
#include <iostream>
#include <array>
#include <vector>
#include <map>
#include <chrono>
#include "BaseCommunicator.h"

auto now() { return std::chrono::high_resolution_clock::now(); }

using timepoint = decltype(now());
double elapsed(timepoint t0, timepoint t1) { return std::chrono::duration<double>{t1-t0}.count();}

class SubD {
public:
	SubD(std::array<int,3> subSize, std::array<int,3> nsubPerDirection, std::array<int, 3> thisIndex);

	std::map<int, std::vector<long> > getNodeSharing(const BaseCommunicator &communicator) const;
private:
	std::vector<long> nodeIndices;
	long totalNumNodes;
public:
	const std::vector<long> &getNodeIndices() const;

	long getTotalNumNodes() const;
	//!< Number of nodes in the whole model.
};

template <typename T>
T product(const std::array<T,3> &d) {
	return d[0]*d[1]*d[2];
}

template <typename T>
T productP1(const std::array<T,3> &d) {
	return (d[0]+1)*(d[1]+1)*(d[2]+1);
}

SubD::SubD(std::array<int, 3> subSize, std::array<int, 3> nsubPerDirection, std::array<int, 3> thisIndex) {
	std::array<long,3> totalSize { subSize[0]*nsubPerDirection[0]+1,
	                               subSize[1]*nsubPerDirection[1]+1,
	                               subSize[2]*nsubPerDirection[2]+1 };

	totalNumNodes = productP1(totalSize);
	nodeIndices.reserve(productP1(subSize));

	std::array<long,3> corner{ subSize[0]*thisIndex[0],
	                           subSize[1]*thisIndex[1],
	                           subSize[2]*thisIndex[2] };

	auto index = [&] (long i, long j, long k) {
		return i+(totalSize[0]*(j+totalSize[1]*k));
	};

	for(long i = 0; i <= subSize[0]; ++i)
		for(long j = 0; j <= subSize[1]; ++j)
			for(long k = 0; k <= subSize[2]; ++k)
				nodeIndices.push_back(index(i+corner[0], j+corner[1], k+corner[2]));
}

/** \brief Obtain the list of processes who are canonic owners of the nodes and for each owner the list of nodes.
 *
 * @param totalNumNodes Maximum index+1 of all the nodes in all the processors.
 * @param glNodes Array of global node numbers.
 * @param communicator Communicator over which the processors taking park communicate.
 * @return A map between process index and the array of global node indices whose canonic owner is that process.
 */
std::map<int, std::vector<long> >
getCanonicalOwners(long totalNumNodes, const std::vector<long> &glNodes, int numProc) {

	auto canonicOwner = [ nPerProc=totalNumNodes/numProc, maxRem =  totalNumNodes%numProc](long node) {
		auto baseProc = node/nPerProc;
		auto rem = node % nPerProc;
		if(rem < maxRem && rem < baseProc-1)
			return baseProc-1;
		return baseProc;
	};

	std::map<int, std::vector<long>> canonicNodePosition;
	for(auto node: glNodes) {
		auto &list = canonicNodePosition[canonicOwner(node)];
		list.push_back(node);
	}
	return canonicNodePosition;
};

std::map<int, std::vector<long> > SubD::getNodeSharing(const BaseCommunicator &communicator) const {
	int numProc = communicator.commSize();

	std::map<int, std::vector<long>> canonicNodePosition = getCanonicalOwners(totalNumNodes,
	                                                                          nodeIndices, communicator.commSize());
	long info[2] = {0,0};
	auto window = communicator.window(info, 2);
	window.open();
	long one = 1;
	for(auto &info : canonicNodePosition)
		window.accumulate(SumHandle, &one, 1, info.first, 0);
	window.close();

	return std::map<int, std::vector<long>>();
}

const std::vector<long> &SubD::getNodeIndices() const {
	return nodeIndices;
}

long SubD::getTotalNumNodes() const {
	return totalNumNodes;
}

void testGlobals(const BaseCommunicator &communicator) {
	communicator.barrier();
	double result;
	auto t1 = now();
	for(int i = 0; i < 100; ++i)
		result += communicator.globalSum(i);
	auto t2 = now();
	if(communicator.rank() == 0)
		std::cout << "100 Global sum took: " << elapsed(t1,t2)*1e6 << "µs. Result is " << result << "." << std::endl;
}


int main(int argc, char *argv[]) {
	int required = MPI_THREAD_MULTIPLE, provided;
	MPI_Init_thread(&argc, &argv, required, &provided);
	BaseCommunicator comm(MPI_COMM_WORLD);
	int myRank = comm.rank();
	int commSize = comm.commSize();
	if(myRank == 0)
		std::cout << "I was requiring " << required << " I got " << provided << std::endl;

	testGlobals(comm);
	int databank[4096];
	int &myData = databank[0];
	myData = 2*myRank;
	int remData;
	int operandData = myRank+3;

	int myDest = (myRank+1)%commSize;

	auto t0 = now();
	{
		auto window = comm.window(&myData, 4096);

		window.sharedLockAll();
		window.fetchAndOp(ProdHandle, &operandData, &remData, myDest, 0);
		window.unlockAll();

	}
	auto dt = elapsed(t0, now());
	if(myRank == 0)
		std::cout << "Time to do the work is " << dt*1e6 << "µs." << std::endl;
	for(int i = 0; i < commSize; ++i) {
		comm.barrier();
		if(myRank == i)
			std::cout << i << " worked with " << myDest << " where I found " << remData << " and I was made to have " << myData << std::endl;
	}
	MPI_Finalize();
}
