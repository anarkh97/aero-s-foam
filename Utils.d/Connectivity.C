#include <iostream>
#include <cstdio>
#include <Utils.d/dbg_alloca.h>
#include <algorithm>
#include <set>
using std::stable_sort;

#include <Utils.d/Connectivity.h>
#include <Utils.d/BinFileHandler.h>
#include <Utils.d/resize_array.h>
#include <Utils.d/dofset.h>
#include <Element.d/Element.h>
#include <Mortar.d/FaceElement.d/FaceElemSet.h>
#include <Mortar.d/FaceElement.d/FaceElement.h>
#include <Element.d/Sommerfeld.d/SommerElement.h>

#include <Utils.d/GlobalInt.h>

//HB: for activating the optimization of renumByComponent without renumbering
//    Currently, one has to use renumAlg = -1 when calling the renumByComponent method
//    to use this optimization
#define HB_RENUMBYCOMP_OPT

Connectivity::Connectivity(int _size, int *_pointer, int *_target, int _removeable, float *_weight)
	: size(_size), pointer(_pointer, _pointer+_size+1),
	  target(_target, _target+_pointer[_size])
{
	size      = _size;
	numtarget = pointer[size];
	if(_removeable) {
		delete[] _pointer;
		delete[] _target;
		delete[] _weight;
	}
	if(_weight != nullptr)
		weight.assign(_weight, _weight+_pointer[_size]);
}

//Connectivity::Connectivity(const Elemset &els)
//{
//	int i;
//
//	size = els.last();
//
//	// Find out the number of targets we will have
//	pointer.resize(size+1) ;
//	int pp = 0;
//	for(i=0; i < size; ++i) {
//		pointer[i] = pp;
//		pp += els[i] ? els[i]->numNodes() : 0;
//	}
//	pointer[size] = pp;
//	numtarget = pp;
//
//	// Create the target array
//	target.resize(pp);
//
//	// Fill it in
//	for(i=0; i < size; ++i) {
//		if(els[i]) els[i]->nodes(target.data()+pointer[i]);
//	}
//}

Connectivity::Connectivity(const Elemset &els, Connectivity *nodeToElem)
{

	size = els.last();

	// locate any nodes that are not connected to rigid or flexible elements
	std::set<int> mpcnodes;
	for(int i=0; i < size; ++i) {
		Element *ele = els[i];
		if(ele->isMpcElement() && !ele->isRigidElement()) {
			int *nodes = els[i]->nodes();
			for(int j=0; j<ele->numNodes()-ele->numInternalNodes(); ++j) {
				bool connected = false;
				for(int k=0; k<nodeToElem->num(nodes[j]); ++k) {
					Element *kele = els[(*nodeToElem)[nodes[j]][k]];
					if(kele->isRigidElement() || !kele->isMpcElement()) {
						connected = true;
						break;
					}
				}
				if(!connected) mpcnodes.insert(nodes[j]);
			}
			delete [] nodes;
		}
	}

	size += mpcnodes.size();

	// Find out the number of targets we will have
	pointer.resize(size+1);
	int pp = 0;
	for(int i = 0; i < size-mpcnodes.size(); ++i) {
		pointer[i] = pp;
		Element *ele = els[i];
		if(ele)  {
			if(ele->isRigidElement() || !ele->isMpcElement()) {
				pp += ele->numNodes();
			}
		}
	}
	for(auto i=size-mpcnodes.size(); i < size; ++i) {
		pointer[i] = pp;
		pp++;
	}

	pointer[size] = pp;
	numtarget = pp;

	// Create the target array
	target.resize(pp);

	// Fill it in
	for(int i = 0; i < size-mpcnodes.size(); ++i) {
		Element *ele = els[i];
		if(ele) {
			if(ele->isRigidElement() || !ele->isMpcElement()) {
				ele->nodes(target.data()+pointer[i]);
			}
		}
	}
	auto i = size-mpcnodes.size();
	for (auto mpcnode : mpcnodes) {
		target[pointer[i++]] = mpcnode;
	}
}

Connectivity::Connectivity(Elemset*els, int numSom, SommerElement **som)
{
	int i, j, k;

	size = numSom;

	pointer.resize(numSom+1);
	for (i = 0 ; i < numSom+1 ; i++)
		pointer[i] = i;
	numtarget = numSom;

	target.resize(numSom);

	int nE = els->last();
	Connectivity *eToN = new Connectivity(els->asSet());
	Connectivity *nToE = eToN->reverse();
	int *eleCount = new int[nE];
	int nno, no, el, nei;

	for (i = 0 ; i < numSom ; i++) {
		el = -1;
		for (j = 0 ; j < nE ; j++) eleCount[j] = 0;
		//find the adjascent neighbour: el
		nno = som[i]->numNodes();
		for (j = 0 ; j < nno ; j++) {
			no = som[i]->getNode(j);
			for (k = 0 ; k < nToE->num(no) ; k++) {
				nei = (*nToE)[no][k];
				eleCount[nei]++;
				if (eleCount[nei]==nno) {
					el = nei;
					continue;
				}
			}
			if (el) continue;
		}
		if (el==-1) fprintf(stderr,"Connectivity.C, No adjascent element\n");
		target[i] = el;
	}
	delete [] eleCount; eleCount = NULL;
}

Connectivity::Connectivity(BinFileHandler& f, bool oldSower)
{
	if(!oldSower) {
		f.read(&size,1);
		pointer.resize(size);
		--size;
		f.read(pointer.data(),size+1);
		f.read(&numtarget,1);
		target.resize(numtarget);
		f.read(target.data(),numtarget);
	}
	else {
		//cerr << " *** WARNING: using OLD_SOWER ::Connectivity(BinFileHandler& f) \n";
		f.read(&size, 1);
		f.read(&numtarget, 1);
		pointer.resize(size+1);
		target.resize(numtarget);
		f.read(pointer.data(), size+1);
		f.read(target.data(), numtarget);
	}
}

size_t Connectivity::write(BinFileHandler& f)
{
	int _size = size+1;
	f.write(&_size,1);
	f.write(pointer.data(),_size);
	f.write(&numtarget,1);
	f.write(target.data(),numtarget);
	return 0;
}

size_t Connectivity::writeg(BinFileHandler& f)
{
	int _size = size+1;
	f.write(&_size,1);
	f.write(pointer.data(),_size);
	f.write(&numtarget,1);
	GlobalInt *gtarget = new GlobalInt[numtarget];
	for(int i = 0; i < numtarget; i++)
		gtarget[i] = target[i];
	f.write(gtarget,numtarget);
	delete [] gtarget;
	return 0;
}

Connectivity::Connectivity(int _size, int *_count)
{
	size    = _size;
	pointer.resize(size+1);
	pointer[0] = 0;
	int i;
	for(i=0; i < _size; ++i)
		pointer[i+1] = pointer[i] + _count[i];
	numtarget = pointer[size];
	target.resize(numtarget);
}

Connectivity::Connectivity(int _size, int count)
{
	size    = _size;
	pointer.resize(size+1);
	pointer[0] = 0;
	for(int i=0; i < _size; ++i)
		pointer[i+1] = pointer[i] + count;
	numtarget = pointer[size];
	target.resize(numtarget);
	for(int i=0; i < numtarget; ++i)
		target[i] = i;
}

Connectivity::Connectivity(FaceElemSet* els, int _size)
{
	int i;

	size = (_size == 0) ? els->last() : _size;

	// Find out the number of targets we will have
	pointer.resize(size+1);
	int pp = 0;
	for(i=0; i < size; ++i) {
		pointer[i] = pp;
		pp += (*els)[i] ? (*els)[i]->nNodes() : 0;
	}
	pointer[size] = pp;
	numtarget = pp;

	// Create the target array
	target.resize(pp);

	// Fill it in
	for(i=0; i < size; ++i) {
		if((*els)[i]) (*els)[i]->GetNodes(target.data()+pointer[i]);
	}
}

Connectivity::~Connectivity() {}

int
Connectivity::offset(int i, int j) const
{
	int ii;
	for(ii = pointer[i]; ii < pointer[i+1]; ++ii)
		if(target[ii] == j) return ii;

	return -1; // We didn't find a connection between i and j
}

int
Connectivity::cOffset(int i, int j) const
{
	int ii;
	for(ii = pointer[i]; ii < pointer[i+1]; ++ii)
		if(target[ii] == j) return (ii - pointer[i]);

	return -1; // We didn't find a connection between i and j
}

Connectivity*
Connectivity::transconOne( Connectivity* tc)
// D. Rixen :for every pointer
// of tc only one entry in target of tc is considered
// (used to associate a mpc term for a interface node
//  to only one sub, with a preference for already included sub)
{
	int i,j,k;

	// First find the biggest target so we can size arrays correctly
	int tgmax=-1;

	for(i =0; i < tc->numtarget; ++i)
		if(tc->target[i] > tgmax) tgmax=tc->target[i];
	tgmax++; // Important adjustment

	// Now we can size the array that flags if a target has been visited
	int *flags = (int *) dbg_alloca(sizeof(int)*tgmax);

	// For every pointer, build the number of occurrence of targets
	// and choose one target per tc->pointer (max occurrence)
	// At the same time build new pointers np

	int cp = 0;
	std::vector<int> np(size+1);
	int ii;
	for(i = 0; i < size; ++i) {
		np[i] = cp;
		for(ii = 0; ii < tgmax; ++ii)
			flags[ii] = 0;
		//-- store number of occurence in flag
		for(j = pointer[i]; j < pointer[i+1]; ++j){
			int intermed = target[j];
			for (k = 0; k < tc->num(intermed); ++k)
				flags[(*tc)[intermed][k]]++;
		}
		//-- set pointer and size target
		for(j = pointer[i]; j < pointer[i+1]; ++j){
			int intermed = target[j];
			//-- target with max occ.
			int targMaxOcc;
			int maxOcc = 0;
			for (k = 0; k < tc->num(intermed); ++k){
				if(flags[(*tc)[intermed][k]]==-1){
					maxOcc=0;
					break;
				}
				else if(flags[(*tc)[intermed][k]]> maxOcc){
					maxOcc = flags[(*tc)[intermed][k]];
					targMaxOcc = (*tc)[intermed][k];
				}
			}
			if(maxOcc!=0){cp++;flags[targMaxOcc]=-1;}
		}
	}
	np[size] = cp;
	// Now allocate and fill new target
	std::vector<int> ntg(cp);
	cp = 0;
	for(i = 0; i < size; ++i) {
		for(ii = 0; ii < tgmax; ++ii)
			flags[ii] = 0;
		//-- store number of occurence in flag
		for(j = pointer[i]; j < pointer[i+1]; ++j){
			int intermed = target[j];
			for (k = 0; k < tc->num(intermed); ++k)
				flags[(*tc)[intermed][k]]++;
		}
		//-- fill target
		for(j = pointer[i]; j < pointer[i+1]; ++j){
			int intermed = target[j];
			//-- target with max occ.
			int targMaxOcc;
			int maxOcc = 0;
			for (k = 0; k < tc->num(intermed); ++k){
				if(flags[(*tc)[intermed][k]]==-1){
					maxOcc=0;
					break;
				}
				else if(flags[(*tc)[intermed][k]]> maxOcc){
					maxOcc = flags[(*tc)[intermed][k]];
					targMaxOcc = (*tc)[intermed][k];
				}
			}
			if(maxOcc!=0){
				ntg[cp]=targMaxOcc;cp++;flags[targMaxOcc]=-1;}
		}
	}

	Connectivity *res = new Connectivity();
	res->size      = size;
	res->numtarget = np[size];
	res->pointer   = std::move(np);
	res->target    = std::move(ntg);
	return res;
}

int
Connectivity::num(int nd, int *mask) const
{
	int res=0;
	int jstrt = pointer[nd];
	int jstop = pointer[nd+1];
	int j;
	for(j = jstrt; j < jstop; ++j)
		if(mask[target[j]]) res++;
	return res;
}


void
Connectivity::findPseudoDiam(int *s, int *e, int *mask)
{
	int i,k,nw;
	// Select the node with the lowest connectivity
	int cmin = numtarget+1;
	int cmax = 0;
	for(i = 0; i < size; ++i) {
		if(mask[i]) {
			if((nw = num(i,mask)) < cmin) {
				cmin=nw;
				*s = i;
			}
			if(nw > cmax) cmax = nw;
		}
	}

	int *ls  =  new int[numtarget];
	int *xls =  new int[numtarget+1];

	// Created rooted level structure
	int w;
	int h = rootLS(*s, xls, ls, w, mask);
	int subconsize = xls[h];
	int *sorted = (int *) dbg_alloca((cmax+1)*sizeof(int));
	*e = ls[xls[h-1]]; // Give a default end point in case h == subconsize.
	while(h < subconsize) {
		for(k=0; k <= cmax; ++k) sorted[k] = -1;
		// Find in the last level the node with the minimum connectivity
		//int maxweight = subconsize;
		int kstrt = xls[h-1];
		int kstop = xls[h];
		for(k=kstrt; k < kstop; ++k)
		{
			sorted[num(ls[k],mask)] = ls[k];
		}
		int w_e = subconsize;
		for(k = 0; k <= cmax; ++k)
			if(sorted[k] >= 0) {
				int nh = rootLS(sorted[k], xls, ls, w, mask);
				if(w < w_e) {
					if(nh > h) {
						*s = sorted[k];
						h = nh;
						break;
					}
					*e = sorted[k];
					w_e = w;
				}
			}
		if(k > cmax) break;
	}
	delete [] xls;
	delete [] ls;
	return;
}

int
Connectivity::rootLS(int root, int *xls, int *ls, int &w, int *mask)
{
	int i, j;
	w = 0;

	int *locMask = new int[size];

	if(mask)
		for(i = 0; i < size; ++i) locMask[i] = mask[i];
	else
		for(i = 0; i < size; ++i) locMask[i] = 1;

	locMask[root] = 0;
	ls[0] = root;
	xls[0] = 0;
	int nlvl = 1;
	int nf = 1;
	while(nf > xls[nlvl-1]) {
		xls[nlvl] = nf;
		int lbegin = xls[nlvl-1];
		int lend   = xls[nlvl];
		for (i = lbegin; i <lend; ++i) {
			int n1 = ls[i];
			int jstart = pointer[n1];
			int jstop  = pointer[n1+1];
			for(j=jstart; j < jstop; ++j) {
				int n2 = target[j];
				if(locMask[n2]) {
					locMask[n2] = 0;
					ls[nf++] = n2;
				}
			}
		}
		if(nf-xls[nlvl] > w) w = nf-xls[nlvl];
		nlvl++;
	}

	delete [] locMask;

	return nlvl-1;
}

compStruct
Connectivity::renumByComponent(int renumAlg)
{
	// size = total number of nodes
	int *globalRenum = new int[size];
	int *mark = new int[size];
	int *ls   = new int[size];
	int *xls  = new int[size+1];
	ResizeArray<int> xcomp(0,2);
	// Initialize mark to zero, accounting for missing node #s
	// Initialize globalMask

	int inode;
	for(inode = 0; inode < size; ++inode) {
		mark[inode] = (num(inode) != 0) ? 1 : 0;
		globalRenum[inode] = -1;
	}

	// Loop over nodes checking which ones are marked
	// and belong to the same component.
	int j, k, nextNum = 0, count = 0;
	int *locMask = new int[size];
	for(inode = 0; inode < size; ++inode)
		locMask[inode] = 0;

	int currentNumber = 0;
	xcomp[0] = currentNumber;

	int *lrenum = 0;
	if(renumAlg>0) { lrenum = new int[size]; }

	for(inode = 0; inode < size; ++inode)
	{
		if(mark[inode] == 1) {
			// Find all neighbors of inode

			int w;
			int h = rootLS(inode, xls, ls, w);

			// Declare and set local mask

			for(j=0; j<xls[h]; ++j) {
				locMask[ls[j]] = 1;
				mark[ls[j]] = 0;
			}

			// call renumbering for local mask
			switch(renumAlg) {
				case 1:
					renumSloan(locMask, nextNum, lrenum);
					break;
				case 2:
					renumRCM(locMask, nextNum, lrenum);
					break;
				default:
					break;
			}

			// Assemble local mask into global mask
			if(renumAlg>0) {
				for(j=0; j<xls[h]; ++j) {
					k = ls[j];
					globalRenum[k] = lrenum[k];
				}
			}
			else {
#ifdef HB_RENUMBYCOMP_OPT
				if(renumAlg==-1){
					//HB: an optimization would be to loop over ls (0...xls[h])
					// but to recover the same answer as below, one has to sort
					// the ls array. so this may be not worthy ...
					// anyway, do we really need to sort the ls array ??? Not sorting
					// it will give us an output with a different ordering but we can
					// probably take care of it by providing an associated globalRenum array ???
					std::sort(ls,ls+xls[h]);
					for(j=0; j<xls[h]; ++j) {
						k = ls[j];
						globalRenum[k] = currentNumber + j;
					}
				}else{
#endif
					int cnt = 0;
					for(j=0; j<size; ++j)
						if(locMask[j])
							globalRenum[j] = currentNumber+cnt++;
#ifdef HB_RENUMBYCOMP_OPT
				}
#endif
			}
			currentNumber += xls[h];
			count += 1;
			xcomp[count] = currentNumber;
			// reset locMask to zero
			for(j=0; j<xls[h]; ++j)
				locMask[ls[j]] = 0;
		}
	}

	compStruct ret;
	ret.numComp = count;
	ret.xcomp   = xcomp.yield();
	ret.renum   = globalRenum;
	delete [] mark; delete [] ls; delete [] xls;
	delete [] locMask;
	if(lrenum) { delete [] lrenum; }
	return ret;
}

int *
Connectivity::renumRCM(int *mask, int &nextNum, int *renum)
{
	int i,j,k;
	if(mask == 0)
	{
		mask = (int *) dbg_alloca(sizeof(int)*size);
		for(i=0; i < size; ++i)
			mask[i] = (num(i)) ? 1 : 0;
	}

	int s_node, e_node;

	findPseudoDiam( &s_node, &e_node, mask);

	if(renum == 0) renum = new int[size];
	int *order  = (int *) dbg_alloca(sizeof(int)*size);
	int *degree = (int *) dbg_alloca(sizeof(int)*size);

	// get the degree of all the nodes
	for(i =0; i < size; ++i)
		degree[i] = num(i,mask);

	order[0] =e_node;
	mask[e_node] = -mask[e_node]; // mark we have seen this node
	int lastNode=1; // number of nodes which have been assigned a number
	for(i = 0; i < lastNode; ++i) {
		int curNode = order[i];
		int firstNeighb = lastNode;
		// Look at the neighbors of this node and add that to the list
		for(j = pointer[curNode]; j < pointer[curNode+1]; ++j) {
			int neighbNode = target[j];
			if(mask[neighbNode] > 0) {
				order[lastNode] = neighbNode;
				mask[neighbNode] = -mask[neighbNode];
				lastNode += 1;
			}
		}
		// now sort the added nodes by degree
		if(firstNeighb < lastNode-1) {
			for(j = firstNeighb; j < lastNode-1; ++j)
				for(k = j+1; k < lastNode; ++k)
					if(degree[order[k]] < degree[order[j]]) {
						int tmp = order[k];
						order[k] = order[j];
						order[k] = tmp;
					}
		}
	}
	// now reverse the order
	for(i = lastNode; i--; ) {
		renum[order[i]] = (nextNum++);
	}
	return renum;
}


int *
Connectivity::renumSloan(int *mask, int &nextNum, int *renum)
{
	int i,j,k;
	int s_node, e_node;
	float w1=1.0, w2=2.0;

	if(mask == 0)
	{
		mask = (int *) dbg_alloca(sizeof(int)*size);
		for(i=0; i < size; ++i)
			mask[i] = (num(i)) ? 1 : 0;
	}

	findPseudoDiam( &s_node, &e_node, mask);

	int *ls  = new int[numtarget];
	int *xls = new int[numtarget+1];

	int w;
	int h = rootLS(e_node, xls,ls,w,mask);
	// now give a distance to each point

	if(renum == 0) renum = new int[size];
	int *distance = renum;
	int *status = distance;

	for(i = 0; i < size; ++i)
		distance[i] = -1;

	for(i = 0; i < h; ++i) {
		for(j = xls[i]; j < xls[i+1]; ++j)
			distance[ls[j]] = i;
	}
	delete [] xls;
	delete [] ls;

	float *priority = (float *) dbg_alloca(sizeof(float)*size);
	// initialize the priority values
	for(i = 0; i < size; ++i) {
		if(mask[i]) {
			priority[i] = w1*distance[i] - w2*num(i);
			distance[i] = -2; // status and distance are stored in the same array
		}
	}

	// initalize the queue with the starting point
	int *q = (int *) dbg_alloca(sizeof(int)*size);
	// maximum size is all the points in the q
	q[0] = s_node;
	int nn=1;
	status[s_node] = -1;

	// While the queue is not empty
	while(nn > 0) {
		// Find in the queue the point with maximum priority
		int next = 0;
		float maxpri = priority[q[next]];
		for(i = 1; i < nn; ++i) {
			if(priority[q[i]] > maxpri) {
				maxpri = priority[q[i]];
				next = i;
			}
		}
		// Remove the next node numbered from the queue
		int nextnode = q[next];
		q[next] = q[nn-1];
		nn = nn-1;
		int istart = pointer[nextnode], istop = pointer[nextnode+1];
		if(status[nextnode] == -1) {
			status[nextnode] = 0; // So we don't step on our feet.
			// Preactive node. Examine its neighbors
			for(i = istart; i < istop; ++i) {
				int neighbor = target[i];
				priority[neighbor] += w2; // update the priority
				if(status[neighbor] == -2) { // if neighbor was inactive, it becomes pre-active
					// NB the next loop will set the same nodes as active
					q[nn] = neighbor;
					nn++;
					status[neighbor] = -1;
				}
			}
		}
		status[nextnode] = nextNum;
		nextNum++;
		// Now scan the preactive neighbors, make them active and their neighbors become
		// preactive if they were inactive
		for(i = istart; i < istop; ++i) {
			int neighbor = target[i];
			if(status[neighbor] == -1) {
				status[neighbor] = 0;
				priority[neighbor] += w2;
				int kstart = pointer[neighbor], kstop = pointer[neighbor+1];
				for(k = kstart; k < kstop ; ++k) {
					int kn = target[k];
					priority[kn] += w2;
					if(status[kn] == -2) { // This node is inactive. must become preactive
						status[kn] = -1;
						q[nn] = kn;
						nn++;
					}
				}
			}
		}
	}
	return status;
}

void
Connectivity::print(FILE *f, int node)
{
	if(node == -1) {
		int i;
		int maxdist = 0;
		fprintf(f, "Size %d\n",size);
		fflush(f);
		for(i = 0; i < size; ++i) {
			fprintf(f, "%d ->", i+1);
			int j;
			for(j = pointer[i]; j < pointer[i+1]; ++j) {
				fprintf(f, " %d", target[j]+1);
				if(i-target[j] > maxdist ) maxdist = i-target[j];
			}
			fprintf(f,"\n");
		}
		fprintf(f, "Max dist %d\n", maxdist);
		fflush(f);
	} else {
		fprintf(f,"%d ->",node+1);
		int j;
		for(j=pointer[node]; j<pointer[node+1]; ++j) {
			fprintf(f," %d", target[j]+1);
		}
		fprintf(f,"\n");
	}
}


int
Connectivity::findMaxDist(int *renum)
{
	int count = numNonZeroP();
	int i;
	int maxdist = 0;
	int localDist;
	double avgDist=0;
	for(i = 0; i < size; ++i) {
		localDist = 0;
		int j;
		for(j = pointer[i]; j < pointer[i+1]; ++j) {
			if(renum[i]-renum[target[j]] > maxdist )
				maxdist = renum[i]-renum[target[j]];
			if(renum[i]-renum[target[j]] > localDist)
				localDist = renum[i]-renum[target[j]];
		}
		avgDist += localDist;
	}
	fprintf(stderr, "Total   distance: %e\n", avgDist);
	fprintf(stderr, "Average distance: %e\n", avgDist/count);
	fprintf(stderr, "Maximum distance: %d\n", maxdist);
	return maxdist;
}

int
Connectivity::findProfileSize(EqNumberer *eqn, int unroll)
{
	int i;
	int profileSize = 0;
	int localSize   = 0;
	for(i = 0; i < size; ++i)
	{
		localSize  = 0;
		int fdof       = eqn->firstdof(i);
		int localWidth = eqn->weight(i);

		int j;
		for(j = pointer[i]; j < pointer[i+1]; ++j)
		{
			int fjdof = eqn->firstdof(target[j]);
			if(fjdof < 0) continue;

			// This is for skyline unrolling of level unroll
			fjdof = fjdof - (fjdof % unroll);

			if(fdof - fjdof > localSize)
				localSize = fdof - fjdof;
		}
		localSize   *= localWidth;
		localSize   += ((localWidth*(localWidth+1))/2);
		profileSize += localSize;
	}
	// fprintf(stderr,"---------------------\n");
	fprintf(stderr, "Profile Size: %d\n", profileSize);

	return profileSize;

}

Connectivity*
Connectivity::merge(Connectivity *con2)
{
	int size1 = csize();
	int size2 = con2->csize();

	int i;
	int *cp = new int[size1 + size2+1];
	int fp = 0;
	for(i = 0; i < size1; ++i) {
		cp[i] = fp;
		fp += num(i);
	}
	for(i = 0; i < size2; ++i) {
		cp[i+size1] = fp;
		fp += con2->num(i);
	}
	cp[size1 + size2] = fp;

	int *ct = new int[fp];
	fp = 0;
	for(i=0; i<size1; ++i) {
		int j;
		for(j =0; j < num(i); ++j)
			ct[fp++] = (*this)[i][j];
	}
	for(i = 0; i < size2; ++i) {
		int j;
		for(j =0; j < con2->num(i); ++j)
			ct[fp++] = (*con2)[i][j];
	}

	return new Connectivity(size1+size2, cp, ct);
}

void
Connectivity::combine(Connectivity *con2, int *&cmap, int *cmap2)
{
	// add con2 to this and make combined cmap
	int i, j;
	int size1 = csize();
	int size2 = con2->csize();
	int *tmp1 = (int *) dbg_alloca(sizeof(int)*size1);
	for(i=0; i<size1; ++i) tmp1[i] = -1;
	int *tmp2 = (int *) dbg_alloca(sizeof(int)*size2);
	for(i=0; i<size2; ++i) tmp2[i] = -1;

	std::vector<int> cp(size1 + size2+1);
	int fp = 0;
	for(i = 0; i < size1; ++i) {
		cp[i] = fp;
		fp += num(i);
		for(j = 0; j < size2; j++) {
			if(cmap2[j] == cmap[i]) {
				tmp2[j] = i;
				tmp1[i] = j;
				fp += con2->num(j);
				break;
			}
		}
	}
	int count = 0;
	for(i = 0; i < size2; ++i) {
		if(tmp2[i] == -1) {
			cp[count+size1] = fp;
			fp += con2->num(i);
			count++;
		}
	}
	cp[size1 + count] = fp;
	int *new_cmap = new int[size1 + count];

	std::vector<int> ct(fp);
	fp = 0;
	for(i=0; i<size1; ++i) {
		new_cmap[i] = cmap[i];
		for(j =0; j < num(i); ++j)
			ct[fp++] = (*this)[i][j];
		if(tmp1[i] != -1) {
			for(j =0; j < con2->num(tmp1[i]); ++j)
				ct[fp++] = -1 - (*con2)[tmp1[i]][j];  // merged are -ve
		}
	}
	count = 0;
	for(i = 0; i < size2; ++i) {
		if(tmp2[i] == -1) {
			new_cmap[size1+count] = cmap2[i];
			count++;
			for(j =0; j < con2->num(i); ++j)
				ct[fp++] = -1 - (*con2)[i][j]; // merged are -ve
		}
	}
	delete [] cmap;
	cmap = new_cmap;


	size      = size1+count;
	pointer   = cp;
	target    = ct;
	numtarget = cp[size1+count];
}

int
Connectivity::numNonZeroP() const
{
	int count = 0;
	int i;
	for(i=0; i < size; ++i)
		if(pointer[i+1] != pointer[i]) ++count;
	return count;
}

// Collapse creates a connecivity that represents the grouping
// of a reflexive connectivity into subconnected groups 
Connectivity *
Connectivity::collapse() {
	int i,j;
	int count = 0;
	bool *flag = new bool[size];
	int grIndex;
	int *group = new int[size];
	for(i=0; i < size; ++i)
		flag[i] = false;
	grIndex = 0;
	for(i = 0; i < size; ++i)
		if(flag[i] == false) {
			count++;
			int follow = grIndex;
			flag[i] = true;
			group[grIndex++] = i;
			while(follow < grIndex) {
				int s = group[follow];
				for(j = 0; j < num(s); ++j)
					if((*this)[s][j] >= 0 && flag[ (*this)[s][j] ] == false) {
						flag[ (*this)[s][j] ] = true;
						group[grIndex++] = (*this)[s][j];
					}
				follow++;
			}
		}
	int *ptr = new int[count+1];
	count = 0;
	grIndex = 0;
	for(i=0; i < size; ++i)
		flag[i] = false;
	for(i = 0; i < size; ++i)
		if(flag[i] == false) {
			ptr[count++] = grIndex;
			int follow = grIndex;
			flag[i] = true;
			group[grIndex++] = i;
			while(follow < grIndex) {
				int s = group[follow];
				for(j = 0; j < num(s); ++j)
					if((*this)[s][j] >= 0 && (*this)[s][j] <size && flag[ (*this)[s][j] ] == false) {
						flag[ (*this)[s][j] ] = true;
						group[grIndex++] = (*this)[s][j];
					}
					else if((*this)[s][j] >= size)
						fprintf(stderr, "error connectivity is not reflexive\n");
				follow++;
			}
		}
	ptr[count] = grIndex;
	delete [] flag;
	return new Connectivity(count, ptr, group);
}

Connectivity *
Connectivity::copy()
{
	if(size==0) return new Connectivity();
	int *nptr = new int[size+1];
	int *ntrg = new int[numConnect()];
	int i;
	for(i = 0; i <= size; ++i)
		nptr[i] = pointer[i];
	for(i = 0; i < numConnect(); ++i)
		ntrg[i] = target[i];
	return new Connectivity(size,nptr,ntrg);
}

// Trim a connectivity to remove targets that have only one corresponding
// target in the marker connectivity
Connectivity *
Connectivity::trim(Connectivity *marker)
{
	int i,j;
	int count;
	count = 0;
	for(i = 0; i < numConnect(); ++i)
		if(marker->num( target[i] ) > 1)
			count++;
	std::vector<int> ntrg(count);
	std::vector<int> nptr(size+1);
	count = 0;
	for(i = 0; i < size; ++i) {
		nptr[i] = count;
		for(j = 0; j < num(i); ++j)
			if( marker->num( (*this)[i][j] ) > 1)
				ntrg[count++] = (*this)[i][j];
	}
	nptr[size] = count;
	return new Connectivity(size, nptr, ntrg, std::vector<float>());
}

// PJSA: modify a connectivity to insert connection between each target and itself 
// needed to make bodyToBody, since some bodys may have no MPCs
Connectivity *
Connectivity::modify()
{
	int i,j;
	int count;
	count = 0;
	for(i = 0; i < size; ++i)
		if(num(i) == 0) count++;
	int *ntrg = new int[numtarget+count];
	int *nptr = new int[size+1];
	count = 0;
	for(i = 0; i < size; ++i) {
		nptr[i] = count;
		if(num(i) == 0) ntrg[count++] = i;
		else {
			for(j = 0; j < num(i); ++j)
				ntrg[count++] = (*this)[i][j];
		}
	}
	nptr[size] = count;
	return new Connectivity(size, nptr, ntrg);
}

// this one is to remove connection with self
Connectivity *
Connectivity::modifyAlt()
{
	int i,j;
	int count;
	count = 0;
	for(i = 0; i < size; ++i)
		for(j = 0; j < num(i); ++j)
			if((*this)[i][j] == i) count++;
	int *ntrg = new int[numtarget-count];
	int *nptr = new int[size+1];
	count = 0;
	for(i = 0; i < size; ++i) {
		nptr[i] = count;
		for(j = 0; j < num(i); ++j)
			if((*this)[i][j] != i) ntrg[count++] = (*this)[i][j];
	}
	nptr[size] = count;
	return new Connectivity(size, nptr, ntrg);
}

Connectivity *
Connectivity::combineAll(int addSize, int* cmap)
{
	int i,j,k;
	int count;
	bool foundi, foundj;
	//count = 0;
	//for(i = 0; i < size; ++i)
	//if(num(i) == 0) count++;
	int *ntrg = new int[numtarget+(addSize)*(addSize)];
	int *nptr = new int[size+1];

	count = 0;
	for(i = 0; i < size; ++i) {
		nptr[i] = count;
		//if(num(i) == 0) ntrg[count++] = i;
		//else {
		for(j = 0; j<num(i); ++j)
			ntrg[count++] = (*this)[i][j];
		//}
		int oldNumI = num(i);
		//fprintf(stderr, " ... i = %d, cmap[addCount] = %d ...\n", i,cmap[addCount]);
		foundi = false;
		for(k = 0; k < addSize; ++k)
			if(cmap[k] == i) { foundi = true; break; }
		if(foundi) {
			for(k = 0; k < addSize; ++k) {
				foundj = false;
				for(j = 0; j < oldNumI; ++j) {
					//fprintf(stderr, " ... i = %d, j = %d, k = %d, cmap[k] = %d, (*this)[i][j] = %d ...\n", i,j,k,cmap[k],(*this)[i][j]);
					if (cmap[k] == (*this)[i][j]) {
						foundj = true;
						//fprintf(stderr, " ... found!!! : i = %d, j = %d, k = %d ...\n", i,j,k);
						break;
					}
				}
				if (!foundj) {
					//fprintf(stderr, " ... NOT found!!! : i = %d, j = %d, k = %d ...\n", i,j,k);
					ntrg[count++] = cmap[k];
				}
			}
		}
	}
	nptr[size] = count;
	return new Connectivity(size, nptr, ntrg);
}

void
Connectivity::sortTargets()
{
	for(int i = 0; i < size; ++i)
		stable_sort(target.data()+pointer[i], target.data()+pointer[i+1]);
}

//----------------------------------------------------------------------

void Connectivity::renumberTargets(int *map)  {

	for (int i = 0; i < numtarget; i++)  {
		if (map[target[i]] < 0)
			fprintf(stderr, "target(i) is neg: %d\n", target[i]);
		target[i] = map[target[i]];
	}
}

//----------------------------------------------------------------------

void Connectivity::renumberTargets(std::map<int, int> &map)  {

	for (int i = 0; i < numtarget; i++)  {
		if (map.find(target[i]) == map.end())
			fprintf(stderr, "target(i) does not exist: %d\n", target[i]);
		target[i] = map[target[i]];
	}
}

//----------------------------------------------------------------------

Connectivity *
Connectivity::subSection(bool *select)
{
	int i, j;
	int tgcnt = 0;
	for(i = 0; i < size; ++i)
		if(select[i])
			tgcnt += num(i);
	int *newPtr  = new int[size+1];
	int *newTrg = new int[tgcnt];
	// fprintf(stderr, "Size is %d\n", tgcnt);
	tgcnt = 0;
	for(i = 0; i < size; ++i) {
		newPtr[i] = tgcnt;
		if(select[i])
			for(j = 0; j < num(i); ++j)
				newTrg[tgcnt++] = (*this)[i][j];
	}
	newPtr[size] = tgcnt;
	return new Connectivity(size, newPtr, newTrg);
}

CountedConnectivity::CountedConnectivity(int ns) : Connectivity(ns) {
	cnt = new int[csize()];
}

CountedConnectivity::~CountedConnectivity()
{
	delete [] cnt;
}

void
CountedConnectivity::end_count()
{
	int i;
	for(i=0; i < csize(); ++i)
		cnt[i] = pointer[i];
	Connectivity::end_count();
}
void
Connectivity::end_count()
{
	for(int i=0; i < size; ++i)
		pointer[i+1] += pointer[i];
	numtarget = pointer[size];
	target.resize(numtarget);
}

void
Connectivity::countlink(int from, int)
{
	pointer[from]++;
}

void
Connectivity::addlink(int from, int to)
{
	int t = --pointer[from];
	target[t] = to;
}

void
CountedConnectivity::remove(int from, int to)
{
	cnt[from]--;
	int i, ib, ie;
	ib = pointer[from];
	ie = ib+cnt[from];
	for(i =ib; i < ie; ++i)
		if(target[i] == to) break;
	for(; i < ie; ++i)
		target[i] = target[i+1];
}

Connectivity::Connectivity(int ns)
{
	size = ns;
	pointer.resize(size+1);
	int i;
	for(i=0; i < size+1; ++i)
		pointer[i] = 0;
	numtarget = 0;
}

//----------------------------------------------------------------------
// HB
//----------------------------------------------------------------------
// WARNING: ASSUME 1 DOF PER NODE !!!
double
Connectivity::estimateComponentCost(EqNumberer *eqn, compStruct &cs, double *cost, double *bandwidth, double coef, int unroll)
{
	//fprintf(stderr," Enter Connectivity::estimateCompCost\n");
	int nComp = cs.numComp;
	int i,j;
	double totalCost = 0.0;
	for(i=0;i<nComp;i++){
		cost[i] = 0.0;
		bandwidth[i] = 0.0;
		int ndof = 0;
		for(j=cs.xcomp[i];j<cs.xcomp[i+1];j++){
			//int nodej = cs.renum[j];
			int nodej = cs.order[j];

			int localSize  = 0;
			int fdof       = eqn->firstdof(nodej);
			int localWidth = eqn->weight(nodej);
			//ndof += localWidth;
			ndof++;
			int k;
			for(k = pointer[nodej]; k < pointer[nodej+1]; ++k) {
				int fjdof = eqn->firstdof(target[k]);
				if(fjdof < 0) continue;

				// This is for skyline unrolling of level unroll
				fjdof = fjdof - (fjdof % unroll);

				//fprintf(stderr," fdof - fjdof = %d\n",fdof - fjdof);
				//fprintf(stderr," fjdof - fdof = %d\n",fjdof - fdof);
				if(fdof - fjdof > localSize)
					localSize = fdof - fjdof;
				//if(fjdof - fdof > localSize)
				//  localSize = fjdof - fdof;
			}
			double localSizefb = (localSize) ? localSize : localSize+1;
			cost[i] += localSize*localSize*localWidth + ((localWidth*(localWidth+1))/2) + coef*localSizefb*localWidth;
			//bandwidth[i] += localSize*localSize;
			bandwidth[i] += localSizefb*localSizefb;
		}
		bandwidth[i] = sqrt(bandwidth[i]/ndof);
		totalCost += cost[i];
	}
	return totalCost;
}

double
Connectivity::estimateCost(EqNumberer *eqn, double &cost, double &bandwidth, double coef, int unroll)
{
	//fprintf(stderr," Enter Connectivity::estimateCost\n");
	int i,j;
	double totalCost = 0.0;
	cost = 0.0;
	bandwidth = 0.0;
	int ndof = 0;
	// int profileSize = 0;
	int localSize   = 0;
	for(i = 0; i < size; ++i) {
		localSize  = 0;
		int fdof       = eqn->firstdof(i);
		int localWidth = eqn->weight(i);
		ndof++;

		for(j = pointer[i]; j < pointer[i+1]; ++j) {
			int fjdof = eqn->firstdof(target[j]);
			if(fjdof < 0) continue;

			// This is for skyline unrolling of level unroll
			fjdof = fjdof - (fjdof % unroll);

			if(fdof - fjdof > localSize)
				localSize = fdof - fjdof;
		}
		double localSizefb = (localSize) ? localSize : localSize+1;
		cost += localSize*localSize*localWidth + ((localWidth*(localWidth+1))/2) + coef*localSizefb*localWidth;
		//bandwidth += localSize*localSize;
		bandwidth += localSizefb*localSizefb;
	}
	bandwidth = sqrt(bandwidth/ndof);
	totalCost += cost;

	return totalCost;
}


Connectivity::Connectivity(SommerElement  **els, int numSommer)
{
	int i;

	size = numSommer;

	// Find out the number of targets we will have
	pointer.resize(size+1);
	int pp = 0;
	for(i=0; i < size; ++i) {
		pointer[i] = pp;
		pp += els[i] ? els[i]->numNodes() : 0;
	}
	pointer[size] = pp;
	numtarget = pp;

	// Create the target array
	target.resize(pp);

	// Fill it in
	for(i=0; i < size; ++i)
		if( els[i]) els[i]->nodes(target.data()+pointer[i]);
}

Connectivity::Connectivity(SommerElement  **els, int numSommer, int max)
{
	int i;

	size = numSommer;

	// Find out the number of targets we will have
	pointer.resize(size+1);
	int pp = 0;
	int nmax = 0;
	for(i=0; i < size; ++i) {
		pointer[i] = pp;
		int nn = els[i] ? els[i]->numNodes() : 0 ;
		if (nn > nmax) nmax = nn;
		pp += (nn > max) ? max : nn;
	}
	pointer[size] = pp;
	numtarget = pp;

	// Create the target array
	target.resize(pp);

	int *buffer = (int*) dbg_alloca(nmax*sizeof(int));

	// Fill it in
	for(i=0; i < size; ++i)
		if(els[i]) {
			els[i]->nodes(buffer);
			for(int j=0;j<pointer[i+1]-pointer[i] ;j++)
				target[pointer[i]+j] = buffer[j];
		}
}


Connectivity::Connectivity(FILE *f, int n)
{
	numtarget = n;
	target.reserve(numtarget);

	int r1 = fscanf(f,"%d",&size);

	pointer.reserve(size+1);

	for (int i = 0; i < size; ++i) {
		int m;
		int r2 = fscanf(f,"%d",&m);
		pointer.push_back(target.size());
		if(n != 0 && target.size() + m > numtarget) {
			fprintf(stderr," *** ERROR: Connectivity has too many connections\n");
			exit(1);
		}
		for (int j = 0; j < m; ++j) {
			int t;
			int r3 = fscanf(f,"%d",&t);
			target.push_back(t-1);
		}
	}
	pointer.push_back(target.size());
}

Connectivity::Connectivity(int _size, std::vector<int> _pointer, std::vector<int> _target,
                           std::vector<float> _weight) : size(_size), numtarget(_target.size()),
                                                         pointer(std::move(_pointer)),
                                                         target(std::move(_target)), weight(std::move(_weight)){

}
