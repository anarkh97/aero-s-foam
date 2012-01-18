#include <Element.d/Element.h>
#include <iostream>
#include <set>
#include <queue>
#include <Utils.d/Connectivity.h>

void Elemset::deleteElem(int i)
{
  elem[i]->~Element();
  elem[i] = 0;
}

void Elemset::deleteElems() 
{
  if (!elem) {
    return;
  }

  if (myData) { 
    for (int i = 0; i < size(); ++i) {
      if (elem[i]) elem[i]->~Element();
    }
  }
  
  delete[] elem;
  elem = 0;
  emax = 0;
}

// Last returns the true last defined element
int Elemset::last() const
{
 int last = size();
 while(--last >= 0)
    if(elem[last] != 0) break;

 return last+1;
}


void Elemset::list()
{
  std::cerr << "Size of Elemset : " << emax << std::endl;
  for(int i =0;i<emax ; i++)
    {
      if(elem[i])
	{
	  std::cerr << i <<": ";
	  int num = elem[i]->numNodes();
	  int* nodes = elem[i]->nodes();
	  std::cerr << " nodes : ";
	  for(int j=0;j<num;++j)
      	    std::cerr << nodes[j]<< ",";
	  delete [] nodes;
#ifdef SOWER_DEBUG
          std::cerr << " pressure : " << elem[i]->getPressure();
          std::cerr << " preload : " << elem[i]->getPreLoad();
#endif
          std::cerr << std::endl;
	}
      else
	{
	  std::cerr<<".";
	  std::cerr << "empty element" << std::endl;
	}
    }
  std::cerr << std::endl;
}

typedef Element *ElemP;

Elemset::Elemset(int initsize) : ba(4096, initsize)
{
  emax = initsize;
  if(emax > 0) {
    elem = new ElemP[emax];
    for(int i = 0; i < emax; ++i) elem[i] = 0;
  }
  else elem = 0;
  myData = false;
}

Elemset::Elemset(Elemset &globalset, int numlocal, int *localToGlobalMap)
{
  emax = numlocal;
  if(emax > 0) {
    elem = new ElemP[emax];
    for(int i = 0; i < emax; ++i) elem[i] = globalset[localToGlobalMap[i]];
  }
  else elem = 0;
  myData = false;
}

void
Elemset::elemadd(int num, Element *el)
{
  if(num >= emax) // resize elem[]
   {
    int newsize = ((num+1)*3)/2;
    ElemP *np = new ElemP[newsize];
    int i;
    for(i= 0; i < emax; ++i)
     np[i] = elem[i];
    for(emax = newsize; i < emax; ++i)
     np[i] = 0;
    if(elem) delete [] elem;
    elem = np;
   }
  if(elem[num]) cerr << " *** WARNING: found repeated ELEMENT# " << num+1 << endl; 
  elem[num] = el;
}

Element *& Elemset::operator[] (int num)
{
  if(num >= emax) // resize elem[]
   {
    int newsize = ((num+1)*3)/2;
    ElemP *np = new ElemP[newsize];
    int i;
    for(i= 0; i < emax; ++i)
     np[i] = elem[i];
    for(emax = newsize; i < emax; ++i)
     np[i] = 0;
    if(elem) delete [] elem;
    elem = np;
   }
  return elem[num];
}

template<>
class SetAccess<std::pair<Elemset, DofSet> > 
{
    Elemset &els;
    DofSet ds;
  public:
    SetAccess(Elemset &eset, DofSet dofset) : els(eset), ds(dofset) {}
    int size() {
      return els.last();
    }
    int numNodes(int i) {
      if(els[i] == 0)
        return 0;
      return els[i]->numNodes() - els[i]->numInternalNodes();
    }
    void getNodes(int i, int *nd) {
      if(els[i]) {
        int *p = new int[els[i]->numNodes()];
        els[i]->nodes(p);
        for(int j = 0; j < els[i]->numNodes() - els[i]->numInternalNodes(); ++j)
          nd[j] = p[j];
        delete [] p;
      }
    }
};

#include <Element.d/Rigid.d/RigidBeam.h>
void
Elemset::collapseRigid6(std::set<int> &blockedNodes)
{
  SetAccess<std::pair<Elemset, DofSet> > sa(*this, DofSet(DofSet::DispAndRot));
  Connectivity eToN(sa);
  Connectivity *nToE = eToN.reverse();
  Connectivity *nToNRigid = nToE->transcon(&eToN);

  Element **newSet = new Element *[last()];
  int iEl = 0;

  // Examine all the nodes to form groups and for each group, collapse the elements
  // that form it into a single group.
  std::set<int> visitedNodes;
  int _last = last();
  for(int i = 0; i < nToNRigid->csize(); ++i) {
    // Skip nodes that have no connection or that have already been treated
    if(nToNRigid->num(i) == 0 || visitedNodes.find(i) != visitedNodes.end())
      continue;
    std::set<int> nodeGroup, elementGroup;
    std::queue<int> fifo;
    for(fifo.push(i); !fifo.empty(); fifo.pop()) {
      int j = fifo.front();
      if(visitedNodes.insert(j).second == false)
        continue; // We've already examined all the neighbors of this node
      //std::cerr << "Starting at " << j+1 << ": ";

      for(int ik = 0; ik < nToNRigid->num(j); ++ik) {
        int k = (*nToNRigid)[j][ik];
        // If this node has not been seen before, then we can
        // add it to the queue.
        if(nodeGroup.insert(k).second)
          fifo.push(k);
        //  std::cerr << " " << k+1 << "(" << nodeGroup.size()<<")";
      }
      // std::cerr << std::endl;
      for(int el = 0; el < nToE->num(j); ++el)
        elementGroup.insert((*nToE)[j][el]);
    }
    // We now have the group of nodes that form a complete rigid body.
    // And the group of elements. 

    // choose master node
    int masterNode = *nodeGroup.begin();
    for(std::set<int>::iterator it = blockedNodes.begin(); it != blockedNodes.end(); ++it) {
      if(nodeGroup.find(*it) != nodeGroup.end()) { 
        masterNode = *it;
        break; 
      }
    }

    // Erase the involved elements from the original set;
    // TODO Check if we should keep track of the existence of the element
    for(std::set<int>::iterator it = elementGroup.begin(); it != elementGroup.end(); ++it)
      elem[*it] = 0;

    // Now create the new elements
    int *nodes = new (ba) int[nodeGroup.size()+1];
    nodes[1] = masterNode;
    for(std::set<int>::iterator it = nodeGroup.begin(); it != nodeGroup.end(); ++it) {
      if(*it != nodes[1]) {
        nodes[0] = *it; // slave
        newSet[iEl++] = new (ba) RigidBeam(nodes);
      }
    }
  }
  // Finish the new set
  for(int i = 0; i < emax; ++i)
    if(elem[i] != 0)
      newSet[iEl++] = elem[i];

  delete [] elem;
  elem = newSet;
  emax = iEl;
}

#ifdef SALINAS
#include <Element.d/CtcVirtualElt.d/CtcVirtualElt.h>
void
Elemset::elemadd(int num, int etype, int nnodes, int*n)
{
   Element *ele;

   switch(etype)
   {
     case 64: // CKT virtual Elements for contact =64
       ele = new (ba) CtcVirtualElt(nnodes,n);
       break;
     default:
       std::cerr << "Element Type " << etype << " is not Supported." << std::endl;
       exit(-1);
       return;
   }

 elemadd(num, ele);

}
#endif

