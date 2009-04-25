#include <stdio.h>
#include <stdlib.h>
#include <Element.d/Element.h>

typedef Node *NodeP;

CoordSet::CoordSet(int initialSize) : ba(4096, initialSize)
{
  nmax  = initialSize;
  if(nmax > 0) {
    nodes = new NodeP[nmax];
    for(int i = 0; i < nmax ; ++i) nodes[i] = 0;
  }
  else nodes = 0;
}

CoordSet::CoordSet(const CoordSet &cs, int n, int *ndn)
{
 nmax  = n;
 if(nmax > 0) {
   nodes = new NodeP[nmax];
   for(int i = 0; i < nmax ; ++i) nodes[i] = cs.nodes[ndn[i]];
 }
 else nodes = 0;
}

void CoordSet::nodeadd(int n, double *xyz)
{
 if(n<0) { return ; }
 if(n >= nmax) // resize nodes[]
  {
    int newsize = ((n+1)*3)/2;
    NodeP *np = new NodeP[newsize] ;
    int i;
    for(i= 0; i < nmax ; ++i)
     np[i] = nodes[i] ;
    for(nmax = newsize; i < nmax; ++i)
     np[i] = 0 ;
    if(nodes) delete [] nodes ;
    nodes = np ;
  }
 nodes[n] = new (ba) Node(xyz) ;
}

void CoordSet::nodeadd(int n, Node &node)
{
 if(n<0) { return ; }
 if(n >= nmax) // resize nodes[]
  {
    int newsize = ((n+1)*3)/2;
    NodeP *np = new NodeP[newsize] ;
    int i;
    for(i= 0; i < nmax ; ++i)
     np[i] = nodes[i] ;
    for(nmax = newsize; i < nmax; ++i)
     np[i] = 0 ;
    if(nodes) delete [] nodes ;
    nodes = np ;
  }
 nodes[n] = new (ba) Node(node) ;
}

//------------------------------------------------------

void CoordSet::nodeCopy(int numNodes, CoordSet &nd)  {

  nmax = numNodes;
  if(nmax < 0) { 
    fprintf(stderr,"*** ERROR: Number of nodes is less than zero: %d\n",
			nmax);
    return; 
  }

  Node **np = new Node *[numNodes] ;
  int i;
  for(i= 0; i < numNodes; ++i)
    np[i] = nd[i] ;

  if(nodes) delete [] nodes ;
  nodes = np ;
}

//------------------------------------------------------

Node &
CoordSet::getNode(int i)
{
 // KHPXXX: check if the node actually exists before returning
 //         a reference to it!
 //if(nodes[i]==0) {
 //  fprintf(stderr,"*** ERROR: Node %d does not Exist! Exiting\n",i+1);
 //  exit(-1);
 //}
 return *nodes[i];
}

void
CoordSet::getCoordinates(int *nn, int numNodes, 
                         double *xx, double *yy, double *zz)
{
 int iNode;
 for(iNode=0; iNode<numNodes; ++iNode) {
   xx[iNode] = getNode(nn[iNode]).x;
   yy[iNode] = getNode(nn[iNode]).y;
   zz[iNode] = getNode(nn[iNode]).z;
 }
}

CoordSet * CoordSet::copy() {

    CoordSet * newset = new CoordSet;
  
    double * xyz=new double[3];

    int i;
    for(i= 0; i < nmax ; ++i) {
      if (nodes[i]) {
        xyz[0]=nodes[i]->x;
        xyz[1]=nodes[i]->y;
        xyz[2]=nodes[i]->z;
	newset->nodeadd(i,xyz);
      }
    }
    
    return newset;
}    

// Last returns the true last defined node
int CoordSet::last()
{
 int last = size();
 while(--last >= 0)
    if(nodes[last] != 0) break;

 return last+1;
}

int CoordSet::exist(int nodeNumber)
{
 // check if node is larger than the number of nodes
 // stored in the CoordSet
 // if(nodeNumber > size()) return 0;
 if(nodeNumber >= size()) return 0; // PJSA 4-13-05

 // otherwise, check the pointer to the node in the
 // CoordSet to see if the node is used or not.
 if(nodes[nodeNumber]) 
   return 1;
 else
   return 0;
}

