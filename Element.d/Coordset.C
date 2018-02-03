#include <cstdio>
#include <cstdlib>
#include <Driver.d/GeoSource.h>
#include <Element.d/Element.h>

extern GeoSource * geoSource;

typedef Node *NodeP;

CoordSet::CoordSet(int initialSize) : ba(4096, initialSize)
{
  last = 0;
  nmax  = initialSize;
  if(nmax > 0) {
    nodes = new NodeP[nmax];
    for(int i = 0; i < nmax ; ++i) nodes[i] = 0;
  }
  else nodes = 0;
}

CoordSet::~CoordSet() 
{ 
  if (nodes) delete [] nodes;
}

CoordSet &
CoordSet::operator = (const CoordSet & other)
{
  if (this != &other) // protect against invalid self-assignment
    {
      // 1: allocate new memory and copy the elements
      Node ** new_nodes = new Node * [other.nmax];
      std::copy(other.nodes, other.nodes + other.nmax, new_nodes);

      // 2: deallocate old memory
      if(nodes) { delete [] nodes; }

      // 3: assign the new memory to the object
      nodes = new_nodes;
      nmax = other.nmax;
      last = other.last;
    }
  // by convention, always return *this
  return *this;
}

Node *& CoordSet::operator[] (int n)
{
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
 if(n >= last) last = n+1;
 return nodes[n];
}

void CoordSet::nodeadd(int n, double *xyz, int cp, int cd)
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
 if(n >= last) last = n+1;
 nodes[n] = new (ba) Node(xyz, cp, cd) ;
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
 if(n >= last) last = n+1;
 nodes[n] = new (ba) Node(node) ;
}

Node &
CoordSet::getNode(int i)
{
 // KHPXXX: check if the node actually exists before returning
 //         a reference to it!
 if(nodes[i]==0) {
   fprintf(stderr,"*** ERROR: Node %d does not exist! Exiting ...\n",i+1);
   exit(-1);
 }
 return *nodes[i];
}

const Node &
CoordSet::getNode(int i) const
{
    // KHPXXX: check if the node actually exists before returning
    //         a reference to it!
    if(nodes[i]==0) {
        fprintf(stderr,"*** ERROR: Node %d does not exist! Exiting ...\n",i+1);
        exit(-1);
    }
    return *nodes[i];
}

void
CoordSet::getCoordinates(const int *nn, int numNodes,
                         double *xx, double *yy, double *zz) const
{
 int iNode;
 for(iNode=0; iNode<numNodes; ++iNode) {
   xx[iNode] = getNode(nn[iNode]).x;
   yy[iNode] = getNode(nn[iNode]).y;
   zz[iNode] = getNode(nn[iNode]).z;
 }
}

int
CoordSet::size() const
{
  return last;
}

int CoordSet::nnz() const
{
  int ret = 0;
  for(int i = 0; i < nmax; ++i)
    if(nodes[i] != 0) ret++;
  return ret;
}

NFrameData*
CoordSet::dofFrame(int i) const
{
  // get a pointer to the degree-of-freedom frame for the specified i^th node
  // if one is defined and it differs from the basic frame.
  int cd;
  if(i < nmax && nodes[i] && (cd = nodes[i]->cd) > 0) return &(geoSource->getNFrames()[cd]);
  else return nullptr;
}

std::pair<Eigen::Matrix<double, 3, 1>, int> CoordSet::computeSums() const {
    using V3D = Eigen::Matrix<double, 3, 1>;
    V3D coordSum{0.0, 0.0, 0.0};
    int count = 0;
    for(int i = 0; i < nmax; ++i)
        if(nodes[i]) {
            ++count;
            coordSum += V3D{nodes[i]->x, nodes[i]->y, nodes[i]->z};
        }
    return {coordSum, count};
}
