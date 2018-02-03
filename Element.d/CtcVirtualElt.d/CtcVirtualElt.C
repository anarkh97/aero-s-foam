#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
#include <Element.d/CtcVirtualElt.d/CtcVirtualElt.h>
//#include <Corotational.d/BarCorotator.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>
//#include <Utils.d/linkfc.h>

CtcVirtualElt::CtcVirtualElt(int nnodes,int* nodenums)
{
  if (nnodes!=2) {
    cout << "Only virtual Elt with two nodes are supported yet " << endl;
    exit(-1);
  }
        nn[0] = nodenums[0];
        nn[1] = nodenums[1];
}

Element *
CtcVirtualElt::clone()
{
	return new CtcVirtualElt(*this);
}

void
CtcVirtualElt::renum(int *table)
{
	nn[0] = table[nn[0]];
	nn[1] = table[nn[1]];
}

void
CtcVirtualElt::renum(EleRenumMap& table)
{
	nn[0] = table[nn[0]];
	nn[1] = table[nn[1]];
}

int
CtcVirtualElt::numNodes() const
{
        return 2;
}

int *
CtcVirtualElt::nodes(int *p) const
{
        if(p == 0) p = new int[2];
        p[0] = nn[0];
        p[1] = nn[1];
        return p;
}

int
CtcVirtualElt::numDofs() const
{
  return 1;
}

int *
CtcVirtualElt::dofs(DofSetArray &dsa, int *p)
{
  if(p == 0) p = new int[1];
  dsa.number(nn[0],DofSet::Contact, p  ); // the first node is the virtual one
  return p;
}

void
CtcVirtualElt::markDofs( DofSetArray &dsa )
{
  cout << "> Should not have internal ctc dofs " << endl;
  exit (-1);
  // Never called in dofset.C because the virtual elements have been removed from the packed elt set
  //  dsa[nn[0]].zero();
  //dsa[nn[0]].mark(DofSet::Contact);
    
 
}

int
CtcVirtualElt::getTopNumber()
{
 return 0;
}


FullSquareMatrix
CtcVirtualElt::massMatrix(CoordSet &cs, double *mel, int flg) // should never be called
  // The virtual elt has no structural property -> not in packedEltSet
  // Hence, it is in dsa,c_dsa ... but not in connectivity tab allDofs (built from packedEltSet) of Krr
{
  cout << "No mass Matrix for virtual elements" << endl;
  exit(-1);
  
  return *(new FullSquareMatrix); // just to avoid warning

}

FullSquareMatrix
CtcVirtualElt::stiffness(CoordSet &cs,double *k, int flg) // should never be called
{
  cout << "No stiffness Matrix for virtual elements " << endl;
  exit (-1);
 
  return *(new FullSquareMatrix); // just to avoid warning
}





//double
//CtcVirtualElt::getMass(CoordSet& cs)
//{
//  mass=0;
//  return mass;
//}

//void
//CtcVirtualElt::getGravityForce(CoordSet& cs, double *gravityAcceleration, Vector &gravityForce)
//{
//  double massPerNode = 0;
//  gravityForce[0] = 0;
//  gravityForce[1] = 0;
//  gravityForce[2] = 0;
//  gravityForce[3] = 0;
//  gravityForce[4] = 0;
//  gravityForce[5] = 0;
//}







