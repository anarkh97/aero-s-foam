#include <Element.d/Element.h>
#include <iostream>
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
	  int num = elem[i]->totNumNodes();
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

