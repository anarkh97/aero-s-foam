#include <Utils.d/BlockAlloc.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>

void *
BlockAlloc::getMem(size_t nbyte)
{
  if(int(nbyte) > blLen) return new char[nbyte];
  if(block == 0 || (index+int(nbyte)) > blLen) {
     block = new char[blLen];
     allBlocks[nblock++] = block;
     index = 0;
  }
  if(nbyte & 0x7) {
    nbyte = (nbyte+8)-(nbyte & 0x7);
  }
  void *p = block+index;
  index += nbyte;
  return p;
}

BlockAlloc::~BlockAlloc()
{
 int iBlock;
 for(iBlock = 0; iBlock < nblock; ++iBlock)
   delete [] allBlocks[iBlock]; 
}

void * operator new(size_t nbyte, BlockAlloc &block)
{
 return block.getMem(nbyte);
}

void * operator new[](size_t nbyte, BlockAlloc &block)
{
 return block.getMem(nbyte);
}

BlockAlloc global_ba;

void operator delete(void *p, BlockAlloc &block)
{
  // this is only used for exception unwinding
  free(p);
}
