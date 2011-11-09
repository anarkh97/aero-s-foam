// $Id$

#ifndef _ContactFixedSizeAllocator_h_
#define _ContactFixedSizeAllocator_h_

#include "contact_assert.h"
#include "ContactErrors.h"
#include <vector>

// WARNING:  Make sure to use the sizeof() operator when giving the storage
//           bytes to the initialization methods.  This will avoid data
//           alignment problems, which often give a "Bus Error".
//
// NOTE:  The Delete_Frag() function does not set the pointer to Null.

using namespace std;


class ContactFixedSizeAllocator {
public:
  
  ContactFixedSizeAllocator();
  ContactFixedSizeAllocator(
    int storage_bytes,                // Forced to be at least sizeof(char*).
    int block_length,                 // At least 2.
    int initial_length = 0,           // At least 2 (defaults to block_length).
    const char* object_name = 0);     // Name is optional.
 ~ContactFixedSizeAllocator();
  
  // Purges existing memory and resets with new parameters.
  void  Resize(int storage_bytes, int block_length, int initial_length = 0);
  
  void  Set_Name(const char* name);   // attaches a text name to the allocator
  const char* Name() const { return object_name; }
  
  inline void* New_Frag();
  inline void  Delete_Frag(void*);
  
  void  Purge();     // Frees all memory to the std::system (warns if frags exist).
  
private:
  char* object_name;   // Any text.  Used for identifying the allocator.
  
  // Definition variables:
  
  int storage_bytes;   // Bytes of memory in each frag.
  int block_size;      // Number of bytes per allocated block.
  int initial_size;    // Number of bytes allocated for very first block.
  int num_obj;         // Number of objects per block
  
  // State variables:
  
  std::vector<char*> mem_blocks;
  std::vector<void*> free_spots;

  unsigned int num_free;

  // Internal methods:
  
  int   Sanity_Check() const;
  
  // not allowed
  ContactFixedSizeAllocator(const ContactFixedSizeAllocator& src)
    { PRECONDITION(0); }
  const ContactFixedSizeAllocator& operator=(const ContactFixedSizeAllocator&)
    { PRECONDITION(0); return *this; }
};

void ContactFixedSizeAllocator::Delete_Frag(void* p)
{
#ifndef CONTACTALLOCATOR_BYPASS

  if (p)  // Don't want to try to delete a null pointer.
  {                      
    if(free_spots.size() <= num_free) free_spots.resize(num_free + num_obj);
    free_spots[num_free++] = p;
  }

#else

  char* tmp = (char*)p;
  if (tmp) delete [] tmp;

#endif
}

void* ContactFixedSizeAllocator::New_Frag()
{
// Compiling with this defined will bypass the fancy memory allocation stuff.
// Instead, regular new/delete operations will be performed.
#ifndef CONTACTALLOCATOR_BYPASS

  if(num_free == 0) {
    char *new_block = new char[block_size];
    mem_blocks.push_back(new_block);
    num_free = num_obj;
    if(free_spots.size() < num_free) free_spots.resize(num_free);
    for(int i = 0; i < num_obj; ++i) {
      free_spots[i] = &(new_block[i*storage_bytes]);
    }  
  }

  return free_spots[--num_free];

#else
  
  return (void*)new char[storage_bytes];
  
#endif
}





#endif
