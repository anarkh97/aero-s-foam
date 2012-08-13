// This file is part of a modified version of ACME: Algorithms for Contact in
// a Multiphysics Environment, derived from version 2.7f
//
// Copyright (C) 2007 Sandia Corporation
// Copyright (C) 2011 Stanford University 
//
// ACME is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ACME is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with ACME.  If not, see <http://www.gnu.org/licenses/>.


#ifndef ContactDoublyLinkedList_h_
#define ContactDoublyLinkedList_h_

#include <vector>
#include "ContactTopologyEntity.h"
#include "ContactInteractionEntity.h"
#include "contact_assert.h"
#include "Contact_Defines.h"

class ContactTopologyDLL {

 public:
  ContactTopologyDLL();
  ~ContactTopologyDLL();
  int Append( ContactTopologyEntity<Real>* );
  
  inline void IteratorStart() {
    current = data.begin(); 
    while(current < data.end()) {
      if(*current != NULL) return;
      ++current;
    }
  };
  inline ContactTopologyEntity<Real>* IteratorForward() { 
    while(current < data.end()) {
      ContactTopologyEntity<Real> *cur_object = *current;
      current++;
      if(cur_object != NULL) return cur_object;
    }
    return NULL;
  };

  void SetEntity(const int pos, ContactTopologyEntity<Real>* entity) {
    data[pos] = entity;
  }

  void RemoveEntity(const int pos) {
    data[pos] = NULL;
    num_entities--;
  }

  inline int NumEntities() {return num_entities;};
  
  void Display(ContactParOStream&);
  void Clear();
  
  std::vector< ContactTopologyEntity<Real>* > *EntityList() { return &data; };

 private:
  std::vector< ContactTopologyEntity<Real>* > data;
  std::vector< ContactTopologyEntity<Real>* >::iterator current;
  int num_entities;
};

//=================================================================================================================

class ContactInteractionDLL {

 public:
  ContactInteractionDLL();
  ~ContactInteractionDLL();
  void Append( ContactInteractionEntity* );
  
  void DeletePrev() {
    std::vector< ContactInteractionEntity* >::iterator data_to_erase = current; 
    data_to_erase--;   
    data.erase(data_to_erase);
  }
  
  inline void IteratorStart() {current = data.begin();};

  inline ContactInteractionEntity* IteratorForward() { 
    if(current < data.end()) {
      return *(current++);
    } else {
      return NULL;
    }
  };

  inline ContactInteractionEntity* HeadEntity() {return data[0];};

  inline ContactInteractionEntity* CurrentEntity() {return *current;};

  inline int NumEntities() {return data.size();};
  
  void Display(ContactParOStream&);
  void Reset();
  void Clear();

 private:
  std::vector< ContactInteractionEntity* > data;
  std::vector< ContactInteractionEntity* >::iterator current;
};

#endif
