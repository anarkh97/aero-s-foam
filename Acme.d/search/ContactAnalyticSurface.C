// $Id$

#include "ContactAnalyticSurface.h"

ContactAnalyticSurface::ContactAnalyticSurface( 
    ContactSearch::AnalyticSurface_Type Type, int Host_Index, int Entity_Key )
  : ContactTopologyEntity<Real>( 0, Host_Index, NULL, CT_ANALYTIC_SURFACE)
{
  type = Type;
  entity_key = Entity_Key;
  ProcArrayIndex(Host_Index);
  Global_ID(0,Host_Index);
}

ContactAnalyticSurface::ContactAnalyticSurface( 
    ContactSearch::AnalyticSurface_Type Type, int Host_Index)
  : ContactTopologyEntity<Real>( 0, Host_Index, NULL, CT_ANALYTIC_SURFACE)
{
  type = Type;
  entity_key = -1;
  ProcArrayIndex(Host_Index);
  Global_ID(0,Host_Index);
}

ContactAnalyticSurface::ContactAnalyticSurface( 
    ContactSearch::AnalyticSurface_Type Type, 
    ContactAnalyticSurface* surf)
  : ContactTopologyEntity<Real>( 0, surf->host_array_index, NULL, CT_ANALYTIC_SURFACE)
{
  type = Type;
  entity_key = surf->entity_key;
  
  BlockID(surf->block_id);
  HostArrayIndex(surf->host_array_index);
  ProcArrayIndex(surf->proc_array_index);
  Global_ID(surf->global_id.HiInt(),surf->global_id.LoInt());
}

ContactAnalyticSurface::~ContactAnalyticSurface() {}

void
ContactAnalyticSurface::Display(ContactParOStream& postream)
{
  postream<<"ContactEntity: "<<global_id<<"\n";
  postream<<"               entity type:       CT_ANALYTIC_SURFACE\n";
  postream<<"               surface type:      ";
  switch (type) {
  case ContactSearch::PLANE:
    postream<<"PLANE\n";
    break;
  case ContactSearch::SPHERE:
    postream<<"SPHERE\n";
    break;
  case ContactSearch::CYLINDER_INSIDE:
    postream<<"CYLINDER_INSIDE\n";
    break;
  case ContactSearch::CYLINDER_OUTSIDE:
    postream<<"CYLINDER_OUTSIDE\n";
    break;
  default:
    postream<<"\n";
    break;
  }
  postream<<"               surface_id:        "<<proc_array_index<<"\n";
}

