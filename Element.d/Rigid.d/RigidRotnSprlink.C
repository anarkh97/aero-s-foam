#include <Element.d/Rigid.d/RigidRotnSprlink.h>
#include <Element.d/Joint.d/LinearConstraintType1.h>

RigidRotnSprlink::RigidRotnSprlink(int* _nn)
{
  nnodes = 2;
  nn = new int[nnodes];
  for(int i = 0; i < nnodes; ++i) nn[i] = _nn[i];
}

void
RigidRotnSprlink::setProp(StructProp* _prop, bool _myProp)
{
  nSubElems = int(_prop->kx != 0.0) + int(_prop->ky != 0.0) + int(_prop->kz != 0.0);
  subElems = new Element * [nSubElems];
  int indices[2] = { 0, 1 };
  int count = 0;
  if(_prop->kx != 0.0) {
    DofSet xx[2] = { DofSet::Xrot,  DofSet::Xrot };
    subElems[count++] = new LinearConstraintType1(indices, xx);
  }
  if(_prop->ky != 0.0) {
    DofSet yy[2] = { DofSet::Yrot,  DofSet::Yrot };
    subElems[count++] = new LinearConstraintType1(indices, yy);
  }
  if(_prop->kz != 0.0) {
    DofSet zz[2] = { DofSet::Zrot,  DofSet::Zrot };
    subElems[count++] = new LinearConstraintType1(indices, zz);
  }
  initialize(nnodes, nn);
  SuperElement::setProp(_prop, _myProp);
}

