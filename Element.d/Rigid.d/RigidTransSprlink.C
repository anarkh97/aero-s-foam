#include <Element.d/Rigid.d/RigidTransSprlink.h>
#include <Element.d/Joint.d/LinearConstraintType1.h>

RigidTransSprlink::RigidTransSprlink(int* _nn)
 : SuperElement(true)
{ 
  nnodes = 2;
  nn = new int[nnodes];
  for(int i = 0; i < nnodes; ++i) nn[i] = _nn[i];
}

void
RigidTransSprlink::setProp(StructProp* _prop, bool _myProp)
{
  nSubElems = int(_prop->kx != 0.0) + int(_prop->ky != 0.0) + int(_prop->kz != 0.0);
  subElems = new Element * [nSubElems];
  int indices[2] = { 0, 1 };
  int count = 0;
  if(_prop->kx != 0.0) {
    DofSet xx[2] = { DofSet::Xdisp,  DofSet::Xdisp };
    subElems[count++] = new LinearConstraintType1(indices, xx);
  }
  if(_prop->ky != 0.0) {
    DofSet yy[2] = { DofSet::Ydisp,  DofSet::Ydisp };
    subElems[count++] = new LinearConstraintType1(indices, yy);
  }
  if(_prop->kz != 0.0) {
    DofSet zz[2] = { DofSet::Zdisp,  DofSet::Zdisp };
    subElems[count++] = new LinearConstraintType1(indices, zz);
  }
  for(int i = 0; i < nSubElems; ++i) subElems[i]->buildFrame(*css); // since these elements had not been instantiated when
                                                                    // SupereElement::buildFrame was called
  SuperElement::setProp(_prop, _myProp);
}
