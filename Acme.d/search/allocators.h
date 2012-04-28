// $Id$

class ContactFixedSizeAllocator;

template<typename DataType>
  void ContactNode_SizeAllocator(ContactFixedSizeAllocator& alloc);
template<typename DataType>
void ContactLineEdgeL2_SizeAllocator(ContactFixedSizeAllocator& alloc);
void ContactLineEdgeQ3_SizeAllocator(ContactFixedSizeAllocator& alloc);
void ContactLineFaceL2_SizeAllocator(ContactFixedSizeAllocator& alloc);
void ContactLineFaceQ3_SizeAllocator(ContactFixedSizeAllocator& alloc);
template<typename DataType>
  void ContactQuadFaceL4_SizeAllocator(ContactFixedSizeAllocator& alloc);
void ContactQuadFaceQ8_SizeAllocator(ContactFixedSizeAllocator& alloc);
void ContactQuadFaceQ9_SizeAllocator(ContactFixedSizeAllocator& alloc);
template<typename DataType>
  void ContactTriFaceL3_SizeAllocator(ContactFixedSizeAllocator& alloc);
void ContactTriFaceQ6_SizeAllocator(ContactFixedSizeAllocator& alloc);
template<typename DataType>
  void ContactHexElemL8_SizeAllocator(ContactFixedSizeAllocator& alloc);
template<typename DataType>
  void ContactWedgeElemL6_SizeAllocator(ContactFixedSizeAllocator& alloc);
void ContactNodeNodeInteraction_SizeAllocator(ContactFixedSizeAllocator& alloc);
void ContactNodeFaceInteraction_SizeAllocator(ContactFixedSizeAllocator& alloc);
void ContactNodeSurfaceInteraction_SizeAllocator(ContactFixedSizeAllocator& alloc);
void ContactFaceFaceInteraction_SizeAllocator(ContactFixedSizeAllocator& alloc);
void ContactFaceCoverageInteraction_SizeAllocator(ContactFixedSizeAllocator& alloc);
void ContactElementElementInteraction_SizeAllocator(ContactFixedSizeAllocator& alloc);
void ContactDLLTopologyNode_SizeAllocator(ContactFixedSizeAllocator& alloc);
void ContactDLLInteractionNode_SizeAllocator(ContactFixedSizeAllocator& alloc);
void ContactLLNode_SizeAllocator(ContactFixedSizeAllocator& alloc);
void ContactLLTopologyNode_SizeAllocator(ContactFixedSizeAllocator& alloc);
void ContactLLNodeData_SizeAllocator(ContactFixedSizeAllocator& alloc);
void ContactPolyVert_SizeAllocator(ContactFixedSizeAllocator& alloc);
void ContactPolyEdge_SizeAllocator(ContactFixedSizeAllocator& alloc);
void ContactPoly_SizeAllocator(ContactFixedSizeAllocator& alloc);
void ContactFaceFaceGraphNode_SizeAllocator(ContactFixedSizeAllocator& alloc);
void ContactCartesianHexElementL8_SizeAllocator(ContactFixedSizeAllocator& alloc);
void ContactHexElementL8_SizeAllocator(ContactFixedSizeAllocator& alloc);
template<typename DataType>
  void ContactShellQuadFaceL4_SizeAllocator(ContactFixedSizeAllocator& alloc);
template<typename DataType>
  void ContactShellTriFaceL3_SizeAllocator(ContactFixedSizeAllocator& alloc);
void ContactShellNode_SizeAllocator(ContactFixedSizeAllocator& alloc);
