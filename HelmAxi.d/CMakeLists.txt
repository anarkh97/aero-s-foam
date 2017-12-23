ADD_LIBRARY(HelmAxi
coefFourier.C
FourierDescrip.C
FourierHelmBCs.C
FourierProbTyp.C
HelmAxiQuad8.C
HelmAxiQuad.C
HelmAxiTri6.C
HelmAxiTri.C
Line2AxiSommer.C
LineAxiSommer.C
MDAxiData.C
MDAxiDesc.C
MPCData.C
quad8aximas.f
quad8axistif1.f
quad8axistif2.f
quad8shape.f
quadaximas.f
quadaxistif1.f
quadaxistif2.f
ScatterData.C
StaticTimersFetiHAxi.C
tri3aximas.f
tri3axistif1.f
tri3axistif2.f
tri6aximas.f
tri6axistif1.f
tri6axistif2.f
trig3shape.f
trig6shape.f
triggaussq.f 
FetiHAxi.d/DistrComplexVector.C     FetiHAxi.d/GCRC.C    FetiHAxi.d/GMRESDCV.C
FetiHAxi.d/DistrComplexVectorSet.C  FetiHAxi.d/GCRDCV.C  FetiHAxi.d/MPCSolver.C
FetiHAxi.d/FetiHAxi.C               FetiHAxi.d/GMRESC.C)
TARGET_LINK_LIBRARIES(HelmAxi Eigen3::Eigen)

SET(INC_DIR ..)

INCLUDE_DIRECTORIES( ${INC_DIR} )
