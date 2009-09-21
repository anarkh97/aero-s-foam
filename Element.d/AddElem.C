 #include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include <Utils.d/DistHelper.h> //HB

#include <Element.d/Element.h>
#include <Element.d/Truss.d/TwoNodeTruss.h>
#include <Element.d/Truss.d/TwoNodeTrussF.h>
#include <Element.d/Beam.d/EulerBeam.h>
#include <Element.d/Beam.d/TimoshenkoBeam.h>
#include <Element.d/Shell.d/ThreeNodeShell.h>
#include <Element.d/Shell.d/FourNodeShell.h>
#include <Element.d/Shell.d/Therm3NoShell.h>
#include <Element.d/Shell.d/Therm4NoShell.h>
#include <Element.d/Quad4.d/FourNodeQuad.h>
#include <Element.d/Quad4.d/Quad.h>
#include <Element.d/Brick.d/EightNodeBrick.h>
#include <Element.d/Tetra.d/Tetrahedral.h>
#include <Element.d/Tetra10.d/TenNodeTetrahedral.h>
#include <Element.d/Penta.d/Pentahedral.h>
#include <Element.d/Membrane.d/Membrane.h>
#include <Element.d/Spring.d/TorSpring.h>
#include <Element.d/Spring.d/TransSprlink.h>
#include <Element.d/Spring.d/RotnSprlink.h>
#include <Element.d/CompShell.d/Compo3NodeShell.h>
#include <Element.d/CompShell.d/Compo4NodeShell.h>
#include <Element.d/Triangle3.d/Triangle3.h>
#include <Element.d/Triangle3.d/ThermTriangle.h>
#include <Element.d/ThermQuad.d/ThermQuadGal.h>
#include <Element.d/ThermQuad.d/Therm3DQuad.h>
#include <Element.d/Brick.d/ThermBrick.h>
#include <Element.d/Tetra.d/ThermIsoParamTetra.h>
#include <Element.d/Truss.d/Therm2NodeBar.h>
#include <Element.d/Helm.d/HelmQuadGal.h>
#include <Element.d/Helm.d/Tetra10HelmGal.h>
#include <Element.d/Helm.d/HelmQuadGls.h>
#include <Element.d/Helm.d/TetraHelmGal.h>
#include <Element.d/Helm.d/TetraHelmGLS.h>
#include <Element.d/Helm.d/HelmTri3Gal.h>
#include <Element.d/Helm.d/HelmTri3Gls.h>
#include <Element.d/Helm.d/HelmBrick.h>
#include <Element.d/Helm.d/HelmTri6Gal.h>
#include <Element.d/Helm.d/HelmQuad8Gal.h>
#include <Element.d/Helm.d/HelmLagQuadGal.h>
#include <Element.d/Helm.d/Tetra10HelmGal.h>
#include <Element.d/Helm.d/HelmBrickGLS.h>
#include <Element.d/Helm.d/HelmPenta.h> //HB

#include <Element.d/Helm.d/HelmIsoParamHexa.h>
#include <Element.d/Helm.d/HelmSpectralIsoParamHexa.h>
#include <Element.d/Helm.d/HelmIsoParamTetra.h>
#include <Element.d/Helm.d/HelmIsoParamQuad.h>
#include <Element.d/Helm.d/HelmSpectralIsoParamQuad.h>
#include <Element.d/Helm.d/HelmIsoParamTri.h>
#include <Element.d/Helm.d/LEIsoParamQuad.h>
#include <Element.d/Helm.d/LEIsoParamTri.h>
#include <Element.d/Helm.d/LEIsoParamHexa.h>
#include <Element.d/Helm.d/LEIsoParamTetra.h>

#include <Element.d/Shear.d/ShearPanel.h>
#include <Element.d/Convection.d/BarConvec.h>
#include <Element.d/Convection.d/QuadConvec.h>
#include <Element.d/Convection.d/TriangleConvec.h>
#include <Element.d/Radiation.d/BarRadiation.h>
#include <Element.d/Radiation.d/TriangleRadiation.h>
#include <Element.d/Radiation.d/QuadRadiation.h>
#include <Element.d/ContactResistance.d/QuadContact.h>
#include <Element.d/ContactResistance.d/BrickContact.h>
#include <Element.d/ContactResistance.d/PentaContact.h>
#include <Element.d/BulkFluid.d/TriangleBulk.h>
#include <Element.d/BulkFluid.d/TetraBulk.h>
#include <Element.d/BulkFluid.d/PentaBulk.h>
#include <Element.d/CtcVirtualElt.d/CtcVirtualElt.h>
#include <Element.d/Truss.d/TwoNodeTrussRigid.h>
#include <Element.d/Beam.d/RigidBeam.h>
#include <Element.d/Spring.d/RigidSpring.h>
#include <Element.d/Spring.d/RigidSpringTr.h>
#include <Element.d/Spring.d/RigidSpringRo.h>
#include <Element.d/Brick.d/EightNodeBrickRigid.h>
#include <Element.d/Solid.d/RigidSolid.h>
#include <Element.d/Solid.d/RigidSolid6Dof.h>
#include <Element.d/Brick20.d/Brick20.h>
#include <Element.d/Shell.d/RigidThreeNodeShell.h>
#include <Element.d/NonLinearity.d/NLHexahedral.h>
#include <Element.d/NonLinearity.d/NLMembrane.h>
#include <Element.d/Shell.d/ConnectedTri.h>

#include <HelmAxi.d/HelmAxiQuad.h>
#include <HelmAxi.d/HelmAxiTri.h>
#include <HelmAxi.d/HelmAxiQuad8.h>
#include <HelmAxi.d/HelmAxiTri6.h>

#include <Element.d/MpcElement.d/MpcElement.h>
#include <Element.d/MpcElement.d/FsiElement.h>

#include <map>
extern map<int,double > weightList;

#include <Element.d/Truss.d/TwoNodeTrussRigidMpc.h>
#include <Element.d/Beam.d/RigidMpcBeam.h>
#include <Element.d/Spring.d/RigidMpcSpring.h>
#include <Element.d/Spring.d/RigidMpcSpringTr.h>
#include <Element.d/Spring.d/RigidMpcSpringRo.h>
#include <Element.d/Solid.d/RigidMpcSolid6Dof.h>
#include <Element.d/Solid.d/RigidMpcSolid.h>
#include <Element.d/Rigid.d/RBE2Mpc.h>
#include <Element.d/Rigid.d/RBE2.h>

#include <Element.d/Brick32.d/Brick32.h> 
#include <Element.d/Penta26.d/Penta26.h> 
#include <Element.d/Helm.d/HelmBrick32.h> 
#include <Element.d/Helm.d/HelmPenta26.h> 
#include <Element.d/Penta15.d/Penta15.h> 

#include <Element.d/DEM.d/DEMHelm2d.h>
#include <Element.d/DEM.d/DEMHelm3d.h>
#include <Element.d/DEM.d/DEMLE2d.h>
#include <Element.d/DEM.d/DEMLE3d.h>

//ADDED FOR SLOSHING PROBLEM, EC, 20070713
#include <Element.d/FluidQuad.d/SloshQuadGal.h> 
#include <Element.d/FluidQuad.d/BarSloshFS.h>
#include <Element.d/FluidTetra.d/SloshTetra.h> 
#include <Element.d/FluidTriangle3.d/SloshTriangleFS.h> 

//ADDED FOR HEV PROBLEM, EC, 20070815
#include <Element.d/FluidQuad.d/HEVibQuadGal.h> 
#include <Element.d/FluidTetra.d/HEVibTetra.h> 

#include <Driver.d/Domain.h>
extern Domain *domain;
extern std::auto_ptr<ElementFactory> elemFact;

void
Elemset::elemadd(int num, int etype, int nnodes, int*n)
{
  Element *ele = elemFact->elemadd(num, etype, nnodes, n, ba);
  elemadd(num, ele);
  
  map<int, double >::iterator it = weightList.find(etype);
  if(it == weightList.end())
    {
      ele->setWeight(1.0);
      ele->setTrueWeight(1.0);
    }
  else
    {
      ele->setWeight(it->second);
      ele->setTrueWeight(it->second);
    }
  return;
}

Element*
ElementFactory::elemadd(int num, int etype, int nnodes, int*n, BlockAlloc& ba)
{
   Element *ele;
   bool grbmeig = (domain->solInfo().probType == SolverInfo::Modal && domain->solInfo().rbmflg == 1);
   bool rigidmpc = ((domain->solInfo().type == 2) || geoSource->getDirectMPC()
                   || (grbmeig && (etype != 66) && (etype != 73) && (etype != 74)) // safe to leave these as rigid for GRBM
                   || (domain->solInfo().isNonLin() && etype != 65));
   switch(etype) 
   {
     case 1:
       ele = new (ba) TwoNodeTruss(n);
       break;
     case 2:
       ele = new (ba) FourNodeQuad(n);
       // KHP: new Quad element that will work
       //      for 4, 8, 12 nodes per element
       //ele = new (ba) Quad(4, n);
       break;
     case 3:
       ele = new (ba) Therm3DQuad(n);
       break;
     case 4:
       ele = new (ba) Triangle3(n);
       break;
     case 6: {
       int nn[3] = { n[0], n[1], -1 };
       if(nnodes > 2) nn[2] = n[2];
       ele = new (ba) EulerBeam(nn);
       break;
       }
     case 7: {
       int nn[3] = { n[0], n[1], -1 };
       if(nnodes > 2) nn[2] = n[2];
       ele = new (ba) TimoshenkoBeam(nn);
       break;
       }
     case 8:
       ele = new (ba) ThreeNodeShell(n);
       break;
     case 9:
       ele = new (ba) Therm2NodeBar(n);
       break;
     case 10:
       ele = new (ba) ThermQuadGal(n);
       break;
     case 11:
       ele = new (ba) TorSpring(n);
       break;
     case 17:
       ele = new (ba) EightNodeBrick(n);
       break;
     case 18:
       ele = new (ba) ShearPanel(n);
       break;
     case 19:
       ele = new (ba) Membrane(n);
       break;
     case 20:
       ele = new (ba) Compo3NodeShell(n);
       break;
     case 21:
       ele = new (ba) TransSprlink(n);
       break;
     case 22:
       ele = new (ba) RotnSprlink(n);
       break;
     case 23:
       ele = new (ba) Tetrahedral(n);
       break;
     case 24:
       ele = new (ba) Pentahedral(n);
       break;
     case 25:
       ele = new (ba) TenNodeTetrahedral(n);
       break;
     case 30:
       ele = new (ba) HelmQuadGal(n);
       break;
     case 31:
       ele = new (ba) HelmQuadGls(n);
       break;
     case 32:
       ele = new (ba) HelmQuad8Gal(n);
       break;
     case 35:
       ele = new (ba) HelmTri3Gal(n);
       break;
     case 36:
       ele = new (ba) HelmTri3Gls(n);
       break;
     case 38:
       ele = new (ba) HelmTri6Gal(n);
       break;
     case 40:
       ele = new (ba) TetraHelmGal(n);
       break;
     case 41:
       ele = new (ba) TetraHelmGLS(n);
       break;
     case 42:
       ele = new (ba) Tetra10HelmGal(n);
	 break;
     case 43:
       ele = new (ba) HelmLagQuadGal(nnodes,n);
       break;
     case 44:
       ele = new (ba) HelmBrickGLS(n);
       break;
     case 45:
       ele = new (ba) HelmBrick(n);
       break;
     case 46:
       ele = new (ba) Therm3NoShell(n);
       break;
     case 47:
       ele = new (ba) BarConvec(n);
       break;
     case 48:
       ele = new (ba) QuadConvec(n);
       break;
     case 49:
       ele = new (ba) TriangleConvec(n);
       break;
     case 50:
       ele = new (ba) ThermIsoParamTetra(nnodes,n);
       break;     
     case 51:
       ele = new (ba) ThermBrick(n);
       break;
     case 53:
       ele = new (ba) ThermTriangle(n);
       break;
     case 56:
       ele = new (ba) BarRadiation(n);
       break;
     case 57:
       ele = new (ba) TriangleRadiation(n);
       break;
     case 58:
       ele = new (ba) QuadRadiation(n);
       break;
     case 59:
       ele = new (ba) HelmAxiTri6(n);
       break;
     case 60:
       ele = new (ba) HelmAxiQuad(n);
       break;
     case 61:
       ele = new (ba) HelmAxiTri(n);
       break;
     case 62:
       ele = new (ba) HelmAxiQuad8(n);
       break;
     case 63:
       ele = new (ba) HelmLagQuadGal(nnodes,n);
       break;
     case 64: // CKT virtual Elements for contact =64
       ele = new (ba) CtcVirtualElt(nnodes,n);
       break;
     case 65:
       if(rigidmpc) ele = new (ba) TwoNodeTrussRigidMpc(n); // FETI
       else ele = new (ba) TwoNodeTrussRigid(n); // direct
       break;
     case 66:
       if(rigidmpc) ele = new (ba) RigidMpcBeam(n);
       else ele = new (ba) RigidBeam(n);
       break;
     case 67:
       if(rigidmpc) ele = new (ba) RigidMpcSpring(n);
       else ele = new (ba) RigidSpring(n);
       break;
     case 68:
       if(rigidmpc) ele = new (ba) RigidMpcSpringTr(n);
       else ele = new (ba) RigidSpringTr(n);
       break;
     case 69:
       if(rigidmpc) ele = new (ba) RigidMpcSpringRo(n);
       else ele = new (ba) RigidSpringRo(n);
       break;
     case 70:
       if(rigidmpc) ele = new (ba) RigidMpcSolid(8,n);
       else ele = new (ba) EightNodeBrickRigid(n);
       break;
     case 71:
       if(rigidmpc) ele = new (ba) RigidMpcSolid(nnodes,n);
       else ele = new (ba) RigidSolid(nnodes,n);
       break;
     case 72:
       ele = new (ba) Brick20(n);
       break;
     case 73:
       if(rigidmpc) ele = new (ba) RigidMpcSolid6Dof(3,n);
       else ele = new (ba) RigidThreeNodeShell(n);
       break;
     case 74:
       if(rigidmpc) ele = new (ba) RigidMpcSolid6Dof(nnodes,n);
       else ele = new (ba) RigidSolid6Dof(nnodes,n);
       break;
     case 75:  // PJSA 3-30-05
       if(rigidmpc) ele = new (ba) RBE2Mpc(nnodes,n);
       else ele = new (ba) RBE2(nnodes,n);
       break;
     case 80:
       ele = new (ba) ConnectedTri(n);
       break;
     case 81:
       ele = new (ba) QuadContact(n);  // Interface element
       break;
     case 82:
       ele = new (ba) BrickContact(n);  // Interface element
       break;
     case 83:
       ele = new (ba) PentaContact(n);  // Interface element
       break;
     case 84:
       ele = new (ba) TriangleBulk(n);  // 2D Bulk Fluid Element
       break;
     case 85:
       ele = new (ba) TetraBulk(n);  // 3D Bulk Fluid Element
       break;
     case 86:
       ele = new (ba) PentaBulk(n);  // 3D Bulk Fluid Element
       break;
     case 88:
       {
         // PJSA: first count number of unique nodes & if there are only 3 make a tri
         int i, j;
         int count = 0;
         bool isUnique[4] = { true, true, true, true };
         for(i=0; i<4; ++i)
           for(j=0; j<i; ++j)
             if((n[i] == n[j]) && isUnique[j]) { count++; isUnique[i] = false; }
         if(count > 1) {
           filePrint(stderr," *** ERROR: FourNodeShell type 88, number %d, has %d duplicate nodes \n",num+1,count);
           exit(-1);
         }
         else if(count == 1) {
           int new_n[3];
           j=0;
           for(i=0; i<4; ++i) if(isUnique[i]) new_n[j++] = n[i];
           //filePrint(stderr," *** WARNING: converting degenerate FourNodeShell to ThreeNodeShell, ele = %d,  nodes: %6d %6d %6d\n",
           //         num+1,new_n[0]+1,new_n[1]+1,new_n[2]+1);
           ele = new (ba) ThreeNodeShell(new_n);
         }
         else {
           ele = new (ba) FourNodeShell(n); // superelement comprising two ThreeNodeShells
         }
       }
       break;
     case 90:
       ele = new (ba) HelmPenta(n); 
       break;
     case 91:
       //filePrint(stderr," *** add a Brick32 element\n");
       ele = new (ba) Brick32(n); 
       break;
     case 92:
       //filePrint(stderr," *** add a Penta26 element\n");
       ele = new (ba) Penta26(n); 
       break;
     case 93:
       //filePrint(stderr," *** add a HelmBrick32 element\n");
       ele = new (ba) HelmBrick32(n); 
       break;
     case 94:
       //filePrint(stderr," *** add a HelmPenta26 element\n");
       ele = new (ba) HelmPenta26(n); 
       break;
     case 95:
       ele = new (ba) HelmIsoParamHexa(nnodes,n);
       break;
     case 96:
       ele = new (ba) HelmIsoParamTetra(nnodes,n);
       break;
     case 97:
       //filePrint(stderr," *** add a Penta15 element\n");
       ele = new (ba) Penta15(n);
       break;
     case 98:
       ele = new (ba) HelmIsoParamQuad(nnodes,n);
       break;
     case 99:
       ele = new (ba) HelmIsoParamTri(nnodes,n);
       break;
     case 100:
       ele = new (ba) LEIsoParamQuad(nnodes,n);
       break;
     case 101:
       ele = new (ba) LEIsoParamTri(nnodes,n);
       break;
     case 102:
       ele = new (ba) LEIsoParamHexa(nnodes,n);
       break;
     case 103:
       ele = new (ba) LEIsoParamTetra(nnodes,n);
       break;
     case 105:
       ele = new (ba) HelmSpectralIsoParamHexa(nnodes,n);
       break;
     case 108:
       ele = new (ba) HelmSpectralIsoParamQuad(nnodes,n);
       break;
     case 111:
       ele = new (ba) TwoNodeTrussF(n);
       break;
     case 201:
       ele = new (ba) NLHexahedral(n,true); // linear kinematics
       break;
     case 202:
       ele = new (ba) NLHexahedral(n,false); // non-linear kinematics
       break;
     case 203:
       ele = new (ba) NLMembrane(n,false);
       break;
     case 204:
       ele = new (ba) NLMembrane(n, true);
       break;
     case 2020:
       ele = new (ba) Compo4NodeShell(n);
       break;
     case 4646:
       ele = new (ba) Therm4NoShell(n);
       break; 

     case 1100:
       ele = new (ba) DGMHelm2d_4(nnodes,n);
       break;
     case 1110:
       ele = new (ba) DGMHelm2d_4t(nnodes,n);
       break;
     case 1101:
       ele = new (ba) DGMHelm2d_8(nnodes,n);
       break;
     case 1111:
       ele = new (ba) DGMHelm2d_8t(nnodes,n);
       break;
     case 1102:
       ele = new (ba) DGMHelm2d_16(nnodes,n);
       break;
     case 1103:
       ele = new (ba) DGMHelm2d_32(nnodes,n);
       break;
     case 1104:
       ele = new (ba) DGMHelm2d_Eva2_8(nnodes,n);
       break;
     case 1120:
       ele = new (ba) DEMHelm2d_4(nnodes,n);
       break;
     case 1130:
       ele = new (ba) DEMHelm2d_4t(nnodes,n);
       break;
     case 1121:
       ele = new (ba) DEMHelm2d_8(nnodes,n);
       break;
     case 1131:
       ele = new (ba) DEMHelm2d_8t(nnodes,n);
       break;
     case 1122:
       ele = new (ba) DEMHelm2d_16(nnodes,n);
       break;
     case 1123:
       ele = new (ba) DEMHelm2d_32(nnodes,n);
       break;

     case 1200:
       ele = new (ba) DGMLE2d_4(nnodes,n);
       break;
     case 1201:
       ele = new (ba) DGMLE2d_16(nnodes,n);
       break;
     case 1220:
       ele = new (ba) DEMLE2d_4(nnodes,n);
       break;

     case 1150:
       ele = new (ba) DGMHelm3d_6(nnodes,n);
       break;
     case 1160:
       ele = new (ba) DGMHelm3d_6t(nnodes,n);
       break;
     case 1151:
       ele = new (ba) DGMHelm3d_26(nnodes,n);
       break;
     case 1161:
       ele = new (ba) DGMHelm3d_26t(nnodes,n);
       break;
     case 1152:
       ele = new (ba) DGMHelm3d_56(nnodes,n);
       break;
     case 1162:
       ele = new (ba) DGMHelm3d_56t(nnodes,n);
       break;
     case 1153:
       ele = new (ba) DGMHelm3d_98(nnodes,n);
       break;
     case 1170:
       ele = new (ba) DEMHelm3d_6(nnodes,n);
       break;
     case 1171:
       ele = new (ba) DEMHelm3d_26(nnodes,n);
       break;
     case 1172:
       ele = new (ba) DEMHelm3d_56(nnodes,n);
       break;
     case 1173:
       ele = new (ba) DEMHelm3d_98(nnodes,n);
       break;

     case 1250:
       ele = new (ba) DGMLE3d_6(nnodes,n);
       break;
     case 1251:
       ele = new (ba) DGMLE3d_26(nnodes,n);
       break;
     case 1252:
       ele = new (ba) DGMLE3d_50(nnodes,n);
       break;

// ADDED FOR SLOSHING ELEMENTS, EC, 20070713
     //case 5001:
     case 301:
       ele = new (ba) SloshQuadGal(n);
       break;
     //case 5002:
     case 302:
       ele = new (ba) BarSloshFS(n);
       break;
     //case 5011:
     case 311:
       ele = new (ba) SloshTetra(n);
       break;
     //case 5012:
     case 312:
       ele = new (ba) SloshTriangleFS(n);
       break;
     //case 5101:
     case 321:
       ele = new (ba) HEVibQuadGal(n);
       break;
     //case 5111:
     case 331:
       ele = new (ba) HEVibTetra(n);
       break;
//-----------------------------------------
     default:
       std::cerr << "Element Type " << etype << " is not Supported." << std::endl;
       assert(0);
       break;
   }
   ele->setElementType(etype);
   return ele;
}

void
Elemset::mpcelemadd(int num, LMPCons *mpc)
{
   Element *ele;
   ele = new (ba) MpcElement(mpc);
   ele->setElementType(1001);
   elemadd(num, ele);

   map<int, double >::iterator it = weightList.find(1001);
   if(it == weightList.end()) {
     ele->setWeight(1.0);
     ele->setTrueWeight(1.0);
   } else {
     ele->setWeight(it->second);
     ele->setTrueWeight(it->second);
   }
}

void
Elemset::fsielemadd(int num, LMPCons *fsi)
{
   Element *ele;
   ele = new (ba) FsiElement(fsi);
   ele->setElementType(1002);
   elemadd(num, ele);

   map<int, double >::iterator it = weightList.find(1002);
   if(it == weightList.end()) {
     ele->setWeight(1.0);
     ele->setTrueWeight(1.0);
   } else {
     ele->setWeight(it->second);
     ele->setTrueWeight(it->second);
   }
}

