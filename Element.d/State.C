#include <Math.d/Vector.h>
#include <Element.d/State.h>

Vector nillVec;

void
State::getDV(int node, double xyz[3], double v[3])
{
 int loc1 = DSA->locate(node,DofSet::Xdisp);
 int loc  = dsa->locate(node,DofSet::Xdisp);

 if(loc >= 0) { 	// free
   xyz[0] = disp[loc];
     v[0] = veloc[loc];
 } else if(loc1 >= 0) { // prescribed
   xyz[0] = bcx[loc1];
     v[0] = vcx[loc1];
 } else {		// not defined
   xyz[0] = 0.0;
     v[0] = 0.0;
   }

 loc1 = DSA->locate(node,DofSet::Ydisp);
 loc  = dsa->locate(node,DofSet::Ydisp);

 if(loc >= 0) {         // free
   xyz[1] = disp[loc];
     v[1] = veloc[loc];
 } else if(loc1 >= 0) { // prescribed
   xyz[1] = bcx[loc1];
     v[1] = vcx[loc1];
 } else {               // not defined
   xyz[1] = 0.0;
     v[1] = 0.0;
   }


 loc1 = DSA->locate(node,DofSet::Zdisp);
 loc  = dsa->locate(node,DofSet::Zdisp);

 if(loc >= 0) {         // free
   xyz[2] = disp[loc];
     v[2] = veloc[loc];
 } else if(loc1 >= 0) { // prescribed
   xyz[2] = bcx[loc1];
     v[2] = vcx[loc1];
 } else {               // not defined
   xyz[2] = 0.0;
     v[2] = 0.0;
   }

}

void
State::getTemp(int node, double Temp[1], double dTempdt[1])
{
 int tloc1 = DSA->locate(node,DofSet::Temp);
 int tloc  = dsa->locate(node,DofSet::Temp);

 if(tloc >= 0) {         // free
//   Temp = d_n[tloc];
//   dTempdt = v_n[tloc];
   Temp[0] = disp[tloc];
   dTempdt[0] = veloc[tloc];

 } else if(tloc1 >= 0) { // prescribed(constrained)
   Temp[0] = bcx[tloc1];
   dTempdt[0] = 0.0; //vcx[tloc1];
 } else {               // not defined
   Temp[0] = 0.0;
   dTempdt[0] = 0.0;
   }
}


void
State::getDVRot(int node, double xyz[6], double v[6])
{
 int loc;
 loc = dsa->locate(node,DofSet::Xdisp);
 xyz[0] = (loc >= 0) ? disp[loc]  : 0;
   v[0] = (loc >= 0) ? veloc[loc] : 0;
 loc = dsa->locate(node,DofSet::Ydisp);
 xyz[1] = (loc >= 0) ? disp[loc]  : 0;
   v[1] = (loc >= 0) ? veloc[loc] : 0;
 loc = dsa->locate(node,DofSet::Zdisp);
 xyz[2] = (loc >= 0) ? disp[loc]  : 0;
   v[2] = (loc >= 0) ? veloc[loc] : 0;
 loc = dsa->locate(node,DofSet::Xrot);
 xyz[3] = (loc >= 0) ? disp[loc]  : 0;
   v[3] = (loc >= 0) ? veloc[loc] : 0;
 loc = dsa->locate(node,DofSet::Yrot);
 xyz[4] = (loc >= 0) ? disp[loc]  : 0;
   v[4] = (loc >= 0) ? veloc[loc] : 0;
 loc = dsa->locate(node,DofSet::Zrot);
 xyz[5] = (loc >= 0) ? disp[loc]  : 0;
   v[5] = (loc >= 0) ? veloc[loc] : 0;
}

//------------------------------------------------------------------

