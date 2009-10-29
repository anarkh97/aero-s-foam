#include <Element.d/Solid.d/RigidMpcSolid.h>
#include <Utils.d/dofset.h>
#include <Driver.d/Mpc.h>

RigidMpcSolid::RigidMpcSolid(int _nNodes, int* _nn )
 : RigidMpcElement(_nNodes, 3, DofSet::XYZdisp, 3*(_nNodes-2), _nn)
{
  if(_nNodes < 3) {
    cerr << " *** ERROR: minimum number of nodes is 3 for rigid solid. Exiting...\n";
    exit(-1);
  }
}

void
RigidMpcSolid::computeMPCs(CoordSet &cs)
{
  for(int i = 0; i < nMpcs; ++i)
    mpcs[i] = new LMPCons(0, 0.0);

  double x[4], y[4], z[4];
  int best[3];

  // Finding 3 nodes that make the largest triangle
  //-----------------------------------------------

  // Node 1

  Node &Nd1 = cs.getNode(nn[0]);
  x[0] = Nd1.x; y[0] = Nd1.y; z[0] = Nd1.z;

  double DisTemp1 = 0.0;

  for (int i1 = 1; i1 < nNodes; ++i1){

    Node &nd2 = cs.getNode(nn[i1]);
    x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;

    double dx = x[1] - x[0];
    double dy = y[1] - y[0];
    double dz = z[1] - z[0];
    double length = sqrt(dx * dx + dy * dy + dz * dz);

     if (length >= DisTemp1){
       DisTemp1 = length;
       best[0] = i1;
    }
  }

  Node &nd1 = cs.getNode(nn[best[0]]);
  x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;


  // Node 2

  double DisTemp2 = 0.0;

  for (int i2 = 0; i2 < nNodes; ++i2){

    if (i2 == best[0]) {continue;}

    Node &Nd2 = cs.getNode(nn[i2]);
    x[1] = Nd2.x; y[1] = Nd2.y; z[1] = Nd2.z;

    double dx = x[1] - x[0];
    double dy = y[1] - y[0];
    double dz = z[1] - z[0];
    double length = sqrt(dx * dx + dy * dy + dz * dz);

    if (length >= DisTemp2){
      DisTemp2 = length;
      best[1] = i2;
    }
  }

  Node &nd2 = cs.getNode(nn[best[1]]);
  x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;

  double dx_12 = x[1] - x[0];
  double dy_12 = y[1] - y[0];
  double dz_12 = z[1] - z[0];


  // Node 3

  double LargestArea = 0.0;

  for (int i3 = 0; i3 < nNodes; ++i3) {
    if (i3 == best[0] || i3 == best[1]) { continue; }

    Node &Nd3 = cs.getNode(nn[i3]);
    x[2] = Nd3.x; y[2] = Nd3.y; z[2] = Nd3.z;

    double dx_13 = x[2] - x[0];
    double dy_13 = y[2] - y[0];
    double dz_13 = z[2] - z[0];

    double diff1 = dy_12 * dz_13 - dy_13 * dz_12;
    double diff2 = dx_13 * dz_12 - dx_12 * dz_13;
    double diff3 = dx_12 * dy_13 - dx_13 * dy_12;

    double Area = sqrt( diff1 * diff1 + diff2 * diff2 + diff3 * diff3 );
    if (Area >= LargestArea) {
      LargestArea = Area;
      best[2] = i3;
    }
  }

  Node &nd3 = cs.getNode(nn[best[2]]);
  x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

  double dx_13 = x[2] - x[0];
  double dy_13 = y[2] - y[0];
  double dz_13 = z[2] - z[0];


  // cross product of d_12 and d_13 

  double cross[3];

  cross[0] = dy_12 * dz_13 - dy_13 * dz_12;
  cross[1] = dx_13 * dz_12 - dx_12 * dz_13;
  cross[2] = dx_12 * dy_13 - dx_13 * dy_12;
  double crossNorm = sqrt( cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]);


  // restraining the triangle        

  double c[3];
  int index = 0;

  for (int j = 0; j < 3; ++j){
    int n1 = best[j];
    for (int jj = j + 1; jj < 3; ++jj) {
      int n2 = best[jj];
      double dx = x[jj] - x[j];
      double dy = y[jj] - y[j];
      double dz = z[jj] - z[j];
      double length = sqrt( dx * dx + dy * dy + dz * dz);

      c[0] = dx/length;
      c[1] = dy/length;
      c[2] = dz/length;

      mpcs[index]->addterm(new LMPCTerm(nn[n1], 0, c[0])); 
      mpcs[index]->addterm(new LMPCTerm(nn[n1], 1, c[1]));
      mpcs[index]->addterm(new LMPCTerm(nn[n1], 2, c[2]));
      mpcs[index]->addterm(new LMPCTerm(nn[n2], 0, -c[0])); 
      mpcs[index]->addterm(new LMPCTerm(nn[n2], 1, -c[1])); 
      mpcs[index]->addterm(new LMPCTerm(nn[n2], 2, -c[2])); 

      ++index;
    }
  }


  // Restraining the rest of the nodes

  double dot_1212 = dx_12 * dx_12 + dy_12 * dy_12 + dz_12 * dz_12;
  double dot_1213 = dx_12 * dx_13 + dy_12 * dy_13 + dz_12 * dz_13;
  double dot_1313 = dx_13 * dx_13 + dy_13 * dy_13 + dz_13 * dz_13;

  int N1 = best[0];
  int N2 = best[1];
  int N3 = best[2];


  for (int i = 0; i < nNodes; ++i) {

    if (i == best[0] || i == best[1] || i == best[2]) { continue; }

    Node &nd4 = cs.getNode(nn[i]);
    x[3] = nd4.x; y[3] = nd4.y; z[3] = nd4.z;

    double dx_14 = x[3] - x[0];
    double dy_14 = y[3] - y[0];
    double dz_14 = z[3] - z[0];

    // finding gama for the constraint equation

    double mult = dx_14 * cross[0] + dy_14 * cross[1] + dz_14 * cross[2];
    double gama = mult / (crossNorm * crossNorm);

    // finding alpha and beta for the constraint equation

    double dot_1412 = dx_14 * dx_12 + dy_14 * dy_12 + dz_14 * dz_12;
    double dot_1413 = dx_14 * dx_13 + dy_14 * dy_13 + dz_14 * dz_13;

    double betaNum = dot_1413 - (dot_1412 * dot_1213) / dot_1212;
    double betaDenom = dot_1313 - (dot_1213 * dot_1213) / dot_1212;

    // BETA
    double beta = betaNum / betaDenom;

    // ALPHA
    double alpha = dot_1412 / dot_1212 - beta * ( dot_1213 / dot_1212);

    // Stiffness

    int N4 = i;

    mpcs[index]->addterm(new LMPCTerm(nn[N1], 0, - (1 - alpha - beta)));
    mpcs[index]->addterm(new LMPCTerm(nn[N1], 1, gama * (dz_13 - dz_12)));
    mpcs[index]->addterm(new LMPCTerm(nn[N1], 2, gama * (dy_12 - dy_13)));

    mpcs[index]->addterm(new LMPCTerm(nn[N2], 0, - alpha));
    mpcs[index]->addterm(new LMPCTerm(nn[N2], 1, - gama * dz_13));
    mpcs[index]->addterm(new LMPCTerm(nn[N2], 2, gama * dy_13));

    mpcs[index]->addterm(new LMPCTerm(nn[N3], 0, - beta));
    mpcs[index]->addterm(new LMPCTerm(nn[N3], 1, gama * dz_12));
    mpcs[index]->addterm(new LMPCTerm(nn[N3], 2, - gama * dy_12));

    mpcs[index]->addterm(new LMPCTerm(nn[N4], 0, 1.0));

    //

    mpcs[index + 1]->addterm(new LMPCTerm(nn[N1], 0, gama * (dz_12 - dz_13)));
    mpcs[index + 1]->addterm(new LMPCTerm(nn[N1], 1, -(1 - alpha - beta)));
    mpcs[index + 1]->addterm(new LMPCTerm(nn[N1], 2, gama * (dx_13 - dx_12)));

    mpcs[index + 1]->addterm(new LMPCTerm(nn[N2], 0, gama * dz_13));
    mpcs[index + 1]->addterm(new LMPCTerm(nn[N2], 1, - alpha));
    mpcs[index + 1]->addterm(new LMPCTerm(nn[N2], 2, - gama * dx_13));

    mpcs[index + 1]->addterm(new LMPCTerm(nn[N3], 0, - gama * dz_12));
    mpcs[index + 1]->addterm(new LMPCTerm(nn[N3], 1, - beta));
    mpcs[index + 1]->addterm(new LMPCTerm(nn[N3], 2, gama * dx_12));

    mpcs[index + 1]->addterm(new LMPCTerm(nn[N4], 1, 1.0));

    // 

    mpcs[index + 2]->addterm(new LMPCTerm(nn[N1], 0, gama * (dy_13 - dy_12)));
    mpcs[index + 2]->addterm(new LMPCTerm(nn[N1], 1, gama * (dx_12 - dx_13)));
    mpcs[index + 2]->addterm(new LMPCTerm(nn[N1], 2, - (1 - alpha - beta)));

    mpcs[index + 2]->addterm(new LMPCTerm(nn[N2], 0, - gama * dy_13));
    mpcs[index + 2]->addterm(new LMPCTerm(nn[N2], 1, gama * dx_13));
    mpcs[index + 2]->addterm(new LMPCTerm(nn[N2], 2, - alpha));

    mpcs[index + 2]->addterm(new LMPCTerm(nn[N3], 0, gama * dy_12));
    mpcs[index + 2]->addterm(new LMPCTerm(nn[N3], 1, - gama * dx_12));
    mpcs[index + 2]->addterm(new LMPCTerm(nn[N3], 2, - beta));

    mpcs[index + 2]->addterm(new LMPCTerm(nn[N4], 2, 1.0));

    index += 3;
  }
}

int 
RigidMpcSolid::getTopNumber()
{ 
  switch(nNodes) {
    case 4 : return 123; // 4-node tetra
    case 6 : return 124; // 6-node penta
    case 8 : return 117; // 8-node hexa
    case 10 : return 125; // 10-node tetra
    case 15 : return 197; // 15-node hexa
    case 20 : return 172; // 20-node brick
    case 26 : return 192; // 26-node hexa 
    case 32 : return 191; // 32-node wedge
    default : return 101;
  }  
}

int 
RigidMpcSolid::numTopNodes() 
{ 
  switch(nNodes) {
    case 4 : return 4; // 4-node tetra
    case 6 : return 6; // 6-node penta
    case 8 : return 8; // 8-node hexa
    case 10 : return 10; // 10-node tetra
    case 15 : return 15; // 15-node hexa
    case 20 : return 20; // 20-node brick
    case 26 : return 26; // 26-node hexa 
    case 32 : return 32; // 32-node wedge
    default : return 2;
  } 
}
