#include <Element.d/Beam.d/RigidMpcBeam.h>
#include <Utils.d/dofset.h>
#include <Driver.d/Mpc.h>

RigidMpcBeam::RigidMpcBeam(int* _nn)
 : RigidMpcElement(2, 6, DofSet::XYZdisp | DofSet::XYZrot, 6, _nn)
{
  /* do nothing */
}

void 
RigidMpcBeam::computeMPCs(CoordSet& cs)
{
  for(int i = 0; i < nMpcs; ++i)
    mpcs[i] = new LMPCons(0, 0.0);

  Node &nd1 = cs.getNode(nn[0]);
  Node &nd2 = cs.getNode(nn[1]);

  double dx = nd2.x - nd1.x;
  double dy = nd2.y - nd1.y;
  double dz = nd2.z - nd1.z;

  double length = sqrt( dx*dx + dy*dy + dz*dz );

  if(length == 0.0) {
    cerr << " *** WARNING: Rigid beam has zero length, nodes: " << nn[0]+1 << " " << nn[1]+1 << endl;
    for(int i = 0; i < 6; ++i) {
      mpcs[i]->addterm(new LMPCTerm(nn[0], i, 1.0));
      mpcs[i]->addterm(new LMPCTerm(nn[1], i, -1.0));
    }
  }
  else {
    double c1[3], c2[3], c3[3];
    c1[0] = dx/length;
    c1[1] = dy/length;
    c1[2] = dz/length;

    // translation in x, y, z (1 constraint equation)
    for(int i = 0; i < 3; ++i) {
      mpcs[0]->addterm(new LMPCTerm(nn[0], i, c1[i]));
      mpcs[0]->addterm(new LMPCTerm(nn[1], i, -c1[i]));
    }

    // rotation in x, y, z (3 constraint equations)
    for(int i = 0; i < 3; ++i) {
      mpcs[i+1]->addterm(new LMPCTerm(nn[0], i+3, 1.0));
      mpcs[i+1]->addterm(new LMPCTerm(nn[1], i+3, -1.0));
    }

    // relative rotation in y and z (2 constraint equations)
    double N1 = sqrt( c1[0]*c1[0] + c1[1]*c1[1] );
    double N2 = sqrt( c1[0]*c1[0] + c1[2]*c1[2] );

    if (N1 > N2) {
      c2[0] = -c1[1]/N1;
      c2[1] = c1[0]/N1;
      c2[2] = 0.0;
    }
    else {
      c2[0] = c1[2]/N2;
      c2[1] = 0.0;
      c2[2] = -c1[0]/N2;
    }

    mpcs[4]->addterm(new LMPCTerm(nn[0], 0, c2[0]));
    mpcs[4]->addterm(new LMPCTerm(nn[0], 1, c2[1]));
    mpcs[4]->addterm(new LMPCTerm(nn[0], 2, c2[2]));
    mpcs[4]->addterm(new LMPCTerm(nn[0], 3, -c2[1] * dz + c2[2] * dy));
    mpcs[4]->addterm(new LMPCTerm(nn[0], 4, c2[0] * dz - c2[2] * dx));
    mpcs[4]->addterm(new LMPCTerm(nn[0], 5, -c2[0] * dy + c2[1] * dx));
    mpcs[4]->addterm(new LMPCTerm(nn[1], 0, -c2[0]));
    mpcs[4]->addterm(new LMPCTerm(nn[1], 1, -c2[1]));
    mpcs[4]->addterm(new LMPCTerm(nn[1], 2, -c2[2]));

    c3[0] = c1[1] * c2[2] - c1[2] * c2[1];
    c3[1] = c1[2] * c2[0] - c1[0] * c2[2];
    c3[2] = c1[0] * c2[1] - c1[1] * c2[0];

    mpcs[5]->addterm(new LMPCTerm(nn[0], 0, c3[0]));
    mpcs[5]->addterm(new LMPCTerm(nn[0], 1, c3[1]));
    mpcs[5]->addterm(new LMPCTerm(nn[0], 2, c3[2]));
    mpcs[5]->addterm(new LMPCTerm(nn[0], 3, -c3[1] * dz + c3[2] * dy));
    mpcs[5]->addterm(new LMPCTerm(nn[0], 4, c3[0] * dz - c3[2] * dx));
    mpcs[5]->addterm(new LMPCTerm(nn[0], 5, -c3[0] * dy + c3[1] * dx));
    mpcs[5]->addterm(new LMPCTerm(nn[1], 0, -c3[0]));
    mpcs[5]->addterm(new LMPCTerm(nn[1], 1, -c3[1]));
    mpcs[5]->addterm(new LMPCTerm(nn[1], 2, -c3[2]));
  }
}

int 
RigidMpcBeam::getTopNumber() 
{ 
  return 106; 
}

