#ifdef STRUCTOPT

#include <Structopt.d/Element_opt.d/Element_opt.h>
#include <Structopt.d/Element_opt.d/Shell.d/ThreeNodeShell_opt.h>
#include <Structopt.d/Element_opt.d/Quad4.d/FourNodeQuad_opt.h>
#include <Structopt.d/Element_opt.d/Brick.d/EightNodeBrick_opt.h>
#include <Structopt.d/Element_opt.d/ViscDamp.d/ViscDamp2D.h>
#include <Structopt.d/Element_opt.d/ViscDamp.d/ViscDamp3D.h>
#include <Structopt.d/Element_opt.d/PMLWave.d/PMLWaveElement2D.h>
#include <Structopt.d/Element_opt.d/PMLWave.d/PMLWaveElement3D.h>
#include <Structopt.d/Element_opt.d/PiezoElements.d/PiezoBrick8.h>
#include <Structopt.d/Element_opt.d/PiezoElements.d/PiezoQuad4.h>
#include <Structopt.d/Element_opt.d/Spring.d/TransSpring_opt.h>
#include <Structopt.d/Element_opt.d/Truss.d/TwoNodeTruss_opt.h>

Element*
ElementFactory_opt::elemadd(int num, int etype, int nnodes, int*n,
			    BlockAlloc& ba)
{
  Element *ele = 0;
  switch(etype) 
    {
     case 1:
       ele = new (ba) TwoNodeTruss_opt(n);
       break;
    case 2:
      ele = new (ba) FourNodeQuad_opt(n);
      break;
    case 8:
      ele = new (ba) ThreeNodeShell_opt(n);
      break;
    case 17:
      ele = new (ba) EightNodeBrick_opt(n);
      break;
     case 21:
       ele = new (ba) TransSprlink_opt(n);
       break;
    case 110:
      ele = new (ba) ViscDamp2D(n);
      break;
    case 111:
      ele = new (ba) ViscDamp3D<4>(n);
      break;
    case 112:
      ele = new (ba) PMLWaveElement2D(n);
      break;
    case 113:
      ele = new (ba) PMLWaveElement3D(n);
      break;
    case 114:
      ele = new (ba) PiezoBrick8(n);
      break;
    case 115:
      ele = new (ba) PiezoQuad4(n);
      break;
    default:
      ele = ElementFactory::elemadd(num, etype, nnodes, n, ba);
    }
  return ele;
}

#endif
