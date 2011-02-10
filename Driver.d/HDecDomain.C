#include <cstdlib>
#include <cmath>
#include <Utils.d/dbg_alloca.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include <Utils.d/linkfc.h>
#include <Utils.d/dofset.h>
#include <Utils.d/Memory.h>
#include <Element.d/Sommerfeld.d/SommerElement.h>
#include <Threads.d/PHelper.h>
#include <Math.d/Skyline.d/SkyMatrix.h>
#include <Math.d/SparseMatrix.h>
#include <Solvers.d/Solver.h>
#include <Math.d/matrix.h>
#include <Utils.d/Connectivity.h>
#include <Driver.d/SubDomain.h>
#include <Threads.d/Paral.h>
#include <Driver.d/DomainOp.h>
#include <Solvers.d/Rbm.h>
#include <Driver.d/PolygonSet.h>
#include <Utils.d/OutputInfo.h>

#include <Feti.d/DistrVector.h>
#include <Feti.d/Feti.h>


extern const char* problemTypeMessage[];

#ifdef DISTRIBUTED
extern FILE *debugFile;
#endif

static void cross(double *a, double *b, double *c) 
{
   c[0] = a[1]*b[2]-a[2]*b[1];
   c[1] = a[2]*b[0]-a[0]*b[2];
   c[2] = a[0]*b[1]-a[1]*b[0];
}

