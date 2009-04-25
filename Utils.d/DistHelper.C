#include <Utils.d/DistHelper.h>

#ifdef DISTRIBUTED
class SysCom;
extern SysCom *syscom;
extern SysCom *scom;
#include <Comm.d/Communicator.h>
extern Communicator *structCom;
#endif

extern int verboseFlag;
extern int salinasFlag;

void
filePrint(FILE *file, const char *format, ...)
{
 if(!salinasFlag || (salinasFlag && verboseFlag)) {
  va_list args;
  va_start(args, format);
#ifdef DISTRIBUTED 
  // cerr << "structCom->myID() = " << structCom->myID() << endl;
  if( !structCom || structCom->myID() == 0 )
#endif
  vfprintf(file, format, args);
  va_end(args);
 }
}

void
filePrint2(FILE *file, const char *format, ...)
{
#ifdef DISTRIBUTED
 if(!salinasFlag || (salinasFlag && verboseFlag)) {
  va_list args;
  va_start(args, format);
  if( !structCom || structCom->myID() == 0 )
  vfprintf(file, format, args);
  va_end(args);
 }
#endif
}

