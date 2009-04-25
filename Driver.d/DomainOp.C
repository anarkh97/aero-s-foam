#include <Driver.d/SubDomain.h>
#include <Driver.d/DomainOp.h>

ThreadManager *threadManager;

DomainTask::DomainTask(Domain *parent,int ns,Connectivity *cn,Connectivity *sToN)
{
 g       = parent;
 subn    = ns;
 connect = cn;
 op      = 0;
 subToN  = sToN;
}

void
DomainTask::run()
{
 if(op)
   (sd->*op)();
 else
   sd = new SubDomain(*g,subn,*connect, *subToN);
}


