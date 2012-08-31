#ifdef STRUCTOPT

#include <cstdlib>
#include <cstdio>
#include <cstring>

#include <Driver.d/Domain.h>
#include <Structopt.d/Driver_opt.d/Domain_opt.h>
#include <Hetero.d/FlExchange.h>

#include <Structopt.d/Optvar.h>
#include <Structopt.d/Optpro.h>
#include <Structopt.d/Relsol.h>
#include <Structopt.d/Structopt_sd.h>

extern Domain * domain;
extern int yyrelparse(void);


//------------------------------------------------------------------------------

void Domain_opt::reliabilitySolve() {

  relpro.reset(new Optpro(1));                 //Define Reliability Problem
  structrel.reset(new Structopt_sd(1));        //Define Structural Rel. Interface

  reliabilityInput();                       // Reading Reliability Input File

  dynamic_cast<Structopt_sd*>(structrel.get())
    ->build(dynamic_cast<Domain_opt*>(domain),relpro.get());           //Build Reliability Problem    
                                             
  relpro->solve(structrel.get());                  //Solve Reliability Problem

  fclose(relunitout);                        //Closing Reliability Output
  fclose(relprotout);

#ifdef AEROELASTIC
#ifdef STRUCTMECH
  extern Communicator *structToelectroCom;
  if (structToelectroCom) structrel->sndOptpar(-1,-1);
#endif
  if (structrel->dynpros)  structrel->dynpros[0]->cmdCom(-1);
#endif

}

//------------------------------------------------------------------------------

void Domain_opt::reloptInitialize(Structopt_sd* structopt) 
{

  relpro.reset(new Optpro(1));                       //Define Reliability Problem
  structrel.reset(new Structopt_sd(1));              //Define Structural Rel. Interface

  reliabilityInput();                              //Reading Reliability Input File

  dynamic_cast<Structopt_sd*>(structrel.get())->
    buildInopt(dynamic_cast<Domain_opt*>(domain),relpro.get(),structopt);  //Build Reliability Problem
                                                   //within Structural Optimization    
}

//------------------------------------------------------------------------------

void Domain_opt::reliabilityInput() {

  //Reading Reliability Analysis Input-File

  FILE *relin = freopen(relinputfile,"r",stdin);
  if (relin == 0) {
    fprintf(stderr,"\n  *** ERROR: Could not open input file: %s\n",
            relinputfile);
    exit(-1);
  }

  //Checking for Errors during Reading
    
  int relerror = yyrelparse();

  if (relerror) {
  fprintf(stderr,
  "\n *** ERROR: Reliability Analysis-Input file contained errors. Aborting ... \n");
  exit(relerror);
  }
  
  //Open Outputfiles for Reliability Analysis
  
  char * basename = getBasename(relinputfile); 
  int fnamesize   = strlen(basename)+4;
  char * relfile  = static_cast<char*>(malloc(sizeof(char*)*(fnamesize)));
  strcpy(relfile,basename);
  strcat(relfile,".rel"); 
  relunitout = fopen(relfile,"w"); 

  relpro->relsol->fsize       = fnamesize;
  relpro->relsol->relprotfile = static_cast<char*>(malloc(sizeof(char*)*(fnamesize)));
  strcpy(relpro->relsol->relprotfile,basename);
  strcat(relpro->relsol->relprotfile,".rpo");   
  relprotout = fopen(relpro->relsol->relprotfile,"w"); 

  // set Outputfile

  relpro->setOutput();

  //Printing Problem

  relpro->print();                          
}

#endif
