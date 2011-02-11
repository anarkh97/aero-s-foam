#ifdef STRUCTOPT

#include <cstdlib>
#include <cstdio>
#include <string.h>
#include <Utils.d/dbg_alloca.h>

#include <Driver.d/Domain.h>
#include <Structopt.d/Driver_opt.d/Domain_opt.h>
#include <Structopt.d/Driver_opt.d/GeoSource_opt.h>
#include <Hetero.d/FlExchange.h>

#include <Structopt.d/Optvar.h>
#include <Structopt.d/Optpro.h>
#include <Structopt.d/Optsol.h>
#include <Structopt.d/Structopt_sd.h>

#include <Structopt.d/Element_opt.d/Element_opt.h>

extern Domain * domain;
extern int yyoptparse(void);

//------------------------------------------------------------------------------
void Domain_opt::structoptSolve() {

  optpro.reset(new Optpro(0));                 //Define Optimization Problem
  structopt.reset(new Structopt_sd(0, this, optpro.get()));        //Define Structural Opt. Interface

  structoptInput();                          // Reading Optimiztion Input File 
  static_cast<Structopt_sd*>(structopt.get())
    ->build(dynamic_cast<Domain_opt*>(domain), optpro.get());           //Build Structural Optimization    
                                             
  optpro->solve(structopt.get());                  //Solve Optimization problem

  fclose(optunitout);                        //Closing Optimization Output
  fclose(optprotout);

#ifdef AEROELASTIC
#ifdef STRUCTMECH
  extern Communicator *structToelectroCom;
  if (structToelectroCom) structopt->sndOptpar(-1,-1);
#endif
  if (sinfo.thermoeFlag >=0){
     structopt->sndOptpar(-1,-1);
     }
  if (structopt->dynpros)   structopt->dynpros[0]->cmdCom(-1);
  if (structopt->NLdynpros) structopt->NLdynpros[0]->cmdCom(-1);
#endif
}

//------------------------------------------------------------------------------
void Domain_opt::structoptInput() 
{
  //Reading Optimization Input-File
  FILE *optin = freopen(this->optinputfile,"r",stdin);
  if (optin == 0) {
    fprintf(stderr,"\n  *** ERROR: Could not open input file: %s\n",
            this->optinputfile);
    exit(-1);
  }

  //Checking for Errors during Reading
    
  int opterror = yyoptparse();

  if (opterror) {
    fprintf(stderr,
	    "\n *** ERROR: Optimization-Input file contained errors. Aborting ... \n");
    exit(opterror);
  }
  
  //Open Outputfiles for Optimization
  char * basename = getBasename(optinputfile); 
  int fnamesize   = strlen(basename)+5;
  boost::scoped_array<char> optfile(new char[fnamesize]);
  strcpy(optfile.get(),basename);
  strcat(optfile.get(),".opt"); 
  optunitout = fopen(optfile.get(),"w"); 

  optpro->optsol->fsize       = fnamesize;
  optpro->optsol->optprotfile = static_cast<char*>(malloc(sizeof(char*)*(fnamesize)));
  strcpy(optpro->optsol->optprotfile,basename);
  strcat(optpro->optsol->optprotfile,".nlp");   
  optprotout = fopen(optpro->optsol->optprotfile,"w"); 

  // set Outputfile

  optpro->setOutput();

  //Printing Problem

  optpro->print();                          
}

//------------------------------------------------------------------------------
void
Domain_opt::postProcessing(Vector &sol, double *bcx, Vector& force, 
			   double time)
{
  OutputInfo *oinfo = geoSource->getOutputInfo();
  int numOutInfo = geoSource->getNumOutInfo();
  
  //double oldtime=sinfo.dt;
  //sinfo.dt=time;
  
   int exactNumNodes  = 0;
   int lastNodeNumber = nodes.size();

   int *nodeTable = new int[lastNodeNumber];

   int inode;
   for(inode=0; inode<lastNodeNumber; ++inode)
    nodeTable[inode] = -1;

   for(inode=0; inode<lastNodeNumber; ++inode) {
     if(nodes[inode] == 0) continue;
     exactNumNodes += 1;
     nodeTable[inode] = exactNumNodes;
   }
   
   Node oldnod;

 // Wext = external energy
 // Wela = elastic energy
 // Wkin = kinetic energy
 // Wdmp = damping energy
 // Total Energy = Wext+Wela+Wkin+Wdmp 

 double Wext=0.0,Wela=0.0,Wkin=0.0,Wdmp=0.0;
 double zero = 0.0;

 int idLC = actLC;

  // organize displacements

  int i;

  assert(dynamic_cast<Domain_opt*>(domain) != 0);
  double (*xyz)[11] = new double[dynamic_cast<Domain_opt*>(domain)->numnodes][11];

  for (i = 0; i < numnodes; ++i) {
    xyz[i][0] = xyz[i][1] = xyz[i][2] = xyz[i][3] = xyz[i][4] = xyz[i][5] = 0.0;
    xyz[i][6] = xyz[i][7] = xyz[i][8] = xyz[i][9] = xyz[i][10] = 0.0;
  }

  this->mergeDistributedDisp(xyz, sol.data(), bcx); 

  // loop over all output data

 for(i=0; i<numOutInfo; ++i) 
 {
   int w = oinfo[i].width;
   int p = oinfo[i].precision;

   // check if output file is assigned to current loadcase
   if( oinfo[i].loadcase != idLC && oinfo[i].loadcase > -1 ) 
     { continue; }

   if(oinfo[i].interval == 1) {
     switch(oinfo[i].type)
     {
       case OutputInfo::Displacement:
     	 if (oinfo[i].nodeNumber == -1) {

           double x,y,z;
           fprintf(oinfo[i].filptr,"  %f\n",time);

           for(inode=0; inode<lastNodeNumber; ++inode) {

             if(nodes[inode] == 0) continue;

             int xloc  = c_dsa->locate( inode, DofSet::Xdisp);
             int xloc1 =   dsa->locate( inode, DofSet::Xdisp);

	     if(xloc >= 0) 
               x = sol[xloc];
 	     else if(xloc1 >= 0)
               x = bcx[xloc1];
             else
               x = 0.0;

             int yloc  = c_dsa->locate( inode, DofSet::Ydisp);
             int yloc1 =   dsa->locate( inode, DofSet::Ydisp);

	     if(yloc >= 0) 
               y = sol[yloc];
 	     else if(xloc1 >= 0)
               y = bcx[yloc1];
             else
               y = 0.0;

             int zloc  = c_dsa->locate( inode, DofSet::Zdisp);
             int zloc1 =   dsa->locate( inode, DofSet::Zdisp);

	     if(zloc >= 0) 
               z = sol[zloc];
 	     else if(zloc1 >= 0)
               z = bcx[zloc1];
             else
               z = 0.0;

             // ... considering the variation of shape
   	     oldnod=nodescopy->getNode(inode);
             double vx = nodes[inode]->x - oldnod.x ;
             double vy = nodes[inode]->y - oldnod.y ;
             double vz = nodes[inode]->z - oldnod.z ;
	    
	     x=x+vx;   y=y+vy;    z=z+vz;

             fprintf(oinfo[i].filptr,"% *.*E\t% *.*E\t% *.*E\n"
                                    ,w,p,x,w,p,y,w,p,z);
           }
           fflush(oinfo[i].filptr);
         }
       	 else { 
           int iNode = oinfo[i].nodeNumber;
           fprintf(stderr," .... WATCH: node variation not conidered in output\n");
            geoSource->outputNodeVectors(i, xyz+(7*iNode), 1, zero);
          }
          break;
       case OutputInfo::Disp6DOF:
          if (oinfo[i].nodeNumber == -1) {
            fprintf(oinfo[i].filptr,"  %f\n",time);
            for(inode=0; inode<lastNodeNumber; ++inode) {
   
               if(nodes[inode] == 0) continue;
   
               int xloc  = c_dsa->locate( inode, DofSet::Xdisp);
               int xloc1 =   dsa->locate( inode, DofSet::Xdisp);
   
               double x,y,z,xr,yr,zr;
   
               if(xloc >= 0)	// dof exists and is free
                 x = sol[xloc];
               else if(xloc1 >= 0)	// dof exists and is constrained
                 x = bcx[xloc1];
               else	        // dof does not exist
                 x = 0.0;
   
               int yloc  = c_dsa->locate( inode, DofSet::Ydisp);
               int yloc1 =   dsa->locate( inode, DofSet::Ydisp);
   
               if(yloc >= 0)
                 y = sol[yloc];
               else if(xloc1 >= 0)
                 y = bcx[yloc1];
               else
                 y = 0.0;
   
               int zloc  = c_dsa->locate( inode, DofSet::Zdisp);
               int zloc1 =   dsa->locate( inode, DofSet::Zdisp);
   
   	       if(zloc >= 0) 
                 z = sol[zloc];
    	       else if(zloc1 >= 0)
                 z = bcx[zloc1];
               else
                 z = 0.0;

               int xrot  = c_dsa->locate( inode, DofSet::Xrot);
               int xrot1 =   dsa->locate( inode, DofSet::Xrot);
   
               if(xrot >= 0)
                 xr = sol[xrot];
               else if(xrot1 >= 0)
                 xr = bcx[xrot1];
               else
                 xr = 0.0;
   
               int yrot  = c_dsa->locate( inode, DofSet::Yrot);
               int yrot1 =   dsa->locate( inode, DofSet::Yrot);
   
               if(yrot >= 0)
                 yr = sol[yrot];
               else if(yrot1 >= 0)
                 yr = bcx[yrot1];
               else
                 yr = 0.0;
   
               int zrot  = c_dsa->locate( inode, DofSet::Zrot);
               int zrot1 =   dsa->locate( inode, DofSet::Zrot);
   
               if(zrot >= 0)
                 zr = sol[zrot];
               else if(zrot1 >= 0)
                 zr = bcx[zrot1];
               else
                 zr = 0.0;
   
               // ... considering the variation of shape
      	       oldnod=nodescopy->getNode(inode);
               double vx = nodes[inode]->x - oldnod.x ;
               double vy = nodes[inode]->y - oldnod.y ;
               double vz = nodes[inode]->z - oldnod.z ;
	    
   	       x=x+vx;   y=y+vy;    z=z+vz;

               fprintf(oinfo[i].filptr,
                 "%d % *.*E\t% *.*E\t% *.*E % *.*E\t % *.*E\t % *.*E\n"
                 ,nodeTable[inode], w,p,x,w,p,y,w,p,z,w,p,xr,w,p,yr,w,p,zr);
            }
            fflush(oinfo[i].filptr);
          }
          else  {
            int iNode = oinfo[i].nodeNumber;
            geoSource->outputNodeVectors6(i, xyz+(7*iNode), 1, zero);
          }
          break;
       case OutputInfo::Temperature:
         fprintf(oinfo[i].filptr,"  %f\n",time);
          if(oinfo[i].nodeNumber == -1) {
            geoSource->outputNodeScalars(i, xyz[0], 0, zero);
            int inode;
            for (inode = 0; inode < numnodes; ++inode)
              geoSource->outputNodeScalars(i, xyz[inode]+6, 1);
          } 
  	  else  {
            // Only one node was requested for output
            int iNode = oinfo[i].nodeNumber;
 	    geoSource->outputNodeScalars(i, xyz[iNode]+6, 1, zero);
          }
          break;
       case OutputInfo::StressXX:
         getStressStrain(sol,bcx,i,SXX);
         break;
       case OutputInfo::StressYY:
         getStressStrain(sol,bcx,i,SYY);
         break;
       case OutputInfo::StressZZ:
         getStressStrain(sol,bcx,i,SZZ);
         break;
       case OutputInfo::StressXY:
         getStressStrain(sol,bcx,i,SXY);
         break;
       case OutputInfo::StressYZ:
         getStressStrain(sol,bcx,i,SYZ);
         break;
       case OutputInfo::StressXZ:
         getStressStrain(sol,bcx,i,SXZ);
         break;
       case OutputInfo::StrainXX:
         getStressStrain(sol,bcx,i,EXX);
         break;
       case OutputInfo::StrainYY:
         getStressStrain(sol,bcx,i,EYY);
         break;
       case OutputInfo::StrainZZ:
         getStressStrain(sol,bcx,i,EZZ);
         break;
       case OutputInfo::StrainXY:
         getStressStrain(sol,bcx,i,EXY);
         break;
       case OutputInfo::StrainYZ:
         getStressStrain(sol,bcx,i,EYZ);
         break;
       case OutputInfo::StrainXZ:
         getStressStrain(sol,bcx,i,EXZ);
         break;
       case OutputInfo::StressVM:
         getStressStrain(sol,bcx,i,VON);
         break;
       case OutputInfo::StressPR1:
         getPrincipalStress(sol,bcx,i,PSTRESS1);
         break;
       case OutputInfo::StressPR2:
         getPrincipalStress(sol,bcx,i,PSTRESS2);
         break;
       case OutputInfo::StressPR3:
         getPrincipalStress(sol,bcx,i,PSTRESS3);
         break;
       case OutputInfo::StrainPR1:
         getPrincipalStress(sol,bcx,i,PSTRAIN1);
         break;
       case OutputInfo::StrainPR2:
         getPrincipalStress(sol,bcx,i,PSTRAIN2);
         break;
       case OutputInfo::StrainPR3:
         getPrincipalStress(sol,bcx,i,PSTRAIN3);
         break;
       case OutputInfo::InXForce:
         getElementForces(sol, bcx, i, INX);
         break;
       case OutputInfo::InYForce:
         getElementForces(sol, bcx, i, INY);
         break;
       case OutputInfo::InZForce:
         getElementForces(sol, bcx, i, INZ);
         break;
       case OutputInfo::AXMoment:
         getElementForces(sol, bcx, i, AXM);
         break;
       case OutputInfo::AYMoment:
         getElementForces(sol, bcx, i, AYM);
         break;
       case OutputInfo::AZMoment:
         getElementForces(sol, bcx, i, AZM);
         break;
       case OutputInfo::Energies: {
          Wext = force *  sol;   // Wext = external energy
          Wela =   0.5 * Wext;   // Wela = elastic energy 
          double error = Wext+Wela+Wkin+Wdmp;
          static_cast<GeoSource_opt*>(geoSource)->outputEnergies(i, zero, Wext, Wela, Wkin, Wdmp, error); }
         break;
       case OutputInfo::StrainVM:
	 getStressStrain(sol,bcx,i,STRAINVON);
         break;
       case OutputInfo::YModulus:
         getElementAttr(i,YOUNG, time);
         break;
       case OutputInfo::MDensity:
         getElementAttr(i,MDENS, time);
         break;
       case OutputInfo::Thicknes:
         getElementAttr(i,THICK, time);
         break;
       case OutputInfo::Composit:
         getCompositeData(i,time);
         break;
	 /*
       case OutputInfo::Ctexp:
         getElementAttr(i,CTEXP);
         break;
	 */
       case OutputInfo::Helmholtz:
         // ... PRINT (REAL) HELMHOLTZ SOLUTION
         for(inode=0; inode<numnodes; ++inode) {
              int loc  = c_dsa->locate( inode, DofSet::Helm);
              int loc1 =   dsa->locate( inode, DofSet::Helm);

              double xHelm;
              if(loc >= 0)        // dof exists and is free
                xHelm = sol[loc];
              else if(loc1 >= 0)  // dof exists and is constrained
                xHelm = bcx[loc1];
              else                // dof does not exist
                xHelm = 0.0;

              fprintf(oinfo[i].filptr,"% *.*E\n",w,p,xHelm);
         }
         break;
       case OutputInfo::ShapeAtt:
         fprintf(oinfo[i].filptr,"  %f\n",time);
	 
         for(inode=0; inode<lastNodeNumber; ++inode) {

           if(nodes[inode] == 0) continue;

	   oldnod=nodescopy->getNode(inode);
           double x = nodes[inode]->x - oldnod.x ;
           double y = nodes[inode]->y - oldnod.y ;
           double z = nodes[inode]->z - oldnod.z ;

           fprintf(oinfo[i].filptr,"% *.*E\t% *.*E\t% *.*E\n"
                               ,w,p,x,w,p,y,w,p,z);
	 }
	 break;
       case OutputInfo::ShapeStc:
         fprintf(oinfo[i].filptr,"  %f\n",time);
	 
         for(inode=0; inode<lastNodeNumber; ++inode) {

           if(nodes[inode] == 0) continue;

 	   oldnod=nodescopy->getNode(inode);
           double x = nodes[inode]->x - oldnod.x ;
           double y = nodes[inode]->y - oldnod.y ;
           double z = nodes[inode]->z - oldnod.z ;

           fprintf(oinfo[i].filptr,"% *.*E\t% *.*E\t% *.*E\n"
                               ,w,p,x,w,p,y,w,p,z);
	 }
	 break;
     default:
       fprintf(stderr, "---------------> %d\n", oinfo[i].type);
       assert(0);       
       break;
     }
    }
  } 

 //sinfo.dt=oldtime;

 // delete temporary arrays
 delete [] nodeTable;
 delete [] xyz;

}

/*------------------------------------------------------------------------------
void
Domain_opt::eigenOutput(Vector& eigenValues, VectorSet& eigenVectors,
                    double time)
{
  OutputInfo *oinfo = geoSource->getOutputInfo();
  int numOutInfo = geoSource->getNumOutInfo();

  int i,inode;
  const double pi = 3.141592653589793;

  enum {YOUNG,MDENS,THICK};

  Node oldnod;

  int iInfo;
  for(iInfo = 0; iInfo < numOutInfo; ++iInfo)
  {
    if(oinfo[iInfo].interval == 1) {
      int w = oinfo[iInfo].width;
      int p = oinfo[iInfo].precision;

      switch(oinfo[iInfo].type) {
      
      case OutputInfo::EigenPair:
      case OutputInfo::Disp6DOF:
 
        // ... Print eigenvalues and eigenvectors
        int imode,firstmode,maxmode;
      
        if (oinfo[iInfo].nodeNumber == -1 ) {
          firstmode=0;  
	  maxmode =sinfo.nEig; }
        else {
          firstmode=oinfo[iInfo].nodeNumber;
	  maxmode=firstmode+1; }
	  
        for(imode=firstmode; imode < maxmode; ++imode) {

           fprintf(oinfo[iInfo].filptr,"  %f\n",sqrt(eigenValues[imode])/(2.0*pi));

           for(i=0; i<numnodes; ++i) {

             // If you want to skip nodes that are not defined.
             if( dsa->firstdof(i) == -1 ) continue;

               int xloc = c_dsa->locate( i, DofSet::Xdisp);
               double x = (xloc >= 0) ? eigenVectors[imode][xloc] : 0;
 
               int yloc = c_dsa->locate( i, DofSet::Ydisp);
               double y = (yloc >= 0) ? eigenVectors[imode][yloc] : 0;

               int zloc = c_dsa->locate( i, DofSet::Zdisp);
               double z = (zloc >= 0) ? eigenVectors[imode][zloc] : 0;

               // ... considering the variation of shape
   	       oldnod=nodescopy->getNode(i);
               double vx = nodes[i]->x - oldnod.x ;
               double vy = nodes[i]->y - oldnod.y ;
               double vz = nodes[i]->z - oldnod.z ;
	   
	       x=x+vx;   y=y+vy;    z=z+vz;

               if(oinfo[iInfo].type == OutputInfo::Disp6DOF) {
                  xloc  = c_dsa->locate( i, DofSet::Xrot);
                  double xrot = (xloc >= 0) ? eigenVectors[imode][xloc] : 0;

                  yloc  = c_dsa->locate( i, DofSet::Yrot);
                  double yrot = (yloc >= 0) ? eigenVectors[imode][yloc] : 0;

                  zloc  = c_dsa->locate( i, DofSet::Zrot);
                  double zrot = (zloc >= 0) ? eigenVectors[imode][zloc] : 0;

                  fprintf(oinfo[iInfo].filptr,
                         "%d % *.*E\t% *.*E\t% *.*E\t% *.*E\t% *.*E\t% *.*E\n",
                         i+1,w,p,x,w,p,y,w,p,z,w,p,xrot,w,p,yrot,w,p,zrot);
               } else
                  fprintf(oinfo[iInfo].filptr,"% *.*E\t% *.*E\t% *.*E\n",
                                         w,p,x,w,p,y,w,p,z);
             }
          }
          fflush(oinfo[iInfo].filptr);
          break;
        case OutputInfo::ShapeAtt:
          fprintf(oinfo[iInfo].filptr,"  %f\n",time);

          for(inode=0; inode<numnodes; ++inode) {
	    oldnod=nodescopy->getNode(inode);
            double x = nodes[inode]->x - oldnod.x ;
            double y = nodes[inode]->y - oldnod.y ;
            double z = nodes[inode]->z - oldnod.z ;

            fprintf(oinfo[iInfo].filptr,"% *.*E\t% *.*E\t% *.*E\n"
                                ,w,p,x,w,p,y,w,p,z);
	  }
	  break;
        case OutputInfo::ShapeStc:
          fprintf(oinfo[iInfo].filptr,"  %f\n",time);

          for(inode=0; inode<numnodes; ++inode) {
	    oldnod=nodescopy->getNode(inode);
            double x = nodes[inode]->x - oldnod.x ;
            double y = nodes[inode]->y - oldnod.y ;
            double z = nodes[inode]->z - oldnod.z ;

            fprintf(oinfo[iInfo].filptr,"% *.*E\t% *.*E\t% *.*E\n"
                               ,w,p,x,w,p,y,w,p,z);
	  }
	  break;
        case OutputInfo::YModulus:
          getElementAttr(iInfo,YOUNG);
          break;
        case OutputInfo::MDensity:
          getElementAttr(iInfo,MDENS);
          break;
        case OutputInfo::Thicknes:
          getElementAttr(iInfo,THICK);
          break;
       }
     }  
   }

 // Print omega and frequency values to screen
 fprintf(stderr," Mode\tOmega^2\t\tFrequency\n");
 fprintf(stderr," --------------------------------------\n");
 int imode;
 for(imode=0; imode<sinfo.nEig; ++imode)
   fprintf(stderr," %d\t%e\t%e\n",imode+1,eigenValues[imode],
           sqrt(eigenValues[imode])/(2.0*pi));

}
*/

//------------------------------------------------------------------------------

void
Domain_opt::buildOptGrad() 
{
   // .... check if arrays for analytical SA exist

   if ( optgradFlag ) return;

   // .... variable position of fe-nodes

   gradnodes.reset(new CoordSet(numnodes));

   double *xyz = static_cast<double*>(dbg_alloca(sizeof(double)*3));
   
   xyz[0]=0.0; xyz[1]=0.0; xyz[2]=0.0; 
   
   int i;
   for (i=0;i<numnodes;i++) { gradnodes->nodeadd(i,xyz) ; }

   // .... allocate memory for SA element routines 
 
  int maxNodesPerElement = 20;

   if(elDisp    == 0)  elDisp    = new Vector(maxNumDOFs);
   if(elForce.get() == 0)  elForce.reset(new Vector(maxNumDOFs));
   if(elstress == 0)  elstress = new Vector(maxNodesPerElement);
   if(elweight == 0)  elweight = new Vector(maxNodesPerElement);
   /* if(elNodTemp == 0)  elNodTemp = new Vector(maxNumNodes,0.0); */
 
   elKelD.reset(new Vector(maxNumDOFs));
   elAdj.reset(new Vector(maxNumDOFs));
   elGrad.reset(new Vector(maxNumDOFs));
   elRHS.reset(new Vector(maxNumDOFs));	 

   kelData.reset(new double[maxNumDOFs*maxNumDOFs]); 
   gkelData.reset(new double[maxNumDOFs*maxNumDOFs]);

   gradelstress.reset(new Vector(maxNodesPerElement));

   c_elDisp.reset(new ComplexVector(maxNumDOFs));
   c_elKelD.reset(new ComplexVector(maxNumDOFs));
   c_elAdj.reset(new ComplexVector(maxNumDOFs));
   c_elGrad.reset(new ComplexVector(maxNumDOFs));
   c_elRHS.reset(new ComplexVector(maxNumDOFs));

   /*
     gradelNodTemp = new Vector(maxNumNodes);
   */

   // .... variable load (Neumann boundary conditons)

   gradnbcLC.reset(new boost::scoped_array<BCond>[numLC]);
   gradcnbcLC.reset(new boost::scoped_array<ComplexBCond>[numLC]);
   
   int ilc;
   for (ilc = 0;ilc<numLC;ilc++) {
     numNeuman      = numNmLC[ilc];
     nbc            = nbcLC[ilc];
     
     if (numNeuman) { 
       gradnbcLC[ilc].reset(new BCond[numNeuman]);
       gradnbc      = gradnbcLC[ilc].get();

       for (i=0;i<numNeuman;i++) {
         gradnbc[i].nnum   = nbc[i].nnum;
         gradnbc[i].dofnum = nbc[i].dofnum;
         gradnbc[i].val    = 0; 
       }
     }

     if (numComplexNeuman) { 
       gradcnbcLC[ilc].reset(new ComplexBCond[numComplexNeuman]);
       gradcnbc        = gradcnbcLC[ilc].get();

       for (i=0;i<numComplexNeuman;i++) {
         gradcnbc[i].nnum   = cnbc[i].nnum;
         gradcnbc[i].dofnum = cnbc[i].dofnum;
         gradcnbc[i].reval   = 0; 
         gradcnbc[i].imval   = 0; 
       }
     }

   }

   // .... variable material properties   
   static_cast<GeoSource_opt*>(geoSource)->buildOptGradProp();
   
   // set flag that all arrays for analytical SA exist

   optgradFlag=1;
}

//------------------------------------------------------------------------------

OptActInfo* Domain_opt::buildOptInf()
{

  int iele;
  int vtyp,vtypset = 0;
  int eltyp,fltyp,thtyp;
    
  actInf = new OptActInfo;
  
  for(iele=0; iele<numele; ++iele) {
    vtyp = dynamic_cast<Element_opt*>(packedEset[iele])->chkOptInf(*gradnodes);
    if (vtyp)  {
      actInf->add(iele,vtyp);
      vtypset = 1;
    }
  } 

#ifdef AEROELASTIC

  if (!flExchanger) {
    fprintf(stderr,"ERROR: flExchanger does not exist in buildOptInf\n");
  }
  else {

    // check for fluid attributes
    fltyp = flExchanger->chkOptInfFL();

    // check for electrostatic attributes
    eltyp = structopt->chkOptInfEL();

    // check for thermal attributes
    thtyp = structopt->chkOptInfTH();
  
    int mptyp = fltyp + eltyp + thtyp;

    // add appropriate variable flag
 
    if (eltyp && thtyp) 
      actInf->add(-1,Optvar::mpattr);
    else if (eltyp && 
      (vtypset == Optvar::attribute || vtypset == Optvar::composite) ) 
      actInf->add(-1,Optvar::mpattr);
    else if (thtyp && 
      (vtypset == Optvar::attribute || vtypset == Optvar::composite) ) 
      actInf->add(-1,Optvar::mpattr);
    else if (fltyp)
      actInf->add(-1,fltyp);
    else if (eltyp) 
      actInf->add(-1,eltyp);
    else if (thtyp)
      actInf->add(-1,thtyp);
    else if ( mptyp != 0 ) {
      fprintf(stderr,"\n Error in Domain_opt::buildOptInf\n");
      fprintf(stderr,"   Combination of STCVAR is not supported. \n");
      exit(-1);
    }

  }
#endif // AEROELASTIC  
  return actInf;
}

//------------------------------------------------------------------------------

int Domain_opt::getOptvarTyp()
{ 
  return actInf->getTyp(); 
}

//------------------------------------------------------------------------------

void
Domain_opt::zeroGrad()
{
   int i;
   for (i=0;i<numnodes;i++) {
      (*gradnodes)[i]->x = 0.0 ; 
      (*gradnodes)[i]->y = 0.0 ;     
      (*gradnodes)[i]->z = 0.0 ;
   } 

   int ilc;
   for (ilc = 0;ilc<numLC;ilc++) {
     numNeuman = numNmLC[ilc];
     nbc       = nbcLC[ilc];
     gradnbc   = gradnbcLC[ilc].get();
     gradcnbc  = gradcnbcLC[ilc].get();
     for (i=0;i<numNeuman;i++) {
       gradnbc[i].nnum   = nbc[i].nnum;
       gradnbc[i].dofnum = nbc[i].dofnum;
       gradnbc[i].val    = 0; 
     }
     for (i=0;i<numComplexNeuman;i++) {
       gradcnbc[i].nnum   = cnbc[i].nnum;
       gradcnbc[i].dofnum = cnbc[i].dofnum;
       gradcnbc[i].reval  = 0; 
       gradcnbc[i].imval  = 0; 
     }
   }

   static_cast<GeoSource_opt*>(geoSource)->zeroGradProp();

#ifdef AEROELASTIC
   flExchanger->zerogradFluidAttr();
#endif
}   


//------------------------------------------------------------------------------
void
Domain_opt::setGradActiveLC(int lc)
{
  // watch: lc is global load case number
  if (gradnbcLC) gradnbc = gradnbcLC[lc].get();
  if (gradcnbcLC) gradcnbc = gradcnbcLC[lc].get();
} 

//------------------------------------------------------------------------------
#ifdef AEROELASTIC
//------------------------------------------------------------------------------

void Domain_opt::addGradExtConstPseudoLoad(Vector &constPseudo)
{

  int neav = 0;

  if (structopt) neav = structopt->numCurElecAttrVar;
  if (structrel) neav = structrel->numCurElecAttrVar;

  if (neav) flExchanger->addGradExtConstPseudoLoad(nodes,constPseudo);

}

//------------------------------------------------------------------------------
#endif
//------------------------------------------------------------------------------

/*
void
Domain_opt::computeTimePseudoLoad(Vector & transPseudo, int tIndex, double t,
                              int saType, Vector * gradAf)
{

  // add derivative of external domain (electrostatic, meshmotion, fluid)

#ifdef AEROELASTIC
  
     if (sinfo.aeroFlag >= 0 && sinfo.thermoeFlag < 0 ) { 

       double tFluid = flExchanger->getGradExtDomainLoad(nodes,transPseudo,
                                                         tIndex,t,saType);

       // ... save derivative of fluid load
       if ( gradAf ) (*gradAf) = transPseudo;
     }

  if (sinfo.thermoeFlag >=0){
     if (saType == 1){ //direct
        flExchanger->getGradStrucTemp(gradtemprcvd);
        fprintf(stderr," ... GradTempLoad done receiving [S]\n");
     
        //fprintf(stderr,"[S] .. the grad nodal temperatures are \n");
        //int i;
        //for(i=0;i<numnodes;i++)
        //   fprintf(stderr," %e \n",gradtemprcvd[i]);
        buildDThermalForceDs(transPseudo, 0,1);
     }
     else if(saType == 2){ //adjoint

        double* full_load = new double[6*numnodes];
        int i;
	for (i=0;i<6*numnodes;i++) full_load[i] = 0.0;
        flExchanger->getTempElecAdjLoad(full_load, numnodes,sinfo.econdThermoeFlag);
     
        //now break this info down in to active dofs and add to transPseudo
        int cn;
        for (i=0;i<numnodes;i++){     
           cn = c_dsa->locate(i,DofSet::Xdisp);
              if (cn >=0) transPseudo[cn] += full_load[i*6   ];
           cn = c_dsa->locate(i,DofSet::Ydisp);
              if (cn >=0) transPseudo[cn] += full_load[i*6 + 1];
           cn = c_dsa->locate(i,DofSet::Zdisp);
              if (cn >=0) transPseudo[cn] += full_load[i*6 + 2];
           cn = c_dsa->locate(i,DofSet::Xrot);
              if (cn >=0) transPseudo[cn] += full_load[i*6 + 3];
           cn = c_dsa->locate(i,DofSet::Yrot);
              if (cn >=0) transPseudo[cn] += full_load[i*6 + 4];
           cn = c_dsa->locate(i,DofSet::Zrot);
              if (cn >=0) transPseudo[cn] += full_load[i*6 + 5];
        }
        delete [] full_load;
     }
  }
#endif
}
*/

//------------------------------------------------------------------------------

void
Domain_opt::getMidPointMass(double& totmas, double* midPoint)
{

  totmas = 0;

  midPoint[0] = 0;  midPoint[1] = 0;  midPoint[2] = 0;

  int numActelm = actInf->size();
  
  int iele,actelm;
  for (actelm=0; actelm < numActelm; ++actelm) {

    iele = actInf->getEntry(actelm);

    totmas      += packedEset[iele]->getMass(nodes);
    double *xyz  = packedEset[iele]->getMidPoint(nodes);
   
    midPoint[0] += xyz[0];
    midPoint[1] += xyz[1];
    midPoint[2] += xyz[2];
    
    delete [] xyz;

  }

  totmas = totmas / static_cast<double>(numActelm); 

  midPoint[0] = midPoint[0]/ static_cast<double>(numActelm);
  midPoint[1] = midPoint[1]/ static_cast<double>(numActelm);
  midPoint[2] = midPoint[2]/ static_cast<double>(numActelm);

}


//------------------------------------------------------------------------------
void
Domain_opt::postProcessing(ComplexVector &sol, DComplex *bcx, ComplexVector& force, 
			   double time)
{
  OutputInfo *oinfo = geoSource->getOutputInfo();
  int numOutInfo = geoSource->getNumOutInfo();
  
  //double oldtime=sinfo.dt;
  //sinfo.dt=time;
  
   int exactNumNodes  = 0;
   int lastNodeNumber = nodes.size();

   int *nodeTable = new int[lastNodeNumber];

   int inode;
   for(inode=0; inode<lastNodeNumber; ++inode)
    nodeTable[inode] = -1;

   for(inode=0; inode<lastNodeNumber; ++inode) {
     if(nodes[inode] == 0) continue;
     exactNumNodes += 1;
     nodeTable[inode] = exactNumNodes;
   }
   
   Node oldnod;

 // Wext = external energy
 // Wela = elastic energy
 // Wkin = kinetic energy
 // Wdmp = damping energy
 // Total Energy = Wext+Wela+Wkin+Wdmp 

 double Wext=0.0,Wela=0.0,Wkin=0.0,Wdmp=0.0;
 double zero = 0.0;

 int idLC = actLC;

  // organize displacements

  int i;

  assert(dynamic_cast<Domain_opt*>(domain) != 0);
  DComplex (*xyz)[11] = new DComplex[dynamic_cast<Domain_opt*>(domain)->numnodes][11];

  for (i = 0; i < numnodes; ++i) {
    xyz[i][0] = xyz[i][1] = xyz[i][2] = xyz[i][3] = xyz[i][4] = xyz[i][5] = 0.0;
    xyz[i][6] = xyz[i][7] = xyz[i][8] = xyz[i][0] = xyz[i][10] = 0.0;
  }

  mergeDistributedDisp(xyz, sol.data(), bcx); 

  // loop over all output data

 for(i=0; i<numOutInfo; ++i) 
 {
   int w = oinfo[i].width;
   int p = oinfo[i].precision;

   // check if output file is assigned to current loadcase
   if( oinfo[i].loadcase != idLC && oinfo[i].loadcase > -1 ) 
     { continue; }

   if(oinfo[i].interval == 1) {
     switch(oinfo[i].type)
     {
       case OutputInfo::Displacement:
     	 if (oinfo[i].nodeNumber == -1) {

           DComplex x,y,z;
           fprintf(oinfo[i].filptr,"  %f\n",time);

           for(inode=0; inode<lastNodeNumber; ++inode) {

             if(nodes[inode] == 0) continue;

             int xloc  = c_dsa->locate( inode, DofSet::Xdisp);
             int xloc1 =   dsa->locate( inode, DofSet::Xdisp);

	     if(xloc >= 0) 
               x = sol[xloc];
 	     else if(xloc1 >= 0)
               x = bcx[xloc1];
             else
               x = 0.0;

             int yloc  = c_dsa->locate( inode, DofSet::Ydisp);
             int yloc1 =   dsa->locate( inode, DofSet::Ydisp);

	     if(yloc >= 0) 
               y = sol[yloc];
 	     else if(xloc1 >= 0)
               y = bcx[yloc1];
             else
               y = 0.0;

             int zloc  = c_dsa->locate( inode, DofSet::Zdisp);
             int zloc1 =   dsa->locate( inode, DofSet::Zdisp);

	     if(zloc >= 0) 
               z = sol[zloc];
 	     else if(zloc1 >= 0)
               z = bcx[zloc1];
             else
               z = 0.0;

             // ... considering the variation of shape
   	     oldnod=nodescopy->getNode(inode);
             double vx = nodes[inode]->x - oldnod.x ;
             double vy = nodes[inode]->y - oldnod.y ;
             double vz = nodes[inode]->z - oldnod.z ;
	    
	     x=x+vx;   y=y+vy;    z=z+vz;

             fprintf(oinfo[i].filptr,"% *.*E\t% *.*E\t% *.*E\n"
		     ,w,p,ScalarTypes::norm(x),
		     w,p,ScalarTypes::norm(y),
		     w,p,ScalarTypes::norm(z));
           }
           fflush(oinfo[i].filptr);
         }
       	 else { 
           int iNode = oinfo[i].nodeNumber;
           fprintf(stderr," .... WATCH: node variation not conidered in output\n");
            geoSource->outputNodeVectors(i, xyz+(7*iNode), 1, zero);
          }
          break;
       case OutputInfo::Disp6DOF:
          if (oinfo[i].nodeNumber == -1) {
            fprintf(oinfo[i].filptr,"  %f\n",time);
            for(inode=0; inode<lastNodeNumber; ++inode) {
   
               if(nodes[inode] == 0) continue;
   
               int xloc  = c_dsa->locate( inode, DofSet::Xdisp);
               int xloc1 =   dsa->locate( inode, DofSet::Xdisp);
   
               DComplex x,y,z,xr,yr,zr;
   
               if(xloc >= 0)	// dof exists and is free
                 x = sol[xloc];
               else if(xloc1 >= 0)	// dof exists and is constrained
                 x = bcx[xloc1];
               else	        // dof does not exist
                 x = 0.0;
   
               int yloc  = c_dsa->locate( inode, DofSet::Ydisp);
               int yloc1 =   dsa->locate( inode, DofSet::Ydisp);
   
               if(yloc >= 0)
                 y = sol[yloc];
               else if(xloc1 >= 0)
                 y = bcx[yloc1];
               else
                 y = 0.0;
   
               int zloc  = c_dsa->locate( inode, DofSet::Zdisp);
               int zloc1 =   dsa->locate( inode, DofSet::Zdisp);
   
   	       if(zloc >= 0) 
                 z = sol[zloc];
    	       else if(zloc1 >= 0)
                 z = bcx[zloc1];
               else
                 z = 0.0;

               int xrot  = c_dsa->locate( inode, DofSet::Xrot);
               int xrot1 =   dsa->locate( inode, DofSet::Xrot);
   
               if(xrot >= 0)
                 xr = sol[xrot];
               else if(xrot1 >= 0)
                 xr = bcx[xrot1];
               else
                 xr = 0.0;
   
               int yrot  = c_dsa->locate( inode, DofSet::Yrot);
               int yrot1 =   dsa->locate( inode, DofSet::Yrot);
   
               if(yrot >= 0)
                 yr = sol[yrot];
               else if(yrot1 >= 0)
                 yr = bcx[yrot1];
               else
                 yr = 0.0;
   
               int zrot  = c_dsa->locate( inode, DofSet::Zrot);
               int zrot1 =   dsa->locate( inode, DofSet::Zrot);
   
               if(zrot >= 0)
                 zr = sol[zrot];
               else if(zrot1 >= 0)
                 zr = bcx[zrot1];
               else
                 zr = 0.0;
   
               // ... considering the variation of shape
      	       oldnod=nodescopy->getNode(inode);
               double vx = nodes[inode]->x - oldnod.x ;
               double vy = nodes[inode]->y - oldnod.y ;
               double vz = nodes[inode]->z - oldnod.z ;
	    
   	       x=x+vx;   y=y+vy;    z=z+vz;

               fprintf(oinfo[i].filptr,
                 "%d % *.*E\t% *.*E\t% *.*E % *.*E\t % *.*E\t % *.*E\n"
		       ,nodeTable[inode], w,p,ScalarTypes::norm(x),
		       w,p,ScalarTypes::norm(y),w,p,ScalarTypes::norm(z),
		       w,p,ScalarTypes::norm(xr),w,p,ScalarTypes::norm(yr),
		       w,p,ScalarTypes::norm(zr));
            }
            fflush(oinfo[i].filptr);
          }
          else  {
            int iNode = oinfo[i].nodeNumber;
            geoSource->outputNodeVectors6(i, xyz+(7*iNode), 1, zero);
          }
          break;
       case OutputInfo::Temperature:
         fprintf(oinfo[i].filptr,"  %f\n",time);
          if(oinfo[i].nodeNumber == -1) {
            geoSource->outputNodeScalars(i, xyz[0], 0, zero);
            int inode;
            for (inode = 0; inode < numnodes; ++inode)
              geoSource->outputNodeScalars(i, xyz[inode]+6, 1);
          } 
  	  else  {
            // Only one node was requested for output
            int iNode = oinfo[i].nodeNumber;
 	    geoSource->outputNodeScalars(i, xyz[iNode]+6, 1, zero);
          }
          break;
       case OutputInfo::StressXX:
         getStressStrain(sol,bcx,i,SXX);
         break;
       case OutputInfo::StressYY:
         getStressStrain(sol,bcx,i,SYY);
         break;
       case OutputInfo::StressZZ:
         getStressStrain(sol,bcx,i,SZZ);
         break;
       case OutputInfo::StressXY:
         getStressStrain(sol,bcx,i,SXY);
         break;
       case OutputInfo::StressYZ:
         getStressStrain(sol,bcx,i,SYZ);
         break;
       case OutputInfo::StressXZ:
         getStressStrain(sol,bcx,i,SXZ);
         break;
       case OutputInfo::StrainXX:
         getStressStrain(sol,bcx,i,EXX);
         break;
       case OutputInfo::StrainYY:
         getStressStrain(sol,bcx,i,EYY);
         break;
       case OutputInfo::StrainZZ:
         getStressStrain(sol,bcx,i,EZZ);
         break;
       case OutputInfo::StrainXY:
         getStressStrain(sol,bcx,i,EXY);
         break;
       case OutputInfo::StrainYZ:
         getStressStrain(sol,bcx,i,EYZ);
         break;
       case OutputInfo::StrainXZ:
         getStressStrain(sol,bcx,i,EXZ);
         break;
       case OutputInfo::StressVM:
         getStressStrain(sol,bcx,i,VON);
         break;
       case OutputInfo::StressPR1:
         getPrincipalStress(sol,bcx,i,PSTRESS1);
         break;
       case OutputInfo::StressPR2:
         getPrincipalStress(sol,bcx,i,PSTRESS2);
         break;
       case OutputInfo::StressPR3:
         getPrincipalStress(sol,bcx,i,PSTRESS3);
         break;
       case OutputInfo::StrainPR1:
         getPrincipalStress(sol,bcx,i,PSTRAIN1);
         break;
       case OutputInfo::StrainPR2:
         getPrincipalStress(sol,bcx,i,PSTRAIN2);
         break;
       case OutputInfo::StrainPR3:
         getPrincipalStress(sol,bcx,i,PSTRAIN3);
         break;
       case OutputInfo::InXForce:
         getElementForces(sol, bcx, i, INX);
         break;
       case OutputInfo::InYForce:
         getElementForces(sol, bcx, i, INY);
         break;
       case OutputInfo::InZForce:
         getElementForces(sol, bcx, i, INZ);
         break;
       case OutputInfo::AXMoment:
         getElementForces(sol, bcx, i, AXM);
         break;
       case OutputInfo::AYMoment:
         getElementForces(sol, bcx, i, AYM);
         break;
       case OutputInfo::AZMoment:
         getElementForces(sol, bcx, i, AZM);
         break;
	 /*
       case OutputInfo::Energies: {
          Wext = force *  sol;   // Wext = external energy
          Wela =   0.5 * Wext;   // Wela = elastic energy 
          double error = Wext+Wela+Wkin+Wdmp;
	  assert(static_cast<GeoSource_opt*>(geoSource) != 0);
          static_cast<GeoSource_opt*>(geoSource)->outputEnergies(i, zero, Wext, Wela, Wkin, Wdmp, error); }
         break;
	 */
       case OutputInfo::StrainVM:
	 getStressStrain(sol,bcx,i,STRAINVON);
         break;
       case OutputInfo::YModulus:
         getElementAttr(i,YOUNG, time);
         break;
       case OutputInfo::MDensity:
         getElementAttr(i,MDENS, time);
         break;
       case OutputInfo::Thicknes:
         getElementAttr(i,THICK, time);
         break;
       case OutputInfo::Composit:
         getCompositeData(i,time);
         break;
	 /*
       case OutputInfo::Ctexp:
         getElementAttr(i,CTEXP);
         break;
       case OutputInfo::Helmholtz:
         // ... PRINT (REAL) HELMHOLTZ SOLUTION
         for(inode=0; inode<numnodes; ++inode) {
              int loc  = c_dsa->locate( inode, DofSet::Helm);
              int loc1 =   dsa->locate( inode, DofSet::Helm);

              double xHelm;
              if(loc >= 0)        // dof exists and is free
                xHelm = sol[loc];
              else if(loc1 >= 0)  // dof exists and is constrained
                xHelm = bcx[loc1];
              else                // dof does not exist
                xHelm = 0.0;

              fprintf(oinfo[i].filptr,"% *.*E\n",w,p,xHelm);
         }
         break;
	 */
       case OutputInfo::ShapeAtt:
         fprintf(oinfo[i].filptr,"  %f\n",time);
	 
         for(inode=0; inode<lastNodeNumber; ++inode) {

           if(nodes[inode] == 0) continue;

	   oldnod=nodescopy->getNode(inode);
           double x = nodes[inode]->x - oldnod.x ;
           double y = nodes[inode]->y - oldnod.y ;
           double z = nodes[inode]->z - oldnod.z ;

           fprintf(oinfo[i].filptr,"% *.*E\t% *.*E\t% *.*E\n"
                               ,w,p,x,w,p,y,w,p,z);
	 }
	 break;
       case OutputInfo::ShapeStc:
         fprintf(oinfo[i].filptr,"  %f\n",time);
	 
         for(inode=0; inode<lastNodeNumber; ++inode) {

           if(nodes[inode] == 0) continue;

 	   oldnod=nodescopy->getNode(inode);
           double x = nodes[inode]->x - oldnod.x ;
           double y = nodes[inode]->y - oldnod.y ;
           double z = nodes[inode]->z - oldnod.z ;

           fprintf(oinfo[i].filptr,"% *.*E\t% *.*E\t% *.*E\n"
                               ,w,p,x,w,p,y,w,p,z);
	 }
	 break;
     default:
       fprintf(stderr, "---------------> %d\n", oinfo[i].type);
       assert(0);       
       break;
     }
    }
  } 

 //sinfo.dt=oldtime;

 // delete temporary arrays
 delete [] nodeTable;
 delete [] xyz;

}

#endif
