#ifdef STRUCTOPT

#include <Structopt.d/Structopt_base.h>
#include <Structopt.d/Optcrit.h>

// global singletons
std::auto_ptr<Structopt> structopt;
std::auto_ptr<Structopt> structrel;

std::auto_ptr<Optpro> optpro;
std::auto_ptr<Optpro> relpro;

Structopt::~Structopt() {}

//------------------------------------------------------------------------------
double* Structopt::getptreleattr( int& loc1, int& loc2 ) {

  double *p;
	  
  SPropContainer& sprop = geoSource->getStructProps();
	   	      
  switch (loc2) {
  case 0: 
    p=&(sprop[loc1].A); 
    break;
  case 1: 
    p=&(sprop[loc1].E); 
    break;
  case 2: 
    p=&(sprop[loc1].nu); 
    break;
  case 3: 
    p=&(sprop[loc1].rho); 
    break;
  case 5: 
    p=&(sprop[loc1].k); 
    break;
  case 6: 
    p=&(sprop[loc1].eh); 
    break;
  case 10: 
    p=&(sprop[loc1].W); 
    break;
  default: 
    p=0; 
  }
  return p;
}

//------------------------------------------------------------------------------
double* Structopt::getptrgradeleattr( int& loc1, int& loc2 ) {

  double *p;
	  
  SPropContainer& gradsprop = static_cast<GeoSource_opt*>(geoSource)->getGradStructProps();
	   	      
  switch (loc2) {
  case 0: 
    p=&(gradsprop[loc1].A); 
    break;
  case 1: 
    p=&(gradsprop[loc1].E); 
    break;
  case 2: 
    p=&(gradsprop[loc1].nu); 
    break;
  case 3: 
    p=&(gradsprop[loc1].rho); 
    break;
  case 5: 
    p=&(gradsprop[loc1].k); 
    break;
  case 6: 
    p=&(gradsprop[loc1].eh); 
    break;
  case 10: 
    p=&(gradsprop[loc1].W); 
    break;
  default: 
    p=0; 
  }
  return p;
}

//------------------------------------------------------------------------------
double* Structopt::getptrcomposite(int& loc1, int& loc2)
{
  double* p=0;
  if(loc2 >= 0)
    {
      // layers
      int layer = loc2/100000;
      int typ   = loc2 - 100000*layer;
      
      LayInfo* layinfo = geoSource->getLayerInfo(loc1);
      
      int numlayer = layinfo->nLayers();
      
      if (numlayer < layer) 
	{
	  fprintf(stderr," *** ERROR: Layer %d to optimize does not exist\n",
		  layer);
	  exit(-1);
	}
      
      switch (typ) 
	{      
	case 1: 
	  p=&(layinfo->data[layer-1][0]); 
	  break;
	case 2: 
	  p=&(layinfo->data[layer-1][1]); 
	  break;
	case 4: 
	  p=&(layinfo->data[layer-1][3]); 
	  break;
	case 7: 
	  p=&(layinfo->data[layer-1][6]); 
	  break;
	case 8: 
	  p=&(layinfo->data[layer-1][7]); 
	  break;
	case 9: 
	  p=&(layinfo->data[layer-1][8]); 
	  break;
	default: 
	  p=0; 
	}
    }
  else
    {
      // coefficients
      int i = (-loc2-1)/6;
      int j = (-loc2-1)%6;
      assert(0<=i && i<6 && 0<=j && j<6);
      CoefData* pCfDta = geoSource->getCoefData(loc1);
      p = pCfDta->values() + (std::min(i,j)+6*std::max(i,j));
    }
  
  return p;
}

//------------------------------------------------------------------------------
double * Structopt::getptrgradcomposite( int& loc1, int& loc2 ) 
{
  double* p = 0;
  if(loc2 >= 0)
    {
      // layers
      int layer = loc2/100000;
      int typ   = loc2 - 100000*layer;
      
      LayInfo * layinfo = geoSource->getLayerInfo(loc1);      
      switch (typ) 
	{	   	   
	case 1: 
	  p=&(layinfo->grad[layer-1][0]); 
	  break;
	case 2: 
	  p=&(layinfo->grad[layer-1][1]); 
	  break;
	case 4: 
	  p=&(layinfo->grad[layer-1][3]); 
	  break;
	case 7: 
	  p=&(layinfo->grad[layer-1][6]); 
	  break;
	case 8: 
	  p=&(layinfo->grad[layer-1][7]); 
	  break;
	case 9: 
	  p=&(layinfo->grad[layer-1][8]); 
	  break;
	default: 
	  p=0; 
	}
    }
  else
    {
      // coefficients
      int i = (-loc2-1)/6;
      int j = (-loc2-1)%6;
      assert(0<=i && i<6 && 0<=j && j<6);
      CoefData* pCfDta = static_cast<GeoSource_opt*>(geoSource)->getGradCoefData(loc1);
      p = pCfDta->values() + (std::min(i,j)+6*std::max(i,j));
    }
  return p;
}


//------------------------------------------------------------------------------
void Structopt::gradcrit(int ivar, double** gc, int numAC, int* listAC )
{
  int icrit;

  double oldval;  

  if (!numAC) {
  
    for (icrit=0;icrit<optpro->numcrit;icrit++) {
  
      oldval= optpro->opc[icrit]->grad;

      optpro->opc[icrit]->gradcrit(this,ivar);

      gc[icrit][ivar] = optpro->opc[icrit]->grad;

      optpro->opc[icrit]->grad = oldval;    
    }
  }
  else {

    for (icrit=0;icrit<numAC;icrit++) {

      int ac=listAC[icrit];

      oldval= optpro->opc[ac]->grad;

      optpro->opc[ac]->gradcrit(this,ivar);

      gc[icrit][ivar] = optpro->opc[ac]->grad;

      optpro->opc[ac]->grad = oldval;
    }
  }
}

//------------------------------------------------------------------------------
double Structopt::getgradpartdisplevel(int aFlag, int quoFlag, int size, 
				       int * nodeid, int*distyp, double* refVal,                                        
				       double* difVal,double& powFac,
				       double& time,int anaId)
{
  return 0;
}


#endif
