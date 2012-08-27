#include <Structopt.d/Driver_opt.d/GeoSource_opt.h>
#include <Structopt.d/Element_opt.d/Element_opt.h>

#include <cassert>

//-----------------------------------------------------------------------
#ifdef STRUCTOPT

void GeoSource_opt::buildOptGradProp()
{
   if (!numProps) {
     fprintf(stderr," ERROR: need at least one material in optimization\n");
     exit(-1);
   }

   for (int i=0; i<na; ++i) 
     {
       Element_opt *ele = dynamic_cast<Element_opt*>(elemSet[ attrib[i].nele ]);
       assert(ele != 0);
       if (attrib[i].attr >= 0)
	 { ele->setGradProp(&(sgradprops[attrib[i].attr])); }
     }

   // .... variable composite properties
   for (int i=0;i<numLayInfo;i++) {
     layInfo[i]->setGrad();
     layInfo[i]->zeroGrad();
   }

   for (int i=0; i < na; ++i) 
     {
       Element_opt *ele = dynamic_cast<Element_opt*>(elemSet[ attrib[i].nele ]);
       assert(ele != 0);
       if (attrib[i].cmp_attr >= 0) 
	 {
	   if(coefData[attrib[i].cmp_attr] != 0) 	     
	     { ele->setDCCoefs(getGradCoefData(attrib[i].cmp_attr)->values()); }
	   if(layInfo[attrib[i].cmp_attr] != 0)
	     { ele->setCompositeGrad(layInfo[attrib[i].cmp_attr]->gradval()); }
	 }
     }


   // Set up a derivative for each existing element frame

   for (int iFrame = 0; iFrame < numEframes; iFrame++)  
     {
       Element_opt *ele = dynamic_cast<Element_opt*>(elemSet[efd[iFrame].elnum]);
       assert(ele != 0);
       if (ele) 
	 {
	   double* dframe = new double[9];
	   for(int i =0; i < 9; ++i) { dframe[i] = 0; }
	   ele->setdFrame((EFrame*)dframe);
	 }
     }
}

//-----------------------------------------------------------------------

void GeoSource_opt::zeroGradProp()
{
   int i;
   for (i=0;i<numProps;i++) {
      sgradprops[i].E         = 0.0; 
      sgradprops[i].A         = 0.0;
      sgradprops[i].nu        = 0.0;
      sgradprops[i].rho       = 0.0;
      sgradprops[i].eh        = 0.0;
      sgradprops[i].Ixx       = 0.0;
      sgradprops[i].Iyy       = 0.0;
      sgradprops[i].Izz       = 0.0;
      sgradprops[i].c         = 0.0;
      sgradprops[i].k         = 0.0;
      sgradprops[i].Q         = 0.0;
      sgradprops[i].W         = 0.0;
      sgradprops[i].P         = 0.0;
      sgradprops[i].Ta        = 0.0;
      sgradprops[i].kappaHelm = 0.0;
   }

   for (i=0;i<numLayInfo;i++) 
     { layInfo[i]->zeroGrad(); }   

   for (int iFrame = 0; iFrame < numEframes; iFrame++)  
     {
       Element_opt *ele = dynamic_cast<Element_opt*>(elemSet[efd[iFrame].elnum]);       
       assert(ele != 0);
       if (!ele) 
	 { ele->zerodFrame(); }
     }
}

void GeoSource_opt::outputEnergies(int fileNum, double time, double Wext, 
			       double Wela, double Wkin, double Wdmp,
                               double error)  {


  fprintf(oinfo[fileNum].filptr,"%20.12e %20.12e %20.12e %20.12e %20.12e %20.12e\n",time,Wext,Wela,Wkin,Wdmp,error);

  fflush(oinfo[fileNum].filptr);
}


CoefData* GeoSource_opt::getGradCoefData(int i)
{
  assert(i>=0 && i<numCoefData);
  CoefData* ret = 0;
  std::map<int, CoefData>::iterator lb = d_coefData.lower_bound(i);
  if(lb != d_coefData.end() && lb->first == i)
    { ret = &(lb->second); }
  else
    {
      CoefData n; n.zero();
      std::map<int, CoefData>::iterator pos = d_coefData.insert(lb, std::make_pair(i,n));
      assert(pos->first == i);
      ret = &(pos->second);
    }
  return ret;
}



//-----------------------------------------------------------------------
#endif // STRUCTOPT
//-----------------------------------------------------------------------
