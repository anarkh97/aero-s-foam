#ifndef _HEADER_H_
#define _HEADER_H_

static const char*header[] = {
"Vector DISP under %s for %s\n%d\n",
"Scalar TEMP under %s for %s\n%d\n",
"Scalar SIGMAXX under %s for %s\n%d\n",
"Scalar SIGMAYY under %s for %s\n%d\n",
"Scalar SIGMAZZ under %s for %s\n%d\n",
"Scalar SIGMAXY under %s for %s\n%d\n",
"Scalar SIGMAYZ under %s for %s\n%d\n",
"Scalar SIGMAXZ under %s for %s\n%d\n",
"Scalar STRAINXX under %s for %s\n%d\n",
"Scalar STRAINYY under %s for %s\n%d\n",
"Scalar STRAINZZ under %s for %s\n%d\n",
"Scalar STRAINXY under %s for %s\n%d\n",
"Scalar STRAINYZ under %s for %s\n%d\n",
"Scalar STRAINXZ under %s for %s\n%d\n",
"Scalar HEATFLXX under %s for %s\n%d\n",
"Scalar HEATFLXY under %s for %s\n%d\n",
"Scalar HEATFLXZ under %s for %s\n%d\n",
"Scalar GRDTEMPX under %s for %s\n%d\n",
"Scalar GRDTEMPY under %s for %s\n%d\n",
"Scalar GRDTEMPZ under %s for %s\n%d\n",
"Scalar VONMISES under %s for %s\n%d\n",
"Scalar STRESSP1 under %s for %s\n%d\n",
"Scalar STRESSP2 under %s for %s\n%d\n",
"Scalar STRESSP3 under %s for %s\n%d\n",
"Scalar STRAINP1 under %s for %s\n%d\n",
"Scalar STRAINP2 under %s for %s\n%d\n",
"Scalar STRAINP3 under %s for %s\n%d\n",
"Scalar INXFORCE under %s for %s\n%d\n",
"Scalar INYFORCE under %s for %s\n%d\n",
"Scalar INZFORCE under %s for %s\n%d\n",
"Scalar AXMOMENT under %s for %s\n%d\n",
"Scalar AYMOMENT under %s for %s\n%d\n",
"Scalar AZMOMENT under %s for %s\n%d\n",
"#   Time         Wext         Waer         Wela         Wkin         Wdmp         Error\n    ____         ____         ____         ____         ____         ____          ____\n\n",
"Vector AEROFORCE under %s for %s\n%d\n",
"Vector MODE under %s for %s\n%d\n",
"Scalar VONMISES_STRAIN under %s for %s\n%d\n",
"Scalar PRESSURE under %s for %s\n%d\n",//PRESSURE is for HELMHOLTZ also
"Vector DISP6DOF under %s for %s\n%d\n",
"Vector MODE6 under %s for %s\n%d\n",
"Scalar AEROFORX under %s for %s\n%d\n",
"Scalar PRESSURE under %s for %s\n%d\n",//PRESSURE is for HELMHOLTZ also
"Scalar AEROFORY under %s for %s\n%d\n",
"Scalar AEROFORZ under %s for %s\n%d\n",
"Scalar AEROMOMX under %s for %s\n%d\n",
"Scalar AEROMOMY under %s for %s\n%d\n",
"Scalar AEROMOMZ under %s for %s\n%d\n",
"Vector VELOCITY under %s for %s\n%d\n",
"Vector ACCELERATION under %s for %s\n%d\n",
"Scalar YMODULUS under %s for %s\n%d\n",
"Scalar MDENSITY under %s for %s\n%d\n",
"Scalar THICKNES under %s for %s\n%d\n",
"Vector SHAPEATT under %s for %s\n%d\n",
"Vector SHAPESTC under %s for %s\n%d\n",
"Vector COMPOSIT under %s for %s\n%d\n",
"Scalar DX under %s for %s\n%d\n",
"Scalar DY under %s for %s\n%d\n",
"Scalar DZ under %s for %s\n%d\n",
"Scalar RX under %s for %s\n%d\n",
"Scalar RY under %s for %s\n%d\n",
"Scalar RZ under %s for %s\n%d\n",
"Scalar MD under %s for %s\n%d\n",
"Scalar MR under %s for %s\n%d\n",
"Scalar MT under %s for %s\n%d\n",
"Vector MODE under %s for %s\n%d\n",
"Element to Node Table\n",
"Node to Element Table\n",
"Node to Node Table\n",
"# Time      HeatFlux n+1        HeatFlux n+1/2\n",
"Vector HEATFLX under %s for %s\n%d\n",
"Vector GRDTEMP under %s for %s\n%d\n",
"Vector VELOCITY6 under %s for %s\n%d\n",
"Vector ACCELERATION6 under %s for %s\n%d\n",
"# Time    Alpha(s)\n",
"# Time    Error\n",
"Vector REACTIONS under %s for %s\n%d\n",
"",
"",
"Scalar CONPRESS under %s for %s\n%d\n",
"",
"",
"",
"",
"Scalar PRESSURE under %s for %s\n%d\n",
"Vector MODE under %s for %s\n%d\n",
"Scalar MODE under %s for %s\n%d\n",
"Vector STRESSP1DIRECTION under %s for %s\n%d\n",
"Vector STRESSP2DIRECTION under %s for %s\n%d\n",
"Vector STRESSP3DIRECTION under %s for %s\n%d\n",
"Vector STRAINP1DIRECTION under %s for %s\n%d\n",
"Vector STRAINP2DIRECTION under %s for %s\n%d\n",
"Vector STRAINP3DIRECTION under %s for %s\n%d\n",
"Scalar POTENTIAL under %s for %s\n%d\n",
"Scalar SLOSHING_X_DISPLACEMENT under %s for %s\n%d\n",
"Scalar SLOSHING_Y_DISPLACEMENT under %s for %s\n%d\n",
"Scalar SLOSHING_Z_DISPLACEMENT under %s for %s\n%d\n",
"Scalar SLOSHING_DISPLACEMENTS under %s for %s\n%d\n",
"Scalar TDENFORCEMENT under %s for %s\n%d\n",
"Scalar DAMAGE under %s for %s\n%d\n",
"Scalar EQPLSTRN under %s for %s\n%d\n",
"Scalar TEMPERATUREVEL under %s for %s\n%d\n",
"Scalar PRESSUREVEL under %s for %s\n%d\n",
"Scalar PRESSUREACC under %s for %s\n%d\n",
"Scalar HEATREACTIONS under %s for %s\n%d\n",
"Vector REACTIONS6 under %s for %s\n%d\n",
"",
"",
"",
"",
"",
"",
"",
"Matrix ROTATIONS under %s for %s\n%d\n",
"Scalar EXTFORCEX under %s for %s\n%d\n",
"Scalar EXTFORCEY under %s for %s\n%d\n",
"Scalar EXTFORCEZ under %s for %s\n%d\n",
"Scalar EXTMOMENTX under %s for %s\n%d\n",
"Scalar EXTMOMENTY under %s for %s\n%d\n",
"Scalar EXTMOMENTZ under %s for %s\n%d\n",
"",
""
};

static const char * ele_header[] = {
"",
"",
"ElemScalar ELE_SIGMAXX under %s using %s_pattern\n",
"ElemScalar ELE_SIGMAYY under %s using %s_pattern\n",
"ElemScalar ELE_SIGMAZZ under %s using %s_pattern\n",
"ElemScalar ELE_SIGMAXY under %s using %s_pattern\n",
"ElemScalar ELE_SIGMAYZ under %s using %s_pattern\n",
"ElemScalar ELE_SIGMAXZ under %s using %s_pattern\n",
"ElemScalar ELE_STRAINXX under %s using %s_pattern\n",
"ElemScalar ELE_STRAINYY under %s using %s_pattern\n",
"ElemScalar ELE_STRAINZZ under %s using %s_pattern\n",
"ElemScalar ELE_STRAINXY under %s using %s_pattern\n",
"ElemScalar ELE_STRAINYZ under %s using %s_pattern\n",
"ElemScalar ELE_STRAINXZ under %s using %s_pattern\n",
"ElemScalar ELE_HEATFLXX under %s using %s_pattern\n",
"ElemScalar ELE_HEATFLXY under %s using %s_pattern\n",
"ElemScalar ELE_HEATFLXZ under %s using %s_pattern\n",
"ElemScalar ELE_GRDTEMPX under %s using %s_pattern\n",
"ElemScalar ELE_GRDTEMPY under %s using %s_pattern\n",
"ElemScalar ELE_GRDTEMPZ under %s using %s_pattern\n",
"ElemScalar ELE_VONMISES under %s using %s_pattern\n",
"ElemScalar ELE_STRESSP1 under %s using %s_pattern\n",
"ElemScalar ELE_STRESSP2 under %s using %s_pattern\n",
"ElemScalar ELE_STRESSP3 under %s using %s_pattern\n",
"ElemScalar ELE_STRAINP1 under %s using %s_pattern\n",
"ElemScalar ELE_STRAINP2 under %s using %s_pattern\n",
"ElemScalar ELE_STRAINP3 under %s using %s_pattern\n",
"ElemScalar ELE_INXFORCE under %s using %s_pattern\n",
"ElemScalar ELE_INYFORCE under %s using %s_pattern\n",
"ElemScalar ELE_INZFORCE under %s using %s_pattern\n",
"ElemScalar ELE_AXMOMENT under %s using %s_pattern\n",
"ElemScalar ELE_AYMOMENT under %s using %s_pattern\n",
"ElemScalar ELE_AZMOMENT under %s using %s_pattern\n",
"",
"",
"",
"ElemScalar ELE_VONMISES_STRAIN under %s using %s_pattern\n",
"",
"",
"ElemScalar ELE_AEROFORX under %s using %s_pattern\n",
"ElemScalar ELE_AEROFORY under %s using %s_pattern\n",
"ElemScalar ELE_AEROFORZ under %s using %s_pattern\n",
"ElemScalar ELE_AEROMOMX under %s using %s_pattern\n",
"ElemScalar ELE_AEROMOMY under %s using %s_pattern\n",
"ElemScalar ELE_AEROMOMZ under %s using %s_pattern\n",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"ElemScalar ELE_DAMAGE under %s using %s_pattern\n"
};

#endif
