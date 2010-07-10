%x incl

%{
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <Parser.d/AuxDefs.h>
#include <Driver.d/Domain.h>
#include <Parser.d/p.tab.H>
#include <stack>
//#include <Mortar.d/FaceElement.d/SurfaceEntity.h>
int line_num = 1;

// used for include
stack<YY_BUFFER_STATE> includedFiles;
stack<int> lineNumber;	// for fancy error message

%}

%e 40000
%p 20000
%a 40000
%n 20000
%k 40000
%o 30000


DIGIT [0-9]
HEXDIGIT [0-9A-Fa-f]
SPACES [ \t]+
NAME [a-zA-Z_][a-zA-Z_0-9]*
TRAIL [a-zA-Z]*
%%
"INCLUDE"|"include" 		{ BEGIN(incl); 	} // go to include state 
<incl>[ \t]* 
<incl>[^ \t\n]+		{
				// clean up !
				int c = yyinput();
				while(c == '\t' || c == ' ')
					c = yyinput();
				if(c != '\n')
					yyunput(c,yytext_ptr);
				string fname = string(yytext);
				if(fname.length() > 0)
				{
					if(fname[0]=='"' && fname[ fname.length()-1 ]=='"' && fname.length()> 2)
						fname = fname.substr(1,fname.length()-2);
						
					FILE * f = fopen(fname.c_str(),"r");
					if(f!=NULL)
					{	// open the file and substitute it as the current buffer
						//std::cout << " ... Including file " << fname << endl;
						includedFiles.push(YY_CURRENT_BUFFER);
						lineNumber.push(line_num);
						line_num = 1;
						yy_switch_to_buffer( yy_create_buffer(f,YY_BUF_SIZE) );
					}
					else
					{
						std::cout << " **************************************************** "<< endl 
							  << " ** ERROR :  Could'nt open file " << fname << endl
							  << " **************************************************** "<< endl 
;
					}
					BEGIN(INITIAL); // go back to parsing
				}	
		    }
<<EOF>>	    {
			if(includedFiles.size()==0)
			{	// end of everything
				yyterminate();
			}
			else
			{	// switch to previous buffer
		    	   yy_delete_buffer( YY_CURRENT_BUFFER );
            		   yy_switch_to_buffer( includedFiles.top());
			   line_num = lineNumber.top()+1;
			   lineNumber.pop();
			   includedFiles.pop();
			}
	    }

"NODE"[a-zA-Z]*|"node"[a-zA-Z]*    { return NODE;    }

"ZERO"[a-zA-Z]*|"zero"[a-zA-Z]*    { return ZERO;    }

"READMODE"[a-zA-Z]*|"readmode"[a-zA-Z]* { return READMODE; }

"PRES"[a-zA-Z]*|"pres"[a-zA-Z]*    { return PRESSURE; }

"PREL"[a-zA-Z]*|"prel"[a-zA-Z]*    { return PRELOAD; }

"BUCK"[a-zA-Z]*|"buck"[a-zA-Z]*    { return BUCKLE; }

"CONTROL"[a-zA-Z]*|"control"[a-zA-Z]*    { return CONTROL;  }

"STRU"[a-zA-Z]*|"stru"[a-zA-Z]*    { return OPTIMIZATION;  }

"TOPO"[a-zA-Z]*|"topo"|"topology"  { return TOPOLOGY; }

"SOWE"[a-zA-Z]*|"sowe"|"topology"  { return SOWER; }

"BINARY"|"binary"                  { return BINARY; }

"BINARYINPUT"|"binaryinput"        { return BINARYINPUT; }

"BINARYOUTPUT"|"binaryoutput"      { return BINARYOUTPUT; }

"ON"|"on"|"On"                     { yylval.ival = 1; return SWITCH; }
"OFF"|"off"|"Off"                  { yylval.ival = 0; return SWITCH; }
"YES"|"yes"|"Yes"                  { yylval.ival = 1; return SWITCH; }
"NO"|"no"|"No"                     { yylval.ival = 0; return SWITCH; }
"TRUE"|"true"|"True"               { yylval.ival = 1; return SWITCH; }
"FALSE"|"false"|"False"            { yylval.ival = 0; return SWITCH; }

"DECO"[a-zA-Z]{0,7}|"deco"[a-zA-Z]{0,7} { return DECOMPOSE; }

"nsubs"|"NSUBS"                    { return NSUBS; }

"outfile"|"OUTFILE"                { return DECOMPFILE; }

"exit"|"EXIT"                      { return EXITAFTERDEC; }

"skip"|"SKIP"                      { return SKIPDECCALL; }

"outload"|"OUTLOAD"                { return OUTPUTWEIGHT; }

"outmem"|"OUTMEM"                  { return OUTPUTMEMORY; }

"deter"|"DETER"                    { return DETER; }

"WEIGHTS"|"weights"                { return WEIGHTLIST; }

"decomposition"|"DECOMPOSITION"    { return DECOMPOSITION; }

"global"|"GLOBAL"                  { return GLOBAL; }

"matcher"|"MATCHER"|"Matcher"      { return MATCHER; }

"geometry"|"GEOMETRY"              { return GEOMETRY; }

"cpumap"|"CPUMAP"                  { return CPUMAP; }

"TOPF"[a-zA-Z]*|"topf"[a-zA-Z]*    { return TOPFILE;  }

"ATTR"[a-zA-Z]*|"attr"[a-zA-Z]*    { return ATTRIBUTES; }

"THETA"|"theta"                    { return THETA; }

"IDENTITY"|"identity"              { return IDENTITY; }

"EFRA"[a-zA-Z]*|"efra"[a-zA-Z]*    { return EFRAMES; }

"THIRDNODE"|"thirdnode"            { return THIRDNODE; }

"BOFF"[a-zA-Z]*|"boff"[a-zA-Z]*    { return BOFFSET; }

"idisp6pita"|"IDISP6PITA"          { return PITADISP6; }

"ivel6pita"|"IVEL6PITA"            { return PITAVEL6; }         

"IDISP"[a-zA-Z]*|"idisp"[a-zA-Z]*  { return IDIS;  }

"IDISP6"[a-zA-Z]*|"idisp6"[a-zA-Z]* { return IDIS6; }

"GEPS"[a-zA-Z]*|"geps"[a-zA-Z]*    { return GEPS; }

"IVEL"[a-zA-Z]*|"ivel"[a-zA-Z]*    { return IVEL; }

"IACC"[a-zA-Z]*|"iacc"[a-zA-Z]*    { return IACC; }

"ITEM"[a-zA-Z]*|"item"[a-zA-Z]*    { return ITEMP; }

"TEMP"[a-zA-Z]*|"temp"[a-zA-Z]*    { return TEMP; }

"PDIR"[a-zA-Z]*|"pdir"[a-zA-Z]*    { return PDIR; }

"HEFR"[a-zA-Z]*|"hefr"[a-zA-Z]*    { return HEFRS; }

"FLUX"[a-zA-Z]*|"flux"[a-zA-Z]*    { return FLUX; } 

"END"|"end"                        { return END;  }

"LOAD"|"load"                      { return LOAD; }

"USE"|"use"                        { return USE;  }

"SENS"[a-zA-Z]*|"sens"[a-zA-Z]*    { return SENSORS; }

"ACTU"[a-zA-Z]*|"actu"[a-zA-Z]*    { return ACTUATORS; }

"MATE"[a-zA-Z]*|"mate"[a-zA-Z]*    { return MATERIALS; }

"DISP"[a-zA-Z]*|"disp"[a-zA-Z]*    { return DISP; }

"THRU"[a-zA-Z]*|"thru"[a-zA-Z]*    { return THRU; }

"STEP"[a-zA-Z]*|"step"[a-zA-Z]*    { return STEP; }

"LMPC"[a-zA-Z]*|"lmpc"[a-zA-Z]*    { return LMPC; }

"HLMPC"[a-zA-Z]*|"hlmpc"[a-zA-Z]*  { return HLMPC; }

"CRHS"|"crhs"                      { return CRHS; }

"CTYPE"|"ctype"                    { return CTYPE; }

"FORC"[a-zA-Z]*|"forc"[a-zA-Z]*    { return FORCE; }

"COND"[a-zA-Z]*|"cond"[a-zA-Z]*    { return CONDITION; }

"MASS"[a-zA-Z]*|"mass"[a-zA-Z]*    { return MASS; }

"QLOA"[a-zA-Z]*|"qloa"[a-zA-Z]*    { return QLOAD; }

"PLOA"[a-zA-Z]*|"ploa"[a-zA-Z]*    { return PLOAD; }

"STAT"[a-zA-Z]*|"stat"[a-zA-Z]*    { return STATS; }

"AXIHDIR"|"axihdir"                { return AXIHDIR; }

"AXIHNEU"|"axihneu"                { return AXIHNEU; }

"AXINUMMODES"|"axinummodes"        { return AXINUMMODES; }

"AXINUMSLICES"|"axinumslices"      { return AXINUMSLICES; }

"AXIHSOMMER"|"axihsommer"          { return AXIHSOMMER; }

"AXIMPC"|"aximpc"                  { return AXIMPC; }

"ELSCATT"[a-zA-Z]*|"elscatt"[a-zA-Z]* { return ELSCATTERER; }

"FFP"|"ffp"                        {return FFP; }

"FFPDIR"|"ffpdir"                  {return FFPDIR; }

"HSOM"[a-zA-Z]*|"hsom"[a-zA-Z]*    { return HSOMMERFELD; }

"SCATT"[a-zA-Z]*|"scatt"[a-zA-Z]*  { return SCATTERER; }

"IHDI"[a-zA-Z]*|"ihdi"[a-zA-Z]*    { return IHDIRICHLET; }

"IHNE"[a-zA-Z]*|"ihne"[a-zA-Z]*    { return IHNEUMANN; }

"IHDS"[a-zA-Z]*|"ihds"[a-zA-Z]*    { return IHDSWEEP; }

"FACO"[a-zA-Z]*|"faco"[a-zA-Z]*    { return FACOUSTICS; }

"ELHSOM"[a-zA-Z]*|"elhsom"[a-zA-Z]* { return ELHSOMMERFELD; }

"HPNEU"[a-zA-Z]*|"hpneu"[a-zA-Z]*  { return HNEUMAN; }

"TURKEL"|"turkel"                  { return TURKEL; }

"SHIFT"|"shift"                    { return SHIFT; }

"IMPE"[a-zA-Z]*|"impe"[a-zA-Z]*    { return IMPE; }

"freq"|"FREQ"|"frequency"|"FREQUENCY" { return FREQ; }

"freqsweep"|"FREQSWEEP"            { return FREQSWEEP; }

"freqsweep1"|"FREQSWEEP1"          { return FREQSWEEP1; }

"freqsweep2"|"FREQSWEEP2"          { return FREQSWEEP2; }

"RECONS"[a-zA-Z]*|"recons"[a-zA-Z]* { return RECONS; }

"Taylor"|"taylor"|"TAYLOR"         { yylval.ival = 0; return RECONSALG; }

"Pade1"|"pade1"|"PADE1"            { yylval.ival = 1; return RECONSALG; }

"Pade"|"pade"|"PADE"               { yylval.ival = 2; return RECONSALG; }

"padepivot"|"PADEPIVOT"            { return PADEPIVOT; }

"Fourier"|"fourier"|"FOURIER"      { yylval.ival = 3; return RECONSALG; }

"PadeLanczos"|"padelanczos"|"PADELANCZOS" { yylval.ival = 4; return RECONSALG; }

"GalProjection"|"galprojection"|"GALPROJECTION" { yylval.ival = 5; return RECONSALG; }

"KrylovGalProjection"|"krylovgalprojection"|"KRYLOVGALPROJECTION" { yylval.ival= 6; return RECONSALG; }

"QRGalProjection"|"qrgalprojection"|"QRGALPROJECTION" { yylval.ival = 7; return RECONSALG; }

"padepoles"|"PADEPOLES"            { return PADEPOLES; }

"modefilter"|"MODEFILTER"          { return MODEFILTER; }

"test_ulrich"                      { return TESTULRICH; }

"numcgm"|"NUMCGM"                  { return NUMCGM; }

"interfacelumped"|"INTERFACELUMPED" { return INTERFACELUMPED; }

"MODE"|"mode"                      { return MODE; }

"CG"|"cg"|"PCG"|"pcg"              { yylval.ival = 0; return ITERTYPE; }

"GMRES"|"gmres"                    { yylval.ival = 1; return ITERTYPE; }

"GCR"|"gcr"                        { yylval.ival = 2; return ITERTYPE; }

"CGAL"|"cgal"                      { yylval.ival = 3; return ITERTYPE; }

"BCG"|"bcg"                        { yylval.ival = 4; return ITERTYPE; }

"CR"|"cr"                          { yylval.ival = 5; return ITERTYPE; }

"FETIH"|"fetih"                    { return FETIH; }

"HFETI"|"hfeti"                    { return HFETI; }

"tolcgm"|"TOLCGM"                  { return TOLCGM; }

"spacedim"|"SPACEDIM"              { return SPACEDIMENSION; }

"krylovtype"|"KRYLOVTYPE"          { return KRYLOVTYPE; }

"savememcoarse"|"SAVEMEMCOARSE"    { return SAVEMEMCOARSE; }

"HELMHOLTZ"|"helmholtz"            { return HELMHOLTZ; }

"DEM"                              { return DEM; }

"k"|"K"                            { return WAVENUMBER; }

"bgtl"|"BGTL"                      { return BGTL; }

"incidence"|"INCIDENCE"            { return INCIDENCE; }

"pointsource"|"POINTSOURCE"        { return POINTSOURCE; }

"planewave"|"PLANEWAVE"            { return PLANEWAVE; }

"ksweep"|"KSWEEP"                  { return HELMSWEEP; }

"ksweep1"|"KSWEEP1"                { return HELMSWEEP1; }

"ksweep2"|"KSWEEP2"                { return HELMSWEEP2; }

"HDIR"[a-zA-Z]*|"hdir"[a-zA-Z]*    { return HDIRICHLET; }

"HNEU"[a-zA-Z]*|"hneu"[a-zA-Z]*    { return HNEUMAN; }

"HNBO"[a-zA-Z]*|"hnbo"[a-zA-Z]*    { return HNBO; }

"HDNB"[a-zA-Z]*|"hdnb"[a-zA-Z]*    { return HNBO; }

"HARB"[a-zA-Z]*|"harb"[a-zA-Z]*    { return ELHSOMMERFELD; }

"HWIB"[a-zA-Z]*|"hwib"[a-zA-Z]*    { return HWIBO; }

"HELMMF"|"helmmf"                  { return HELMMF; }

"HELMSO"[a-zA-Z]*|"helmso"[a-zA-Z]* { return HELMSO; }

"HSCB"[a-zA-Z]*|"hscb"[a-zA-Z]*    { return HSCBO; }

"ATDARB"[a-zA-Z]*|"atdarb"[a-zA-Z]* { return ATDARB; }

"ATDNEU"[a-zA-Z]*|"atdneu"[a-zA-Z]* { return ATDNEU; }

"ATDDNB"[a-zA-Z]*|"atddnb"[a-zA-Z]* { return ATDDNB; }

"ATDDIR"[a-zA-Z]*|"atddir"[a-zA-Z]* { return ATDDIR; }

"ATDROB"[a-zA-Z]*|"atdrob"[a-zA-Z]* { return ATDROB; }

"BLOCKDIAG" 			   { return BLOCKDIAG; }

"FETI"|"feti"                      { return FETI; }

"OLD"|"old"                        { yylval.ival = 0; return FETI2TYPE; }

"NEW"|"new"                        { yylval.ival = 1; return FETI2TYPE; }

"DP"|"dp"                          { return DP; }

"DPH"|"dph"                        { return DPH; }

"noco"[a-zA-Z]*|"NOCO"[a-zA-Z]*    { return NOCOARSE; }

"proj"[a-zA-Z]*|"PROJ"[a-zA-Z]*    { return PROJ; }

"outerloop"|"OUTERLOOP"            { return OUTERLOOP; }
"interf_solver"|"INTERF_SOLVER"    { return OUTERLOOP; }

"hermitian"                        { return HERMITIAN; }

"local_solver"|"LOCAL_SOLVER"      { return LOCALSOLVER; }

"coarse_solver"|"COARSE_SOLVER"    { return COARSESOLVER; }

"aux_coarse_solver" { return AUXCOARSESOLVER; }

"cct_solver"|"CCT_SOLVER"          { return CCTSOLVER; }

"rebuildcct"|"REBUILDCCT"          { return REBUILDCCT; }

"skyl"[a-zA-Z]*|"SKYL"[a-zA-Z]*    { yylval.ival = 0; return SOLVERTYPE; }

"spar"[a-zA-Z]*|"SPAR"[a-zA-Z]*    { yylval.ival = 1; return SOLVERTYPE; }

"sgisparse"[a-zA-Z]*|"SGISPARSE"[a-zA-Z]* { yylval.ival = 2; return SOLVERTYPE; }

"sgisky"[a-zA-Z]*|"SGISKY"[a-zA-Z]* { yylval.ival = 3; return SOLVERTYPE; }

"fron"[a-zA-Z]*|"FRON"[a-zA-Z]*    { yylval.ival = 5; return SOLVERTYPE; }

"bloc"[a-zA-Z]*|"BLOC"[a-zA-Z]*    { yylval.ival = 6; return SOLVERTYPE; }

"para"[a-zA-Z]*|"PARA"[a-zA-Z]*    { yylval.ival = 7; return SOLVERTYPE; }

"spooles"[a-zA-Z]*|"SPOOLES"[a-zA-Z]* { yylval.ival = 8; return SOLVERTYPE; }

"mumps"|"MUMPS" { yylval.ival = 9; return SOLVERTYPE; }

"purelydiag"[a-zA-Z]*|"PURELYDIAG"[a-zA-Z]* { yylval.ival = 10; return SOLVERTYPE; }

"mumps_icntl"|"MUMPS_ICNTL" { return MUMPSICNTL; }

"mumps_cntl"|"MUMPS_CNTL" { return MUMPSCNTL; }

"pivot"[a-zA-Z]*|"PIVOT"[a-zA-Z]*  { return PIVOT; }

"scaled"[a-zA-Z]*|"SCALED"[a-zA-Z]* { return SCALED; }

"spooles_tau"|"SPOOLES_TAU"        { return SPOOLESTAU; }

"spooles_maxsize"|"SPOOLES_MAXSIZE" { return SPOOLESMAXSIZE; }

"spooles_maxdomainsize"|"SPOOLES_MAXDOMAINSIZE" { return SPOOLESMAXDOMAINSIZE; }

"spooles_seed"|"SPOOLES_SEED"      { return SPOOLESSEED; }

"spooles_maxzeros"|"SPOOLES_MAXZEROS" { return SPOOLESMAXZEROS; }

"spooles_msglvl"|"SPOOLES_MSGLVL" { return SPOOLESMSGLVL; }

"spooles_scale"|"SPOOLES_SCALE" { return SPOOLESSCALE; }

"spooles_pivot"|"SPOOLES_PIVOT" { return SPOOLESPIVOT; }

"spooles_renum"|"SPOOLES_RENUM" { return SPOOLESRENUM; }

"sparse_maxsup"|"SPARSE_MAXSUP"    { return SPARSEMAXSUP; }

"sparse_defblk"|"SPARSE_DEFBLK"    { return SPARSEDEFBLK; }

"kryp"[a-zA-Z]*|"KRYP"[a-zA-Z]*    { return NLPREC; }

"tolf"[a-zA-Z]*|"TOLF"[a-zA-Z]     { return TOLFETI; }

"stag"[a-zA-Z_]*|"STAG"[a-zA-Z_]*    { return STAGTOL; }

"global_rbm_tol"|"GLOBAL_RBM_TOL"  { return GLOBALTOL; }

"global_cor_rbm_tol"|"GLOBAL_COR_RBM_TOL" { return GLOBALCRBMTOL; }

"cct_tol"|"CCT_TOL"                { return CCTTOL; }

"maxor"[a-zA-Z]*|"MAXOR"[a-zA-Z]*  { return MAXORTHO; }

"vers"[a-zA-Z]*|"VERS"[a-zA-Z]*    { return VERSION; }

"USDF"[a-zA-Z]*|"usdf"[a-zA-Z]*    { return USERDEFINEFORCE; }

"USDD"[a-zA-Z]*|"usdd"[a-zA-Z]*    { return USERDEFINEDISP; }

"printMat"[a-zA-Z]*|printmat[a-zA-Z]*|PRINTMAT[a-zA-Z]* { return PRINTMATLAB; }

"u_proj"[a-zA-Z]*|"U_PROJ"[a-zA-Z]*  { return UPROJ; }

"nulls"[a-zA-Z]*|"NULLS"[a-zA-Z]*  { return NULLSPACE; }

"mpc_type"|"MPC_TYPE"              { return MPCTYPE; }

"dual"                             { yylval.ival = 1; return MPCTYPEID; }

"primal"|"PRIMAL"                  { yylval.ival = 2; return MPCTYPEID; }

"mpc_precn"[a-zA-Z]*|"MPC__PRECN"[a-zA-Z]* { return MPCPRECNO; }

"diag"[a-zA-Z]*|"DIAG"[a-zA-Z]*    { yylval.ival = 1; return MPCPRECNOID; }

"full"[a-zA-Z]*|"FULL"[a-zA-Z]*    { yylval.ival = 2; return MPCPRECNOID; }

"tblock"[a-zA-Z]*|"TBLOCK"[a-zA-Z]* { yylval.ival = 3; return MPCPRECNOID; }

"sblock"[a-zA-Z]*|"SBLOCK"[a-zA-Z]* { yylval.ival = 4; return MPCPRECNOID; }

"mblock"[a-zA-Z]*|"MBLOCK"[a-zA-Z]* { yylval.ival = 5; return MPCPRECNOID; }

"blk_overlap"[a-zA-Z]*|"BLK_OVERLAP"[a-zA-Z]* { return MPCBLK_OVERLAP; }

"c_normalize" { return CNORM; }

"mpc_check" { return MPCCHECK; }

"gmres_residual"|"GMRES_RESIDUAL"  { return GMRESRESIDUAL; }

"pick_any_corner" { return PICKANYCORNER; }

"DIRE"[a-zA-Z]*|"dire"[a-zA-Z]*    { return DIRECT; }

"MULT"[a-zA-Z]*|"mult"[a-zA-Z]*    { return MULTIPLIERS; }

"CONSTRAINTS"|"constraints"        { return CONSTRAINTS; }

"PENA"[a-zA-Z]*|"pena"[a-zA-Z]*    { return PENALTY; }

"RENU"[a-zA-Z]*|"renu"[a-zA-Z]*    { return RENUM; }

"DYNA"[a-zA-Z]*|"dyna"[a-zA-Z]*    { return DYNAM; }

"STAB"[a-zA-Z]*|"stab"[a-zA-Z]*    { return STABLE; }

"rbmf"[a-zA-Z]*|"RBMF"[a-zA-Z]*    { return RBMFILTER; }

"RESTART"|"restart"                { return RESTART; }

"CONV"[a-zA-Z]*|"conv"[a-zA-Z]*    { return CONVECTION; }

"RADI"[a-zA-Z]*|"radi"[a-zA-Z]*    { return RADIATION; }

"newm"[a-zA-Z]*|"NEWM"[a-zA-Z]*    { return NEWMARK; }

"QSTA"[a-zA-Z]*|"qsta"[a-zA-Z]*    { return QSTATIC; }

"modal"|"MODAL"                    { return MODAL; }

"precn"o[a-zA-Z]*|"PRECN"O[a-zA-Z]* { return PRECNO; }

"prectype"|"PRECTYPE"              { return PRECTYPE; }

"nonshifted"|"NONSHIFTED"          { yylval.ival = 0; return PRECTYPEID; }

"shifted"|"SHIFTED"                { yylval.ival = 1; return PRECTYPEID; }

"tolp"[a-zA-Z]*|"TOLP"[a-zA-Z]*    { return TOLPCG; }

"tolb"[a-zA-Z]*|"TOLB"[a-zA-Z]*    { return TOLPCG; }

"tolcr"|"TOLCR"    { return TOLPCG; }

"tolg"[a-zA-Z]*|"TOLG"[a-zA-Z]*    { return TOLPCG; }

"subt"[a-zA-Z]*|"SUBT"[a-zA-Z]*    { return SUBTYPE; }

"maxit"[a-zA-Z]*|"MAXIT"[a-zA-Z]*  { return MAXITR; }

"propo"[a-zA-Z]*|"PROPO"[a-zA-Z]*  { return PROPORTIONING; }

"lines"[a-zA-Z]*|"LINES"[a-zA-Z]*  { return LINESEARCH; }

"debug_icntl"                      { return DEBUGICNTL; }

"debug_cntl"                       { return DEBUGCNTL; }

"ptol"[a-zA-Z]*|"PTOL"[a-zA-Z]*    { return PTOL; }

"pmax"[a-zA-Z]*|"PMAX"[a-zA-Z]*    { return PMAXIT; }

"plan"[a-zA-Z]*|"PLAN"[a-zA-Z]*    { return PLANTOL; }

"nsbs"[a-zA-Z]*|"NSBS"[a-zA-Z]*    { return NSBSPV; }

"neig"[a-zA-Z]*|"NEIG"[a-zA-Z]*    { return NEIGPA; }

"toleig"[a-zA-Z]*|"TOLEIG"[a-zA-Z]* { return TOLEIG; }

"tolj"[a-zA-Z]*|"TOLJ"[a-zA-Z]*    { return TOLJAC; }

"expl"[a-zA-Z]*|"EXPL"[a-zA-Z]*    { return EXPLICIT; }

"lobpcg"[a-zA-Z]*|"LOBPCG"[a-zA-Z]* { return LOBPCG; }

"arpack"|"ARPACK"                  { return ARPACK; }

"filtereig"|"FILTEREIG"            { return FILTEREIG; }

"jacobi"|"JACOBI"                  { return JACOBI; }

"damp"[a-zA-Z]*|"DAMP"[a-zA-Z]*    { return DAMPING; }

"time"|"TIME"                      { return TIME; }

"pita"|"PITA"|"Pita"               { return PITA; }

"mdpita"|"MDPITA"|"MdPita"         { return MDPITA; }

"NoForce"|"noforce"|"NOFORCE"      { return NOFORCE; }

"ConstForce"|"constforce"|"CONSTFORCE" { return CONSTFORCE; }

"CkCoarse"|"ckcoarse"|"CKCOARSE"   { return CKCOARSE; }

"localbases"|"LocalBases"|"LOCALBASES" { return LOCALBASES; }

"remotecoarse"|"REMOTECOARSE"      { return REMOTECOARSE; }

"ORTHOPROJTOL"|"orthoprojtol"      { return ORTHOPROJTOL; }

"newimpl"|"NEWIMPL"|"OOPITA"       { return NEWIMPL; }

"heat"|"HEAT"                      { return HEAT; }

"mech"|"MECH"                      { return MECH; }

"acou"[a-zA-Z]*|"ACOU"[a-zA-Z]*    { return ACOU; }

"OUTPUT"|"output"                  { return OUTPUT; }

"OUTPUT6"|"outpu6"                 { return OUTPUT6; }

"EIGE"[a-zA-Z]*|"eige"[a-zA-Z]*    { return EIGEN; }

"addedmass"                        { return ADDEDMASS; }

"SLOS"[a-zA-Z]*|"slos"[a-zA-Z]*    { return SLOSH; }

"SLGR"[a-zA-Z]*|"slgr"[a-zA-Z]*    { return SLGRAV; }

"HEFS"[a-zA-Z]*|"hefs"[a-zA-Z]*    { return HEFSB; }

"HEIN"[a-zA-Z]*|"hein"[a-zA-Z]*    { return HEINTERFACE; }

"subs"[a-zA-Z]*|"subs"[a-zA-Z]*    { return SUBSPACE;  }

"NONL"[a-zA-Z]*|"nonl"[a-zA-Z]*    { return NL; }

"nlmat"[a-zA-Z]*|"nlmat"[a-zA-Z]*  { return NLMAT; }

"FITA"[a-zA-Z]*|"fita"[a-zA-Z]*    { return FITALG; }

"nlto"[a-zA-Z]*|"NLTO"[a-zA-Z]*    { return NLTOL; }

"maxv"[a-zA-Z]*|"MAXV"[a-zA-Z]*    { return MAXVEC; }

"arcl"[a-zA-Z]*|"ARCL"[a-zA-Z]*    { return ARCLENGTH; }

"dlam"[a-zA-Z]*|"DLAM"[a-zA-Z]*    { return DLAMBDA; }

"lfac"[a-zA-Z]*|"LFAC"[a-zA-Z]*    { return LFACTOR; }

"unsy"[a-zA-Z]*|"UNSY"[a-zA-Z]*    { return UNSYMMETRIC; }

"rebuild"|"REBUILD"                { return REBUILD; }

"tang"[a-zA-Z]*|"TANG"[a-zA-Z]*    { return TANGENT; }

"preco"[a-zA-Z]*|"PRECO"[a-zA-Z]*  { return PRECONDITIONER; }

"DIMA"[a-zA-Z]*|"dima"[a-zA-Z]*    { return DIMASS; }

"COMP"[a-zA-Z]*|"composite"        { return COMPOSITE; }

"COEF"|"coef"                      { return COEF; }

"CFRA"[a-zA-Z]*|"cfra"[a-zA-Z]*    { return CFRAMES; }

"LAYC"|"layc"                      { return LAYC; }

"LAYN"|"layn"                      { return LAYN; }

"LAYD"|"layd"                      { return LAYD; }

"LAYO"|"layo"                      { return LAYO; }

"LAYMAT"|"laymat"                  { return LAYMAT; }

"TRBM"[a-zA-Z]*|"trbm"[a-zA-Z]*    { return TRBM; }

"GRBM"[a-zA-Z]*|"grbm"[a-zA-Z]*    { return GRBM; }

"HZEM"|"hzem"                      { return HZEM; }

"HZEMF"[a-zA-Z]*|"hzemf"[a-zA-Z]*  { return HZEMFILTER; }

"SZEM"[a-zA-Z]*|"szem"[a-zA-Z]*    { return SLZEM; }

"SZFI"[a-zA-Z]*|"szfi"[a-zA-Z]*    { return SLZEMFILTER; }

"GRAV"[a-zA-Z]*|"grav"[a-zA-Z]*    { return GRAVITY; }

"MFTT"|"mftt"                      { return MFTT; }

"MPTT"|"mptt"                      { return MPTT; }

"HFTT"|"hftt"                      { return HFTT; }

"YMTT"|"ymtt"                      { return YMTT; }

"curve"|"CURVE"                    { return CURVE; }

"TETT"|"tett"                      { return TETT; }

"REOR"[a-zA-Z]*|"reor"[a-zA-Z]*    { return REORTHO; }

"scal"[a-zA-Z]*|"SCAL"[a-zA-Z]*    { return SCALING; }

"mpc_scal"[a-zA-Z]*|"MPC_SCAL"[a-zA-Z]* { return MPCSCALING; }

"fsi_scal"[a-zA-Z]*|"FSI_SCAL"[a-zA-Z]* { return FSISCALING; }

"mpc_ele"[a-zA-Z]*|"MPC_ELE"[a-zA-Z]* { return MPCELEMENT; }

"fsi_ele"[a-zA-Z]*|"FSI_ELE"[a-zA-Z]* { return FSIELEMENT; }

"fsi_cor"[a-zA-Z]*|"FSI_COR"[a-zA-Z]* { return FSICORNER; }

"coupled_scale"      { return COUPLEDSCALE; }

"wet_corners"        { return WETCORNERS; }

"noLocalFsiSpliting" { return NOLOCALFSISPLITING; }

"stif"[a-zA-Z]*|"STIF"[a-zA-Z]*    { yylval.ival = 1; return SCALINGTYPE; }

"topological" { yylval.ival = 2; return SCALINGTYPE; }

"sloan"|"SLOAN"    { yylval.ival = 1; return RENUMBERID; }

"rcm"|"RCM"        { yylval.ival = 2; return RENUMBERID; }

"mmd"|"MMD"        { yylval.ival = 3; return RENUMBERID; }

"nd"|"ND"          { yylval.ival = 4; return RENUMBERID; }

"ms"|"MS"          { yylval.ival = 5; return RENUMBERID; }

"ndms"|"NDMS"      { yylval.ival = 6; return RENUMBERID; }

"esmo"[a-zA-Z]*|"ESMO"[a-zA-Z]*    { yylval.ival = 3; return RENUMBERID; }
"meti"[a-zA-Z]*|"METI"[a-zA-Z]*    { yylval.ival = 4; return RENUMBERID; }

"gtgs"[a-zA-Z]*|"GTGS"[a-zA-Z]*    { yylval.ival = FetiInfo::sparse; return GTGSOLVER; }

"gtgp"[a-zA-Z]*|"GTGP"[a-zA-Z]*    { yylval.ival = FetiInfo::pcg;    return GTGSOLVER; }

"corn"[a-zA-Z]*|"CORN"[a-zA-Z]*    { return CORNER; }

"sixdof"|"SIXDOF"                  { yylval.ival = 0; return CORNERTYPE; }

"threedof"|"THREEDOF"              { yylval.ival = 1; return CORNERTYPE; }

"noEnd6"|"noend6"|"NOEND6"         { yylval.ival = 2; return CORNERTYPE; }

"noEnd3"|"noend3"|"NOEND3"         { yylval.ival = 3; return CORNERTYPE; }

"be6"|"BE6"                        { yylval.ival = 0; return CORNERTYPE; }

"be3"|"BE3"                        { yylval.ival = 1; return CORNERTYPE; }

"cp6"|"CP6"                        { yylval.ival = 2; return CORNERTYPE; }

"cp3"|"CP3"                        { yylval.ival = 3; return CORNERTYPE; }

"augment"|"AUGMENT"                { return AUGMENT; }

"none"|"NONE"                      { yylval.ival = 0; return AUGMENTTYPE;}

"Gs"|"gs"|"GS"                     { yylval.ival = 1; return AUGMENTTYPE;}

"EdgeGs"|"edgegs"|"EDGEGS"         { yylval.ival = 2; return AUGMENTTYPE;}

"WeightedEdgeGs"|"weightededgegs"|"WEIGHTEDEDGEGS" { yylval.ival = 3; return AUGMENTTYPE;}

"EdgeWs"|"edgews"|"EDGEWS"         { return EDGEWS; }

"solid"|"SOLID"                    { yylval.ival = 0; return WAVETYPE;}

"shell"|"SHELL"                    { yylval.ival = 1; return WAVETYPE;}

"fluid"|"FLUID"                    { yylval.ival = 2; return WAVETYPE;}

"FLUMAT"|"flumat"                    { return FLUMAT;}
"FABMAT"			{ return FABMAT;}
"THERMMAT"|"thermmat"                    { return THERMMAT; }

"any"|"ANY"                        { yylval.ival = 3; return WAVETYPE; }

"averageK"|"averagek"|"AVERAGEK"   { yylval.ival = 0; return WAVEMETHOD;}

"averageMat"|"averagemat"|"AVERAGEMAT" { yylval.ival = 1; return WAVEMETHOD;}

"uniform"|"UNIFORM"                { yylval.ival = 2; return WAVEMETHOD;}

"orthotol"|"ORTHOTOL"              { return ORTHOTOL; } 

"trans"[a-zA-Z]*|"TRANS"[a-zA-Z]*  { yylval.ival = FetiInfo::translation; return RBMSET; }

"rotat"[a-zA-Z]*|"ROTAT"[a-zA-Z]*  { yylval.ival = FetiInfo::rotation;    return RBMSET; }

"all"[a-zA-Z]*|"ALL"[a-zA-Z]*      { yylval.ival = FetiInfo::all;         return RBMSET; }

"averageTrans"[a-zA-Z]*|"averagetrans"[a-zA-Z]*|"AVERAGETRANS"[a-zA-Z]* { yylval.ival = FetiInfo::averageTran; return RBMSET; }

"averageRot"[a-zA-Z]*|"averagerot"[a-zA-Z]*|"AVERAGEROT"[a-zA-Z]* { yylval.ival = FetiInfo::averageRot;  return RBMSET; }

"averageAll"[a-zA-Z]*|"averageall"[a-zA-Z]*|"AVERAGEALL"[a-zA-Z]* { yylval.ival = FetiInfo::averageAll;  return RBMSET; }

"noprec"[a-zA-Z]*|"NOPREC"[a-zA-Z]* { yylval.ival = FetiInfo::noPrec;    return FETIPREC; }

"lump"[a-zA-Z]*|"LUMP"[a-zA-Z]*    { return LUMPED; }

"diri"[a-zA-Z]*|"DIRI"[a-zA-Z]*    { yylval.ival = FetiInfo::dirichlet; return FETIPREC; }

"idprec"[a-zA-Z]*|"IDPREC"[a-zA-Z]*  { yylval.ival = FetiInfo::identity; return FETIPREC; }

"mrhs"|"MRHS"                      { return MRHS; }

"AERO"|"aero"                      { return AERO; }

"AEROH"|"aeroh"                    { return AEROH; }

"THERMOH"|"thermoh"                { return THERMOH; }

"THERMOE"|"thermoe"                { return THERMOE; }

"COLL"[a-zA-Z]*|"coll"[a-zA-Z]*    { yylval.ival = 1; return COLLOCATEDTYPE; }

"NONC"[a-zA-Z]*|"nonc"[a-zA-Z]*    { yylval.ival = 0; return COLLOCATEDTYPE; }

"A0"|"a0"                          { yylval.ival = 0; return AEROTYPE; }

"PP"|"pp"                          { yylval.ival = 1; return AEROTYPE; }

"A4"|"a4"                          { yylval.ival = 4; return AEROTYPE; }

"A5"|"a5"                          { yylval.ival = 5; return AEROTYPE; }

"A6"|"a6"                          { yylval.ival = 6; return AEROTYPE; }

"A7"|"a7"                          { yylval.ival = 7; return AEROTYPE; }

"MPP"|"mpp"                        { yylval.ival = 8; return AEROTYPE; }

"B0"|"b0"                          { yylval.ival = 10; return AEROTYPE; }

"C0"|"c0"                          { yylval.ival = 20; return AEROTYPE; }

"DX"|"dx"                          { yylval.ival = OutputInfo::DispX;  return DOFTYPE; }

"DY"|"dy"                          { yylval.ival = OutputInfo::DispY;  return DOFTYPE; }

"DZ"|"dz"                          { yylval.ival = OutputInfo::DispZ;  return DOFTYPE; }

"RX"|"rx"                          { yylval.ival = OutputInfo::RotX;   return DOFTYPE; }

"RY"|"ry"                          { yylval.ival = OutputInfo::RotY;   return DOFTYPE; }

"RZ"|"rz"                          { yylval.ival = OutputInfo::RotZ;   return DOFTYPE; }

"DM"|"dm"                          { yylval.ival = OutputInfo::DispMod;return DOFTYPE; }

"RM"|"rm"                          { yylval.ival = OutputInfo::RotMod; return DOFTYPE; }

"TM"|"tm"                          { yylval.ival = OutputInfo::TotMod; return DOFTYPE; }

"elemental"[a-zA-Z]*|"ELEMENTAL"[a-zA-Z]* { yylval.ival =  0; return AVERAGED; }

"nodal"|"NODAL"|"nodalfull"|"NODALFULL"  { yylval.ival =  1; return AVERAGED; }

"nodalpartial"|"NODALPARTIAL"    { yylval.ival =  2; return AVERAGED; }

"upper"[a-zA-Z]*|"UPPER"[a-zA-Z]*  { yylval.ival = 1; return SURFACE; }

"median"[a-zA-Z]*|"MEDIAN"[a-zA-Z]* { yylval.ival = 2; return SURFACE; }

"lower"[a-zA-Z]*|"LOWER"[a-zA-Z]*  { yylval.ival = 3; return SURFACE; }

"gdisplac"|"GDISPLAC"              { yylval.ival = OutputInfo::Displacement; return STRESSID; }

"disp6dof"|"DISP6DOF"              { yylval.ival = OutputInfo::Disp6DOF;     return STRESSID; }

"modaldsp"|"MODALDSP"              { yylval.ival = OutputInfo::ModalDsp;     return STRESSID; }

"modalexf"|"MODALEXF"              { yylval.ival = OutputInfo::ModalExF;     return STRESSID; }

"gtempera"|"GTEMPERA"              { yylval.ival = OutputInfo::Temperature;  return STRESSID; }

"stressxx"|"STRESSXX"              { yylval.ival = OutputInfo::StressXX;     return STRESSID; }

"stressyy"|"STRESSYY"              { yylval.ival = OutputInfo::StressYY;     return STRESSID; }

"stresszz"|"STRESSZZ"              { yylval.ival = OutputInfo::StressZZ;     return STRESSID; }

"stressxy"|"STRESSXY"              { yylval.ival = OutputInfo::StressXY;     return STRESSID; }

"stressxz"|"STRESSXZ"              { yylval.ival = OutputInfo::StressXZ;     return STRESSID; }

"stresszx"|"STRESSZX"              { yylval.ival = OutputInfo::StressXZ;     return STRESSID; }

"stressyz"|"STRESSYZ"              { yylval.ival = OutputInfo::StressYZ;     return STRESSID; }

"strainxx"|"STRAINXX"              { yylval.ival = OutputInfo::StrainXX;     return STRESSID; }

"strainyy"|"STRAINYY"              { yylval.ival = OutputInfo::StrainYY;     return STRESSID; }

"strainzz"|"STRAINZZ"              { yylval.ival = OutputInfo::StrainZZ;     return STRESSID; }

"strainxy"|"STRAINXY"              { yylval.ival = OutputInfo::StrainXY;     return STRESSID; }

"strainzx"|"STRAINZX"              { yylval.ival = OutputInfo::StrainXZ;     return STRESSID; }

"strainxz"|"STRAINXZ"              { yylval.ival = OutputInfo::StrainXZ;     return STRESSID; }

"strainyz"|"STRAINYZ"              { yylval.ival = OutputInfo::StrainYZ;     return STRESSID; }

"heatflxx"|"HEATFLXX"              { yylval.ival = OutputInfo::HeatFlXX;     return STRESSID; }

"heatflxy"|"HEATFLXY"              { yylval.ival = OutputInfo::HeatFlXY;     return STRESSID; }

"heatflxz"|"HEATFLXZ"              { yylval.ival = OutputInfo::HeatFlXZ;     return STRESSID; }

"grdtempx"|"GRDTEMPX"              { yylval.ival = OutputInfo::GrdTempX;     return STRESSID; }

"grdtempy"|"GRDTEMPY"              { yylval.ival = OutputInfo::GrdTempY;     return STRESSID; }

"grdtempz"|"GRDTEMPZ"              { yylval.ival = OutputInfo::GrdTempZ;     return STRESSID; }

"sdisplac"|"SDISPLAC"              { yylval.ival = OutputInfo::SloshDisplacement;     return STRESSID; }

"slshdspx"|"SLSHDSPX"              { yylval.ival = OutputInfo::SloshDispX;     return STRESSID; }

"slshdspy"|"SLSHDSPY"              { yylval.ival = OutputInfo::SloshDispY;     return STRESSID; }

"slshdspz"|"SLSHDSPZ"              { yylval.ival = OutputInfo::SloshDispZ;     return STRESSID; }

"stressvm"|"STRESSVM"              { yylval.ival = OutputInfo::StressVM;     return STRESSID; }

"stressp1"|"STRESSP1"              { yylval.ival = OutputInfo::StressPR1;    return STRESSID; }

"stressp2"|"STRESSP2"              { yylval.ival = OutputInfo::StressPR2;    return STRESSID; }

"stressp3"|"STRESSP3"              { yylval.ival = OutputInfo::StressPR3;    return STRESSID; }

"strainp1"|"STRAINP1"              { yylval.ival = OutputInfo::StrainPR1;    return STRESSID; }

"strainp2"|"STRAINP2"              { yylval.ival = OutputInfo::StrainPR2;    return STRESSID; }

"strainp3"|"STRAINP3"              { yylval.ival = OutputInfo::StrainPR3;    return STRESSID; }

"inxforce"|"INXFORCE"              { yylval.ival = OutputInfo::InXForce;     return STRESSID; }

"inyforce"|"INYFORCE"              { yylval.ival = OutputInfo::InYForce;     return STRESSID; }

"inzforce"|"INZFORCE"              { yylval.ival = OutputInfo::InZForce;     return STRESSID; }

"axmoment"|"AXMOMENT"              { yylval.ival = OutputInfo::AXMoment;     return STRESSID; }

"aymoment"|"AYMOMENT"              { yylval.ival = OutputInfo::AYMoment;     return STRESSID; }

"azmoment"|"AZMOMENT"              { yylval.ival = OutputInfo::AZMoment;     return STRESSID; }

"energies"|"ENERGIES"              { yylval.ival = OutputInfo::Energies;     return STRESSID; }

"raerofor"|"RAEROFOR"              { yylval.ival = OutputInfo::AeroForce;    return STRESSID; }

"geigenpa"|"GEIGENPA"              { yylval.ival = OutputInfo::EigenPair;    return STRESSID; }

"geig6dof"|"GEIG6DOF"              { yylval.ival = OutputInfo::Disp6DOF;     return STRESSID; }

"strainvm"|"STRAINVM"              { yylval.ival = OutputInfo::StrainVM;     return STRESSID; }

"ghelmhol"|"gacoustd"|"pressure"|"gacoupre"   { yylval.ival = OutputInfo::Helmholtz;    return STRESSID; }

"GHELMHOL"|"GACOUSTD"|"PRESSURE"|"GACOUPRE"   { yylval.ival = OutputInfo::Helmholtz;    return STRESSID; }

"aeroforx"|"AEROFORX"              { yylval.ival = OutputInfo::AeroXForce;   return STRESSID; }

"aerofory"|"AEROFORY"              { yylval.ival = OutputInfo::AeroYForce;   return STRESSID; }

"aeroforz"|"AEROFORZ"              { yylval.ival = OutputInfo::AeroZForce;   return STRESSID; }

"aeromomx"|"AEROMOMX"              { yylval.ival = OutputInfo::AeroXMom;     return STRESSID; }

"aeromomy"|"AEROMOMY"              { yylval.ival = OutputInfo::AeroYMom;     return STRESSID; }

"aeromomz"|"AEROMOMZ"              { yylval.ival = OutputInfo::AeroZMom;     return STRESSID; }

"gvelocit"|"GVELOCIT"              { yylval.ival = OutputInfo::Velocity;     return STRESSID; }

"gacceler"|"GACCELER"              { yylval.ival = OutputInfo::Acceleration; return STRESSID; }

"gvel6dof"|"GVEL6DOF"              { yylval.ival = OutputInfo::Velocity6;     return STRESSID; }

"gacc6dof"|"GACC6DOF"              { yylval.ival = OutputInfo::Accel6; return STRESSID; }

"gtempvel"|"GTEMPVEL"              { yylval.ival = OutputInfo::TemperatureFirstTimeDerivative; return STRESSID; }

"gacouvel"|"GACOUVEL"              { yylval.ival = OutputInfo::PressureFirstTimeDerivative; return STRESSID; }

"gacouacc"|"GACOUACC"              { yylval.ival = OutputInfo::PressureSecondTimeDerivative; return STRESSID; }

"ymodulus"|"YMODULUS"              { yylval.ival = OutputInfo::YModulus;     return STRESSID; }

"mdensity"|"MDENSITY"              { yylval.ival = OutputInfo::MDensity;     return STRESSID; }

"thicknes"|"THICKNES"              { yylval.ival = OutputInfo::Thicknes;     return STRESSID; }

"shapeatt"|"SHAPEATT"              { yylval.ival = OutputInfo::ShapeAtt;     return STRESSID; }

"shapestc"|"SHAPESTC"              { yylval.ival = OutputInfo::ShapeStc;     return STRESSID; }

"rigid"|"RIGID"                    { yylval.ival = OutputInfo::Rigid;        return STRESSID; }

"composit"|"COMPOSIT"              { yylval.ival = OutputInfo::Composit;     return STRESSID; }

"eletonod"|"ELASTONOD"             { yylval.ival = OutputInfo::ElemToNode;   return STRESSID; }

"nodtoele"|"NODTOELE"              { yylval.ival = OutputInfo::NodeToElem;   return STRESSID; }

"nodtonod"|"NODTONOD"              { yylval.ival = OutputInfo::NodeToNode;   return STRESSID; }

"raerotfl"|"RAEROTFL"              { yylval.ival = OutputInfo::AeroFlux;     return STRESSID; }

"heatflx"|"HEATFLX"                { yylval.ival = OutputInfo::HeatFlX;       return STRESSID; }

"grdtemp"|"GRDTEMP"                { yylval.ival = OutputInfo::GrdTemp;       return STRESSID; }

"velo6dof"|"VELO6DOF"              { yylval.ival = OutputInfo::Velocity6;    return STRESSID; }

"accl6dof"|"ACCL6DOF"              { yylval.ival = OutputInfo::Accel6;       return STRESSID; }

"gencoord"|"GENCOORD"              { yylval.ival = OutputInfo::ModeAlpha;    return STRESSID; }

"moderror"|"MODERROR"              { yylval.ival = OutputInfo::ModeError;    return STRESSID; }

"reaction"[s]?|"REACTION"[s]?      { yylval.ival = OutputInfo::Reactions; return STRESSID; }

"conpress"|"CONPRESS"              { yylval.ival = OutputInfo::ContactPressure; return STRESSID; }

"farfield"|"FARFIELD"              { yylval.ival = OutputInfo::Farfield;     return STRESSID; }

"geigpres"|"GEIGPRES"              { yylval.ival = OutputInfo::EigenPressure; return STRESSID; }

"geigslsh"|"GEIGSLSH"              { yylval.ival = OutputInfo::EigenSlosh; return STRESSID; }

"gfrmodes"|"GFRMODES"              { yylval.ival = OutputInfo::FreqRespModes; return STRESSID; }

"ghhmodes"|"GHHMODES"              { yylval.ival = OutputInfo::HelmholtzModes; return STRESSID; }

"sp1direc"|"SP1DIREC"              { yylval.ival = OutputInfo::StressPR1Direc; return STRESSID; }

"sp2direc"|"SP2DIREC"              { yylval.ival = OutputInfo::StressPR2Direc; return STRESSID; }

"sp3direc"|"SP3DIREC"              { yylval.ival = OutputInfo::StressPR3Direc; return STRESSID; }

"ep1direc"|"EP1DIREC"              { yylval.ival = OutputInfo::StrainPR1Direc; return STRESSID; }

"ep2direc"|"EP2DIREC"              { yylval.ival = OutputInfo::StrainPR2Direc; return STRESSID; }

"ep3direc"|"EP3DIREC"              { yylval.ival = OutputInfo::StrainPR3Direc; return STRESSID; }

"realimag"|"REALIMAG"              { yylval.ival = OutputInfo::realimag; return COMPLEXOUTTYPE; }

"modphase"|"MODPHASE"              { yylval.ival = OutputInfo::modulusphase; return COMPLEXOUTTYPE; }

"anim"[a-zA-Z]*|"ANIM"[a-zA-Z]*    { yylval.ival = OutputInfo::animate; return COMPLEXOUTTYPE; }

"conface"|"CONFACE"                                  { yylval.ival = 1; return TDENFORC; }
"normal_force_mag"|"NORMAL_FORCE_MAG"                { yylval.ival = 2; return TDENFORC; }
"normal_traction_mag"|"NORMAL_TRACTION_MAG"          { yylval.ival = 3; return TDENFORC; }
"tangential_force_mag"|"TANGENTIAL_FORCE_MAG"        { yylval.ival = 4; return TDENFORC; }
"tangential_traction_mag"|"TANGENTIAL_TRACTION_MAG"  { yylval.ival = 5; return TDENFORC; }
"cdirnorx"|"CDIRNORX"                                { yylval.ival = 6; return TDENFORC; }
"cdirnory"|"CDIRNORY"                                { yylval.ival = 7; return TDENFORC; }
"cdirnorz"|"CDIRNORZ"                                { yylval.ival = 8; return TDENFORC; }
"cdirtanx"|"CDIRTANX"                                { yylval.ival = 9; return TDENFORC; }
"cdirtany"|"CDIRTANY"                                { yylval.ival = 10; return TDENFORC; }
"cdirtanz"|"CDIRTANZ"                                { yylval.ival = 11; return TDENFORC; }
"slip_mag"|"SLIP_MAG"                                { yylval.ival = 12; return TDENFORC; }
"nodal_dissipation"|"NODAL_DISSIPATION"              { yylval.ival = 13; return TDENFORC; }
"contact_area"|"CONTACT_AREA"                        { yylval.ival = 14; return TDENFORC; }
"gap_cur"|"GAP_CUR"                                  { yylval.ival = 15; return TDENFORC; }
"gap_old"|"GAP_OLD"                                  { yylval.ival = 16; return TDENFORC; }

"damage"|"DAMAGE"                  { yylval.ival = OutputInfo::Damage; return STRESSID; }

"MATSPEC"|"matspec"                { return MATSPEC; }

"MATLAW"|"matlaw"                  { return MATSPEC; }

"MATUSAGE"|"matusage"              { return MATUSAGE; }

"bilinplast"|"BILINPLAST"          { return BILINPLAST; }

"linearplanestress"|"LINEARPLANESTRESS" { return LINPLSTRESS; }

"linear"|"LINEAR"                  { return LINEAR; }

"neohookean"|"NEOHOOKEAN"          { return NEOHOOKEAN; }

"simple"[a-zA-Z]*|"SIMPLE"[a-zA-Z]* { return SIMPLE; }

"read"|"READ"                      { return READ; }

"NODALCONTACT"|"nodalcontact"      { return NODALCONTACT; }

"BMPC"|"bmpc"                      { return BMPC; }
"DMPC"|"dmpc"                      { return DMPC; }
"CMPC"|"cmpc"                      { return CMPC; }

"FRIC"[a-zA-Z]*|"fric"[a-zA-Z]*    { return FRIC; }

"GAP"|"gap"                        { return GAP; }

"SURFTOPO"|"SURFACETOPO"           { return SURFACETOPOLOGY; }

"surftopo"|"surfacetopo"           { return SURFACETOPOLOGY; }

"MORTARTIED"|"mortartied"          { return MORTARTIED; }

"SEARCH_TOL"|"search_tol"          { return SEARCHTOL; }

"STD"|"std"                        { return STDMORTAR; }

"DUAL"|"dual"                      { return DUALMORTAR; }

"WETINTERFACE"|"wetinterface"      { return WETINTERFACE; }

"FSINTERFACE"|"fsinterface"        { return FSINTERFACE; }

"TIEDSURFACES"|"tiedsurfaces"      { return TIEDSURFACES; }

"CONTACTSURFACES"|"contactsurfaces" { return CONTACTSURFACES; }

"NO_SECONDARY"|"no_secondary"  { return NOSECONDARY; }

"SDISP"[a-zA-Z]*|"sdisp"[a-zA-Z]*  { return SDISP; }
"SFORC"[a-zA-Z]*|"sforc"[a-zA-Z]*  { return SFORCE; }
"SPRES"[a-zA-Z]*|"spres"[a-zA-Z]*  { return SPRESSURE; }
"ACME"[a-zA-Z]*|"acme"[a-zA-Z]*    { return ACMECNTL; }
"ffi_debug"                        { return FFIDEBUG; }

"NONINPC" { return NONINPC; }
"INPC" { return INPC; }
"GROUP"|"group" { return GROUP; }
"RAND"[a-zA-Z]*|"rand"[a-zA-Z]* { return RANDOM; }
"E"|"e" { yylval.rprop = GeoSource::E; return RPROP; }
"KX"|"kx"|"Kx" { yylval.rprop = GeoSource::KX; return RPROP; }
"KY"|"ky"|"Ky" { yylval.rprop = GeoSource::KY; return RPROP; }
"KZ"|"kz"|"Kz" { yylval.rprop = GeoSource::KZ; return RPROP; }
"MEAN"|"mean" { yylval.ival = 1; return NDTYPE; }
"STDV"|"stdv" { yylval.ival = 2; return NDTYPE; }
"PDF"|"pdf" { yylval.ival = 3; return NDTYPE; }

[A-Za-z][A-Za-z0-9_.]*  { yylval.strval = strdup(yytext); return FNAME; }

"\""[A-Za-z0-9_./]*"\"" { yytext[yyleng-1] = 0; yylval.strval = strdup(yytext+1);
                           return FNAME; }

"\n" { line_num++; /* Count line numbers */ return NewLine; }

{SPACES} /* Just discard spaces */

[+-]?{DIGIT}+ {yylval.ival = atoi(yytext) ; return IntConstant ; }

[+-]?{DIGIT}+"."{DIGIT}* |
[+-]?"."{DIGIT}+ |
[+-]?{DIGIT}+"."{DIGIT}*[eE][+-]?{DIGIT}+ |
[+-]?{DIGIT}+[eE][+-]?{DIGIT}+ |
[+-]?"."{DIGIT}+[eE][+-]?{DIGIT}+ {
       yylval.fval = atof(yytext) ; return DblConstant ; }

[+-]?{DIGIT}+"."{DIGIT}*[dD][+-]?{DIGIT}+ |
[+-]?{DIGIT}+[dD][+-]?{DIGIT}+ |
[+-]?"."{DIGIT}+[dD][+-]?{DIGIT}+ {
       int il;
       for(il = 0; yytext[il] != 'd' && yytext[il] != 'D';) ++il;
       yytext[il] = 'e'; yylval.fval = atof(yytext) ; return DblConstant ; }

"*" { while(yyinput() != '\n') ; line_num++; }



%%
void
yyerror(const char  *)
{
 	fprintf(stderr," ** Syntax Error near line %d : %s\n", line_num,yytext);
}

int yywrap()
{
	//std::cout << "yywrap" << endl;
 	return 1;
}
