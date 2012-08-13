#ifdef USE_EIGEN3
#include <cmath>
#include <stdexcept>
#include <Element.d/FelippaShell.d/ShellMaterial.hpp>

template<typename doublereal>
ShellMaterialTypes2And3<doublereal>::ShellMaterialTypes2And3(
  int _nlayer, doublereal *_mtlayer, bool _couple, doublereal *_aframe)
  : nlayer(_nlayer), mtlayer(_mtlayer,12,_nlayer), couple(_couple), aframe(_aframe)
{
// .....COMPUTE THE THICKNESS FOR TYPE-2 AND TYPE-3 CONSTITUTIVE LAWS 
// .....IT IS ASSUMED CONSTANT AND EQUAL TO THE SUM OF EACH LAYER'S THICKNESS 

    thick = 0;

    for (int ilayer = 0; ilayer < nlayer; ++ilayer) {
      thick += mtlayer(7, ilayer);
    }

// .....COMPUTE THE AREA DENSITY FOR TYPE-2 AND TYPE-3 CONSTITUTIVE LAWS 

    rhoh = 0;

    for (int ilayer = 0; ilayer < nlayer; ++ilayer) {
      rhoh += mtlayer(6, ilayer)*mtlayer(7, ilayer);
    }

}

template<typename doublereal>
void
ShellMaterialTypes2And3<doublereal>::GetConstitutiveResponse(doublereal *_Upsilon, doublereal *_Sigma, doublereal *_D,
                                                             doublereal *eframe, int gp)
{
    // Initialized data 
    doublereal one = 1.;

    // Builtin functions 
    using std::acos;

    // Local variables 
    int i, j, k;
    doublereal c, s, e1, e2, z0, g12, s11, s12, s13, s22, pi, s23, s33, mu1, mu2, nu12;
    int irot, ilayer;
    doublereal dets, thetaf, zsup, zinf, intthick;
    doublereal *data = (_D == NULL) ? new doublereal[36] : _D;
    Eigen::Map<Eigen::Matrix<doublereal,6,1> > Upsilon(_Upsilon), Sigma(_Sigma);
    Eigen::Map<Eigen::Matrix<doublereal,6,6> > D(data);
    Eigen::Block< Eigen::Map<Eigen::Matrix<doublereal,6,6> > >
        Dm = D.topLeftCorner(3,3),     Dmb = D.topRightCorner(3,3),
        Dbm = D.bottomLeftCorner(3,3), Db = D.bottomRightCorner(3,3);
    Eigen::Matrix<doublereal,3,3> C, invT;

// ==================================================================== 
//                                                                      
//     Perform =   Assembles the 6 by 6 Constitutive Matrix According   
//     ---------   to the Type of Constitutive Law Requested.           
//                                                                      
//                                                                      
//     Input/Output =                                                   
//     --------------                                                   
//     NLAYER  <input>  number of layers of the composite element       
//     MTLAYER <input>  material properties of each layer               
//     D       <output> 6 by 6 constitutive matrix                      
//     COUPLE  <input>  type of constitutive law                        
//     EFRAME  <input>  element level 3x3 frame                         
//     AFRAME  <input>  arbitrary 3x3 frame of the constitutive law     
//                                                                      
//                                                                      
//     Computations =                                                   
//     --------------                                                   
//                             [ [D_m]   [D_mb] ]                       
//     [Constitutive_Matrix] = [                ]                       
//            6 by 6           [ [D_bm]  [D_b]  ]                       
//                                                                      
//     where "b" and "m" stand for bending and membrane, respectively:  
//                                                                      
//     The constitutive matrix [D_b] relates the element's curvatures  
//     to the bending moments and the constitutive matrix [D_m]        
//     relates the element's normal efforts to the strains. Similarly,  
//     the constitutive matrices [D_bm] and [D_mb] couple the bending   
//     and membrane effects:                                            
//                                                                      
//     [ M_x  ]            [ k_x  ]                                     
//     [ M_y  ] = [D_b]  * [ k_y  ]                                     
//     [ M_xy ]            [ k_xy ]                                     
//                                                                      
//     [ N_x  ]            [ e_x  ]                                     
//     [ N_y  ] = [D_m]  * [ e_y  ]                                     
//     [ N_xy ]            [ e_xy ]                                     
//                                                                      
//     [ M_x  ]            [ e_x  ]                                     
//     [ M_y  ] = [D_bm] * [ e_y  ]                                     
//     [ M_xy ]            [ e_xy ]                                     
//                                                                      
//     [ N_x  ]            [ k_x  ]                                     
//     [ N_y  ] = [D_mb] * [ k_y  ]                                     
//     [ N_xy ]            [ k_xy ]                                     
//                                                                      
//     The symmetry constraints are:                                    
//                                                                      
//     [D_b]^T  = [D_b]                                                
//     [D_m]^T  = [D_m]                                                
//     [D_bm]^T = [D_mb]                                                
//                                                                      
//     The assembly of these matrices is performed according to the     
//     type of constitutive law.     
//     In the following, the assembly of each one of the four matrices  
//     is briefly summarized.                                           
//                                                                      
//     3.  Constitutive Law of type-2:                                  
//     - - - - - - - - - - - - - - - -                                  
//     The material properties of each layer are known and integration  
//     through the thickness of the composite material is performed.    
//     It is assumed that there is NO coupling between bending and      
//     membrane effects (even though these terms are found non-zero     
//     when numerical integration through the thickness is performed).  
//                                                                      
//     3.1  Pure Bending:                                        
//     - - - - - - - - - - - - -                                        
//                                                                      
//     3.2 Pure Membrane:                                        
//     - - - - - - - - - - - - -                                        
//                                                                      
//     3.3 Coupling Bending-Membrane:                            
//     - - - - - - - - - - - - - - - - - - -                            
//     [D_bm] = zero (assumed)                                          
//                                                                      
//     3.4 Coupling Membrane-Bending:                            
//     - - - - - - - - - - - - - - - - - - -                            
//     [D_mb] = zero (assumed)                                          
//                                                                      
//     4.  Constitutive Law of type-3:                                  
//     - - - - - - - - - - - - - - - -                                  
//     The material properties of each layer are known and integration  
//     through the thickness of the composite material is performed.    
//                                                                      
//     4.1  Pure Bending:                                        
//     - - - - - - - - - - - - -                                        
//                                                                      
//     4.2 Pure Membrane:                                        
//     - - - - - - - - - - - - -                                        
//                                                                      
//     4.3 Coupling Bending-Membrane:                            
//     - - - - - - - - - - - - - - - - - - -                            
//                                                                      
//     4.4 Coupling Membrane-Bending:                            
//     - - - - - - - - - - - - - - - - - - -                            
//                                                                      
//                                                                      
//     Caution =   It is assumed that the element has a constant        
//     ---------   thickness so that no numerical interpolation is      
//                 required.                                            
//                                                                      
// ==================================================================== 
// Author  = Francois M. Hemez                                          
// Date    = June 9th, 1995                                             
// Version = 2.0                                                        
// Comment =                                                            
// ==================================================================== 

    pi = acos(-one);

// .....CLEAR THE 3 BY 3 CONSTITUTIVE MATRIX 

    D.setZero();

// .....CALCULATE THE TOTAL HALF-HEIGHT OF THE LAYER 

    z0 = -0.5 * GetShellThickness();

// .....LOOP ON LAYERS OF THE COMPOSITE SHELL ELEMENT 

    for (ilayer = 0; ilayer < nlayer; ++ilayer) {

// .....INITIALIZE THE LAYER MATERIAL PROPERTIES 

        e1     = mtlayer(0, ilayer);
        e2     = mtlayer(1, ilayer);
        nu12   = mtlayer(2, ilayer);
        g12    = mtlayer(3, ilayer);
        mu1    = mtlayer(4, ilayer);
        mu2    = mtlayer(5, ilayer);
        thetaf = mtlayer(8, ilayer);

// .....CALCULATE THE COMPLIANCE MATRIX [S] WHICH RELATES THE STRESSES 
// .....[s1], [s2] AND [s12] TO THE STRAINS [e1], [e2] AND [e12] IN 
// .....THE COORDINATE SYSTEM {1;2} ASSOCIATED WITH THE FIBER ORIENTATION 

        s11 = 1 / e1;
        s12 = -nu12 / e1;
        s13 = mu1 / g12;
        s22 = 1 / e2;
        s23 = mu2 / g12;
        s33 = 1 / g12;

// .....CALCULATE THE DETERMINANT OF THE COMPLIANCE MATRIX 

        dets = s33 * (s11 * s22 - s12 * s12) - s11 * s23 * s23 - s22 *
                s13 * s13;
        dets += s12 * 2. * s13 * s23;

        if (dets == 0) {
            throw std::runtime_error(
                "*** FATAL ERROR in Routine COMPCST    ***\n"
                "*** The Compliance Matrix is Singular ***\n"
                "*** ... Check Material Properties ... ***\n"
                "*** STOP ALL TREATMENTS RIGHT HERE    ***\n");
            break;
        }

// .....CALCULATE THE INVERSE OF THE COMPLIANCE MATRIX (i.e. THE ELASTICITY 
// .....STIFFNESS MATRIX) WHICH RELATES THE STRAINS [e1], [e2] AND [e12] 
// .....TO THE STRESSES [s1], [s2] AND [s12] IN THE COORDINATE SYSTEM {1;2}
// .....OF THE FIBER ORIENTATION 

        C(0, 0) = s22 * s33 - s23 * s23;
        C(0, 1) = s13 * s23 - s12 * s33;
        C(0, 2) = s12 * s23 - s13 * s22;

        C(1, 0) = C(0, 1);
        C(1, 1) = s11 * s33 - s13 * s13;
        C(1, 2) = s12 * s13 - s11 * s23;

        C(2, 0) = C(0, 2);
        C(2, 1) = C(1, 2);
        C(2, 2) = s11 * s22 - s12 * s12;

        C /= dets;

// .....TRANSFORM ANGLE IN THE RANGE BETWEEN 0-360 

        if (thetaf < 0 || thetaf > 360.) {
            irot = (int) (thetaf / 360.);
            thetaf -= (doublereal) irot * 360.;
            if (thetaf < 0.) {
                thetaf += 360.;
            }
        }

// .....TRANSFORM FROM DEGREE TO RADIAN THE ANGLE BETWEEN THE 
// .....REFERENCE ORIENTATION VECTOR AND THE ORIENTATION OF THE FIBERS 

        thetaf = pi * thetaf / 180.;

// .....CALCULATE THE INVERSE OF THE COMPLIANCE MATRIX (i.e. THE ELASTICITY
// .....STIFFNESS MATRIX WHICH RELATES THE STRAINS [ex], [ey] AND [exy] 
// .....TO THE STRESSES [sx], [sy] AND [sxy] IN THE TRIANGULAR COORDINATE 
// .....SYSTEM {x;y}: 
// .....[C'] = [invT] * [C] * [invT]^t 

        invT = this->andesinvt(eframe, aframe, thetaf);
        C = invT * C * invT.transpose();

// .....INITIALIZE THE UPPER AND LOWER [z] COORDINATES FOR THE LAYER 

        zinf = z0;
        zsup = z0;

        for (i = 0; i < nlayer; ++i) {
            if (i < ilayer) {
                zinf += mtlayer(7, i);
            }
            if (i <= ilayer) {
                zsup += mtlayer(7, i);
            }
        }

//     -------------------------------------------------- 
//       (NUMERICAL INTEGRATION THROUGH THE THICKNESS)    
//     -------------------------------------------------- 

// .....ASSEMBLE THE CONSTITUTIVE MATRIX FOR PURE BENDING 

        Db += C * (zsup * zsup * zsup - zinf * zinf * zinf) / 3.;

// .....ASSEMBLE THE CONSTITUTIVE MATRIX FOR PURE MEMBRANE 

        Dm += C * (zsup - zinf);

        if(couple) {

// .....ASSEMBLE THE CONSTITUTIVE MATRIX FOR COUPLING BENDING-MEMBRANE 

            Dbm += C.transpose() * (zsup * zsup - zinf * zinf) * .5;

// .....ASSEMBLE THE CONSTITUTIVE MATRIX FOR COUPLING MEMBRANE-BENDING 

            Dmb += C * (zsup * zsup - zinf * zinf) * .5;
        }

    }

// .....COMPUTE THE GENERALIZED "STRESSES"

    Sigma = D*Upsilon; // this could alternatively be formed by integrating the stresses through the thickness

    if(_D == NULL) delete [] data;
}

template<typename doublereal>
void
ShellMaterialTypes2And3<doublereal>
::GetLocalConstitutiveResponse(doublereal *_Upsilon, doublereal *_sigma, doublereal z,
                               doublereal *eframe, int)
{
    // Initialized data 
    doublereal one = 1.;

    // Builtin functions 
    using std::acos;

    // Local variables 
    int i, irot, ilayer;
    doublereal c, s, e1, e2, z0, g12, s11, s12, s13, s22, pi, s23, s33, mu1, mu2, nu12;
    doublereal dets, thetaf, zsup, zinf;
    Eigen::Matrix<doublereal,3,3> invT;
    Eigen::Matrix<doublereal,3,1> epsilon;
    Eigen::Matrix<doublereal,3,3> C; // plane stress elasticity stiffness matrix
    Eigen::Map<Eigen::Matrix<doublereal,6,1> > Upsilon(_Upsilon);
    Eigen::Map<Eigen::Matrix<doublereal,3,1> > sigma(_sigma);

    // Some convenient definitions 
    Eigen::VectorBlock< Eigen::Map< Eigen::Matrix<doublereal,6,1> > >
        e = Upsilon.head(3), chi = Upsilon.tail(3);

// .....COMPUTE THE LOCAL STRAINS ON THE SPECIFIED SURFACE

    epsilon = e + z * chi;

    pi = acos(-one);

// .....LOOP ON LAYERS OF THE COMPOSITE SHELL ELEMENT 

    for (ilayer = 0; ilayer < nlayer; ++ilayer) {

// .....INITIALIZE THE UPPER AND LOWER [z] COORDINATES FOR THE LAYER 

      zinf = -0.5 * GetShellThickness();
      zsup = -0.5 * GetShellThickness();

      for (i = 0; i < nlayer; ++i) {
        if (i < ilayer) {
          zinf += mtlayer(7, i);
        }
        if (i <= ilayer) {
          zsup += mtlayer(7, i);
        }
      }

      zinf = (ilayer == 0) ? -std::numeric_limits<double>::infinity() : zinf;
      zsup = (ilayer == nlayer-1) ? std::numeric_limits<double>::infinity() : zsup;

// .....CHECK IF THIS IS THE LAYER CONTAINING z
      // TODO this doesn't deal with the case when z is exactly on the boundary between 2 layers
      // which can happen for example on the median surface

      if(!(z >= zinf && z <= zsup)) continue;

// .....INITIALIZE THE LAYER MATERIAL PROPERTIES 

      e1     = mtlayer(0, ilayer);
      e2     = mtlayer(1, ilayer);
      nu12   = mtlayer(2, ilayer);
      g12    = mtlayer(3, ilayer);
      mu1    = mtlayer(4, ilayer);
      mu2    = mtlayer(5, ilayer);
      thetaf = mtlayer(8, ilayer);

// .....CALCULATE THE COMPLIANCE MATRIX [S] WHICH RELATES THE STRESSES 
// .....[s1], [s2] AND [s12] TO THE STRAINS [e1], [e2] AND [e12] IN 
// .....THE COORDINATE SYSTEM {1;2} ASSOCIATED WITH THE FIBER ORIENTATION 

      s11 = 1 / e1;
      s12 = -nu12 / e1;
      s13 = mu1 / g12;
      s22 = 1 / e2;
      s23 = mu2 / g12;
      s33 = 1 / g12;

// .....CALCULATE THE DETERMINANT OF THE COMPLIANCE MATRIX 

      dets = s33 * (s11 * s22 - s12 * s12) - s11 * s23 * s23 - s22 *
              s13 * s13;
      dets += s12 * 2. * s13 * s23;

      if (dets == 0) {
        throw std::runtime_error(
          "*** FATAL ERROR in Routine COMPCST    ***\n"
          "*** The Compliance Matrix is Singular ***\n"
          "*** ... Check Material Properties ... ***\n"
          "*** STOP ALL TREATMENTS RIGHT HERE    ***\n");
        break;
      }

// .....CALCULATE THE INVERSE OF THE COMPLIANCE MATRIX (i.e. THE ELASTICITY 
// .....STIFFNESS MATRIX) WHICH RELATES THE STRAINS [e1], [e2] AND [e12] 
// .....TO THE STRESSES [s1], [s2] AND [s12] IN THE COORDINATE SYSTEM {1;2}
// .....OF THE FIBER ORIENTATION 

      C(0, 0) = s22 * s33 - s23 * s23;
      C(0, 1) = s13 * s23 - s12 * s33;
      C(0, 2) = s12 * s23 - s13 * s22;

      C(1, 0) = C(0, 1);
      C(1, 1) = s11 * s33 - s13 * s13;
      C(1, 2) = s12 * s13 - s11 * s23;

      C(2, 0) = C(0, 2);
      C(2, 1) = C(1, 2);
      C(2, 2) = s11 * s22 - s12 * s12;

      C /= dets;

// .....TRANSFORM ANGLE IN THE RANGE BETWEEN 0-360 

      if (thetaf < 0 || thetaf > 360.) {
        irot = (int) (thetaf / 360.);
        thetaf -= (doublereal) irot * 360.;
        if (thetaf < 0.) {
          thetaf += 360.;
        }
      }

// .....TRANSFORM FROM DEGREE TO RADIAN THE ANGLE BETWEEN THE 
// .....REFERENCE ORIENTATION VECTOR AND THE ORIENTATION OF THE FIBERS 

      thetaf = pi * thetaf / 180.;

// .....CALCULATE THE INVERSE OF THE COMPLIANCE MATRIX (i.e. THE ELASTICITY
// .....STIFFNESS MATRIX WHICH RELATES THE STRAINS [ex], [ey] AND [exy] 
// .....TO THE STRESSES [sx], [sy] AND [sxy] IN THE TRIANGULAR COORDINATE 
// .....SYSTEM {x;y}: 
// .....[C'] = [invT] * [C] * [invT]^t 

      invT = this->andesinvt(eframe, aframe, thetaf);
      C = invT * C * invT.transpose();

      sigma = C*epsilon;
      break;
    }
}

template
ShellMaterialTypes2And3<double>::ShellMaterialTypes2And3(
  int nlayer, double *mtlayer, bool couple, double *aframe);

template
void
ShellMaterialTypes2And3<double>::GetConstitutiveResponse(double *Upsilon, double *Sigma, double *D,
                                                         double *eframe, int gp);

template
void
ShellMaterialTypes2And3<double>
::GetLocalConstitutiveResponse(double *Upsilon, double *sigma, double z, double *eframe, int);
#endif
