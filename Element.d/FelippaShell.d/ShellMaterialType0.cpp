#ifdef USE_EIGEN3
#include <cmath>
#include <stdexcept>
#include <Element.d/FelippaShell.d/ShellMaterial.hpp>

template<typename doublereal>
void
ShellMaterialType0<doublereal>::GetConstitutiveResponse(doublereal *_Upsilon, doublereal *_Sigma, doublereal *_D,
                                                        doublereal *, int)
{
  // Initialized data 
  doublereal zero = 0.;
  doublereal one = 1.;

  // Local variables 
  doublereal *data = (_D == NULL) ? new doublereal[36] : _D;
  Eigen::Map<Eigen::Matrix<doublereal,6,1> > Upsilon(_Upsilon), Sigma(_Sigma);
  Eigen::Map<Eigen::Matrix<doublereal,6,6> > D(data); 
  Eigen::Block< Eigen::Map<Eigen::Matrix<doublereal,6,6> > >
    Dm = D.topLeftCorner(3,3),     Dmb = D.topRightCorner(3,3),
    Dbm = D.bottomLeftCorner(3,3), Db = D.bottomRightCorner(3,3);

// ==================================================================== 
//                                                                      
//     Perform =   Assembles the 6 by 6 Constitutive Matrix According   
//     ---------   to the Type of Constitutive Law Requested.           
//                                                                      
//                                                                      
//     Input/Output =                                                   
//     --------------                                                   
//     E       <input>  Young modulus                                   
//     thick   <input>  thickness (assumed constant over the element)   
//     nu      <input>  Poisson's ratio                                 
//     D       <output> 6 by 6 constitutive matrix                      
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
//     1.  Constitutive Law of type-0:                                  
//     - - - - - - - - - - - - - - - -                                  
//     The material is assumed isotropic and known via the Young        
//     modulus [E], the Poisson's ratio [nu] and thickness [thick]:     
//                                                                      
//     1.1  Pure Bending:                                        
//     - - - - - - - - - - - - -                                        
//                                                                      
//             [ d_11  d_12   0   ]                                    
//     [D_b] = [ d_12  d_11   0   ]                                    
//             [  0     0    d_33 ]                                    
//                                                                      
//     with:                                                            
//                                                                      
//              [E]*[thick]^3                                           
//     [d_11] = -------------                                           
//              12*(1-[nu]^2)                                           
//                                                                      
//              [nu]*[E]*[thick]^3                                      
//     [d_12] = ------------------                                      
//                12*(1-[nu]^2)                                         
//                                                                      
//              [E]*[thick]^3                                           
//     [d_33] = -------------                                           
//               24*(1+[nu])                                            
//                                                                      
//     1.2 Pure Membrane:                                        
//     - - - - - - - - - - - - -                                        
//                                                                      
//             [ d_44  d_45   0   ]                                    
//     [D_m] = [ d_45  d_55   0   ]                                    
//             [  0     0    d_66 ]                                    
//                                                                      
//     with:                                                            
//                                                                      
//              [E]*[thick]                                             
//     [d_44] = -----------                                             
//              (1-[nu]^2)                                              
//                                                                      
//              [nu]*[E]*[thick]                                        
//     [d_45] = ----------------                                        
//                 (1-[nu]^2)                                           
//                                                                      
//              [E]*[thick]                                             
//     [d_66] = ------------                                            
//               2*(1+[nu])                                             
//                                                                      
//     1.3 Coupling Bending-Membrane:                            
//     - - - - - - - - - - - - - - - - - - -                            
//     [D_bm] = zero (no coupling for isotropic material)               
//                                                                      
//     1.4 Coupling Membrane-Bending:                            
//     - - - - - - - - - - - - - - - - - - -                            
//     [D_mb] = zero (no coupling for isotropic material)               
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

// .....ASSEMBLE THE CONSTITUTIVE MATRIX FOR PURE BENDING

    Db(0, 0) = E * (thick * thick * thick) / ((one - nu * nu) * 12.);
    Db(0, 1) = nu * E * (thick * thick * thick) / ((one - nu * nu) * 12.);
    Db(0, 2) = zero;
    Db(1, 0) = nu * E * (thick * thick * thick) / ((one - nu * nu) * 12.);
    Db(1, 1) = E * (thick * thick * thick) / ((one - nu * nu) * 12.);
    Db(1, 2) = zero;
    Db(2, 0) = zero;
    Db(2, 1) = zero;
    Db(2, 2) = E * (thick * thick * thick) / ((one + nu) * 24.);

// .....ASSEMBLE THE CONSTITUTIVE MATRIX FOR PURE MEMBRANE 
    Dm(0, 0) = E * thick / (one - nu * nu);
    Dm(0, 1) = nu * E * thick / (one - nu * nu);
    Dm(0, 2) = zero;
    Dm(1, 0) = nu * E * thick / (one - nu * nu);
    Dm(1, 1) = E * thick / (one - nu * nu);
    Dm(1, 2) = zero;
    Dm(2, 0) = zero;
    Dm(2, 1) = zero;
    Dm(2, 2) = E * thick / ((one + nu) * 2.);

// .....ASSEMBLE THE CONSTITUTIVE MATRIX FOR COUPLING BENDING-MEMBRANE

    Dbm = Eigen::Matrix<doublereal,3,3>::Zero();

// .....ASSEMBLE THE CONSTITUTIVE MATRIX FOR COUPLING MEMBRANE-BENDING

    Dmb = Eigen::Matrix<doublereal,3,3>::Zero();

// .....COMPUTE THE GENERALIZED "STRESSES"

    Sigma = D*Upsilon;

    if(_D == NULL) delete [] data;
}

template<typename doublereal>
void
ShellMaterialType0<doublereal>
::GetLocalConstitutiveResponse(doublereal *_Upsilon, doublereal *_sigma, doublereal z,
                               doublereal *, int)
{
    // Local variables
    Eigen::Matrix<doublereal,3,1> epsilon;
    Eigen::Matrix<doublereal,3,3> C;
    Eigen::Map<Eigen::Matrix<doublereal,6,1> > Upsilon(_Upsilon);
    Eigen::Map<Eigen::Matrix<doublereal,3,1> > sigma(_sigma);

    // Some convenient definitions 
    Eigen::VectorBlock< Eigen::Map< Eigen::Matrix<doublereal,6,1> > >
        e = Upsilon.head(3), chi = Upsilon.tail(3);

// .....COMPUTE THE LOCAL STRAINS [epsilon] = {epsilonxx,epsilonyy,gammaxy} ON THE SPECIFIED SURFACE

    epsilon = e + z * chi;

// .....GET THE PLANE STRESS ELASTICITY STIFFNESS MATRIX

    doublereal v = E/(1-nu*nu);
    C << v,    v*nu, 0,
         v*nu, v,    0,
         0,    0.,   v*(1-nu)/2;

// .....COMPUTE THE LOCAL STRESSES [sigma] = {sigmaxx,sigmayy,sigmaxy} ON THE SPECIFIED SURFACE

    sigma = C*epsilon;
}

template
void
ShellMaterialType0<double>
::GetConstitutiveResponse(double *Upsilon, double *Sigma, double *D, double *, int);

template
void
ShellMaterialType0<double>
::GetLocalConstitutiveResponse(double *Upsilon, double *sigma, double z, double *, int);
#endif