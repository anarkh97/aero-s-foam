#ifdef USE_EIGEN3
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <Element.d/FelippaShell.d/ShellElementTemplate.hpp>
#include <Element.d/FelippaShell.d/ShellMaterial.hpp>
#include <Eigen/Dense>
#include <Eigen/Geometry>

template<typename doublereal, template<typename> class Membrane, template<typename> class Bending>
void
ShellElementTemplate<doublereal,Membrane,Bending>
::andescrd(int elm, doublereal *_x, doublereal *_y, doublereal *_z,
           doublereal *_eframe, doublereal *_xlp, 
           doublereal *_ylp, doublereal *_zlp, doublereal &area)
{
  // Builtin functions 
  using std::sqrt;
  using std::abs;

  // Local variables 
  int i, j;
  doublereal x21, y21, z21, x13, y13, z13, x32, y32, z32, signedarea,
         projection, xcg, ycg, zcg, side21length, 
         side32length, lengthy, lengthz;
  Eigen::Map<Eigen::Matrix<doublereal,3,1> > x(_x), y(_y), z(_z), xlp(_xlp), ylp(_ylp), zlp(_zlp);
  Eigen::Matrix<doublereal,3,3> xyz; xyz << x[0], x[1], x[2],
                                            y[0], y[1], y[2],
                                            z[0], z[1], z[2];
  Eigen::Matrix<doublereal,3,1> side21, side13, side32, cg, xp, yp, zp;
  Eigen::Map<Eigen::Matrix<doublereal,3,3> > eframe(_eframe);

// ==================================================================== 
//                                                                      
//     Perform =    This subroutine computes basic quantities needed    
//     ---------    for the assembly of the basic and higher order      
//                  composite stiffness matrices.                       
//                                                                      
//                                                                      
//     Inputs/Outputs =                                                 
//     ----------------                                                 
//     ELM    <input>   finite element number                           
//     X      <input>   nodal coordinates in the X-direction            
//     Y      <input>   nodal coordinates in the Y-direction            
//     Z      <input>   nodal coordinates in the Z-direction            
//     EFRAME <output>  element frame                                   
//     XLP    <output>  triangular coordinates in the X-direction       
//     YLP    <output>  triangular coordinates in the Y-direction       
//     ZLP    <output>  triangular coordinates in the Z-direction       
//                                                                      
// ==================================================================== 
// Author  = Francois M. Hemez                                          
// Date    = June 9th, 1994                                             
// Version = 1.0                                                        
// Comment =                                                            
// ==================================================================== 

// .....COMPUTE THE NODAL COORDINATE DIFFERENCES 

    side21 = xyz.col(1)-xyz.col(0);
    side13 = xyz.col(0)-xyz.col(2);
    side32 = xyz.col(2)-xyz.col(1);

// .....COMPUTE THE LENGTH OF SIDE 2-1 

    side21length = side21.norm();

// .....CHECK IF LENGTH 2-1 IS DIFFERENT FROM ZERO 

    if (side21length == 0) {
        throw std::runtime_error(
          "*** FATAL ERROR in ShellElementTemplate::andescrd ***\n"
          "*** Side Between Nodes 1 and 2 Has 0-Length       ***\n"
          "*** Check Coordinates and FE Topology             ***\n");
    }

// .....COMPUTE THE LENGTH OF SIDE 3-2 

    side32length = side32.norm();

// .....CHECK IF LENGTH 3-2 IS DIFFERENT FROM ZERO 

    if (side32length == 0) {
        throw std::runtime_error(
          "*** FATAL ERROR in ShellElementTemplate::andescrd ***\n"
          "*** Side Between Nodes 2 and 3 Has 0-Length       ***\n"
          "*** Check Coordinates and FE Topology             ***\n");
    }

// .....COMPUTE THE DISTANCE OF THE OPPOSING NODE 3 TO SIDE 2-1 

    projection = abs(side21.dot(side32))/side21length;

// .....GET THE AREA OF THE TRIANGLE 

    signedarea = side32length * side32length - projection * projection;

    if (signedarea <= 0) {
        throw std::runtime_error(
          "*** FATAL ERROR in ShellElementTemplate::andescrd ***\n"
          "*** The Area is Negative or Zero                  ***\n"
          "*** Check Coordinates and FE Topology             ***\n");
    }

    area = side21length * .5 * sqrt(signedarea);

// .....COMPUTE THE DIRECTION COSINES OF THE LOCAL SYSTEM 
// .....DIRECTION [X] IS DIRECTED PARALLEL TO THE SIDE 2-1 
// .....DIRECTION [Z] IS THE EXTERNAL NORMAL (COUNTERCLOCKWISE) 
// .....DIRECTION [Y] IS COMPUTED AS [Z] x [X] (TENSORIAL PRODUCT) 

    xp = side21.normalized();
    zp = (side21.cross(side32)).normalized();
    yp = (zp.cross(xp)).normalized();

// .....COMPUTE THE COORDINATES FOR THE CENTER OF GRAVITY 

    cg = xyz.rowwise().sum()/3;

// .....CONSTRUCT THE ELEMENT FRAME 

    eframe << xp, yp, zp;

// .....COMPUTE THE LOCAL COORDINATES 

    Eigen::Matrix<doublereal,3,3> lp = (xyz-cg.replicate(1,3)).transpose()*eframe;
    xlp = lp.col(0);
    ylp = lp.col(1);
    zlp = lp.col(2); 
}

template<typename doublereal, template<typename> class Membrane, template<typename> class Bending>
void
ShellElementTemplate<doublereal,Membrane,Bending>
::andesms(int elm, doublereal *x, doublereal *y, doublereal *z, 
          doublereal *_emass, doublereal *gamma, doublereal *grvfor, 
          bool grvflg, doublereal &totmas, bool masflg)
{
  // Reference:
  // "Parametrized variational principles in dynamics applied
  //  to the optimization of dynamic models of plates",
  // F. J. Brito Castro, C. Militello, C. A. Felippa
  // Computational Mechanics 20 (1997) 285±294
  // Notes: this is the LMM (lumped mass matrix) which I think is formed by row lumping of the BCMM
  // TODO: implement BCMM (boundary consistent mass matrix)
  //       and CMM (consistent mass matrix)


  // Builtin functions 
  using std::sqrt;
  using std::abs;

  // Local variables 
  int i, j, i1, i2, i3;
  doublereal twicearea2, x21, x13, y13, z13, x32, y32, z32, y21, z21,
    ix, iy, iz, rlb, bpr, rlr, area, rhoh, dist[3], mass0, mass1, 
    mass2, mass3, thick;

  Eigen::Map<Eigen::Matrix<doublereal,18,18,Eigen::RowMajor> > emass(_emass); 

// ==================================================================== 
//                                                                      
//     Performs =   This subroutine will form the elemental mass        
//     ----------   matrix of the 3D 3-node ANDES composite shell.      
//                  Lumping is assumed here.                            
//                                                                      
//                                                                      
//     Inputs/Outputs =                                                 
//     ----------------                                                 
//     X        <input>   nodal coordinates in the X-direction          
//     Y        <input>   nodal coordinates in the Y-direction          
//     Z        <input>   nodal coordinates in the Z-direction          
//                                                                      
//     Computations =                                                   
//     ---------------                                                  
//                                                                      
//     The lumped mass matrix [M] is equal to:                          
//                                                                      
//           [ mt 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 ]   
//           [ 0  mt 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 ]   
//           [ 0  0  mt 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 ]   
//           [ 0  0  0  m1 0  0  0  0  0  0  0  0  0  0  0  0  0  0 ]   
//           [ 0  0  0  0  m1 0  0  0  0  0  0  0  0  0  0  0  0  0 ]   
//           [ 0  0  0  0  0  m1 0  0  0  0  0  0  0  0  0  0  0  0 ]   
//           [ 0  0  0  0  0  0  mt 0  0  0  0  0  0  0  0  0  0  0 ]   
//           [ 0  0  0  0  0  0  0  mt 0  0  0  0  0  0  0  0  0  0 ]   
//           [ 0  0  0  0  0  0  0  0  mt 0  0  0  0  0  0  0  0  0 ]   
//     [M] = [ 0  0  0  0  0  0  0  0  0  m2 0  0  0  0  0  0  0  0 ]   
//           [ 0  0  0  0  0  0  0  0  0  0  m2 0  0  0  0  0  0  0 ]   
//           [ 0  0  0  0  0  0  0  0  0  0  0  m2 0  0  0  0  0  0 ]   
//           [ 0  0  0  0  0  0  0  0  0  0  0  0  mt 0  0  0  0  0 ]   
//           [ 0  0  0  0  0  0  0  0  0  0  0  0  0  mt 0  0  0  0 ]   
//           [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  mt 0  0  0 ]   
//           [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  m3 0  0 ]   
//           [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  m3 0 ]   
//           [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  m3]   
//                                                                      
//     with the following ordering of local degrees of freedom:         
//                                                                      
//                                [     U_x1 ]                          
//                                [     U_y1 ]                          
//                                [     U_z1 ]                          
//                                [ theta_x1 ]                          
//                                [ theta_y1 ]                          
//                                [ theta_z1 ]                          
//                                [     U_x2 ]                          
//                                [     U_y2 ]                          
//                                [     U_z2 ]                          
//     Ordering_of_Local_DOFs  =  [ theta_x2 ]                          
//                                [ theta_y2 ]                          
//                                [ theta_z2 ]                          
//                                [     U_x3 ]                          
//                                [     U_y3 ]                          
//                                [     U_z3 ]                          
//                                [ theta_x3 ]                          
//                                [ theta_y3 ]                          
//                                [ theta_z3 ]                          
//                                                                      
//     No rotation of local-to-global basis is implemented since the    
//     mass matrix [M] is formed of 3 by 3 blocks proportional to the   
//     identity. The lumping factors are equal to:                      
//                                                                      
//     [mt] = [rhoh] * [A]        /    3.0                              
//     [m1] = [rhoh] * [A] * [Ix] / 1260.0                              
//     [m2] = [rhoh] * [A] * [Iy] / 1260.0                              
//     [m3] = [rhoh] * [A] * [Iz] / 1260.0                              
//                                                                      
//     where:                                                           
//                                                                      
//     [rhoh]                area density                               
//     [A]                   area                                       
//     [Ix], [Iy] and [Iz]   equivalent pseudo-moments of inertia       
//                                                                      
//     Caution =                                                        
//     ---------                                                        
//     The finite element is assumed to have constant thickness so      
//     that no numerical interpolation is required.                     
//                                                                      
// ==================================================================== 
// Author  = Francois M. Hemez                                          
// Date    = June 9, 1995                                               
// Version = 2.0                                                        
// ==================================================================== 

// .....CLEAR THE OUTPUT MASS MATRIX

    emass.setZero();

//     -------------------------- 
//     CHECKS AND INITIALIZATIONS 
//     -------------------------- 

// .....COMPUTE THE DISTANCE BETWEEN X-, Y- AND Z- NODAL COORDINATES

    x21 = x[1] - x[0];
    y21 = y[1] - y[0];
    z21 = z[1] - z[0];

    x32 = x[2] - x[1];
    y32 = y[2] - y[1];
    z32 = z[2] - z[1];

    x13 = x[0] - x[2];
    y13 = y[0] - y[2];
    z13 = z[0] - z[2];

// .....COMPUTE THE DISTANCE BETWEEN NODES 1-2, 2-3 AND 3-1

    dist[0] = sqrt(x21 * x21 + y21 * y21 + z21 * z21);
    dist[1] = sqrt(x32 * x32 + y32 * y32 + z32 * z32);
    dist[2] = sqrt(x13 * x13 + y13 * y13 + z13 * z13);

// .....COMPUTE THE LENGTH OF SIDE 1-2

    rlr = sqrt(x21 * x21 + y21 * y21 + z21 * z21);

// .....CHECK FOR ZERO-SIDE LENGTH

    if (rlr == 0) {
        throw std::runtime_error(
          "*** FATAL ERROR in ShellElementTemplate::andesms ***\n"
          "*** The Side 1-2 has Zero Length                 ***\n"
          "*** Check Coordinates and FE Topology            ***\n");
    }

// .....COMPUTE THE DISTANCE OF THE OPPOSING NODE (3) TO THAT SIDE (1-2)

    rlb = sqrt(x32 * x32 + y32 * y32 + z32 * z32);
    bpr = abs(x21 * x32 + y21 * y32 + z21 * z32) / rlr;

// .....COMPUTE THE SQUARE OF TWICE THE TRIANGLE'S AREA

    twicearea2 = rlb * rlb - bpr * bpr;

// .....CHECK IF THE TRIANGLE'S AREA IS POSITIVE

    if (twicearea2 <= 0) {
        throw std::runtime_error(
          "*** FATAL ERROR in ShellElementTemplate::andesms ***\n"
          "*** The Area is Negative or Zero                 ***\n"
          "*** Check Coordinates and FE Topology            ***\n");
    }

// .....COMPUTE THE AREA OF THE TRIANGLE

    area = rlr * .5 * sqrt(twicearea2);

// .....COMPUTE THE THREE PSEUDO MOMENTS OF INERTIA

    ix = dist[0] * dist[0] + dist[2] * dist[2];
    iy = dist[0] * dist[0] + dist[1] * dist[1];
    iz = dist[1] * dist[1] + dist[2] * dist[2];

// .....FORM THE MASS COEFFICIENTS PER DEGREE OF FREEDOM

    rhoh = gpmat->GetAreaDensity();
    mass0 = rhoh * area / 3.;
    mass1 = rhoh * area * ix / 1260.;
    mass2 = rhoh * area * iy / 1260.;
    mass3 = rhoh * area * iz / 1260.;

//     ------------------------------------- 
//     ASSEMBLY OF THE ELEMENTAL MASS MATRIX 
//     ------------------------------------- 

// .....FORM THE LUMPED ELEMENT MASS MATRIX 

    for (i = 0; i < 3; ++i) {
        i2 = i + 6;
        i3 = i + 12;
        emass(i, i) = mass0;
        emass(i2, i2) = mass0;
        emass(i3, i3) = mass0;
    }

    for (i = 0; i < 3; ++i) {
        i1 = i + 3;
        i2 = i + 9;
        i3 = i + 15;
        emass(i1, i1) = mass1;
        emass(i2, i2) = mass2;
        emass(i3, i3) = mass3;
    }

// ..... COMPUTE THE BODY FORCE DUE TO GRAVITY IF NEEDED

    if (grvflg) {
        grvfor[0] = mass0 * 3. * gamma[0];
        grvfor[1] = mass0 * 3. * gamma[1];
        grvfor[2] = mass0 * 3. * gamma[2];
    }

// .... ACCUMULATE THE SUBDOMAIN MASS

    if (masflg) {
        totmas += mass0 * 3.;
    }

}


template<typename doublereal, template<typename> class Membrane, template<typename> class Bending>
void
ShellElementTemplate<doublereal,Membrane,Bending>
::andesstf(int elm, doublereal *_estiff, doublereal *_fint, doublereal nu,
           doublereal *x, doublereal *y, doublereal *z, doublereal *_v,
           int ctyp, int flag)
{
  // Initialized data 
  bool debug = false;
  doublereal clr = 0;
  doublereal cqr = 1;
  doublereal betab = 1;
  doublereal alpha = 1.5;
  doublereal betam = .32;

  // Local variables 
  int i, j;
  doublereal xlp[3], ylp[3], zlp[3];
  doublereal area;

  Eigen::Matrix<doublereal,9,3> Lb, Lm;
  Eigen::Matrix<doublereal,3,9> Bb, Bm;

  Eigen::Map< Eigen::Matrix<doublereal,18,18,Eigen::RowMajor> > K(_estiff);
  Eigen::Map< Eigen::Matrix<doublereal,18,1> > F(_fint);
  Eigen::Matrix<doublereal,6,1> Upsilon, Sigma;
  Eigen::Matrix<doublereal,6,6> *D = (_estiff) ? new Eigen::Matrix<doublereal,6,6> : NULL;
  Eigen::Matrix<doublereal,3,3> eframe;
  Eigen::Matrix<doublereal,18,1> vd;
  Eigen::Map<Eigen::Matrix<doublereal,18,1> > v(_v);

  // Some convenient definitions 
  Eigen::Block< Eigen::Map<Eigen::Matrix<doublereal,18,18,Eigen::RowMajor> >,9,9 >
    Km = K.template topLeftCorner<9,9>(),     Kmb = K.template topRightCorner<9,9>(),
    Kbm = K.template bottomLeftCorner<9,9>(), Kb = K.template bottomRightCorner<9,9>();
  Eigen::VectorBlock< Eigen::Map<Eigen::Matrix<doublereal,18,1> >,9 >
    Fm = F.template head<9>(), Fb = F.template tail<9>();
  Eigen::Block< Eigen::Matrix<doublereal,6,6>,3,3 > 
    Dm = D->template topLeftCorner<3,3>(),     Dmb = D->template topRightCorner<3,3>(),
    Dbm = D->template bottomLeftCorner<3,3>(), Db = D->template bottomRightCorner<3,3>();
  Eigen::VectorBlock< Eigen::Matrix<doublereal,6,1>,3 >
    e = Upsilon.template head<3>(), chi = Upsilon.template tail<3>();
  Eigen::VectorBlock< Eigen::Matrix<doublereal,6,1>,3 >
    N = Sigma.template head<3>(), M = Sigma.template tail<3>();

// ================================================================== 
//                                                                    
//     Perform =    This subroutine will form the element stiffness   
//     ---------    matrix for the 3D 3-node ANDES-EFF shell element.  
//                                                                    
//                                                                    
//     Inputs/Outputs =                                               
//     ----------------                                               
//     ELM     <input>  finite element number                         
//     ESTIFF  <output> element stiffness matrix                      
//     NU      <input>  Poisson's ratio                               
//     X       <input>  nodal coordinates in the X-direction          
//     Y       <input>  nodal coordinates in the Y-direction          
//     Z       <input>  nodal coordinates in the Z-direction          
//     CTYP    <input>  composite attribute number                    
//     FLAG    <input>  int specifying whether to return              
//                      transformed element stiffness matrix or       
//                      global element stiffness matrix               
//                                                                    
//                                                                    
//     Computations =                                                 
//     --------------                                                 
//     This subroutine evaluates the stiffness matrix for the 18      
//     degrees of freedom 3-node composite triangle obtained as a     
//     combination of the Assumed Quadratic Rotations bending         
//     triangle plus the EFF membrane with driling degrees of freedom     
//     developed by Militello, Felippa et al. For documentation, see     
//     "The first ANDES elements: 9-dof plate bending triangles",
//      Militello & Felippa, Comput. Methods Appl. Mech. Engrg,
//      Vol. 93 (1991) pp. 217-246
//     "Membrane triangles with corner drilling freedoms I. The EFF element"
//      Alvin, de la Fuente, Haugen & Felippa,
//      Finite Elements in Analysis and Design Vol. 12 (1992) pp. 163-187
//     
//                                                                    
//     The original version of the ANS shell element has been         
//     generalized here to the case of a composite shell element      
//     with complete bending-membrane coupling and/or layered
//     plastic shell with basic-higher order coupling. Five types of         
//     constitutive laws have been used:                    
//        type-0: isotropic linear                              
//        type-1: constitutive coefficients are given                 
//        type-2: properties of each layer are given and no coupling  
//                   is assumed between bending and membrane          
//        type-3: properties of each layer are given and coupling     
//                   between bending and membrane is assumed          
//        type-4: isotropic nonlinear
//                                                                    
//     The stiffness matrix [K] is assembled as the combination of    
//     the basic stiffness and the higher order stiffness matrices, see
//     "Plastic buckling and collapse of thin shell structures,
//      using layered plastic modeling and co-rotational ANDES
//      finite element", Dal Cortivo, Felippa, Bavestrello & Silva,
//      Comput. Methods Appl. Mech. Engrg, Vol. 198 (2009) pp. 785-798.
//                                                                    
//                                                                    
//     Caution =                                                      
//     ---------                                                      
//     The finite element is assumed to have a constant thickness     
//     so that no numerical interpolation is required.
//                                                                    
// ================================================================== 
// Authors = Francois M. Hemez and Philip J. S. Avery                                       
// Date    = September 13, 2011                                             
// Version = 3.0                                                      
// ================================================================== 

// .....GET THE ELEMENT TRIANGULAR COORDINATES 
// .....GET THE ELEMENT LEVEL FRAME

    andescrd(elm, x, y, z, eframe.data(), xlp, ylp, zlp, area);

// .....CONSTRUCT THE PERMUTATION MATRIX 

    Eigen::Matrix<int,18,1> indices;
    indices << 0, 1, 6, 7, 12, 13, 5, 11, 17, // M indices
               2, 3, 4, 8, 9, 10, 14, 15, 16; // B indices
    Eigen::PermutationMatrix<18,18,int> PI(indices);

//  .....ROTATE THE NODAL DISPLACEMENTS TO THE LOCAL 
//       FRAME SYSTEM (LOCAL TO THE SHELL ELEMENT) 
//       AND APPLY PERMUTATION to {M,B} ordering

    if(flag == 1) {

        for(i = 0; i < 18; i += 3)
            vd.segment(i,3) = eframe.transpose()*v.segment(i,3);

        vd = PI.transpose()*vd;
    }
    else {

        // for Corotational Formulation, v is already local
        vd = PI.transpose()*v; 

    }

    // Note: betam = 0.32 is max(1/2*(1-4*nu^2),0.01) assuming nu to be 0.3
    // Reference: "Membrane triangles with corner drilling freedoms III. Implementation and performance evaluation"
    //             Carlos A. Felippa and Scott Alexander
    // It's not clear what to do when nu is not available (e.g. for type 1 composite) or varies through the thickness
    // (e.g. for types 2 and 3 composite)

// .....FORM THE BASIC CURVATURE-TO-DISPLACEMENT MATRIX
 
    Lb = Bending<doublereal>::L(xlp, ylp, clr, cqr);

// .....FORM THE BASIC EXTENSION-TO-DISPLACEMENT MATRIX

    Lm = Membrane<doublereal>::L(xlp, ylp, alpha);

    if(_estiff) K.setZero();
    if(_fint) F.setZero();
    doublereal zeta[3][3] = { { 0.,.5,.5 }, { .5,0.,.5 }, { .5,.5,0. } }; // triangular coordinates of gauss integration points
    doublereal weight[3] = { 1/3., 1/3., 1/3. };
    for(i = 0; i < 3; ++i) {

// .....FORM HIGHER ORDER INTEGRATED CURVATURE-TO-DISPLACEMENT MATRIX
//      AND THE ELEMENT CURVATURE VECTOR 
        Bb = 1/area*Lb.transpose() + Bending<doublereal>::Bd(xlp, ylp, betab, zeta[i]);
        chi = Bb*vd.tail(9);

// .....FORM THE HIGHER ORDER INTEGRATED EXTENSION-TO-DISPLACEMENT MATRIX
//      AND THE ELEMENT EXTENSION VECTOR

        Bm = 1/area*Lm.transpose() + Membrane<doublereal>::Bd(xlp, ylp, betam, zeta[i]);
        e = Bm*vd.head(9);

// .....GET THE TANGENT CONSTITUTIVE MATRIX [D = {Dm,Dmb;Dbm,Db}] 
//      AND THE GENERALIZED STRESSES [Sigma = {N,M}]

        if(i == 0 || ctyp == 4) {
          doublereal *_D = (_estiff) ? D->data() : NULL;
          gpmat->GetConstitutiveResponse(Upsilon.data(), Sigma.data(), _D, eframe.data(), i);
        }

        if(_estiff) {

// .....FORM STIFFNESS FOR PURE BENDING

            Kb.noalias() += (area*weight[i])*(Bb.transpose()*Db*Bb);

// .....FORM STIFFNESS FOR PURE MEMBRANE 

            Km.noalias() += (area*weight[i])*(Bm.transpose()*Dm*Bm);

            if(ctyp != 0 && ctyp != 2) {

// .....FORM STIFFNESS FOR BENDING-MEMBRANE COUPLING

              Kbm.noalias() += (area*weight[i])*(Bb.transpose()*Dbm*Bm);

// .....FORM STIFFNESS FOR MEMBRANE-BENDING COUPLING

              Kmb.noalias() += (area*weight[i])*(Bm.transpose()*Dmb*Bb);

            }

// .....CHECK THE CONSTITUTIVE MATRIX (USED FOR DEBUGGING ONLY)

            if(debug && (i == 0 || ctyp == 4)) {
                std::cerr << "Here are the eigenvalues of the constitutive matrix (element " << elm 
                          << ", gauss point " << i << ") :\n"
                          << D->eigenvalues().transpose() << std::endl;
            }

        }

        if(_fint) {

// ..... FORM THE INTERNAL FORCE FOR BENDING

            Fb.noalias() += (area*weight[i])*Bb.transpose()*M;

// ..... FORM THE INTERNAL FORCE FOR MEMBRANE

            Fm.noalias() += (area*weight[i])*Bm.transpose()*N;

        }

    }

// .....APPLY PERMUTATION 

    if(_estiff) K = PI*K*PI.transpose();
    if(_fint)   F = PI*F;

// .....ROTATE ELEMENT STIFFNESS AND/OR FORCE TO GLOBAL COORDINATES 
//      (only if flag equals 1)

    if (flag == 1 && _estiff) {
        for(i = 0; i < 18; i += 3)
            for(j = 0; j < 18; j += 3)
                 K.template block<3,3>(i,j) = eframe*K.template block<3,3>(i,j)*eframe.transpose();
    }

    if (flag == 1 && _fint) {
        for(i = 0; i < 18; i += 3)
            F.template segment<3>(i) = eframe*F.template segment<3>(i);
    }

// .....CHECK THE POSITIVITY OF THE OUTPUT STIFFNESS MATRIX 
// .....(USED FOR DEBUGGING ONLY) 

    if(debug && _estiff) {
        std::cerr << "Here are the eigenvalues of the stiffness matrix (element " << elm << ") :\n"
                  << K.eigenvalues().transpose() << std::endl;
    }

    if(_estiff) delete D;
}


template<typename doublereal, template<typename> class Membrane, template<typename> class Bending>
void
ShellElementTemplate<doublereal,Membrane,Bending>
::andesvms(int elm, int maxstr, doublereal nu, 
           doublereal *X, doublereal *Y, doublereal *Z,
           doublereal *_v, doublereal *_stress,
	   int ctyp, int strainflg, int surface)
{
  // Initialized data 
  doublereal clr = 0;
  doublereal cqr = 1;
  doublereal betab = 1;
  doublereal alpha = 1.5;
  doublereal betam = .32;

  // Local variables 
  int i, j;
  doublereal xlp[3], ylp[3], zlp[3];
  doublereal area, thick;
  doublereal str[6];
  doublereal z, epszz;

  Eigen::Matrix<doublereal,9,3> Lb, Lm;
  Eigen::Matrix<doublereal,3,9> Bb, Bm;
  Eigen::Matrix<doublereal,3,3> eframe, gframe = Eigen::Matrix<doublereal,3,3>::Identity();
  Eigen::Matrix<doublereal,18,1> vd;
  Eigen::Map<Eigen::Matrix<doublereal,18,1> > v(_v);
  Eigen::Map<Eigen::Matrix<doublereal,Eigen::Dynamic,Eigen::Dynamic> > stress(_stress,maxstr,3);
  Eigen::Matrix<doublereal,3,1> sigma, epsilon;
  Eigen::Matrix<doublereal,6,1> Upsilon, Sigma;

  // Some convenient definitions 
  Eigen::VectorBlock< Eigen::Matrix<doublereal,6,1> >
    e = Upsilon.head(3), chi = Upsilon.tail(3);
  Eigen::VectorBlock< Eigen::Matrix<doublereal,6,1> >
    N = Sigma.head(3), M = Sigma.tail(3);


// ==================================================================== 
//                                                                      
//     -----------------                                                
//     V A R I A B L E S                                                
//     -----------------                                                
//                                                                      
//     elm      <input>   Finite Element Number                         
//     maxstr   <input>   Maximum Number of Stresses                    
//     nu       <input>   Poisson's Ratio (for an Isotropic Element)    
//     globalX  <input>   X- Nodal Coordinates                          
//     globalY  <input>   Y- Nodal Coordinates                          
//     globalZ  <input>   Z- Nodal Coordinates                          
//     globalU  <input>   Global Displacements at the Nodal Joints      
//     stress   <output>  Stresses (Von Mises Stress) of the Element    
//     ctyp     <input>   Type of Constitutive Law (0, 1, 2, or 3)      
//                                                                      
// ==================================================================== 
// Author   = Francois M. Hemez                                         
// Date     = June 10th, 1995                                           
// Version  = 2.0                                                       
// Modified = K. H. Pierson                                             
// Date     = April 11, 1997                                            
// Reason   = Added stress calculations for sigmaxx, sigmayy, sigmaxy   
//            and von mises stress at top, median and bottom surfaces   
//            Also added strain calculations for epsilonxx, epsilonyy,  
//            epsilonzz, epsilonxy and an equivalent strain at top,     
//            median and bottom surfaces.                               
// ==================================================================== 

    thick = nmat->GetShellThickness();

//     ---------------------------------- 
//     STEP 1                             
//     COMPUTE THE TRIANGULAR COORDINATES 
//     ---------------------------------- 

// .....GET THE ELEMENT TRIANGULAR COORDINATES 
// .....GET THE ELEMENT LEVEL FRAME

    andescrd(elm, X, Y, Z, eframe.data(), xlp, ylp, zlp, area);

//     --------------------------------------------------- 
//     STEP 2                                              
//     COMPUTE THE INTEGRATED CURVATURE-NODAL DISPLACEMENT 
//     MATRIX FOR PURE BENDING (BASIC STIFFNESS MATRIX)    
//     --------------------------------------------------- 

    Lb = Bending<doublereal>::L(xlp, ylp, clr, cqr);

//     ------------------------------------------------- 
//     STEP 3                                            
//     COMPUTE THE INTEGRATED STRAIN-NODAL DISPLACEMENT  
//     MATRIX FOR PURE MEMBRANE (BASIC STIFFNESS MATRIX) 
//     ------------------------------------------------- 

    Lm = Membrane<doublereal>::L(xlp, ylp, alpha);

//     ------------------------------------------- 
//     STEP 4                                      
//     ROTATE THE NODAL DISPLACEMENTS TO THE LOCAL 
//     FRAME SYSTEM (LOCAL TO THE SHELL ELEMENT)   
//     ------------------------------------------- 

    // TODO would be better to store vd and v and 6x3 matrices instead of 18x1 vector
    // then this would simply be a matrix matrix product
    for(i = 0; i < 18; i += 3)
        vd.segment(i,3) = eframe.transpose()*v.segment(i,3);

//     ----------------------------------------------------- 
//     STEP 5                                                
//     COMPUTE THE ELEMENTAL EXTENSION AND CURVATURE VECTORS 
//     ----------------------------------------------------- 

    Eigen::Matrix<int,18,1> indices;
    indices << 0, 1, 6, 7, 12, 13, 5, 11, 17, // M indices
               2, 3, 4, 8, 9, 10, 14, 15, 16; // B indices
    Eigen::PermutationMatrix<18,18,int> PI(indices);

    vd = PI.transpose()*vd;

// .....COMPUTE THE Z- COORDINATE OF THE SELECTED SURFACE

    if(surface == 1) z = thick/2; // upper surface
    else if(surface == 2) z = 0;  // median surface
    else z = -thick/2;            // lower surface

    // compute stresses and strains at the nodes
    doublereal zeta[3][3] = { { 1.,0.,0. }, { 0.,1.,0. }, { 0.,0.,1. } }; // triangular coordinates of nodes
    for(i = 0; i < 3; ++i) {

#ifdef COMPATABILITY_MODE
// .....ELEMENTAL CURVATURE COMPUTATION

        chi = (1/area)*Lb.transpose()*vd.tail(9);

// .....ELEMENTAL EXTENSION COMPUTATION

        e = (1/area)*Lm.transpose()*vd.head(9);
#else
// .....ELEMENTAL CURVATURE COMPUTATION (including now the higher order contribution)

        Bb = (1/area)*Lb.transpose() + Bending<doublereal>::Bd(xlp, ylp, betab, zeta[i]);
        chi = Bb*vd.tail(9);

// .....ELEMENTAL EXTENSION COMPUTATION (including now the higher order contribution)

        Bm = (1/area)*Lm.transpose() +  Membrane<doublereal>::Bd(xlp, ylp, betam, zeta[i]);
        e = Bm*vd.head(9);
#endif

//     -------------------------------------------------
//     STEP 7
//     COMPUTE THE STRESS OR STRAIN
//     -------------------------------------------------

        if (strainflg == 1) {

// .....COMPUTE THE LOCAL STRAINS ON THE SPECIFIED SURFACE

            epsilon = e + z * chi;

            // TODO nu is not necessarily defined, and should be obtained from the material
            // furthermore this equation is only correct for isotropic plane stress material
            epszz = -nu / (1. - nu) * (epsilon[0] + epsilon[1]);

// .....CALCULATE VON MISES EQUIVALENT STRAIN

            stress(6, i) = equivstr(epsilon[0], epsilon[1], epszz, 0.5*epsilon[2]);

// .....ROTATE LOCAL STRAINS TO GLOBAL AND CONVERT SHEAR STRAINS TO ENGINEERING SHEAR STRAINS

            str[0] = epsilon[0];
            str[1] = epsilon[1];
            str[2] = epszz;
            str[3] = 0.5*epsilon[2];
            str[4] = 0.;
            str[5] = 0.;
            transform(eframe.data(), gframe.data(), str);
            for (j = 0; j < 3; ++j)
                stress(j, i) = str[j];
            for (j = 3; j < 6; ++j) 
                stress(j, i) = 2*str[j];

        }
        else {

#ifdef COMPATABILITY_MODE
// .....COMPUTE THE GENERALIZED STRESSES [Sigma = {N,M}] WHICH ARE
// .....FORCE AND MOMENT PER UNIT LENGTH

            if(i == 0 || ctyp == 4)
                nmat->GetConstitutiveResponse(Upsilon.data(), Sigma.data(), NULL, eframe.data(), i);

            if (surface == 1) {

// .....ESTIMATE THE LOCAL STRESSES ON THE UPPER SURFACE

                sigma = N/thick + 6*M/(thick*thick); 

            }

            else if (surface == 2) {

// .....ESTIMATE THE LOCAL STRESSES ON THE MEDIAN SURFACE

                sigma = N/thick;

            }

            else if (surface == 3) {

// .....ESTIMATE THE LOCAL STRESSES ON THE LOWER SURFACE

                sigma = N/thick - 6*M/(thick*thick);

            }
#else

// .....COMPUTE THE LOCAL STRESSES ON THE SPECIFIED SURFACE

            nmat->GetLocalConstitutiveResponse(Upsilon.data(), sigma.data(), z, eframe.data(), i);

#endif

// .....CALCULATE VON MISES EQUIVALENT STRESS

            stress(6, i) = equivstr(sigma[0], sigma[1], 0, sigma[2]);

// .....ROTATE LOCAL STRESSES TO GLOBAL

            str[0] = sigma[0];
            str[1] = sigma[1];
            str[2] = 0.;
            str[3] = sigma[2];
            str[4] = 0.;
            str[5] = 0.;
            transform(eframe.data(), gframe.data(), str);
            for (j = 0; j < 6; ++j)
                stress(j, i) = str[j];

        }
    }

}


template<typename doublereal, template<typename> class Membrane, template<typename> class Bending>
doublereal
ShellElementTemplate<doublereal,Membrane,Bending>
::equivstr(doublereal sxx, doublereal syy, doublereal szz, doublereal sxy)
{
    // Builtin functions 
    using std::sqrt;

    // Local variables 
    doublereal s0, dsxx, dsyy, dszz, eq;

// ... COMPUTE MEAN HYDROSTATIC STRESS OR STRAIN 

    s0 = (sxx + syy + szz) / 3.;

// ... COMPUTE DEVIATORIC STRESSES OR STRAINS 

    dsxx = sxx - s0;
    dsyy = syy - s0;
    dszz = szz - s0;

// ... COMPUTE EQUIVALENT STRESS OR STRAIN 

    eq = (dsxx * dsxx + dsyy * dsyy + dszz * dszz) / 2. + sxy * sxy;
    eq = sqrt(eq * 3.);

    return eq;
}

template<typename doublereal, template<typename> class Membrane, template<typename> class Bending>
void
ShellElementTemplate<doublereal,Membrane,Bending>
::transform(doublereal *_lframe, doublereal *_gframe, doublereal *_str)
{
    // Local variables
    doublereal l1, l2, l3, m1, m2, m3, n1, n2, n3;
    Eigen::Map<Eigen::Matrix<doublereal,3,3> > lframe(_lframe), gframe(_gframe);
    Eigen::Map<Eigen::Matrix<doublereal,6,1> > str(_str);
    Eigen::Matrix<doublereal,6,6> T;


// **********************************************************************

// Purpose: to form a transformation matrix from local coordinates to 
//          global coordinates 

// input variables: 
//      xl = x local unit vector = first column of lframe
//      yl = y local unit vector = second column of lframe
//      zl = z local unit vector = third column of lframe
//      xg = x global unit vector = first column of gframe
//      yg = y global unit vector = second column of gframe
//      zg = z global unit vector = third column of gframe
//      str = stress/strain 6x1 vector 
//            sigmaxx, sigmayy, sigmazz, sigma12, sigma23, sigma13 

// local variables: 
//      l1 = direction cosine between xl and xg 
//      l2 = direction cosine between xl and yg 
//      l3 = direction cosine between xl and zg 
//      m1 = direction cosine between yl and xg 
//      m2 = direction cosine between yl and yg 
//      m3 = direction cosine between yl and zg 
//      n1 = direction cosine between zl and xg 
//      n2 = direction cosine between zl and yg 
//      n3 = direction cosine between zl and zg 
//      T  = transformation matrix from local to global 

// **********************************************************************

// Compute direction cosines 
/*
    l1 = xg[0] * xl[0] + xg[1] * xl[1] + xg[2] * xl[2];
    l2 = yg[0] * xl[0] + yg[1] * xl[1] + yg[2] * xl[2];
    l3 = zg[0] * xl[0] + zg[1] * xl[1] + zg[2] * xl[2];
    m1 = xg[0] * yl[0] + xg[1] * yl[1] + xg[2] * yl[2];
    m2 = yg[0] * yl[0] + yg[1] * yl[1] + yg[2] * yl[2];
    m3 = zg[0] * yl[0] + zg[1] * yl[1] + zg[2] * yl[2];
    n1 = xg[0] * zl[0] + xg[1] * zl[1] + xg[2] * zl[2];
    n2 = yg[0] * zl[0] + yg[1] * zl[1] + yg[2] * zl[2];
    n3 = zg[0] * zl[0] + zg[1] * zl[1] + zg[2] * zl[2];
*/
    Eigen::Matrix<doublereal,3,3> dc = gframe.transpose()*lframe;

    l1 = dc(0,0);
    l2 = dc(1,0);
    l3 = dc(2,0);
    m1 = dc(0,1);
    m2 = dc(1,1);
    m3 = dc(2,1);
    n1 = dc(0,2);
    n2 = dc(1,2);
    n3 = dc(2,2);

// Construct the 6x6 transformation matrix 

    T(0, 0) = l1 * l1;
    T(0, 1) = m1 * m1;
    T(0, 2) = n1 * n1;
    T(0, 3) = l1 * 2 * m1;
    T(0, 4) = m1 * 2 * n1;
    T(0, 5) = n1 * 2 * l1;

    T(1, 0) = l2 * l2;
    T(1, 1) = m2 * m2;
    T(1, 2) = n2 * n2;
    T(1, 3) = l2 * 2 * m2;
    T(1, 4) = m2 * 2 * n2;
    T(1, 5) = n2 * 2 * l2;

    T(2, 0) = l3 * l3;
    T(2, 1) = m3 * m3;
    T(2, 2) = n3 * n3;
    T(2, 3) = l3 * 2 * m3;
    T(2, 4) = m3 * 2 * n3;
    T(2, 5) = n3 * 2 * l3;

    T(3, 0) = l1 * l2;
    T(3, 1) = m1 * m2;
    T(3, 2) = n1 * n2;
    T(3, 3) = l1 * m2 + l2 * m1;
    T(3, 4) = m1 * n2 + m2 * n1;
    T(3, 5) = n1 * l2 + n2 * l1;

    T(4, 0) = l2 * l3;
    T(4, 1) = m2 * m3;
    T(4, 2) = n2 * n3;
    T(4, 3) = l2 * m3 + l3 * m2;
    T(4, 4) = m2 * n3 + m3 * n2;
    T(4, 5) = n2 * l3 + n3 * l2;

    T(5, 0) = l3 * l1;
    T(5, 1) = m3 * m1;
    T(5, 2) = n3 * n1;
    T(5, 3) = l3 * m1 + l1 * m3;
    T(5, 4) = m3 * n1 + m1 * n3;
    T(5, 5) = n3 * l1 + n1 * l3;

// Perform the multiplication {str'} = T{str} 

    str = T*str;
}

template<typename doublereal, template<typename> class Membrane, template<typename> class Bending>
void
ShellElementTemplate<doublereal,Membrane,Bending>
::andesups(int elm, doublereal *state, doublereal *X, doublereal *Y, doublereal *Z, doublereal *_v)
{
  // Initialized data 
  doublereal clr = 0;
  doublereal cqr = 1;
  doublereal betab = 1;
  doublereal alpha = 1.5;
  doublereal betam = .32;

  // Local variables 
  int i;
  doublereal xlp[3], ylp[3], zlp[3];
  doublereal area;

  Eigen::Matrix<doublereal,9,3> Lb, Lm;
  Eigen::Matrix<doublereal,3,9> Bb, Bm;
  Eigen::Matrix<doublereal,3,3> eframe;
  Eigen::Matrix<doublereal,18,1> vd;
  Eigen::Map<Eigen::Matrix<doublereal,18,1> > v(_v);
  Eigen::Matrix<doublereal,6,1> Upsilon;

  // Some convenient definitions 
  Eigen::VectorBlock< Eigen::Matrix<doublereal,6,1> >
    e = Upsilon.head(3), chi = Upsilon.tail(3);

// ==================================================================== 
//                                                                      
//     -----------------                                                
//     V A R I A B L E S                                                
//     -----------------                                                
//                                                                      
//     elm      <input>   Finite Element Number                         
//     maxstr   <input>   Maximum Number of Stresses                    
//     nu       <input>   Poisson's Ratio (for an Isotropic Element)    
//     globalX  <input>   X- Nodal Coordinates                          
//     globalY  <input>   Y- Nodal Coordinates                          
//     globalZ  <input>   Z- Nodal Coordinates                          
//     globalU  <input>   Global Displacements at the Nodal Joints      
//     stress   <output>  Stresses (Von Mises Stress) of the Element    
//     ctyp     <input>   Type of Constitutive Law (0, 1, 2, or 3)      
//                                                                      
// ==================================================================== 
// Author   = Francois M. Hemez                                         
// Date     = June 10th, 1995                                           
// Version  = 2.0                                                       
// Modified = K. H. Pierson                                             
// Date     = April 11, 1997                                            
// Reason   = Added stress calculations for sigmaxx, sigmayy, sigmaxy   
//            and von mises stress at top, median and bottom surfaces   
//            Also added strain calculations for epsilonxx, epsilonyy,  
//            epsilonzz, epsilonxy and an equivalent strain at top,     
//            median and bottom surfaces.                               
// ==================================================================== 

//     ---------------------------------- 
//     STEP 1                             
//     COMPUTE THE TRIANGULAR COORDINATES 
//     ---------------------------------- 

// .....GET THE ELEMENT TRIANGULAR COORDINATES 
// .....GET THE ELEMENT LEVEL FRAME

    andescrd(elm, X, Y, Z, eframe.data(), xlp, ylp, zlp, area);

//     --------------------------------------------------- 
//     STEP 2                                              
//     COMPUTE THE INTEGRATED CURVATURE-NODAL DISPLACEMENT 
//     MATRIX FOR PURE BENDING (BASIC STIFFNESS MATRIX)    
//     --------------------------------------------------- 

    Lb = Bending<doublereal>::L(xlp, ylp, clr, cqr);

//     ------------------------------------------------- 
//     STEP 3                                            
//     COMPUTE THE INTEGRATED STRAIN-NODAL DISPLACEMENT  
//     MATRIX FOR PURE MEMBRANE (BASIC STIFFNESS MATRIX) 
//     ------------------------------------------------- 

    Lm = Membrane<doublereal>::L(xlp, ylp, alpha);

//     ----------------------------------------------------- 
//     STEP 5                                                
//     COMPUTE THE ELEMENTAL EXTENSION AND CURVATURE VECTORS 
//     ----------------------------------------------------- 

    Eigen::Matrix<int,18,1> indices;
    indices << 0, 1, 6, 7, 12, 13, 5, 11, 17, // M indices
               2, 3, 4, 8, 9, 10, 14, 15, 16; // B indices
    Eigen::PermutationMatrix<18,18,int> PI(indices);

    vd = PI.transpose()*v; // note: v is already local in this case

    // compute updated material state at the gauss points
    doublereal zeta[3][3] = { { 0.,.5,.5 }, { .5,0.,.5 }, { .5,.5,0. } }; // triangular coordinates of gauss integration points
    for(i = 0; i < 3; ++i) {

// .....ELEMENTAL CURVATURE COMPUTATION

        Bb = (1/area)*Lb.transpose() + Bending<doublereal>::Bd(xlp, ylp, betab, zeta[i]);
        chi = Bb*vd.tail(9);

// .....ELEMENTAL EXTENSION COMPUTATION

        Bm = (1/area)*Lm.transpose() +  Membrane<doublereal>::Bd(xlp, ylp, betam, zeta[i]);
        e = Bm*vd.head(9);

// .....COMPUTE THE UPDATED MATERIAL STATE AT GAUSS POINT i

        gpmat->UpdateState(Upsilon.data(), state, i);
        state += 5*7; // 5 is the number of layers per gauss point, 7 is the number of states per material point

    }

    // compute updated material state at the nodes
    doublereal eta[3][3] = { { 1.,0.,0. }, { 0.,1.,0. }, { 0.,0.,1. } }; // triangular coordinates of nodes
    for(i = 0; i < 3; ++i) {

#ifdef COMPATABILITY_MODE
// .....ELEMENTAL CURVATURE COMPUTATION

        chi = (1/area)*Lb.transpose()*vd.tail(9);

// .....ELEMENTAL EXTENSION COMPUTATION

        e = (1/area)*Lm.transpose()*vd.head(9);
#else
// .....ELEMENTAL CURVATURE COMPUTATION (including now the higher order contribution)

        Bb = (1/area)*Lb.transpose() + Bending<doublereal>::Bd(xlp, ylp, betab, eta[i]);
        chi = Bb*vd.tail(9);

// .....ELEMENTAL EXTENSION COMPUTATION (including now the higher order contribution)

        Bm = (1/area)*Lm.transpose() +  Membrane<doublereal>::Bd(xlp, ylp, betam, eta[i]);
        e = Bm*vd.head(9);
#endif

// .....COMPUTE THE UPDATED MATERIAL STATE AT NODE i

        nmat->UpdateState(Upsilon.data(), state, i);
        state += 3*7; // 3 is the number of layers per node, 7 is the number of states per material point

    }

}
#endif