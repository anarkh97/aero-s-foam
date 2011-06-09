ADD_LIBRARY(Corotational  Corotator.C 
       BeamCorotator.C inc_rottensor.C 
       Shell3Corotator.C BarCorotator.C 
       BarThermalCorotator.C 
       rightmult_rotvar.C mat_to_vec.C 
       pseudorot_2var.C 
       GeomState.C leftmult_rotvar.C 
       pseudorot_var.C vec_to_quat.C 
       form_rottensor.C mat_to_quat.C 
       fitalg3_3nodnew.C vec_to_mat.C 
       quat_to_mat.C crossprod.C  
       tran_fsl.C mat_mult_mat.C  
       GeomNLSolver.C orthonorm3.C 
       DistrGeomState.C normalize.C tran_force.C
       tran_stiff.C SpringCorotator.C 
       TetCorotator.C 
       BrickCorotator.C 
       SuperCorotator.C 
       PentaCorotator.C BarFCorotator.C
       QuadThermalCorotator.C TriangleThermalCorotator.C
       PhantomCorotator.C)


SET(INC_DIR ..)

INCLUDE_DIRECTORIES( ${INC_DIR} )
