#ifndef contact_tolerances_h_
#define contact_tolerances_h_

//---------------------------------------------------------------
//
//        ACME's G L O B A L   T O L E R A N C E S
//
//---------------------------------------------------------------
// see also: search_parameters.par & Contact_Defines.h
// (as well as the individual source files)

// These are the MAX & MIN tolerances for time to contact
// They are used in Interaction_Definition() to determine
// if a moving interaction is valid.
//static Real time_to_contact_tolerance_min = -0.55;
//static Real time_to_contact_tolerance_max =  1.05;
#define time_to_contact_tolerance_min  -0.55
#define time_to_contact_tolerance_max   1.05

// inflation tolerance on bounding boxes: related to JAS's tol_msl & tol_snm
//static Real BOX_INFLATION_FACTOR = 0.01;
#define BOX_INFLATION_FACTOR  0.01

// tolerance that inflates the max remaining gap to get a "reasonable gap" 
//static Real GAP_INFLATION_FACTOR = 0.01;
#define GAP_INFLATION_FACTOR  0.01

// The computed areas in cnodetriange_cpproj for a point exactly on a
// corner (or possibly on an edge) may yield a very small negative number.
// REL_TANG_TOL is used in an if-test to avoid bad consequences from this
// possibility.  A better approach to handle this is desirable. (RMS)
// the tolerance on subtriangle to triangle area in cnodetriangle_cpproj, etc
// NOTE: this was "CTOLT1"
//static Real REL_TANG_TOL  = 1.e-12;
#define REL_TANG_TOL   1.e-12

#endif
