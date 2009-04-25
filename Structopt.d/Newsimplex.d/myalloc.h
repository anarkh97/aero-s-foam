/*********************************************************************/
/***    Copyright (c) Robert J. Vanderbei, 1994                    ***/
/***    All Rights Reserved                                        ***/
/*********************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <string.h>
/*
#include <dmalloc.h>
*/
#include <malloc.h>

#undef MALLOC
#define	MALLOC(name,len,type) {	\
	if ( (name) != NULL )  printf( "MALLOC:" #name " not initialized\n"); \
        if ( (len) > 0 ) { \
	  if ( ( (name) = (type *)malloc( (len) * sizeof(type) ) ) == NULL \
	   	  && (len)>0) { \
	      printf( "MALLOC(" #name "," #len "=%d," #type ")" \
	  	    ": cannot allocate space", len); \
	      exit(1); \
        }} else { \
         (name) = NULL; \
	} \
}

#undef CALLOC
#define	CALLOC(name,len,type) { \
       if ( (name) != NULL )  printf( "CALLOC:" #name " not initialized\n"); \
       if ( (len) > 0 ) { \
	if ( ( (name) = (type *)calloc( (len) , sizeof(type) ) ) == NULL \
	        && (len)>0) { \
	    printf( "CALLOC(" #name "," #len "=%d," #type ")" \
		  ": cannot allocate space", len); \
	    exit(1); \
        }} else { \
         (name) = NULL; \
	} \
}

#undef REALLOC
#define	REALLOC(name,len,type) { \
        if ( (len) > 0 ) { \
	if (((name) = (type *)realloc((name),(len)*sizeof(type))) == NULL \
	       && (len)>0) { \
	    printf( "REALLOC(" #name "," #len "=%d," #type ")" \
		  ": cannot reallocate space", len); \
	    exit(1); \
        } } else { \
         FREE(name) \
        }\
}

#undef FREE
#define	FREE(name) { \
	if ( (name) != NULL ) free( (name) ); \
	(name) = NULL; \
}
