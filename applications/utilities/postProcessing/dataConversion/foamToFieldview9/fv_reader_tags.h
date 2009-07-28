#ifndef FV_READER_TAGS_H
#define FV_READER_TAGS_H

/* Numeric tags (codes) for FIELDVIEW binary file format. */

#define FV_MAGIC	0x00010203	/* decimal 66051 */

/* Content of the file (grid only, results only or combined). */
#define FV_GRIDS_FILE           1
#define FV_RESULTS_FILE         2
#define FV_COMBINED_FILE        3

#define FV_NODES        	1001
#define FV_FACES        	1002
#define FV_ELEMENTS     	1003
#define FV_VARIABLES    	1004
#define FV_BNDRY_VARS   	1006
#define FV_ARB_POLY_FACES       1007
#define FV_ARB_POLY_ELEMENTS    1008
#define FV_ARB_POLY_BNDRY_VARS  1009

#define FV_TET_ELEM_ID          1
#define FV_HEX_ELEM_ID          2
#define FV_PRISM_ELEM_ID        3
#define FV_PYRA_ELEM_ID         4
#define FV_ARB_POLY_ELEM_ID     5

/* Values for "wall_info" array (see comments in fv_encode_elem_header). */
#ifdef __STDC__
#define A_WALL         (07u)
#define NOT_A_WALL     (0u)
#else
#define A_WALL         (07)
#define NOT_A_WALL     (0)
#endif

#endif /* FV_READER_TAGS_H */
