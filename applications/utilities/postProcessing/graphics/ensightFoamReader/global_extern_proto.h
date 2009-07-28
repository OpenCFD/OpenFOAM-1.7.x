/*--------------------------------------------------------------*/
/* Prototype Header file for EnSight External Reader            */
/* DSO Library Routines                                         */
/*                                                              */
/* intended to be included from global_extern.h only            */
/*--------------------------------------------------------------*/
/*  *************************************************************
 *   Copyright 1998 Computational Engineering International, Inc.
 *   All Rights Reserved.
 *
 *        Restricted Rights Legend
 *
 *   Use, duplication, or disclosure of this
 *   software and its documentation by the
 *   Government is subject to restrictions as
 *   set forth in subdivision [(b)(3)(ii)] of
 *   the Rights in Technical Data and Computer
 *   Software clause at 52.227-7013.
 *  *************************************************************
 */
#ifndef GLOBAL_EXTERN_PROTO_H
#define GLOBAL_EXTERN_PROTO_H

#ifdef WIN32
#define W32IMPORT __declspec( dllimport )
#define W32EXPORT __declspec( dllexport )
#else
#define W32IMPORT extern
#define W32EXPORT extern
#endif

/*----------------------
 * Same in All Versions
 *----------------------*/
W32IMPORT int
USERD_get_number_of_model_parts( void );

W32IMPORT int
USERD_get_block_coords_by_component(int block_number,
                                    int which_component,
                                    float *coord_array);

W32IMPORT int
USERD_get_block_iblanking(int block_number,
                          int *iblank_array);

W32IMPORT int
USERD_get_block_scalar_values(int block_number,
                              int which_scalar,
                              float *scalar_array);

W32IMPORT int
USERD_get_block_vector_values_by_component(int block_number,
                                           int which_vector,
                                           int which_component,
                                           float *vector_array);

W32IMPORT int
USERD_get_name_of_reader(char reader_name[Z_MAX_USERD_NAME],
                         int *two_fields);
     
W32IMPORT int
USERD_get_reader_descrip(char descrip[Z_MAXFILENP]);

W32IMPORT int
USERD_set_filenames(char filename_1[],
                    char filename_2[],
                    char the_path[],
                    int swapbytes);
     
W32IMPORT int
USERD_get_number_of_files_in_dataset( void );

W32IMPORT int
USERD_get_dataset_query_file_info(Z_QFILES *qfiles);
     
W32IMPORT int
USERD_get_changing_geometry_status( void );

W32IMPORT int
USERD_get_node_label_status( void );

W32IMPORT int
USERD_get_element_label_status( void );

W32IMPORT int
USERD_get_number_of_variables( void );

W32IMPORT void
USERD_stop_part_building( void );

W32IMPORT int
USERD_bkup(FILE *archive_file,
           int backup_type);



/*-----------------------
 * For Version 1.000 Only
 *-----------------------*/
#if defined USERD_API_100

W32IMPORT int
USERD_get_number_of_global_nodes( void );
     
W32IMPORT int
USERD_get_global_coords(CRD *coord_array);
     
W32IMPORT int
USERD_get_global_node_ids(int *nodeid_array);
     
W32IMPORT int
USERD_get_element_connectivities_for_part(int part_number,
                                          int **conn_array[Z_MAXTYPE]);

W32IMPORT int
USERD_get_element_ids_for_part(int part_number,
                               int *elemid_array[Z_MAXTYPE]);

W32IMPORT int
USERD_get_vector_values(int which_vector,
                        int which_part,
                        int which_type,
                        float *vector_array);
     
W32IMPORT int
USERD_get_part_build_info(int *part_id,
                          int *part_types,
                          char *part_descriptions[Z_BUFL],
                          int *number_of_elements[Z_MAXTYPE],
                          int *ijk_dimensions[3],
                          int *iblanking_options[6]);

W32IMPORT int
USERD_get_scalar_values(int which_scalar,
                        int which_part,
                        int which_type,
                        float *scalar_array);

W32IMPORT int
USERD_get_variable_info(char **var_description,
                        char **var_filename,
                        int *var_type,
                        int *var_classify);

W32IMPORT int
USERD_get_description_lines(int which_type,
                            int which_var,
                            char line1[Z_BUFL],
                            char line2[Z_BUFL]);

W32IMPORT int
USERD_get_variable_value_at_specific(int which_var,
                                     int which_node_or_elem,
                                     int which_part,
                                     int which_elem_type,
                                     int time_step,
                                     float values[3]);

W32IMPORT float
USERD_get_constant_value(int which_var);

W32IMPORT int
USERD_get_solution_times(float *solution_times);
W32IMPORT void
USERD_set_time_step(int time_step);

W32IMPORT int
USERD_get_number_of_time_steps(void);

#endif


/*----------------------
 * New For Version 2.000
 *----------------------*/
#if !defined USERD_API_100

W32IMPORT int
USERD_get_part_coords(int part_number,
                      float **coord_array);

W32IMPORT int
USERD_get_part_node_ids(int part_number,
                        int *nodeid_array);
     
W32IMPORT int
USERD_get_part_elements_by_type(int part_number,
                                int element_type,
                                int **conn_array);
W32IMPORT int
USERD_get_part_element_ids_by_type(int part_number,
                                   int element_type,
                                   int *elemid_array);
     
W32IMPORT int
USERD_get_reader_version(char version_number[Z_MAX_USERD_NAME]);

W32IMPORT int
USERD_get_reader_release(char version_number[Z_MAX_USERD_NAME]);

W32IMPORT int
USERD_get_var_by_component(int which_variable,
                           int which_part,
                           int var_type,
                           int which_type,
                           int complex,
                           int component,
                           float *var_array);

W32IMPORT int
USERD_get_maxsize_info(int *max_number_of_nodes,
                       int *max_number_of_elements[Z_MAXTYPE],
                       int *max_ijk_dimensions[3]);

W32IMPORT void
USERD_exit_routine( void );

W32IMPORT int
USERD_get_gold_variable_info(char **var_description,
                             char **var_filename,
                             int *var_type,
                             int *var_classify,
                             int *var_complex,
                             char **var_ifilename,
                             float *var_freq,
                             int *var_contran,
                             int *var_timeset);
W32IMPORT int
USERD_get_model_extents( float extents[6] );

W32IMPORT int
USERD_get_descrip_lines(int which_type,
                        int which_var,
                        int imag_data,
                        char line1[Z_BUFL],
                        char line2[Z_BUFL]);

W32IMPORT int
USERD_get_var_value_at_specific(int which_var,
                                int which_node_or_elem,
                                int which_part,
                                int which_elem_type,
                                int time_step,
                                float values[3],
                                int imag_data);

W32IMPORT float
USERD_get_constant_val(int which_var, int imag_data);

W32IMPORT int
USERD_get_geom_timeset_number(void);

W32IMPORT int
USERD_get_number_of_timesets(void);

W32IMPORT int
USERD_get_timeset_description(int timeset_number,
                              char timeset_description[Z_BUFL]);

W32IMPORT int
USERD_get_sol_times(int timeset_number,
                    float *solution_times);
W32IMPORT void
USERD_set_time_set_and_step(int timeset_number,
                            int time_step);
W32IMPORT int
USERD_get_num_of_time_steps(int timeset_number);

W32IMPORT int
USERD_get_border_availability(int part_number,
                              int number_of_elements[Z_MAXTYPE]);

W32IMPORT int
USERD_get_border_elements_by_type(int part_number,
                                  int element_type,
                                  int **conn_array,
                                  short *parent_element_type,
                                  int *parent_element_num);

W32IMPORT void
USERD_set_server_number(int serv_num,
                        int tot_servs);

#endif


/*----------------------
 * New For Version 2.010
 *----------------------*/
#if defined USERD_API_201 || defined USERD_API_202 || defined USERD_API_203
W32IMPORT int
USERD_get_ghosts_in_model_flag( void );

W32IMPORT int
USERD_get_ghosts_in_block_flag(int block_number);

W32IMPORT int
USERD_get_block_ghost_flags(int block_number,
                            int *ghost_flags);
#endif

/*--------------------------
 * Modified at Version 2.030
 *--------------------------*/
#if defined USERD_API_201 || defined USERD_API_202

W32IMPORT int
USERD_get_gold_part_build_info(int *part_id,
                               int *part_types,
                               char *part_descriptions[Z_BUFL],
                               int *number_of_nodes,
                               int *number_of_elements[Z_MAXTYPE],
                               int *ijk_dimensions[3],
                               int *iblanking_options[6]);
#endif

#if defined USERD_API_203
W32IMPORT int
USERD_get_gold_part_build_info(int *part_id,
                               int *part_types,
                               char *part_descriptions[Z_BUFL],
                               int *number_of_nodes,
                               int *number_of_elements[Z_MAXTYPE],
                               int *ijk_dimensions[9],
                               int *iblanking_options[6]);
#endif


/*----------------------
 * New For Version 2.030
 *----------------------*/
#if defined USERD_API_203
W32IMPORT int
USERD_get_number_of_material_sets( void );

W32IMPORT int
USERD_get_matf_set_info(int *mat_set_ids,
                        char **mat_set_name);

W32IMPORT int
USERD_get_number_of_materials( int set_index );

W32IMPORT int
USERD_get_matf_var_info(int set_index,
                        int *mat_ids,
                        char **mat_desc);

W32IMPORT int
USERD_size_matf_data(int set_index,
                     int part_id,
                     int wtyp,
                     int mat_type,
                     int *matf_size );

W32IMPORT int
USERD_load_matf_data( int set_index,
                      int part_id,
                      int wtyp,
                      int mat_type,
                      int *ids_list,
                      float *val_list );

W32IMPORT int
USERD_get_nsided_conn( int part_number,
                       int *nsided_conn_array );

W32IMPORT int
USERD_get_nfaced_nodes_per_face( int part_number,
                                 int *nfaced_npf_array );

W32IMPORT int
USERD_get_nfaced_conn( int part_number,
                       int *nfaced_conn_array );


#endif

     
/*--------------------------------------------------------------------*/
#endif /*GLOBAL_EXTERN_PROTO_H*/

