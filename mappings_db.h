
#ifndef MAPPINGS_DB_H
#define MAPPINGS_DB_H

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "sqlite3.h"

#include "bam_commons.h"
#include "commons.h"
#include "log.h"

/**
*  @brief Creates a database for storing mapping info
*  @param db pointer to the sqlite database
*  @param db_folder output directory where .db file will be place
*  @param memory flag indicating if table will be stored in memory or not (0: disk, 1: memory)
*  @return int error code (0: ok, otherwise: error)
*  
*  Creates a database for storing mapping info (only id_seq and paired_end)
*/
int mapping_db_create_short_mappings_database(sqlite3** db, const char* db_folder, int memory);

/**
*  @brief Creates a database for storing mapping info
*  @param db pointer to the sqlite database
*  @param db_folder output directory where .db file will be place
*  @param memory flag indicating if table will be stored in memory or not (0: disk, 1: memory)
*  @return int error code (0: ok, otherwise: error)
*  
*  Creates a database for storing mapping info with whole fields (id_seq, paired_end, tid, mtid, isize, start_coordinate, seq_length)
*/
int mapping_db_create_complete_mappings_database(sqlite3** db, const char* db_folder, int memory);

/**
*  @brief Creates indexes in the mapping database
*  @param db pointer to the sqlite database
*  @return int error code
*  
*  Creates indexes in the mapping database for id_seq field  
*/
int mapping_db_create_indexes(sqlite3** db);

/**
*  @brief Performs calculations of mapping histogram and mean paired end distance
*  @param db pointer to the sqlite database
*  @param max_insert_size max size between paired ends to consider mappings in calculations
*  @param[in,out] num_mappings_histogram pointer to the mapping histogram
*  @param[in,out] mean_paired_end_distance pointer to the mean paired end distance value
*  @return void
*  
*  Performs calculations of mapping histogram and mean paired end distance over the database data
*/
void mapping_db_perform_calculations(sqlite3* db, int max_insert_size, unsigned int* num_mappings_histogram, unsigned long* mean_paired_end_distance);

/**
*  @brief Creates and prepare an insert statement
*  @param db pointer to the sqlite database
*  @return sqlite3_stmt pointer to the prepared statements
*  
*  Creates and prepare an insert statement
*/
sqlite3_stmt* mapping_db_prepare_insert_complete(sqlite3* db);

/**
*  @brief Performs calculations of mapping histogram and mean paired end distance
*  @param db pointer to the sqlite database
*  @param id_seq sequence id
*  @param pairend_end number of paired end (1 or 2)
*  @param insert_stmt pointer to the insert statement
*  @return void
*  
*  Performs calculations of mapping histogram and mean paired end distance over the database data
*/
int mapping_db_insert_short(sqlite3* db, char* id_seq, short int paired_end, sqlite3_stmt* insert_stmt);

/**
*  @brief Performs calculations of mapping histogram and mean paired end distance
*  @param db pointer to the sqlite database
*  @param id_seq sequence id 
*  @param pairend_end number of paired end (1 or 2)
*  @param tid chromosome
*  @param mtid mate chromosome
*  @param start_coordinate start coordinate of the alignment
*  @param seq_length sequence id length
*  @param insert_stmt pointer to the insert statement
*  @return void
*  
*  Performs calculations of mapping histogram and mean paired end distance over the database data
*/
int mapping_db_insert_complete(sqlite3* db, char* id_seq, short int paired_end, int tid, int mtid, int isize, int start_coordinate, int seq_length, sqlite3_stmt* insert_stmt);

/**
*  @brief Begins a transaction in the sqlite database
*  @param db pointer to the sqlite database
*  @return void
*  
*  Begins a transaction in the sqlite database
*/
void mapping_db_begin_transaction(sqlite3* db);

/**
*  @brief Ends a transaction in the sqlite database
*  @param db pointer to the sqlite database
*  @return void
*  
*  Ends a transaction in the sqlite database
*/
void mapping_db_end_transaction(sqlite3* db);

/**
*  @brief Closes the given sqlite database
*  @param db pointer to the sqlite database
*  @return int error code
*  
*  Closes the given sqlite database
*/
int mapping_db_close(sqlite3* db);

#endif
