#ifndef MAPPINGS_DB_H
#define MAPPINGS_DB_H

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "sqlite3.h"

#include "log.h"

int create_mappings_database(sqlite3** db, const char* db_folder);

int create_indexes(sqlite3** db);

int* get_num_mappings_histogram(sqlite3* db, int* num_mappings_histogram);

int insert_mapping(sqlite3* db, char* seq_id, short int paired_end, sqlite3_stmt *insert_stmt);

#endif
