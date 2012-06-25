#include "mappings_db.h"

int create_mappings_database(sqlite3** db, const char* db_folder, int memory) {
    char db_name[strlen(db_folder) + 14];
    sprintf(db_name, "%s/mappings.db", db_folder);
    
    char* error_msg;
    int ret_code;
    
    if (memory) {
        ret_code = sqlite3_open(":memory:", db);
    } else {
        ret_code = sqlite3_open(db_name, db);
    }
    
    // If database could be open, create its tables
    if (ret_code == SQLITE_OK) {
        ret_code = sqlite3_exec(*db, "BEGIN TRANSACTION", NULL, NULL, &error_msg);
        
        ret_code += sqlite3_exec(*db, 
                                 "CREATE TABLE mappings (\
                                 id_seq VARCHAR(32), \
                                 paired_end INTEGER)", 
                                 NULL, 0, &error_msg);
        
        sqlite3_exec(*db, "END TRANSACTION", NULL, NULL, &error_msg);
    }
    
    // Return 0 on success
    return ret_code - SQLITE_OK;
}

int create_complete_mappings_database(sqlite3** db, const char* db_folder, int memory) {
    char db_name[strlen(db_folder) + 14];
    sprintf(db_name, "%s/mappings.db", db_folder);
    
    char* error_msg;
    int ret_code;
    
    if (memory) {
        ret_code = sqlite3_open(":memory:", db);
    } else {
        ret_code = sqlite3_open(db_name, db);
    }
    
    // If database could be open, create its tables
    if (ret_code == SQLITE_OK) {
        ret_code = sqlite3_exec(*db, "BEGIN TRANSACTION", NULL, NULL, &error_msg);
        
        ret_code += sqlite3_exec(*db, 
                                 "CREATE TABLE mappings (\
                                 id_seq VARCHAR(32), \
                                 paired_end INTEGER, \
                                 tid INTEGER, \
                                 start_coordinate INTEGER, \
                                 seq_length INTEGER)",
                                 NULL, 0, &error_msg);
        
        sqlite3_exec(*db, "END TRANSACTION", NULL, NULL, &error_msg);
    }
    
    // Return 0 on success
    return ret_code - SQLITE_OK;
}

int create_indexes(sqlite3** db) {
    char *error_msg;
    int ret_code = 0;
        
    ret_code += sqlite3_exec(*db, 
                             "CREATE INDEX seq_id ON mappings (seq_id)", 
                             NULL, 0, &error_msg);
        
    // Return 0 on success
    return ret_code - SQLITE_OK;
}

int* get_num_mappings_histogram(sqlite3* db, int* num_mappings_histogram) {
    int num_mappings, sqlite_code;
    char query_buf[] = "SELECT COUNT(*), id_seq FROM mappings GROUP BY id_seq";
    sqlite3_stmt* *count_query_stmt;
    
    // ...and retrieve the ID assigned to it
    sqlite3_prepare_v2(db, query_buf, strlen(query_buf), &count_query_stmt, NULL);
    
    while (1) {
        sqlite_code = sqlite3_step(count_query_stmt);
        
        if (sqlite_code == SQLITE_ROW) {  //more rows
            num_mappings = sqlite3_column_int(count_query_stmt, 0);

            if (num_mappings > (MAX_MAPPING_COUNT_IN_HISTOGRAM + 1)) {
                num_mappings = MAX_MAPPING_COUNT_IN_HISTOGRAM + 1;  //histogram counts to 10 or more mappings
            }
            num_mappings_histogram[num_mappings - 1]++;
        } else if (sqlite_code == SQLITE_DONE) {  //no more rows, end reached
            break;
        } else {
            LOG_FATAL("Different code from SQLITE_ROW and SQLITE_DONE received while iterating results of statement. Aborting program\n");
        }
    }

    sqlite3_finalize(count_query_stmt);
    
    return num_mappings_histogram;
}

int insert_mapping(sqlite3* db, char* id_seq, short int paired_end, sqlite3_stmt* insert_stmt) {
    // insert mapping
    sqlite3_bind_text(insert_stmt, 1, id_seq, strlen(id_seq), SQLITE_TRANSIENT);
    sqlite3_bind_int(insert_stmt, 2, paired_end);

    if (sqlite3_step(insert_stmt) != SQLITE_DONE) {
        char log_message[200];
        sprintf("Could not insert mapping with seq_id: %s (paired-end %i)\n", id_seq, paired_end);
        LOG_FATAL(log_message);
    }
    
    sqlite3_reset(insert_stmt);
    
    return 0;
}

int insert_complete_mapping(sqlite3* db, char* id_seq, short int paired_end, int tid, int start_coordinate, int seq_length, sqlite3_stmt* insert_stmt) {
    // insert mapping
    sqlite3_bind_text(insert_stmt, 1, id_seq, strlen(id_seq), SQLITE_TRANSIENT);
    sqlite3_bind_int(insert_stmt, 2, paired_end);
    sqlite3_bind_int(insert_stmt, 3, tid);
    sqlite3_bind_int(insert_stmt, 4, start_coordinate);
    sqlite3_bind_int(insert_stmt, 5, seq_length);

    if (sqlite3_step(insert_stmt) != SQLITE_DONE) {
        printf("error description: %s\n", sqlite3_errmsg(db));
        char log_message[200];
        //sprintf("Could not insert mapping with seq_id: %s (paired-end %i)\n", id_seq, paired_end);
        //LOG_FATAL(log_message);
    }
    
    sqlite3_reset(insert_stmt);
    
    return 0;    
}