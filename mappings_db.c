#include "mappings_db.h"

int create_sorting_database(sqlite3** db, const char* db_folder) {
    char db_name[strlen(db_folder) + 14];
    sprintf(db_name, "%s/mappings.db", db_folder);
    
    char* error_msg;
    int ret_code = sqlite3_open(db_name, db);
    // If database could be open, create its tables
    if (ret_code == SQLITE_OK) {
        ret_code = sqlite3_exec(*db, "BEGIN TRANSACTION", NULL, NULL, &error_msg);
        
        ret_code += sqlite3_exec(*db, 
                                 "CREATE TABLE mappings (\
                                 seq_id VARCHAR(32), \
                                 paired_end INTEGER)", 
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
            num_mappings = sqlite3_column_int(count_query_stmt, 1);
            num_mappings_histogram[num_mappings]++;
        } else if (sqlite_code == SQLITE_DONE) {  //no more rows, end reached
            break;
        } else {
            LOG_FATAL("Different code from SQLITE_ROW and SQLITE_DONE received while iterating results of statement. Aborting program\n");
        }
    }

    sqlite3_finalize(count_query_stmt);
    
    return num_mappings_histogram;
}

int insert_mapping(sqlite3* db, char* seq_id, short int paired_end, sqlite3_stmt *insert_stmt) {
    // insert mapping
    sqlite3_bind_text(insert_stmt, 1, seq_id, strlen(seq_id), SQLITE_TRANSIENT);
    sqlite3_bind_int(insert_stmt, 2, paired_end);
    
    if (sqlite3_step(insert_stmt) != SQLITE_DONE) {
        char log_message[200];
        sprintf("Could not insert mapping with seq_id: %s (paired-end %i)\n", seq_id, paired_end);
        LOG_FATAL(log_message);
    }
    
    sqlite3_reset(insert_stmt);
    
    return 0;
}