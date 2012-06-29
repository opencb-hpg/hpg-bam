#include "mappings_db.h"

int mapping_db_create_short_mappings_database(sqlite3** db, const char* db_folder, int memory) {
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
                                 id_seq VARCHAR(64), \
                                 paired_end INTEGER)", 
                                 NULL, 0, &error_msg);

        //create_indexes(db);  //index creation before insertion
        
        sqlite3_exec(*db, "END TRANSACTION", NULL, NULL, &error_msg);
    }
    
    // Return 0 on success
    return ret_code - SQLITE_OK;
}

// int mapping_db_create_complete_mappings_database(sqlite3** db, const char* db_folder, int memory) {
//     char db_name[strlen(db_folder) + 14];
//     sprintf(db_name, "%s/mappings.db", db_folder);
//     
//     char* error_msg;
//     int ret_code;
//     
//     if (memory) {
//         ret_code = sqlite3_open(":memory:", db);
//     } else {
//         ret_code = sqlite3_open(db_name, db);
//     }
//     
//     // If database could be open, create its tables
//     if (ret_code == SQLITE_OK) {
//         ret_code = sqlite3_exec(*db, "BEGIN TRANSACTION", NULL, NULL, &error_msg);
//         
//         ret_code += sqlite3_exec(*db, 
//                                  "CREATE TABLE mappings (\
//                                  id_seq VARCHAR(64), \
//                                  paired_end INTEGER, \
//                                  tid INTEGER, \
//                                  mtid INTEGER, \
//                                  isize INTEGER, \
//                                  start_coordinate INTEGER, \
//                                  seq_length INTEGER)",
//                                  NULL, 0, &error_msg);        
//         
//         sqlite3_exec(*db, "END TRANSACTION", NULL, NULL, &error_msg);
//     }
//     
//     // Return 0 on success
//     return ret_code - SQLITE_OK;
//}

int mapping_db_create_complete_mappings_database(sqlite3** db, const char* db_folder, int memory) {
    char* error_msg;
    int ret_code;
  
    char db_file[strlen(db_folder) + 14];
    sprintf(db_file, "%s/mappings.db", db_folder);
    remove(db_file);
    
    if (memory) {
        ret_code = sqlite3_open(":memory:", db);
    } else {
        ret_code = sqlite3_open(db_file, db);
    }
    
    // If database could be open, create its tables
    if (ret_code == SQLITE_OK) {
        ret_code = sqlite3_exec(*db, "BEGIN TRANSACTION", NULL, NULL, &error_msg);
        
        ret_code += sqlite3_exec(*db, 
                                 "CREATE TABLE mappings (\
                                 id_seq VARCHAR(64), \
                                 paired_end INTEGER, \
                                 tid INTEGER, \
                                 mtid INTEGER, \
                                 isize INTEGER, \
                                 start_coordinate INTEGER, \
                                 seq_length INTEGER)",
                                 NULL, 0, &error_msg);        
        
        sqlite3_exec(*db, "END TRANSACTION", NULL, NULL, &error_msg);
    }
    
    if (ret_code != SQLITE_OK) {
        printf("ret_code = %d, reason: %s\n", ret_code, sqlite3_errmsg(db));
        LOG_FATAL("Can't create temporary database for calculating mapping histogram. Reason:\n");
        sqlite3_close(db);
    }    
    
    // Return 0 on success
    return ret_code - SQLITE_OK;
}

int mapping_db_create_indexes(sqlite3** db) {
    char *error_msg;
    int ret_code = 0;
     
    ret_code += sqlite3_exec(*db, 
                             "CREATE INDEX id_seq ON mappings (id_seq)", 
                             NULL, 0, &error_msg);

    // Return 0 on success
    return ret_code - SQLITE_OK;
}

void mapping_db_perform_calculations(sqlite3* db, int max_insert_size, unsigned int* num_mappings_histogram, unsigned long* mean_paired_end_distance) {
    int num_mappings, sqlite_code, num_paired_ends;
    unsigned long mean_paired_end_distance_acc = 0;
    char* query_buf = (char*) calloc(100, sizeof(char));
    sprintf(query_buf, "%s %i %s", "SELECT COUNT(*), SUM(isize), id_seq FROM mappings WHERE tid = mtid AND isize <", max_insert_size, "GROUP BY id_seq");
    sqlite3_stmt** count_query_stmt;
    
    // ...and retrieve the ID assigned to it
    sqlite3_prepare_v2(db, query_buf, strlen(query_buf), &count_query_stmt, NULL);
    
    while (1) {
        sqlite_code = sqlite3_step(count_query_stmt);
        
        if (sqlite_code == SQLITE_ROW) {  //more rows
            num_mappings = sqlite3_column_int(count_query_stmt, 0);
            mean_paired_end_distance_acc += sqlite3_column_int(count_query_stmt, 1);

            if (num_mappings > (MAX_MAPPING_COUNT_IN_HISTOGRAM + 1)) {
                num_mappings = MAX_MAPPING_COUNT_IN_HISTOGRAM + 1;  //histogram counts to 10 or more mappings
            }
            num_mappings_histogram[num_mappings - 1]++;
            num_paired_ends += (num_mappings - 1);        
        } else if (sqlite_code == SQLITE_DONE) {  //no more rows, end reached
            break;
        } else {
            char log_message[500];
            sprintf(log_message, "%s%s\n", "Different code from SQLITE_ROW and SQLITE_DONE received while iterating results of statement. Aborting program. Error: \n", sqlite3_errmsg(db));
            LOG_FATAL(log_message);
        }
    }

    sqlite3_finalize(count_query_stmt);
    *mean_paired_end_distance = (mean_paired_end_distance_acc / num_paired_ends);
    
    num_mappings_histogram[0] = num_alignments - 2 * num_mappings_histogram[1] - 3 * num_mappings_histogram[2] - 4 * num_mappings_histogram[3] - 
                                                 5 * num_mappings_histogram[4] - 6 * num_mappings_histogram[5] - 7 * num_mappings_histogram[6] - 
                                                 8 * num_mappings_histogram[7] - 9 * num_mappings_histogram[8] - 10 * num_mappings_histogram[9] - 
                                                 11 * num_mappings_histogram[10];
    
    return;
}

sqlite3_stmt* mapping_db_prepare_insert_complete(sqlite3* db) {
    sqlite3_stmt* insert_complete_mapping_stmt;
    char insert_complete_mapping_buffer[] = "INSERT INTO mappings VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7)";
    
    sqlite3_prepare_v2(db, insert_complete_mapping_buffer, -1, &insert_complete_mapping_stmt, NULL);

    return insert_complete_mapping_stmt;
}

int mapping_db_insert_short(sqlite3* db, char* id_seq, short int paired_end, sqlite3_stmt* insert_stmt) {
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

int mapping_db_insert_complete(sqlite3* db, char* id_seq, short int paired_end, int tid, int mtid, int isize, int start_coordinate, int seq_length, sqlite3_stmt* insert_stmt) {
    // insert mapping
    sqlite3_bind_text(insert_stmt, 1, id_seq, strlen(id_seq), SQLITE_TRANSIENT);
    sqlite3_bind_int(insert_stmt, 2, paired_end);
    sqlite3_bind_int(insert_stmt, 3, tid);
    sqlite3_bind_int(insert_stmt, 4, mtid);
    sqlite3_bind_int(insert_stmt, 5, isize);
    sqlite3_bind_int(insert_stmt, 6, start_coordinate);
    sqlite3_bind_int(insert_stmt, 7, seq_length);

    if (sqlite3_step(insert_stmt) != SQLITE_DONE) {
        printf("error description: %s\n", sqlite3_errmsg(db));
        char log_message[200];
        //sprintf("Could not insert mapping with seq_id: %s (paired-end %i)\n", id_seq, paired_end);
        //LOG_FATAL(log_message);
    }
    
    sqlite3_reset(insert_stmt);
    
    return 0;    
}

void mapping_db_begin_transaction(sqlite3* db) {
    char* zErrMsg;
    sqlite3_exec(db, "PRAGMA cache_size = 100000", NULL, NULL, &zErrMsg);
    sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &zErrMsg);  
}

void mapping_db_end_transaction(sqlite3* db) {
    char* zErrMsg;
    sqlite3_exec(db, "END TRANSACTION", NULL, NULL, &zErrMsg);
}

int mapping_db_close(sqlite3* db) {
    int ret_code = sqlite3_close(db);

    if (ret_code != SQLITE_OK) {
        LOG_WARN("Temporary database could not be closed");
        printf("closing reason: %s\n", sqlite3_errmsg(db));
    }
    
}