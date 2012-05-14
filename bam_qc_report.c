
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "bam_qc_report.h"
#include "qc.h"

/* **********************************************
 *       	Private functions      		*
 * *********************************************/

void print_text_report_file(bam_qc_report_t bam_qc_report, int base_quality, char *inputfilename, char *outfilename);
void generate_map_errors_histogram_datafile(bam_qc_report_t* bam_qc_report_p, char* map_errors_histogram_filename);
void generate_num_mappings_histogram_datafile(bam_qc_report_t* bam_qc_report_p, char* num_mappings_histogram_filename);
void plot_map_histograms(char* map_errors_histogram_filename, char* graph_filename, char* title, char* xlabel, int column, int x_size);
qc_graph_t* set_qc_graph_default_values(qc_graph_t* qc_graph);
void generate_gnuplot_image(qc_graph_t qc_graph, char* data_filename, char* graph_filename);
void generate_html_file_with_images(char* html_filename, char* report_directory, char* in_shortname, int valid);

/* **************************************************************
 *       	Public function implementations      		*
 * *************************************************************/

void generate_report(bam_qc_report_t bam_qc_report, char* inputfilename, int base_quality, char* report_directory, int valid) {
    char in_shortname[MAX_FULL_PATH_LENGTH];
    char processed_filename[MAX_FULL_PATH_LENGTH];
    char data_filename[MAX_FULL_PATH_LENGTH];
    char map_errors_histogram_filename[MAX_FULL_PATH_LENGTH];
    char num_mappings_histogram_filename[MAX_FULL_PATH_LENGTH];
    char graph_filename[MAX_FULL_PATH_LENGTH];
    char html_filename[MAX_FULL_PATH_LENGTH];
    char* str_valid_suffix = "";

    get_filename_from_path(inputfilename, in_shortname);

    str_valid_suffix = (valid) ? (char*) VALID_ALIGNMENT_FILE_SUFFIX : (char*)  INVALID_ALIGNMENT_FILE_SUFFIX;

    // print text report
    sprintf(data_filename, "%s/%s%s%s", report_directory, in_shortname, QC_SUFFIX, str_valid_suffix);
    sprintf(processed_filename, "%s%s", in_shortname, str_valid_suffix);
    print_text_report_file(bam_qc_report, base_quality, processed_filename, data_filename);

    sprintf(map_errors_histogram_filename, "%s/%s%s%s", report_directory, in_shortname, MAP_ERRORS_HISTOGRAM_SUFFIX, str_valid_suffix);
    sprintf(graph_filename, "%s/%s.map_errors_histogram%s.png", report_directory, in_shortname, str_valid_suffix);
    generate_map_errors_histogram_datafile(&bam_qc_report, map_errors_histogram_filename);
    plot_map_histograms(map_errors_histogram_filename, graph_filename, "Histogram of map errors", "Num. of errors", MAP_ERRORS_COLUMN, MAX_MAP_ERRORS_IN_HISTOGRAM);

    sprintf(graph_filename, "%s/%s.map_deletions_histogram%s.png", report_directory, in_shortname, str_valid_suffix);
    generate_map_errors_histogram_datafile(&bam_qc_report, map_errors_histogram_filename);
    plot_map_histograms(map_errors_histogram_filename, graph_filename, "Histogram of deletions", "Num. of deletions", DELETION_COLUMN, MAX_MAP_ERRORS_IN_HISTOGRAM);

    sprintf(graph_filename, "%s/%s.map_insertions_histogram%s.png", report_directory, in_shortname, str_valid_suffix);
    generate_map_errors_histogram_datafile(&bam_qc_report, map_errors_histogram_filename);
    plot_map_histograms(map_errors_histogram_filename, graph_filename, "Histogram of insertions", "Num. of insertions", INSERTION_COLUMN, MAX_MAP_ERRORS_IN_HISTOGRAM);

    sprintf(graph_filename, "%s/%s.map_matchings_histogram%s.png", report_directory, in_shortname, str_valid_suffix);
    generate_map_errors_histogram_datafile(&bam_qc_report, map_errors_histogram_filename);
    plot_map_histograms(map_errors_histogram_filename, graph_filename, "Histogram of matches", "Num. of matches", MATCHING_COLUMN, MAX_MAP_ERRORS_IN_HISTOGRAM);

    sprintf(num_mappings_histogram_filename, "%s/%s%s%s", report_directory, in_shortname, NUM_MAPPINGS_HISTOGRAM_SUFFIX, str_valid_suffix);
    sprintf(graph_filename, "%s/%s.num_mappings_histogram%s.png", report_directory, in_shortname, str_valid_suffix);
    generate_num_mappings_histogram_datafile(&bam_qc_report, num_mappings_histogram_filename);
    plot_map_histograms(num_mappings_histogram_filename, graph_filename, "Histogram of alignments per read", "Num. of aligments", NUM_ALIGNMENTS_COLUMN, MAX_MAPPING_COUNT_IN_HISTOGRAM);

    sprintf(html_filename, "%s/%s%s%s", report_directory, in_shortname, str_valid_suffix, HTML_FILE_SUFFIX);
    generate_html_file_with_images(html_filename, report_directory, in_shortname, valid);
}

/* **************************************************************
 *       	Private function implementations      		*
 * *************************************************************/

void print_text_report_file(bam_qc_report_t bam_qc_report, int base_quality, char *inputfilename, char *outfilename) {
    FILE* fd = (outfilename == NULL ? stdout : fopen(outfilename, "w"));

    fprintf(fd, "\n----------------------------------------------\n");
    fprintf(fd, "         B A M     Q C     R E P O R T            ");
    fprintf(fd, "\n----------------------------------------------\n");
    fprintf(fd, "\nProcessed file  : %s\n", inputfilename);
    fprintf(fd, "\nNumber of alignments  : %li\n", bam_qc_report.num_alignments);
    fprintf(fd, "\nMean alignment length  : %li\n", bam_qc_report.mean_alignment_length);
    fprintf(fd, "\nMean alignment quality: %i\n", bam_qc_report.mean_map_quality);
    fprintf(fd, "\nStrand 0/1  : %3.2f/%3.2f\n", 100.0 * (bam_qc_report.num_alignments - bam_qc_report.strand_counter) / bam_qc_report.num_alignments, 100.0 * bam_qc_report.strand_counter / bam_qc_report.num_alignments);
    fprintf(fd, "\nMean distance between paired ends: %ld\n", bam_qc_report.mean_paired_end_distance);
    fprintf(fd, "\nNumber of covered nts: %ld\n", nts_with_coverage);   //Global variable
    fprintf(fd, "\nMean coverage: %3.2f\n\n", mean_coverage / nts_with_coverage);  //Global variable

    if (outfilename != NULL)  fclose(fd);
}

void generate_map_errors_histogram_datafile(bam_qc_report_t* bam_qc_report_p, char* map_errors_histogram_filename) {
    FILE* fd = (map_errors_histogram_filename == NULL ? stdout : fopen(map_errors_histogram_filename, "w"));

    for (int i = 0; i <= (MAX_MAP_ERRORS_IN_HISTOGRAM + 1); i++) {
        fprintf(fd, "%i\t%i\t%i\t%i\t%i\n", i, bam_qc_report_p->map_error_histogram[i], bam_qc_report_p->map_deletion_histogram[i], bam_qc_report_p->map_insertion_histogram[i], bam_qc_report_p->map_matching_histogram[i]);
    }

    if (map_errors_histogram_filename != NULL)  fclose(fd);
}

void generate_num_mappings_histogram_datafile(bam_qc_report_t* bam_qc_report_p, char* num_mappings_histogram_filename) {
    FILE* fd = (num_mappings_histogram_filename == NULL ? stdout : fopen(num_mappings_histogram_filename, "w"));

    for (int i = 0; i <= (MAX_MAPPING_COUNT_IN_HISTOGRAM + 1); i++) {
        fprintf(fd, "%i\t%u\n", i, bam_qc_report_p->num_mappings_histogram[i]);
    }

    if (num_mappings_histogram_filename != NULL)  fclose(fd);
}

void plot_map_histograms(char* map_errors_histogram_filename, char* graph_filename, char* title, char* xlabel, int column, int x_size) {
    qc_graph_t qc_graph;
    set_qc_graph_default_values(&qc_graph);

    qc_graph.title = title;
    qc_graph.xlabel = xlabel;
    qc_graph.ylabel = "Num. of alignments";
    qc_graph.type = "boxes";
    qc_graph.x_autoscale = 0;
    qc_graph.x_start = 0;
    qc_graph.x_end = (x_size + 1);
    qc_graph.y_autoscale = 1;
    qc_graph.x_column = POS_COLUMN;
    qc_graph.num_y_columns = 1;
    qc_graph.y_columns[0] = column;

    generate_gnuplot_image(qc_graph, map_errors_histogram_filename, graph_filename);
}

qc_graph_t* set_qc_graph_default_values(qc_graph_t* qc_graph) {
    qc_graph->x_autoscale = 1;
    qc_graph->x_start = 1;
    qc_graph->x_end = 100;
    qc_graph->y_autoscale = 1;
    qc_graph->y_start = 0;
    qc_graph->y_end = 100;
    qc_graph->lmargin = 10;
    qc_graph->rmargin = 4;
    qc_graph->tmargin = 3;
    qc_graph->bmargin = 4;
    qc_graph->title = "Title";
    qc_graph->xlabel = "X axis";
    qc_graph->ylabel = "Y axis";
    qc_graph->type = "lines";
    qc_graph->x_column = 0;
    qc_graph->num_y_columns = 1;
    qc_graph->y_columns[0] = 1;
    qc_graph->y_titles[0] = "";
    qc_graph->y_titles[1] = "";
    qc_graph->y_titles[2] = "";
    qc_graph->y_titles[3] = "";
    qc_graph->y_titles[4] = "";
    qc_graph->y_titles[5] = "";
    qc_graph->y_titles[6] = "";
    qc_graph->y_titles[7] = "";
    qc_graph->y_titles[8] = "";
    qc_graph->y_titles[9] = "";

    return qc_graph;
}

void generate_gnuplot_image(qc_graph_t qc_graph, char* data_filename, char* graph_filename) {
    // lines specifying input data and output graph are declared and filled
    char line[MAX_FULL_PATH_LENGTH];

    char gnuplot_filename[1024];
    sprintf(gnuplot_filename, "%s.gnuplot", graph_filename);

    // open the file for writing gnuplot lines
    FILE* graph_fd = fopen(gnuplot_filename, "w");

    if (graph_fd == NULL) {
        LOG_FATAL("Opening of file descriptor for gnuplot execution failed\n");
        return;
    }

    sprintf(line, "set output '%s'\n", graph_filename);
    fprintf(graph_fd, line);
    fprintf(graph_fd, "set terminal png nocrop enhanced font arial 10 size 640,360\n");
    sprintf(line, "set ylabel '%s'\n", qc_graph.ylabel);
    fprintf(graph_fd, line);
    sprintf(line, "set xlabel '%s'\n", qc_graph.xlabel);
    fprintf(graph_fd, line);
    fprintf(graph_fd, "set ytics border in scale 1,0.5 mirror norotate  offset character 0, 0, 0\n");
    sprintf(line, "set title '%s'\n", qc_graph.title);
    fprintf(graph_fd, line);

    if (qc_graph.x_autoscale == 1) {
        fprintf(graph_fd, "set autoscale x\n");
    } else {
        sprintf(line, "set xrange [ %i : %i ] noreverse nowriteback\n", qc_graph.x_start, qc_graph.x_end);
        fprintf(graph_fd, line);
    }

    if (qc_graph.y_autoscale == 1) {
        fprintf(graph_fd, "set autoscale y\n");
    } else {
        sprintf(line, "set yrange [ %i : %i ] noreverse nowriteback\n", qc_graph.x_start, qc_graph.x_end);
        fprintf(graph_fd, line);
    }

    sprintf(line, "set lmargin '%i'\n", qc_graph.lmargin);
    fprintf(graph_fd, line);
    sprintf(line, "set rmargin '%i'\n", qc_graph.rmargin);
    fprintf(graph_fd, line);
    sprintf(line, "set tmargin '%i'\n", qc_graph.tmargin);
    fprintf(graph_fd, line);
    sprintf(line, "set bmargin '%i'\n", qc_graph.bmargin);
    fprintf(graph_fd, line);

    sprintf(line, "plot ");

    for (int i = 0; i < qc_graph.num_y_columns; i++) {
        sprintf(line, "%s%s '%s' using %i:%i title '%s' with %s", line, (i == 0 ? "" : ", "), data_filename, qc_graph.x_column, qc_graph.y_columns[i], qc_graph.y_titles[i], qc_graph.type);
    }
    fprintf(graph_fd, line);

    fclose(graph_fd);

    // build the command line by calling gnuplot followed by is instruction file
    char cmd[1024];
    sprintf(cmd, "gnuplot %s;", gnuplot_filename);

    //execute command line: gnuplot filename.gnuplot
    system(cmd);
}

void generate_html_file_with_images(char* html_filename, char* report_directory, char* in_shortname, int valid) {
    char line[MAX_FULL_PATH_LENGTH];
    char* str_valid_suffix = NULL;

    str_valid_suffix = (valid) ? (char*) VALID_FILE_SUFFIX : (char*) INVALID_FILE_SUFFIX;

    FILE* fd = fopen(html_filename, "w");

    fprintf(fd, "<HTML>\n");
    fprintf(fd, "<HEAD>\n");
    fprintf(fd, "<TITLE>BAM Quality Control Results</TITLE>\n");
    fprintf(fd, "</HEAD>\n");
    fprintf(fd, "<BODY>\n");
    fprintf(fd, "<TABLE WIDTH='70%'>\n");
    fprintf(fd, "<TR><TD ALIGN='center'><U><B>BAM QUALITY CONTROL RESULTS</B></U></TD></TR>\n");
    fprintf(fd, "<TR><TD ALIGN='center'></TD></TR>\n");
    fprintf(fd, "<TR><TD ALIGN='center'></TD></TR>\n");

    // print general statistics from file
    struct stat buf;

    sprintf(line, "%s/%s%s%s", report_directory, in_shortname, QC_SUFFIX, str_valid_suffix);
    stat(line, &buf);

    char* buf_p = (char*) calloc(1, buf.st_size * sizeof(char));

    FILE* aux_fd = fopen(line, "r");

    fread(buf_p, 1, buf.st_size, aux_fd);
    fclose(aux_fd);

    fprintf(fd, "<TR><TD ALIGN='center'><PRE>\n");
    fprintf(fd, buf_p);
    fprintf(fd, "</PRE></TD></TR>\n");

    // print graphics
    sprintf(line, "<TR><TD ALIGN='center'><IMG src='%s/%s.map_errors_histogram%s.png'</TD></TR>\n", report_directory, in_shortname, str_valid_suffix);
    fprintf(fd, line);

    sprintf(line, "<TR><TD ALIGN='center'><IMG src='%s/%s.map_deletions_histogram%s.png'</TD></TR>\n", report_directory, in_shortname, str_valid_suffix);
    fprintf(fd, line);

    sprintf(line, "<TR><TD ALIGN='center'><IMG src='%s/%s.map_insertions_histogram%s.png'</TD></TR>\n", report_directory, in_shortname, str_valid_suffix);
    fprintf(fd, line);

    sprintf(line, "<TR><TD ALIGN='center'><IMG src='%s/%s.map_matchings_histogram%s.png'</TD></TR>\n", report_directory, in_shortname, str_valid_suffix);
    fprintf(fd, line);

    sprintf(line, "<TR><TD ALIGN='center'><IMG src='%s/%s.num_mappings_histogram%s.png'</TD></TR>\n", report_directory, in_shortname, str_valid_suffix);
    fprintf(fd, line);

    fprintf(fd, "</TABLE>\n");
    fprintf(fd, "</BODY>\n");
    fprintf(fd, "</HTML>\n");

    //free dynamic memory and close file
    free(buf_p);
    fclose(fd);
}
