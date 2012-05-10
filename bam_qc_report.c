
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "bam_qc_report.h"
#include "qc.h"



// private functions for datafile and graph generation
//
void print_text_report_file(bam_qc_report_t bam_qc_report, int base_quality, char *inputfilename, char *outfilename);
void generate_map_errors_histogram_datafile(bam_qc_report_t* bam_qc_report_p, char* map_errors_histogram_filename);
void generate_num_mappings_histogram_datafile(bam_qc_report_t* bam_qc_report_p, char* num_mappings_histogram_filename);
void plot_map_histograms(char* map_errors_histogram_filename, char* graph_filename, char* title, char* xlabel, int column, int x_size);
// void generate_position_datafile(qc_report_t qc_report, char* filename);
// void generate_read_quality_datafile(qc_report_t* qc_report_p, char* quality_filename);
// void generate_kmers_datafile(qc_report_t* qc_report_p, char* kmers_filename);
// void generate_kmers_position_datafile(qc_report_t* qc_report_p, char* kmers_position_filename);
// void plot_nt_quality(qc_report_t* qc_report_p, char* datafilename, char* graphfilename);
// void plot_read_quality(char* quality_filename, char* graph_filename);
// void plot_sequence_length_distribution(char* data_filename, char* graph_filename);
// void plot_base_sequence_content(qc_report_t* qc_report_p, char* data_filename, char* graph_filename);
// void plot_base_gc_content(qc_report_t* qc_report_p, char* data_filename, char* graph_filename);
// void generate_gc_histogram_datafile(qc_report_t* qc_report_p, char* gc_histogram_filename);
// void plot_sequence_gc_content(char* data_filename, char* graph_filename);
// void plot_n_content(qc_report_t* qc_report_p, char* data_filename, char* graph_filename);
qc_graph_t* set_qc_graph_default_values(qc_graph_t* qc_graph);
void generate_gnuplot_image(qc_graph_t qc_graph, char* data_filename, char* graph_filename);
void generate_html_file_with_images(char* html_filename, char* report_directory, char* in_shortname, int valid);


//-----------------------------------------------------
// generate_report (public)
//-----------------------------------------------------

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
	//
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


//-----------------------------------------------------
// print_text_report_file
//-----------------------------------------------------

void print_text_report_file(bam_qc_report_t bam_qc_report, int base_quality, char *inputfilename, char *outfilename) {

  FILE* fd = (outfilename==NULL ? stdout : fopen(outfilename, "w"));
  
  fprintf(fd, "\n----------------------------------------------\n");
  fprintf(fd, "         B A M     Q C     R E P O R T            ");
  fprintf(fd, "\n----------------------------------------------\n");
  fprintf(fd, "\nProcessed file  : %s\n", inputfilename);
  fprintf(fd, "\nNumber of alignments  : %li\n", bam_qc_report.num_alignments);
  fprintf(fd, "\nMean alignment length  : %li\n", bam_qc_report.mean_alignment_length);
  //fprintf(fd, "\nMean alignment quality: %i\n", (bam_qc_report.mean_map_quality - base_quality));
  fprintf(fd, "\nMean alignment quality: %i\n", bam_qc_report.mean_map_quality);
  fprintf(fd, "\nStrand 0/1  : %3.2f/%3.2f\n", 100.0 * (bam_qc_report.num_alignments - bam_qc_report.strand_counter) / bam_qc_report.num_alignments, 100.0 * bam_qc_report.strand_counter / bam_qc_report.num_alignments);  
  fprintf(fd, "\nMean distance between paired ends: %ld\n", bam_qc_report.mean_paired_end_distance);
  fprintf(fd, "\nNumber of covered nts: %ld\n", nts_with_coverage);			//Global variable
  fprintf(fd, "\nMean coverage: %3.2f\n\n", mean_coverage / nts_with_coverage);		//Global variable
  
  if (outfilename!=NULL)  fclose(fd);
}


//-----------------------------------------------------
// generate_map_errors_histogram_datafile
//-----------------------------------------------------

void generate_map_errors_histogram_datafile(bam_qc_report_t* bam_qc_report_p, char* map_errors_histogram_filename) {

  FILE* fd = (map_errors_histogram_filename==NULL ? stdout : fopen(map_errors_histogram_filename, "w"));
  
  for (int i=0; i<=(MAX_MAP_ERRORS_IN_HISTOGRAM + 1); i++) {
    fprintf(fd, "%i\t%i\t%i\t%i\t%i\n", i, bam_qc_report_p->map_error_histogram[i], bam_qc_report_p->map_deletion_histogram[i], bam_qc_report_p->map_insertion_histogram[i], bam_qc_report_p->map_matching_histogram[i]);
  }
  
  if (map_errors_histogram_filename!=NULL)  fclose(fd);  
}

//-----------------------------------------------------
// generate_map_errors_histogram_datafile
//-----------------------------------------------------

void generate_num_mappings_histogram_datafile(bam_qc_report_t* bam_qc_report_p, char* num_mappings_histogram_filename) {

  FILE* fd = (num_mappings_histogram_filename==NULL ? stdout : fopen(num_mappings_histogram_filename, "w"));
  
  for (int i=0; i<=(MAX_MAPPING_COUNT_IN_HISTOGRAM + 1); i++) {
    fprintf(fd, "%i\t%u\n", i, bam_qc_report_p->num_mappings_histogram[i]);
  }
  
  if (num_mappings_histogram_filename!=NULL)  fclose(fd);  
}


//-----------------------------------------------------
// plot_sequence_gc_content
//-----------------------------------------------------

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



/*
//-----------------------------------------------------
// generate positon datafile
//-----------------------------------------------------

void generate_position_datafile(qc_report_t qc_report, char* filename) {

	FILE* fd = (filename==NULL ? stdout : fopen(filename, "w"));
	
	int nt_last_position = qc_report.nt_counter[0];
	int length;

	for (int i=1; i<=qc_report.max_read_length; i++) {
	  
	      if (nt_last_position != qc_report.nt_counter[i]) {
		length = nt_last_position - qc_report.nt_counter[i];
		nt_last_position = qc_report.nt_counter[i];
	      } else {
		length = 0;
	      }   
	      
	      fprintf(fd, "%i\t%i\t%i\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\n", i, qc_report.mean_nt_quality[i-1], length, (double) (100.0 * qc_report.nt_type_counter[i-1][A])/qc_report.nt_counter[i-1], (double) (100.0 * qc_report.nt_type_counter[i-1][C])/qc_report.nt_counter[i-1], (double) (100.0 * qc_report.nt_type_counter[i-1][G])/qc_report.nt_counter[i-1], (double) (100.0 * qc_report.nt_type_counter[i-1][T])/qc_report.nt_counter[i-1], (double) (100.0 * qc_report.nt_type_counter[i-1][N])/qc_report.nt_counter[i-1], ((100.0 * (double) qc_report.nt_type_counter[i-1][C])/qc_report.nt_counter[i-1] + (100.0 * (double) qc_report.nt_type_counter[i-1][G])/qc_report.nt_counter[i-1]));
	}
	
	if (fd!=NULL)  fclose(fd);
}

//-----------------------------------------------------
// plot_nt_quality
//-----------------------------------------------------

void plot_nt_quality(qc_report_t* qc_report_p, char* data_filename, char* graph_filename) {
    
    qc_graph_t qc_graph;
    set_qc_graph_default_values(&qc_graph);

    qc_graph.title = "Quality per nucleotide position";
    qc_graph.xlabel = "nt position";
    qc_graph.ylabel = "Quality (normalized Prhed scale)";
    qc_graph.type = "lines";
    qc_graph.x_autoscale = 0;
    qc_graph.x_start = 0;
    qc_graph.x_end = qc_report_p->max_read_length;
    qc_graph.y_autoscale = 1;
    qc_graph.x_column = POS_COLUMN;
    qc_graph.num_y_columns = 1;
    qc_graph.y_columns[0] = NT_QUALITY_COLUMN;

    generate_gnuplot_image(qc_graph, data_filename, graph_filename);
    
}

//-----------------------------------------------------
// plot_sequence_length_distribution
//-----------------------------------------------------

void plot_sequence_length_distribution(char* data_filename, char* graph_filename) {

    qc_graph_t qc_graph;
    set_qc_graph_default_values(&qc_graph);

    qc_graph.title = "Sequence length distribution";
    qc_graph.xlabel = "Sequence length";
    qc_graph.ylabel = "Num. of reads";
    qc_graph.type = "boxes";
    qc_graph.x_autoscale = 1;
    qc_graph.y_autoscale = 1;
    qc_graph.x_column = POS_COLUMN;
    qc_graph.num_y_columns = 1;
    qc_graph.y_columns[0] = NUM_READ_LENGTH_COLUMN;    
 
    generate_gnuplot_image(qc_graph, data_filename, graph_filename);  
}

//-----------------------------------------------------
// plot_base_sequence_content
//-----------------------------------------------------

void plot_base_sequence_content(qc_report_t* qc_report_p, char* data_filename, char* graph_filename) {
    
    qc_graph_t qc_graph;
    set_qc_graph_default_values(&qc_graph);

    qc_graph.title = "Per base sequence content";
    qc_graph.xlabel = "nt position";
    qc_graph.ylabel = "Percentage of nt (%)";
    qc_graph.type = "lines";
    qc_graph.x_autoscale = 0;
    qc_graph.x_start = 0;
    qc_graph.x_end = qc_report_p->max_read_length;
    qc_graph.y_autoscale = 1;
    qc_graph.x_column = POS_COLUMN;
    qc_graph.num_y_columns = 4;
    qc_graph.y_columns[0] = A_COLUMN;
    qc_graph.y_columns[1] = C_COLUMN;
    qc_graph.y_columns[2] = G_COLUMN;
    qc_graph.y_columns[3] = T_COLUMN;
    qc_graph.y_titles[0] = "A %";
    qc_graph.y_titles[1] = "C %";
    qc_graph.y_titles[2] = "G %";
    qc_graph.y_titles[3] = "T %";

    generate_gnuplot_image(qc_graph, data_filename, graph_filename);
    
}

//-----------------------------------------------------
// plot_base_gc_content
//-----------------------------------------------------

void plot_base_gc_content(qc_report_t* qc_report_p, char* data_filename, char* graph_filename) {
    
    qc_graph_t qc_graph;
    set_qc_graph_default_values(&qc_graph);

    qc_graph.title = "Per base GC content";
    qc_graph.xlabel = "nt position";
    qc_graph.ylabel = "GC %";
    qc_graph.type = "lines";
    qc_graph.x_autoscale = 0;
    qc_graph.x_start = 0;
    qc_graph.x_end = qc_report_p->max_read_length;
    qc_graph.y_autoscale = 1;
    qc_graph.x_column = POS_COLUMN;
    qc_graph.num_y_columns = 1;
    qc_graph.y_columns[0] = GC_COLUMN;

    generate_gnuplot_image(qc_graph, data_filename, graph_filename);
    
}

//-----------------------------------------------------
// generate_gc_histogram_datafile
//-----------------------------------------------------

void generate_gc_histogram_datafile(qc_report_t* qc_report_p, char* gc_histogram_filename) {

	FILE* fd = (gc_histogram_filename==NULL ? stdout : fopen(gc_histogram_filename, "w"));
	
	int i;
	for (i=0; i<=100; i++) {
	        fprintf(fd, "%i\t%i\n", i, qc_report_p->gc_histogram[i]);
	}
	
	if (gc_histogram_filename!=NULL)  fclose(fd);    
  
}

//-----------------------------------------------------
// plot_sequence_gc_content
//-----------------------------------------------------

void plot_sequence_gc_content(char* gc_histogram_filename, char* graph_filename) {

    qc_graph_t qc_graph;
    set_qc_graph_default_values(&qc_graph);

    qc_graph.title = "Per sequence GC content";
    qc_graph.xlabel = "GC %";
    qc_graph.ylabel = "Num. of reads";
    qc_graph.type = "boxes";
    qc_graph.x_autoscale = 0;
    qc_graph.x_start = 0;
    qc_graph.x_end = 99;
    qc_graph.y_autoscale = 1;
    qc_graph.x_column = POS_COLUMN;
    qc_graph.num_y_columns = 1;
    qc_graph.y_columns[0] = GC_HISTOGRAM_COLUMN;

    generate_gnuplot_image(qc_graph, gc_histogram_filename, graph_filename);
}

//-----------------------------------------------------
// plot_n_content
//-----------------------------------------------------

void plot_n_content(qc_report_t* qc_report_p, char* data_filename, char* graph_filename) {
    
    qc_graph_t qc_graph;
    set_qc_graph_default_values(&qc_graph);

    qc_graph.title = "Per base sequence content";
    qc_graph.xlabel = "nt position";
    qc_graph.ylabel = "Percentage of N (%)";
    qc_graph.type = "lines";
    qc_graph.x_autoscale = 0;
    qc_graph.x_start = 0;
    qc_graph.x_end = qc_report_p->max_read_length;
    qc_graph.y_autoscale = 1;
    qc_graph.x_column = POS_COLUMN;
    qc_graph.num_y_columns = 1;
    qc_graph.y_columns[0] = N_COLUMN;

    generate_gnuplot_image(qc_graph, data_filename, graph_filename);
    
}

//-----------------------------------------------------
// plot_kmers_position
//-----------------------------------------------------

void plot_kmers_position(qc_report_t* qc_report_p, char* data_filename, char* graph_filename) {

    qc_graph_t qc_graph;
    set_qc_graph_default_values(&qc_graph);

    qc_graph.title = "Relative enrichment over read length";
    qc_graph.xlabel = "nt position";
    qc_graph.ylabel = "Num. of kmers";
    qc_graph.type = "lines";
    qc_graph.x_autoscale = 0;
    qc_graph.x_start = 0;
    qc_graph.x_end = qc_report_p->max_read_length;
    qc_graph.y_autoscale = 1;
    qc_graph.x_column = KMER_POS_COLUMN;
    qc_graph.num_y_columns = 5;
    qc_graph.y_columns[0] = KMER_1_COLUMN;
    qc_graph.y_columns[1] = KMER_2_COLUMN;
    qc_graph.y_columns[2] = KMER_3_COLUMN;
    qc_graph.y_columns[3] = KMER_4_COLUMN;
    qc_graph.y_columns[4] = KMER_5_COLUMN;
    qc_graph.y_titles[0] = qc_report_p->qc_kmers[0].kmer;
    qc_graph.y_titles[1] = qc_report_p->qc_kmers[1].kmer;
    qc_graph.y_titles[2] = qc_report_p->qc_kmers[2].kmer;
    qc_graph.y_titles[3] = qc_report_p->qc_kmers[3].kmer;
    qc_graph.y_titles[4] = qc_report_p->qc_kmers[4].kmer;

    generate_gnuplot_image(qc_graph, data_filename, graph_filename);
  
}

//-----------------------------------------------------
// get_reads_per_length_from_nt_counter
//-----------------------------------------------------

void get_reads_per_length_from_nt_counter(int nt_counter[], int* length_p) {
	
	int nt_last_position = nt_counter[0];
	  
	for (int j=0; j<MAX_LINE_LENGTH ; j++) {
	      if (nt_last_position != nt_counter[j]) {
		length_p[j] = nt_last_position - nt_counter[j];
		nt_last_position = nt_counter[j];
	      } else {
		length_p[j] = 0;
	      }
	}

	for (int j=0; j<MAX_LINE_LENGTH ; j++) {
	    printf("length_p[%i]: %i\n", j, length_p[j]);
	}

}
*/

//-----------------------------------------------------
// set_qc_graph_default_values
//-----------------------------------------------------

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

//-----------------------------------------------------
// generate_gnuplot_image
//-----------------------------------------------------

// void generate_gnuplot_image(qc_graph_t qc_graph, char* data_filename, char* graph_filename) {
// 
// 	// lines specifying input data and output graph are declared and filled
// 	//
// 	char line[MAX_FULL_PATH_LENGTH];
// 	
// 	// graph is parametrized based on qc_graph options and plotted
// 	//
// 	FILE *graph_fd = popen("gnuplot -persist","w");
// 
// 	if (graph_fd == NULL) {
// 	  printf("ERROR opening pipe\n");
// 	}
// 
// 	sprintf(line, "set output '%s'\n", graph_filename);
// 	fprintf(graph_fd, line);
// 	fprintf(graph_fd, "set terminal png nocrop enhanced font arial 10 size 640,360\n");
// 	sprintf(line, "set ylabel '%s'\n", qc_graph.ylabel);
// 	fprintf(graph_fd, line);
// 	sprintf(line, "set xlabel '%s'\n", qc_graph.xlabel);
// 	fprintf(graph_fd, line);
// 	fprintf(graph_fd, "set ytics border in scale 1,0.5 mirror norotate  offset character 0, 0, 0\n");
// 	sprintf(line, "set title '%s'\n", qc_graph.title);
// 	fprintf(graph_fd, line);	
// 
// 	
// 	if (qc_graph.x_autoscale == 1) {
// 	    fprintf(graph_fd, "set autoscale x\n");    
// 	} else {
// 	    sprintf(line, "set xrange [ %i : %i ] noreverse nowriteback\n", qc_graph.x_start, qc_graph.x_end);
// 	    fprintf(graph_fd, line);
// 	}
// 
// 	if (qc_graph.y_autoscale == 1) {
// 	    fprintf(graph_fd, "set autoscale y\n");    
// 	} else {
// 	    sprintf(line, "set yrange [ %i : %i ] noreverse nowriteback\n", qc_graph.x_start, qc_graph.x_end);
// 	    fprintf(graph_fd, line);
// 	}
// 
// 	sprintf(line, "set lmargin '%i'\n", qc_graph.lmargin);
// 	fprintf(graph_fd, line);
// 	sprintf(line, "set rmargin '%i'\n", qc_graph.rmargin);
// 	fprintf(graph_fd, line);
// 	sprintf(line, "set tmargin '%i'\n", qc_graph.tmargin);
// 	fprintf(graph_fd, line);
// 	sprintf(line, "set bmargin '%i'\n", qc_graph.bmargin);
// 	fprintf(graph_fd, line);
// 
// 	sprintf(line, "plot ");
// 	for (int i=0; i<qc_graph.num_y_columns; i++) {
// 	    sprintf(line, "%s%s '%s' using %i:%i title '%s' with %s", line, (i==0 ? "" : ", "), data_filename, qc_graph.x_column, qc_graph.y_columns[i], qc_graph.y_titles[i], qc_graph.type);
// //    	    sprintf(line, "%s%s '%s' using %i:%i notitle with %s", line, (i==0 ? "" : ", "), data_filename, qc_graph.x_column, qc_graph.y_columns[i], qc_graph.type);
// // 	    sprintf(line, "%s%s '%s' using %i:%i %s with %s", line, (i==0 ? "" : ", "), data_filename, qc_graph.x_column, qc_graph.y_columns[i], (strcmp(qc_graph.y_titles[i], "")!=0 ? "title 'TITULO'" : "notitle"), qc_graph.type);
// 	}	
// 	fprintf(graph_fd, line);
// 
// 	pclose(graph_fd); 
// }

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

//-----------------------------------------------------
// generate_html_file_with_images
//-----------------------------------------------------

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

	//sprintf(line, "<TR><TD ALIGN='center'><a href='%s/%s.qc' target='blank'>Summary results</a></TD></TR>\n", report_directory, in_shortname);
	//fprintf(fd, line);	

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