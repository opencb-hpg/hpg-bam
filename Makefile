BIN = ./bin

COMMONS_LIB = libs/commons
COMMONS_CUDA_LIB = libs/commons-cuda
CONTAINERS_LIB = libs/containers
REGIONS_LIB = libs/bioformats/features/region

BAM_LIB = libs/bioformats/bam-sam
XML_LIB = /usr/include/libxml2

LIB = libs


ALL = hpg-bam

CC = gcc
CXX = g++ -fopenmp 
CFLAGS =  -O3 -Wall -std=c99
#CFLAGS = -DVERBOSE_DBG -Wall -std=c99
#CFLAGS = -Wall -std=c99

CINCLUDES = -I. -I/opt/cuda/include -I./samtools-0.1.18/ -I$(BAM_LIB) -I$(COMMONS_LIB) -I$(COMMONS_CUDA_LIB) -I$(CONTAINERS_LIB) -I$(REGIONS_LIB) -lcprops
CUINCLUDES = -I. -I./samtools-0.1.18/ -I$(BAM_LIB) -I$(COMMONS_LIB) -I$(COMMONS_CUDA_LIB) -I$(CONTAINERS_LIB)

NVCC = nvcc
NVCCFLAGS = -g -G -Xptxas -v  -arch=sm_20 -Xcompiler " -fopenmp"
#NVCCFLAGS = -g -G -Xptxas -v  -arch=sm_12

all: $(ALL)

hpg-bam: commons-objects bam_hpc_tools_main.o hpg-bam-objects containers-objects region-objects qc.o sort.o sort_thrust.o convert.o qc_kernel_cuda.o cuda_commons.o
	$(NVCC) $(NVCCFLAGS) $(COMMONS_LIB)/*.o bam_hpc_tools_main.o aligner_dataset.o aligner_dataset_file.o bam_file.o bam_reader.o bam_writer.o \
        bam_data_batch_list.o bam_qc_batch.o bam_qc_report.o bam_coverage.o chrom_alignments.o convert.o qc.o qc_hash.o qc_hash_list.o list.o sort.o \
        cuda_commons.o bam_data_batch.o sort_thrust.o qc_kernel_cuda.o qc_kernel_omp.o gff_data.o gff_reader.o alignment.o region_table.o region.o \
        GeneralHashFunctions.o -o $(BIN)/hpg-bam -L$(LIB) -lbam -lz -I$(XML_LIB) -lxml2 -lcurl -lcprops

commons-objects: 
	(cd $(COMMONS_LIB); $(MAKE) file_utils.o system_utils.o string_utils.o log.o)

qc.o: qc.cu qc.h *.h
	$(NVCC) $(NVCCFLAGS) -c qc.cu $(CUINCLUDES)

sort.o: sort.cu sort.h *.h
	$(NVCC) $(NVCCFLAGS) -c sort.cu $(CUINCLUDES)

convert.o: convert.c convert.h *.h
	$(CC) $(CFLAGS) -c convert.c $(CINCLUDES)

sort_thrust.o: sort_thrust.cu *.h
	$(NVCC) $(NVCCFLAGS) -c sort_thrust.cu $(CUINCLUDES)

qc_kernel_cuda.o: qc_kernel_cuda.cu *.h
	$(NVCC) $(NVCCFLAGS) -c qc_kernel_cuda.cu $(CINCLUDES)

cuda_commons.o: $(COMMONS_CUDA_LIB)/cuda_commons.cu $(COMMONS_CUDA_LIB)/cuda_commons.h *.h
	$(NVCC) $(NVCCFLAGS) -c $(COMMONS_CUDA_LIB)/cuda_commons.cu

bam_hpc_tools_main.o: bam_hpc_tools_main.c *.h
	$(CC) $(CFLAGS) -c bam_hpc_tools_main.c $(CINCLUDES)

containers-objects:
	 $(CC) $(CFLAGS) -c ${CONTAINERS_LIB}/*.c $(CINCLUDES) -lcprops

region-objects:
	$(CC) $(CFLAGS) -c ${REGIONS_LIB}/region.c $(CINCLUDES) -lcprops 

bam-objects:
	$(CC) $(CFLAGS) -c $(BAM_LIB)/*.c $(CINCLUDES)  -fopenmp

hpg-bam-objects: bam-objects
	$(CC) $(CFLAGS) -c bam_qc_batch.c aligner_dataset.c aligner_dataset_file.c bam_qc_report.c bam_coverage.c chrom_alignments.c \
gff_data.c gff_reader.c qc_hash.c qc_hash_list.c qc_kernel_omp.c GeneralHashFunctions.c $(CINCLUDES)  -fopenmp


###################################################################

clean-commons-objects: 
	(cd $(COMMONS_LIB); $(MAKE) clean)

clean-hpg-bam:
	-rm -f *~ \#*\# $(LIB)/*.o $(BIN)/hpg-bam

clean:
	-rm -f *~ \#*\# *.o $(BIN)/*

