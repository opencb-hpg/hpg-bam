BIN = ./bin

COMMONS_LIB = libs/commons
COMMONS_CUDA_LIB = libs/commons-cuda
CONTAINERS_LIB = libs/containers
REGIONS_LIB = libs/bioformats/features/region
SQLITE_LIB = sqlite

SAMTOOLS_LIB = samtools-0.1.18
BAM_LIB = libs/bioformats/bam-sam
XML_LIB = /usr/include/libxml2

LIB = libs

NVCC_DISCOVER = $(shell expr `which nvcc | wc -l` \> 10)

ALL = hpg-bam

CC = gcc
#CC = /home/vrequena/third_party/llvm-build/Release+Asserts/bin/clang -faddress-sanitizer -O1 -fno-omit-frame-pointer -g
CXX = g++ -fopenmp 
#CFLAGS = -O3 -Wall -std=c99 -D_XOPEN_SOURCE=500 
CFLAGS = -O0 -g -Wall -std=c99 -D_XOPEN_SOURCE=500
#CFLAGS = -DVERBOSE_DBG -Wall -std=c99 -D_XOPEN_SOURCE=500

CINCLUDES = -I. -I/opt/cuda/include -I$(SAMTOOLS_LIB) -I$(BAM_LIB) -I$(COMMONS_LIB) -I$(COMMONS_CUDA_LIB) -I$(CONTAINERS_LIB) -I$(REGIONS_LIB) -I$(SQLITE_LIB) -lcprops
CUINCLUDES = -I. -I$(SAMTOOLS_LIB) -I$(BAM_LIB) -I$(COMMONS_LIB) -I$(COMMONS_CUDA_LIB) -I$(CONTAINERS_LIB) -I$(SQLITE_LIB)

NVCC = nvcc
NVCCFLAGS = -g -G -Xptxas -v  -arch=sm_20 -Xcompiler " -fopenmp"
#NVCCFLAGS = -g -G -Xptxas -v  -arch=sm_12 -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_OMP -lgomp

all: $(ALL)

SQLITE_OBJECTS = sqlite3.o

ifeq ($(NVCC_DISCOVER), 1)
CUDA_OBJECTS = qc_kernel_cuda.o cuda_commons.o
endif

ifeq ($(NVCC_DISCOVER), 1)

hpg-bam: commons-objects bam_hpc_tools_main.o hpg-bam-objects containers-objects region-objects sqlite-objects samtools-objects qc.o sort.o sort_thrust.o convert.o filter.o map_validation.o $(CUDA_OBJECTS)
	$(NVCC) $(NVCCFLAGS) $(COMMONS_LIB)/*.o bam_hpc_tools_main.o aligner_dataset.o aligner_dataset_file.o bam_file.o bam_reader.o bam_writer.o map_validation.o\
        bam_qc_batch.o bam_qc_report.o bam_coverage.o chrom_alignments.o convert.o qc.o qc_hash.o qc_hash_list.o list.o sort.o bam_trie.o mappings_db.o\
        bam_data_batch.o sort_thrust.o qc_kernel_omp.o gff_data.o gff_reader.o alignment.o region_table.o region.o dna_map_region.o filter.o\
        GeneralHashFunctions.o $(CUDA_OBJECTS) $(SQLITE_OBJECTS) -o $(BIN)/hpg-bam -L$(LIB) -lbam -lz -I$(XML_LIB) -lxml2 -lcurl -lcprops

else

hpg-bam: commons-objects bam_hpc_tools_main.o hpg-bam-objects containers-objects region-objects sqlite-objects samtools-objects qc.o sort.o sort_thrust.o convert.o filter.o map_validation.o $(CUDA_OBJECTS)
	$(CXX) $(CFLAGS) $(COMMONS_LIB)/*.o bam_hpc_tools_main.o aligner_dataset.o aligner_dataset_file.o bam_file.o bam_reader.o bam_writer.o map_validation.o\
        bam_qc_batch.o bam_qc_report.o bam_coverage.o chrom_alignments.o convert.o qc.o qc_hash.o qc_hash_list.o list.o sort.o bam_trie.o mappings_db.o\
        bam_data_batch.o sort_thrust.o qc_kernel_omp.o gff_data.o gff_reader.o alignment.o region_table.o region.o dna_map_region.o filter.o\
        GeneralHashFunctions.o $(SQLITE_OBJECTS) -o $(BIN)/hpg-bam -L$(LIB) -lbam -lz -I$(XML_LIB) -lxml2 -lcurl -lcprops

endif

commons-objects: clean-commons-objects
	(cd $(COMMONS_LIB); $(MAKE) file_utils.o system_utils.o string_utils.o log.o)

convert.o: convert.c convert.h *.h
	$(CC) $(CFLAGS) -c convert.c $(CINCLUDES)

filter.o: filter.c filter.h *.h
	$(CC) $(CFLAGS) -c filter.c $(CINCLUDES)

ifeq ($(NVCC_DISCOVER), 1)

qc.o: qc.cu qc.h *.h
	$(NVCC) $(NVCCFLAGS) -DCUDA_VERSION -c qc.cu $(CUINCLUDES)

sort.o: sort.cu sort.h *.h
	$(NVCC) $(NVCCFLAGS) -DCUDA_VERSION -c sort.cu $(CUINCLUDES)

sort_thrust.o: sort_thrust.cu *.h
	$(NVCC) $(NVCCFLAGS) -DCUDA_VERSION -c sort_thrust.cu $(CUINCLUDES)

qc_kernel_cuda.o: qc_kernel_cuda.cu *.h
	$(NVCC) $(NVCCFLAGS) -DCUDA_VERSION -c qc_kernel_cuda.cu $(CINCLUDES)

cuda_commons.o: $(COMMONS_CUDA_LIB)/cuda_commons.cu $(COMMONS_CUDA_LIB)/cuda_commons.h *.h
	$(NVCC) $(NVCCFLAGS) -DCUDA_VERSION -c $(COMMONS_CUDA_LIB)/cuda_commons.cu

else

qc.o: qc.c qc.h *.h
	$(CC) $(CFLAGS) -c qc.c $(CINCLUDES)

sort.o: sort.c sort.h *.h
	$(CC) $(CFLAGS) -c sort.c $(CINCLUDES)

sort_thrust.o: sort_thrust.c *.h
	$(CXX) $(CFLAGS) -c sort_thrust.c -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_OMP $(CINCLUDES)

endif

map_validation.o: map_validation.c map_validation.h *.h
	$(CC) $(CFLAGS) -c map_validation.c $(CINCLUDES)

bam_hpc_tools_main.o: bam_hpc_tools_main.c *.h
	$(CC) $(CFLAGS) -c bam_hpc_tools_main.c $(CINCLUDES)

containers-objects:
	 $(CC) $(CFLAGS) -c ${CONTAINERS_LIB}/*.c $(CINCLUDES) -lcprops

region-objects:
	$(CC) $(CFLAGS) -c ${REGIONS_LIB}/*.c $(CINCLUDES) -lcprops 

bam-objects:
	$(CC) $(CFLAGS) -c $(BAM_LIB)/*.c $(CINCLUDES)  -fopenmp

sqlite-objects:
	$(CC) $(CFLAGS) -c $(SQLITE_LIB)/*.c $(CINCLUDES)  -fopenmp

samtools-objects: clean-samtools-objects
	(cd $(SAMTOOLS_LIB); $(MAKE))

hpg-bam-objects: bam-objects
	$(CC) $(CFLAGS) -c bam_trie.c bam_qc_batch.c aligner_dataset.c aligner_dataset_file.c bam_qc_report.c bam_coverage.c chrom_alignments.c \
mappings_db.c gff_data.c gff_reader.c qc_hash.c qc_hash_list.c qc_kernel_omp.c GeneralHashFunctions.c $(CINCLUDES) -fopenmp


###################################################################

clean-commons-objects: 
	(cd $(COMMONS_LIB); $(MAKE) clean)

clean-samtools-objects: 
	(cd $(SAMTOOLS_LIB); $(MAKE) clean)

clean-hpg-bam:
	-rm -f *~ \#*\# $(LIB)/*.o $(SAMTOOLS_LIB)/*.o $(BIN)/hpg-bam

clean:
	-rm -f *~ \#*\# *.o $(SAMTOOLS_LIB)/*.o $(BIN)/*

