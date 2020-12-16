# the compiler: gcc for C program
CC = gcc

#compiler flags:
# -g adds debugging information to the executable file
# -Wall turns on most, but not all, compiler warnings
CFLAGS = -w -pg
DBGCFLAGS = -g -w -fopenmp
# -lm links the math library
#LIBS = -lm -lpthread -lz
LIBS = -lm -pthread -lz
OPTIMIZATION = -O3 -march=native
#sources
SOURCES = NJcluster.c options.c math.c opt.c
NEEDLEMANWUNSCH = needleman_wunsch.c alignment.c alignment_scoring.c
HASHMAP = hashmap.c
KALIGN = kalign/run_kalign.c kalign/tlmisc.c kalign/tldevel.c kalign/parameters.c kalign/rwalign.c kalign/alignment_parameters.c kalign/idata.c kalign/aln_task.c kalign/bisectingKmeans.c kalign/esl_stopwatch.c kalign/aln_run.c kalign/alphabet.c kalign/pick_anchor.c kalign/sequence_distance.c kalign/euclidean_dist.c kalign/aln_mem.c kalign/tlrng.c kalign/aln_setup.c kalign/aln_controller.c kalign/weave_alignment.c kalign/aln_seqseq.c kalign/aln_profileprofile.c kalign/aln_seqprofile.c
OBJECTS = (SOURCES: .c = .o)
# the build target executable:
TARGET = NJcluster

all: $(TARGET)
$(TARGET): $(TARGET).c
	$(CC) $(CFLAGS) $(OPTIMIZATION) -o $(TARGET) $(NEEDLEMANWUNSCH) $(HASHMAP) $(SOURCES) $(LIBS)
debug: $(TARGET).c
	$(CC) $(DBGCFLAGS) -o $(TARGET) $(NEEDLEMANWUNSCH) $(HASHMAP) $(KALIGN) $(SOURCES) $(LIBS)

clean:
	$(RM) $(TARGET)

