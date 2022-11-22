# the compiler: gcc for C program
CC = gcc

#compiler flags:
# -g adds debugging information to the executable file
# -Wall turns on most, but not all, compiler warnings
CFLAGS = -w -pg
DBGCFLAGS = -g -w -fopenmp
# -lm links the math library
#LIBS = -lm -lpthread -lz
LIBS = -lm -pthread -lz -std=gnu99
OPENMP = -fopenmp
OPTIMIZATION = -O3 -march=native
#sources
SOURCES = ancestralclust.c options.c math.c opt.c
NEEDLEMANWUNSCH = needleman_wunsch.c alignment.c alignment_scoring.c
HASHMAP = hashmap.c
KALIGN = kalign/run_kalign.c kalign/tlmisc.c kalign/tldevel.c kalign/parameters.c kalign/rwalign.c kalign/alignment_parameters.c kalign/idata.c kalign/aln_task.c kalign/bisectingKmeans.c kalign/esl_stopwatch.c kalign/aln_run.c kalign/alphabet.c kalign/pick_anchor.c kalign/sequence_distance.c kalign/euclidean_dist.c kalign/aln_mem.c kalign/tlrng.c kalign/aln_setup.c kalign/aln_controller.c kalign/weave_alignment.c kalign/aln_seqseq.c kalign/aln_profileprofile.c kalign/aln_seqprofile.c
WFA2 = WFA2/wavefront_aligner.c WFA2/wavefront_bialigner.c WFA2/mm_allocator.c WFA2/wavefront_align.c WFA2/wavefront_attributes.c WFA2/wavefront_debug.c WFA2/profiler_timer.c WFA2/profiler_counter.c WFA2/cigar.c WFA2/wavefront_compute.c WFA2/wavefront_slab.c WFA2/wavefront.c WFA2/vector.c WFA2/wavefront_components.c WFA2/bitmap.c WFA2/wavefront_backtrace_buffer.c WFA2/wavefront_pcigar.c WFA2/wavefront_heuristic.c WFA2/wavefront_penalties.c WFA2/wavefront_unialign.c WFA2/wavefront_plot.c WFA2/heatmap.c WFA2/wavefront_backtrace.c WFA2/wavefront_extend.c WFA2/wavefront_compute_affine2p.c WFA2/wavefront_backtrace_offload.c WFA2/wavefront_compute_affine.c WFA2/wavefront_compute_linear.c WFA2/wavefront_compute_edit.c WFA2/string_padded.c WFA2/wavefront_bialign.c
OBJECTS = (SOURCES: .c = .o)
# the build target executable:
TARGET = ancestralclust

all: $(TARGET)
$(TARGET): $(TARGET).c
	$(CC) $(OPENMP) -o $(TARGET) $(NEEDLEMANWUNSCH) $(HASHMAP) $(KALIGN) $(WFA2) $(SOURCES) $(LIBS)
debug: $(TARGET).c
	$(CC) $(DBGCFLAGS) -o $(TARGET) $(NEEDLEMANWUNSCH) $(HASHMAP) $(KALIGN) $(WFA2) $(SOURCES) $(LIBS)

clean:
	$(RM) $(TARGET)

