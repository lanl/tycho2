# Kokkos stuff
KOKKOS_PATH = /home/ckgarrett/projects/kokkos
KOKKOS_SRC_PATH = ${KOKKOS_PATH}

default: sweep.x
	echo "Start Build"

CXX = ${KOKKOS_PATH}/bin/nvcc_wrapper
CXXFLAGS = -O3 -Wall -Wextra -Wpedantic -DASSERT_ON=0 -Isrc

KOKKOS_DEVICES = "Cuda,OpenMP"
KOKKOS_ARCH = "SNB,Kepler35"
KOKKOS_CUDA_OPTIONS += "enable_lambda"

DEPFLAGS = -M

# List of sources, header files, and object files
SOURCE = $(wildcard src/*.cc)
HEADERS = $(wildcard src/*.hh)
OBJECTS = $(patsubst src/%.cc, build/%.o, $(SOURCE))

include $(KOKKOS_PATH)/Makefile.kokkos


# Settings for local environment
#MPICC = ${KOKKOS_PATH}/bin/nvcc_wrapper -O3 -std=c++11 -fopenmp -Wall -Wextra -Wpedantic


# Link object files
sweep.x: $(OBJECTS) $(KOKKOS_LINK_DEPENDS)
	@echo Linking $@
	$(CXX) $(KOKKOS_LDFLAGS) $(OBJECTS) -o sweep.x $(KOKKOS_LIBS) $(LIBS) -ldl


# Make object files
build/%.o: src/%.cc $(HEADERS) make.inc $(KOKKOS_CPP_DEPENDS)
	@echo Making $@
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) -c $< -o $@


# Delete object files
.PHONY: clean
clean:
	@echo Delete object files
	rm -f build/*.o sweep.x *.o *.a *.h *.tmp

