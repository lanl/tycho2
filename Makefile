# Settings for local environment
include make.inc


# Add assert level and PETSc flag to compile command
MPICC += -DASSERT_ON=$(ASSERT_ON)


# Source and library info
INC += -Isrc
INC += -I$(KOKKOS_DIR)/include
LIBS = $(KOKKOS_DIR)/lib/libkokkos.a


# List of sources, header files, and object files
SOURCE = $(wildcard src/*.cc)
HEADERS = $(wildcard src/*.hh)
OBJECTS = $(patsubst src/%.cc, build/%.o, $(SOURCE))


# Link object files
sweep.x: $(OBJECTS) 
	@echo Linking $@
	$(MPICC) $(OBJECTS) -o sweep.x $(LIBS) -ldl


# Make object files
build/%.o: src/%.cc $(HEADERS) make.inc
	@echo Making $@
	$(MPICC) $(INC) -c $< -o $@


# Delete object files
.PHONY: clean
clean:
	@echo Delete object files
	rm build/*.o sweep.x

