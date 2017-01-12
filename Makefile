# Settings for local environment
include make.inc


# Add assert level and PETSc flag to compile command
MPICC += -DASSERT_ON=$(ASSERT_ON) -DUSE_PETSC=$(USE_PETSC)


# Include source directory
INC += -Isrc


# Add PETSC include directory and library command
ifeq ($(USE_PETSC), 1)
	INC += $(PETSC_INC)
	LIBS += $(PETSC_LIB)
endif



# List of sources, header files, and object files
SOURCE = $(wildcard src/*.cc)
HEADERS = $(wildcard src/*.hh)
OBJECTS = $(patsubst src%.cc, build%.o, $(SOURCE))


# Link object files
sweep.x: $(OBJECTS)
	@echo Linking $@
	$(MPICC) $(OBJECTS) -o sweep.x ${LIBS}

# Make object files
build/%.o: src/%.cc $(HEADERS) make.inc
	@echo Making $@
	$(MPICC) $(INC) -c $< -o $@

# Delete object files
.PHONY: clean
clean:
	@echo Delete object files
	rm build/*.o

