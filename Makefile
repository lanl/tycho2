# Settings for local environment
include make.inc


# Add assert level and PETSc flag to compile command
MPICC += -DASSERT_ON=$(ASSERT_ON) -DUSE_PETSC=$(USE_PETSC)


# Include source directory
SRC_DIR = $(TOP_DIR)/src
INC += -I$(SRC_DIR)


# Add PETSC include directory and library command
ifeq ($(strip $(USE_PETSC)), 1)
	INC += $(PETSC_INC)
	LIBS += $(PETSC_LIB)
endif



# List of sources, header files, and object files
SOURCE = $(wildcard $(SRC_DIR)/*.cc)
HEADERS = $(wildcard $(SRC_DIR)/*.hh)
OBJECTS = $(patsubst $(SRC_DIR)/%.cc, %.o, $(SOURCE))


# Link object files
sweep.x: $(OBJECTS)
	@echo Linking $@
	$(MPICC) $(OBJECTS) -o sweep.x ${LIBS}


# Make object files
%.o: $(SRC_DIR)/%.cc $(HEADERS) make.inc
	@echo Making $@
	$(MPICC) $(INC) -c $< -o $@


# Delete object files
.PHONY: clean
clean:
	@echo Delete object files
	rm *.o *.x

