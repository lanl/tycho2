# Settings for local environment
include make.inc

# Add assert level and PETSc flag to compile command
MPICC += -DASSERT_ON=$(ASSERT_ON) -DUSE_PETSC=$(USE_PETSC)
MPICXX += -DASSERT_ON=$(ASSERT_ON) -DUSE_PETSC=$(USE_PETSC)

# Include source directory
INC += -Isrc

# Add PETSC include directory and library command
ifeq ($(strip $(USE_PETSC)), 1)
	INC += $(PETSC_INC)
	LIBS += $(PETSC_LIB)
endif

# List of sources, header files, and object files
SOURCE = $(wildcard src/*.cc)
OBJS = $(patsubst src%, $(BUILD_DIR)%, $(SOURCE)) 
OBJECTS = $(OBJS:%.cc=%.o)

CSOURCE = $(wildcard src/*.c)
COBJS = $(patsubst src%, $(BUILD_DIR)%, $(CSOURCE)) 
OBJECTS += $(COBJS:%.c=%.o)

HEADERS = $(wildcard src/*.hh)
HEADERS += $(wildcard src/*.h)

# Link object files
sweep.x: $(OBJECTS) 
	@echo Linking $@
	$(MPICXX) $(OBJECTS) -o sweep.x ${LIBS}

$(BUILD_DIR)/%.o: src/%.c $(HEADERS)
	@echo Making $@
	$(MPICC) $(INC) -c $< -o $@

# Make object files
$(BUILD_DIR)/%.o: src/%.cc $(HEADERS)
	@echo Making $@
	$(MPICXX) $(INC) -c $< -o $@

# Delete object files
.PHONY: clean
clean:
	@echo Delete object files
	rm -f $(BUILD_DIR)/*.o *.x

print-%  : ; @echo $* = $($*)

