##### Settings for local environment #####
include make.inc


##### Variables #####
CPP += -DPARALLEL_MESH_READ_SERIAL_ONLY

SRC_MESH = ../src/SerialMesh.cc \
           ../src/ParallelMesh.cc

INC = -I../src \
      -I$(MOAB_DIR)/include \
      -I$(METIS_DIR)/include

LIB_MOAB = $(MOAB_DIR)/lib/libMOAB.a \
       $(MOAB_DIR)/lib/libiMesh.a \
       $(CGNS_DIR)/lib/libcgns.a \
       $(HDF5_DIR)/lib/libhdf5.a \
       -ldl -lz 

LIB_METIS = $(METIS_DIR)/lib/libmetis.a


##### Applications #####
APPS = SerialToParallelMesh \
       PartitionColumns \
       MoabToSerialMesh \
       RefineSerialMesh \
       SerialMeshToMoab \
       ParallelMeshToMoab
       


##### Make Commands #####
# all
.PHONY: all
all: $(APPS)

# clean
.PHONY: clean
clean: 
	@echo Delete executables
	rm *.x
	@echo " "

# SerialToParallelMesh
.PHONY: SerialToParallelMesh
SerialToParallelMesh:
	@echo Making SerialToParallelMesh
	$(CPP) $(INC) $(SRC_MESH) SerialToParallelMesh.cc -o SerialToParallelMesh.x $(LIB_METIS)
	@echo " "

# PartitionColumns
.PHONY: PartitionColumns
PartitionColumns:
	@echo Making PartitionColumns
	$(CPP) $(INC) $(SRC_MESH) PartitionColumns.cc -o PartitionColumns.x
	@echo " "

# MoabToSerialMesh
.PHONY: MoabToSerialMesh
MoabToSerialMesh:
	@echo Making MoabToSerialMesh
	$(CPP) $(INC) $(SRC_MESH) MoabToSerialMesh.cc -o MoabToSerialMesh.x $(LIB_MOAB)
	@echo " "

# RefineSerialMesh
.PHONY: RefineSerialMesh
RefineSerialMesh:
	@echo Making RefineSerialMesh
	$(CPP) $(SRC_MESH) RefineSerialMesh.cc -o RefineSerialMesh.x
	@echo " "

# SerialMeshToMoab
.PHONY: SerialMeshToMoab
SerialMeshToMoab:
	@echo Making SerialMeshToMoab
	$(CPP) $(INC) $(SRC_MESH) SerialMeshToMoab.cc -o SerialMeshToMoab.x $(LIB_MOAB)
	@echo " "

# ParallelMeshToMoab
.PHONY: ParallelMeshToMoab
ParallelMeshToMoab:
	@echo Making ParallelMeshToMoab
	$(CPP) $(INC) $(SRC_MESH) ParallelMeshToMoab.cc -o ParallelMeshToMoab.x $(LIB_MOAB)
	@echo " "

