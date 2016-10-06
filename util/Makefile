##### Settings for local environment #####
include make.inc


##### Variables #####
CPP1 = $(CPP) -DPARALLEL_MESH_READ_SERIAL_ONLY

SRC_MESH = ../src/SerialMesh.cc \
           ../src/ParallelMesh.cc

SRC_MESH2 = ../src/SerialMesh.cc \
            ../src/ParallelMesh.cc \
            ../src/Comm.cc \
            ../src/Assert.cc

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
APPS = PartitionMetis \
       SerialMeshInfo \
       ParallelMeshInfo \
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

# PartitionMetis
.PHONY: PartitionMetis
PartitionMetis:
	@echo Making PartitionMetis
	$(CPP1) $(INC) $(SRC_MESH) PartitionMetis.cc -o PartitionMetis.x $(LIB_METIS)
	@echo " "

# SerialMeshInfo
.PHONY: SerialMeshInfo
SerialMeshInfo:
	@echo Making SerialMeshInfo
	$(CPP1) $(INC) $(SRC_MESH) SerialMeshInfo.cc -o SerialMeshInfo.x
	@echo " "

# ParallelMeshInfo
.PHONY: ParallelMeshInfo
ParallelMeshInfo:
	@echo Making ParallelMeshInfo
	$(CPP) $(INC) $(SRC_MESH2) ParallelMeshInfo.cc -o ParallelMeshInfo.x
	@echo " "

# PartitionColumns
.PHONY: PartitionColumns
PartitionColumns:
	@echo Making PartitionColumns
	$(CPP1) $(INC) $(SRC_MESH) PartitionColumns.cc -o PartitionColumns.x
	@echo " "

# MoabToSerialMesh
.PHONY: MoabToSerialMesh
MoabToSerialMesh:
	@echo Making MoabToSerialMesh
	$(CPP1) $(INC) $(SRC_MESH) MoabToSerialMesh.cc -o MoabToSerialMesh.x $(LIB_MOAB)
	@echo " "

# RefineSerialMesh
.PHONY: RefineSerialMesh
RefineSerialMesh:
	@echo Making RefineSerialMesh
	$(CPP1) $(SRC_MESH) RefineSerialMesh.cc -o RefineSerialMesh.x
	@echo " "

# SerialMeshToMoab
.PHONY: SerialMeshToMoab
SerialMeshToMoab:
	@echo Making SerialMeshToMoab
	$(CPP1) $(INC) $(SRC_MESH) SerialMeshToMoab.cc -o SerialMeshToMoab.x $(LIB_MOAB)
	@echo " "

# ParallelMeshToMoab
.PHONY: ParallelMeshToMoab
ParallelMeshToMoab:
	@echo Making ParallelMeshToMoab
	$(CPP1) $(INC) $(SRC_MESH) ParallelMeshToMoab.cc -o ParallelMeshToMoab.x $(LIB_MOAB)
	@echo " "

