import struct

UINT64_MAX = 18446744073709551615

with open('tet.smesh', 'wb') as f:
    f.write(struct.pack('Q', 1))    # num cells
    f.write(struct.pack('Q', 4))    # num faces
    f.write(struct.pack('Q', 4))    # num nodes
    
    # Cell Data
    # bounding faces
    f.write(struct.pack('Q', 0))
    f.write(struct.pack('Q', 1))
    f.write(struct.pack('Q', 2))
    f.write(struct.pack('Q', 3))
    
    # bounding nodes
    f.write(struct.pack('Q', 0))
    f.write(struct.pack('Q', 1))
    f.write(struct.pack('Q', 2))
    f.write(struct.pack('Q', 3))
    
    
    # Face Data
    # bounding cells (0)
    f.write(struct.pack('Q', 0))
    f.write(struct.pack('Q', UINT64_MAX))
    
    # bounding nodes (0)
    f.write(struct.pack('Q', 1))
    f.write(struct.pack('Q', 2))
    f.write(struct.pack('Q', 3))
    
    # bounding cells (1)
    f.write(struct.pack('Q', 0))
    f.write(struct.pack('Q', UINT64_MAX))
    
    # bounding nodes (1)
    f.write(struct.pack('Q', 0))
    f.write(struct.pack('Q', 2))
    f.write(struct.pack('Q', 3))
    
    # bounding cells (2)
    f.write(struct.pack('Q', 0))
    f.write(struct.pack('Q', UINT64_MAX))
    
    # bounding nodes (2)
    f.write(struct.pack('Q', 0))
    f.write(struct.pack('Q', 1))
    f.write(struct.pack('Q', 3))
    
    # bounding cells (3)
    f.write(struct.pack('Q', 0))
    f.write(struct.pack('Q', UINT64_MAX))
    
    # bounding nodes (3)
    f.write(struct.pack('Q', 0))
    f.write(struct.pack('Q', 1))
    f.write(struct.pack('Q', 2))
    
    
    # Node Data
    # coords (0)
    f.write(struct.pack('d', 0.0))
    f.write(struct.pack('d', 0.0))
    f.write(struct.pack('d', 0.0))
    
    # coords (1)
    f.write(struct.pack('d', 1.0))
    f.write(struct.pack('d', 0.0))
    f.write(struct.pack('d', 0.0))
    
    # coords (2)
    f.write(struct.pack('d', 0.5))
    f.write(struct.pack('d', 1.0))
    f.write(struct.pack('d', 0.0))
    
    # coords (3)
    f.write(struct.pack('d', 0.5))
    f.write(struct.pack('d', 0.5))
    f.write(struct.pack('d', 1.0))