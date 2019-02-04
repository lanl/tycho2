template<typename CellType, typname... Ts>
class Cell{

  // a cell seems like a bunch of facets...
  // i.e. spatial and structural information.
  // Not sure if we should store material
  // iformation because it could be redundant.
  // Probably want things like volumes and
  // and surface normals cached if they are used
  // along with the spatial information when sweeping.
  
};
