namespace mesh{
  
template<typename MeshTag, typename... Ts>
class Cell;

template<typename MeshTag, typename... Args>
auto make(Args&&... args){
  return make(MeshTag{}, std::forward<Args>(args)...);
};

#include "mesh/Tetrahedral.hpp"

// free functions
template<typename PointT>
constexpr static auto cross_product(PointT&& a, PointT&& b)
{
  return make_array (a[1] * b[2] - a[2] * b[1],
		     a[2] * b[0] - a[0] * b[2],
		     a[0] * b[1] - a[1] * b[0]);
}

} // namespace mesh
