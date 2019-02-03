template<typename MeshTag>
class Face;

template<typename MeshTag, typename Args...>
class Type;

#include "mesh/Tetrahedral.hpp"

// free helper functions
template<typename PointT>
constexpr static auto cross_product(PointT&& a, PointT&& b)
{
  return details::make_array (a[1] * b[2] - a[2] * b[1],
			      a[2] * b[0] - a[0] * b[2],
			      a[0] * b[1] - a[1] * b[0]);
}

