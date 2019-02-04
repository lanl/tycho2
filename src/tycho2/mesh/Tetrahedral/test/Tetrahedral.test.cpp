#define CATCH_CONFIG_MAIN
#include "catch2/catch.hpp"
#include "tycho2.hpp"

SCENARIO("GaussLegendre can return points and weights "
	 "of the corresponding order", "[GaussLegendre]"){
  namespace mesh = tycho2::mesh;
  auto myMesh = mesh::make<mesh::Tetrahedral>();
}
