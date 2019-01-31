#define CATCH_CONFIG_MAIN
#include "catch2/catch.hpp"
#include "tycho2.hpp"

SCENARIO("GaussLegendre can return points and weights "
	 "of the corresponding order", "[GaussLegendre]"){
  namespace quad = tycho2::quadrature;

  GIVEN("A Gauss-Legendre Rule of second order"){

    auto gauss_legendre = quad::GaussLegendre<2>{};
    
    WHEN("Querying for points and weights"){

      auto ref_pts = std::array<double,2>{-0.5773502691896257,
					  0.5773502691896257};
      auto pts = gauss_legendre.points();
	
      REQUIRE(pts == ref_pts);
      
      auto ref_wts = std::array<double,2>{1.0,1.0};
      
      auto wts = gauss_legendre.weights();

      REQUIRE(wts == ref_wts);

      //AND_THEN("using the points and weights to compute the "
      //	       "integral of f(x)=x^2 on (-1,1)"){
	
      //}

    }
    
  }
}
