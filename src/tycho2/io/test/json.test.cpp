#define CATCH_CONFIG_MAIN
#include "catch2/catch.hpp"
#include "nlohmann/json.hpp"
#include <fstream>

SCENARIO("Integration testing for json library" "[json]"){
using json = nlohmann::json;
 
 GIVEN("a json input file for the gold standard problem"){

   std::ifstream input("gold-standard.json");
   json json_object;
   input >> json_object;
   
   REQUIRE(json_object["snOrder"] == 8);
   REQUIRE(json_object["iterMax"] == 100);
   REQUIRE(json_object["errMax"] == 1e-10);
   REQUIRE(json_object["maxCellsPerStep"] == 100);
   REQUIRE(json_object["intraAngleP"] == 3);
   REQUIRE(json_object["interAngleP"] == 1);
   REQUIRE(json_object["nGroups"] == 2);
   REQUIRE(json_object["sigmaT1"] == 10);
   REQUIRE(json_object["sigmaS1"] == 5);
   REQUIRE(json_object["sigmaT2"] == 10);
   REQUIRE(json_object["sigmaS2"] == 5);
   REQUIRE(json_object["OutputFile"] == true);
   REQUIRE(json_object["OutputFilename"] == "out.psi");
   REQUIRE(json_object["SourceIteration"] == true);
   REQUIRE(json_object["MPIType"] == "TychoTwoSided");
   REQUIRE(json_object["DD_IterMax"] == 100);
   REQUIRE(json_object["DD_ErrMax"] == 1e-10);
   REQUIRE(json_object["SweepType"] == "TraverseGraph");
   REQUIRE(json_object["GaussElim"] == "NoPivot");   
 }
}
