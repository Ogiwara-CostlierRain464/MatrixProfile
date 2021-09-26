#include <gtest/gtest.h>
#include "../src/stomp.h"

using namespace std;

class Unit: public ::testing::Test{};

TEST_F(Unit, SlidingDotProduct){
  std::vector<double> a{1,1,1,1};
  std::vector<double> b{1,1,1,1};
  CArray out;
  SlidingDotProduct(a, b , out);
  out[0].real();
}

int main(int argc, char **argv){
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}