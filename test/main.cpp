#include <gtest/gtest.h>
#include "../src/stomp.h"

using namespace std;

class Unit: public ::testing::Test{};

TEST_F(Unit, SlidingDotProduct){
  std::vector<double> Q{1,1};
  std::vector<double> T{1,1,1,1};
  std::vector<double> QT;
  SlidingDotProduct(Q, T , QT);

}

int main(int argc, char **argv){
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}