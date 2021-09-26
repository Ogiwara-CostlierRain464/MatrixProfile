#include <iostream>
#include <vector>
#include <cfloat>
#include <cassert>
#include <cmath>
#include "stomp.h"

int main(int argc, char* argv[]){
  if(argc < 5){
    printf("Usage: STOMP <window_len> <input_file> <profile output file> <index output file>\n");
    exit(0);
  }

  size_t m = atoi(argv[1]);
  std::vector<double> T;
  readFile(argv[2], T, "%lf");
  size_t n = T.size();
  size_t l = n - m + 1;
  std::vector<float> profile(l, CC_MIN);
  std::vector<uint32_t> profile_idx(l, 0);
  std::cout << "Starting STOMP" << std::endl;


  return 0;
}