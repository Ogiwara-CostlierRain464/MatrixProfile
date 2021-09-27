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
  read_file(argv[2], T, "%lf");
  size_t n = T.size();
  size_t l = n - m + 1;
  std::vector<double> P(l, -FLT_MAX);
  std::vector<uint32_t> I(l, 0);
  std::cout << "Starting STOMP" << std::endl;
  STOMP(T, P, I, m);
  std::cout << "Now writing result to files" << std::endl;
  FILE* f1 = fopen( argv[3], "w");
  FILE* f2 = fopen( argv[4], "w");
  for(size_t i = 0; i < P.size(); ++i){
    fprintf(f1, "%f\n", P[i]);
    fprintf(f2, "%u\n", I[i]);
  }
  fclose(f1); fclose(f2);
  std::cout << "Done" << std::endl;

  return 0;
}