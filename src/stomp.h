#ifndef MATRIXPROFILE_STOMP_H
#define MATRIXPROFILE_STOMP_H

#include <cfloat>
#include <cstddef>
#include "fft.h"

#define CC_MIN -FLT_MAX

union mp_entry{
  float floats[2];
  uint32_t ints[2];
  uint64_t ulong;

  void from(float distance, uint32_t idx){
    floats[0] = distance;
    ints[1] = idx;
  }
};


inline void ComputeSum(const std::vector<double> &T,
                       std::vector<double> &out_sum){
  out_sum = std::vector<double>(T.size(), 0);
  for(size_t i = 0; i < out_sum.size(); i++){
    out_sum[i] = out_sum[i] + T[i];
  }
}


inline double mean(const std::vector<double> &sum,
                   size_t i,
                   size_t m){
  // μ_i is the mean of T_i,m
  double coeff = 1.0 / (double) m;
  size_t i_plus_m = i + m;
  size_t n = sum.size();

  if(i == 0){
    return sum[m - 1] * coeff;
  }else{
    assert(i < n - 1);
    return (sum[i_plus_m-1] - sum[i-1]) * coeff;
  }
}

inline double Std(const std::vector<double> &sum,
                  size_t i,
                  size_t m,
                  double mean){
  // σ_i is the standard deviation of T_i,m
  double coeff = 1.0 / (double) m;
  size_t i_plus_m = i + m;
  size_t n = sum.size();
  if(i == 0){
    return 1 / sqrt((sum[m-1] * coeff) - mean * mean);
  }else{
    assert(i_plus_m < n + m);
    return 1 / sqrt(((sum[i_plus_m-1] - sum[i-1]) * coeff) - (mean * mean));
  }
}

void ElementwiseMultiplication(const CArray &Q_raf, const CArray &T_af, CArray &out){
  assert(Q_raf.size() == T_af.size());
  out = CArray(Q_raf.size());
  for(size_t i = 0; i < Q_raf.size(); i++){
    out[i] = Q_raf[i] * T_af[i];
  }
}

/**
 * @param Q: A query Q.
 * @param T: user provided time series T.
 */
void SlidingDotProduct(std::vector<double> &Q, std::vector<double> &T, CArray &out){
  size_t n = T.size(), m = Q.size();
  // copy
  std::vector<double> T_a(T.begin(), T.end());
  T_a.resize(2 * n, 0);
  std::vector<double> Q_ra(Q.rbegin(), Q.rend());
  Q_ra.resize(2 * n, 0);
  CArray Q_raf, T_af;
  from_double_vec(Q_ra, Q_raf);
  from_double_vec(T_a, T_af);
  fft(Q_raf); fft(T_af);
  ElementwiseMultiplication(Q_raf, T_af, out);
  ifft(out);
}

void STOMP(std::vector<double> &T,
           std::vector<float> &profile,
           std::vector<uint32_t> &profile_idx,
           size_t m){
  size_t n = T.size();
  size_t l = n - m + 1;
  // mu, sigma
  // どこかで必ず計算O(1)

}

void readFile(const char *filename,
              std::vector<double>& v,
              const char *format_str){
  FILE *f = fopen(filename, "r");
  if(f == nullptr){
    std::cerr << "Unable to open:" << filename << std::endl;
    exit(-1);
  }
  double num;
  while (!feof(f)){
    fscanf(f, format_str, &num);
    v.push_back(num);
  }
  v.pop_back();
  fclose(f);
}


#endif //MATRIXPROFILE_STOMP_H
