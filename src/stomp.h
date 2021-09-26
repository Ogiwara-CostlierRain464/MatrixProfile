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


inline void compute_sum(const std::vector<double> &T,
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
void SlidingDotProduct(const std::vector<double> &Q,
                       const std::vector<double> &T,
                       std::vector<double> &QT){
  assert(T.size() > Q.size());
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
  CArray tmp;
  ElementwiseMultiplication(Q_raf, T_af, tmp);
  ifft(tmp);
  QT = std::vector<double>(n - m + 1);
  for(size_t i = m; i < n; i++){
    QT[i - m] = tmp[i].real();
  }
}

void CalculateDistanceProfile(const std::vector<double> &QT,
                              const std::vector<double> &sum,
                              size_t m,
                              size_t i,
                              std::vector<double> &out_D){
  // at here, we don't skip trivial match
  // return D[i]
  double mu_i = mean(sum, i, m);
  double std_i = Std(sum, i, m, mu_i);

  out_D = std::vector<double>(QT.size());

  for(size_t j = 0; j < QT.size(); j++){
    double mu_j = mean(sum, i, m);
    double std_j = Std(sum, i, m, mu_j);

    double d_i_j = sqrt(2 * m * (
        1 - ((QT[j] - m * mu_i * mu_j) / (m * std_i * std_j))
      ));
    out_D[j] = d_i_j;
  }
}

void ElementWiseMin(std::vector<double> &P,
                    std::vector<uint32_t> &I,
                    // d_i_0, d_i_1, ...
                    const std::vector<double> &D,
                    size_t i,
                    size_t m){
  // if D it is smaller than P, update P
  // expect it is within trivial match zone.
  assert(P.size() == I.size());
  assert(P.size() == D.size());
  for(size_t j = 0; j < P.size() ; j++){
    if((i - (m/2) < j) and (j < i + (m/2))){
      // Skip trivial match.
      continue;
    }
    if(D[j] < P[j]){
      P[j] = D[j];
      I[j] = i;
    }
  }
}

void STOMP(std::vector<double> &T,
           std::vector<double> &P,
           std::vector<uint32_t> &I,
           size_t m){
  size_t n = T.size(), l = n - m + 1;
  // mu, sigma
  std::vector<double> sum;
  compute_sum(T, sum);
  std::vector<double> QT;
  SlidingDotProduct(
    std::vector<double>(T.begin(), T.begin()+m),
      T, QT);
  std::vector<double> QT_first(QT.begin(), QT.end());
  std::vector<double> D;
  CalculateDistanceProfile(QT, sum, m, 0, D);
  P = D; I = std::vector<uint32_t>(l, 0);
  for(size_t i = 1 ; i < l ; i++){
    for(size_t j = l-1; j >= 1; j--){
      QT[j] = QT[j-1] - T[j-1]*T[i-1] + T[j+m-1]*T[i+m-1];
    }
    QT[0] = QT_first[i];
    CalculateDistanceProfile(QT, sum, m, i, D);
    ElementWiseMin(P, I, D, i, m);
  }
}

void read_file(const char *filename,
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
