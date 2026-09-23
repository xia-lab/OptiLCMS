#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List lipids_build_dot_buffer_cpp(
    NumericVector mass1,
    NumericVector intensity1,
    NumericVector mass2,
    NumericVector intensity2,
    double bin,
    double min_mz,
    double max_mz,
    bool drive_reference = false) {

  const int n1 = mass1.size();
  const int n2 = mass2.size();

  if (n1 == 0 || n2 == 0) {
    return List::create(
      _["mz"] = NumericVector(0),
      _["m"] = NumericVector(0),
      _["r"] = NumericVector(0),
      _["baseM"] = 0.0,
      _["baseR"] = 0.0,
      _["size"] = 0
    );
  }

  const double last_mass1 = mass1[n1 - 1];
  const double last_mass2 = mass2[n2 - 1];
  const double last_mass = (last_mass1 > last_mass2) ? last_mass1 : last_mass2;

  double focused_mz = min_mz;
  int remain_m = 0;
  int remain_l = 0;

  NumericVector buf_mz(n1 + n2);
  NumericVector buf_m(n1 + n2);
  NumericVector buf_r(n1 + n2);
  int size = 0;

  double base_m = -std::numeric_limits<double>::max();
  double base_r = -std::numeric_limits<double>::max();

  while (focused_mz <= max_mz) {
    double sum_l = 0.0;
    int i = remain_l;
    while (i < n2) {
      const double mz = mass2[i];
      if (mz < focused_mz - bin) {
        i++;
        continue;
      }
      if (mz < focused_mz + bin) {
        sum_l += intensity2[i];
        i++;
      } else {
        break;
      }
    }
    remain_l = i;

    double sum_m = 0.0;
    i = remain_m;
    while (i < n1) {
      const double mz = mass1[i];
      if (mz < focused_mz - bin) {
        i++;
        continue;
      }
      if (mz < focused_mz + bin) {
        sum_m += intensity1[i];
        i++;
      } else {
        break;
      }
    }
    remain_m = i;

    buf_mz[size] = focused_mz;
    buf_m[size] = sum_m;
    buf_r[size] = sum_l;
    if (sum_m > base_m) base_m = sum_m;
    if (sum_l > base_r) base_r = sum_l;
    size++;

    if (drive_reference) {
      if (focused_mz + bin > last_mass2) break;
      if (remain_l < n2) {
        focused_mz = mass2[remain_l];
      } else {
        break;
      }
    } else {
      if (focused_mz + bin > last_mass) break;

      const double next_m = (remain_m < n1) ? mass1[remain_m] : R_PosInf;
      const double next_l = (remain_l < n2) ? mass2[remain_l] : R_PosInf;

      if (next_l > focused_mz + bin && next_m <= focused_mz + bin) {
        focused_mz = next_m;
      } else if (next_l <= focused_mz + bin && next_m > focused_mz + bin) {
        focused_mz = next_l;
      } else {
        focused_mz = (next_m < next_l) ? next_m : next_l;
      }
    }
  }

  if (size == 0) {
    return List::create(
      _["mz"] = NumericVector(0),
      _["m"] = NumericVector(0),
      _["r"] = NumericVector(0),
      _["baseM"] = 0.0,
      _["baseR"] = 0.0,
      _["size"] = 0
    );
  }

  return List::create(
    _["mz"] = buf_mz[Range(0, size - 1)],
    _["m"] = buf_m[Range(0, size - 1)],
    _["r"] = buf_r[Range(0, size - 1)],
    _["baseM"] = base_m,
    _["baseR"] = base_r,
    _["size"] = size
  );
}
