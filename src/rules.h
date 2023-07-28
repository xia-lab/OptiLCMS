#ifndef RULES_H
#define RULES_H

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace std;

// Names of all rules
vector<string> static getBio_rules_name() {
  const vector<string> bio_rules_name = { // 18
    "O-NH", "NH3-O", "H2", "CH2", "NH", 
    "O", "NH3", "H2O", "CO", "C2H4", 
    "CH2O","S", "H2S", "C2H2O", "CO2", 
    "C5H8", "SO3", "HPO3"
  };
  return bio_rules_name;
}

vector<string> static getAdc_rules_name() {
  const vector<string> adc_rules_name = { // 22
    "Na-H", "K-H", "HCOOH", "HCOONa", 
    "HCOOK", "C2H4O2", "CH3COONa", "CH3COOK",  "H2SO4", 
    "HNO3", "NaNO3", "H2CO3", "NaHCO3", "KHCO3", 
    "H3PO4", "2Na-2H", "2K-2H", "HCN", "CH3OH", 
    "NaOH", "KOH", "CH3CN"
  };
  return adc_rules_name;
}

vector<string> static getFrg_rules_name() {
  const vector<string> frg_rules_name = { // 5
    "CO2", "CH20", "H20", "NH3", "H"
  };
  return frg_rules_name;
}

vector<string> static getAll_rules_name() {
  const vector<string> all_rules_name = { // 45
    "O-NH", "NH3-O", "H2", "CH2", "NH",
    "O", "NH3", "H2O", "CO", "C2H4",
    "CH2O","S", "H2S", "C2H2O", "CO2",
    "C5H8", "SO3", "HPO3",
    "Na-H", "K-H", "HCOOH", "HCOONa",
    "HCOOK", "C2H4O2", "CH3COONa", "CH3COOK", "H2SO4", 
    "HNO3", "NaNO3", "H2CO3", "NaHCO3", "KHCO3", 
    "H3PO4", "2Na-2H", "2K-2H", "HCN", "CH3OH",
    "NaOH", "KOH", "CH3CN",
    "CO2", "CH20", "H20", "NH3", "H"
  };
  return all_rules_name;
}

// formula vector of all rules
static vector<int> * getBio_formulas(){
  static vector<int> bio_formulas[18] = {
    // C H O N P S Na K
    {0,-1,1,-1,0,0,0,0}, {0,3,-1,1,0,0,0,0}, {0,2,0,0,0,0,0,0}, {1,2,0,0,0,0,0,0}, {0,1,0,1,0,0,0,0}, 
    {0,0,1,0,0,0,0,0}, {0,3,0,1,0,0,0,0}, {0,2,1,0,0,0,0,0}, {1,0,1,0,0,0,0,0}, {2,0,4,0,0,0,0,0}, 
    {1,2,1,0,0,0,0,0}, {0,0,0,0,0,1,0,0}, {0,2,0,0,0,1,0,0}, {2,2,1,0,0,0,0,0}, {1,0,2,0,0,0,0,0}, 
    {5,8,0,0,0,0,0,0}, {0,0,3,0,0,1,0,0}, {0,1,3,0,1,0,0,0}
  };
  return bio_formulas;
}

static vector<int> * getAdc_formulas(){
  static vector<int> adc_formulas[22] = {
    // C H O N P S Na K
    {0,-1,0,0,0,0,1,0}, {0,-1,0,0,0,0,0,1}, {1,2,2,0,0,0,0,0}, {1,1,2,0,0,0,1,0}, 
    {1,1,2,0,0,0,0,1}, {2,4,2,0,0,0,0,0}, {2,3,2,0,0,0,1,0}, {2,3,2,0,0,0,0,1}, {0,2,4,0,0,1,0,0}, 
    {0,1,3,1,0,0,0,0}, {0,0,3,1,0,0,1,0}, {1,2,3,0,0,0,0,0}, {1,1,3,0,0,0,1,0}, {1,1,3,0,0,0,0,1}, 
    {0,3,4,0,1,0,0,0}, {0,-2,0,0,0,0,2,0}, {0,-2,0,0,0,0,0,2}, {1,1,0,1,0,0,0,0}, {1,4,1,0,0,0,0,0}, 
    {0,1,1,0,0,0,1,0}, {0,1,1,0,0,0,0,1}, {2,3,0,1,0,0,0,0}
  };
  return adc_formulas;
}

static vector<int> * getFrg_formulas(){
  static vector<int> frg_formulas[5] = {
    // C H O N P S Na K
    {1,0,2,0,0,0,0,0}, {1,2,1,0,0,0,0,0}, {0,2,1,0,0,0,0,0}, {0,3,0,1,0,0,0,0}, {1,0,0,0,0,0,0,0}
  };
  return frg_formulas;
}

static vector<int> * getAll_formulas(){
  static vector<int> all_formulas[45] = {
    // C H O N P S Na K
    {0,-1,1,-1,0,0,0,0}, {0,3,-1,1,0,0,0,0}, {0,2,0,0,0,0,0,0}, {1,2,0,0,0,0,0,0}, {0,1,0,1,0,0,0,0}, 
    {0,0,1,0,0,0,0,0}, {0,3,0,1,0,0,0,0}, {0,2,1,0,0,0,0,0}, {1,0,1,0,0,0,0,0}, {2,0,4,0,0,0,0,0}, 
    {1,2,1,0,0,0,0,0}, {0,0,0,0,0,1,0,0}, {0,2,0,0,0,1,0,0}, {2,2,1,0,0,0,0,0}, {1,0,2,0,0,0,0,0}, 
    {5,8,0,0,0,0,0,0}, {0,0,3,0,0,1,0,0}, {0,1,3,0,1,0,0,0},
    {0,-1,0,0,0,0,1,0}, {0,-1,0,0,0,0,0,1}, {1,2,2,0,0,0,0,0}, {1,1,2,0,0,0,1,0}, 
    {1,1,2,0,0,0,0,1}, {2,4,2,0,0,0,0,0}, {2,3,2,0,0,0,1,0}, {2,3,2,0,0,0,0,1}, {0,2,4,0,0,1,0,0}, 
    {0,1,3,1,0,0,0,0}, {0,0,3,1,0,0,1,0}, {1,2,3,0,0,0,0,0}, {1,1,3,0,0,0,1,0}, {1,1,3,0,0,0,0,1}, 
    {0,3,4,0,1,0,0,0}, {0,-2,0,0,0,0,2,0}, {0,-2,0,0,0,0,0,2}, {1,1,0,1,0,0,0,0}, {1,4,1,0,0,0,0,0}, 
    {0,1,1,0,0,0,1,0}, {0,1,1,0,0,0,0,1}, {2,3,0,1,0,0,0,0},
    {1,0,2,0,0,0,0,0}, {1,2,1,0,0,0,0,0}, {0,2,1,0,0,0,0,0}, {0,3,0,1,0,0,0,0}, {1,0,0,0,0,0,0,0}
  };
  return all_formulas;
}

// mass different of all rules
vector<double> static getBio_ms_change (){
  const vector<double> bio_ms_change {
    0.9840, 1.0316, 2.0157, 14.0157, 15.0109,
    15.9949, 17.0265, 18.0106, 27.9949, 28.0313, 
    30.0106, 31.9720, 33.9877, 42.0106, 43.9898,
    68.0626, 79.9568, 79.9663
  };
  return bio_ms_change;
}

vector<double> static getAdc_ms_change (){
  const vector<double> adc_ms_change = {
    21.9819, 37.9559, 46.0055, 67.9874, 83.9614, 
    60.0211, 82.0031, 97.9770, 97.9674, 62.9956, 
    84.9776, 62.0004, 83.9823, 99.9562, 97.9769, 
    43.9639, 75.9118, 27.0109, 32.0262, 39.9925, 
    55.9664, 41.0265
  };
  
  return adc_ms_change;
}

vector<double> static getFrg_ms_change (){
  const vector<double> frg_ms_change = {
    43.9898, 30.0106, 18.0106, 17.0265, 1.0078
  };
  return frg_ms_change;
}

vector<double> static getAll_ms_change (){
  const vector<double> all_ms_change = {
    0.9840, 1.0316, 2.0157, 14.0157, 15.0109,
    15.9949, 17.0265, 18.0106, 27.9949, 28.0313, 
    30.0106, 31.9720, 33.9877, 42.0106, 43.9898,
    68.0626, 79.9568, 79.9663,
    21.9819, 37.9559, 46.0055, 67.9874, 83.9614, 
    60.0211, 82.0031, 97.9770, 97.9674, 62.9956, 
    84.9776, 62.0004, 83.9823, 99.9562, 97.9769, 
    43.9639, 75.9118, 27.0109, 32.0262, 39.9925, 
    55.9664, 41.0265,
    43.9898, 30.0106, 18.0106, 17.0265, 1.0078
  };
  return all_ms_change;
}

// change direction of all rules
vector<int> static getBio_dir_change(){
  const vector<int> bio_dir_change = { // 0, add + minus; 1, add only; -1, minus only
    0,0,0,0,0,
    0,0,0,0,0,
    1,0,0
  };
  return bio_dir_change;
}

vector<int> static getAdc_dir_change(){
  const vector<int> adc_dir_change = {
    1,1,1,1,1,
    1,1,1,1,1,
    1,1,1,1,1,
    1,1,1
  };
  return adc_dir_change;
}

vector<int> static getFrg_dir_change(){
  const vector<int> frg_dir_change = {
    -1, -1, -1, -1, -1
  };
  return frg_dir_change;
}

vector<int> static getAll_dir_change(){
  const vector<int> all_dir_change = {
    0,0,0,0,0,
    0,0,0,0,0,
    1,0,0,
    1,1,1,1,1,
    1,1,1,1,1,
    1,1,1,1,1,
    1,1,1,
    -1, -1, -1, -1, -1
  };
  return all_dir_change;
}

#endif