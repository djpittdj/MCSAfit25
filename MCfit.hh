#include <boost/algorithm/string.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/foreach.hpp>
#include <boost/program_options.hpp>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include "Vector.hh"

const double COULOMB = 332.0636;
const double BOLTZMANN = 0.0019872041;
const int    QMconfig_num = 1331;
const int    HSE_num = 12;
const int    HSP_num = 13;
string errChangeGen("Total charge generation step reached, but still no hit. Terminate");
const Vector dx(-0.2, 0., 0.), dy(0., 0.2, 0.), dz(0., 0., 0.2);

double SumVector(vector<double> vec);

class MCfit {
private:
  double switchdist;
  double cutoff;
  double initTemperature, finalTemperature;
  double temperatureScalingFac;
  double chargePrecision;
  int maxChargeGenStep;
  int cycleMCStep;
  int chargeOutputFreq;
  vector<Vector> r0_HSE, r0_HSP;
  map<string, double> dict_epsilon, dict_sigma;
  vector<string> hse_name, hsp_name;
  map<string, string> hse_name2type, hsp_name2type;
  string fname_RMSD;
  vector<double> EQM, Evdw;
  vector<double> q0, q;
  // Internal variable to keep track of RMSD
  double RMSD;
  
public:
  MCfit(double i_switchdist, double i_cutoff, double i_initTemperature,
	double i_finalTemperature, double i_tScalingFac,
        double i_chargePrecision, int i_maxChargeGenStep, int i_cycleMCStep,
	int i_chargeOutputFreq,
	vector<Vector> i_r0_HSE, vector<Vector> i_r0_HSP,
	map<string, double> i_dict_epsilon, map<string, double> i_dict_sigma,
	vector<string> i_hse_name, vector<string> i_hsp_name,
	map<string, string> i_hse_name2type, map<string, string> i_hsp_name2type,
	string i_fname_RMSD,
	vector<double> i_EQM, vector<double> i_Evdw,
	vector<double> i_q0);
  double ComputeVDW(double epsilon_i, double sigma_i, Vector r_i,
                     double epsilon_j, double sigma_j, Vector r_j);
  double ComputeElec(double q_i, Vector r_i, double q_j, Vector r_j);
  // Debug output
  void Show_variables();
  void Show_r0_HSE();
  void Show_r0_HSP();
  void Show_dict_epsilon();
  void Show_dict_sigma();
  void Show_hse_name();
  void Show_hsp_name();
  void Show_hse_name2type();
  void Show_hsp_name2type();
  void Show_EQM();
  void Show_Evdw();
  void Show_q0();
  void Show_q();
  void Show_q(vector<double> q);
  void DebugOutput();
  void GenConstrChargeSet();
  double ComputeRMSD();
  double ComputeRMSD(vector<double> qq);
  void MC();
};
