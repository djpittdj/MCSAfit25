#include "MCfit.hh"
using namespace std;
using namespace boost;
namespace po = boost::program_options;
typedef vector<string> split_vector_type;
typedef boost::mt19937 base_generator_type;

// A helper function to simplify the main part.
template<class T>
ostream& operator<<(ostream& os, const vector<T>& v) {
  copy(v.begin(), v.end(), ostream_iterator<T>(cout, " "));
  return os;
}

// Sum components of a vector
double SumVector(vector<double> vec) {
  double sum=0.;
  for (vector<double>::iterator iter=vec.begin(); iter!=vec.end(); ++iter)
    sum+=(*iter);

  return sum;
}

// Initializer:
// Set values, read related files and store them in vectors.
MCfit::MCfit(double i_switchdist, double i_cutoff, double i_initTemperature,
	     double i_finalTemperature, double i_tScalingFac,
             double i_chargePrecision, int i_maxChargeGenStep, int i_cycleMCStep,
	     int i_chargeOutputFreq,
	     vector<Vector> i_r0_HSE, vector<Vector> i_r0_HSP,
	     map<string, double> i_dict_epsilon, map<string, double> i_dict_sigma,
	     vector<string> i_hse_name, vector<string> i_hsp_name,
	     map<string, string> i_hse_name2type, map<string, string> i_hsp_name2type,
	     string i_fname_RMSD,
	     vector<double> i_EQM, vector<double> i_Evdw,
	     vector<double> i_q0) {
  switchdist = i_switchdist;
  cutoff = i_cutoff;
  initTemperature = i_initTemperature;
  finalTemperature = i_finalTemperature;
  temperatureScalingFac = i_tScalingFac;
  chargePrecision = i_chargePrecision;
  maxChargeGenStep = i_maxChargeGenStep;
  cycleMCStep = i_cycleMCStep;
  chargeOutputFreq = i_chargeOutputFreq;
  fname_RMSD = i_fname_RMSD;
  
  r0_HSE.reserve(HSE_num*sizeof(Vector));
  r0_HSE = i_r0_HSE;

  r0_HSP.reserve(HSP_num*sizeof(Vector));
  r0_HSP = i_r0_HSP;

  dict_epsilon = i_dict_epsilon;
  dict_sigma = i_dict_sigma;

  hse_name.reserve(HSE_num*sizeof(string));
  hse_name = i_hse_name;
  hsp_name.reserve(HSP_num*sizeof(string));
  hsp_name = i_hsp_name;

  hse_name2type = i_hse_name2type;
  hsp_name2type = i_hsp_name2type;

  EQM.reserve(QMconfig_num*sizeof(double));
  EQM = i_EQM;

  Evdw.reserve(QMconfig_num*sizeof(double));
  Evdw = i_Evdw;

  q0.reserve((HSE_num+HSP_num)*sizeof(double));
  q0 = i_q0;

  q.reserve((HSE_num+HSP_num)*sizeof(double));

  // Initialize q with q0
  q = q0;
}

// The code is adapted from mindy/mdenergy source code
// But actually in this version, the vdW force is not calulated but
// imported from the vdW file calculated separately from a Python
// program getMMvdw25.py
double MCfit::ComputeVDW(double epsilon_i, double sigma_i, Vector r_i,
                         double epsilon_j, double sigma_j, Vector r_j) {
  double switch2 = pow(switchdist, 2.);
  double cut2 = pow(cutoff, 2.);

  double c1 = 1./(cut2 - switch2);
  c1 = pow(c1, 3.);
  double c3 = 4.*c1;

  Vector vr_ij = r_i - r_j;
  double dist2 = vr_ij.length2();

  // If distance is larger than the cutoff, don't calculate and return 0.
  if (dist2 > cut2) return 0.;

  double r = sqrt(dist2);
  double r_1 = 1./r;
  double r_2 = pow(r_1, 2.);
  double r_6 = pow(r_1, 6.);
  double r_12 = pow(r_1, 12.);
  double switchVal = 1.;
  if (dist2 > switch2) {
    double c2 = cut2 - dist2;
    double c4 = c2*(cut2 + 2*dist2 - 3.*switch2);
    switchVal = c2*c4*c1;
  }

  epsilon_i = -1.*epsilon_i;
  sigma_i = pow(2., 5./6.)*sigma_i;
  epsilon_j = -1.*epsilon_j;
  sigma_j = pow(2., 5./6.)*sigma_j;
  double epsilon_ij = sqrt(epsilon_i*epsilon_j);
  double sigma_ij = 0.5*(sigma_i+sigma_j);
  sigma_ij = pow(sigma_ij, 6.);

  double B = 4.*sigma_ij*epsilon_ij;
  double A = B*sigma_ij;
  double AmBterm = (A*r_6-B)*r_6;
  double vdwenergy = switchVal*AmBterm;

  return vdwenergy;
}

double MCfit::ComputeElec(double q_i, Vector r_i,
                          double q_j, Vector r_j) {
  double cut2 = pow(cutoff, 2.);
  double kqq = COULOMB*q_i*q_j;
  Vector vr_ij = r_i - r_j;
  double dist2 = vr_ij.length2();
  double r_ij = sqrt(dist2);
  double r_1 = 1./r_ij;
  double efac = 1. - dist2/cut2;
  double prefac = kqq * r_1 * efac;
  double elecenergy = prefac * efac;

  return elecenergy;
}

void MCfit::Show_variables() {
  cout<<"#Show variable"<<endl;
  cout<<"#switchdist = "<<switchdist<<endl;
  cout<<"#cutoff = "<<cutoff<<endl;
  cout<<"#initTemperature = "<<initTemperature<<endl;
  cout<<"#finalTemperature = "<<finalTemperature<<endl;
  cout<<"#temperatureScalingFac = "<<temperatureScalingFac<<endl;
  cout<<"#chargePrecision = "<<chargePrecision<<endl;
  cout<<"#maxChargeGenStep = "<<maxChargeGenStep<<endl;
  cout<<"#cycleMCStep = "<<cycleMCStep<<endl;
  cout<<"#chargeOutputFreq = "<<chargeOutputFreq<<endl;
  cout<<"#fname_RMSD "<<fname_RMSD<<endl;
}

void MCfit::Show_r0_HSE() {
  cout<<"Show r0_HSE"<<endl;
  BOOST_FOREACH(Vector xyz, r0_HSE) cout<<xyz<<" ";
  cout<<endl;
}

void MCfit::Show_r0_HSP() {
  cout<<"Show r0_HSP"<<endl;
  BOOST_FOREACH(Vector xyz, r0_HSP) cout<<xyz<<" ";
  cout<<endl;
}

void MCfit::Show_dict_epsilon() {
  cout<<"Show dict_epsilon"<<endl;
  for (map<string, double>::iterator iter=dict_epsilon.begin();
       iter!=dict_epsilon.end();
       ++iter) {
    cout<<(*iter).first<<"\t"<<(*iter).second<<endl;
  }
}

void MCfit::Show_dict_sigma() {
  cout<<"Show dict_sigma"<<endl;
  for (map<string, double>::iterator iter=dict_sigma.begin();
       iter!=dict_sigma.end();
       ++iter) {
    cout<<(*iter).first<<"\t"<<(*iter).second<<endl;
  }
}

void MCfit::Show_hse_name() {
  cout<<"Show hse_name:"<<endl;
  BOOST_FOREACH (string name, hse_name) cout<<name<<endl;
  cout<<endl;
}

void MCfit::Show_hsp_name() {
  cout<<"Show hsp_name:"<<endl;
  BOOST_FOREACH (string name, hsp_name) cout<<name<<endl;
  cout<<endl;
}

void MCfit::Show_hse_name2type() {
  cout<<"Show hse_name2type:"<<endl;
  for (map<string, string>::iterator iter=hse_name2type.begin();
       iter!=hse_name2type.end();
       ++iter) {
    cout<<(*iter).first<<"\t"<<(*iter).second<<endl;
  }
  cout<<endl;
}

void MCfit::Show_hsp_name2type() {
  cout<<"Show hsp_name2type:"<<endl;
  for (map<string, string>::iterator iter=hsp_name2type.begin();
       iter!=hsp_name2type.end();
       ++iter) {
    cout<<(*iter).first<<"\t"<<(*iter).second<<endl;
  }
  cout<<endl;
}

void MCfit::Show_EQM() {
  cout<<"Show EQM:"<<endl;
  BOOST_FOREACH (double eqm, EQM) cout<<eqm<<endl;
}

void MCfit::Show_Evdw() {
  cout<<"Show Evdw:"<<endl;
  BOOST_FOREACH (double evdw, Evdw) cout<<evdw<<endl;
}

void MCfit::Show_q0() {
  cout<<"#Show q0:"<<endl;
  BOOST_FOREACH(double q_i, q0) cout<<q_i<<" ";
  cout<<endl;
}

void MCfit::Show_q() {
  cout<<"Show q:"<<endl;
  BOOST_FOREACH(double q_i, q) cout<<q_i<<" ";
  cout<<endl;
}

void MCfit::Show_q(vector<double> q) {
  cout<<"Show q:"<<endl;
  BOOST_FOREACH(double q_i, q) cout<<q_i<<" ";
  cout<<endl;
}

void MCfit::DebugOutput() {
  Show_variables(); Show_r0_HSE(); Show_r0_HSP();
  Show_dict_epsilon(); Show_dict_sigma(); Show_hse_name(); Show_hsp_name();
  Show_hse_name2type(); Show_hsp_name2type();
  Show_EQM(); Show_Evdw(); Show_q0(); Show_q();
}

void MCfit::GenConstrChargeSet() {
  base_generator_type generator(42);
  generator.seed(static_cast<unsigned int>(std::time(0)));
  //variate_generator combines the generator and the distribution to make a random number generator
  boost::uniform_real<> uni_dist(0, 1);
  boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni_real(generator, uni_dist);

  boost::uniform_int<> uni_dist40(-40, 40);
  boost::variate_generator<base_generator_type&, boost::uniform_int<> > uni_int40(generator, uni_dist40);

  boost::uniform_int<> uni_dist100(-100, 100);
  boost::variate_generator<base_generator_type&, boost::uniform_int<> > uni_int100(generator, uni_dist100);

  double q_1, q_2, q_3, q_4, q_5, q_6, q_7, q_8,
          q_9, q_10, q_11, q_12,
          q_13, q_14, q_15, q_16, q_17, q_18, q_19,
	  q_20, q_21, q_22, q_23, q_24, q_25;

  vector<double> q_HSE2_8(7), q_HSE9_12(4), q_HSP14_19(6), q_HSP20_25(6);
  // q0: 0.09, -0.08, 0.09, 0.09, -0.7, 0.22, 0.25, 0.13,
  //     -0.36, 0.32, -0.05, 0.09,
  //     0.09, -0.05, 0.09, 0.09, 0.19, 0.13, 0.19,
  //     -0.51, 0.44, -0.51, 0.44, 0.32, 0.18

  for (unsigned chargeGenStep=0; chargeGenStep<maxChargeGenStep; ++chargeGenStep) {
    q_2 = q0.at(1) + chargePrecision*uni_int40();  q_HSE2_8.at(0)=q_2;
    q_3 = q0.at(2) + chargePrecision*uni_int40();  q_HSE2_8.at(1)=q_3;
    q_1 = q_3;
    q_4 = q_3;                                     q_HSE2_8.at(2)=q_4;
    q_5 = q0.at(4) + chargePrecision*uni_int100(); q_HSE2_8.at(3)=q_5;
    q_6 = q0.at(5) + chargePrecision*uni_int100(); q_HSE2_8.at(4)=q_6;
    q_7 = q0.at(6) + chargePrecision*uni_int100(); q_HSE2_8.at(5)=q_7;
    q_8 = q0.at(7) + chargePrecision*uni_int100(); q_HSE2_8.at(6)=q_8;
    if (chargeGenStep == maxChargeGenStep) {
      cout<<errChangeGen<<endl;
      exit(1);
    }

    if (SumVector(q_HSE2_8) == 0.) {
      break;
    }
  }

  for (unsigned chargeGenStep=0; chargeGenStep<maxChargeGenStep; ++chargeGenStep) {
    q_9 = q0.at(8) + chargePrecision*uni_int100();  q_HSE9_12.at(0)=q_9;
    q_10 = q0.at(9) + chargePrecision*uni_int100(); q_HSE9_12.at(1)=q_10;
    q_11 = q0.at(10) + chargePrecision*uni_int40(); q_HSE9_12.at(2)=q_11;
    q_12 = q0.at(11) + chargePrecision*uni_int40(); q_HSE9_12.at(3)=q_12;
    if (chargeGenStep == maxChargeGenStep) {
      cout<<errChangeGen<<endl;
      exit(1);
    }

    if (SumVector(q_HSE9_12) == 0.) {
      break;
    }
  }

  for (unsigned chargeGenStep=0; chargeGenStep<maxChargeGenStep; ++chargeGenStep) {
    q_14 = q0.at(13) + chargePrecision*uni_int40();  q_HSP14_19.at(0)=q_14;
    q_15 = q0.at(14) + chargePrecision*uni_int40();  q_HSP14_19.at(1)=q_15;
    q_13 = q_15;
    q_16 = q_15;                                     q_HSP14_19.at(2)=q_16;
    q_17 = q0.at(16) + chargePrecision*uni_int100(); q_HSP14_19.at(3)=q_17;
    q_18 = q0.at(17) + chargePrecision*uni_int100(); q_HSP14_19.at(4)=q_18;
    q_19 = q0.at(18) + chargePrecision*uni_int100(); q_HSP14_19.at(5)=q_19;
    if (chargeGenStep == maxChargeGenStep) {
      cout<<errChangeGen<<endl;
      exit(1);
    }

    if (SumVector(q_HSP14_19) == 0.64) {
      break;
    }
  }

  for (unsigned chargeGenStep=0; chargeGenStep<maxChargeGenStep; ++chargeGenStep) {
    q_20 = q0.at(19) + chargePrecision*uni_int100(); q_HSP20_25.at(0)=q_20;
    q_21 = q0.at(20) + chargePrecision*uni_int100(); q_HSP20_25.at(1)=q_21;
    q_22 = q_20;                                     q_HSP20_25.at(2)=q_22;
    q_23 = q_21;                                     q_HSP20_25.at(3)=q_23;
    q_24 = q0.at(23) + chargePrecision*uni_int100(); q_HSP20_25.at(4)=q_24;
    q_25 = q0.at(24) + chargePrecision*uni_int100(); q_HSP20_25.at(5)=q_25;
    if (chargeGenStep == maxChargeGenStep) {
      cout<<errChangeGen<<endl;
      exit(1);
    }

    if (SumVector(q_HSP20_25) == 0.36) {
      break;
    }
  }

  double q_array[] = {q_1, q_2, q_3, q_4, q_5, q_6, q_7, q_8,
                       q_9, q_10, q_11, q_12,
		       q_13, q_14, q_15, q_16, q_17, q_18, q_19,
                       q_20, q_21, q_22, q_23, q_24, q_25};

  vector<double> tmp_q(q_array, q_array+sizeof(q_array)/sizeof(double));
  q = tmp_q;
}

double MCfit::ComputeRMSD() {
  unsigned configcounter = 0;
  double RMSD = 0.;
  for (int x=0; x<11; ++x) {
    for (int y=-5; y<6; ++y) {
      for (int z=-5; z<6; ++z) {
	double EMM = 0.;
	double Eelec = 0.;
	for (unsigned i=0; i<hse_name.size(); ++i) {
	  Vector r_i = r0_HSE.at(i);
	  for (unsigned j=0; j<hsp_name.size(); ++j) {
	    Vector r_j = r0_HSP.at(j) + x*dx + y*dy + z*dz;

	    double e_elec = ComputeElec(q.at(i), r_i, q.at(12+j), r_j);
	    Eelec += e_elec;
	  }
	}

	EMM = Evdw.at(configcounter) + Eelec;
	RMSD += pow((EQM.at(configcounter) - EMM), 2.);
	++configcounter;
      }
    }
  }
  RMSD = sqrt(RMSD/double(QMconfig_num));
  return RMSD;
}

double MCfit::ComputeRMSD(vector<double> qq) {
  unsigned configcounter = 0;
  double RMSD = 0.;
  for (int x=0; x<11; ++x) {
    for (int y=-5; y<6; ++y) {
      for (int z=-5; z<6; ++z) {
	double EMM = 0.;
	double Eelec = 0.;
	for (unsigned i=0; i<hse_name.size(); ++i) {
	  Vector r_i = r0_HSE.at(i);
	  for (unsigned j=0; j<hsp_name.size(); ++j) {
	    Vector r_j = r0_HSP.at(j) + x*dx + y*dy + z*dz;

	    double e_elec = ComputeElec(qq.at(i), r_i, qq.at(12+j), r_j);
	    Eelec += e_elec;
	  }
	}
	EMM = Evdw.at(configcounter) + Eelec;
	RMSD += pow((EQM.at(configcounter) - EMM), 2.);
	++configcounter;
      }
    }
  }
  RMSD = sqrt(RMSD/double(QMconfig_num));
  return RMSD;
}

void MCfit::MC() {
  base_generator_type generator(42);
  generator.seed(static_cast<unsigned int>(std::time(0)));
  //variate_generator combines the generator and the distribution to make a random number generator
  boost::uniform_real<> uni_dist(0, 1);
  boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni_real(generator, uni_dist);
  
  vector<double> q_old, q_best;
  q_old = q;
  q_best = q;
  cout<<"#Initial charge set:"<<endl;
  Show_q();
  
  double RMSD_old, RMSD_new, RMSD_best;
  RMSD = ComputeRMSD();
  RMSD_old = RMSD;
  RMSD_best = RMSD;
  cout<<"#Initial RMSD "<<RMSD_old<<endl;

  double temperature = initTemperature;

  ofstream rmsdfile(fname_RMSD.c_str());
  while (temperature > finalTemperature) {
    cout<<"Current temperature is :"<<temperature<<endl;
    for (unsigned MCstep=0; MCstep<cycleMCStep; ++MCstep) {
      RMSD_old = RMSD;
      q_old = q;
      
      // Generate a new charge set and calculate its RMSD with QM data
      GenConstrChargeSet();
      RMSD_new = ComputeRMSD();
      
      // Metropolis condition
      bool ACCEPT = false;
      if (RMSD_new <= RMSD_old) {
	ACCEPT = true;
      } else {
	double BoltzmannFac = exp(-1.*(RMSD_new - RMSD_old)/BOLTZMANN/temperature);
	double a = uni_real();
	if (a <= BoltzmannFac) {
          ACCEPT = true;
	} else {
          ACCEPT = false;
	}
      }

      if (ACCEPT == true) {
	RMSD = RMSD_new;
      } else {
	q = q_old;
	RMSD = RMSD_old;
      }
      
      if (RMSD < RMSD_best) {
        RMSD_best = RMSD;
	q_best = q;
      }
      
      // Output
      if (MCstep % chargeOutputFreq == 0) {
	cout<<"#"<<MCstep<<endl;
	cout<<"Best charge set so far\t";
	Show_q(q_best);
	cout<<"RMSD_best "<<RMSD_best<<endl;
        cout<<"New charge set\t";
	Show_q();
	cout<<"New RMSD: "<<RMSD_new<<endl;
	cout<<endl;
      }
      rmsdfile<<MCstep<<"\t"<<RMSD_old<<endl;
    }
    temperature *= temperatureScalingFac;
  }

  rmsdfile.close();
}

int main(int argc, char* argv[]) {
  time_t t_begin=time(NULL);

  try {
    // Name for the input config file
    // and for the output RMSD file
    string config_fname, fname_RMSD;
    // For reading lines from a file
    string line;
    double switchdist, cutoff, initTemperature, finalTemperature,
           temperatureScalingFac, chargePrecision;
    int maxChargeGenStep, cycleMCStep, chargeOutputFreq;

    // Declare a group of options that will be
    // allowed only on command line
    po::options_description generic("Generic options");
    generic.add_options()
      ("help", "produce help messages")
      ("config, c", po::value<string>(&config_fname)->default_value("MCSAfit25.config"),
                    "name of the input running parameter file.")
    ;

    // Declare a group of options that will be
    // allowed both on command line and in
    // config file
    po::options_description config("Running parameters");
    config.add_options()
      ("switchdist", po::value<double>(&switchdist)->default_value(15.),
		      "Switch distance<=cutoff")
      ("cutoff", po::value<double>(&cutoff)->default_value(17.),
		      "cutoff distance")
      ("initTemperature", po::value<double>(&initTemperature)->default_value(500.),
		      "SA initial temperature")
      ("finalTemperature", po::value<double>(&finalTemperature)->default_value(500.),
		      "SA final temperature")
      ("temperatureScalingFac", po::value<double>(&temperatureScalingFac)->default_value(0.9),
                      "Temperature scaling factor for SA")
      ("chargePrecision", po::value<double>(&chargePrecision)->default_value(0.001),
		      "Charge precision")
      ("maxChargeGenStep", po::value<int>(&maxChargeGenStep)->default_value(10000),
		      "Max number of steps to generate a constrained charge set")
      ("cycleMCStep", po::value<int>(&cycleMCStep)->default_value(100),
		      "Number of MC steps in each temperature cycle")
      ("chargeOutputFreq", po::value<int>(&chargeOutputFreq)->default_value(5),
                      "Charge output frequency")
      ("fname_RMSD", po::value<string>(&fname_RMSD)->default_value(string("RMSD25.dat")),
                      "The name for the output RMSD file")
      ("r0fname", po::value<string>()->default_value(string("s1_0_0_0.gjf")),
		      "Initial coordinate for HSE and HSP, from Gaussian input file")
      ("epsilon-sigma", po::value<string>()->default_value(string("par_sub.prm")),
		      "Read from parameter file to build dictionaries of epsilon and sigma")
      ("hse-rtf", po::value<string>()->default_value(string("HSE12_QM.rtf")),
		      "Read from HSE12.rtf to get atom name, type and charge")
      ("hsp-rtf", po::value<string>()->default_value(string("HSP13_QM.rtf")),
		      "Read from HSP13.rtf to get atom name, type and charge")
      ("fnameEQM", po::value<string>()->default_value(string("QM.dat")),
		      "Read EQM data")
      ("fnameEvdw", po::value<string>()->default_value(string("MMvdw25.dat")),
		      "Read Evdw data from a pre-calculated file")
      ("fnameq0", po::value<string>()->default_value(string("CHARMMq0.dat")),
		      "Read q0 data, the default value is from CHARMM")
    ;

    po::options_description cmdline_options;
    cmdline_options.add(generic).add(config);

    po::options_description config_file_options;
    config_file_options.add(config);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
	      options(cmdline_options).run(), vm);
    po::notify(vm);

    ifstream ifs(config_fname.c_str());
    if (!ifs) {
      cout<<"Can not open config file: "<<config_fname<<endl;
      cout<<"May not be in the current directory"<<endl;
      return 0;
    } else {
      store(parse_config_file(ifs, config_file_options), vm);
      notify(vm);
    }

    if (chargeOutputFreq > cycleMCStep) {
      cout<<"chargeOutputFreq should be less than cycleMCStep"<<endl;
      return 0;
    }

    if (vm.count("help")) {
      cout<<"Usage: MCfit25 [options]\n";
      cout<<config;
      return 0;
    }

    vector<Vector> r0_HSE, r0_HSP;

    if (vm.count("r0fname")) {
      string r0fname = vm["r0fname"].as<string>();
      // Build r0_HSE and r0_HSP
      vector<Vector> r0(25);
      ifstream r0_file(r0fname.c_str());
      unsigned i = 0;
      while (getline(r0_file, line)) {
	if (boost::starts_with(line, "H") or boost::starts_with(line, "N") or boost::starts_with(line, "C")) {
	  split_vector_type tokens;
	  split(tokens, line, is_any_of("\t "), token_compress_on);
	  Vector xyz(atof(tokens.at(1).c_str()),
		    atof(tokens.at(2).c_str()),
		    atof(tokens.at(3).c_str()));
	  r0.at(i++)=xyz;
	}
      }
      r0_file.close();

      for (int i=0; i<HSE_num; ++i) r0_HSE.push_back(r0.at(i));
      for (int i=0; i<HSP_num; ++i) r0_HSP.push_back(r0.at(i+HSE_num));
    }

    // Build dictionaries for epsilon and sigma
    map<string, double> dict_epsilon, dict_sigma;

    if (vm.count("epsilon-sigma")) {
      string dictfname = vm["epsilon-sigma"].as<string>();
      // Build dict_epsilon and dict_sigma
      ifstream par_sub_file(dictfname.c_str());
      while (getline(par_sub_file, line)) {
	if (!boost::starts_with(line, "!")) {
	  split_vector_type tokens;
	  split(tokens, line, is_any_of("\t "), token_compress_on);

	  dict_epsilon[tokens.at(0)] = atof(tokens.at(1).c_str());
	  dict_sigma[tokens.at(0)] = atof(tokens.at(2).c_str());
	}
      }
      par_sub_file.close();
    }

    vector<string> hse_name, hsp_name;
    map<string, string> hse_name2type, hsp_name2type;

    if (vm.count("hse-rtf")) {
      string hse_rtf_fname = vm["hse-rtf"].as<string>();
      ifstream hse_qm_file(hse_rtf_fname.c_str());
      while (getline(hse_qm_file, line)) {
	if (boost::starts_with(line, "ATOM")) {
	  split_vector_type tokens;
	  split(tokens, line, is_any_of(" "), token_compress_on);

	  hse_name.push_back(tokens.at(1));
	  hse_name2type[tokens.at(1)] = tokens.at(2);
	}
      }
      hse_qm_file.close();
    }

    if (vm.count("hsp-rtf")) {
      string hsp_rtf_fname = vm["hsp-rtf"].as<string>();
      ifstream hsp_qm_file(hsp_rtf_fname.c_str());
      while (getline(hsp_qm_file, line)) {
	if (boost::starts_with(line, "ATOM")) {
	  split_vector_type tokens;
	  split(tokens, line, is_any_of(" "), token_compress_on);

	  hsp_name.push_back(tokens.at(1));
	  hsp_name2type[tokens.at(1)] = tokens.at(2);
	}
      }
      hsp_qm_file.close();
    }

    vector<double> EQM, Evdw, q0;

    if (vm.count("fnameEQM")) {
      string fnameEQM = vm["fnameEQM"].as<string>();
      ifstream qm_file(fnameEQM.c_str());
      while (getline(qm_file, line)) {
	split_vector_type tokens;
	split(tokens, line, is_any_of(" \t"), token_compress_on);
	EQM.push_back(atof(tokens.at(3).c_str()));
      }
      qm_file.close();
    }

    if (vm.count("fnameEvdw")) {
      string fnameEvdw = vm["fnameEvdw"].as<string>();
      ifstream vdw_file(fnameEvdw.c_str());
      while (getline(vdw_file, line)) {
	split_vector_type tokens;
	split(tokens, line, is_any_of(" \t"), token_compress_on);
	Evdw.push_back(atof(tokens.at(3).c_str()));
      }
      vdw_file.close();
    }

    if (vm.count("fnameq0")) {
      string fnameq0 = vm["fnameq0"].as<string>();
      ifstream q0_file(fnameq0.c_str());
      while (getline(q0_file, line)) {
	split_vector_type tokens;
	split(tokens, line, is_any_of(" \t"), token_compress_on);
	q0.push_back(atof(tokens.at(0).c_str()));
      }
      q0_file.close();
    }

    MCfit mcfit(switchdist, cutoff, initTemperature,
		finalTemperature, temperatureScalingFac,
		chargePrecision, maxChargeGenStep, cycleMCStep,
		chargeOutputFreq,
		r0_HSE, r0_HSP,
		dict_epsilon, dict_sigma,
		hse_name, hsp_name,
		hse_name2type, hsp_name2type,
		fname_RMSD,
		EQM, Evdw,
		q0);
    mcfit.Show_variables();
    mcfit.MC();

    //double q_array[23] = {-0.065, 0.071, 0.071, -0.776, 0.226, 0.345, 0.128, -0.396, 0.333, -0.055, 0.118, -0.047, 0.057, 0.057, 0.226, 0.192, 0.155, -0.587, 0.476, -0.587, 0.476, 0.361, 0.221 };
    //vector<double> tmp_q(q_array, q_array+sizeof(q_array)/sizeof(double));
    //cout<<"RMSD "<<mcfit.ComputeRMSD(tmp_q)<<endl;

    cout<<"Program takes "<<time(NULL) - t_begin<<" seconds to finish."<<endl;
  }
  catch(std::exception& e)
  {
    cout << e.what() << "\n";
    return 1;
  }

  return 0;
}
