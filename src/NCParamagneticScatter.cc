#include "NCParamagneticScatter.hh"

//Include various utilities from NCrystal's internal header files:
#include "NCrystal/internal/NCString.hh"
#include "NCrystal/internal/NCRandUtils.hh"
#include "NCrystal/internal/NCMath.hh"

//Components to calculate the magnetic scattering cross sections
double gfunc( double temperature, double D_const, int mag_scat )
{
  //g functions, which accounts for the spin
  //temperature : temperature of material K
  //D_const : zero-field splitting constant, eV
  //mag_scat : magnetic scattering option, int
  //-1, 0, 1 represent respectively down, elastic and up scattering
  //g0=g+
  nc_assert( mag_scat==-1 || mag_scat==0 || mag_scat==1 );
  double exp = NCrystal::exp_negarg_approx(-D_const / (NCrystal::constant_boltzmann * temperature));
  double g = 4. / 3. / (1 + 2 * exp);
  if ( mag_scat == 0 || mag_scat == 1 ) g *= exp;
  return g;
}

double ffunc( double incident_neutron_E, double hwhm, double D_const,
              int mag_scat, double msd )
{
  //f functions, integral of magnetic form factor
  //incident_neutron_E : incident neutron energy, eV
  //hwhm : half width at half maximum, float, Aa^-1
  //D_const : zero-field splitting constant, eV
  //mag_scat : magnetic scattering option, int
  //-1, 0, 1 represent respectively down, elastic and up scattering
  //msd: mean-squared displacement, Aa^2
  nc_assert( mag_scat==-1 || mag_scat==0 || mag_scat==1 );
  double f;
  if ( mag_scat == -1 && incident_neutron_E <= D_const ) {
    f = 0.;
  }
  else {
    double A = 2 * (msd + std::log(2) / (hwhm * hwhm)) / NCrystal::const_hhm; // eV^-1
    double cm = std::sqrt(1 + mag_scat*D_const / incident_neutron_E);
    f  = NCrystal::exp_negarg_approx(-A * incident_neutron_E * (1 - cm) * (1 - cm));
    f -= NCrystal::exp_negarg_approx(-A * incident_neutron_E * (1 + cm) * (1 + cm));
    f *= 1. / (4 * cm * A * incident_neutron_E);
  }
  return f;
}

double kgffunc( double temperature, double incident_neutron_E, double hwhm,
                double D_const, int mag_scat, double msd )
{
  //inelastic and elastic magnetic cross sections
  //The factor 0.5 is used to compute the cross section
  //per atom instead of per paramagnetic center
  //temperature : temperature of material K
  //incident_neutron_E : incident neutron energy, eV
  //hwhm : half width at half maximum, float, Aa^-1
  //D_const : zero-field splitting constant, eV
  //mag_scat : magnetic scattering option, int
  //-1, 0, 1 represent respectively down, elastic and up scattering
  //msd: mean-squared displacement, Aa^2
  double g = gfunc(temperature, D_const, mag_scat);
  double f = ffunc(incident_neutron_E, hwhm, D_const, mag_scat, msd);
  double kgf;
  if ( mag_scat == -1 && incident_neutron_E <= D_const ) {
    kgf = 0.;
  }
  else {
    kgf = g * f;
    kgf *= std::sqrt(1 + mag_scat*D_const / incident_neutron_E);
    kgf *= 0.5;
  }
  return kgf;
}

struct D_constContribution{
	
	double splittingConstant;
	double contribution;
	
};

std::vector<D_constContribution> calcD_constContribution( double temperature,
                                                          double incident_neutron_E,
                                                          double hwhm, 
                                                          double D_const_loc,
                                                          double hwhm_d,
                                                          double eta,
                                                          int mag_scat, double msd )
{
  //magnetic cross section for zero-field splitting
  //which follows a pseudo-Voigt distribution
  //temperature : temperature of material K
  //incident_neutron_E : incident neutron energy, eV
  //hwhm : half width at half maximum, float, Aa^-1
  //D_const_loc : central value of D, eV
  //hwhm_d : phenomenological hwhm of D, eV
  //eta : constant between 0 and 1, 1 represents pure Lorentzian
  //0 represents pure Gaussian
  //mag_scat : magnetic scattering option, int
  //-1, 0, 1 represent respectively down, elastic and up scattering
  //msd: mean-squared displacement, Aa^2
  //To do: improvement of bins of Lorentzian distribution
  //       for accelerating the calculation
  std::vector<D_constContribution> result;

  double D_lower_bound = NCrystal::ncmax( 0.01e-3, D_const_loc - 4 * hwhm_d );
  double D_upper_bound = NCrystal::ncmin( 1, D_const_loc + 4 * hwhm_d );
  unsigned int num_D = 100; //these values can be changed later
  double Sigma = hwhm_d / std::sqrt(2 * std::log(2));
  
  double D_const_weight_norm = 0.; //for normalizing the weight of D_const
  for ( auto D_const : NC::linspace( D_lower_bound, D_upper_bound, num_D ) ) {
    double G = 1. / Sigma / std::sqrt(2 * NCrystal::kPi) * NCrystal::exp_negarg_approx(-0.5 * (D_const - D_const_loc) * (D_const - D_const_loc) / Sigma / Sigma);
    double L = NCrystal::kInvPi * hwhm_d / ((D_const - D_const_loc) * (D_const - D_const_loc) + hwhm_d * hwhm_d);
    D_const_weight_norm += (1 - eta) * G + eta * L;
  }
  
  for ( auto D_const : NC::linspace( D_lower_bound, D_upper_bound, num_D ) ) {
		//D_const_weight obtained from the Gaussian probability distribution function
    double G = 1. / Sigma / std::sqrt(2 * NCrystal::kPi) * NCrystal::exp_negarg_approx(-0.5 * (D_const - D_const_loc) * (D_const - D_const_loc) / Sigma / Sigma);
    double L = NCrystal::kInvPi * hwhm_d / ((D_const - D_const_loc) * (D_const - D_const_loc) + hwhm_d * hwhm_d);
    double D_const_weight = (1 - eta) * G + eta * L;
    D_const_weight /= D_const_weight_norm; //normalization of weight
	
    double xs; //dimensionless cross section, before multiplying sigma_m
    if ( mag_scat==2 ) {
      xs  = kgffunc( temperature, incident_neutron_E, hwhm, D_const, -1, msd );
      xs += kgffunc( temperature, incident_neutron_E, hwhm, D_const,  0, msd );
      xs += kgffunc( temperature, incident_neutron_E, hwhm, D_const,  1, msd );
    }
    else xs = kgffunc( temperature, incident_neutron_E, hwhm, D_const, mag_scat, msd );
  
	  result.emplace_back();
	  result.back().splittingConstant = D_const;
	  result.back().contribution = D_const_weight * xs;
  }
	return result;
}

double mupdf( NC::RNG& rng, double incident_neutron_E, double hwhm,
              double D_const, int mag_scat, double msd )
{
  //angular probability distribution functions
  //incident_neutron_E : incident neutron energy, eV
  //hwhm : half width at half maximum, float, Aa^-1
  //D_const : zero-field splitting constant, eV
  //mag_scat : magnetic scattering option, int
  //-1, 0, 1 represent respectively down, elastic and up scattering
  //msd: mean-squared displacement, Aa^2
  nc_assert( mag_scat==-1 || mag_scat==0 || mag_scat==1 );
  double A = 2 * (msd + std::log(2) / (hwhm * hwhm)) / NCrystal::const_hhm; // eV^-1
  if ( mag_scat == -1 && incident_neutron_E <= D_const ) {
    return 1.0;
  }
  else {
    double B = 2 * A * std::sqrt(incident_neutron_E * (incident_neutron_E + mag_scat*D_const));
    //treatment analog to incoherent elastic scattering
    if ( B < 0.01 ) {
      //Rejection method:
      
      double maxval = NCrystal::exp_smallarg_approx(B);
      while (true) {
        double mu = rng.generate() * 2.0 - 1.0;
        if ( rng.generate() * maxval < NCrystal::exp_smallarg_approx(B * mu) )
          return mu;
      }
      
    } else {
      //Transformation method:
      
      // If f(x)=N*exp(a*x) is a normalised distribution on [-1,1], then
      // N=a/(exp(a)-exp(-a)) and the commulative probability function is F(x)=(
      // exp(a*(x+1)) -1 ) / ( exp(2*a) -1 ). With R a uniformly distributed
      // random number in (0,1], solving R=F(x) yields:
      //
      // x(R) = log( 1 + R * ( exp(2*a)-1 ) ) / a - 1
      //
      // Which can preferably be evaluated with expm1/log1p functions.
      return NCrystal::ncclamp(std::log1p( rng.generate() * std::expm1(2.0 * B) ) / B - 1.0, -1.0, 1.0);
    }
  }
}

double Eppdf( double incident_neutron_E, double D_const, int mag_scat )
{
  //energy probability distribution functions
  //incident_neutron_E : incident neutron energy, eV
  //D_const : zero-field splitting constant, eV
  //mag_scat : magnetic scattering option, int
  //-1, 0, 1 represent respectively down, elastic and up scattering
  nc_assert( mag_scat==-1 || mag_scat==0 || mag_scat==1 );
  double Ep;
  if ( mag_scat == -1 && incident_neutron_E <= D_const ) {
    Ep = incident_neutron_E;
  }
  else {
    Ep = incident_neutron_E + mag_scat*D_const;
  }
  return Ep;
}

bool NCP::ParamagneticScatter::isApplicable( const NC::Info& info )
{
  //Accept if input is NCMAT data with @CUSTOM_<pluginname> section:
  return info.countCustomSections(pluginNameUpperCase()) > 0;
}

NCP::ParamagneticScatter NCP::ParamagneticScatter::createFromInfo( const NC::Info& info )
{
  //Parse the content of our custom section. In case of syntax errors, we should
  //raise BadInput exceptions, to make sure users gets understandable error
  //messages. We should try to avoid other types of exceptions.

  //Get the relevant custom section data (and verify that there are not multiple
  //such sections in the input data):
  if ( info.countCustomSections( pluginNameUpperCase() ) != 1 )
    NCRYSTAL_THROW2(BadInput,"Multiple @CUSTOM_"<<pluginNameUpperCase()<<" sections are not allowed");
  auto data = info.getCustomSection( pluginNameUpperCase() );

  // data is here a vector of lines, and each line is a vector of words. In our
  // case, we want to accept sections of the form (units are barn, angstrom^-1,
  // eV, boolean, as is usual in NCrystal):
  //
  // @CUSTOM_<ourpluginname>
  //    <sigmavalue> <hwhm value> <D_const value> <magnetic scattering option>
  //

  //Verify we have exactly one line and four or five words:
  if ( data.size() != 1 || ( data.at(0).size()!=5
       && data.at(0).size()!=6 ) )
    NCRYSTAL_THROW2(BadInput,"Data in the @CUSTOM_"<<pluginNameUpperCase()
                    <<" section should be four or five numbers on a single line");

  //Parse and validate values:
  double sigma, hwhm, D_const_loc, hwhm_d, eta;
  if ( ! NC::safe_str2dbl( data.at(0).at(0), sigma )
       || ! NC::safe_str2dbl( data.at(0).at(1), hwhm )
       || ! NC::safe_str2dbl( data.at(0).at(2), D_const_loc )
       || ! NC::safe_str2dbl( data.at(0).at(3), hwhm_d )
       || ! NC::safe_str2dbl( data.at(0).at(4), eta )
       || ! (sigma>0.0) || ! (hwhm>0.0) || ! (D_const_loc>0.0) || ! (hwhm_d>=0.0) || ! (eta>=0.0) || ! (eta<=1.0) )
    NCRYSTAL_THROW2( BadInput,"Invalid values specified in the @CUSTOM_"<<pluginNameUpperCase()
                     <<" sigma, hwhm, D_const_loc (should be three strictly positive floating point values), hwhm_d should be positive, eta should be between 0. and 1." );
    
  int mag_scat;
  // inelastic and elastic magnetic cross sections are both considered
  if ( data.at(0).size()==5 ) mag_scat = 2;
  // only one component is considered
  else if ( ! NC::safe_str2int( data.at(0).at(5), mag_scat )
            || ! (mag_scat > -2) || ! (mag_scat < 2) )
    NCRYSTAL_THROW2( BadInput,"Invalid value specified in the @CUSTOM_"<<pluginNameUpperCase()
                     <<" mag_scat (must be -1 (for down-scattering) or 1 (for up-scattering) or 0 (for elastic scattering)" );
  
  //Getting the temperature
  double temperature = info.getTemperature().get();
    
  //Getting the mean-squared displacement (MSD)
  double msd = 0.;
  if ( info.hasAtomMSD() ) msd = info.getAtomInfos().front().msd().value();

  //Parsing done! Create and return our model:
  return ParamagneticScatter(sigma,hwhm,D_const_loc,temperature,mag_scat,msd,hwhm_d,eta);
}

NCP::ParamagneticScatter::ParamagneticScatter( double sigma, double hwhm, double D_const_loc,
                                               double temperature, int mag_scat, double msd,
                                               double hwhm_d, double eta )
  : m_sigma(sigma),
    m_hwhm(hwhm),
    m_D_const_loc(D_const_loc),
    m_temperature(temperature),
    m_mag_scat(mag_scat),
    m_msd(msd),
    m_hwhm_d(hwhm_d),
    m_eta(eta)
{
  //Important note to developers who are using the infrastructure in the
  //testcode/ subdirectory: If you change the number or types of the arguments
  //for the constructor here, you should make sure to perform a corresponding
  //change in three files in the testcode/ directory: _cbindings.py,
  //__init__.py, and NCForPython.cc - that way you can still instantiate your
  //model directly from your python test code).

  nc_assert( m_sigma > 0.0 );
  nc_assert( m_hwhm > 0.0 );
  nc_assert( m_D_const_loc > 0.0 );
  nc_assert( m_temperature > 0.0);
  nc_assert( m_hwhm_d >= 0.0 );
  nc_assert( m_eta >= 0.0 );
  nc_assert( m_eta <= 1.0 );
}

double NCP::ParamagneticScatter::calcCrossSection( double neutron_ekin ) const
{
  double xs;
  
  if ( m_hwhm_d < 1.e-9 ) {
    if ( m_mag_scat==2 ) {
      xs  = kgffunc( m_temperature, neutron_ekin, m_hwhm, m_D_const_loc, -1, m_msd );
      xs += kgffunc( m_temperature, neutron_ekin, m_hwhm, m_D_const_loc,  0, m_msd );
      xs += kgffunc( m_temperature, neutron_ekin, m_hwhm, m_D_const_loc,  1, m_msd );
    }
    else xs = kgffunc( m_temperature, neutron_ekin, m_hwhm, m_D_const_loc, m_mag_scat, m_msd );
  }
  else {
    auto contribs = calcD_constContribution( m_temperature, neutron_ekin, m_hwhm,
                                             m_D_const_loc, m_hwhm_d, m_eta,
                                             m_mag_scat, m_msd );
    NCrystal::StableSum sum;
    for( auto &e:contribs ) {
      sum.add( e.contribution );
    }
    xs = sum.sum();
  }
  return m_sigma * xs;
}

NCP::ParamagneticScatter::ScatEvent NCP::ParamagneticScatter::sampleScatteringEvent( NC::RNG& rng, double neutron_ekin ) const
{
  ScatEvent result;
  double D_const;

  if ( m_hwhm_d < 1.e-9 ) {
    D_const = m_D_const_loc;
  }
  else {
    auto contribs = calcD_constContribution( m_temperature, neutron_ekin, m_hwhm,
                                             m_D_const_loc, m_hwhm_d, m_eta,
                                             m_mag_scat, m_msd );
  
    std::vector<double> v;
    v.reserve( contribs.size() );
  
    NCrystal::StableSum sum;
    for( auto &e:contribs ) {
      sum.add( e.contribution );
      v.push_back( sum.sum() );
    }
  
    auto idx = pickRandIdxByWeight( rng, v );
    D_const = contribs.at( idx ).splittingConstant;
  }

  if ( m_mag_scat == 2 ) {
        
    double rand = rng.generate(); //random number
    double kgf_down = kgffunc( m_temperature, neutron_ekin, m_hwhm, D_const, -1, m_msd );
    double kgf_el   = kgffunc( m_temperature, neutron_ekin, m_hwhm, D_const,  0, m_msd );
    double kgf_up   = kgffunc( m_temperature, neutron_ekin, m_hwhm, D_const,  1, m_msd );
    double kgf_tot  = kgf_down + kgf_el + kgf_up;
    
    if ( rand < kgf_down / kgf_tot ) {
      //down-scattering happens
      if ( neutron_ekin <= D_const ) {
        result.ekin_final = neutron_ekin;
        result.mu = 1.;
      }
      else {
        result.ekin_final = Eppdf( neutron_ekin, D_const, -1 );
        result.mu         = mupdf( rng, neutron_ekin, m_hwhm, D_const, -1, m_msd );
      }
    }
    else if ( rand < (kgf_down + kgf_up) / kgf_tot ) {
      //up-scattering happens
      result.ekin_final = Eppdf( neutron_ekin, D_const, 1 );
      result.mu         = mupdf( rng, neutron_ekin, m_hwhm, D_const, 1, m_msd );
    }
    else {
      //elastic scattering happens
      result.ekin_final = Eppdf( neutron_ekin, D_const, 0 );
      result.mu         = mupdf( rng, neutron_ekin, m_hwhm, D_const, 0, m_msd );
    }
  }
      
  else if ( m_mag_scat == -1 ) {
    if ( neutron_ekin <= D_const ) {
      result.ekin_final = neutron_ekin;
      result.mu = 1.;
    }
    else {
      result.ekin_final = Eppdf( neutron_ekin, D_const, -1 );
      result.mu         = mupdf( rng, neutron_ekin, m_hwhm, D_const, -1, m_msd );
    }
  }
    
  else if ( m_mag_scat == 1 ) {
    result.ekin_final = Eppdf( neutron_ekin, D_const, 1 );
    result.mu         = mupdf( rng, neutron_ekin, m_hwhm, D_const, 1, m_msd );
  }
    
  else {
    result.ekin_final = Eppdf( neutron_ekin, D_const, 0 );
    result.mu         = mupdf( rng, neutron_ekin, m_hwhm, D_const, 0, m_msd );
  }

  return result;
}
