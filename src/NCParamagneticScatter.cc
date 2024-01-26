#include "NCParamagneticScatter.hh"

//Include various utilities from NCrystal's internal header files:
#include "NCrystal/internal/NCString.hh"
#include "NCrystal/internal/NCRandUtils.hh"

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

double hfunc( double incident_neutron_E, double D_const,
              int mag_scat, double hwhm_o )
{
  //h functions, integral of Lorentzian functions 
  //which account for reorientational motions
  //incident_neutron_E : incident neutron energy, eV
  //D_const : zero-field splitting constant, eV
  //mag_scat : magnetic scattering option, int
  //-1, 0, 1 represent respectively down, elastic and up scattering
  //hwhm_o : phenomenological hwhm of oxygen, eV
  nc_assert( mag_scat==-1 || mag_scat==0 || mag_scat==1 );
  double h;
  if ( hwhm_o < 1.e-9 ) {
    h = 1;
  }
  else {
    if ( mag_scat == -1 && incident_neutron_E <= D_const ) {
      h = 0;
    }
    else {
      h = NCrystal::atan_approx((incident_neutron_E + mag_scat*D_const) / hwhm_o);
      h *= NCrystal::kInvPi;
      h += 0.5;
    }
  }
  return h;
}

double kgfhfunc( double temperature, double incident_neutron_E, double hwhm,
                 double D_const, int mag_scat, double msd, double hwhm_o )
{
  //inelastic and elastic magnetic cross sections
  //temperature : temperature of material K
  //incident_neutron_E : incident neutron energy, eV
  //hwhm : half width at half maximum, float, Aa^-1
  //D_const : zero-field splitting constant, eV
  //mag_scat : magnetic scattering option, int
  //-1, 0, 1 represent respectively down, elastic and up scattering
  //msd: mean-squared displacement, Aa^2
  //hwhm_o : phenomenological hwhm of oxygen, eV
  double g = gfunc(temperature, D_const, mag_scat);
  double f = ffunc(incident_neutron_E, hwhm, D_const, mag_scat, msd);
  double h = hfunc(incident_neutron_E, D_const, mag_scat, hwhm_o );
  double kgfh;
  if ( mag_scat == -1 && incident_neutron_E <= D_const ) {
    kgfh = 0.;
  }
  else {
    kgfh = g * f * h;
    kgfh *= std::sqrt(1 + mag_scat*D_const / incident_neutron_E);
  }      
  return kgfh;
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

double Eppdf( NC::RNG& rng, double incident_neutron_E, double hwhm,
              double D_const, int mag_scat, double msd, double hwhm_o )
{
  //energy probability distribution functions
  //incident_neutron_E : incident neutron energy, eV
  //hwhm : half width at half maximum, float, Aa^-1
  //D_const : zero-field splitting constant, eV
  //mag_scat : magnetic scattering option, int
  //-1, 0, 1 represent respectively down, elastic and up scattering
  //msd: mean-squared displacement, Aa^2
  //hwhm_o : phenomenological hwhm of oxygen, eV
  nc_assert( mag_scat==-1 || mag_scat==0 || mag_scat==1 );
  double Ep;
  if ( mag_scat == -1 && incident_neutron_E <= D_const ) {
    Ep = incident_neutron_E;
  }
  else {
    if ( hwhm_o < 1.e-9 ) {
      Ep = incident_neutron_E + mag_scat*D_const;
    }
    else {
      double randEp = rng.generate();
      double h = hfunc(incident_neutron_E, D_const, mag_scat, hwhm_o );
      double atan = NCrystal::atan_approx((incident_neutron_E + mag_scat*D_const) / hwhm_o);
      Ep = incident_neutron_E + mag_scat*D_const - hwhm_o * std::tan(atan - NCrystal::kPi * randEp * h);
    }
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
  if ( data.size() != 1 || ( data.at(0).size()!=4
       && data.at(0).size()!=5 ) )
    NCRYSTAL_THROW2(BadInput,"Data in the @CUSTOM_"<<pluginNameUpperCase()
                    <<" section should be four or five numbers on a single line");

  //Parse and validate values:
  double sigma, hwhm, D_const, hwhm_o;
  if ( ! NC::safe_str2dbl( data.at(0).at(0), sigma )
       || ! NC::safe_str2dbl( data.at(0).at(1), hwhm )
       || ! NC::safe_str2dbl( data.at(0).at(2), D_const )
       || ! NC::safe_str2dbl( data.at(0).at(3), hwhm_o )
       || ! (sigma>0.0) || ! (hwhm>0.0) || ! (D_const>0.0) || ! (hwhm_o>=0.0) )
    NCRYSTAL_THROW2( BadInput,"Invalid values specified in the @CUSTOM_"<<pluginNameUpperCase()
                     <<" sigma, hwhm, D_const (should be three strictly positive floating point values), hwhm_o should be positive" );
    
  int mag_scat;
  // inelastic and elastic magnetic cross sections are both considered
  if ( data.at(0).size()==4 ) mag_scat = 2;
  // only one component is considered
  else if ( ! NC::safe_str2int( data.at(0).at(4), mag_scat )
            || ! (mag_scat > -2) || ! (mag_scat < 2) )
    NCRYSTAL_THROW2( BadInput,"Invalid value specified in the @CUSTOM_"<<pluginNameUpperCase()
                     <<" mag_scat (must be -1 (for down-scattering) or 1 (for up-scattering) or 0 (for elastic scattering)" );
  
  //Getting the temperature
  double temperature = info.getTemperature().get();
    
  //Getting the mean-squared displacement (MSD)
  double msd = 0.;
  if ( info.hasAtomMSD() ) msd = info.getAtomInfos().front().msd().value();

  //Parsing done! Create and return our model:
  return ParamagneticScatter(sigma,hwhm,D_const,temperature,mag_scat,msd,hwhm_o);
}

NCP::ParamagneticScatter::ParamagneticScatter( double sigma, double hwhm, double D_const,
                                               double temperature, int mag_scat, 
                                               double msd, double hwhm_o )
  : m_sigma(sigma),
    m_hwhm(hwhm),
    m_D_const(D_const),
    m_temperature(temperature),
    m_mag_scat(mag_scat),
    m_msd(msd),
    m_hwhm_o(hwhm_o)
{
  //Important note to developers who are using the infrastructure in the
  //testcode/ subdirectory: If you change the number or types of the arguments
  //for the constructor here, you should make sure to perform a corresponding
  //change in three files in the testcode/ directory: _cbindings.py,
  //__init__.py, and NCForPython.cc - that way you can still instantiate your
  //model directly from your python test code).

  nc_assert( m_sigma > 0.0 );
  nc_assert( m_hwhm > 0.0 );
  nc_assert( m_D_const > 0.0 );
  nc_assert( m_temperature > 0.0);
  nc_assert( m_hwhm_o >= 0.0 );
}

double NCP::ParamagneticScatter::calcCrossSection( double neutron_ekin ) const
{
  double xs_in_barn;
  if ( m_mag_scat==2 ) {
    xs_in_barn  = kgfhfunc( m_temperature, neutron_ekin, m_hwhm, m_D_const, -1, m_msd, m_hwhm_o );
    xs_in_barn += kgfhfunc( m_temperature, neutron_ekin, m_hwhm, m_D_const,  0, m_msd, m_hwhm_o );
    xs_in_barn += kgfhfunc( m_temperature, neutron_ekin, m_hwhm, m_D_const,  1, m_msd, m_hwhm_o );
  }
  else xs_in_barn = kgfhfunc( m_temperature, neutron_ekin, m_hwhm, m_D_const, m_mag_scat, m_msd, m_hwhm_o );
  xs_in_barn *= m_sigma;
    
  return xs_in_barn;
}

NCP::ParamagneticScatter::ScatEvent NCP::ParamagneticScatter::sampleScatteringEvent( NC::RNG& rng, double neutron_ekin ) const
{
  ScatEvent result;

  //if ( ! (neutron_ekin > 1.e-3) ) {
    //Special case: We are asked to sample a scattering event for a neutron
    //energy where we have zero cross section! Although in a real simulation we
    //would usually not expect this to happen, users with custom code might
    //still generate such calls. The only consistent thing to do when the cross
    //section is zero is to not change the neutron state parameters, which means:
    //result.ekin_final = neutron_ekin;
    //result.mu = 1.0;
    //return result;
  //}

  //Implement our actual model here. Of course it is trivial for the example
  //model. For a more realistic or complicated model, it might be that
  //additional helper classes or functions should be created and used, in order
  //to keep the code here manageable:

  //result.ekin_final = neutron_ekin;//Elastic
  //result.mu = randIsotropicScatterMu(rng).dbl(); //Isotropic

  if ( m_mag_scat == 2 ) {
        
    double rand = rng.generate(); //random number
    double kgf_down = kgfhfunc( m_temperature, neutron_ekin, m_hwhm, m_D_const, -1, m_msd, m_hwhm_o );
    double kgf_el   = kgfhfunc( m_temperature, neutron_ekin, m_hwhm, m_D_const,  0, m_msd, m_hwhm_o );
    double kgf_up   = kgfhfunc( m_temperature, neutron_ekin, m_hwhm, m_D_const,  1, m_msd, m_hwhm_o );
    double kgf_tot  = kgf_down + kgf_el + kgf_up;
    
    if ( rand < kgf_down / kgf_tot ) {
      //down-scattering happens
      if ( neutron_ekin <= m_D_const ) {
        result.ekin_final = neutron_ekin;
        result.mu = 1.;
      }
      else {
        result.ekin_final = Eppdf( rng, neutron_ekin, m_hwhm, m_D_const, -1, m_msd, m_hwhm_o );
        result.mu         = mupdf( rng, neutron_ekin, m_hwhm, m_D_const, -1, m_msd );
      }
    }
    else if ( rand < (kgf_down + kgf_up) / kgf_tot ) {
      //up-scattering happens
      result.ekin_final = Eppdf( rng, neutron_ekin, m_hwhm, m_D_const, 1, m_msd, m_hwhm_o );
      result.mu         = mupdf( rng, neutron_ekin, m_hwhm, m_D_const, 1, m_msd );
    }
    else {
      //elastic scattering happens
      result.ekin_final = Eppdf( rng, neutron_ekin, m_hwhm, m_D_const, 0, m_msd, m_hwhm_o );
      result.mu         = mupdf( rng, neutron_ekin, m_hwhm, m_D_const, 0, m_msd );
    }
  }
      
  else if ( m_mag_scat == -1 ) {
    if ( neutron_ekin <= m_D_const ) {
      result.ekin_final = neutron_ekin;
      result.mu = 1.;
    }
    else {
      result.ekin_final = Eppdf( rng, neutron_ekin, m_hwhm, m_D_const, -1, m_msd, m_hwhm_o );
      result.mu         = mupdf( rng, neutron_ekin, m_hwhm, m_D_const, -1, m_msd );
    }
  }
    
  else if ( m_mag_scat == 1 ) {
    result.ekin_final = Eppdf( rng, neutron_ekin, m_hwhm, m_D_const, 1, m_msd, m_hwhm_o );
    result.mu         = mupdf( rng, neutron_ekin, m_hwhm, m_D_const, 1, m_msd );
  }
    
  else {
    result.ekin_final = Eppdf( rng, neutron_ekin, m_hwhm, m_D_const, 0, m_msd, m_hwhm_o );
    result.mu         = mupdf( rng, neutron_ekin, m_hwhm, m_D_const, 0, m_msd );
  }

  return result;
}
