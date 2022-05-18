#include "NCParamagneticScatter.hh"

//Include various utilities from NCrystal's internal header files:
#include "NCrystal/internal/NCString.hh"
#include "NCrystal/internal/NCRandUtils.hh"

//Components to calculate the magnetic scattering cross section
double gfunc( double temperature, double D_const, int mag_scat )
{
  //g functions, which accounts for the spin
  //temperature : temperature of material K
  //D_const : zero-field splitting constant, eV
  //mag_scat : magnetic scattering option, int
  //-1, 0, 1 represent respectively down, elastic and up scattering
  //g0=g+
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
  double A, B, f;
  A = 2 * (msd + std::log(2) / (hwhm * hwhm)) / NCrystal::const_hhm; // eV^-1
  if ( mag_scat == -1 ) {
    if (incident_neutron_E > D_const) {
      B = 2 * A * std::sqrt(incident_neutron_E * (incident_neutron_E - D_const));
      f = NCrystal::exp_negarg_approx(-A * (2 * incident_neutron_E - D_const));
      f *= 0.5 * (NCrystal::exp_approx(B) - NCrystal::exp_negarg_approx(-B)) / B;
    }
    else f = 0;
  }
  else if ( mag_scat == 1 ) {
    B = 2 * A * std::sqrt(incident_neutron_E * (incident_neutron_E + D_const));
    f = NCrystal::exp_negarg_approx(-A * (2 * incident_neutron_E + D_const));
    f *= 0.5 * (NCrystal::exp_approx(B) - NCrystal::exp_negarg_approx(-B)) / B;
  }
  else {
    B = 2 * A * incident_neutron_E;
    f = NCrystal::exp_negarg_approx(-B);
    f *= 0.5 * (NCrystal::exp_approx(B) - NCrystal::exp_negarg_approx(-B)) / B;
  }
  return f;
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
  //    <sigmavalue> <hwhm value> <D_const value> <down-scattering option>
  //

  //Verify we have exactly one line and four words:
  if ( data.size() != 1 || data.at(0).size()!=4 )
    NCRYSTAL_THROW2(BadInput,"Data in the @CUSTOM_"<<pluginNameUpperCase()
                    <<" section should be four numbers on a single line");

  //Parse and validate values:
  double sigma, hwhm, D_const;
  if ( ! NC::safe_str2dbl( data.at(0).at(0), sigma )
       || ! NC::safe_str2dbl( data.at(0).at(1), hwhm )
       || ! NC::safe_str2dbl( data.at(0).at(2), D_const )
       || ! (sigma>0.0) || ! (hwhm>0.0) || ! (D_const>0.0) )
    NCRYSTAL_THROW2( BadInput,"Invalid values specified in the @CUSTOM_"<<pluginNameUpperCase()
                     <<" sigma, hwhm, D_const (should be three positive floating point values)" );
    
  int mag_scat;
  if ( ! NC::safe_str2int( data.at(0).at(3), mag_scat )
       || ! (mag_scat > -2) || ! (mag_scat < 2) )
    NCRYSTAL_THROW2( BadInput,"Invalid value specified in the @CUSTOM_"<<pluginNameUpperCase()
                     <<" mag_scat (must be -1 (for down-scattering) or 1 (for up-scattering) or 0 (for elastic scattering)" );
  
  //Getting the temperature
  double temperature = info.getTemperature().get();
    
  //Getting the mean-squared displacement (MSD)
  double msd = 0.;
  if ( info.hasAtomMSD() ) msd = info.getAtomInfos().front().msd().value();

  //Parsing done! Create and return our model:
  return ParamagneticScatter(sigma,hwhm,D_const,temperature,mag_scat,msd);
}

NCP::ParamagneticScatter::ParamagneticScatter( double sigma, double hwhm, double D_const,
                                               double temperature, int mag_scat, double msd)
  : m_sigma(sigma),
    m_hwhm(hwhm),
    m_D_const(D_const),
    m_temperature(temperature),
    m_mag_scat(mag_scat),
    m_msd(msd)
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
}

double NCP::ParamagneticScatter::calcCrossSection( double neutron_ekin ) const
{
  double g = gfunc(m_temperature, m_D_const, m_mag_scat);
  double f = ffunc(neutron_ekin, m_hwhm, m_D_const, m_mag_scat, m_msd);
  double xs_in_barn = m_sigma * g * f;
  if ( m_mag_scat == -1 ) {
    if ( neutron_ekin > m_D_const ) xs_in_barn *= std::sqrt(1 - m_D_const / neutron_ekin);
    else xs_in_barn = 0.;
  }
  else if ( m_mag_scat == 1 ) xs_in_barn *= std::sqrt(1 + m_D_const / neutron_ekin);
    
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
  result.mu = randIsotropicScatterMu(rng).dbl(); //Isotropic
    
  if ( m_mag_scat == -1 ) {
    if ( neutron_ekin > m_D_const ) result.ekin_final = neutron_ekin - m_D_const;
    else result.ekin_final = neutron_ekin;
  }
  else if ( m_mag_scat == 1 ) result.ekin_final = neutron_ekin + m_D_const;
  else result.ekin_final = neutron_ekin;

  return result;
}


