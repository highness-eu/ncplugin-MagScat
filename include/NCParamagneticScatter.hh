#ifndef NCPlugin_ParamagneticScatter_hh
#define NCPlugin_ParamagneticScatter_hh

#include "NCrystal/NCPluginBoilerplate.hh"//Common stuff (includes NCrystal
                                          //public API headers, sets up
                                          //namespaces and aliases)

namespace NCPluginNamespace {

  //We implement the magnetic scattering model in this completely custom C++ helper
  //class. That decouples it from NCrystal interfaces (which is nice in case the
  //NCrystal API changes at some point), and it makes it easy to directly
  //instantiate and test the modelling implementation from standalone C++ code.
  //
  //We mark the class as MoveOnly, to make sure it doesn't get copied around by
  //accident (since it could easily end up having large data members).

  class ParamagneticScatter final : public NC::MoveOnly {
  public:

    //A few static helper functions which can extract relevant data from NCInfo
    //objects (the createFromInfo function will raise BadInput exceptions in
    //case of syntax errors in the @CUSTOM_ section data):

    static bool isApplicable( const NC::Info& );
    static ParamagneticScatter createFromInfo( const NC::Info& );//will raise BadInput in case of syntax errors

    //Paramagnetic scattering is considered. The cross section formule is
    //presented in Zimmer's paper.

    //Constructor gets constant cross section value, half width at half maximum,
    //zero-field splitting constant, temperature and magnetic scattering option:
    ParamagneticScatter( double sigma, double hwhm, double D_const,
                         double temperature, int mag_scat, 
                         double msd, double tau );

    //Provide cross sections for a given neutron:
    double calcCrossSection( double neutron_ekin ) const;

    //Sample scattering event (rng is random number stream). Results are given
    //as the final ekin of the neutron and scat_mu which is cos(scattering_angle).
    struct ScatEvent { double ekin_final, mu; };
    ScatEvent sampleScatteringEvent( NC::RNG& rng, double neutron_ekin ) const;

  private:
    //Data members:
    double m_sigma;
    double m_hwhm;
    double m_D_const;
    double m_temperature;
    int    m_mag_scat;
    double m_msd;
    double m_tau;
  };

}
#endif
