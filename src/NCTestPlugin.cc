
#include "NCTestPlugin.hh"
#include "NCrystal/internal/utils/NCMsg.hh"
#include "NCrystal/internal/utils/NCMath.hh"

void NCP::customPluginTest()
{
  //This function is called by NCrystal after the plugin is loaded, but only if
  //the NCRYSTAL_PLUGIN_RUNTESTS environment variable is set to "1". In case of
  //errors or anything about the plugin not working, simply throw an exception
  //(which is what the nc_assert_always function does below, but feel free to
  //simply throw an exception directly).

  //Note, emit all messages here and elsewhere in plugin code with NCPLUGIN_MSG
  //(or NCPLUGIN_WARN for warnings), never raw usage of std::cout or printf!

  NCPLUGIN_MSG("Testing plugin");

  //Here we add a test which does not really make sense from a physics
  //perspective, but merely test that the model is activated correctly and has
  //an effect on cross sections when the appropriate @CUSTOM_MAGSCAT section
  //is present in NCMAT data.
  //
  //To do so, we take a file (BeO) from stdlib, append a @CUSTOM_MACSCAT
  //section, and note the effect on inelastic cross sections at 1K and 5Aa,
  //which should make it jump from <1e-4barn/atom to around 1 barn/atom (with
  //the given parameters).

  constexpr double custom_incohelas_sigma = 5.0;//barn
  constexpr double custom_incohelas_wl_threshold = 1.2;//angstrom
  const std::string base_data = "stdlib::BeO_sg186.ncmat";
  const std::string cfg_params = ";temp=1K;comp=inelas";

  std::string testdata;
  {
    testdata = NC::FactImpl::createTextData(base_data)->rawDataCopy();
    std::ostringstream customsection;
    customsection << "@CUSTOM_" << pluginNameUpperCase() << "\n"
                  << "3.66  1.5  0.4e-3  0.  1. \n";
    testdata += customsection.str();
  }

  if ( true ) {
    NCPLUGIN_MSG("Test NCMAT data begin:");
    NCRYSTAL_RAWOUT(testdata);
    NCPLUGIN_MSG("Test NCMAT data end.");
  }

  //Let us load and exercise this testdata. We do NOT want to register the
  //testdata as a virtual file, since we want this test to not leave anything
  //registered after it is done running. So we simply put the raw data directly
  //into a MatCfg object (which is the C++ object representing a cfg-string):

  auto cfg = NC::MatCfg::createFromRawData( std::string(testdata),
                                            ";comp=inelas;temp=1K" );
  auto cfg_ref = NC::MatCfg(base_data+cfg_params);

  auto scat = NC::createScatter( cfg );
  auto scat_ref = NC::createScatter( cfg_ref );
  auto wl = NC::NeutronWavelength{ 5.0 };

  auto xs_ref = scat_ref.crossSectionIsotropic( wl );
  auto xs = scat.crossSectionIsotropic( wl );
  NCPLUGIN_MSG("Ref   BeO         XS @ "<<wl<<": "<<xs_ref
               <<" (should be less than 1e-4barn)");
  NCPLUGIN_MSG("Dummy BeO+MagScat XS @ "<<wl<<": "<<xs
               <<" (should be around 1 barn)");

  nc_assert_always( xs_ref.dbl() < 1e-4 );
  nc_assert_always( xs.dbl() > 0.8 && xs.dbl() < 1.2 );

  //Note: We do not test the scatter sampling here, but should probably do so
  //(to test the energy changes are as expected).

  NCPLUGIN_MSG("All tests of plugin were successful!");
}
