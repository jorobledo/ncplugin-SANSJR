#include "NCPhysicsModel.hh"
#include <iostream>

//Include various utilities from NCrystal's internal header files:
#include "NCrystal/internal/NCString.hh"
#include "NCrystal/internal/NCRandUtils.hh"
#include "sasmodels/sas_core_shell_sphere.c"
#include "NCrystal/internal/NCMath.hh"

bool NCP::PhysicsModel::isApplicable( const NC::Info& info )
{
  //Accept if input is NCMAT data with @CUSTOM_<pluginname> section:
  return info.countCustomSections(pluginNameUpperCase()) > 0;
}

NCP::PhysicsModel NCP::PhysicsModel::createFromInfo( const NC::Info& info )
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
  // case, we want to accept sections of the form (units are barn and angstrom as
  // is usual in NCrystal):
  //
  // @CUSTOM_<ourpluginname>
  //    <sigmavalue> <wavelength threshold value>
  //

  //Verify we have exactly one line and two words:
  if ( data.size() != 1 || data.at(0).size()!=2 )
    NCRYSTAL_THROW2(BadInput,"Data in the @CUSTOM_"<<pluginNameUpperCase()
                    <<" section should be two numbers on a single line");

  //Parse and validate values:
  double radius, thickness, sld_core, sld_shell, sld_solvent, result;
  // if ( ! NC::safe_str2dbl( data.at(0).at(0), sigma )
  //      || ! NC::safe_str2dbl( data.at(0).at(1), lambda_cutoff )
  //      || ! (sigma>0.0) || !(lambda_cutoff>=0.0) )
  //   NCRYSTAL_THROW2( BadInput,"Invalid values specified in the @CUSTOM_"<<pluginNameUpperCase()
  //                    <<" section (should be two positive floating point values)" );

  //Parsing done! Create and return our model:
  return PhysicsModel(radius, thickness, sld_core, sld_shell, sld_solvent);
}

NCP::PhysicsModel::PhysicsModel( double radius, double thickness, double sld_core, double sld_shell, double sld_solvent)
  : m_radius(radius),
    m_thickness(thickness),
    m_sld_core(sld_core),
    m_sld_shell(sld_shell),
    m_sld_solvent(sld_solvent)
{
  //Important note to developers who are using the infrastructure in the
  //testcode/ subdirectory: If you change the number or types of the arguments
  //for the constructor here, you should make sure to perform a corresponding
  //change in three files in the testcode/ directory: _cbindings.py,
  //__init__.py, and NCForPython.cc - that way you can still instantiate your
  //model directly from your python test code).

  nc_assert( m_radius > 0.0 );
  nc_assert( m_thickness > 0.0);
}

double NCP::PhysicsModel::calcIQ(double Q) const
{
  double F1=0.0, F2=0.0;
  Fq(Q, &F1, &F2, m_radius, m_thickness, m_sld_core, m_sld_shell, m_sld_solvent);

  return F2;
}

double NCP::PhysicsModel::calcDiffCrossSection( double neutron_ekin ) const
{

  double q_min = 0.0;
  double lastk = 2.*NC::kPi/NC::ekin2wl(neutron_ekin);
  unsigned int n = 100;
  // std::cout << "Energy:" << neutron_ekin << ", last k:" << lastk << std::endl;

  // CHECK FIX
  auto integrand = [this](double Q){return Q * this->calcIQ(Q);};


  double xs = NC::integrateSimpsons(integrand, q_min, 2*lastk, n);

  double vol;
  double ksq = lastk * lastk;
  xs *= form_volume(m_radius, m_thickness)/(ksq);
  return xs;
}

NCP::PhysicsModel::ScatEvent NCP::PhysicsModel::sampleScatteringEvent( NC::RNG& rng, double neutron_ekin ) const
{
  ScatEvent result;

  if ( ! (neutron_ekin > 0) ) {
    //Special case: We are asked to sample a scattering event for a neutron
    //energy where we have zero cross section! Although in a real simulation we
    //would usually not expect this to happen, users with custom code might
    //still generate such calls. The only consistent thing to do when the cross
    //section is zero is to not change the neutron state parameters, which means:
    result.ekin_final = neutron_ekin;
    result.mu = 1.0;
    return result;
  }

  //Implement our actual model here. Of course it is trivial for the example
  //model. For a more realistic or complicated model, it might be that
  //additional helper classes or functions should be created and used, in order
  //to keep the code here manageable:

  result.ekin_final = neutron_ekin;//Elastic
  result.mu = randIsotropicScatterMu(rng).dbl();

  return result;
}


