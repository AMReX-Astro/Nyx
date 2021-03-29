#ifndef _constants_cosmo_
#define _constants_cosmo_

#ifndef AMREX_USE_FLOAT

// Relation of our code units & CGS 
static constexpr amrex::Real  M_unit = 1.98848e33;    // M_sun
static constexpr amrex::Real  L_unit = 3.0856776e24;  // Mpc
static constexpr amrex::Real  V_unit = 1.e5;          // km/s
static constexpr amrex::Real  T_unit = L_unit/V_unit; // time unit

  
//   Fundamental constants
  
   static constexpr amrex::Real  Gconst = 6.67408e-8
           * M_unit*T_unit*T_unit/(L_unit*L_unit*L_unit);// Newton [g^-1*s^-2*cm^3]
  
static constexpr amrex::Real  k_B    = 1.38064852e-16          
  * T_unit*T_unit/(M_unit*L_unit*L_unit);// Boltzmann [g*cm^2/s^2*K]
  
//static constexpr amrex::Real  hbar   = 1.0545718e-27           
  //  * T_unit/(M_unit*L_unit**2);// Planck/2pi [g*cm^2/s]
  
static constexpr amrex::Real  n_A    = 6.022140857e23          
  * M_unit;// Avogadro's number [mol^-1]

static constexpr amrex::Real  m_nucleon    = 1.0/n_A;
  
static constexpr amrex::Real  m_proton = 1.672621e-24          
  / M_unit;// Proton mass [g]
  
//static constexpr amrex::Real  sigma_T =  6.6524587158e-25         
//  / L_unit**2;// Thomson cross section [cm^2]
  
static constexpr amrex::Real  c_light = 2.99792458e10   
  / V_unit;// Speed of light [cm/s] 

static constexpr amrex::Real  Hubble_const = 100.;    //            Hubble constant / h

  // Useful quantities and conversions
  
  static constexpr amrex::Real  mp_over_kb = m_proton/k_B;

static constexpr amrex::Real  density_to_cgs = M_unit / (L_unit*L_unit*L_unit);

  // For internal energy
  static constexpr amrex::Real  e_to_cgs = V_unit * V_unit;

// For source terms we convert [erg/(s*cm^3) = g/(s^3*cm)] into code units
static constexpr amrex::Real  heat_from_cgs = L_unit*(T_unit*T_unit*T_unit / M_unit);

// For AGN accretion rate
//static constexpr amrex::Real  eddington_const = 4.00*pi * Gconst * m_proton / (sigma_T * c_light);
#else
// Hard-code gravity source type
#define grav_source_type 1

// Relation of our code units & CGS 
#define M_unit ( 1.98848e33)    // M_sun
#define L_unit ( 3.0856776e24)  // Mpc
#define V_unit ( 1.e5)          // km/s
#define T_unit ( L_unit/V_unit) // time unit

  
//   Fundamental constants
  
#define Gconst ( 6.67408e-8 * M_unit*T_unit*T_unit/(L_unit*L_unit*L_unit))// Newton [g^-1*s^-2*cm^3]
  
#define k_B    ( 1.38064852e-16          * T_unit*T_unit/(M_unit*L_unit*L_unit))// Boltzmann [g*cm^2/s^2*K]
  
//#define hbar   ( 1.0545718e-27           
//  * T_unit/(M_unit*L_unit**2))// Planck/2pi [g*cm^2/s]
  
#define n_A    ( 6.022140857e23          * M_unit)// Avogadro's number [mol^-1]

#define m_nucleon    ( 1.0/n_A)
  
#define m_proton ( 1.672621e-24          / M_unit)// Proton mass [g]
  
//#define sigma_T (  6.6524587158e-25         
//  / L_unit**2)// Thomson cross section [cm^2]
  
#define c_light ( 2.99792458e10   / V_unit)// Speed of light [cm/s] 

#define Hubble_const ( 100.)    //            Hubble constant / h

// Useful quantities and conversions
  
#define mp_over_kb ( m_proton/k_B)

#define density_to_cgs ( M_unit / (L_unit*L_unit*L_unit))

// For internal energy
#define e_to_cgs ( V_unit * V_unit)

// For source terms we convert [erg/(s*cm^3) ( g/(s^3*cm)] into code units
#define heat_from_cgs ( L_unit*(T_unit*T_unit*T_unit / M_unit))

// For AGN accretion rate
//#define eddington_const ( 4.00*pi * Gconst * m_proton / (sigma_T * c_light))

#endif
#endif
