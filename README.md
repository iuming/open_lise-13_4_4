# open_lise++-13_4_4
Readme and Last Changes. 

Readme

The program LISE++ is designed to predict the intensity and purity of secondary ion beams for experiments using radioactive beams produced by in-flight separators. The LISE++ program also facilitates tuning experiments where its results can be quickly compared to on-line data. The LISE++ package includes configuration files for most of the existing fragment and recoil separators. Projectile fragmentation, fusion-evaporation, fusion-fission, Coulomb fission, abrasion-fission and two body nuclear reactions models are included in this program and can be used as the production reaction mechanism to simulate experiments at beam energies above the Coulomb barrier.

The LISE++ program contains more features and options than those described in the documentation at the official site. You are strongly encouraged to experiment with them and see the effects they have on the results. A very large amount of physics is incorporated in this program, from projectile fragmentation, fission and other reaction mechanism models, cross section systematics, electron stripping models, energy loss models to beam optics, just to list a few. All the references for the calculations are directly accessible within the program itself (see the various option windows) and you are encouraged to consult them for detailed information.

The LISE++ name is borrowed from the well known evolution of the C programming language, and is meant to indicate that the program is no longer limited to a fixed configuration like it was in the original “LISE” program, but can be configured to match any type of device or add to an existing device using the concept of blocks. 

The program is constantly expanding and evolving using the feedback of users around the world. At the time of this writing, many “satellite” tools have been incorporated into the LISE++ framework, which are accessible with buttons on the main toolbar and include:
·	Physical calculator
·	Kinematics calculator
·	Evaporation calculator
·	Units converter
·	Mathematical calculator
·	The program PACE4 (fusion-evaporation code) by A.Gavron et all.
·	Spectrometric calculator by J. Kantele
·	The program CHARGE by Th. Stöhlker et al. 
·	The program GLOBAL by W. E. Meyerhof et al.
·	The program BI (search for 2-dimentional peaks)
·	MOTER by H.A.Thiessen et al.: raytracing code with optimization capabilities operating under MS Windows

Citations of this program please use the following: 
O.B.Tarasov and D.Bazin, NIM B (2008) 4657-4664

Official LISE++ web site: http://lise.nscl.msu.edu

The authors wish to thank Prof. B.M. Sherrill, Prof. D.J. Morrissey, and Dr. H. Weick for fruitful discussions. The authors thank Drs. M. Lewitowicz and O. Sorlin for their help and support in the development of the first versions of the LISE code, Drs. M. Portillo and M. Hausmann for extensive testing the LISE++ code, their remarks, and requests. This work was supported by DOE #DE-FG02-06ER41413, NSF #PHY-01-10253, -06-06007, and -11-02511 grants.


Copyright © 1985-2015 by 
LISE++ group development (Oleg B. Tarasov & Daniel Bazin) 
National Superconducting Cyclotron Laboratory @ Michigan State University.
All rights reserved.

_________________________________________________

Last Changes
(see complete list @ http://lise.nscl.msu.edu/changes.html)

______________________________________________
9.8.4	01/08/17

* Modification of EPAX 3.1 for secondary reactions
* Corrections in the Kinematics calculator dialog for the case A_lab=0
______________________________________________

9.7 	11/01/13    

* New official version 9.7

* MARS spectrometer extended configuration in LISE++ package 
* Plotting Envelopes with rotation blocks:
* neutron and gamma induced reactions in the Kinematic Calculator
* Gas Cell utility modification
* Correction in Monte Carlo Eloss and Range plots
* New block :  "Shift”
* E-dipole, E-quad: Analytical & MC solutions for Electrostatical dipoles
* New address for dynamical menu. New indexation!!!
* Correction in Multiple Reactions use  (when AF and PF togehter, and AF is the first)
* Calculation backward transmission (realistic)
* User DifCS in MC transmission
* Utility menu : Differential cross sections utility, LAB<->CM converter 
* Input of ions rays after target in MC
* several locations for output MC file
* A,Z,Q filters in MC
______________________________________________
version 9.5        06.03.13     

* new official version

______________________________________________

* lise2013.dbf -- update
* User AME2012.lme
* new DBF time string "private"
______________________________________________
* RF_Buncher 
- configurations revision
- analytical calculations
- Graph_goodies (ellipse plot)

* Pairing energies (P_x) plots
* lise2013.dbf, user_mass_excess_2012.lme  -- new default database and mass file

* A1900_2013.lcn, A1900_2013.lopt -- new default files
______________________________________________

* Acculinna2 extended configuration and example
* Angular straggling plot modification
______________________________________________
* RESOLUT configurations update
* Optics Quad Setup : corrections for second order labels
* Resolut 3gap configurations
  Buncher block proerty mdifcation

* Solenoid block dialog
  - Help link
  - changes in the length message
  - LengthWienFilter modification

* Corrections for Plot25 due to abundance values in other modes 
* Fixed: Bug in the case charge state & stripper 
______________________________________________

version  9.4.80    21.01.13     

* RF buncher
 - RESOLUT.lcn in the package (also GXPF1B & GXPF1B5 masses)
 - RF buncher in MC calculations
 - RF buncher : tuning method

______________________________________________

version  9.4.74   03.01.13

* Abundance value of stable isotopes in T1/2 plot
* Database half-life time scale has been modified (now from 10^-21sec)

______________________________________________

version  9.4.66   07.12.12 

* modification  in "decay analysis": evap+res -> Parent  and  Init Cross section
* modification : interpolation  @ distribution2 
* correction in Abrasion-Ablation for Breakup  -> interpolation2
* fission NP increased for n-evaporation
* modification in thermalization subroutine for new axis
* corrections for _o_mass.cpp  
* RGD - plot -- more digits in output file
______________________________________________

version  9.4.56   06.11.12

* New Dn and Dp database calculated values
* Read/Write in/to file correction for EPAX mode
_________________________________________________________________________ 

version  9.4.54   30.10.12    

* Modifications in the "Global" code and in LISE_global.dll
* Secondary reactions : accounting reduced cross sections in the 1D-mode
______________________________________________

version  9.4.52   26.10.12     

* NRV - partner site : links from some LISE++ dialogs and menu
* New user ME files: AME2011 & AME2011+GXPF1 (Z=16-22)
* AA correction for the Breakup channel
* No FTP server more. Download from HTTP
* Correction for DB plots (extrapolation for User ME mass correction)
* Update: Analysis of prefragment plots
______________________________________________

version  9.4.46  12.09.12    

* Modification of iso-files (table2012, discovery, discovery_lab)
* thermalization of excitation energy based on NPA531 (1991) 709
______________________________________________
version  9.4.40  28.08.12     

* Choice of EPAX model for secondary reactions  @ SR dialog
* Fix of bug for reaction mechanism parameter in the case of AF
______________________________________________
version  9.4.34  27.08.12     

* Utility menu: Plots of NSCL and FRIB rates
______________________________________________
version  9.4.32  14.08.12

* S_d & S_t plots in the Database menu
* Possibility to turn on/off a,p,n channels
* Prefragment search utility modification
* Modification for Initial prefgraments plot for AA @ the Evaporation Calculator

______________________________________________
version  9.4.28  30.07.12

*  Initial prefgraments plot for AA @ the Evaporation Calculator
______________________________________________
version  9.4.24  18.07.12

* new block "RF_buncher" (still under construction) 
  - dialog
  - RF-buncher plot
  - RF-separator correction :  instead r->Gtv use r->Gtv_real
  - no RF resolution in plot options, new bunch length in the beam dialog

* Reaction mechanism
  - new spontaneous fission formula (Karpov,Zagrebaev)
  - probability for compound nucleus formation P_CN (E,l) V.Zagrebaev & W.Greiner, PRC78, 034610 (2008)
  - new excitation plot with PCN option for Fusion-Residueal and Fusion-Fission

______________________________________________
version 9.4.15   25.06.12  

* new option : Zp_power for Leon's model might be changed 
* modification for non-equilibrium case for Global code:  Q_INIT = max(Z - 28, Q_INIT)

* new default option file is A1900_2012.lopt
* new option file A1900_2012.lopt
* new option file GSI_FRS_2012.lopt

______________________________________________
version 9.4.12   18.06.12  

* new GSI configuration  "FRS - TA2 - CaveC (2012).lcn"		
______________________________________________
version 9.4.11     11.05.12     

* EPAX 3.01
* Abrasion-Ablation revision
______________________________________________
version 9.4.04   beta    17.04.12     

* New S800 extended configuration
* all A1900 configurations have been changed: 
   a. wedge defect thickness = 0.3%
   b. extended configurations : beam emittance, tuning dipole after target, I2-wedge 

______________________________________________
Version 9.4.03   beta    16.03.12

* New LISE++ site : lise.nscl.msu.edu
______________________________________________
Version 9.4.03   beta    16.03.12

* New LISE++ site : lise.nscl.msu.edu

______________________________________________
version 9.4.02 & 9.3.12   16.03.12     

* Fixed: calculated transmission has been cleared during reading file from disk.
 !!! Important for versions  9.3.10-11  and 9.4.1 beta!!!!

______________________________________________
Version 9.4.01   beta    15 .03.12

* Asymmetrical momentum distributions
______________________________________________
Version 9.3.11 13.03.12

* Checkbox for Faraday cup block to make it disable from the Setup window
* Speed Optimization for Global matrices recalculation
* Convolution Momentum distribution analysis plots from the Universal
  parameterization (Convolution) momentum distribution dialog
______________________________________________
Version 9.3.8    02.03.12

* Correction in the Dispersion method for secondary target
* Address of the LISE++ site ("groups/lise") has been changed.
  Avoiding warning message, when LISE++ is checking for new version
______________________________________________
Version 9.3.6    22.02.12

* Apertures for Rotation Blocks are disable by default
* LowMatrixLimit set to 1e-10 for XX,YY,TT,PP component. 
  This condition is not applied to Rotation Blocks
* Monte Carlo:  Passing Aperture bug for X-use -- fixed
* If Optical Block Length equal to 0, then bounds are used just for one point (like slits)
* Envelope plots: axis title has been changed
______________________________________________
Version 9.3.5    16.02.12

* New properties of identification and label fonts for 2D-plots
* Message for old version calculations

______________________________________________
Version 9.3.3    23.01.12     

* The Preference dialog: choice of the LISE++ working directory 
   for users with administrative privileges 
* Initialization of the density array has been moved at beginning
* New GSI configuration files
______________________________________________
Version 9.3      06.01.12     

* New official version
______________________________________________
Version 9.2.169  28.12.11     

* New names for elements 114 & 116
______________________________________________
Version 9.2.168  27.12.11     

* New isotope discovery data
* Modifications in "discovery history" to read isotopes of light elements  (Z<10)
* Changes in the "secondary reactions" calculations (case of zero or negative primary cross sections)
* Changes in the "secondary reactions" dialog
______________________________________________
Version 9.2.164  13.12.11     

* modification in the evaporation cascade subroutines
* Abrasion-Fission 3EER : using max values instead mean values for E.E.R.
* fix some insignificant bugs 
  => Evaporation calculator : chose nucleus
  => Decay analysis dialog : absent window
______________________________________________
Version 9.2.155  08.12.11     

* element, rate and background colors have been changed
* new order of decays. New file "table2012.iso" with decay information
* Option of element, rate and background colors editing in the Isotope dialog

______________________________________________
Version 9.2.151  06.12.11     

* Revision of decays and half-lives
* Correction in half-life calculation in the database and isotope dialogs

* Extended version of the isotope dialog : "Decay Analysis"
   => Call the "Decay analysis" dialog from the "Utility" menu & statistics window
   => Double click left mouse button - show only database information
   => "conflict" analysis of decay branch 

* New decay spontaneous fission formula : maximum of all them
* revision of constant masses
* corrections for 1D plot legends

*  New decays : 
    => b- & a
    => p & b+
    => p & a
    => SF & b+
    => SF & b-
______________________________________________
Version 9.2.134 30.11.11     

* "inter"  modification for Distribution class  (important! impact to transmission calculation)
* "inter2" modification for Distribution2 class (important! impact to transmission calculation)
* Show discovery history availability in the chart of nuclides
* New functions in LISE for EXCEL (inter, A1900 dipoles and so on)
* Update of Isotope Discoveries for elements Ca,In,Sn,Pt
* Fix of the new link for of JAEA (Japan Atomic Nuclear Agency Chart) 2010 Nuclear Chart
* Call LISE++ database dialog from Transmission results window
* Beta+ & Beta- results in Transmission results window

______________________________________________
Version 9.2.126  11.11.11     

<http://groups.nscl.msu.edu/lise/9_2/9_2_126/9_2_126.pdf>

* Spontaneous fission dialog & plots
* Yield Plot gated for downstream block 
* Message (Gauge) how many blocks have remained at reading file 

* Options for PID plot :
 1) white color is available
 2) "+4,+5,+6" sizes are available
 3) option using "s" or "sec" is available.

* Correction for Electromagnetic Excitation energy plot
______________________________________________
Version 9.2.118  04.11.11	

* Correction: Z=112  -> "Cn" name
* Fixed for the potential (fission): Inverse values for shells
* Spontaneous fission dialog (template)
* Modification of "Utilities" menu
______________________________________________
Version 9.2.113  29.07.11

* Modification for taulise_angle ("wedge_curiosity" with monoeneregtic wedge)
* Writing energy deposition 2D plot to file
______________________________________________
version 9.2.111  26.07.11

* Modification for GFS block
* Quadrupole dialog: corrections to view a charge stte of setting fragment
* PACE IV  v.4.19.2
______________________________________________
version 9.2.107  07.06.11

* new NSCL high order optics configuration 
"A1900_extended_COSY_S800BL_LISE.lcn", where
A1900 based on COSY maps, and transfer hall and S800BL are based on LISE++ 2nd order calculations
______________________________________________
version 9.2.106  27.05.11

* Corrections in dipole matrix calculations
* Corrections in high order global matrix calculations
* NSCL high order optics configuration update
______________________________________________
version 9.2.102  25.05.11

* New extended NSCL configurations (A1900, S800BL, A1900+S800BL)
* "Link" sign in Quad-Edit utility
* New isotopes and new element names update
______________________________________________
version 9.2.99   20.05.11

* Correction for Ellipse Aperture MC transmission
* Icons for Quads
______________________________________________
version 9.2.97   12.05.11

* The Slits dialog modification :
 - Angular acceptance for each block ellipse or rectangle
 - Aperture option:  use for MC calculation
 - Aperture option:  rectangle & ellipse

* Options "Bounds On/Off" for All MC calculations
* Corrections for Envelope transmission statistics
* Correction for read/write slits characteristics in/from file
______________________________________________
version 9.2.88   03.05.11

* Beam emittance shapes
* The "Setup" dialog -- new information cell: number of blocks value 
* The MC dialog - new button: the same block as X or as Y
______________________________________________
version 9.2.85   26.04.11

 Submenu "Optics" in "Calculations" menu

* Quad & Dipoles settings: EDIT
* Quad & Dipoles settings: View & Print
* First order matrix elements: view, print
* First order matrix elements: PLOT

______________________________________________
version 9.2.78 beta  19.03.11

* Correction in contour analysis for Mode=4 (spectrum from file)
* PACE4 (v.4.18.4)    Gamma and Gate for A,Z for detailed analysis
* InitRandom(void)  in LISE++
* PACE4 (v.4.18.3)  new random generator (taken from LISE++)
______________________________________________
version 9.2.75 beta  06.03.11

* PACE 4.18 -- Detailed analysis of emitted particles
* Defect thickness can be equal to 0
* Absorbed dose dialog from the target dialog
______________________________________________
version 9.2.71 beta  04.02.11

* Ranges of MC_tab, CrossSec and SecReact increased up to 15 000
* User DB file from the Production mechanism
* Correction for Abrasion-Fission regions restore procedure
______________________________________________
version 9.2.68 beta  01.02.11

* modifications for Abrasion-Fission 3-EER dialogue
* wedge curiosity revision
______________________________________________
version 9.2.66 beta  26.01.11

* DF4 distributions class has been updated:
  + New dX_dP & dY_dP components
  + Debug plots have been modified: two windows 
  + new classes "sigma_dP" in BLOCK, and "sigma_dS" in BLOCK_Component
  + modification in transmission calculations due to this update

* Monte Carlo :
  + corrections in transmission statistics : reactions in material
  + corrections in energy deposition : start point
  + modifications in the transmission statistics window
  + modifications in material passing subroutine
______________________________________________
version 9.2.57 beta  21.01.11
* MC rays generator: new option "Range"
* Modifications in the MC and MC rays generator dialogs 
* Range of  Momentum distribution plot for the Convolution model has been increased

* Customizable Chart of the Nuclides:
 - USing user iso and isolit files
 - Axis, Grid and Border options 
 - Identification in 2D plot mode=25,35 (CS,T12 and so on)
 - Option "Seconds" or "Days,Years" for T_1/2 2D-plots
 - Plot.c virtual size has been changed from 1000 to 10000
 - Graph.c virtual size has been renamed to VirtSizeG, and still equal to 1000
 - subscript text in Graphics
______________________________________________
version 9.2.43 beta  07.01.11

* Thickness defect for Target and Stripper
* New option for stripper : use  for charge state calculations
______________________________________________
version 9.2.38 beta  23.12.10	

* Stripper foil utility :  Rotating target + Pulsing beam
* Stripper foil utility :  Initial temperature
* Stripper foil utility :  Rotating target -> using X_radius instead reduced one,
  link for LISE documentation, correction for gaussian beam density distribution

* Correction due to "superlong" configurations 
______________________________________________
version 9.2.33 beta  10.12.10	

* Dipole (dispersive block): Transport solution 
  of 1st and 2-nd orders including  fringing fields 

* New utility dialog: "The First- and Second-Order Matrix Elements for an Ideal Magnet"
* Increasing possible number of blocks from 94 to 194
* Cuts Edge effect option for transmission calculation  (the "Option" dialog)
* correction in link with cosy matrices
* correction in the transmission statistics dialog
* correction: indexation in dynamical submenu has been changed
______________________________________________
version 9.2.24 beta  19.11.10

* new option for quadrupoles : keep matrices, calculate fields
* Update COSY links for optical blocks
* Calculate automatically Drift-block matrices
* Analyzing 2D-historgam ROOT files by the BI code
* Modification in plot drawing
* Fixed: Crash in ShowCalc (Increasing number of rows and string size)
* Analyzing 1D-historgam ROOT files by the BI code
* Fixed : Pseudo MC plot: crash in rotation blocks
______________________________________________
version 9.2.16 beta  23.09.10

* Monte Carlo : List of isotopes. Modifications to get initial yields
______________________________________________
version 9.2.15 beta  15.09.10

* Corrections for TKE-calculations in pseudo MC plot
* Discovery information has been implemented for Z=23,26,27,36,47,48,56
* The Line of observed neutron-rich isotopes for Z=63,72-89 has been updated
______________________________________________
version 9.2.13 beta  07.09.10

Monte Carlo simulation
* Normalization in Energy deposition for group of isotopes
* Modification for Energy loss in detector (range), where particle is stopped
* increase of the limit of the number of isotopes (1e4) for the MC group
* increase of the limit of the number of rays (1e6) to write in output file
______________________________________________
version 9.2.8 beta  27.08.10

* Correction in MC transmission calculation against crash 
 in the reaction place subroutine at negative energy   
______________________________________________
version 9.2.7 beta  17.08.10

* modifications in angular transmission for the fission case
* Celement::name ->correction for Abrasion-Fission reactions 
* c_plot.cpp  : corrections for Envelope mode -- show all in the case of fission
* modifications for the case : fission & beam inclination (Distribution mode)
* modifications for the case : fission & low energy  (MC mode)
______________________________________________

* new version PACE4 (4.17) from Jenna Smith (smithj@nscl.msu.edu):
   1. GILBERT & CAMERON parameter for LITTLE-A is FACLA==0
   2. Result output modifications
______________________________________________
version 9.2.3 beta  12.08.10

* modifications in MakeAngOptic: no shifts
* modifications in BlockOptics : local coefficients changed
______________________________________________
version 9.2.2 beta  11.08.10

* taulise anglular :: TR_A -- changes for mom_store and its gate
* correction for Calc_list in Faraday cup case
* correction for MakeAngOptic. Substitution for Extrapolation function
* 1st attempt for rotation blocks
______________________________________________
version 9.2.1 beta  10.08.10

* DF4 class:  new AX nad Ay distributions
* new parameters: i_sigma_mrad[]
* DF4 debug for angular distributions
* DF4prev in taulise (anglular and optical)
* DF4 debug plot: current and previous modes
* defocusing in taulise_optic
______________________________________________
version 9.1.19      09.08.10

* Corrections: Target-wedge, Wedge-wedge and Kicker optimizations crash 
______________________________________________
version 9.1.18      05.08.10

* Changing Local to Global matrices for accumulation uncertainty coefficient "mm_%"
______________________________________________
version 9.1.17      04.08.10

* New NSCL configuration and option files
______________________________________________
version 9.1.16      30.07.10

* Correction to plot calculated primary beam charge states after opening file
* Charge state plot. Non-equilibrium case has been excluded
* Brho calculations: modifications in the get_spline subroutine
* Fixed : loading LPP-files with AF production mechanism
______________________________________________
version 9.1.12      20.07.10

* New 45 isotopes from the RIKEN experiment have been implemented
* Charge state after reaction in the Production mechanism dialog
* Correction for Pseudo MC plot detector resolution
______________________________________________
version 9.1.9       30.05.10

* Corrections for Primary beam transmission without target
______________________________________________
version 9.1.8       28.05.10

 * Corrections for Fission reactions in MC calculations - Isotope group option
______________________________________________
version 9.1.5       21.05.10

 * Evaporation calculator: output CS files format has been changed
______________________________________________
version 9.1.4       19.05.10
 
 * Correction for the MC mode in the case of Moyal shape of energy straggling
______________________________________________
version 9.1.3       12.05.10
 
 * List of isotopes in MC mode
 * Escape and gauge in MC plot accumulation
______________________________________________
version 9.1          30.04.10

 * New official version
______________________________________________
version 9.0.43      26.04.10

 * The MC dialog has been changed
 * Transmission modification in case of fission without angular acceptance
______________________________________________
version 9.0.39      19.04.10

* New configuration "A1900_expanded"
* Statistics by block for Monte Carlo Envelope mode
* Using physical aperture to calculate transmission inside of blocks for 
  MC envelope mode
______________________________________________
version 9.0.34      15.04.10

* Quadrupole and Sextupole modes in the Drift block
  with possibilities to calculate optical matrices (1st and 2nd order) 

* Modifications in LISE_Excel : corrections for GLOBAL calculations
______________________________________________
version 9.0.28      12.04.10

* Secondary reaction area has been increased for fragmentation
* New MARS spectrometer configuration, The Discovery button in the Isotopes dialog is disable
______________________________________________
version 9.0.26      04.04.10

* Corection for fission case in the new version of transmission calculations
* Modifications in the optical block transmission
______________________________________________
version 9.0.23      30.03.10

* Modifications in the wedge transmission subroutine
______________________________________________
version 9.0.22      25.03.10

 * The Gas Density dialog is available from the  MCP144 utility
 * Setup window timer update
______________________________________________
version 9.0.20      23.03.10

 * Modification of the MC gate dialog
 * Correction for the case: Wedge thickness =0, angle!=0
______________________________________________
version 9.0.18      20.03.10

 * Modifications in the Discovery database. As,W,Au have been added
 * Several gates for MC calculations
 * Block position for the Envelope mode in MC calculation
______________________________________________
version 9.0.14      15.03.10

 * The “Envelope” mode in Monte Carlo calculations
______________________________________________
version 9.0.10      10.03.10

 * Transmission statistics block by block for Monte Carlo calculations
______________________________________________
version 9.0.08      05.03.10

 * Preparation for 64-bit Windows 7
 * Preparation for new installer
 * new main_init arrays
 * new "version" class
 * Update User fies menu
 * GlobalPermissionSetupRedraw
 * test_ESC was killer in AbrasionFission (Low) in o_Manage_AA.cpp

______________________________________________
version 8.5.51      01.03.10

 * New LISE fragment separator configurations files
______________________________________________
version 8.5.50      28.02.10

 * MC modifications: ray generator to file. + Energy loss and Cross sections
 * MC modifications for counters and speed calcualtion: no timer, no division
______________________________________________
version 8.5.47      18.02.10

 * Modifications of ToF calculations in MC mode
______________________________________________
version 8.5.46      16.02.10

 * Option on/off for a material to be used in Q-state calculations 
   (dialoges "Material", "Wedge", "Target")
 * Corrections for Q-state calculations in MC mode
______________________________________________
version 8.5.44      09.02.10

 * Corrections for MC high order calculations
______________________________________________
version 8.5.43      04.02.10

 Wien filter revision:
 * "floating" matrix solution: dispersion is calculated 
    for each ion for MC & distribution methods
 * Pseudo MC and ellipses plot were corrected as well
 * Wien filter dialog modification (lengths,disp.coefficient)
  http://groups.nscl.msu.edu/lise/8_5/8_5_040_Wien.pdf

 * Ellipse plot: possibility to write data in file
______________________________________________
version 8.5.38      02.02.10

 * Access to Discovery of isotopes 
 * Element names in the Isotope, Database dialog and Show Statistics window
 * Monte Carlo transmission: Isotope group calculation  -> secondary target
______________________________________________
version 8.5.34      28.01.10

 * Monte Carlo transmission: Isotope group calculation  
______________________________________________
version 8.5.28      15.01.10

 * Custom shape degrader:  polynomial solution
   http://groups.nscl.msu.edu/lise/8_5/8_5_028_CustomShapeDegrader.pdf

 * Database plot: option "decay mode"
 * Plotting data : "NO LINE" option through the Method drawing dialog
 * Plotting errors : option "turn  off"
 * Plotting Legend or caption : option "turn  off"
 * "BI" code: corrections for 1-dimensional file 
 * Solenoid block in the "Setup" dialog: modifications 
 * Brho scanning utility: modifications for fission case

______________________________________________
version 8.5.17      08.01.10

 * New option "Custom shape" for a energy-loss degrader ("wedge")
   Available trhough the "Wedge" dialog
______________________________________________
version 8.5.11      15.12.09

 * Several modifications incuding the correction of E=f(B,v) function for Wien-filter 
______________________________________________
version 8.5.8       05.11.09

 * New version of ISOL catcher includes fusion-residual, user excitation files,
   materials fron and behind the target
 * Option : Y-title turn on 180 degrees for screen (Mac's case)
 * Erasing space in COSY log-files for high order optics. MP's request
 * Magnification optic coefficients should be not equal to zero. LISE++ checks and corrects.
 * Modification for Abrasion-Fission reaction if User Cross section & Secondary reactions options set
 * Changes to get new option and configuration files for users without admin.privileges
 * Correction in Monte Carlo transmission calculations for Length 
______________________________________________
version 8.4.1       17.09.09

 * Official version 8.4
