/**
 * @file montecarlo.C
 * @brief ROOT macro to be used with CAFAna executable to run the selection.
 * @author justin.mueller@colostate.edu
*/
#include "include/analysis.h"
#include "include/container.h"
#include "include/csv_maker.h"
#include "sbnana/CAFAna/Core/Binning.h"

using namespace ana;

/**
 * The main function of the selection. Creates a container for the CAFAna
 * Spectrum objects and populates it with a variety of variables that define
 * the selection.
 * @return none.
*/

void mc_analyzer()
{

  //SpecContainer spectra("/pnfs/icarus/persistent/users/lkashur/numi_nu_cosmic_v09_89_01p01/output_flat/*.flat.root", "spectra_numi_mc_test.root", -1, -1);
  SpecContainer spectra("/pnfs/icarus/persistent/users/lkashur/bnb_nu_cosmic_v09_84_00_01/output_flat/*.flat.root", "spectra_bnb_mc_test.root", -1, -1);
  spectra.add_spectrum1d("sInfoVar", Binning::Simple(1, 0, 2), kInfoVar); // Analysis
  //  spectra.add_spectrum1d("sSignalVar", Binning::Simple(1, 0, 2), kSignalVar); // Signal/Efficiency
  //spectra.add_spectrum1d("sSelectedVar", Binning::Simple(1, 0, 2), kSelectedVar); // Selected/Purity
  spectra.run();
  
}
