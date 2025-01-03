/**
 * @file csv_maker.h
 * @brief Header file defining a dummy SpillMultiVar for dumping particle info.
 * @author justin.mueller@colostate.edu
*/
#ifndef CSV_MAKER_H
#define CSV_MAKER_H

#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include "cuts.h"
#include "variables.h"
#include "sbnana/CAFAna/Core/MultiVar.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

/**
 * MC
 */
//std::ofstream output("output_mc_numi_test_mu100ph50.log"); // analysis
//std::ofstream output("output_mc_numi_1mu0pi1pi0_masscut.log");
std::ofstream output("output_mc_bnb_1mu0pi1pi0_masscut_TEST.log");
//std::ofstream output("output_mc_signal.log"); // signal/efficiency
//std::ofstream output("output_mc_selected.log"); // selected/purity

/**
 * DATA
 */
//std::ofstream output("output_data_numi_test.log"); // analysis

#define GUARD(VAL) std::isinf(VAL) ? -9999 : VAL
#define OUT(STREAM,TAG) STREAM << std::fixed << TAG << ","
#define CSV(VAL) VAL << ","

/**
 * Writes information about a failed containment cut.
 * @param sr is an SRSpillProxy that attaches to the StandardRecord of the
 * current spill.
 * @param i the truth interaction (signal)
 * @return None.
*/
void write_file_info(const caf::SRSpillProxy* sr, const caf::SRInteractionTruthDLPProxy& i)
{
    output  << CSV(sr->hdr.run) << CSV(sr->hdr.evt) << CSV(sr->hdr.subrun)
            << CSV(i.nu_id) <<CSV(i.momentum[0]) << CSV(vars::id(i))
            << CSV(std::string(sr->hdr.sourceName))
            << std::endl;
}

/**
 * Writes reconstructed variables (truth and reco) for selected/signal
 * interactions.
 * @param sr is an SRSpillProxy that attaches to the StandardRecord of the
 * current spill.
 * @param i the truth interaction (signal)
 * @param j the reco interaction (selected).
 * @return None.
 */
void write_pair(const caf::SRSpillProxy* sr, const caf::SRInteractionTruthDLPProxy& i, const caf::SRInteractionDLPProxy& j)
{
  output  << CSV(sr->hdr.run) << CSV(sr->hdr.evt) << CSV(sr->hdr.subrun) // USE FOR EVERYTHING ELSE
            //<< CSV(i.nu_energy_init + i.nu_position[2]) << CSV(sr->hdr.evt) << CSV(sr->hdr.subrun) // USE FOR DETECTOR SYSTEMATICS
	  << CSV(std::string(sr->hdr.sourceName))
          << CSV(i.nu_id) << CSV(vars::id(i)) << CSV(vars::id(j)) << CSV(sr->hdr.triggerinfo.global_trigger_det_time)
	  << CSV(vars::category(i))
	  << CSV(vars::category_topology(i))
	  << CSV(vars::category_interaction_mode(i))
	  << CSV(vars::muon_momentum_mag(i))
	  << CSV(vars::muon_momentum_mag(j))
	  << CSV(vars::muon_beam_costheta(i))
	  << CSV(vars::muon_beam_costheta(j))
	  << CSV(vars::pi0_leading_photon_energy(i))
	  << CSV(vars::pi0_leading_photon_energy(j))
	  << CSV(vars::pi0_leading_photon_conv_dist(i))
          << CSV(vars::pi0_leading_photon_conv_dist(j))
	  << CSV(vars::pi0_leading_photon_cosphi(i))
	  << CSV(vars::pi0_leading_photon_cosphi(j))
	  << CSV(vars::pi0_subleading_photon_energy(i))
	  << CSV(vars::pi0_subleading_photon_energy(j))
	  << CSV(vars::pi0_subleading_photon_conv_dist(i))
          << CSV(vars::pi0_subleading_photon_conv_dist(j))
	  << CSV(vars::pi0_subleading_photon_cosphi(i))
	  << CSV(vars::pi0_subleading_photon_cosphi(j))
	  << CSV(vars::pi0_costheta(i))
	  << CSV(vars::pi0_costheta(j))
	  << CSV(vars::pi0_mass(i))
	  << CSV(vars::pi0_mass(j))
	  << CSV(vars::pi0_momentum_mag(i))
	  << CSV(vars::pi0_momentum_mag(j))
	  << CSV(vars::pi0_beam_costheta(i))
	  << CSV(vars::pi0_beam_costheta(j))
	  << CSV(vars::transverse_momentum_mag(i))
	  << CSV(vars::transverse_momentum_mag(j))
	  << CSV(cuts::all_1mu0pi2gamma_cut(j))
	  << std::endl;
}


void write_reco(const caf::SRSpillProxy* sr, const caf::SRInteractionDLPProxy& j)
{
  output  << CSV(sr->hdr.run) << CSV(sr->hdr.evt) << CSV(sr->hdr.subrun)
	  << CSV(std::string(sr->hdr.sourceName))
	  << CSV(vars::id(j))
    //<< CSV(j.vertex[0])
    //<< CSV(j.flash_time)
    //<< CSV(vars::muon_momentum_mag(j))
    //<< CSV(vars::muon_costheta_z(j))
	  << CSV(vars::pi0_leading_photon_energy(j))
	  << CSV(vars::pi0_subleading_photon_energy(j))
    //<< CSV(vars::pi0_costheta(j))
    //<< CSV(vars::pi0_leading_photon_start_to_vertex(j))
    //<< CSV(vars::pi0_subleading_photon_start_to_vertex(j))
	  << CSV(vars::pi0_mass(j))
    //<< CSV(vars::pi0_momentum_mag(j))
    //<< CSV(vars::pi0_costheta_z(j))
	  << CSV(cuts::all_1mu0pi2gamma_cut(j))
	  << std::endl;
}


/**
 * Writes reconstructed variables (reco) for selected interactions.
 * @param sr is an SRSpillProxy that attaches to the StandardRecord of the
 * current spill.
 * @param j the reco interaction (selected).
 * @return None.
 */
const SpillMultiVar kDataInfo([](const caf::SRSpillProxy* sr)
{
  for(auto const & i : sr->dlp)
    {
      if(cuts::all_1mu0pi2gamma_cut(i))
	{
	  OUT(output,"DATA");
	  write_reco(sr, i);
	}
    }  
  return std::vector<double>{1};
});


const SpillMultiVar kSelectedVar([](const caf::SRSpillProxy* sr)
{
  /**
   * Loop over all selected interactions for purity study.
   */
  for(auto const & i : sr->dlp)
    {
      if(cuts::all_1mu0pi2gamma_cut(i))
	{
	  if(cuts::matched(i))
	    {
	      OUT(output, "SELECTED");
	      const auto & t = sr->dlp_true[i.match_ids[0]];
	      output << CSV(sr->hdr.run) << CSV(sr->hdr.evt) << CSV(sr->hdr.subrun) 
		     << CSV(std::string(sr->hdr.sourceName))
		     << CSV(vars::id(i))
		     << CSV(vars::category_topology(t))
		//<< CSV(vars::category_interaction_mode(t))
		//<< CSV(vars::muon_ke(i))
		//<< CSV(vars::muon_ke(t))
		//<< CSV(vars::pi0_subleading_photon_energy(i))
		//<< CSV(vars::pi0_leading_photon_energy(i))
		//<< CSV(vars::pi0_costheta(i))
		//<< CSV(vars::pi0_mass(i))
		     << std::endl;
	    } // Matched to true interaction
	} // Selection cuts
    }
  return std::vector<double>{1};
});

const SpillMultiVar kSignalVar([](const caf::SRSpillProxy* sr)
{
  /**
   * Loop over all truth interactions for efficiency study.
   */
  for(auto const & i : sr->dlp_true)
    {
      if(cuts::neutrino(i))
	{
	  int category(vars::category(i));
	  if(category == 0)
	    {
	      if(cuts::matched(i))
		{
		  OUT(output, "SIGNAL");
		  const auto & r = sr->dlp[i.match_ids[0]];
		  output << CSV(std::string(sr->hdr.sourceName))
		         << CSV(vars::id(i))
		         << CSV(vars::category_interaction_mode(i))
			 << CSV(i.nu_id)
			 << CSV(i.energy_init)
			 << CSV(vars::pi0_subleading_photon_energy(i))
			 << CSV(vars::pi0_leading_photon_energy(i))
		    //<< CSV(vars::pi0_costheta(i))
			 << CSV(vars::pi0_mass(i))
			 << CSV(cuts::all_1mu0pi2gamma_cut(r))
			 << std::endl;
		    
		} // Matched to reco interaction
	    } // 1mu1pi0 signal
	} // Neutrinos (not cosmics)
    } // True interactions

  return std::vector<double>{1};

});
			 

const SpillMultiVar kInfoVar([](const caf::SRSpillProxy* sr)
{
    /**
     * Loop over truth interactions for efficiency metrics and for signal-level
     * variables of interest.
    */
    for(auto const & i : sr->dlp_true)
    {
        if(cuts::neutrino(i))
        {
	  
	  // TEST
	  /*
	    cout << i.id << endl;
	    //for(const auto & j : cuts::count_primaries(i))
	    //  {
	    //cout << j << endl;
	    //}
	    //cout << cuts::count_primaries(i) << endl;
	    cout << cuts::fiducial_muon_containment_cut(i) << endl;
	    //cout << cuts::topological_1mu2gamma_cut(i) << endl;
	    cout << category << endl;
	    for(auto const  & p : i.particles)
	      {
		cout << p.pid << " " << p.pdg_code <<  " " << p.ke << " " <<  " " << p.is_contained << endl;
	      }
	  */

	  int category(vars::category(i));
	  if(category == 0)
            {
                if(cuts::matched(i))
                {
                    OUT(output, "SIGNAL");
                    const auto & r = sr->dlp[i.match_ids[0]];
                    write_pair(sr, i, r);
                }
            }
        }
    }

    /**
     * Loop over reconstructed interactions for purity metrics and for
     * reconstructed variables of interest.
    */
    for(auto const & i : sr->dlp)
    {
      if(cuts::all_1mu0pi2gamma_cut(i))
        {
            if(cuts::matched(i))
            {
                const auto & t = sr->dlp_true[i.match_ids[0]];
                OUT(output, "SELECTED");
                write_pair(sr, t, i);
            }
        }
    }

    return std::vector<double>{1};
});

#endif
