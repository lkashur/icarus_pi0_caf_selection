/**
 * @file cuts.h
 * @brief Header file for definitions of selection cuts.
 * @author justin.mueller@colostate.edu 
 */

#include <functional>
#include <vector>
#include <TVector3.h>
#include <string>
#include <sstream>
#include <numeric>
#include <iostream>

#ifndef CUTS_H
#define CUTS_H

#define MIN_PHOTON_ENERGY_BNB 25.0
#define MIN_MUON_ENERGY_BNB 143.425
#define MIN_PION_ENERGY_BNB 25.0

#define MIN_PHOTON_ENERGY_NUMI 50.0
#define MIN_MUON_ENERGY_NUMI 250.345
#define MIN_PION_ENERGY_NUMI 25.0

struct truth_inter {
  int num_primary_muons;
  int num_primary_muons_thresh;
  int num_primary_pions;
  int num_primary_pions_thresh;
  int num_primary_pi0s;
  int num_primary_pi0s_thresh;
  int num_nonprimary_pi0s;
  bool is_fiducial;
  bool has_contained_tracks;
  bool is_neutrino;
  bool is_cc;
  double muon_momentum_mag;
  double muon_beam_costheta;
  double pi0_leading_photon_energy;
  double pi0_leading_photon_conv_dist;
  double pi0_subleading_photon_energy;
  double pi0_subleading_photon_conv_dist;
  double pi0_costheta;
  double pi0_mass;
  double pi0_momentum_mag;
  double pi0_beam_costheta;
};

struct reco_pi0 {
  double muon_momentum_mag;
  double muon_beam_costheta;
  double leading_photon_energy;
  double leading_photon_cosphi;
  double leading_photon_conv_dist;
  double subleading_photon_energy;
  size_t subleading_photon_cosphi;
  double subleading_photon_conv_dist;
  double costheta;
  double mass;
  double momentum_mag;
  double beam_costheta;
};

namespace cuts
{
  /**
   * Apply a cut on whether a match exists.
   * @tparam T the object type (true or reco, interaction or particle).
   * @param obj the object to select on.
   * @return true if the object is matched.
   */
    template<class T>
      bool matched(const T & obj) { return obj.match_ids.size() > 0; }

    /**
     * Apply a cut on the validity of the flash match.
     * @tparam T the type of interaction (true or reco).
     * @param interaction on which to place the flash validity cut.
     * @return true if the interaction is flash matched and the time is valid.
     */
    template<class T>
      bool valid_flashmatch(const T & interaction) { return !std::isnan(interaction.flash_time) && interaction.is_flash_matched == 1; }

    /**
     * Check if the particle meets final state signal requirements.
     * Particles must be primary and have an energy above threshold.
     * Muons must have a length of at least 50 cm (143.425 MeV), protons
     * must have an energy above 50 MeV, and all other particles must have
     * an energy above 25 MeV.
     * @tparam T the type of particle (true or reco).
     * @param particle to check.
     * @return true if the particle is a final state signal particle.
     */
    template<class T> bool final_state_signal(const T & p)
      {
	bool passes(false);
	if(p.is_primary)
	  {
	    double energy(p.pid > 1 ? p.csda_ke : p.calo_ke);
	    if constexpr (std::is_same_v<T, caf::SRParticleTruthDLPProxy>)
			   energy = p.ke;
	        
	    // Photons have minimum energy requirement (50 MeV)
	    if(p.pid == 0 && energy >= MIN_PHOTON_ENERGY_BNB)
	      {
		passes = true;
	      }
	    // Electrons have no energy requirement
	    if(p.pid == 1)
	      {
		passes = true;
	      }
	    // Muons have minimum energy requirement (250.345 MeV)
	    if(p.pid == 2 && energy >= MIN_MUON_ENERGY_BNB)
	      {
		passes = true;
	      }
	    // Pions have minimum energy requirement (25 MeV)
	    if(p.pid == 3 && energy >= MIN_PION_ENERGY_BNB)
	      {
		passes = true;
	      }
	    // Protons have no minimum energy requirement
	    if(p.pid == 4)
	      {
		passes = true;
	      }
	    // Kaons have no minimum energy requirement
	    if(p.pid == 5)
	      {
		passes = true;
	      }
	  }
	return passes;
      }

    /**
     * Count the primaries of the interaction with cuts applied to each particle.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to find the topology of.
     * @return the count of primaries of each particle type within the
     * interaction.
     */
    template<class T>
      std::vector<uint32_t> count_primaries(const T & interaction)
      {
	std::vector<uint32_t> counts(5, 0);
	for(auto &p : interaction.particles)
	  {
	    if(final_state_signal(p))
	      ++counts[p.pid];
	  }
	return counts;
      }

    /**
     * Find the topology of the interaction with cuts applied to each particle.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to find the topology of.
     * @return the topology of the interaction as a string (e.g 0ph0e1mu0pi1p).
     */
    template<class T>
      std::string topology(const T & interaction)
      {
	std::vector<uint32_t> counts(count_primaries(interaction));
	std::stringstream ss;
            ss  << counts[0] << "ph"
                << counts[1] << "e"
                << counts[2] << "mu"
                << counts[3] << "pi"
                << counts[4] << "p";
            return ss.str();
      }

    /**
     * Apply selection for 1mu + 2gamma + 0pi + X topology.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction has 1mu + 1pi0 + 0pi + X topology.
     */
    template <class T> bool topological_1mu0pi2gamma_cut(const T & interaction)
      {
	std::vector<uint32_t> c(count_primaries(interaction));
	return c[2] == 1 && c[0] == 2 && c[3] == 0;
      }

    /**
     * Obtain all relevant reco information about pi0 candidate.
     * @tparam T the type of interaction (reco).
     * @param interaction to select on.
     * @return stuff.
     */
    template<class T> reco_pi0 reco_pi0_info(const T & interaction)
      {
	reco_pi0 s;

	// Loop over particles
	size_t leading_muon_index(0);
	size_t leading_photon_index(0);
	size_t subleading_photon_index(0);
	double max_csda_ke(-99999);
	double max_calo_ke0(-99999);
	double max_calo_ke1(-99999);
	for(size_t i(0); i < interaction.particles.size(); ++i)
	  {
	    const auto & p = interaction.particles[i];

	    // Primary muons
	    if(p.is_primary && p.pid == 2 && p.csda_ke >= MIN_MUON_ENERGY_BNB)
	      {
		if(p.csda_ke > max_csda_ke)
		  {
		    max_csda_ke = p.csda_ke;
		    leading_muon_index = i;
		  }
	      }

	    // Primary photons
	    if(p.is_primary && p.pid == 0 && p.calo_ke >= MIN_PHOTON_ENERGY_BNB)
	      {
		if(p.calo_ke > max_calo_ke0)
		  {
		    max_calo_ke0 = p.calo_ke;
		    leading_photon_index = i;
		  }
	      } // end photon loop
	  } // end particle loop

	for(size_t i(0); i < interaction.particles.size(); ++i)
	  {
	    const auto & p = interaction.particles[i];
	    // Priary photons
	    if(p.is_primary && p.pid == 0 && p.calo_ke >= MIN_PHOTON_ENERGY_BNB)
	      {
		if(p.calo_ke > max_calo_ke1 && p.calo_ke < max_calo_ke0)
		  {
		    max_calo_ke1 = p.calo_ke;
		    subleading_photon_index = i;
		  }
	      }
	  }

	// Get info about pi0 from selected photons
	TVector3 beamdir(0, 0, 1); // BNB
	//TVector3 beamdir(0.39431672, 0.04210058, 0.91800973); // NuMI
	TVector3 vertex(interaction.vertex[0], interaction.vertex[1], interaction.vertex[2]);
	
	TVector3 muon_momentum;
	muon_momentum.SetX(interaction.particles[leading_muon_index].momentum[0]);
	muon_momentum.SetY(interaction.particles[leading_muon_index].momentum[1]);
	muon_momentum.SetZ(interaction.particles[leading_muon_index].momentum[2]);
	double muon_momentum_mag = muon_momentum.Mag();
	double muon_beam_costheta = muon_momentum.Unit().Dot(beamdir);
	
	double leading_photon_energy = interaction.particles[leading_photon_index].calo_ke;
	TVector3 leading_photon_start_point;
	leading_photon_start_point.SetX(interaction.particles[leading_photon_index].start_point[0]);
	leading_photon_start_point.SetY(interaction.particles[leading_photon_index].start_point[1]);
	leading_photon_start_point.SetZ(interaction.particles[leading_photon_index].start_point[2]);
	double leading_photon_conv_dist = (vertex - leading_photon_start_point).Mag();
	TVector3 leading_photon_dir;
	leading_photon_dir.SetX(leading_photon_start_point[0] - vertex[0]);
	leading_photon_dir.SetY(leading_photon_start_point[1] - vertex[1]);
	leading_photon_dir.SetZ(leading_photon_start_point[2] - vertex[2]);
	leading_photon_dir = leading_photon_dir.Unit();
	TVector3 leading_photon_start_dir;
	leading_photon_start_dir.SetX(interaction.particles[leading_photon_index].start_dir[0]);
	leading_photon_start_dir.SetY(interaction.particles[leading_photon_index].start_dir[1]);
	leading_photon_start_dir.SetZ(interaction.particles[leading_photon_index].start_dir[2]);
	leading_photon_start_dir = leading_photon_start_dir.Unit();
	double leading_photon_cosphi = leading_photon_dir.Dot(leading_photon_start_dir);
	TVector3 leading_photon_momentum = leading_photon_energy * leading_photon_dir;

	double subleading_photon_energy = interaction.particles[subleading_photon_index].calo_ke;
        TVector3 subleading_photon_start_point;
        subleading_photon_start_point.SetX(interaction.particles[subleading_photon_index].start_point[0]);
        subleading_photon_start_point.SetY(interaction.particles[subleading_photon_index].start_point[1]);
        subleading_photon_start_point.SetZ(interaction.particles[subleading_photon_index].start_point[2]);
	double subleading_photon_conv_dist = (vertex - subleading_photon_start_point).Mag();
	TVector3 subleading_photon_dir;
	subleading_photon_dir.SetX(subleading_photon_start_point[0] - vertex[0]);
	subleading_photon_dir.SetY(subleading_photon_start_point[1] - vertex[1]);
	subleading_photon_dir.SetZ(subleading_photon_start_point[2] - vertex[2]);
	subleading_photon_dir = subleading_photon_dir.Unit();
	TVector3 subleading_photon_start_dir;
	subleading_photon_start_dir.SetX(interaction.particles[subleading_photon_index].start_dir[0]);
	subleading_photon_start_dir.SetY(interaction.particles[subleading_photon_index].start_dir[1]);
	subleading_photon_start_dir.SetZ(interaction.particles[subleading_photon_index].start_dir[2]);
	subleading_photon_start_dir = subleading_photon_start_dir.Unit();
	double subleading_photon_cosphi = subleading_photon_dir.Dot(subleading_photon_start_dir);
	TVector3 subleading_photon_momentum = subleading_photon_energy * subleading_photon_dir;

	double cos_opening_angle = leading_photon_dir.Dot(subleading_photon_dir);
	double mass = sqrt(2*leading_photon_energy*subleading_photon_energy*(1-cos_opening_angle));
	TVector3 momentum = leading_photon_momentum + subleading_photon_momentum;
	double momentum_mag = momentum.Mag();
	double beam_costheta = momentum.Unit().Dot(beamdir);

	// Output
	s.muon_momentum_mag = muon_momentum_mag;
	s.muon_beam_costheta = muon_beam_costheta;
	s.leading_photon_energy = leading_photon_energy;
	s.leading_photon_cosphi = leading_photon_cosphi;
	s.leading_photon_conv_dist = leading_photon_conv_dist;
	s.subleading_photon_energy = subleading_photon_energy;
	s.subleading_photon_cosphi = subleading_photon_cosphi;
	s.subleading_photon_conv_dist = subleading_photon_conv_dist;
	s.costheta = cos_opening_angle;
	s.mass = mass;
	s.momentum_mag = momentum_mag;
	s.beam_costheta = beam_costheta;

	return s;
      }

    /**
     * Apply a fiducial volume cut. Interaction vertex must be within 25 cm of
     * x and y detector faces, 50 cm of downstream (+) z face, and 30 cm of
     * upstream (-) z face.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the vertex is in the fiducial volume.
     */
    template<class T>
      bool fiducial_cut(const T & interaction)
      {
	return interaction.is_fiducial && !(interaction.vertex[0] > 210.215 && interaction.vertex[1] > 60 && (interaction.vertex[2] > 290 && interaction.vertex[2] < 390));
      }
    
    /**
     * Apply a containment volume cut. All points within the interaction must be
     * at least 5 cm from the detector boundaries.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction is contained.
     */
    template<class T>
      bool containment_cut(const T & interaction) { return interaction.is_contained; }

    /**
     * Apply a containment cut to all primary tracks in interaction.
     * All points within tracks must be at least 5 cm from detector boundaries.
     * @tparam T the type o interaction (true or reco).
     * @param interaction to select on.
     * @return true if interactin is contained.
     */
    template<class T> bool track_containment_cut(const T & interaction)
      {
	bool passes(true);
	for(auto & p : interaction.particles)
	  {
	    if(p.is_primary && p.pid > 1 && !p.is_contained)
	      {
		passes = false;
	      }
	  }
	return passes;
      }   

    /**
     * Obtain all relevant truth information about interaction.
     * @tparam T the type of interaction (true).
     * @param interaction to select on.
     * @return ture if interaction is contained.
     */
    template<class T> truth_inter true_interaction_info(const T & interaction)
      {
	truth_inter s;

	int primary_muon_count(0);
	int primary_muon_count_thresh(0);
	int primary_pion_count(0);
	int primary_pion_count_thresh(0);
	int primary_pi0_count(0);
	int primary_pi0_count_thresh(0);
	int nonprimary_pi0_count(0);
	bool is_neutrino(false);
	bool is_cc(false);
	//unordered_map< int, vector<double> > primary_pi0_map;
	unordered_map<int, vector<pair<size_t, double>> > primary_pi0_map;
	unordered_map< int, vector<double> > nonprimary_pi0_map;
	
	TVector3 vertex(interaction.vertex[0], interaction.vertex[1], interaction.vertex[2]);

	// Loop over particles
	size_t leading_muon_index(0);
	double max_csda_ke(-99999);
	for(size_t i(0); i < interaction.particles.size(); ++i)
	  {
	    const auto & p = interaction.particles[i];

	    // Primary muons
	    if(p.pid == 2 && p.is_primary)
	      {
		primary_muon_count++;
		if(p.ke >= MIN_MUON_ENERGY_BNB) primary_muon_count_thresh++;
		
		if(p.ke > max_csda_ke)
		  {
		    max_csda_ke = p.ke;
		    leading_muon_index = i;
		  }
		
	      }
	    // Primary pions
	    if(p.pid == 3 && p.is_primary)
	      {
		primary_pion_count++;
		if(p.ke >= MIN_PION_ENERGY_BNB) primary_pion_count_thresh++;
	      }
	    // Primary photons
	    if(p.pdg_code == 22 && p.is_primary && p.parent_pdg_code == 111)
              {
                primary_pi0_map[p.parent_track_id].push_back({i, p.ke});
              }
	    // Nonprimary photons
	    if(p.pdg_code == 22 && !p.is_primary && p.parent_pdg_code == 111)
              {
                nonprimary_pi0_map[p.parent_track_id].push_back(p.ke);
              }
	  } // end particle loop

	// Loop over primary pi0s
	for(auto const & pi0 : primary_pi0_map)
          {
            int num_primary_photon_daughters(0);
	    int num_primary_photon_daughters_thresh(0);
            for(auto & pair : pi0.second)
              {
		num_primary_photon_daughters++;
                if(pair.second >= MIN_PHOTON_ENERGY_BNB) num_primary_photon_daughters_thresh++;
              }
            if(num_primary_photon_daughters == 2) primary_pi0_count++;
	    if(num_primary_photon_daughters_thresh == 2) primary_pi0_count_thresh++;
          }
	
	// Nonprimary pi0s
	nonprimary_pi0_count = nonprimary_pi0_map.size();

	// Output
	s.num_primary_muons = primary_muon_count;
	s.num_primary_muons_thresh = primary_muon_count_thresh;
	s.num_primary_pions = primary_pion_count;
	s.num_primary_pions_thresh = primary_pion_count_thresh;
	s.num_primary_pi0s = primary_pi0_count;
	s.num_primary_pi0s_thresh = primary_pi0_count_thresh;
	s.num_nonprimary_pi0s = nonprimary_pi0_count;
	s.is_fiducial = fiducial_cut<T>(interaction);
	s.has_contained_tracks = track_containment_cut<T>(interaction);
	if(interaction.nu_id > -1) is_neutrino = true;
	s.is_neutrino = is_neutrino;
	if(interaction.current_type == 0) is_cc = true;
	s.is_cc = is_cc;
	
	// Add true pi0 info, if it exists
	TVector3 beamdir(0, 0, 1); // BNB
        //TVector3 beamdir(0.39431672, 0.04210058, 0.91800973); // NuMI

	double muon_momentum_mag;
	double muon_beam_costheta;
	double pi0_leading_photon_energy;
	TVector3 pi0_leading_photon_dir;
	double pi0_leading_photon_conv_dist;

	double pi0_subleading_photon_energy;
	TVector3 pi0_subleading_photon_dir;
	double pi0_subleading_photon_conv_dist;

	double pi0_costheta;
	double pi0_mass;
	double pi0_momentum_mag;
	double pi0_beam_costheta;

	vector<size_t> pi0_photon_indices;
	if(primary_pi0_count_thresh == 1)
	  {
	    for(auto const & pi0 : primary_pi0_map)
	      {
		for(auto pair : pi0.second)
		  {
		    if(pair.second >= MIN_PHOTON_ENERGY_BNB)
		      {
			pi0_photon_indices.push_back(pair.first);
		      }
		  }
	      }

	    const auto & muon = interaction.particles[leading_muon_index];
	    TVector3 muon_momentum;
	    muon_momentum.SetX(muon.momentum[0]);
	    muon_momentum.SetY(muon.momentum[1]);
	    muon_momentum.SetZ(muon.momentum[2]);
	    double muon_momentum_mag = muon_momentum.Mag();
	    double muon_beam_costheta = muon_momentum.Unit().Dot(beamdir);
	    

	    const auto & pi0_photon0 = interaction.particles[pi0_photon_indices[0]];
	    const auto & pi0_photon1 = interaction.particles[pi0_photon_indices[1]];
	    size_t leading_photon_index;
	    size_t subleading_photon_index;
	    if(pi0_photon0.ke > pi0_photon1.ke)
	      {
		leading_photon_index = pi0_photon_indices[0];
		subleading_photon_index = pi0_photon_indices[1];
	      }
	    else
	      {
		leading_photon_index = pi0_photon_indices[1];
		subleading_photon_index = pi0_photon_indices[0];
	      }
	    const auto & pi0_leading_photon = interaction.particles[leading_photon_index];
	    const auto & pi0_subleading_photon = interaction.particles[subleading_photon_index];
	        
	    pi0_leading_photon_energy = pi0_leading_photon.ke;
	    TVector3 pi0_leading_photon_start_point(pi0_leading_photon.start_point[0], pi0_leading_photon.start_point[1], pi0_leading_photon.start_point[2]);
	    TVector3 pi0_leading_photon_dir(pi0_leading_photon.momentum[0], pi0_leading_photon.momentum[1], pi0_leading_photon.momentum[2]);
	    pi0_leading_photon_dir = pi0_leading_photon_dir.Unit();
	    pi0_leading_photon_conv_dist = (vertex - pi0_leading_photon_start_point).Mag();
	    TVector3 pi0_leading_photon_momentum(pi0_leading_photon.momentum[0], pi0_leading_photon.momentum[1], pi0_leading_photon.momentum[2]);

	    pi0_subleading_photon_energy = pi0_subleading_photon.ke;
	    TVector3 pi0_subleading_photon_start_point(pi0_subleading_photon.start_point[0], pi0_subleading_photon.start_point[1], pi0_subleading_photon.start_point[2]);
	    TVector3 pi0_subleading_photon_dir(pi0_subleading_photon.momentum[0], pi0_subleading_photon.momentum[1], pi0_subleading_photon.momentum[2]);
            pi0_subleading_photon_dir = pi0_subleading_photon_dir.Unit();
	    pi0_subleading_photon_conv_dist = (vertex - pi0_subleading_photon_start_point).Mag();
	    TVector3 pi0_subleading_photon_momentum(pi0_subleading_photon.momentum[0], pi0_subleading_photon.momentum[1], pi0_subleading_photon.momentum[2]);
	    
	    TVector3 pi0_momentum = pi0_leading_photon_momentum + pi0_subleading_photon_momentum;
	    pi0_beam_costheta = pi0_momentum.Unit().Dot(beamdir);

	    pi0_costheta = pi0_leading_photon_dir.Dot(pi0_subleading_photon_dir);
	    pi0_mass = sqrt(2*pi0_leading_photon_energy*pi0_subleading_photon_energy*(1-pi0_costheta));

	    s.muon_momentum_mag = muon_momentum_mag;
	    s.muon_beam_costheta = muon_beam_costheta;
	    s.pi0_leading_photon_energy = pi0_leading_photon_energy;
	    s.pi0_leading_photon_conv_dist = pi0_leading_photon_conv_dist;
	    s.pi0_subleading_photon_energy = pi0_subleading_photon_energy;
	    s.pi0_subleading_photon_conv_dist = pi0_subleading_photon_conv_dist;
	    s.pi0_costheta = pi0_costheta;
	    s.pi0_mass = pi0_mass;
	    s.pi0_momentum_mag = pi0_momentum.Mag();
	    s.pi0_beam_costheta = pi0_beam_costheta;
	  }
	
	return s;
      }
    
    /**
     * Apply a flash time cut.  The interaction must be matched to an in-time
     * flash.  The in-time definition is valid for BNB simulation.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if interaction has been matched to an in-time flash.
     */
    template<class T> bool flash_cut_bnb(const T & interaction)
      {
	if(!valid_flashmatch(interaction))
          {
            return false;
          }
        else
	  {
	    return (interaction.flash_time >= 0) && (interaction.flash_time <= 1.6);
	  }
      }

    /**
     * Apply a flash time cut.  The interaction must be matched to an in-time
     * flash.  The in-time definition is valid for NuMI simulation.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction has been matched to an in-time flash.
     */
    template<class T> bool flash_cut_numi(const T & interaction)
      {
	if(!valid_flashmatch(interaction))
	  {
	    return false;
	  }
	else
	  {
	    return (interaction.flash_time >= 0) && (interaction.flash_time <= 9.6);
	  }
      }
    
    /**
     * Apply pi0 mass cut (eta meson rejection).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction passes above cuts.
     */
    template<class T> bool pi0_mass_cut(const T & interaction)
      {
	reco_pi0 s = reco_pi0_info(interaction);
	return s.mass < 400;
      }

    /**
     * Apply fiducial, track containment, topological (1mu + 2gamma + 0pi), and flahs time (NuMI) cut
     * (logical "and" of each).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction passes above cuts.
     */
    template<class T> bool all_1mu0pi2gamma_cut(const T & interaction)
      {
	return topological_1mu0pi2gamma_cut<T>(interaction) && fiducial_cut<T>(interaction) && track_containment_cut<T>(interaction) && flash_cut_bnb<T>(interaction) && pi0_mass_cut<T>(interaction);
      }
    
    /**
     * Defined the true neutrino interaction classification.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction is a neutrino interaction.
     */
    template<class T>
      bool neutrino(const T & interaction) { return interaction.nu_id > -1; }

    /**
     * Define the true cosmic interaction classification.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction is a cosmic.
     */
    template<class T>
      bool cosmic(const T & interaction) { return interaction.nu_id == -1; }

    /**
     * Define the true neutrino interaction classification.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction is a neutrino interaction.
     */
    template<class T>
      bool matched_neutrino(const T & interaction) { return interaction.match_ids.size() > 0 && neutrino(interaction); }

    /**
     * Define the true neutrino interaction classification.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction is a neutrino interaction.
     */
    template<class T>
      bool matched_cosmic(const T & interaction) { return interaction.match_ids.size() > 0 && cosmic(interaction); }

    /**
     * Define the true 1mu1pi0 interaction classification.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction is a 1mu1p neutrino interaction.
     */
    template<class T> bool signal_1mu0pi1pi0(const T & interaction)
      {
	truth_inter s = true_interaction_info(interaction);
	//cout << s.num_primary_pi0s << s.num_primary_pi0s_thresh << endl;
	return s.num_primary_muons_thresh == 1 && s.num_primary_pions_thresh == 0 && s.num_primary_pi0s_thresh == 1 && s.is_cc && s.is_neutrino;
      }
    
    template<class T> bool other_nu_1mu0pi1pi0(const T & interaction)
      {
	truth_inter s = true_interaction_info(interaction);
	return !(s.num_primary_muons_thresh == 1 && s.num_primary_pions_thresh == 0 && s.num_primary_pi0s_thresh == 1 && s.num_nonprimary_pi0s == 0 && s.is_cc) && s.is_neutrino;
      }
}
#endif
