/**
 * @file cuts.h
 * @brief Header file for definitions of selection cuts.
 * @author justin.mueller@colostate.edu 
*/
#ifndef CUTS_H
#define CUTS_H

#define MIN_PHOTON_ENERGY 50.0
#define MIN_MUON_ENERGY 250.345
#define MIN_PION_ENERGY 25

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
};

#include <functional>
#include <vector>
#include <TVector3.h>
#include <string>
#include <sstream>
#include <numeric>
#include <iostream>

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
	    if(p.pid == 0 && energy >= MIN_PHOTON_ENERGY)
 	      {
		passes = true;
	      }
	    // Electrons have no energy requirement
	    if(p.pid == 1)
	      {
		passes = true;
	      }
	    // Muons have minimum energy requirement (250.345 MeV)
	    if(p.pid == 2 && energy >= MIN_MUON_ENERGY)
	      {
		passes = true;
	      }
	    // Pions have minimum energy requirement (25 MeV)
	    if(p.pid == 3 && energy >= MIN_PION_ENERGY)
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
     * Apply pi0 mass cut (to reject eta mesons).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if interaction passes pi0 mass cut.
     */
    template<class T> bool pi0_mass_cut(const T & interaction)
      {
	bool mass_pass = false;
	if(topological_1mu2gamma_cut(interaction))
	  {
	    TVector3 vertex(interaction.vertex[0], interaction.vertex[1], interaction.vertex[2]);
	    vector<double> photon_energies;
	    vector<TVector3> photon_dirs;

	    // Loop over particles
	    for(auto & p : interaction.particles)
	      {
		// primary photons
		if(p.is_primary && p.pid == 0 && p.calo_ke > MIN_PHOTON_ENERGY)
		  {
		    // photon energies
		    photon_energies.push_back(p.calo_ke);

		    // photon directions
		    TVector3 ph_dir;
		    ph_dir.SetX(p.start_point[0] - vertex[0]);
		    ph_dir.SetY(p.start_point[1] - vertex[1]);
		    ph_dir.SetZ(p.start_point[2] - vertex[2]);
		    photon_dirs.push_back(ph_dir.Unit());
		  } // end primary photon loop
	      } // end particle loop

	    double cos_theta = photon_dirs[0].Dot(photon_dirs[1]);
	    double pi0_mass = sqrt(2*photon_energies[0] * photon_energies[1] * (1-cos_theta));

	    if(pi0_mass < 400) mass_pass = true;
	  }

	return mass_pass;
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
	unordered_map< int, vector<double> > primary_pi0_map;
	unordered_map< int, vector<double> > nonprimary_pi0_map;
	
	// Loop over particles
	for(auto & p : interaction.particles)
	  {
	    // Primary muons
	    if(p.pid == 2 && p.is_primary)
	      {
		primary_muon_count++;
		if(p.ke >= MIN_MUON_ENERGY) primary_muon_count_thresh++;
	      }
	    // Primary pions
	    if(p.pid == 3 && p.is_primary)
	      {
		primary_pion_count++;
		if(p.ke >= MIN_PION_ENERGY) primary_pion_count_thresh++;
	      }
	    // Primary photons
	    if(p.pdg_code == 22 && p.is_primary && p.parent_pdg_code == 111)
              {
                primary_pi0_map[p.parent_track_id].push_back(p.ke);
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
            for(auto & e : pi0.second)
              {
		num_primary_photon_daughters++;
                if(e >= MIN_PHOTON_ENERGY) num_primary_photon_daughters_thresh++;
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
	
	return s;
      }
    
    /**
     * Apply a flash tie cut.  The interaction must be matched ot an in-time
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
     * Apply fiducial, track containment, topological (1mu + 2gamma + 0pi), and flahs time (NuMI) cut
     * (logical "and" of each).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction passes above cuts.
     */
    template<class T> bool all_1mu0pi2gamma_cut(const T & interaction)
      {
	return topological_1mu0pi2gamma_cut<T>(interaction) && fiducial_cut<T>(interaction) && track_containment_cut<T>(interaction) && flash_cut_numi<T>(interaction);
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
	return s.num_primary_muons_thresh == 1 && s.num_primary_pions_thresh == 0 && s.num_primary_pi0s_thresh == 1 && s.num_nonprimary_pi0s == 0 && s.is_cc && s.is_neutrino;
      }
    
    template<class T> bool other_nu_1mu0pi1pi0(const T & interaction)
      {
	truth_inter s = true_interaction_info(interaction);
	return !(s.num_primary_muons_thresh == 1 && s.num_primary_pions_thresh == 0 && s.num_primary_pi0s_thresh == 1 && s.num_nonprimary_pi0s == 0 && s.is_cc) && s.is_neutrino;
      }
}
#endif
