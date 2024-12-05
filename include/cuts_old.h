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
     * Apply a 1mu2gamma topological cut.  The interaction must have a topology
     * matching 1mu2gamma as defined by the conditions in the count_primaries() function.
     * @tparm T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction has a 1mu2gamma topology.
     */
    template<class T>
      bool topological_1mu2gamma_cut(const T & interaction)
      {
	std::vector<uint32_t> c(count_primaries(interaction));
        return c[2] == 1 && c[0] == 2;
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
     * Apply a containment cut to the primary muon.
     * All points in muon must be at least 5 cm from
     * the detector boundaries.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the muon is contained.
     */
    template<class T>
      bool muon_containment_cut(const T & interaction)
      {
	bool primary_muon_contained = false;
	std::vector<uint32_t> counts(count_primaries(interaction));
	if(counts[2] >= 1)
	  {
	    double leading_ke(0);
	    size_t index(0);
	    for(size_t i(0); i < interaction.particles.size(); ++i)
	      {
		const auto & p = interaction.particles[i];
		double energy = p.csda_ke;
		if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
			   energy = p.ke;
		if(p.pid == 2 && p.is_primary && energy > leading_ke)
		  {
		    leading_ke = energy;
		    index = i;
		  }
	      }

	    if(interaction.particles[index].is_contained)
	      {
		primary_muon_contained = true;
	      }
	  }
	
	return primary_muon_contained;
      }

    /**
     * Define the true 1mu0pi1pi0 interaction topology.
     * @tparam T the type of interactino (true or reco).
     * @param interaction to select on.
     * @return true if interaction is a 1mu0pi1pi0 neutrino interaction.
     */
    template<class T>
      bool true_1mu0pi1pi0_topology_cut(const T & interaction)
      {
        int primary_muon_count = 0;
	int primary_pion_count = 0;
        int primary_pi0_count = 0;
	int nonprimary_pi0_count = 0;
        unordered_map< int, vector<double> > primary_pi0_map;
	unordered_map< int, vector<double> > nonprimary_pi0_map;
        // Loop over particles
        for(auto & p : interaction.particles)
          {
            // Primary muon
            if(p.pid == 2 && p.is_primary && p.ke >= MIN_MUON_ENERGY)
              {
                ++primary_muon_count;
              }
	    // Primary pion
	    if(p.pid == 3 && p.is_primary && p.ke >= MIN_PION_ENERGY)
	      {
		++primary_pion_count;
	      }
            // Group primary pi0 photons
            if(p.pdg_code == 22 && p.is_primary && p.parent_pdg_code == 111)
              {
                primary_pi0_map[p.parent_track_id].push_back(p.ke);
              }
	    // Group nonprimary pi0 photons
	    if(p.pdg_code == 22 && !p.is_primary && p.parent_pdg_code == 111)
	      {
		nonprimary_pi0_map[p.parent_track_id].push_back(p.ke);
	      }
          }

        // Loop over primary pi0s 
        for (auto const& pi0 : primary_pi0_map)
          {
            int num_primary_photon_daughters = 0;
            for(auto & e : pi0.second)
              {
                if(e >= MIN_PHOTON_ENERGY) num_primary_photon_daughters++;
              }

            if(num_primary_photon_daughters == 2)
              {
                primary_pi0_count++;
              }
          }

	// Loop over nonprimary pi0s
	for(auto const& pi0 : nonprimary_pi0_map)
	  {
	    int num_nonprimary_photon_daughters = 0;
	    for(auto e : pi0.second)
	      {
		num_nonprimary_photon_daughters++;
	      }
	    nonprimary_pi0_count++;
	  }

	if (primary_muon_count == 1 && primary_pion_count == 0 && primary_pi0_count == 1 && nonprimary_pi0_count == 0 && interaction.current_type == 0)
          {
            return true;
          }
        else
          {
            return false;
          }
      }

    
    /**
     * Apply a flash time cut. The interaction must be matched to an in-time
     * flash. The in-time definition is valid for BNB simulation.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction has been matched to an in-time flash.
     */
    template<class T>
        bool flash_cut_bnb(const T & interaction)
        {
            if(!valid_flashmatch(interaction))
                return false;
            else
                return (interaction.flash_time >= 0) && (interaction.flash_time <= 1.6);
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
     * Apply a flash time cut. The interaction must be matched to an in-time
     * flash. The in-time definition is valid for BNB data.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction has been matched to an in-time flash.
     */
    template<class T>
        bool flash_cut_data(const T & interaction)
        {
            if(!valid_flashmatch(interaction))
                return false;
            else
                return (interaction.flash_time >= -0.5) && (interaction.flash_time <= 1.4);
        }

    /**
     * Apply a fiducial and containment cut (logical "and" of both).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction passes the fiducial and containment cut.
     */
    template<class T>
        bool fiducial_containment_cut(const T & interaction) {return fiducial_cut<T>(interaction) && containment_cut<T>(interaction);}

    /**
     * Apply a fiducial and muon containment cut (logical "and" of both).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction passes the fiducial and muon containment cut.
     */
    template<class T>
      bool fiducial_muon_containment_cut(const T & interaction) {return fiducial_cut<T>(interaction) && muon_containment_cut(interaction);}

    /**
     *
     *
     *
     */
    template<class T> bool fiducial_track_containment_cut(const T & interaction) {return fiducial_cut<T>(interaction) && track_containment_cut(interaction);}

    /**
     * Apply a fiducial, containment, and topological (1mu2gamma) cut (logical
     * "and" of each).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction passes the fiducial, containment, and
     * topological cut.
     */
    template<class T>
      bool fiducial_muon_containment_topological_1mu2gamma_cut(const T & interaction) {return fiducial_cut<T>(interaction) && muon_containment_cut<T>(interaction) && topological_1mu2gamma_cut<T>(interaction);}

    /**
     * Apply a fiducial, containment, topological (1mu2gamma), and flash time cut
     * (logical "and" of each).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction passes the fiducial, containment,
     * topological, and flash time cut.
     */
    template<class T>
        bool all_1mu2gamma_bnb_cut(const T & interaction) {return topological_1mu2gamma_cut<T>(interaction) && fiducial_cut<T>(interaction) && flash_cut_bnb<T>(interaction) && muon_containment_cut<T>(interaction);}

    /**
     * Apply fiducial, containment, topological (1mu2gamma), and flash time(NuMI) cut
     * (logical "and" of each).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction passes the fiducial, containment,
     * topological, and flash time cut.
     */
    template<class T> bool all_1mu2gamma_numi_cut(const T & interaction)
      {
	return topological_1mu2gamma_cut<T>(interaction) && fiducial_cut<T>(interaction) && muon_containment_cut<T>(interaction) && flash_cut_numi<T>(interaction);
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
	return topological_1mu0pi2gamma_cut<T>(interaction) && fiducial_track_containment_cut<T>(interaction) && flash_cut_numi<T>(interaction);
      }

    /**
     * Apply a fiducial, containment, topological (1mu2gamma), and flash time cut
     * (logical "and" of each).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction passes the fiducial, containment,
     * topological, and flash time cut.
     */
    template<class T>
        bool all_1mu2gamma_data_cut(const T & interaction) { return topological_1mu2gamma_cut<T>(interaction) && fiducial_cut<T>(interaction) && flash_cut_data<T>(interaction) && muon_containment_cut<T>(interaction); }
    
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
    template<class T>
        bool signal_1mu0pi1pi0(const T & interaction) {return true_1mu0pi1pi0_topology_cut<T>(interaction) && neutrino(interaction);}
    
    /**
     * Define the true "other neutrino" interaction classification.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction is an "other neutrino" interaction.
     */
    template<class T>
      bool other_nu_1mu0pi1pi0(const T & interaction) {return !true_1mu0pi1pi0_topology_cut<T>(interaction) && neutrino(interaction);}
    
}
#endif
