/**
 * @file variables.h
 * @brief Header file for definitions of selection variables.
 * @author justin.mueller@colostate.edu
*/
#ifndef VARIABLES_H
#define VARIABLES_H

#define ELECTRON_MASS 0.5109989461
#define MUON_MASS 105.6583745
#define PION_MASS 139.57039
#define PROTON_MASS 938.2720813

#include <algorithm>
#include <iostream>
#include <TVector3.h>
#include <string>
#include <iostream>

namespace vars
{

    /**
     * Variable for counting interactions/particles.
     * @tparam T the type of object (true or reco, interaction or particle).
     * @param interaction/particle to apply the variable on.
     * @return 1.0 (always).
     */
    template<class T>
        double count(const T & obj) { return 1.0; }
   
    /**
     * Variable for id (unique identifier for the object).
     * @tparam T the type of object (true or reco, interaction or particle).
     * @param interaction/particle to apply the variable on.
     * @return the id of the interaction/particle.
    */
    template<class T>
        double id(const T & obj) { return obj.id; }

    /**
     * Variable for enumerating interaction categories.  This is a basic
     * categorization using only signal, neutrino background, and cosmic
     * background as the three categories.
     * 0: 1mu1pi0 (contained and fiducial)
     * 1: 1mu1pi0 (not contained or fiducial)               
     * 2: Other nu                
     * 3: cosmic
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the enumerated category of the interaction.   
     */
    template<class T>
      double category(const T & interaction)
      {
        // Cosmic background               
        double cat(3);

        // Signal
        if(cuts::signal_1mu0pi1pi0(interaction))
          {
	    if(cuts::fiducial_cut(interaction) && cuts::track_containment_cut(interaction))
	      {
		cat = 0;
	      }
	    else cat = 1;
          }
        // Neutrino Background                         
        else if(cuts::other_nu_1mu0pi1pi0(interaction))
          {
            cat = 2;
          }
        return cat;
      }

    /**
     * Variable for enumerating interaction categories. This classifies the          
     * interactions based on the visible final states.
     * 0: 1mu1pi0, 1: 1muNpi0, 2: 1muCex, 3: NC, 4: Other, 5: Cosmic        
     * @tparam T the type of interaction (true or reco).             
     * @param interaction to apply the variable on.      
     * @return the enumerated category of the interaction.                                               
     */
    template<class T> double category_topology(const T & interaction)
      {

	truth_inter s = cuts::true_interaction_info(interaction);
	
	// Cosmic
	uint16_t cat(7);

	// Neutrino
	if(s.is_neutrino)
	  {
	    // 1mu0pi1pi0 (in-phase)
	    if(s.num_primary_muons_thresh == 1 && s.num_primary_pions_thresh == 0 && s.num_primary_pi0s_thresh == 1 && s.is_cc && s.is_fiducial && s.has_contained_tracks)
	      {
		cat = 0;
	      }
	    // 1mu0pi1pi0 (out-of-phase)
	    else if( (s.num_primary_muons == 1 && s.num_primary_pions == 0 && s.num_primary_pi0s == 1 && s.is_cc) && (s.num_primary_muons_thresh != 1 || s.num_primary_pions_thresh != 0 || s.num_primary_pi0s_thresh != 1 || !s.is_fiducial || !s.has_contained_tracks) )
	      {
		cat = 1;
	      }
	    // 1muNpi1pi0
	    else if(s.num_primary_muons == 1 && s.num_primary_pions > 0 && s.num_primary_pi0s == 1 && s.is_cc)
	      {
		cat = 2;
	      }
	    // 1muNpi0pi0
	    else if(s.num_primary_muons == 1 && s.num_primary_pions > 0 && s.num_primary_pi0s == 0 && s.is_cc)
	      {
		cat = 3;
	      }
	    // 1muNpi0
	    else if(s.num_primary_muons == 1 && s.num_primary_pi0s > 1 && s.is_cc)
	      {
		cat = 4;
	      }
	    // NC 1pi0
	    else if(s.num_primary_muons == 0 && s.num_primary_pi0s == 1 && !s.is_cc)
	      {
		cat = 5;
	      }
	    // Other
	    else
	      {
		cat = 6;
	      }
	  }
	return cat;
      }

    /**
     * Variable for enumerating interaction categories. This categorization
     * uses the interaction type (generator truth) classify the interactions
     * 0: nu_mu CC QE, 1: nu_mu CC Res, 2: nu_mu CC MEC, 3: nu_mu CC DIS, 4: nu_mu CC Coh, 5: nu_e CC, 6: NC, 7: Cosmic
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the enumerated category of the interaction.
    */
    template<class T>
        double category_interaction_mode(const T & interaction)
        {
            double cat(7);

            if(interaction.nu_id > -1)
            {
                if(interaction.current_type == 0)
                {
                    if(abs(interaction.pdg_code) == 14)
                    {
                        if(interaction.interaction_mode == 0) cat = 0;
                        else if(interaction.interaction_mode == 1) cat = 1;
                        else if(interaction.interaction_mode == 10) cat = 2;
                        else if(interaction.interaction_mode == 2) cat = 3;
                        else if(interaction.interaction_mode == 3) cat = 4;
                        else cat = 8;
                    }
                    else cat = 5;
                }
                else cat = 6;
            }

            return cat;
        }

    /**
     * Find muon momentum (magnitude).
     * Assumes 1mu2gamma cut has been made.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to operate on.
     * @return the muon momentum.
     */
    template<class T> double muon_momentum_mag(const T & interaction)
      {
	if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                       {
                         truth_inter s = cuts::true_interaction_info(interaction);
                         return s.muon_momentum_mag;
                       }
	else
	  {
	    reco_pi0 s = cuts::reco_pi0_info(interaction);
            return s.muon_momentum_mag;
	  }
      }

    /**
     * Find muon-beam angle (cosine).
     * Assumes 1mu2gamma cut has been made.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to operate on.
     * @return the muon-beam angle.
     */
    template<class T> double muon_beam_costheta(const T & interaction)
      {
        if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                       {
                         truth_inter s = cuts::true_interaction_info(interaction);
                         return s.muon_beam_costheta;
                       }
        else
          {
            reco_pi0 s = cuts::reco_pi0_info(interaction);
            return s.muon_beam_costheta;
          }
      }

    /**
     * Find pi0 leading photon energy.
     * Assumes 1mu2gamma cut has been made.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to operate on.
     * @return the pi0 leading photon energy.
     */
    template<class T> double pi0_leading_photon_energy(const T & interaction)
      {
	if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
		       {
			 truth_inter s = cuts::true_interaction_info(interaction);
			 return s.pi0_leading_photon_energy;
		       }
	else
	  {
	    reco_pi0 s = cuts::reco_pi0_info(interaction);
	    return s.leading_photon_energy;
	  }
      }

    
    /**
     * Find pi0 leading photon conversion distance.
     * Assumes 2gamma cut has been made.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to operate on.
     * @return the pi0 leading photon conversion distance.
     */
    template<class T> double pi0_leading_photon_conv_dist(const T & interaction)
      {
	if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                       {
                         truth_inter s = cuts::true_interaction_info(interaction);
                         return s.pi0_leading_photon_conv_dist;
                       }
	else
	  {
	    reco_pi0 s = cuts::reco_pi0_info(interaction);
	    return s.leading_photon_conv_dist;
	  }
      }

    /**
     * Find pi0 leading photon cos(phi).
     * Assumes 2gamma cut has been made.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to operate on.
     * @return the pi0 leading photon cos(phi).
     */
    template<class T> double pi0_leading_photon_cosphi(const T & interaction)
      {
	if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                       {
                         truth_inter s = cuts::true_interaction_info(interaction);
                         return 0;
                       }
	else
	  {
	    reco_pi0 s = cuts::reco_pi0_info(interaction);
	    return s.leading_photon_cosphi;
	  }
      }
    
    /**
     * Find pi0 subleading photon energy.
     * Assumes 2gamma cut has been made.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to operate on.
     * @return the pi0 subleading photon energy.
     */
    template<class T> double pi0_subleading_photon_energy(const T & interaction)
      {
	if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                       {
                         truth_inter s = cuts::true_interaction_info(interaction);
                         return s.pi0_subleading_photon_energy;
                       }
        else
          {
            reco_pi0 s = cuts::reco_pi0_info(interaction);
            return s.subleading_photon_energy;
          }
      }

    /**
     * Find pi0 subleading photon conversion distance.
     * Assumes 2gamma cut has been made.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to operate on.
     * @return the pi0 subleading photon conversion distance.
     */
    template<class T> double pi0_subleading_photon_conv_dist(const T & interaction)
      {
	if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                       {
                         truth_inter s = cuts::true_interaction_info(interaction);
                         return s.pi0_subleading_photon_conv_dist;
                       }
	else
	  {
	    reco_pi0 s = cuts::reco_pi0_info(interaction);
	    return s.subleading_photon_conv_dist;
	  }
      }
    
    /**
     * Find pi0 subleading photon cos(phi).
     * Assumes 2gamma cut has been made.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to operate on.
     * @return the pi0 subleading photon cos(phi).
     */
    template<class T> double pi0_subleading_photon_cosphi(const T & interaction)
      {
	if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                       {
                         truth_inter s = cuts::true_interaction_info(interaction);
                         return 0;
                       }
	else
          {
            reco_pi0 s = cuts::reco_pi0_info(interaction);
            return s.subleading_photon_cosphi;
          }
      }

    /**
     * Find pi0 opening angle (cosine).
     * Assumes 1mu2gamma cut has been made.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to operate on.
     * @return the pi0 opening angle.
     */
    template<class T> double pi0_costheta(const T & interaction)
      {
	if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                       {
                         truth_inter s = cuts::true_interaction_info(interaction);
                         return s.pi0_costheta;
                       }
	else
	  {
	    reco_pi0 s = cuts::reco_pi0_info(interaction);
            return s.costheta;
	  }
      }
    
    /**
     * Find pi0 mass.
     * Assues 1mu2gamma cut has been made.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to operate on.
     * @return the pi0 mass.
     */
    template<class T> double pi0_mass(const T & interaction)
      {
	if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                       {
                         truth_inter s = cuts::true_interaction_info(interaction);
                         return s.pi0_mass;
                       }
	else
	  {
	    reco_pi0 s = cuts::reco_pi0_info(interaction);
	    return s.mass;
	  }
      }

    /**
     * Find pi0 momentum (magnitude).
     * Assumes 1mu2gamma cut has been made.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to operate on.
     * @return the pi0 momentum.
     */
    template<class T> double pi0_momentum_mag(const T & interaction)
      {
	if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                       {
                         truth_inter s = cuts::true_interaction_info(interaction);
                         return s.pi0_momentum_mag;
                       }
	else
	  {
	    reco_pi0 s = cuts::reco_pi0_info(interaction);
	    return s.momentum_mag;
	  }

      }
    
    /**
     * Find pi0-beam angle (cosine).
     * Assumes 1mu2gamma cut has been made.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to operate on.
     * @return the pi0-beam angle.
     */
    template<class T> double pi0_beam_costheta(const T & interaction)
      {
	if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                       {
                         truth_inter s = cuts::true_interaction_info(interaction);
                         return s.pi0_beam_costheta;
                       }
	else
	  {
	    reco_pi0 s = cuts::reco_pi0_info(interaction);
            return s.beam_costheta;
	  }
      }

    /**
     *
     *
     *
     *
     */
    template<class T> double transverse_momentum_mag(const T & interaction)
      {
	if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                       {
                         truth_inter s = cuts::true_interaction_info(interaction);
                         return s.transverse_momentum_mag;
                       }
	else
	  {
	    reco_pi0 s = cuts::reco_pi0_info(interaction);
	    return s.transverse_momentum_mag;
	  }
	
      }
    
}

#endif
