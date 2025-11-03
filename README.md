# My [T-Userdefined] Tallies

This repository is simply a collection of user-defined [PHITS](https://phits.jaea.go.jp/) tallies, for [T-Userdefined], that I have written for my own research work and have shared on the corresponding [PHITS Wiki topic](https://meteor.nucl.kyushu-u.ac.jp/phitswiki/usercodes).

Each directory contains its own `README.md` file providing more detail on the `usrtally.f` code accomanying it.

## Multi-coincidence neutron and gamma-ray detection

This tally allows event-by-event "list-mode" collection of events with arbitrarily-set minimum number of required interaction coincidences (as inputted and then evaluated by a [Counter]) each satisfying a set minimum energy deposition threshold; 
these requirements can be specified individually for both neutron interactions and gamma-ray interactions. An arbitrarily long list of cells/regions treated within the coincidence logic can be provided.

This tally was developed with a "camera-style" detector in mind that requires, for imaging, neutrons to undergo two interactions and gamma rays three interactions, each in a different region (scintillator bar) of the geometry (detector array). 
It outputs, for each coincidence event satisfying scoring criteria, the region number, energy deposition, and interction coordinates of each interaction, along with additional identifying information and (optionally) extra information about the nature of the interactions. 
In this way, it minimally emulates the output of a real detector system (interaction energy deposition and location), though without the associated experimental systematic uncertainties, necessary for particle imaging 
but also allows a much more detailed look at what exact interactions contribute to the energy deposition in each interaction only available in simulation.

## Ray-tracing through materials

This simple user-defined tally is used to identify the path lengths through distinct materials traversed by particles. 
If all source particles are neutrinos (whose reaction physics are disabled by default), this means all of these paths fall on a straight line, or "ray", traced by that source particle.
