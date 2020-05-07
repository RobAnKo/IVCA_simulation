# IVCA_simulation
A set of functions to model the growth of colonies in vitro after differential irradiation

colony_size_sim_head_script initializes options and runs colony_size_sim_main.

colony_size_sim_main calls colony_size_sim_per_dose for every given dose

colony_size_sim_per_dose calls colony_size_sim_per_well for every given well

colony_size_sim_per_well_arrayed does the actual simulation for one well

concatenation of the results happens within colony_size_sim_main.
