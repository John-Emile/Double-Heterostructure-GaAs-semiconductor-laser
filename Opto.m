%% Implementing the constants and givens

%Givens
scattering_loss_coefficient = 2500; %per meter
radiative_recombination_lifetime = 2e-9; %2ns

%Constants
electron_charge = 1.602e-19;
boltzmann_constant = 1.38e-23;
speed_of_light = 3e8;

%Refrective index for double-heterostructure GaAs semiconductor
x = 0.3; %For Al0.3Ga0.7As
refractive_index1 = 3.5;
refractive_index2 = 3.5-0.62*x; 

%Recombination lifetime
internal_efficiency = 0.5;
recombination_lifetime = internal_efficiency * radiative_recombination_lifetime;

%Planks and reduced Planks constants
planks_constant = 6.626e-34;
reduced_planks_constant = planks_constant / (2 * pi);

%Electron and holes effective mass
electron_mass = 9.109e-31;
effective_hole_mass = 0.5 * electron_mass;
effective_electron_mass = 0.067 * electron_mass;
reduced_effective_mass = (effective_hole_mass * effective_electron_mass) / (effective_electron_mass + effective_hole_mass);
