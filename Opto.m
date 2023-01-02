close all;
clc;
clear all;

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

%Energy bandgap and wavelength
band_gap = 1.424; 
wavelength_at_band_gap = 1.24e-6 / band_gap;

%Valence band and conduction band postions
valence_band_edge = 0;
conduction_band_edge = band_gap - valence_band_edge;

%y = logspace(a,b,n) generates n points between decades 10^a and 10^b.
carrier_concentration = logspace(0,30,5000); 

%% Quasi-Fermi levels and the electron/hole concentration using zero-Kelvin approximation
fermi_conduction_band_edge = conduction_band_edge + (3*pi*pi)^(2/3) * (reduced_planks_constant^2) * (carrier_concentration.^(2/3)) / (2*effective_electron_mass*electron_charge); %Efc
fermi_valence_band_edge = valence_band_edge -(3*pi*pi)^(2/3) * (reduced_planks_constant^2) * (carrier_concentration.^(2/3)) / (2*effective_hole_mass*electron_charge); %Efv

%For debugging, q_debugging should be zero
quasi_fermi_levels_difference = band_gap + (3*pi*pi)^(2/3) * (reduced_planks_constant^2) * (carrier_concentration.^(2/3)) / (2*reduced_effective_mass*electron_charge); %Efc-Efv
q_debugging = (fermi_conduction_band_edge-fermi_valence_band_edge)-quasi_fermi_levels_difference;

%Plotting
figCounter = 1; %For the ease of plotting
figure(figCounter)
figCounter = figCounter +1;

%semilogy(X,Y) plots x- and y-coordinates using a linear scale on the x-axis and a base-10 logarithmic scale on the y-axis.
semilogy(fermi_valence_band_edge, carrier_concentration) 
hold on %To plot both of them on the same graph
semilogy(fermi_conduction_band_edge, carrier_concentration)
xlim([-3, 3]) %Limits choosed based on trial and error
grid on
title("Quasi-Fermi levels vs electron/hole concentration using 0-k approximation")
ylabel("Electron/Hole concentration")
xlabel("Quasi-Fermi levels")

%% Quasi-Fermi levels and the electron/hole concentration exact
conduction_band_edge1 = -3:0.05:3; %Limits choosed based on trial and error

%Constant multiplied by the integration
const_electron_denisty = sqrt(2) * (effective_electron_mass^(3/2)) / ((pi^2) * (reduced_planks_constant^3)); 

%Function 
electron_concentration_dE = @(E) sqrt((E - band_gap) * electron_charge) ./ (1+exp((E - conduction_band_edge1) .* electron_charge / (boltzmann_constant*300)));

%q = integral(fun,xmin,xmax,Name,Value)numerically integrates function fun
%from xmin to xmax using global adaptive quadrature and default error tolerances with one or more Name,Value pair arguments
electron_concentration = const_electron_denisty * integral(electron_concentration_dE, band_gap, inf, 'ArrayValued', true) .* electron_charge;

fermi_valence_band_edge1 = -3:0.005:3;
const_hole_denisty = sqrt(2) * (effective_hole_mass ^ (3/2)) / ((pi ^ 2) * (reduced_planks_constant ^ 3));
hole_concentration_dE = @(E) const_hole_denisty * sqrt(-E * electron_charge) ./ (1 + exp((-E + fermi_valence_band_edge1) .* electron_charge / (boltzmann_constant * 300)));
hole_concentration = integral(hole_concentration_dE, -inf, 0, 'ArrayValued', true) .* electron_charge;

%Plotting
figure(figCounter)
figCounter = figCounter +1;
semilogy(conduction_band_edge1, electron_concentration)
hold on
semilogy(fermi_valence_band_edge1, hole_concentration)
xlim([-3, 3])
grid on
title("Quasi-Fermi levels vs electron/hole concentration exact")
ylabel("Electron/Hole concentration")
xlabel("Quasi-Fermi levels")


%% Optical Joint Denisty of states

%Dimensions
thickness = 0.1e-6;
len = 1000e-6;
width = 10e-6;

%Current density assumptions to get electrons concentration
%current_density = 0:10e3:700e4; %per meter squared
current_density = [200e4 400e4 600e4 800e4];
frequency = logspace(14,15,1000);
%frequency = 1e14:9e11:1e15;

%Initial value for optical joint denisty of states
optical_joint_density_of_states = zeros(1, 1000);

for x = 1:length(frequency)
    %Make sure that photon energy > Energy gap
    if planks_constant * frequency(x) >= band_gap * electron_charge
        optical_joint_density_of_states(x) = (2*sqrt(2) * (reduced_effective_mass ^ (3/2)) * sqrt(planks_constant * frequency(x) - band_gap * electron_charge)) / (pi * reduced_planks_constant ^ (2));
    else
        optical_joint_density_of_states(x) = 0;
    end
end
figure(figCounter)
figCounter = figCounter + 1;
plot(frequency, optical_joint_density_of_states)
title("Optical joint density of states vs frequency")
ylabel("Optical joint denisty of states")
xlabel("Frequency")
grid on

%% Gain & Fermi inversion factor
%Initial value, 1000 is the frequency length
optical_absorption_coefficient = zeros(1, 1000);
gain_coefficient = zeros(1, 1000);
gain_const = zeros(1, 1000);
gain = zeros(1, 1000);
maximum_gain = zeros(1, length(current_density));


for i = 1:length(current_density)
    carrier_density = current_density(i) * recombination_lifetime/ (thickness * electron_charge);
    
    quasi_fermi_conduction_band_diff = (3*pi*pi) ^ (2/3) * (reduced_planks_constant^2) * (carrier_density^(2/3)) / (2*effective_electron_mass*electron_charge); %Efc - Ec
    quasi_fermi_valence_band_diff = (3*pi*pi) ^ (2/3) * (reduced_planks_constant^2) * (carrier_density ^(2/3)) / (2*reduced_effective_mass*electron_charge); %Ev - Efv
 
    for c = 1:length(frequency)
        %Get (reduced plank const * k)^2 alone --> Same method as Sheet 6 pb 2
        hk_squared = (planks_constant*frequency(c) - (band_gap * electron_charge)) * 2 * reduced_effective_mass; %(h_bar*k)^2 = (E-Eg)*2mr
        
        %Use it in E1-Ev = -(hk_squared)/2mh
        valence_band_hk_const = -hk_squared /(2*effective_hole_mass*electron_charge);
        
        %Use it in E2-Ec = (hk_squared)/2me
        conduction_band_hk_const = hk_squared / (2 * effective_electron_mass * electron_charge);
        
        %Get E2-Efc & E1-Efv to replace them and get Fermi inversion factor
        conduction_band_edge_energy_2 = (conduction_band_hk_const) - (quasi_fermi_conduction_band_diff); %Efc-E2
        valence_band_edge_energy_2 = (valence_band_hk_const) + (quasi_fermi_valence_band_diff); %E1-Efv
        
        %Get exponentiels to get Fermi inversion factor
        conduction_band_filling_factor_energy_2 = 1 / (1 + exp(conduction_band_edge_energy_2 * electron_charge / (boltzmann_constant * 300)));
        valence_band_filling_factor_energy_1 = 1 / (1 + exp(valence_band_edge_energy_2 * electron_charge / (boltzmann_constant * 300)));
        
        finv(c) = conduction_band_filling_factor_energy_2 - valence_band_filling_factor_energy_1;
        
        gain_const(c) = ((speed_of_light/ frequency(c)) ^ 2) / (8 * pi * radiative_recombination_lifetime * (refractive_index1 ^ 2)); 
        
        gain_coeff(c) = gain_const(c) * optical_joint_density_of_states(c) * finv(c);
    end
    maximum_gain(i) = max(gain_coeff);

%Plotting the fermi inversion factor
figure(figCounter)
hold on
plot(frequency,finv)
xlim([2e14 6e14])
title("Fermi Inversion Factor")
ylabel("Fermi inversion factor")
xlabel("Frequency")
grid on
%figCounter = figCounter +1;

%Plotting the gain with the current denisty
figure(figCounter +1)
hold on
%figCounter = figCounter +1;
plot(frequency,gain_coeff)
xlim([3e14 4e14])
ylim([-2e4 4e4])
title("Gain coefficient")
xlabel("Frequency")
ylabel("Gain coefficient")
grid on 

end

%Plotting the max gain with the current denisty
figure(figCounter +2)
figCounter = figCounter +3;
plot(current_density,maximum_gain)
%xlim([0 10e6])
xlabel("injection current")
ylabel("Peak gain coefficient")
grid on

%% Finding the transparency current
Jtr = 230;
scattering_loss_coefficient = 25;
J_transparency = Jtr +(scattering_loss_coefficient/internal_efficiency);
Jth = 1.5 * J_transparency;
R1=1;
gamma_0 = internal_efficiency*(Jth-J_transparency);
R2 = (1/R1)*exp(-2*len*(gamma_0 - scattering_loss_coefficient));
disp('power reflectivity of the other mirror is');
disp(R2);



