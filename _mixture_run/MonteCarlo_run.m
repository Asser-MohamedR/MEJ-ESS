% MEJ-ESS: An Enhanced Monte Carlo Simulator for Electron Transport in Gaseous Detectors
%
% This program is a derivative work based on METHES, which is copyrighted by Mohamed Rabie.
%
% Original METHES Copyright and License:
%     Copyright (c) 2015, Mohamed Rabie
%     All rights reserved.
%     
%     Redistribution and use in source and binary forms, with or without
%     modification, are permitted provided that the following conditions are
%     met:
%     - Redistributions of source code must retain the above copyright notice, 
%     this list of conditions and the following disclaimer.
%     - Redistributions in binary form must reproduce the above copyright notice, 
%     this list of conditions and the following disclaimer in the documentation 
%     and/or other materials provided with the distribution
%     - Proper reference is made in publications reporting results obtained 
%     using this software. At present, the preferred way to reference METHES is as follows: 
%     M. Rabie and C.M. Franck, "A Monte Carlo collision code for electron transport in low temperature plasmas", 
%     to be submitted to J. Phys. D: Appl. Phys. (2015) 
%     as well as M. Rabie and C.M. Franck, METHES code, downloaded from www.lxcat.net/download/METHES/, 2015
%         
%     THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%     AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%     IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
%     ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
%     LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
%     CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
%     SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
%     INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
%     CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
%     ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
%     POSSIBILITY OF SUCH DAMAGE.
%
% MEJ-ESS Modifications Copyright:
%     Copyright (c) 2024-2025, [Asser Mohamed]
%     All rights reserved.
%
%     The modifications and enhancements in MEJ-ESS are licensed under the
%     same terms as the original METHES software.
%
% For more information, please refer to the README.md file in the repository.

clear all
close all
clc
pkg load statistics;
 kB = 1.3806488 *10^(-23);
 p = 1e5;
 Temp = 300  

%A=100:200:1000;
%B=2000:2000:50000;

%E_values = [10000:10000:50000]

N= p/(kB*Temp);

%EN_values = E_values.*1e21./N*100
EN_values = 300;

E_mean_results = zeros(size(EN_values));
drift_velocity_results = zeros(size(EN_values));
ion_coeff= zeros(size(EN_values));
att_coeff= zeros(size(EN_values));
long_diffusion=zeros(size(EN_values));
trans_diffusion=zeros(size(EN_values));
Trans_yield=zeros(size(EN_values));

Stucture_Size = 1;  %number of elements
all_results(Stucture_Size) = struct();


%% INPUT
functionsDir = '../_functions';
addpath(functionsDir);

% directories of cross sections
gasDir = {'../_Xsection/CO2_Biagi'};
%"Freon",96.2,"C4H10", 3.5,"SF6",0.3
% sumforumla gases
gas = {'CO2'};


% E/N in Townsend
 
 for k=2
    % mixing ratio
    mix = [1];
 for i = 1:length(EN_values)



% pressure in Pascal
p = 1e5;
% temperature in Kelvin
Temp = 300;
% start electron number
N0 = 1e4;
% maximum allowed electron number
Ne_max = 1e6;
% energy sharing factor in interval [0,1]
W = 0.5;
% error tolerance for flux drift velocity
w_err = 0.01;
% error tolerance for flux diffusion constant
DN_err = 0.01;
% minimum number of collisions before steady-state
col_equ = 1e7;
% maximum number of collisions of simulation
col_max =1.5e7;
% conserve electron number (1) or not(0)
conserve = 1;
% plot data (1) or not(0)
interactive  =0;
%%
% Penning transfer parameters
penning_enabled = 0; % Flag to enable/disable Penning effect (1/0)
r_p = 0.42 ; %penning probability

%% import cross sections
%%
Xsec = importLXcat;
Xsec.dir = gasDir;
Xsec.interactive = interactive;
Xsec = Xsec.importXsections();
Xsec = Xsec.fillThresholds();
Xsec = Xsec.getEnergy();
Xsec = Xsec.prepareForFit1();
Xsec = Xsec.momentum2elastic();
Xsec = Xsec.effective2elastic();
Xsec = Xsec.prepareForFit2();
Xsec = Xsec.totalXsection();

%% PLOTS
textsize = 14;
xscale = 'log' ; yscale = 'log';
Xsec.plotXsections(textsize,xscale,yscale);
%%

 


EN = EN_values(i);
sig = MonteCarlo;
sig.Xsec = Xsec;
sig.gas = gas;
% mixing ratio of gases
sig.mix = mix;
% molecular mass of gases in kg
sig = sig.mass_in_kg();

% conditions:
sig.N0 = N0 ; % number of initial electrons
sig.Ne_max = Ne_max ; % maximum number of electrons
sig.p = p; % Druck in Pascal
sig.Temp = Temp; % Temperatur in Kelvin
sig.EN = EN; % in Td

% numerics
sig.N_energy = 2000; % length of energy vector
sig.E_max = 1e6;
sig.E.energy = 0;
sig.W = W;
sig.iso = 1;
sig.w_err = w_err;
sig.DN_err = DN_err;
sig.col_equ = col_equ;
sig.col_max = col_max;
sig.conserve = conserve; % conserve (1) electron number or not (0)
sig.interactive = interactive; % plot (1) or do not plot data (0)
%%

% initial elctrons
sig.sigma_xyz = [0 0 0];
sig.pos_xyz = [0 0 0];

% check mixture
sig = sig.checkFractionSum();
% gas density in m^-3
sig = sig.gasNumberDensity();
% total cross section
sig = sig.totalXsection(); 
% electric field at positions x
sig = sig.solvePoisson_3D(100,0);
% set initial electron position and velocity 
sig = sig.initialParticles();

sig.counter = 0;
sig.flux.v_int_sum = [0 0 0];
sig.flux.D_sum = 0;
sig.flux.N = 0;
sig.E.E_sum = 0; 
sig.E.EEPF_sum = 0;
tic
while sig.End == 0

    % perform a flight for all electrons without collisions
    sig = sig.freeFlight();
    % calcultes the collective data of electron swarm
    sig = sig.collectMeanData();
    % get mean energy and EEDF
    sig = sig.energyData();
    % get flux transport coefficients 
    sig = sig.fluxData();
    % get bulk transport coefficients
    sig = sig.bulkData();
    % get reaction rates from electron number
    sig = sig.rateDataCount();
    % get reaction rates from convolution with EEDF
    sig = sig.rateDataConv();
    % decides which type of collision will happen
    sig = sig.collisionMatrix();
    % perform eleastic collisions
    sig = sig.elasticCollision();
    % perform ineleastic collisions
    sig = sig.inelasticCollision();
    % perform ionization collisions
    sig = sig.ionCollision();
    % perform attachment collisions
    sig = sig.attachCollision();
    % plot collective data
    sig = sig.plotMeanData();
    % check if equilibrium is reached 
    sig = sig.checkSST();
    % print results in command window
    sig = sig.printOnScreen();
    % end simulation
    sig = sig.endSimulation();
       
end

toc

 results = sig.getResults();

 
 E_mean_results(i) = results.E.E_mean
 drift_velocity_results(i) = results.flux.w(3)
 ion_coeff(i) = results.rates.count.ion_coeff/1e2
 att_coeff(i) = results.rates.count.att_coeff/1e2
 long_diffusion(i) = sqrt( results.bulk.DN(3) * 2 * 1e4 *1/(drift_velocity_results(i)*100*N))*1e4
 Trans_yield(i) = sqrt(results.bulk.DN(1) * results.bulk.DN(2));
 trans_diffusion(i) = sqrt( Trans_yield(i) * 2 * 1e4 *1/(drift_velocity_results(i)*100*N))*1e4;
close all
%clear results
    % Optional: Save individual results:
   
 %save(['results_EN_' num2str(EN)], 'results');

     end

drift_velocity_results_cm_ut=drift_velocity_results./1e4




all_results(k+1).diffusion = long_diffusion; 
all_results(k+1).trans_diffusion = trans_diffusion; 
all_results(k+1).drift_velocity_results = drift_velocity_results;
all_results(k+1).EN_values = EN_values;
all_results(k+1).E_mean_results = E_mean_results;
all_results(k+1).Ion_Coeff=ion_coeff;
all_results(k+1).Att_Coeff=att_coeff;
all_results(k+1).drift_velocity_cm_ut = drift_velocity_results_cm_ut;
tosave=all_results(k+1);
save('all_results_octave.mat', 'tosave');
  end


%% Plot and save combined results
