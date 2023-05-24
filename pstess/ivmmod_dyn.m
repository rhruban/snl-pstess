function [dc,Ec,ddelta_states,dE_states,delta_statesIni,E_statesIni] = ivmmod_dyn(delta_states,E_states,bus,Time,kSim,Flag)
% Syntax: [dc,Ec,ddelta_states,dE_states,delta_statesIni,E_statesIni] = ...
%         ivmmod_dyn(delta_states,E_states,bus,Time,kSim,Flag)
%
% Purpose: Implement state or output variables to model power injection
%
% Inputs:
%   Time = vector of simulation time
%   kSim = current simulation index.  Current time = Time(kSim).
%   Flag:
%       If Flag==0, Initialize d_statesIni, E_statesIni at t = 0.
%       If Flag==1, Calculate d, E at Time(kSim)
%       If Flag==2, Calculate dd_states, dE_states at Time(kSim)
%   d_state = cell array of mac_ang injection states.  Used to set d.
%       d_state{k} = column vector of states corresponding to the kth IVM in
%       mac_con.
%   E_state = cell array of E injection states.  Same format as d_states.
%   bus = initial bus matrix from the solved power flow.
%
% Outputs:
%   d_statesIni == cell array of initial of d_states
%   E_statesIni == cell array of initial of E_states
%   dd_state = cell array of d/dt of d_states.
%   dE_state = cell array of d/dt of E_states.
%   dc = n_ivm by 1 column vector of ivmmod_d_sig commands at t = Time(kSim).
%       ivmmod_d_sig set the mac_ang for an IVM generator in mac_ivm.m
%       ivmmod_d_sig(k,kSim) = d(k).  k corresponds to the kth IVM in mac_con.
%   Ec = n_ivm by 1 column vector of ivmmod_e_sig commands at t = Time(kSim).
%       ivmmod_e_sig set the edprime for an IVM generator in mac_ivm.m
%       ivmmod_e_sig(k,kSim) = E(k).  k corresponds to the kth IVM in mac_con.
%
% Global:
%       g.mac.ivmmod_data = general variable for storing data when necessary

%-----------------------------------------------------------------------------%
% Version history
%
% Version: 1.0 (initial version)
% Author:  D. Trudnowski
% Date:    2020
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

% g.mac.ivmmod_d_sig -- mac_ang command signal
% g.mac.ivmmod_e_sig -- edprime command signal

% bus numbers where ivm's are connected
busnum = g.bus.bus_int(g.mac.mac_con(g.mac.mac_ivm_idx,2));

%% Initialize output variables
dc = zeros(g.mac.n_ivm,1);
Ec = zeros(g.mac.n_ivm,1);
ddelta_states = cell(g.mac.n_ivm,1);
dE_states = cell(g.mac.n_ivm,1);
delta_statesIni = cell(g.mac.n_ivm,1);
E_statesIni = cell(g.mac.n_ivm,1);

%% Define and initialize state derivatives at t = 0.
if (Flag == 0)
    % Note: there are no differential equations, just make the DFE zero
    g.mac.ivmmod_data = nan;  % not used
    for k = 1:g.mac.n_ivm
        if  (g.playin.n_playin~=0)&&(any(g.mac.mac_ivm_idx(k)==g.playin.playin_ivm_idx))
            % play_in gen dynamics (do not change)
            delta_statesIni{k}(1,1) = 0;
            E_statesIni{k}(1,1) = 0;
        else
            % other IVM dynamics
            delta_statesIni{k}(1,1) = 0;
            E_statesIni{k}(1,1) = 0;
        end

    end

%% Calculate dc and Ec
elseif (Flag == 1)
    for k = 1:g.mac.n_ivm
        if  (g.playin.n_playin~=0)&&(any(g.mac.mac_ivm_idx(k)==g.playin.playin_ivm_idx))
            % play_in gen dynamics (do not change)
            Ec(k) = abs(g.playin.PMU.v(kSim));
            dc(k) = angle(g.playin.PMU.v(kSim));
        else
            % other IVM dynamics
            if (Time(kSim) > 1)
                dE = 0.1;
            else
                dE = 0;
            end

            % edprime(,1) is initial internal voltage E
            Ec(k) = g.mac.edprime(g.mac.mac_ivm_idx(k),1) + dE;
            if (Time(kSim) > 4)
                da = 0.2;
            else
                da = 0;
            end

            % mac_ang(,1) is initial internal voltage d
            dc(k) = g.mac.mac_ang(g.mac.mac_ivm_idx(k),1) + da;
        end
    end

%% Calculate derivatives
elseif (Flag == 2)
    % No differential equations for this example, set derivatives to zero
    for k = 1:g.mac.n_ivm
        if  (g.playin.n_playin~=0)&&(any(g.mac.mac_ivm_idx(k)==g.playin.playin_ivm_idx))
            % play_in gen dynamics (do not change)
            ddelta_states{k}(1,1) = 0;
            dE_states{k}(1,1) = 0;
        else
            % other IVM dynamics
            ddelta_states{k}(1,1) = 0;
            dE_states{k}(1,1) = 0;
        end
    end
end

end  % function end

% eof
