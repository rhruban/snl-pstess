function smppi(i,k,flag)
% Syntax: exc_exac8b(i,k,flag)
%
% Purpose: exac8b excitation system, (exc_con(i,1) = ?4?)
%          with vectorized computation option
%
% Input:   i - generator number
%              if 0 - vectorized computation
%          k - integer time
%          flag - 0 - initialization
%                 1 - network interface computation
%                 2 - generator dynamics computation
%
% See Also: smpexc, exc_dc12, exc_st3

% g.exc.exc_con
% Standard for DC12                           |For exac8b
% 1)  Exciter Type                            |4
% 2)  Machine Number                          |Mac #
% 3)  Input Filter Time Constant, TR          |TR
% 4)  Voltage Regulator Gain, KA              |Kvp
% 5)  Voltage Regulator Time Constant, TA     |TA
% 6)  Volatge Regulator Time Constant, TB     |Kvi
% 7)  Voltage Regulator Time Constant, TC     |Kc
% 8)  Max Voltage Regulator Output, VRmax     |VRmax
% 9)  Min Voltage Regulator Output, VRmin     |VRmin
% 10) Exciter Constant, KE                    |KE
% 11) Exciter Time Contstant, TE              |TE
% 12) E1                                      |E1
% 13) SE(E1)                                  |SE(E1)
% 14) E2                                      |E2
% 15) SE(E2)                                  |SE(E2)
% 16) Stabilizing Gain, KF                    |Kvd
% 17) Stabilizing Time Constant, TF           |Tvd
% 18)                                         |Kd
% 19)                                         |Vimax
% 20)                                         |limflg

% default exc_con
% exc_con = [
% % Type Mac TR    KVP TA    KVI KC   VRmax VRmin KE   TE   E1   SE1    E2   SE2   KVD   TVD   KD   Vimax  limflg
%   4    1   0.02  4   0.02  4   0.1  10    -10   1.0  0.2  5.0  0.001  7.0  0.01  20.0  0.05  1.0  0.1    0.0;];

%-----------------------------------------------------------------------------%
% Version History
%
% Version: 1.0
% Date:    January 2023
% Author:  Ronald Hruban
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

if (i ~= 0)
    error('exc_xac8b: all calculations must be vectorized.');
end

% vectorized computation
if (flag == 0)   % initialization
    % n           -- index of machine(s) with simple exac8b exciters
    % err         -- summing junction error
    % - states -
    % Efd         -- Actually V_E
    % R_f         -- State for Derivative Control
    % V_R         -- Voltage Regulator Output V_R
    % V_As        -- integrator state
    % V_TR        -- input filter state
    % exc_pot(,3) -- reference voltage
    if (g.exc.n_smppi ~= 0)
        n = g.mac.mac_int(g.exc.exc_con(g.exc.smppi_idx,2)); % machine indexer

        % Initialize V_E
        % Assume 0 < I_N < 0.433
        %g.exc.Efd(g.exc.smppi_idx,1) = g.mac.vex(n,1)./g.mac.mac_spd(n,1)+0.5774.*g.exc.exc_con(g.exc.smppi_idx,7).*g.mac.fldcur(n,1); % Calc V_E
        g.exc.Efd(g.exc.smppi_idx,1) = g.mac.vex(n,1)./g.mac.mac_spd(n,1);

        % % 0.433 < I_N < 0.75
        % I_N = g.exc.exc_con(g.exc.smppi_idx,7).*g.mac.fldcur(n,1)./g.exc.Efd(g.exc.smppi_idx,1);
        % FEX0 = find(I_N > 0.433);
        % if ~isempty(FEX0)
        %     g.exc.Efd(g.exc.smppi_idx(FEX0),1) = sqrt(((g.mac.vex(n(FEX0),1)./g.mac.mac_spd(n(FEX0),1)).^2+(g.exc.exc_con(g.exc.smppi_idx(FEX0),7).*g.mac.fldcur(n(FEX0),1)).^2)./0.75);
        % end
        % 
        % % 0.75 < I_N < 1
        % I_N = g.exc.exc_con(g.exc.smppi_idx,7).*g.mac.fldcur(n,1)./g.exc.Efd(g.exc.smppi_idx,1);
        % FEX1 = find(I_N > 0.75);
        % if ~isempty(FEX1)
        %     g.exc.Efd(g.exc.smppi_idx(FEX1),1) = g.mac.vex(n(FEX1),1)./g.mac.mac_spd(n(FEX1),1)./1.732+g.exc.exc_con(g.exc.smppi_idx(FEX1),7).*g.mac.fldcur(n(FEX1),1);
        % end
        % 
        % % I_N > 1
        % I_N = g.exc.exc_con(g.exc.smppi_idx,7).*g.mac.fldcur(n,1)./g.exc.Efd(g.exc.smppi_idx,1);
        % FEX2 = find(I_N > 1);
        % if ~isempty(FEX2)
        %     error('exc_exac8b: F_EX is zero because In is greater than 1');
        % end
        % 
        % % I_N < 0
        % FEX3 = find(I_N < 0);
        % if ~isempty(FEX3)
        %     error('exc_exac8b: I_N is negative');
        % end

        % check E_fd >= 0 !!!


        % Saturation Initialization
        % add if else for providing A and B vice Se and E
        % A
        g.exc.exc_pot(g.exc.smppi_idx,1) = exp((log(g.exc.exc_con(g.exc.smppi_idx,13)).*g.exc.exc_con(g.exc.smppi_idx,14)-log(g.exc.exc_con(g.exc.smppi_idx,15)).*g.exc.exc_con(g.exc.smppi_idx,12))...
            ./(g.exc.exc_con(g.exc.smppi_idx,14)-g.exc.exc_con(g.exc.smppi_idx,12)));
        % B
        g.exc.exc_pot(g.exc.smppi_idx,2) = (-log(g.exc.exc_con(g.exc.smppi_idx,13))+log(g.exc.exc_con(g.exc.smppi_idx,15)))...
            ./(g.exc.exc_con(g.exc.smppi_idx,14)-g.exc.exc_con(g.exc.smppi_idx,12));
        % Se
        Se = 0; %g.exc.exc_pot(g.exc.smppi_idx,1).*exp(g.exc.exc_pot(g.exc.smppi_idx,2).*g.exc.Efd(g.exc.smppi_idx,1));

        % need if then for no Saturation !!!

        % need if then V_E < A !!!


        % output of transducer, state
        g.exc.V_TR(g.exc.smppi_idx,1) = g.mac.eterm(n,1);
        % need v_comp calc !!!

        % output of Ta, state
        g.exc.V_R(g.exc.smppi_idx,1) = g.exc.exc_con(g.exc.smppi_idx,18).*g.mac.fldcur(n,1) + (g.exc.exc_con(g.exc.smppi_idx,10) + Se).*g.exc.Efd(g.exc.smppi_idx,1);
        
        
        % if then for PID vs PD vs PI vs P
        con_type = 2;
        switch con_type
            case 1 %PID

            case 2 %PD
                % output of PID integrator, state
                % g.exc.V_As(g.exc.smppi_idx,1) = g.exc.V_R(g.exc.smppi_idx,1);
                g.exc.V_As(g.exc.smppi_idx,1) = zeros(g.exc.n_smppi,1);

                % V_ref
                g.exc.exc_pot(g.exc.smppi_idx,3) = g.exc.V_R(g.exc.smppi_idx,1)./g.exc.exc_con(g.exc.smppi_idx,4)+g.exc.V_TR(g.exc.smppi_idx,1);% - v_sig

                % output of PID 'derivative', state
                % g.exc.R_f(g.exc.smppi_idx,k) = zeros(g.exc.n_smppi,1);
                g.exc.R_f(g.exc.smppi_idx,k) = g.exc.exc_con(g.exc.smppi_idx,16)./g.exc.exc_con(g.exc.smppi_idx,17).*(g.exc.exc_pot(g.exc.smppi_idx,3)-g.exc.V_TR(g.exc.smppi_idx,1)); % add V_sig between vref and V_TR

            case 3 %PI

            case 4 %P

        end

        % should I check intit within limiters??????????????????????????

    end
end

if (flag == 1)   % network interface computation
    % vectorized computation
    if (g.exc.n_smppi ~= 0)                 % check for simple PI exciters
        n = g.mac.mac_int(g.exc.exc_con(g.exc.smppi_idx,2));
        I_N = g.exc.exc_con(g.exc.smppi_idx,7).*g.mac.fldcur(n,k)./g.exc.Efd(g.exc.smppi_idx,k);
        Fex = nan(size(I_N));
        
        %find 0<I_N<0.433
        FEX0 = find(I_N < 0.433);
        if ~isempty(FEX0)
            Fex(FEX0) = 1-0.5771.*I_N(FEX0);
        end

        %find 0.433<I_N<0.75
        FEX1 = find((0.433<=I_N)&(I_N < 0.75));
        if ~isempty(FEX1)
            Fex(FEX1) = sqrt(0.75-I_N(FEX1).^2);
        end

        %find 0.75<I_N<1
        FEX2 = find((0.75<=I_N)&(I_N < 1));
        if ~isempty(FEX2)
            Fex(FEX2) = 1.732.*(1-I_N(FEX2));
        end
        %find 1<In
        FEX3 = find(I_N > 1);
        if ~isempty(FEX3)
            Fex(FEX3) = 0;
        end
        
        % I_N < 0
        FEX4 = find(I_N < 0);
        if ~isempty(FEX4)
            error('exc_exac8b: I_N is negative');
        end

        g.mac.vex(n,k) = g.exc.Efd(g.exc.smppi_idx,k).*g.mac.mac_spd(n,k);%.*Fex;
    end
end

if (flag == 2)   % exciter dynamics calculation
    % vectorized computation

    if (g.exc.n_smppi ~= 0)
        % machine numbers
        n = g.mac.mac_int(g.exc.exc_con(g.exc.smppi_idx,2));
        TR = g.exc.smppi_TR_idx;
        no_TR = g.exc.smppi_noTR_idx;

        % for exciters with zero TR
        if ~isempty(no_TR)
            n_nTR = n(no_TR);
            g.exc.dV_TR(g.exc.smppi_idx(no_TR),k) = zeros(length(no_TR),1);
            g.exc.V_TR(g.exc.smppi_idx(no_TR),k) = g.mac.eterm(n_nTR,k);
        end

        % for exciters with nonzero TR
        if ~isempty(TR)
            n_TR = n(TR);
            g.exc.dV_TR(g.exc.smppi_idx(TR),k) = ...
                (g.mac.eterm(n_TR,k) - g.exc.V_TR(g.exc.smppi_idx(TR),k)) ...
                ./g.exc.exc_con(g.exc.smppi_idx(TR),3);
        end

        % error defined for all simple exciters
        err = g.exc.exc_sig(g.exc.smppi_idx,k) ...
            + g.exc.exc_pot(g.exc.smppi_idx,3) ...
            - g.exc.V_TR(g.exc.smppi_idx,k) ...
            + g.exc.pss_out(g.exc.smppi_idx,k);

        myLim = g.exc.exc_con(g.exc.smppi_idx,19);
        err(err>g.exc.exc_con(g.exc.smppi_idx,19))=myLim(err>g.exc.exc_con(g.exc.smppi_idx,19));
        err(err<-g.exc.exc_con(g.exc.smppi_idx,19))=-myLim(err<-g.exc.exc_con(g.exc.smppi_idx,19));

        % Integral State
        % g.exc.dV_As(g.exc.smppi_idx,k) = err*g.exc.exc_con(g.exc.smppi_idx,4);
        g.exc.dV_As(g.exc.smppi_idx,k) = zeros(g.exc.n_smppi,1);

        % need check if kvd is zero logic??

        g.exc.dR_f(g.exc.smppi_idx,k) = (err.*g.exc.exc_con(g.exc.smppi_idx,16)./g.exc.exc_con(g.exc.smppi_idx,17)-g.exc.R_f(g.exc.smppi_idx,k))./g.exc.exc_con(g.exc.smppi_idx,17);

        g.exc.dV_R(g.exc.smppi_idx,k) = (g.exc.exc_con(g.exc.smppi_idx,4).*err...
            +g.exc.V_As(g.exc.smppi_idx,k)...
            +g.exc.exc_con(g.exc.smppi_idx,16)./g.exc.exc_con(g.exc.smppi_idx,17).*err...
            -g.exc.R_f(g.exc.smppi_idx,k)...
            -g.exc.V_R(g.exc.smppi_idx,k))./g.exc.exc_con(g.exc.smppi_idx,5);

        % anti-windup reset
        maxlmt = find(g.exc.V_R(g.exc.smppi_idx,k) ...
            > g.exc.exc_con(g.exc.smppi_idx,8));

        if ~isempty(maxlmt)
            g.exc.V_R(g.exc.smppi_idx(maxlmt),k) = ...
                g.exc.exc_con(g.exc.smppi_idx(maxlmt),8);

            pos_rate = find(g.exc.dV_R(g.exc.smppi_idx(maxlmt),k) > 0);

            prl = length(pos_rate);
            if (prl ~= 0)
                g.exc.dV_R(g.exc.smppi_idx(maxlmt(pos_rate)),k) = ...
                    zeros(prl,1);
            end
        end

        minlmt = find(g.exc.V_R(g.exc.smppi_idx,k) ...
            < g.exc.exc_con(g.exc.smppi_idx,9));

        if ~isempty(minlmt)
            g.exc.V_R(g.exc.smppi_idx(minlmt),k) = ...
                g.exc.exc_con(g.exc.smppi_idx(minlmt),9);

            neg_rate = find(g.exc.dV_R(g.exc.smppi_idx(minlmt),k) < 0);

            nrl = length(neg_rate);
            if (nrl ~= 0)
                g.exc.dV_R(g.exc.smppi_idx(minlmt(neg_rate)),k) = ...
                    zeros(nrl,1);
            end
        end
        

        %need Se if then statement
        Se = 0; %g.exc.exc_pot(g.exc.smppi_idx,1).*exp(g.exc.exc_pot(g.exc.smppi_idx,2).*g.exc.Efd(g.exc.smppi_idx,k));

        g.exc.dEfd(g.exc.smppi_idx,k) = ...
            (g.exc.V_R(g.exc.smppi_idx,k)-g.exc.exc_con(g.exc.smppi_idx,18).*g.mac.fldcur(n,k)-(g.exc.exc_con(g.exc.smppi_idx,10)+Se).*g.exc.Efd(g.exc.smppi_idx,k)) ...
            ./g.exc.exc_con(g.exc.smppi_idx,5);

        % Te integral limit
        minInt = find(g.exc.Efd(g.exc.smppi_idx,k) ...
            < zeros(size(g.exc.Efd(g.exc.smppi_idx,k))));

        if ~isempty(minInt)
            g.exc.Efd(g.exc.smppi_idx(minInt),k) = ...
                zeros(size(g.exc.Efd(g.exc.smppi_idx(minInt),k)));

            neg_Int = find(g.exc.dEfd(g.exc.smppi_idx(minInt),k) < 0);

            nIl = length(neg_Int);
            if (nIl ~= 0)
                g.exc.dEfd(g.exc.smppi_idx(minInt(neg_Int)),k) = ...
                    zeros(nIl,1);
            end
        end

    end

end

end  % function end

% eof
