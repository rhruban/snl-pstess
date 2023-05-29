function mac_sub(i,k,bus,flag)
% Syntax: mac_sub(i,k,bus,flag)
%
% Purpose: voltage-behind-subtransient-reactance generator
%          model, with vectorized computation option
%
% Note:    state variables include mac_ang, mac_spd, eqprime,
%          psikd, edprime, psikq
%
% Input:   i - generator number
%            - 0 for vectorized computation
%          k - integer time
%          bus - solved loadflow bus data
%          flag - 0 - initialization
%                 1 - network interface computation
%                 2 - generator dynamics computation
%                 3 - state matrix building
%
% See Also: mac_em, mac_tra
%
% Algorithm: PSLF model from John Undrill without saturation

%-----------------------------------------------------------------------------%
% Version history
%
% Version:  2.x?
% Purpose:  Implemented genrou model of PSLF
% Date:     July 2008
% Author:   Dan Trudnowski
%
% Version:  2.1
% Date:     September 1997
% Purpose:  Saturation modified for low Eqprime
%
% Version:  2.0
% Date:     June 1996
% Author:   Graham Rogers
% Purpose:  change to allow multple generator model types
%
% Version:  1.0
% Author:   Joe H. Chow
% Date:     March 1991
% Modified: Correction to saturation GJR May 6 1995
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

if (g.mac.n_sub ~= 0 && i == 0)
    if (flag == 0) % initialization
        % vectorized computation

        % check parameters
        uets_idx = find(g.mac.mac_con(g.mac.mac_sub_idx,8) ...
                        ~= g.mac.mac_con(g.mac.mac_sub_idx,13));

        if ~isempty(uets_idx)
            g.mac.mac_con(g.mac.mac_sub_idx(uets_idx),13) = ...
                g.mac.mac_con(g.mac.mac_sub_idx(uets_idx),8);
            %
            wstr = '\nmac_sub: xqpp made equal to xdpp at generator %0.0f.';
            warning(sprintf(wstr,g.mac.mac_sub_idx(uets_idx)));
        end

        notp_idx = find(g.mac.mac_con(g.mac.mac_sub_idx,14) == 0);
        if ~isempty(notp_idx)
            g.mac.mac_con(g.mac.mac_sub_idx(notp_idx),14) = ...
                999.0*ones(length(notp_idx),1);
        end

        notpp_idx = find(g.mac.mac_con(g.mac.mac_sub_idx,15) == 0);
        if ~isempty(notpp_idx)
            g.mac.mac_con(g.mac.mac_sub_idx(notpp_idx),15) = ...
                999.0*ones(length(notpp_idx),1);

            % set x'q = x"q
            g.mac.mac_con(g.mac.mac_sub_idx(notpp_idx),12) = ...
                g.mac.mac_con(g.mac.mac_sub_idx(notpp_idx),13);
        end

        % busnum -- bus number vector
        busnum = g.bus.bus_int(g.mac.mac_con(g.mac.mac_sub_idx,2));

        % mac_pot(,1) -- scaled MVA base
        % mac_pot(,2) -- base kV
        g.mac.mac_pot(g.mac.mac_sub_idx,1) = g.sys.basmva*ones(g.mac.n_sub,1) ...
                                             ./g.mac.mac_con(g.mac.mac_sub_idx,3);

        g.mac.mac_pot(g.mac.mac_sub_idx,2) = ones(g.mac.n_sub,1);

        g.mac.mac_pot(g.mac.mac_sub_idx,8) = g.mac.mac_con(g.mac.mac_sub_idx,7) ...
                                             - g.mac.mac_con(g.mac.mac_sub_idx,4);

        g.mac.mac_pot(g.mac.mac_sub_idx,9) = ...
            (g.mac.mac_con(g.mac.mac_sub_idx,8) ...
             - g.mac.mac_con(g.mac.mac_sub_idx,4)) ...
            ./g.mac.mac_pot(g.mac.mac_sub_idx,8);

        g.mac.mac_pot(g.mac.mac_sub_idx,7) = g.mac.mac_con(g.mac.mac_sub_idx,6) ...
                                             - g.mac.mac_con(g.mac.mac_sub_idx,7);

        g.mac.mac_pot(g.mac.mac_sub_idx,10) = ...
            (g.mac.mac_con(g.mac.mac_sub_idx,7) ...
             - g.mac.mac_con(g.mac.mac_sub_idx,8)) ...
            ./g.mac.mac_pot(g.mac.mac_sub_idx,8);

        g.mac.mac_pot(g.mac.mac_sub_idx,6) = ...
            g.mac.mac_pot(g.mac.mac_sub_idx,10)./g.mac.mac_pot(g.mac.mac_sub_idx,8);

        g.mac.mac_pot(g.mac.mac_sub_idx,13) = g.mac.mac_con(g.mac.mac_sub_idx,12) ...
                                              - g.mac.mac_con(g.mac.mac_sub_idx,4);

        g.mac.mac_pot(g.mac.mac_sub_idx,14) = ...
            (g.mac.mac_con(g.mac.mac_sub_idx,13) ...
             - g.mac.mac_con(g.mac.mac_sub_idx,4)) ...
            ./g.mac.mac_pot(g.mac.mac_sub_idx,13);

        g.mac.mac_pot(g.mac.mac_sub_idx,12) = g.mac.mac_con(g.mac.mac_sub_idx,11) ...
                                              - g.mac.mac_con(g.mac.mac_sub_idx,12);

        g.mac.mac_pot(g.mac.mac_sub_idx,15) = ...
            (g.mac.mac_con(g.mac.mac_sub_idx,12) ...
             - g.mac.mac_con(g.mac.mac_sub_idx,13)) ...
            ./g.mac.mac_pot(g.mac.mac_sub_idx,13);

        g.mac.mac_pot(g.mac.mac_sub_idx,11) = g.mac.mac_pot(g.mac.mac_sub_idx,15) ...
                                              ./g.mac.mac_pot(g.mac.mac_sub_idx,13);

        % eterm -- terminal bus voltage
        % theta -- bus angle in rad
        g.mac.eterm(g.mac.mac_sub_idx,1) = bus(busnum,2);
        g.bus.theta(busnum,1) = bus(busnum,3)*pi/180;

        % pelect -- electrical output power, active
        % qelect -- electrical output power, reactive
        g.mac.pelect(g.mac.mac_sub_idx,1) = ...
            bus(busnum,4).*g.mac.mac_con(g.mac.mac_sub_idx,22);

        g.mac.qelect(g.mac.mac_sub_idx,1) = ...
            bus(busnum,5).*g.mac.mac_con(g.mac.mac_sub_idx,23);

        % current magnitude on generator base
        curr = sqrt(g.mac.pelect(g.mac.mac_sub_idx,1).^2 ...
                    + g.mac.qelect(g.mac.mac_sub_idx,1).^2) ...
               ./g.mac.eterm(g.mac.mac_sub_idx,1) ...
               .*g.mac.mac_pot(g.mac.mac_sub_idx,1);

        % power factor angle (rad)
        phi = atan2(g.mac.qelect(g.mac.mac_sub_idx,1), ...
                    g.mac.pelect(g.mac.mac_sub_idx,1));

        % complex voltage in system reference frame
        v = g.mac.eterm(g.mac.mac_sub_idx,1).*exp(1j*g.bus.theta(busnum,1));

        % complex current in system reference frame
        curr = curr.*exp(1j*(g.bus.theta(busnum,1) - phi));

        % voltage behind sub-transient reactance in system frame
        ei = v + (g.mac.mac_con(g.mac.mac_sub_idx,5) ...
                  + 1j*g.mac.mac_con(g.mac.mac_sub_idx,11)).*curr;

        % mac_ang -- machine angle (delta)
        % mac_spd -- machine speed at steady state
        g.mac.mac_ang(g.mac.mac_sub_idx,1) = atan2(imag(ei),real(ei));
        g.mac.mac_spd(g.mac.mac_sub_idx,1) = ones(g.mac.n_sub,1);

        % system reference frame rotation to Park's frame
        rot = 1j*exp(-1j*g.mac.mac_ang(g.mac.mac_sub_idx,1));

        % current on generator base in Park's frame
        curr = curr.*rot;
        mcurmag = abs(curr);

        % mechanical power = electrical power + losses on generator base
        g.mac.pmech(g.mac.mac_sub_idx,1) = ...
            g.mac.pelect(g.mac.mac_sub_idx,1).*g.mac.mac_pot(g.mac.mac_sub_idx,1) ...
            + g.mac.mac_con(g.mac.mac_sub_idx,5).*(mcurmag.*mcurmag);

        % curdg -- d axis current on generator base
        % curqg -- q axis current on generator base
        g.mac.curdg(g.mac.mac_sub_idx,1) = real(curr);
        g.mac.curqg(g.mac.mac_sub_idx,1) = imag(curr);

        % curd -- d axis current on system base
        % curq -- q axis current on system base
        g.mac.curd(g.mac.mac_sub_idx,1) = ...
            real(curr)./g.mac.mac_pot(g.mac.mac_sub_idx,1);

        g.mac.curq(g.mac.mac_sub_idx,1) = ...
            imag(curr)./g.mac.mac_pot(g.mac.mac_sub_idx,1);

        % ed -- d axis voltage in Park's frame
        % eq -- q axis voltage in Park's frame
        v = v.*rot;
        g.mac.ed(g.mac.mac_sub_idx,1) = real(v);
        g.mac.eq(g.mac.mac_sub_idx,1) = imag(v);

        % q axis voltage behind resistance
        eqra = g.mac.eq(g.mac.mac_sub_idx,1) ...
               + g.mac.mac_con(g.mac.mac_sub_idx,5) ...
                 .*g.mac.curqg(g.mac.mac_sub_idx,1);

        g.mac.psidpp = ...
            eqra + g.mac.mac_con(g.mac.mac_sub_idx,8) ...
                   .*g.mac.curdg(g.mac.mac_sub_idx,1);

        g.mac.psikd(g.mac.mac_sub_idx,1) = ...
            eqra + g.mac.mac_con(g.mac.mac_sub_idx,4) ...
                   .*g.mac.curdg(g.mac.mac_sub_idx,1);

        g.mac.eqprime(g.mac.mac_sub_idx,1) = ...
            eqra + g.mac.mac_con(g.mac.mac_sub_idx,7) ...
                   .*g.mac.curdg(g.mac.mac_sub_idx,1);

        edra = -g.mac.ed(g.mac.mac_sub_idx,1) ...
               - g.mac.mac_con(g.mac.mac_sub_idx,5) ...
                 .*g.mac.curdg(g.mac.mac_sub_idx,1);

        g.mac.psiqpp = ...
            edra + g.mac.mac_con(g.mac.mac_sub_idx,13) ...
                   .*g.mac.curqg(g.mac.mac_sub_idx,1);

        g.mac.psikq(g.mac.mac_sub_idx,1) = ...
            edra + g.mac.mac_con(g.mac.mac_sub_idx,4) ...
                   .*g.mac.curqg(g.mac.mac_sub_idx,1);

        % this is the negative of Edprime in block diagram
        g.mac.edprime(g.mac.mac_sub_idx,1) = ...
            edra + g.mac.mac_con(g.mac.mac_sub_idx,12) ...
                   .*g.mac.curqg(g.mac.mac_sub_idx,1);

        % compute saturation
        inv_sat = inv([0.64 0.8 1; 1 1 1; 1.44 1.2 1]);
        b = [0.8*ones(g.mac.n_sub,1), ...
             ones(g.mac.n_sub,1) + g.mac.mac_con(g.mac.mac_sub_idx,20), ...
             1.2*(ones(g.mac.n_sub,1) + g.mac.mac_con(g.mac.mac_sub_idx,21))];

        g.mac.mac_pot(g.mac.mac_sub_idx,3) = b*inv_sat(1,:)';
        g.mac.mac_pot(g.mac.mac_sub_idx,4) = b*inv_sat(2,:)';
        g.mac.mac_pot(g.mac.mac_sub_idx,5) = b*inv_sat(3,:)';

        % No saturation for now
        E_Isat = g.mac.eqprime(g.mac.mac_sub_idx,1);

        % E_Isat = g.mac.mac_pot(g.mac.mac_sub_idx,3) ...
        %          .*g.mac.eqprime(g.mac.mac_sub_idx,1).^2 ...
        %          + g.mac.mac_pot(g.mac.mac_sub_idx,4) ...
        %            .*g.mac.eqprime(g.mac.mac_sub_idx,1) ...
        %          + g.mac.mac_pot(g.mac.mac_sub_idx,5);
        %
        % nosat_idx = find(g.mac.eqprime(g.mac.mac_sub_idx,1) < 0.8);
        % if ~isempty(nosat_idx)
        %     E_Isat(nosat_idx) = g.mac.eqprime(g.mac.mac_sub_idx(nosat_idx),1);
        % end

        Eqpe = g.mac.eqprime(g.mac.mac_sub_idx,1) ...
               - g.mac.psikd(g.mac.mac_sub_idx,1) ...
               - g.mac.mac_pot(g.mac.mac_sub_idx,8) ...
                 .*g.mac.curdg(g.mac.mac_sub_idx,1);

        g.mac.vex(g.mac.mac_sub_idx,1) = ...
            g.mac.eqprime(g.mac.mac_sub_idx,1) ...
            + g.mac.mac_pot(g.mac.mac_sub_idx,7) ...
              .*(g.mac.curdg(g.mac.mac_sub_idx,1) ...
                 + g.mac.mac_pot(g.mac.mac_sub_idx,6).*Eqpe);

        g.mac.fldcur(g.mac.mac_sub_idx,1) = g.mac.vex(g.mac.mac_sub_idx,1);

        % psi is in system base and is the voltage behind xpp
        g.mac.psi_re(g.mac.mac_sub_idx,1) = ...
            sin(g.mac.mac_ang(g.mac.mac_sub_idx,1)).*(-g.mac.psiqpp) ...
            + cos(g.mac.mac_ang(g.mac.mac_sub_idx,1)).*g.mac.psidpp;

        g.mac.psi_im(g.mac.mac_sub_idx,1) = ...
            -cos(g.mac.mac_ang(g.mac.mac_sub_idx,1)).*(-g.mac.psiqpp) ...
            + sin(g.mac.mac_ang(g.mac.mac_sub_idx,1)).*g.mac.psidpp;

        % end initialization
    end

    if (flag == 1)  % network interface computation
        % vectorized computation

        % wrt machine reference
        g.mac.mac_ang(g.mac.mac_sub_idx,k) = g.mac.mac_ang(g.mac.mac_sub_idx,k) ...
                                             - g.sys.mach_ref(k)*ones(g.mac.n_sub,1);

        g.mac.psidpp = g.mac.mac_pot(g.mac.mac_sub_idx,9) ...
                       .*g.mac.eqprime(g.mac.mac_sub_idx,k) ...
                       + g.mac.mac_pot(g.mac.mac_sub_idx,10) ...
                         .*g.mac.psikd(g.mac.mac_sub_idx,k);

        g.mac.psiqpp = g.mac.mac_pot(g.mac.mac_sub_idx,14) ...
                       .*g.mac.edprime(g.mac.mac_sub_idx,k) ...
                       + g.mac.mac_pot(g.mac.mac_sub_idx,15) ...
                         .*g.mac.psikq(g.mac.mac_sub_idx,k);

        g.mac.psi_re(g.mac.mac_sub_idx,k) = ...
            g.mac.mac_spd(g.mac.mac_sub_idx,k) ...
            .*(sin(g.mac.mac_ang(g.mac.mac_sub_idx,k)).*(-g.mac.psiqpp) ...
               + cos(g.mac.mac_ang(g.mac.mac_sub_idx,k)).*g.mac.psidpp);

        g.mac.psi_im(g.mac.mac_sub_idx,k) = ...
            g.mac.mac_spd(g.mac.mac_sub_idx,k) ...
            .*(-cos(g.mac.mac_ang(g.mac.mac_sub_idx,k)).*(-g.mac.psiqpp) ...
               + sin(g.mac.mac_ang(g.mac.mac_sub_idx,k)).*g.mac.psidpp);

        % end of interface
    end

    if (flag == 2 || flag == 3)  % generator dynamics calculation
        % vectorized computation

        g.mac.psiqpp = g.mac.mac_pot(g.mac.mac_sub_idx,14) ...
                       .*g.mac.edprime(g.mac.mac_sub_idx,k) ...
                       + g.mac.mac_pot(g.mac.mac_sub_idx,15) ...
                         .*g.mac.psikq(g.mac.mac_sub_idx,k);

        g.mac.psidpp = g.mac.mac_pot(g.mac.mac_sub_idx,9) ...
                       .*g.mac.eqprime(g.mac.mac_sub_idx,k) ...
                       + g.mac.mac_pot(g.mac.mac_sub_idx,10) ...
                         .*g.mac.psikd(g.mac.mac_sub_idx,k);

        g.mac.curd(g.mac.mac_sub_idx,k) = ...
            sin(g.mac.mac_ang(g.mac.mac_sub_idx,k)) ...
            .*g.mac.cur_re(g.mac.mac_sub_idx,k) ...
            - cos(g.mac.mac_ang(g.mac.mac_sub_idx,k)) ...
              .*g.mac.cur_im(g.mac.mac_sub_idx,k);

        g.mac.curq(g.mac.mac_sub_idx,k) = ...
            cos(g.mac.mac_ang(g.mac.mac_sub_idx,k)) ...
            .*g.mac.cur_re(g.mac.mac_sub_idx,k) ...
            + sin(g.mac.mac_ang(g.mac.mac_sub_idx,k)) ...
              .*g.mac.cur_im(g.mac.mac_sub_idx,k);

        g.mac.curdg(g.mac.mac_sub_idx,k) = g.mac.curd(g.mac.mac_sub_idx,k) ...
                                           .*g.mac.mac_pot(g.mac.mac_sub_idx,1);

        g.mac.curqg(g.mac.mac_sub_idx,k) = g.mac.curq(g.mac.mac_sub_idx,k) ...
                                           .*g.mac.mac_pot(g.mac.mac_sub_idx,1);

        mcurmag = abs(g.mac.curdg(g.mac.mac_sub_idx,k) ...
                      + 1j*g.mac.curqg(g.mac.mac_sub_idx,k));

        % No saturation for now
        E_Isat = g.mac.eqprime(g.mac.mac_sub_idx,k);

        % E_Isat = g.mac.mac_pot(g.mac.mac_sub_idx,3) ...
        %          .*g.mac.eqprime(g.mac.mac_sub_idx,k).^2 ...
        %          + g.mac.mac_pot(g.mac.mac_sub_idx,4) ...
        %            .*g.mac.eqprime(g.mac.mac_sub_idx,k) ...
        %          + g.mac.mac_pot(g.mac.mac_sub_idx,5);
        %
        % nosat_idx = find(g.mac.eqprime(g.mac.mac_sub_idx,1) < 0.8);
        % if ~isempty(nosat_idx)
        %     E_Isat(nosat_idx) = g.mac.eqprime(g.mac.mac_sub_idx(nosat_idx),k);
        % end

        Eqpe = g.mac.eqprime(g.mac.mac_sub_idx,k) ...
               - g.mac.psikd(g.mac.mac_sub_idx,k) ...
               - g.mac.mac_pot(g.mac.mac_sub_idx,8) ...
                 .*g.mac.curdg(g.mac.mac_sub_idx,k);

        g.mac.fldcur(g.mac.mac_sub_idx,k) = ...
            g.mac.eqprime(g.mac.mac_sub_idx,k) ...
            + g.mac.mac_pot(g.mac.mac_sub_idx,7) ...
              .*(g.mac.curdg(g.mac.mac_sub_idx,k) ...
                 + g.mac.mac_pot(g.mac.mac_sub_idx,6).*Eqpe);

        g.mac.deqprime(g.mac.mac_sub_idx,k) = ...
            (g.mac.vex(g.mac.mac_sub_idx,k) ...
             - g.mac.fldcur(g.mac.mac_sub_idx,k)) ...
            ./g.mac.mac_con(g.mac.mac_sub_idx,9);

        g.mac.dpsikd(g.mac.mac_sub_idx,k) = ...
            Eqpe./g.mac.mac_con(g.mac.mac_sub_idx,10);

        Edpe = g.mac.edprime(g.mac.mac_sub_idx,k) ...
               - g.mac.psikq(g.mac.mac_sub_idx,k) ...
               - g.mac.mac_pot(g.mac.mac_sub_idx,13) ...
                 .*g.mac.curqg(g.mac.mac_sub_idx,k);

        Hold = -g.mac.edprime(g.mac.mac_sub_idx,k) ...
               + g.mac.mac_pot(g.mac.mac_sub_idx,12) ...
                 .*(-g.mac.curqg(g.mac.mac_sub_idx,k) ...
                    - g.mac.mac_pot(g.mac.mac_sub_idx,11).*Edpe);

        g.mac.dedprime(g.mac.mac_sub_idx,k) = ...
            Hold./g.mac.mac_con(g.mac.mac_sub_idx,14);

        g.mac.dpsikq(g.mac.mac_sub_idx,k) = ...
            Edpe./g.mac.mac_con(g.mac.mac_sub_idx,15);

        % Multiply psiqpp by gen spd
        g.mac.ed(g.mac.mac_sub_idx,k) = ...
            -g.mac.mac_con(g.mac.mac_sub_idx,5).*g.mac.curdg(g.mac.mac_sub_idx,k) ...
            - (g.mac.mac_spd(g.mac.mac_sub_idx,k).*g.mac.psiqpp ...
               - g.mac.mac_con(g.mac.mac_sub_idx,13) ...
                 .*g.mac.curqg(g.mac.mac_sub_idx,k));

        % Multiply psidpp by gen spd
        g.mac.eq(g.mac.mac_sub_idx,k) = ...
            -g.mac.mac_con(g.mac.mac_sub_idx,5).*g.mac.curqg(g.mac.mac_sub_idx,k) ...
            + (g.mac.mac_spd(g.mac.mac_sub_idx,k).*g.mac.psidpp ...
               - g.mac.mac_con(g.mac.mac_sub_idx,8) ...
                 .*g.mac.curdg(g.mac.mac_sub_idx,k));

        g.mac.eterm(g.mac.mac_sub_idx,k) = sqrt(g.mac.ed(g.mac.mac_sub_idx,k).^2 ...
                                                + g.mac.eq(g.mac.mac_sub_idx,k).^2);

        % On system base
        g.mac.pelect(g.mac.mac_sub_idx,k) = ...
            g.mac.eq(g.mac.mac_sub_idx,k).*g.mac.curq(g.mac.mac_sub_idx,k) ...
            + g.mac.ed(g.mac.mac_sub_idx,k).*g.mac.curd(g.mac.mac_sub_idx,k);

        g.mac.qelect(g.mac.mac_sub_idx,k) = ...
            g.mac.eq(g.mac.mac_sub_idx,k).*g.mac.curd(g.mac.mac_sub_idx,k) ...
            - g.mac.ed(g.mac.mac_sub_idx,k).*g.mac.curq(g.mac.mac_sub_idx,k);

        g.mac.dmac_ang(g.mac.mac_sub_idx,k) = ...
            g.sys.basrad*(g.mac.mac_spd(g.mac.mac_sub_idx,k) - ones(g.mac.n_sub,1));

        % electrical torque on generator base
        % Te = g.mac.pelect(g.mac.mac_sub_idx,k) ...
        %      .*g.mac.mac_pot(g.mac.mac_sub_idx,1) ...
        %      + g.mac.mac_con(g.mac.mac_sub_idx,5).*mcurmag.*mcurmag;

        % electrical torque on generator base
        Te = g.mac.psidpp.*g.mac.curqg(g.mac.mac_sub_idx,k) ...
             - g.mac.psiqpp.*g.mac.curdg(g.mac.mac_sub_idx,k);

        g.mac.dmac_spd(g.mac.mac_sub_idx,k) = ...
            (g.mac.pmech(g.mac.mac_sub_idx,k)./g.mac.mac_spd(g.mac.mac_sub_idx,k) ...
             + g.mac.pm_sig(g.mac.mac_sub_idx,k) - Te ...
             - g.mac.mac_con(g.mac.mac_sub_idx,17) ...
               .*(g.mac.mac_spd(g.mac.mac_sub_idx,k) - ones(g.mac.n_sub,1))) ...
            ./(2*g.mac.mac_con(g.mac.mac_sub_idx,16));

        % end rate calculation
    end
elseif (g.mac.n_sub ~= 0 && i ~= 0)
    error('mac_sub: all calculations must be vectorized.');
end

end  % function end

% eof
