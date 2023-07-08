function mac_sub(i,k,bus,flag)
% Syntax: mac_sub(i,k,bus,flag)
%
% Purpose: voltage-behind-subtransient-reactance generator
%          model, with vectorized computation option, and saturation
%
% Note:    state variables include mac_ang, mac_spd, eqprime,
%          psikd, edprime, psikq  (eqprime = psifd, psikd = psikd, edprime
%          = psi1q, psikq = psi2q, the old state names were maintained to 
%          a direct swap with mac_sub.m)
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
% See Also: mac_em, mac_tra, mac_sub
%
% Algorithm: PSLF model from John Undrill with saturation

%-----------------------------------------------------------------------------%
% Version history
%
% Version:  3.0
% Purpose:  Implement genrou model of PSLF with Saturation
% Date:     October 2022
% Author:   Ron Hruban
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

        if ~isempty(uets_idx)   % check that x''d and x''q are equal
            g.mac.mac_con(g.mac.mac_sub_idx(uets_idx),13) = ...
                g.mac.mac_con(g.mac.mac_sub_idx(uets_idx),8);
            %
            wstr = '\nmac_sub: xqpp made equal to xdpp at generator %0.0f.';
            warning(sprintf(wstr,g.mac.mac_sub_idx(uets_idx)));
        end

        notp_idx = find(g.mac.mac_con(g.mac.mac_sub_idx,14) == 0);
        if ~isempty(notp_idx)   % set T'qo to 999 if it is set to 0
            g.mac.mac_con(g.mac.mac_sub_idx(notp_idx),14) = ...
                999.0*ones(length(notp_idx),1);
        end

        notpp_idx = find(g.mac.mac_con(g.mac.mac_sub_idx,15) == 0);
        if ~isempty(notpp_idx)  % set T''qo to 999 if it is set to 0
            g.mac.mac_con(g.mac.mac_sub_idx(notpp_idx),15) = ...
                999.0*ones(length(notpp_idx),1);

            % set x'q = x"q if T''qo is set to zero
            g.mac.mac_con(g.mac.mac_sub_idx(notpp_idx),12) = ...
                g.mac.mac_con(g.mac.mac_sub_idx(notpp_idx),13);
        end

        % mac_indx ensures saturation values col 20/21 are 0 if not given

        % busnum -- bus number vector
        busnum = g.bus.bus_int(g.mac.mac_con(g.mac.mac_sub_idx,2));

        % mac_pot(,1) -- scaled MVA base
        % mac_pot(,2) -- base kV

        % system mva base / machine mva base
        g.mac.mac_pot(g.mac.mac_sub_idx,1) = g.sys.basmva*ones(g.mac.n_sub,1) ...
                                             ./g.mac.mac_con(g.mac.mac_sub_idx,3);
        
        % base kV?
        g.mac.mac_pot(g.mac.mac_sub_idx,2) = ones(g.mac.n_sub,1);

        % x'd - xL
        g.mac.mac_pot(g.mac.mac_sub_idx,8) = g.mac.mac_con(g.mac.mac_sub_idx,7) ...
                                             - g.mac.mac_con(g.mac.mac_sub_idx,4);

        % (x''d - xL)/(x'd - xL)
        g.mac.mac_pot(g.mac.mac_sub_idx,9) = ...
            (g.mac.mac_con(g.mac.mac_sub_idx,8) ...
             - g.mac.mac_con(g.mac.mac_sub_idx,4)) ...
            ./g.mac.mac_pot(g.mac.mac_sub_idx,8);

        % xd - x'd
        g.mac.mac_pot(g.mac.mac_sub_idx,7) = g.mac.mac_con(g.mac.mac_sub_idx,6) ...
                                             - g.mac.mac_con(g.mac.mac_sub_idx,7);
        % (x'd - x''d)
        g.mac.mac_pot(g.mac.mac_sub_idx,17) = ...
            g.mac.mac_con(g.mac.mac_sub_idx,7) ...
             - g.mac.mac_con(g.mac.mac_sub_idx,8);
        
        % (x'd - x''d)/(x'd - xL)
        g.mac.mac_pot(g.mac.mac_sub_idx,10) = ...
            g.mac.mac_pot(g.mac.mac_sub_idx,17) ...
            ./g.mac.mac_pot(g.mac.mac_sub_idx,8);

        % (x'd - x''d)/(x'd - xL)^2
        g.mac.mac_pot(g.mac.mac_sub_idx,6) = ...
            g.mac.mac_pot(g.mac.mac_sub_idx,10)./g.mac.mac_pot(g.mac.mac_sub_idx,8);

        % x'q - xL
        g.mac.mac_pot(g.mac.mac_sub_idx,13) = g.mac.mac_con(g.mac.mac_sub_idx,12) ...
                                              - g.mac.mac_con(g.mac.mac_sub_idx,4);

        % (x''q - xL)/(x'q - xL)
        g.mac.mac_pot(g.mac.mac_sub_idx,14) = ...
            (g.mac.mac_con(g.mac.mac_sub_idx,13) ...
             - g.mac.mac_con(g.mac.mac_sub_idx,4)) ...
            ./g.mac.mac_pot(g.mac.mac_sub_idx,13);

        % xq - x'q
        g.mac.mac_pot(g.mac.mac_sub_idx,12) = g.mac.mac_con(g.mac.mac_sub_idx,11) ...
                                              - g.mac.mac_con(g.mac.mac_sub_idx,12);
        % (x'q - x''q)/(x'q - xL)
        g.mac.mac_pot(g.mac.mac_sub_idx,15) = ...
            (g.mac.mac_con(g.mac.mac_sub_idx,12) ...
             - g.mac.mac_con(g.mac.mac_sub_idx,13)) ...
            ./g.mac.mac_pot(g.mac.mac_sub_idx,13);

        % (x'q - x''q)/(x'q - xL)^2
        g.mac.mac_pot(g.mac.mac_sub_idx,11) = g.mac.mac_pot(g.mac.mac_sub_idx,15) ...
                                              ./g.mac.mac_pot(g.mac.mac_sub_idx,13);
        % (xq - xL)/(xd - xL)
        g.mac.mac_pot(g.mac.mac_sub_idx,16) = ...
            (g.mac.mac_con(g.mac.mac_sub_idx,11) - g.mac.mac_con(g.mac.mac_sub_idx,4))...
            ./(g.mac.mac_con(g.mac.mac_sub_idx,6) - g.mac.mac_con(g.mac.mac_sub_idx,4));

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

        % current magnitude on generator base,  S/V*Vbase?
        curr = sqrt(g.mac.pelect(g.mac.mac_sub_idx,1).^2 ...
                    + g.mac.qelect(g.mac.mac_sub_idx,1).^2) ...
               ./g.mac.eterm(g.mac.mac_sub_idx,1) ...
               .*g.mac.mac_pot(g.mac.mac_sub_idx,1);

        % power factor angle (rad)
        phi = atan2(g.mac.qelect(g.mac.mac_sub_idx,1), ...
                    g.mac.pelect(g.mac.mac_sub_idx,1));

        % complex voltage in system reference frame
        vsys = g.mac.eterm(g.mac.mac_sub_idx,1).*exp(1j*g.bus.theta(busnum,1));

        % complex current in system reference frame
        currsys = curr.*exp(1j*(g.bus.theta(busnum,1) - phi));

        % add ryans off nom init!!!!!!!!!
        g.mac.mac_spd(g.mac.mac_sub_idx,1) = ones(g.mac.n_sub,1);

        % voltage behind sub-transient reactance in system frame
        % ei = v + (ra + j xq) i
        ei = vsys + (g.mac.mac_con(g.mac.mac_sub_idx,5) ...
                  + 1j*g.mac.mac_con(g.mac.mac_sub_idx,11)).*curr;

        % mac_ang -- machine angle (delta)
        % mac_spd -- machine speed at steady state
        g.mac.mac_ang(g.mac.mac_sub_idx,1) = atan2(imag(ei),real(ei));
        

        % system reference frame rotation to Park's frame
        rot = 1j*exp(-1j*g.mac.mac_ang(g.mac.mac_sub_idx,1));

        % current on generator base in Park's frame
        curr = currsys.*rot;
        mcurmag = abs(currsys);

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
        v = vsys.*rot;
        g.mac.ed(g.mac.mac_sub_idx,1) = real(v);
        g.mac.eq(g.mac.mac_sub_idx,1) = imag(v);

        % begin block diagram solve using ed eq curdg curqg for initial
        % conditions withou sat

        % Ed known g.mac.ed(g.mac.mac_sub_idx,1)
        % Eq known g.mac.eq(g.mac.mac_sub_idx,1)
        % Id known g.mac.curdg(g.mac.mac_sub_idx,1)
        % Iq known g.mac.curqg(g.mac.mac_sub_idx,1)
        % psidpp
        g.mac.psidpp(g.mac.mac_sub_idx,1) = 1./g.mac.mac_spd(g.mac.mac_sub_idx,1)...
            .* (g.mac.eq(g.mac.mac_sub_idx,1) + ...
            g.mac.mac_con(g.mac.mac_sub_idx,5) .* g.mac.curqg(g.mac.mac_sub_idx,1) + ...
            g.mac.mac_con(g.mac.mac_sub_idx,8) .* g.mac.curdg(g.mac.mac_sub_idx,1));
        % psifd
        g.mac.eqprime(g.mac.mac_sub_idx,1) = ...
            (g.mac.psidpp(g.mac.mac_sub_idx,1) + ...
            g.mac.mac_pot(g.mac.mac_sub_idx,10).*g.mac.mac_pot(g.mac.mac_sub_idx,8).*g.mac.curdg(g.mac.mac_sub_idx,1))./...
            (g.mac.mac_pot(g.mac.mac_sub_idx,10) + g.mac.mac_pot(g.mac.mac_sub_idx,9));
        % psikd
        g.mac.psikd(g.mac.mac_sub_idx,1) = ...
            g.mac.eqprime(g.mac.mac_sub_idx,1) - ...
            g.mac.mac_pot(g.mac.mac_sub_idx,8).*g.mac.curdg(g.mac.mac_sub_idx,1);
        % Efd 
        g.mac.Efd(g.mac.mac_sub_idx,1) = ...
            g.mac.eqprime(g.mac.mac_sub_idx,1) + ...
            g.mac.mac_pot(g.mac.mac_sub_idx,7).*g.mac.curdg(g.mac.mac_sub_idx,1);
        % \delta known g.mac.mac_ang(g.mac.mac_sub_idx,1)
        % psi1q
        g.mac.edprime(g.mac.mac_sub_idx,1) = ...
            -g.mac.mac_pot(g.mac.mac_sub_idx,12).*g.mac.curqg(g.mac.mac_sub_idx,1);
        % psi2q
        g.mac.psikq(g.mac.mac_sub_idx,1) = ...
            g.mac.edprime(g.mac.mac_sub_idx,1) - ...
            g.mac.mac_pot(g.mac.mac_sub_idx,13).*g.mac.curqg(g.mac.mac_sub_idx,1);
        % psiqpp
        g.mac.psiqpp(g.mac.mac_sub_idx,1) = -1./g.mac.mac_spd(g.mac.mac_sub_idx,1)...
            .* (g.mac.ed(g.mac.mac_sub_idx,1) + ...
            g.mac.mac_con(g.mac.mac_sub_idx,5) .* g.mac.curdg(g.mac.mac_sub_idx,1) - ...
            g.mac.mac_con(g.mac.mac_sub_idx,13) .* g.mac.curqg(g.mac.mac_sub_idx,1));
        % psipp
        g.mac.psipp(g.mac.mac_sub_idx,1) = ....
            sqrt(g.mac.psidpp(g.mac.mac_sub_idx,1).^2 + g.mac.psiqpp(g.mac.mac_sub_idx,1).^2);
        
        % % Se
        % g.mac.Se(g.mac.mac_sub_idx,1) = 0.0; % TODO fix
        
%         % q axis voltage behind resistance
%         eqra = g.mac.eq(g.mac.mac_sub_idx,1) ...
%                + g.mac.mac_con(g.mac.mac_sub_idx,5) ...
%                  .*g.mac.curqg(g.mac.mac_sub_idx,1);
% 
%         g.mac.psidpp = ...
%             eqra + g.mac.mac_con(g.mac.mac_sub_idx,8) ...
%                    .*g.mac.curdg(g.mac.mac_sub_idx,1);
% 
%         g.mac.psikd(g.mac.mac_sub_idx,1) = ...
%             eqra + g.mac.mac_con(g.mac.mac_sub_idx,4) ...
%                    .*g.mac.curdg(g.mac.mac_sub_idx,1);
% 
%         g.mac.eqprime(g.mac.mac_sub_idx,1) = ...
%             eqra + g.mac.mac_con(g.mac.mac_sub_idx,7) ...
%                    .*g.mac.curdg(g.mac.mac_sub_idx,1);
% 
%         edra = -g.mac.ed(g.mac.mac_sub_idx,1) ...
%                - g.mac.mac_con(g.mac.mac_sub_idx,5) ...
%                  .*g.mac.curdg(g.mac.mac_sub_idx,1);
% 
%         g.mac.psiqpp = ...
%             edra + g.mac.mac_con(g.mac.mac_sub_idx,13) ...
%                    .*g.mac.curqg(g.mac.mac_sub_idx,1);
% 
%         g.mac.psikq(g.mac.mac_sub_idx,1) = ...
%             edra + g.mac.mac_con(g.mac.mac_sub_idx,4) ...
%                    .*g.mac.curqg(g.mac.mac_sub_idx,1);
% 
%         % this is the negative of Edprime in block diagram
%         g.mac.edprime(g.mac.mac_sub_idx,1) = ...
%             edra + g.mac.mac_con(g.mac.mac_sub_idx,12) ...
%                    .*g.mac.curqg(g.mac.mac_sub_idx,1);

        % compute saturation

        % solve for A and B
        % Quadratic
        g.mac.mac_pot(g.mac.mac_sub_idx,3) = sqrt(g.mac.mac_con(g.mac.mac_sub_idx,20)./g.mac.mac_con(g.mac.mac_sub_idx,21));
        g.mac.mac_pot(g.mac.mac_sub_idx,4) = (1.2.*g.mac.mac_pot(g.mac.mac_sub_idx,3)-1)./(g.mac.mac_pot(g.mac.mac_sub_idx,3)-1); %A
        g.mac.mac_pot(g.mac.mac_sub_idx,5) = g.mac.mac_con(g.mac.mac_sub_idx,20)./(1-g.mac.mac_pot(g.mac.mac_sub_idx,4)).^2; %B

        % get rid of nan and inf A and B cause by zeros for S1 and S12
        An_indx = find(isnan(g.mac.mac_pot(g.mac.mac_sub_idx,4)));
        g.mac.mac_pot(g.mac.mac_sub_idx(An_indx),4) = zeros(length(An_indx),1);
        Ai_indx = find(isinf(g.mac.mac_pot(g.mac.mac_sub_idx,4)));
        g.mac.mac_pot(g.mac.mac_sub_idx(Ai_indx),4) = zeros(length(Ai_indx),1);

        Bn_indx = find(isnan(g.mac.mac_pot(g.mac.mac_sub_idx,5)));
        g.mac.mac_pot(g.mac.mac_sub_idx(Bn_indx),5) = zeros(length(Bn_indx),1);
        Bi_indx = find(isinf(g.mac.mac_pot(g.mac.mac_sub_idx,5)));
        g.mac.mac_pot(g.mac.mac_sub_idx(Bi_indx),5) = zeros(length(Bi_indx),1);


%         g.mac.mac_pot(g.mac.mac_sub_idx,3) = sqrt(g.mac.mac_con(g.mac.mac_sub_idx,20)./g.mac.mac_con(g.mac.mac_sub_idx,21));
%         g.mac.mac_pot(g.mac.mac_sub_idx,4) = (1.2.*g.mac.mac_pot(g.mac.mac_sub_idx,3)-1)./(g.mac.mac_pot(g.mac.mac_sub_idx,3)-1);
%         g.mac.mac_pot(g.mac.mac_sub_idx,5) = g.mac.mac_con(g.mac.mac_sub_idx,20)./(1-g.mac.mac_pot(g.mac.mac_sub_idx,4)).^2;

        % Scaled Quadratic
%         g.mac.mac_pot(g.mac.mac_sub_idx,3) = sqrt(g.mac.mac_con(g.mac.mac_sub_idx,20)./g.mac.mac_con(g.mac.mac_sub_idx,21)./1.2);
%         g.mac.mac_pot(g.mac.mac_sub_idx,4) = (1.2.*g.mac.mac_pot(g.mac.mac_sub_idx,3)-1)./(g.mac.mac_pot(g.mac.mac_sub_idx,3)-1);
%         g.mac.mac_pot(g.mac.mac_sub_idx,5) = g.mac.mac_con(g.mac.mac_sub_idx,20)./(1-g.mac.mac_pot(g.mac.mac_sub_idx,4)).^2;

        % Exponential
%         g.mac.mac_pot(g.mac.mac_sub_idx,5) = g.mac.mac_con(g.mac.mac_sub_idx,20);
%         g.mac.mac_pot(g.mac.mac_sub_idx,4) = log(g.mac.mac_con(g.mac.mac_sub_idx,21)./g.mac.mac_con(g.mac.mac_sub_idx,20))./log(1.2);
        
        % zero A and B if psipp < A 
        lowPSIPP_indx = find(g.mac.psipp(g.mac.mac_sub_idx,1) < g.mac.mac_pot(g.mac.mac_sub_idx,4));
        g.mac.mac_pot(g.mac.mac_sub_idx(lowPSIPP_indx),5) = zeros(length(lowPSIPP_indx),1);

        % Se
        g.mac.Se(g.mac.mac_sub_idx,1) = g.mac.mac_pot(g.mac.mac_sub_idx,5).*(g.mac.psipp(g.mac.mac_sub_idx,1) - g.mac.mac_pot(g.mac.mac_sub_idx,4)).^2;

        % There needs to be code here that sets A and B to zero to if psipp
        % is less than A. This is just for initalization thought need to
        % find a way to change them back for later in the simulation. Also
        % later in the simulation the Se calculations need an if else
        % statement to capture when psipp is less than A.

        % newton raphson setup
        x0 = [...
            g.mac.ed(g.mac.mac_sub_idx,1);...      %1
            g.mac.eq(g.mac.mac_sub_idx,1);...      %2
            g.mac.curdg(g.mac.mac_sub_idx,1);...   %3
            g.mac.curqg(g.mac.mac_sub_idx,1);...   %4
            g.mac.Efd(g.mac.mac_sub_idx,1);...     %5
            g.mac.mac_ang(g.mac.mac_sub_idx,1);... %6
            g.mac.eqprime(g.mac.mac_sub_idx,1);... %7
            g.mac.psikd(g.mac.mac_sub_idx,1);...   %8
            g.mac.edprime(g.mac.mac_sub_idx,1);... %9
            g.mac.psikq(g.mac.mac_sub_idx,1);...   %10
            g.mac.psidpp(g.mac.mac_sub_idx,1);...  %11
            g.mac.psiqpp(g.mac.mac_sub_idx,1);...  %12
            g.mac.psipp(g.mac.mac_sub_idx,1);...   %13
            g.mac.Se(g.mac.mac_sub_idx,1);...      %14
            ];

        l = length(g.mac.mac_sub_idx);
        rb = @(n) ((n-1)*l+1:(n-1)*l+l);

        % this all needs to be sparse lxl
        zb = sparse(zeros(l));
        ib = sparse(eye(l));
        
        myfun = @(x) [...
            -x(rb(11)) + g.mac.mac_pot(g.mac.mac_sub_idx,9).*x(rb(7)) + g.mac.mac_pot(g.mac.mac_sub_idx,10).*x(rb(8));...
            -x(rb(12)) + g.mac.mac_pot(g.mac.mac_sub_idx,14).*x(rb(9)) + g.mac.mac_pot(g.mac.mac_sub_idx,15).*x(rb(10));...
            -x(rb(2)) + g.mac.mac_spd(g.mac.mac_sub_idx,1).*x(rb(11)) - g.mac.mac_con(g.mac.mac_sub_idx,5).*x(rb(4)) - g.mac.mac_con(g.mac.mac_sub_idx,8).*x(rb(3));...
            -x(rb(1)) - g.mac.mac_spd(g.mac.mac_sub_idx,1).*x(rb(12)) - g.mac.mac_con(g.mac.mac_sub_idx,5).*x(rb(3)) + g.mac.mac_con(g.mac.mac_sub_idx,13).*x(rb(4));...
            -x(rb(7)) - x(rb(14)).*x(rb(11)) - g.mac.mac_pot(g.mac.mac_sub_idx,7).*x(rb(3)) + x(rb(5));...
            -x(rb(8)) + x(rb(7)) - g.mac.mac_pot(g.mac.mac_sub_idx,8).*x(rb(3));...
            -x(rb(9)) - g.mac.mac_pot(g.mac.mac_sub_idx,16).*x(rb(14)).*x(rb(12)) - g.mac.mac_pot(g.mac.mac_sub_idx,12).*x(rb(4));...
            -x(rb(10)) + x(rb(9)) - g.mac.mac_pot(g.mac.mac_sub_idx,13).*x(rb(4));...
            x(rb(13)) - sqrt(x(rb(11)).^2 + x(rb(12)).^2);...
            x(rb(14)) - g.mac.mac_pot(g.mac.mac_sub_idx,5).*(x(rb(13))-g.mac.mac_pot(g.mac.mac_sub_idx,4)).^2;...   %this line can be changed for different saturation functions
            x(rb(1)).*sin(x(rb(6))) + x(rb(2)).*cos(x(rb(6)));...
            x(rb(2)).*sin(x(rb(6))) - x(rb(1)).*cos(x(rb(6)));...
            x(rb(3)).*sin(x(rb(6))) + x(rb(4)).*cos(x(rb(6)));...
            x(rb(4)).*sin(x(rb(6))) - x(rb(3)).*cos(x(rb(6)));...
            ];

        J1_11 = @(x) diag(sin(x(rb(6))));
        J1_12 = @(x) diag(-cos(x(rb(6))));
        J2_11 = @(x) diag(cos(x(rb(6))));
        J2_12 = @(x) diag(sin(x(rb(6))));
        J3_3 = diag(-g.mac.mac_con(g.mac.mac_sub_idx,8));
        J3_4 = diag(-g.mac.mac_con(g.mac.mac_sub_idx,5));
        J3_5 = diag(-g.mac.mac_pot(g.mac.mac_sub_idx,7));
        J3_6 = diag(-g.mac.mac_pot(g.mac.mac_sub_idx,8));
        J3_13 = @(x) diag(sin(x(rb(6))));
        J3_14 = @(x) diag(-cos(x(rb(6))));
        J4_3 = diag(-g.mac.mac_con(g.mac.mac_sub_idx,5));
        J4_4 = diag(g.mac.mac_con(g.mac.mac_sub_idx,13));
        J4_7 = diag(-g.mac.mac_pot(g.mac.mac_sub_idx,12));
        J4_8 = diag(-g.mac.mac_pot(g.mac.mac_sub_idx,13));
        J4_13 = @(x) diag(cos(x(rb(6))));
        J4_14 = @(x) diag(sin(x(rb(6))));
        J6_11 = @(x) diag((x(rb(1)).*cos(x(rb(6))) - x(rb(2)).*sin(x(rb(6)))));
        J6_12 = @(x) diag((x(rb(2)).*cos(x(rb(6))) + x(rb(1)).*sin(x(rb(6)))));
        J6_13 = @(x) diag((x(rb(3)).*cos(x(rb(6))) - x(rb(4)).*sin(x(rb(6)))));
        J6_14 = @(x) diag((x(rb(4)).*cos(x(rb(6))) + x(rb(3)).*sin(x(rb(6)))));
        J7_1 = diag(g.mac.mac_pot(g.mac.mac_sub_idx,9));
        J8_1 = diag(g.mac.mac_pot(g.mac.mac_sub_idx,10));
        J9_2 = diag(g.mac.mac_pot(g.mac.mac_sub_idx,14));
        J10_2 = diag(g.mac.mac_pot(g.mac.mac_sub_idx,15));
        J11_3 = diag(g.mac.mac_spd(g.mac.mac_sub_idx,1));
        J11_5 = @(x) diag(-x(rb(14)));
        J11_9 = @(x) diag((-x(rb(11))./sqrt(x(rb(11)).^2 + x(rb(12)).^2)));
        J12_4 = diag(-g.mac.mac_spd(g.mac.mac_sub_idx,1));
        J12_7 = @(x) diag((-x(rb(14)).*g.mac.mac_pot(g.mac.mac_sub_idx,16)));
        J12_9 = @(x) diag((-x(rb(12))./sqrt(x(rb(11)).^2 + x(rb(12)).^2)));
        J13_10 = @(x) diag((-2.*g.mac.mac_pot(g.mac.mac_sub_idx,5).*(x(rb(13)) - g.mac.mac_pot(g.mac.mac_sub_idx,4))));
        J14_5 = @(x) diag(-x(rb(11)));
        J14_7 = @(x) diag((-x(rb(12)).*g.mac.mac_pot(g.mac.mac_sub_idx,16)));

        J = @(x) [... 
            zb,       zb,       zb,       zb,       zb, zb,       J7_1,  J8_1,  zb,    zb,    -1*ib,    zb,       zb,        zb;...
            zb,       zb,       zb,       zb,       zb, zb,       zb,    zb,    J9_2,  J10_2, zb,       -1*ib,    zb,        zb;...
            zb,       -1*ib,    J3_3,     J4_3,     zb, zb,       zb,    zb,    zb,    zb,    J11_3,    zb,       zb,        zb;...
            -1*ib,    zb,       J3_4,     J4_4,     zb, zb,       zb,    zb,    zb,    zb,    zb,       J12_4,    zb,        zb;...
            zb,       zb,       J3_5,     zb,       ib, zb,       -1*ib, zb,    zb,    zb,    J11_5(x), zb,       zb,        J14_5(x);...
            zb,       zb,       J3_6,     zb,       zb, zb,       ib,    -1*ib, zb,    zb,    zb,       zb,       zb,        zb;...
            zb,       zb,       zb,       J4_7,     zb, zb,       zb,    zb,    -1*ib, zb,    zb,       J12_7(x), zb,        J14_7(x);...
            zb,       zb,       zb,       J4_8,     zb, zb,       zb,    zb,    ib,    -1*ib, zb,       zb,       zb,        zb;...
            zb,       zb,       zb,       zb,       zb, zb,       zb,    zb,    zb,    zb,    J11_9(x), J12_9(x), ib,        zb;...
            zb,       zb,       zb,       zb,       zb, zb,       zb,    zb,    zb,    zb,    zb,       zb,       J13_10(x), ib;...
            J1_11(x), J2_11(x), zb,       zb,       zb, J6_11(x), zb,    zb,    zb,    zb,    zb,       zb,       zb,        zb;...
            J1_12(x), J2_12(x), zb,       zb,       zb, J6_12(x), zb,    zb,    zb,    zb,    zb,       zb,       zb,        zb;...
            zb,       zb,       J3_13(x), J4_13(x), zb, J6_13(x), zb,    zb,    zb,    zb,    zb,       zb,       zb,        zb;...
            zb,       zb,       J3_14(x), J4_14(x), zb, J6_14(x), zb,    zb,    zb,    zb,    zb,       zb,       zb,        zb;...
            ];
        
        y = [...
            zeros(10*l,1);...
            real(vsys);...
            imag(vsys);...
            real(currsys);...
            imag(currsys)...
            ];

        % start iteration process for main Newton_Raphson solution
        iter = 0;  % initialize iteration counter
        itermax = 40;
        Tol = 0.00005;
        xhat = x0;
        conv_flag = 1;

        g.Jac_iter = J(xhat);

        while ((conv_flag == 1) && (iter <= itermax))
            iter = iter + 1;

            % form the Jacobian matrix
            Jac_iter = J(xhat);

            % mismatch vectors
            resid = y-myfun(xhat);

            % solve for increments
            delx = Jac_iter\resid;

            % update 
            xhat = xhat + delx;

            % check convergence
            if norm(delx) < Tol
                conv_flag = 0;
            end
        end

        if (iter > itermax)
            wstr = '\nGenrou: saturation initialization failed to converge after %0.0f ';
            warning(sprintf(wstr,[itermax]));
           
        else
            disp(sprintf('Genrou: saturation initialization iterations = %0.0f',iter));
        end

        % extract from xhat
        g.mac.ed(g.mac.mac_sub_idx,1) = xhat(rb(1),1);
        g.mac.eq(g.mac.mac_sub_idx,1) = xhat(rb(2),1);
        g.mac.curdg(g.mac.mac_sub_idx,1) = xhat(rb(3),1);
        g.mac.curqg(g.mac.mac_sub_idx,1) = xhat(rb(4),1);
        g.mac.Efd(g.mac.mac_sub_idx,1) = xhat(rb(5),1);
        g.mac.mac_ang(g.mac.mac_sub_idx,1) = xhat(rb(6),1);
        g.mac.eqprime(g.mac.mac_sub_idx,1) = xhat(rb(7),1);
        g.mac.psikd(g.mac.mac_sub_idx,1) = xhat(rb(8),1);
        g.mac.edprime(g.mac.mac_sub_idx,1) = xhat(rb(9),1);
        g.mac.psikq(g.mac.mac_sub_idx,1) = xhat(rb(10),1);
        g.mac.psidpp(g.mac.mac_sub_idx,1) = xhat(rb(11),1);
        g.mac.psiqpp(g.mac.mac_sub_idx,1) = xhat(rb(12),1);
        g.mac.psipp(g.mac.mac_sub_idx,1) = xhat(rb(13),1);
        g.mac.Se(g.mac.mac_sub_idx,1) = xhat(rb(14),1);

%         inv_sat = inv([0.64 0.8 1; 1 1 1; 1.44 1.2 1]);
%         b = [0.8*ones(g.mac.n_sub,1), ...
%              ones(g.mac.n_sub,1) + g.mac.mac_con(g.mac.mac_sub_idx,20), ...
%              1.2*(ones(g.mac.n_sub,1) + g.mac.mac_con(g.mac.mac_sub_idx,21))];
% 
%         g.mac.mac_pot(g.mac.mac_sub_idx,3) = b*inv_sat(1,:)';
%         g.mac.mac_pot(g.mac.mac_sub_idx,4) = b*inv_sat(2,:)';
%         g.mac.mac_pot(g.mac.mac_sub_idx,5) = b*inv_sat(3,:)';
% 
%         % No saturation for now
%         E_Isat = g.mac.eqprime(g.mac.mac_sub_idx,1);
% 
%         % E_Isat = g.mac.mac_pot(g.mac.mac_sub_idx,3) ...
%         %          .*g.mac.eqprime(g.mac.mac_sub_idx,1).^2 ...
%         %          + g.mac.mac_pot(g.mac.mac_sub_idx,4) ...
%         %            .*g.mac.eqprime(g.mac.mac_sub_idx,1) ...
%         %          + g.mac.mac_pot(g.mac.mac_sub_idx,5);
%         %
%         % nosat_idx = find(g.mac.eqprime(g.mac.mac_sub_idx,1) < 0.8);
%         % if ~isempty(nosat_idx)
%         %     E_Isat(nosat_idx) = g.mac.eqprime(g.mac.mac_sub_idx(nosat_idx),1);
%         % end
% 
%         Eqpe = g.mac.eqprime(g.mac.mac_sub_idx,1) ...
%                - g.mac.psikd(g.mac.mac_sub_idx,1) ...
%                - g.mac.mac_pot(g.mac.mac_sub_idx,8) ...
%                  .*g.mac.curdg(g.mac.mac_sub_idx,1);

        g.mac.vex(g.mac.mac_sub_idx,1) = ...
            g.mac.Efd(g.mac.mac_sub_idx,1);

        g.mac.fldcur(g.mac.mac_sub_idx,1) = g.mac.vex(g.mac.mac_sub_idx,1);

        % psi is in system base and is the voltage behind xpp (actually
        % psipp)
        g.mac.psi_re(g.mac.mac_sub_idx,1) = ...
            sin(g.mac.mac_ang(g.mac.mac_sub_idx,1)).*(-g.mac.psiqpp(g.mac.mac_sub_idx,1)) ...
            + cos(g.mac.mac_ang(g.mac.mac_sub_idx,1)).*g.mac.psidpp(g.mac.mac_sub_idx,1);

        g.mac.psi_im(g.mac.mac_sub_idx,1) = ...
            -cos(g.mac.mac_ang(g.mac.mac_sub_idx,1)).*(-g.mac.psiqpp(g.mac.mac_sub_idx,1)) ...
            + sin(g.mac.mac_ang(g.mac.mac_sub_idx,1)).*g.mac.psidpp(g.mac.mac_sub_idx,1);
        
        % Resolve for A and B (because some had to be zero'd for low psipp
        % init
        % Quadratic
        g.mac.mac_pot(g.mac.mac_sub_idx,3) = sqrt(g.mac.mac_con(g.mac.mac_sub_idx,20)./g.mac.mac_con(g.mac.mac_sub_idx,21));
        g.mac.mac_pot(g.mac.mac_sub_idx,4) = (1.2.*g.mac.mac_pot(g.mac.mac_sub_idx,3)-1)./(g.mac.mac_pot(g.mac.mac_sub_idx,3)-1); %A
        g.mac.mac_pot(g.mac.mac_sub_idx,5) = g.mac.mac_con(g.mac.mac_sub_idx,20)./(1-g.mac.mac_pot(g.mac.mac_sub_idx,4)).^2; %B

        % get rid of nan and inf A and B cause by zeros for S1 and S12
        An_indx = find(isnan(g.mac.mac_pot(g.mac.mac_sub_idx,4)));
        g.mac.mac_pot(g.mac.mac_sub_idx(An_indx),4) = zeros(length(An_indx),1);
        Ai_indx = find(isinf(g.mac.mac_pot(g.mac.mac_sub_idx,4)));
        g.mac.mac_pot(g.mac.mac_sub_idx(Ai_indx),4) = zeros(length(Ai_indx),1);

        Bn_indx = find(isnan(g.mac.mac_pot(g.mac.mac_sub_idx,5)));
        g.mac.mac_pot(g.mac.mac_sub_idx(Bn_indx),5) = zeros(length(Bn_indx),1);
        Bi_indx = find(isinf(g.mac.mac_pot(g.mac.mac_sub_idx,5)));
        g.mac.mac_pot(g.mac.mac_sub_idx(Bi_indx),5) = zeros(length(Bi_indx),1);

        % end initialization
    end

    if (flag == 1)  % network interface computation
        % vectorized computation

        % wrt machine reference
        g.mac.mac_ang(g.mac.mac_sub_idx,k) = g.mac.mac_ang(g.mac.mac_sub_idx,k) ...
                                             - g.sys.mach_ref(k)*ones(g.mac.n_sub,1);

        g.mac.psidpp(g.mac.mac_sub_idx,k) = g.mac.mac_pot(g.mac.mac_sub_idx,9) ...
                       .*g.mac.eqprime(g.mac.mac_sub_idx,k) ...
                       + g.mac.mac_pot(g.mac.mac_sub_idx,10) ...
                         .*g.mac.psikd(g.mac.mac_sub_idx,k);

        g.mac.psiqpp(g.mac.mac_sub_idx,k) = g.mac.mac_pot(g.mac.mac_sub_idx,14) ...
                       .*g.mac.edprime(g.mac.mac_sub_idx,k) ...
                       + g.mac.mac_pot(g.mac.mac_sub_idx,15) ...
                         .*g.mac.psikq(g.mac.mac_sub_idx,k);
        
        g.mac.psipp(g.mac.mac_sub_idx,k) = ...
            sqrt(g.mac.psidpp(g.mac.mac_sub_idx,k).^2 + g.mac.psiqpp(g.mac.mac_sub_idx,k).^2);
        
        g.mac.Se(g.mac.mac_sub_idx,k) = ...
            g.mac.mac_pot(g.mac.mac_sub_idx,5).*(g.mac.psipp(g.mac.mac_sub_idx,k) - g.mac.mac_pot(g.mac.mac_sub_idx,4)).^2;
        
        % zero Se if psipp < A 
        lowPSIPP_indx = find(g.mac.psipp(g.mac.mac_sub_idx,k) < g.mac.mac_pot(g.mac.mac_sub_idx,4));
        g.mac.Se(g.mac.mac_sub_idx(lowPSIPP_indx),k) = zeros(length(lowPSIPP_indx),1);

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
        
        pdkd = g.mac.eqprime(g.mac.mac_sub_idx,k)...
            - g.mac.psikd(g.mac.mac_sub_idx,k)...
            - g.mac.curdg(g.mac.mac_sub_idx,k).*g.mac.mac_pot(g.mac.mac_sub_idx,8);

        g.mac.fldcur(g.mac.mac_sub_idx,k) = ...
            g.mac.eqprime(g.mac.mac_sub_idx,k) ...
            + g.mac.mac_pot(g.mac.mac_sub_idx,7).*(g.mac.curdg(g.mac.mac_sub_idx,k) + g.mac.mac_pot(g.mac.mac_sub_idx,6).*pdkd)...
            + g.mac.psidpp(g.mac.mac_sub_idx,k).*g.mac.Se(g.mac.mac_sub_idx,k);


        g.mac.psi_re(g.mac.mac_sub_idx,k) = ...
            g.mac.mac_spd(g.mac.mac_sub_idx,k) ...
            .*(sin(g.mac.mac_ang(g.mac.mac_sub_idx,k)).*(-g.mac.psiqpp(g.mac.mac_sub_idx,k)) ...
               + cos(g.mac.mac_ang(g.mac.mac_sub_idx,k)).*g.mac.psidpp(g.mac.mac_sub_idx,k));

        g.mac.psi_im(g.mac.mac_sub_idx,k) = ...
            g.mac.mac_spd(g.mac.mac_sub_idx,k) ...
            .*(-cos(g.mac.mac_ang(g.mac.mac_sub_idx,k)).*(-g.mac.psiqpp(g.mac.mac_sub_idx,k)) ...
               + sin(g.mac.mac_ang(g.mac.mac_sub_idx,k)).*g.mac.psidpp(g.mac.mac_sub_idx,k));

        % end of interface
    end

    if (flag == 2 || flag == 3)  % generator dynamics calculation
        % vectorized computation

        g.mac.psiqpp(g.mac.mac_sub_idx,k) = g.mac.mac_pot(g.mac.mac_sub_idx,14) ...
                       .*g.mac.edprime(g.mac.mac_sub_idx,k) ...
                       + g.mac.mac_pot(g.mac.mac_sub_idx,15) ...
                         .*g.mac.psikq(g.mac.mac_sub_idx,k);

        g.mac.psidpp(g.mac.mac_sub_idx,k) = g.mac.mac_pot(g.mac.mac_sub_idx,9) ...
                       .*g.mac.eqprime(g.mac.mac_sub_idx,k) ...
                       + g.mac.mac_pot(g.mac.mac_sub_idx,10) ...
                         .*g.mac.psikd(g.mac.mac_sub_idx,k);

        g.mac.psipp(g.mac.mac_sub_idx,k) = ...
            sqrt(g.mac.psidpp(g.mac.mac_sub_idx,k).^2 + g.mac.psiqpp(g.mac.mac_sub_idx,k).^2);

        g.mac.Se(g.mac.mac_sub_idx,k) = ...
            g.mac.mac_pot(g.mac.mac_sub_idx,5).*(g.mac.psipp(g.mac.mac_sub_idx,k) - g.mac.mac_pot(g.mac.mac_sub_idx,4)).^2;
        
        % zero Se if psipp < A 
        lowPSIPP_indx = find(g.mac.psipp(g.mac.mac_sub_idx,k) < g.mac.mac_pot(g.mac.mac_sub_idx,4));
        g.mac.Se(g.mac.mac_sub_idx(lowPSIPP_indx),k) = zeros(length(lowPSIPP_indx),1);

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

%         mcurmag = abs(g.mac.curdg(g.mac.mac_sub_idx,k) ...
%                       + 1j*g.mac.curqg(g.mac.mac_sub_idx,k));
% 
%         % No saturation for now
%         E_Isat = g.mac.eqprime(g.mac.mac_sub_idx,k);

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

%         Eqpe = g.mac.eqprime(g.mac.mac_sub_idx,k) ...
%                - g.mac.psikd(g.mac.mac_sub_idx,k) ...
%                - g.mac.mac_pot(g.mac.mac_sub_idx,8) ...
%                  .*g.mac.curdg(g.mac.mac_sub_idx,k);
        
        % pre-dpsikd
        pdkd = g.mac.eqprime(g.mac.mac_sub_idx,k)...
            - g.mac.psikd(g.mac.mac_sub_idx,k)...
            - g.mac.curdg(g.mac.mac_sub_idx,k).*g.mac.mac_pot(g.mac.mac_sub_idx,8);
        
        g.mac.fldcur(g.mac.mac_sub_idx,k) = ...
            g.mac.eqprime(g.mac.mac_sub_idx,k) ...
            + g.mac.mac_pot(g.mac.mac_sub_idx,7).*(g.mac.curdg(g.mac.mac_sub_idx,k) + g.mac.mac_pot(g.mac.mac_sub_idx,6).*pdkd)...
            + g.mac.psidpp(g.mac.mac_sub_idx,k).*g.mac.Se(g.mac.mac_sub_idx,k);

        % dpsifd being called deqprime
        g.mac.deqprime(g.mac.mac_sub_idx,k) = ...
            (g.mac.vex(g.mac.mac_sub_idx,k) ...
             - g.mac.fldcur(g.mac.mac_sub_idx,k)) ...
            ./g.mac.mac_con(g.mac.mac_sub_idx,9);

        % dpsikd 
        g.mac.dpsikd(g.mac.mac_sub_idx,k) = ...
            pdkd./g.mac.mac_con(g.mac.mac_sub_idx,10);

%         Edpe = g.mac.edprime(g.mac.mac_sub_idx,k) ...
%                - g.mac.psikq(g.mac.mac_sub_idx,k) ...
%                - g.mac.mac_pot(g.mac.mac_sub_idx,13) ...
%                  .*g.mac.curqg(g.mac.mac_sub_idx,k);

%         Hold = -g.mac.edprime(g.mac.mac_sub_idx,k) ...
%                + g.mac.mac_pot(g.mac.mac_sub_idx,12) ...
%                  .*(-g.mac.curqg(g.mac.mac_sub_idx,k) ...
%                     - g.mac.mac_pot(g.mac.mac_sub_idx,11).*Edpe);

        % pre-dpsi2q
        pd2q = g.mac.edprime(g.mac.mac_sub_idx,k)...
            - g.mac.psikq(g.mac.mac_sub_idx,k)...
            - g.mac.curqg(g.mac.mac_sub_idx,k).*g.mac.mac_pot(g.mac.mac_sub_idx,13);

        % dpsi1q being called dedprime
        g.mac.dedprime(g.mac.mac_sub_idx,k) = ...
            -(g.mac.edprime(g.mac.mac_sub_idx,k)...
            + g.mac.mac_pot(g.mac.mac_sub_idx,12).*(g.mac.curqg(g.mac.mac_sub_idx,k) + g.mac.mac_pot(g.mac.mac_sub_idx,11).*pd2q)...
            + g.mac.mac_pot(g.mac.mac_sub_idx,16).*g.mac.Se(g.mac.mac_sub_idx,k).*g.mac.psiqpp(g.mac.mac_sub_idx,k))...
            ./g.mac.mac_con(g.mac.mac_sub_idx,14);

        % dpsi2q being called dpsikq
        g.mac.dpsikq(g.mac.mac_sub_idx,k) = ...
            pd2q./g.mac.mac_con(g.mac.mac_sub_idx,15);

        % Ed
        g.mac.ed(g.mac.mac_sub_idx,k) = ...
            -g.mac.mac_con(g.mac.mac_sub_idx,5).*g.mac.curdg(g.mac.mac_sub_idx,k) ...
            - g.mac.mac_spd(g.mac.mac_sub_idx,k).*g.mac.psiqpp(g.mac.mac_sub_idx,k) ...
            + g.mac.mac_con(g.mac.mac_sub_idx,13).*g.mac.curqg(g.mac.mac_sub_idx,k);

        % Eq
        g.mac.eq(g.mac.mac_sub_idx,k) = ...
            -g.mac.mac_con(g.mac.mac_sub_idx,5).*g.mac.curqg(g.mac.mac_sub_idx,k) ...
            + g.mac.mac_spd(g.mac.mac_sub_idx,k).*g.mac.psidpp(g.mac.mac_sub_idx,k) ...
            - g.mac.mac_con(g.mac.mac_sub_idx,8).*g.mac.curdg(g.mac.mac_sub_idx,k);

        g.mac.eterm(g.mac.mac_sub_idx,k) = sqrt(g.mac.ed(g.mac.mac_sub_idx,k).^2 ...
                                                + g.mac.eq(g.mac.mac_sub_idx,k).^2);

        % On system base P and Q
        g.mac.pelect(g.mac.mac_sub_idx,k) = ...
            g.mac.eq(g.mac.mac_sub_idx,k).*g.mac.curq(g.mac.mac_sub_idx,k) ...
            + g.mac.ed(g.mac.mac_sub_idx,k).*g.mac.curd(g.mac.mac_sub_idx,k);

        g.mac.qelect(g.mac.mac_sub_idx,k) = ...
            g.mac.eq(g.mac.mac_sub_idx,k).*g.mac.curd(g.mac.mac_sub_idx,k) ...
            - g.mac.ed(g.mac.mac_sub_idx,k).*g.mac.curq(g.mac.mac_sub_idx,k);

        % delta_dot
        g.mac.dmac_ang(g.mac.mac_sub_idx,k) = ...
            g.sys.basrad*(g.mac.mac_spd(g.mac.mac_sub_idx,k) - ones(g.mac.n_sub,1));

        % electrical torque on generator base
        Te = g.mac.psidpp(g.mac.mac_sub_idx,k).*g.mac.curqg(g.mac.mac_sub_idx,k) ...
             - g.mac.psiqpp(g.mac.mac_sub_idx,k).*g.mac.curdg(g.mac.mac_sub_idx,k);

        % w_dot
        g.mac.dmac_spd(g.mac.mac_sub_idx,k) = ...
            ((g.mac.pmech(g.mac.mac_sub_idx,k) ...
            + g.mac.pm_sig(g.mac.mac_sub_idx,k)  ...
            - g.mac.mac_con(g.mac.mac_sub_idx,17) ...
            .*(g.mac.mac_spd(g.mac.mac_sub_idx,k) - ones(g.mac.n_sub,1))) ...
            ./g.mac.mac_spd(g.mac.mac_sub_idx,k)...
            - Te) ...
            ./(2*g.mac.mac_con(g.mac.mac_sub_idx,16));

        % end rate calculation
    end
elseif (g.mac.n_sub ~= 0 && i ~= 0)
    error('mac_sub: all calculations must be vectorized.');
end

end  % function end

% eof
