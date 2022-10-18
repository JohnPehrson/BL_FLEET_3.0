function [u_plus,u_vd_plus,y_plus,u_plus_c,u_vd_plus_c,y_plus_c] = InnerVariableCalculator(u,y,Mach,T_inf,P_inf,U_inf)
%This function calculates the inner variables given a number of inputs to
%provide Van Driest corrected results. 
%Assumes adiabatic walls, ideal gas, and measurements close enough to the
%wall to directly measure the skin friction.

%% Invariants
Pr = 0.73;
r = Pr^(1/3); %recovery percentage
gamma = 1.4;
R = 287; %J/kgK, gas constant
%sutherland constants
S = 110.4; %K
T_ref = 273.15; %K
mu_ref = 1.716*(10^-5); %kg/ms

%change y to mm
y = y./1000;

%% Wall temperature
T_w = (1+r*((gamma-1)/2).*(Mach^2)).*T_inf;

%% Effective velocity using the Van Driest Scaling
a = sqrt(1-T_inf/T_w);
u_vd = (U_inf/a).*asin((a.*u)./U_inf);

%% Wall viscosities
rho_w = P_inf./(T_w.*R);
mu_w = mu_ref.*((T_w/T_ref).^1.5).*((T_ref+S)/(T_w+S));
v_w = mu_w./rho_w;

%% Wall shear stress and u_tau
    %fit the near-wall velocity gradient
    fitbelow = 0.0005;
    fitpoints = y<fitbelow;
    p = polyfit(y(fitpoints),u_vd(fitpoints),1);
    x_fit = linspace(0,fitbelow);
    y_fit = polyval(p,x_fit);
    dudy_w = p(1);
    t_w = mu_w.*dudy_w;
        t_w_clauser = 17;
    u_tau = sqrt(t_w/rho_w);
    u_tau_c = sqrt(t_w_clauser/rho_w);

    figure;
    plot(y,u_vd);
    xlabel('y (height)');
    ylabel('u (velocity)')
    hold on;
    plot(x_fit,y_fit);

%% Inner Variables
u_plus = u./u_tau;
u_vd_plus = u_vd./u_tau;
y_plus = (y.*u_tau)./v_w;

u_plus_c = u./u_tau_c;
u_vd_plus_c = u_vd./u_tau_c;
y_plus_c = (y.*u_tau_c)./v_w;
end