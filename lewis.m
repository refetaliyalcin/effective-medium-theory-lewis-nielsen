function [n_eff_m,k_eff_m] = lewis(f_v,n_pigment,k_pigment,n_medium,k_medium)
    % Ertas Uslu, M., Yalcin, R. A., Misirlioglu, I. B., & Sendur, K. (2019). Morphology induced spectral reflectance lineshapes in VO2 thin films. Journal of Applied Physics, 125(22), 223103. doi:10.1063/1.5086272 
    %https://doi.org/10.1063/1.5086272
    A_e=0.4; %depolarization factor 0.4 for sphere
    f_max=0.637; %, fmax is the maximum packing fraction that is 0.637 for random close pack
    ev_pigment_e_a_=n_pigment.^2-k_pigment.^2; % eq. 2.30
    ev_pigment_e_a__=n_pigment.*k_pigment*2;
    ev_pigment_e_b_=n_medium.^2-k_medium.^2;
    ev_pigment_e_b__=n_medium.*k_medium*2;
    e_b=ev_pigment_e_b_-ev_pigment_e_b__*1i;
    e_a=ev_pigment_e_a_-ev_pigment_e_a__*1i;
    psi=1+((1-f_max)/f_max^2)*f_v;
    A=1/A_e-1;
    B=(e_a./e_b-1)./(e_a./e_b+A);
    e_mg=e_b.*(1+A*B*f_v)./(1-B*psi*f_v);
    e_mg_=real(e_mg);
    e_mg__=-imag(e_mg);
    n_eff_m=sqrt(0.5*(e_mg_+sqrt(e_mg_.^2+e_mg__.^2)));
    k_eff_m=sqrt(0.5*(-e_mg_+sqrt(e_mg_.^2+e_mg__.^2)));