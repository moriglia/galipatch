function [L_geom, W, EpsilonEff] = rectagularResonantPatchParams(wavelength, EpsilonR, h_sub)
    W = wavelength./2*sqrt(2/(EpsilonR +1));
    EpsilonEff = (EpsilonR + 1)/2 + ((EpsilonR - 1)/2)./sqrt(1 + 12*h_sub./W);
    deltaL = h_sub*0.412*(EpsilonEff+0.3).*(W./h_sub+0.264)/((EpsilonEff-0.258).*(W./h_sub+0.8));
    lambdaEff = wavelength./sqrt(EpsilonEff);
    L_geom = lambdaEff/2 - 2*deltaL;
end