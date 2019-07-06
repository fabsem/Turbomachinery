function [deltaCDdelta] = vavra_LeakageLosses(delta,h,c,s,alpha2,Cl)

  %VAVRA_LEAKAGELOSSES Estimate Leakage losses using Vavra Correlation
  %
  % Input data required:
  %
  % delta --> leakage height
  % h --> blade height
  % c --> chord
  % s --> pitch
  % alpha2 -->
  % Cl --> lift coefficient
  %
  % Example:
  %
  % deltaCDdelta = vavra_LeakageLosses(delta,h,c,s,alpha2,Cl)

  deltaCDdelta = 0.25 * (delta * c * Cl^2) / (h * s * sind(alpha2));

end
