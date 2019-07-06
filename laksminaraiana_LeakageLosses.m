function [deltaCDdelta]=laksminaraiana_LeakageLosses(delta,h,Cl)

  %LAKSMINARAIANA_LEAKAGELOSSES Estimate Leakage losses using Laksminaraiana Correlation
  %
  % Input data required:
  %
  % delta --> leakage height
  % h --> blade height
  % Cl --> lift coefficient
  %
  % Example:
  %
  % deltaCDdelta = laksminaraiana_LeakageLosses(delta,h,Cl)

  deltaCDdelta = 0.7 * delta / h * Cl^2;

end
