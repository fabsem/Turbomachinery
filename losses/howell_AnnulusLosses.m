function [deltaCDa] = howell_AnnulusLosses(s,h)

  %HOWELL_ANNULUSLOSSES Estimate annulus losses using Howell Correlation
  %
  % Input data required:
  %
  % s --> pitch
  % h --> blade height
  %
  % Example:
  %
  % deltaCDa = howell_AnnulusLosses(s,h)

  deltaCDa = 0.020 * s / h;

end
