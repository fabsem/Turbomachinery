function [deltaCDa]=vavra_AnulusLosses(c,h)

  %VAVRA_ANNULUSLOSSES Estimate annulus losses using Howell Correlation
  %
  % Input data required:
  %
  % c --> chord
  % h --> blade height
  %
  % Example:
  %
  % deltaCDa = vavra_AnnulusLosses(c,h)

  deltaCDa = 0.02 * c / h;

end
