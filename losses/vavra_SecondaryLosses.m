function [deltaCDs]=vavra_SecondaryLosses(c,h,Cl)

  %VAVRA_SECONDARYLOSSES Estimate Secondary losses using Howell Correlation
  %
  % Input data required:
  %
  % h --> blade height
  % c --> chord
  % Cl --> lift coefficient
  %
  % Example:
  %
  % deltaCDs = vavra_SecondaryLosses(c,h, Cl)

  deltaCDs = 0.04 * c / h * Cl^2;

end 
