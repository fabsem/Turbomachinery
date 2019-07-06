function [deltaCDs] = howell_SecondaryLosses(lambdaS,Cl)

  %HOWELL_SECONDARYLOSSES Estimate Secondary losses using Howell Correlation
  %
  % Input data required:
  %
  % lambdaS --> fixed parameter, function of Re
  % Cl --> lift coefficient
  %
  % Example:
  %
  % deltaCDa = howell_SecondaryLosses(lambdaS,Cl)

if isempty(lambdaS)
  lambdaS = 0.015 %or 0.02
end

  deltaCDs = lambdaS * Cl^2;

end
