%% Check the loading condition for ROTOR 1
% --> Look at the DIFFUSION FACTOR
%link between the momentum thickness and the chord
DiffFact = @(v1,v2,v1t,v2t,sigma) (v1 - v2)/v1 + abs((v1t - v2t)/(2 * v1 * sigma));

HUBfp.DiffFact_IGV = DiffFact(HUBfp.v1,HUBfp.v2,HUBfp.v1t,HUBfp.v2t,HUBfp.sigma_IGV);
MIDfp.DiffFact_IGV = DiffFact(MIDfp.v1,MIDfp.v2,MIDfp.v1t,MIDfp.v2t,MIDfp.sigma_IGV);
TIPfp.DiffFact_IGV = DiffFact(TIPfp.v1,TIPfp.v2,TIPfp.v1t,TIPfp.v2t,TIPfp.sigma_IGV);

HUBfp.theta_c_IGV = thetaC_Diff(HUBfp.DiffFact_IGV);
MIDfp.theta_c_IGV = thetaC_Diff(MIDfp.DiffFact_IGV);
TIPfp.theta_c_IGV = thetaC_Diff(TIPfp.DiffFact_IGV);
