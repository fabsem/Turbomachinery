%% Check the loading condition for ROTOR 2
% --> Look at the DIFFUSION FACTOR
%link between the momentum thickness and the chord
DiffFact = @(w1,w2,w1t,w2t,sigma) (w1 - w2)/w1 + abs((w1t - w2t)/(2 * w1 * sigma));

HUB.DiffFact2 = DiffFact(HUB.w3,HUB.w4,HUB.w3t,HUB.w4t,HUB.sigma2);
MID.DiffFact2 = DiffFact(MID.w3,MID.w4,MID.w3t,MID.w4t,MID.sigma2);
TIP.DiffFact2 = DiffFact(TIP.w3,TIP.w4,TIP.w3t,TIP.w4t,TIP.sigma2);

HUB.theta_c2 = thetaC_Diff(HUB.DiffFact2);
MID.theta_c2 = thetaC_Diff(MID.DiffFact2);
TIP.theta_c2 = thetaC_Diff(TIP.DiffFact2);

HUB.Yprofile2 = profileLossesLieblein(HUB.theta_c2,HUB.sigma2,HUB.beta3,HUB.beta4);
MID.Yprofile2 = profileLossesLieblein(MID.theta_c2,MID.sigma2,MID.beta3,MID.beta4);
TIP.Yprofile2 = profileLossesLieblein(TIP.theta_c2,TIP.sigma2,TIP.beta3,TIP.beta4);
