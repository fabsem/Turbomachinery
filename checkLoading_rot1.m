%% Check the loading condition for ROTOR 1
% --> Look at the DIFFUSION FACTOR
%link between the momentum thickness and the chord
DiffFact = @(w1,w2,w1t,w2t,sigma) (w1 - w2)/w1 + abs((w1t - w2t)/(2 * w1 * sigma));

HUB.DiffFact1 = DiffFact(HUB.w1,HUB.w2,HUB.w1t,HUB.w2t,HUB.sigma1);
MID.DiffFact1 = DiffFact(MID.w1,MID.w2,MID.w1t,MID.w2t,MID.sigma1);
TIP.DiffFact1 = DiffFact(TIP.w1,TIP.w2,TIP.w1t,TIP.w2t,TIP.sigma1);

HUB.theta_c1 = thetaC_Diff(HUB.DiffFact1);
MID.theta_c1 = thetaC_Diff(MID.DiffFact1);
TIP.theta_c1 = thetaC_Diff(TIP.DiffFact1);
