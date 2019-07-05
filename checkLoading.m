%% Check the loading condition
% --> Look at the DIFFUSION FACTOR
%link between the momentum thickness and the chord
DiffFact = @(w1,w2,w1t,w2t,sigma) (w1 - w2)/w1 + abs((w1t - w2t)/(2 * w1 * sigma));

DiffFact_hub = DiffFact(w1_hub,w2_hub,w1t_hub,w2t_hub,sigma_hub);
DiffFact_mid = DiffFact(w1_mid,w2_mid,w1t_mid,w2t_mid,sigma_mid);
DiffFact_tip = DiffFact(w1_tip,w2_tip,w1t_tip,w2t_tip,sigma_tip);

theta_c_hub = thetaC_Diff(DiffFact_hub);
theta_c_mid = thetaC_Diff(DiffFact_mid);
theta_c_tip = thetaC_Diff(DiffFact_tip);

Yprofile_hub = profileLossesLieblein(theta_c_hub,sigma_hub,beta1_hub,beta2_hub);
Yprofile_mid = profileLossesLieblein(theta_c_mid,sigma_mid,beta1_mid,beta2_mid);
Yprofile_tip = profileLossesLieblein(theta_c_tip,sigma_tip,beta1_tip,beta2_tip);
