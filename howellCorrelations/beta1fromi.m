function [HUB_beta,MID_beta,TIP_beta] = beta1fromi(HUB,MID,TIP,rotor)

options = optimset('Display','off');

if rotor == 1

  % ROT1 - HUB
  fun = @(x) root3ibeta(x,HUB.sigma1);
  x0 = [4,-0.2,40];
  x = fsolve(fun,x0,options);clc;
  HUB_beta = -x(3);

  % ROT1 - MID
  fun = @(x) root3ibeta(x,MID.sigma1);
  x0 = [4,-0.2,40];
  x = fsolve(fun,x0,options);clc;
  MID_beta = -x(3);

  % ROT1 - TIP
  fun = @(x) root3ibeta(x,TIP.sigma1);
  x0 = [4,-0.2,40];
  x = fsolve(fun,x0,options);clc;
  TIP_beta = -x(3);

  elseif rotor == 2

    % ROT2 - HUB
    fun = @(x) root3ibeta(x,HUB.sigma2);
    x0 = [4,-0.2,40];
    x = fsolve(fun,x0,options);clc;
    HUB_beta = x(3);

% ROT2 - HUB
    fun = @(x) root3ibeta(x,MID.sigma2);
    x0 = [4,-0.2,40];
    x = fsolve(fun,x0,options);clc;
    MID_beta = x(3);

    % ROT2 - HUB
    fun = @(x) root3ibeta(x,TIP.sigma2);
    x0 = [4,-0.2,40];
    x = fsolve(fun,x0,options);clc;
    TIP_beta = x(3);

    else
      disp('No rotor selected')
      return
    end

end
