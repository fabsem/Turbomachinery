function [deltaBeta_optimal] = howellCorrelation(beta2, Re, sigma, varargin)

% HOWELLCORRELATION Returns the optimal deflection angle using Howell Correlation by Fabio Semeraro
%
% Input Parameters: beta2, Re and sigma (solidity).
%
% If no Re and no sigma provided, PSI and PHI are set = 1.
% This condition corresponds to Re = 3e5 and sigma = 1.
%
% Activating the flag 'plot', a plot of the current point is shown
% Activating the flag 'savefig', a figure is saved in the working directory
%
% Example:
%
%   [deltaBeta_optimal] = howellCorrelation(beta2,Re,sigma,'plot')


beta2 = abs(beta2);

%Plot Font Size and Style
%set(0,'defaultAxesFontSize',18)
%set(0,'defaultAxesTickLabelInterpreter','latex')
%set(0,'defaultLabelInterpreter','latex')
%set(0,'defaultTextInterpreter','latex')



%% Digital Howel Correlation

% x = [0.0753051155544071;
%      1.72942092962867;
%      3.94313165411582;
%      6.90167056175885;
%      11.2585475634034;
%      15.7967627456072;
%      20.6989526529906;
%      25.4042240110794;
%      29.6412187310655;
%      33.3212152687614;
%      37.0912317147062;
%      40.5846966155977;
%      44.9952393317753;
%      49.5845235003895];
%
% y = [32.8901584004155;
%      32.5149311867048;
%      31.7144897429239;
%      30.3944862806198;
%      28.0351423872587;
%      25.7236215701549;
%      23.4603566173288;
%      21.7179520470873;
%      20.4009781009262;
%      19.4144378083615;
%      18.4755041980438;
%      17.6307019821691;
%      16.6458928416861;
%      15.8036873539340];

% Degree of the interpolant polynomial
% n = 5;

% Coefficients of the interpolant polynomial
% For efficiency they're not calculated everytime!
% The coefficients for n = 5 are the defaul values provided.

%p = polyfit(x,y,n);
p = [2.58221652243903e-07,...
    -3.80609714991578e-05,...
     0.00202908260060230,...
    -0.0422142600384826,...
    -0.161326379175734,...
    32.9049224700194];

plot_check = 0;
if plot_check

x_howell = linspace(0,50);
y_howell = polyval(p,x_howell);

figure(1)
plot(x,y,'*')
hold on
plot(x_howell,y_howell,'linewidth',2)
xlabel('\beta_2 [°]')
ylabel('\Delta\beta^* [°]')
xlim([0,50])
ylim([10,40])
end


if isempty(Re) && isempty(sigma)
    Re = 3e5;
    sigma = 1;

    PSI = 1;
    PHI = 1;

else
    [PHI] = howellCorrectionPHI(Re);
    [PSI] = howellCorrectionPSI(sigma);

end

deltaBeta_optimal_p = polyval(p,beta2);
deltaBeta_optimal = deltaBeta_optimal_p * (PHI * PSI);


if strcmp(varargin,'plot')

    x_howell = linspace(0,50);
    y_howell = polyval(p,x_howell);

    figure(1)
    plot(x_howell,y_howell,'linewidth',2)
    hold on
    plot(ones(1,100) * beta2,linspace(0,deltaBeta_optimal_p),'r','linewidth',2)
    plot(beta2,10,'.r','markersize',30)
    plot(linspace(0,beta2),ones(1,100) * deltaBeta_optimal_p,'r','linewidth',2)
    plot(0,deltaBeta_optimal_p,'.r','markersize',30)
    beta2value = num2str(beta2);
    beta_optimalValue = num2str(round(deltaBeta_optimal_p,2));
    text(beta2+2,12,beta2value,'fontsize',18)
    text(1,deltaBeta_optimal_p + 2,beta_optimalValue,'fontsize',18)
    xlabel('$\beta_2 [^\circ]$')
    ylabel('$\Delta \beta^* [^\circ]$')
    xlim([0,50])
    ylim([10,40])
    grid on

elseif strcmp(varargin,'savefig')

    x_howell = linspace(0,50);
    y_howell = polyval(p,x_howell);

    figure(5)
    plot(x_howell,y_howell,'linewidth',2)
    hold on
    plot(ones(1,100) * beta2,linspace(0,deltaBeta_optimal),'r','linewidth',2)
    plot(beta2,10,'.r','markersize',30)
    plot(linspace(0,beta2),ones(1,100) * deltaBeta_optimal,'r','linewidth',2)
    plot(0,deltaBeta_optimal,'.r','markersize',30)
    beta2value = num2str(beta2);
    beta_optimalValue = num2str(round(deltaBeta_optimal,2));
    text(beta2+2,12,beta2value,'fontsize',18)
    text(1,deltaBeta_optimal + 2,beta_optimalValue,'fontsize',18)
    xlabel('$\beta_2 [^\circ]$')
    ylabel('$\Delta \beta^* [^\circ]$')
    xlim([0,50])
    ylim([10,40])
    grid on

    saveas(figure(5),'howellCorrelation.png');

end




end


function [PHI] = howellCorrectionPHI(Re)

    Re = Re/1e5;

    % Points
    x = [1.01526717557252;
         1.03816793893130;
         1.06870229007634;
         1.10687022900763;
         1.19083969465649;
         1.32824427480916;
         1.56488549618321;
         1.87022900763359;
         2.07633587786260;
         2.39694656488550;
         2.80916030534351;
         3.24427480916031;
         3.71755725190840;
         4.16030534351145;
         4.52671755725191;
         4.74045801526718;
         4.89312977099237;
         4.98473282442748];

    y = [0.842019543973942;
         0.841042345276873;
         0.840553745928339;
         0.841530944625407;
         0.844462540716612;
         0.855211726384365;
         0.881596091205212;
         0.920684039087948;
         0.944136807817590;
         0.973452768729642;
         0.996905537459284;
         1.00863192182410;
         1.01547231270358;
         1.01889250814332;
         1.01840390879479;
         1.01840390879479;
         1.01889250814332;
         1.01938110749186];

% Degree of the interpolant polynomial
% n = 8;

% Coefficients of the interpolant polynomial
% For efficiency they're not calculated everytime!
% The coefficients for n = 8 are the defaul values provided.

%p = polyfit(x,y,n);

    p = [9.42494984400196e-06,...
        -0.000198309453450787,...
         0.00284101529994087,...
        -0.0308406429566715,...
         0.210311660958647,...
        -0.833055689414330,...
         1.81082973735231,...
        -1.86215614529862,...
        1.54441412997659];


    plotcheck = 0;
    if plotcheck

    x_lin = linspace(1,5);
    y_lin = linspace(0.8,1.1);

    corr = polyval(p,x_lin);

    plot(x,y,'*')
    hold on
    plot(x_lin,corr)
    ylim([0.8,1.1])
    xlim([1,5])
    end

    PHI = polyval(p,Re);
end


function [PSI] = howellCorrectionPSI(sigma)


    % Points
%     x = [1.01526717557252;
%          1.03816793893130;
%          1.06870229007634;
%          1.10687022900763;
%          1.19083969465649;
%          1.32824427480916;
%          1.56488549618321;
%          1.87022900763359;
%          2.07633587786260;
%          2.39694656488550;
%          2.80916030534351;
%          3.24427480916031;
%          3.71755725190840;
%          4.16030534351145;
%          4.52671755725191;
%          4.74045801526718;
%          4.89312977099237;
%          4.98473282442748];
%
%     y = [0.842019543973942;
%          0.841042345276873;
%          0.840553745928339;
%          0.841530944625407;
%          0.844462540716612;
%          0.855211726384365;
%          0.881596091205212;
%          0.920684039087948;
%          0.944136807817590;
%          0.973452768729642;
%          0.996905537459284;
%          1.00863192182410;
%          1.01547231270358;
%          1.01889250814332;
%          1.01840390879479;
%          1.01840390879479;
%          1.01889250814332;
%          1.01938110749186];

% Degree of the interpolant polynomial
% n = 3;

% Coefficients of the interpolant polynomial
% For efficiency they're not calculated everytime!
% The coefficients for n = 3 are the defaul values provided.

%p = polyfit(x,y,n);

    p = [0.0135944481257312,...
         0.0565468050967924,...
        -0.637605175475846,...
         1.56991400907626];


    plotcheck = 0;
    if plotcheck

    x_lin = linspace(0.4,1.6);
    y_lin = linspace(0.7,1.3);


    corr = polyval(p,x_lin);

    plot(x,y,'*')
    hold on
    plot(x_lin,corr)
    ylim([0.7,1.3])
    xlim([0.4,1.6])
    end

    PSI = polyval(p,1/sigma);
end
