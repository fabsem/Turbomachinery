function F = root3ibeta(x,sigma)

%sigma = 1;
%x(1) = i0
%x(2) = n
%x(3) = beta1
F(1) = x(2) - 0.025 * sigma + 0.06 + ((x(3)/90)^(1 + 1.2 * sigma))/ ...
       (1.5 + 0.43 * sigma);
F(2) = x(1) - (x(3)^(0.914 + (sigma^3)/160))/(5 + 46 * exp(-2.3 * sigma)) + ...
        0.1 * sigma^3 * exp((x(3) - 70)/4);

end
