function [CHI_Ma] = traupel_Ma_transonic(M1,flag)

if flag == 0

  x_data = [0.00474978795589487; ...
      0.276873056262370; ...
      0.865648854961832;...
      1.39505230421261];

  y_data = [1.01295884733898; ...
      1.08012583072975; ...
      1.41969701849883; ...
      1.83161013946597];

  n = 3;

  p = polyfit(x_data,y_data,n);


  x = linspace(0,1.4);

  y_calc = polyval(p,x);

  plot(x,y_calc)
  hold on
  plot(x_data,y_data,'*')

else
 %for n = 3, p is already available;
   p = [-0.146134430423657, ...
        0.550877936302488, ...
        0.103083718685308, ...
        1.01245680911916];

   CHI_Ma = polyval(p,M1);
   
end

end
