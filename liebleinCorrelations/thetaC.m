X = [8.50020974543508e-05;
     0.0439383569204954;
     0.118662928046276;
     0.225680568741307;
     0.313217274192480;
     0.393521073896629;
     0.460889099860906;
     0.533148610160510;
     0.588817256529706;
     0.622988099706356;
     0.651571986841234;
     0.673687987106175;
     0.694513500982492];
Y = [0.00587817101980438;
     0.00575276532797562;
     0.00592674364692115;
     0.00654405758064160;
     0.00787023749070865;
     0.00962356213984501;
     0.0113715879569323;
     0.0144119400348840;
     0.0180189727625314;
     0.0210437227238941;
     0.0292274744441746;
     0.0374085767484306;
     0.0508938099338382];
 
 n = 10;
 p = polyfit(X,Y,n);
 
 x = linspace(0,0.7);
 thetac = polyval(p,x);
 
 figure(1)
 plot(X,Y,'*')
 hold on
 plot(x,thetac)
 
 %%
 %constrained
x0 = [X(1);X(end)];
y0 = [Y(1);Y(end)];


%x = eta85(1:22,1);
%y = eta85(1:22,2);
x = X;
y = Y;
n = 15;
v(:,n+1) = ones(length(x),1,class(x));

for j = n:-1:1
    v(:,j) = x.*v(:,j+1);
end

c = v;

d = y;

A = [];
b = [];

Aeq = x0.^(n:-1:0);
beq = y0;

p = lsqlin(c,d,A,b,Aeq,beq);

yhat = polyval(p,x);

plot(x,y,'.b-')

hold on

%plot(x0,y0,'gx','linewidth',4)
plot(x,yhat,'r','linewidth',2)
%hold off