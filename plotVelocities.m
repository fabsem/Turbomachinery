function plotVelocities(v1ax,v2ax,v4ax,v1,v2,v4,w1,w2,w3,w4,v1t,v2t,v4t,w1t,w2t,w3t,w4t,U1,U2,station)

if station==[1 0 0]
    location='Hub';
elseif station==[0 1 0]
    location='Midspan';
elseif station==[0 0 1]
    location='Tip';
end
    
V1=[v1t v1ax];
W1=[w1t v1ax];
V2=[v2t v2ax];
W2=[w2t v2ax];
W3=[w3t v2ax];
V4=[v4t v4ax];
W4=[w4t v4ax];

figure
subplot(2,1,1)
vectarrow([0 0],V1)
text(v1t/2,v1ax/2,'V1')
hold on
vectarrow([0 0],W1)
text(w1t/2,v1ax/2,'W1')
hold on
vectarrow(W1,V1)
text((v1t+w1t)/2,v1ax,'U1')
hold on
vectarrow([0 0],V2)
text(v2t/2,v2ax/2,'V2')
hold on
vectarrow([0 0],W2)
text(w2t/2,v2ax/2,'W2')
hold on
vectarrow(W2,V2)
text((v2t+w2t)/2,v2ax,'U1')
%legend('v1','w1','U1','v2','w2','U1')
title(['Rotor 1 ',location])
set(gca, 'YDir','reverse')
axis([-200 300 0 200])

subplot(2,1,2)
vectarrow([0 0],V2)
text(v2t/2,v2ax/2,'V2')
hold on
vectarrow([0 0],W3)
text(w3t/2,v2ax/2,'W3')
hold on
vectarrow(W3,V2)
text((v2t+w3t)/2,v2ax,'U2')
hold on
vectarrow([0 0],V4)
text(v4t/2,v4ax/2,'V4')
hold on
vectarrow([0 0],W4)
text(w4t/2,v4ax/2,'W4')
hold on
vectarrow(W4,V4)
text((v4t+w4t)/2,v4ax,'U2')
%legend('v2','w3','U2','v4','w4','U2')
set(gca, 'YDir','reverse')
title(['Rotor 2 ', location])
axis([-200 300 0 200])

end
