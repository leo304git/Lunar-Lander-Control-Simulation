function yke=ke(m,u) %kinetic energy
L=length(u);
yke=0;
for i=1:L
    yke=yke+0.5*m*(u(i))^2;
end
