t = [1,2,10]';
w = ones(length(t),1)/3;
gamma = 1/2;
[grid, density]= compute_esd_ode(t, w, gamma);

%%
plot(grid,density)