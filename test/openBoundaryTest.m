clear; 
clc;
r_fluid = [1, 1];
r_fixed = [0.1, 0.1];
n_fixed = [-1, 0];
Lb = 3;
r_new = r_fluid - ((r_fluid-r_fixed)*transpose(n_fixed)+Lb)*n_fixed;


scatter(r_fluid(1), r_fluid(2), 100, "LineWidth",2)
hold on
scatter(r_fixed(1), r_fixed(2), 100, "filled")
hold on
scatter(r_new(1), r_new(2), 100, "LineWidth",2)
grid on