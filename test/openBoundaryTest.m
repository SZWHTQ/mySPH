clear; 
clc;
r_fluid = [3.14, -1];
r_fixed = [0, 1];
n_fixed = [1, 0];
Lb = 3;
r_new = r_fluid - ((r_fluid-r_fixed)*transpose(n_fixed)+Lb)*n_fixed;


pic1 = scatter(r_fluid(1), r_fluid(2), 100, "LineWidth",2);
hold on
pic2 = scatter(r_fixed(1), r_fixed(2), 100, "filled");
hold on
pic3 = scatter(r_new(1), r_new(2), 100, "LineWidth",2);
grid on

legend([pic1,pic2,pic3],"Fluid","Fixed","New",'Location','southeast')