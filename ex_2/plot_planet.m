function [] = plot_planet(x0,y0,z0,radius,planet)

% CHANGE THE FILE TO THE APPROPIATE PLANET JPG FILE
planetTexture = imread(planet); 
planetTexture = flipud(planetTexture);

[x, y, z] = sphere;
x = radius*x + x0;
y = radius*y + y0;
z = radius*z + z0;
surf(x,y,z,'FaceColor','texturemap','CData',planetTexture,'EdgeColor','none');

end