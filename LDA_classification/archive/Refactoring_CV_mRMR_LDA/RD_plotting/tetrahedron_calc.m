
sz = 10;

x = 0; 
y = 0; 
z = 0; 
pt = [x, y, z]; 

tetra_X = [x + sz*sqrt(8/9), x - sz*sqrt(2/9), x - sz*sqrt(2/9), x]; 
tetra_Y = [               y, y + sz*sqrt(2/3), y - sz*sqrt(2/3), y];
tetra_Z = [    z - sz*(1/3),     z - sz*(1/3),     z - sz*(1/3), z + sz];

plot3(tetra_X, tetra_Y, tetra_Z)