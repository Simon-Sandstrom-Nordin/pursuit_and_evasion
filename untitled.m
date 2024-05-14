colormap([1 1 1; 0 0 0]);
M = random_walk_map_generator(.25, 302, 25);
figure(1)
imagesc(M)
axis equal

colormap([1 1 1; 0 0 0]);
M = random_walk_map_generator(.5, 302, 25);
figure(2)
imagesc(M)
axis equal

colormap([1 1 1; 0 0 0]);
M = random_walk_map_generator(.75, 302, 25);
figure(3)
imagesc(M)
axis equal
