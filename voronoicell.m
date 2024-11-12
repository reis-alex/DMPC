function [voro] = voronoicell(actual_robot,other_robots)

[~,~,others] = size(other_robots(1,1,:));
interp_x{1} = actual_robot(1,:,1);
interp_y{1} = actual_robot(2,:,1);

for j = 1:others
    interp_x{j+1} = other_robots(1,:,j);
    interp_y{j+1} = other_robots(2,:,j);
    sym{j} = reflect([actual_robot(1,:,1); actual_robot(2,:,1)],[other_robots(1,:,j); other_robots(2,:,j)]);
end

interp = vertcat(horzcat([interp_x{:}]',[interp_y{:}]'),[sym{:}]');
dt      = delaunayTriangulation(interp);
[V,R_v] = voronoiDiagram(dt);
jj = 1;

% generate box
Vbox = [];
Vbox = [Vbox; [actual_robot(1,:,1)  actual_robot(2,:,1)] + 1*[+0 +0.5]];
Vbox = [Vbox; [actual_robot(1,:,1)  actual_robot(2,:,1)] + 1*[-0 -0.5]];
Vbox = [Vbox; [actual_robot(1,:,1)  actual_robot(2,:,1)] + 1*[+0.5 -0]];
Vbox = [Vbox; [actual_robot(1,:,1)  actual_robot(2,:,1)] + 1*[-0.5 -0]];
box = Polyhedron('V',Vbox);

vert = V([R_v{1}(:)],:);
voro = Polyhedron('V',vert(~isinf(vert(:,1)),:));
% voro{i}(jj) = voro{i}(jj).intersect(Xc);
voro = voro.intersect(box);
voro.minHRep();

clear interp_x interp_y interp V R_v dt sym vert j jj Vbox box

end
