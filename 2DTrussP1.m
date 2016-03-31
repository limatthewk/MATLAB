clc;  clear all; hold off;

%**************************************************************************
% This program will analyze a 2D truss structure. It reads an input .txt 
% file about the structure's nodes, elements, boundary conditions, and 
% degrees of freedom. Then it solves for the truss structure's nodal 
% displacements, reactions, and internal forces/stresses.
%**************************************************************************

%Read data from input file

finput = fopen('input_proj1.txt');
%finput = fopen('input.txt');

nnode = fscanf(finput, ['nodes: %d\n %*s %*s\n']);        % # of nodes
size_nc = [2 nnode];
node_coor = fscanf(finput,'%f',size_nc);          % coordinates of each node
node_coor = transpose(node_coor);

nelem = fscanf(finput, ['\n elements: %d\n %*s %*s %*s %*s\n']);% # of elements
size_ne = [4 nelem];
elemdata = fscanf(finput,'%f',size_ne); % node1, node2, Area, Young's Modulus
elemdata = transpose(elemdata);

nforce = fscanf(finput, ['\n force_BCs: %d\n %*s %*s %*s\n']);% # of force BCs
size_nf = [3 nforce];
forcedata = fscanf(finput,'%f',size_nf);    % node, degree of freeedom, value
forcedata = transpose(forcedata);

ndisp = fscanf(finput, ['\n displacement_BCs: %d\n %*s %*s %*s\n']);% number of displacement BCs
size_nf = [3 ndisp];
dispdata = fscanf(finput,'%f',size_nf); % node, degree of freedom, value
dispdata = transpose(dispdata);

% find E, A, L, and angle of each element in radians
for i = 1:1:nelem
    elem1 = elemdata(i,1)
    elem2 = elemdata(i,2)
    pt1x = node_coor(elem1,1)
    pt1y = node_coor(elem1,2)
    pt2x = node_coor(elem2,1)
    pt2y = node_coor(elem2,2)
    length(i) = sqrt((pt1x-pt2x)^2+((pt1y-pt2y)^2))
    angle(i) = atan((pt2y-pt1y)/(pt2x-pt1x))
    area(i) = elemdata(i,3)
    E(i) = elemdata(i,4)
    if pt2x < pt1x
        angle(i) = angle(i) + pi
    end
    if pt2y < pt1y && pt2x >= pt1x
        angle(i) = angle(i) + 2*pi
    end
end

%convert length in feet to inches
length = length*12

%store k matrices for each element
k = cell(nelem,1)
for i = 1:1:nelem
    k{i} = make_k(angle(i))*area(i)*E(i)/length(i)
end

%connectivity
for i = 1:1:nelem
    a(i,1) = [2*elemdata(i,1)-1]
    a(i,2) = [2*elemdata(i,1)]
    a(i,3) = [2*elemdata(i,2)-1]
    a(i,4) = [2*elemdata(i,2)]
end

%build complete element stiffness matrix from connectivity
ktotal = zeros(nnode*2)
for i = 1:1:nelem
    ktotal(a(i,1),a(i,1)) = ktotal(a(i,1),a(i,1))+ k{i}(1,1)
    ktotal(a(i,1),a(i,2)) = ktotal(a(i,1),a(i,2))+ k{i}(1,2)
    ktotal(a(i,1),a(i,3)) = ktotal(a(i,1),a(i,3))+ k{i}(1,3)
    ktotal(a(i,1),a(i,4)) = ktotal(a(i,1),a(i,4))+ k{i}(1,4)
    
    ktotal(a(i,2),a(i,1)) = ktotal(a(i,2),a(i,1))+ k{i}(2,1)
    ktotal(a(i,2),a(i,2)) = ktotal(a(i,2),a(i,2))+ k{i}(2,2)
    ktotal(a(i,2),a(i,3)) = ktotal(a(i,2),a(i,3))+ k{i}(2,3)
    ktotal(a(i,2),a(i,4)) = ktotal(a(i,2),a(i,4))+ k{i}(2,4)
    
    ktotal(a(i,3),a(i,1)) = ktotal(a(i,3),a(i,1))+ k{i}(3,1)
    ktotal(a(i,3),a(i,2)) = ktotal(a(i,3),a(i,2))+ k{i}(3,2)
    ktotal(a(i,3),a(i,3)) = ktotal(a(i,3),a(i,3))+ k{i}(3,3)
    ktotal(a(i,3),a(i,4)) = ktotal(a(i,3),a(i,4))+ k{i}(3,4)
    
    ktotal(a(i,4),a(i,1)) = ktotal(a(i,4),a(i,1))+ k{i}(4,1)
    ktotal(a(i,4),a(i,2)) = ktotal(a(i,4),a(i,2))+ k{i}(4,2)
    ktotal(a(i,4),a(i,3)) = ktotal(a(i,4),a(i,3))+ k{i}(4,3)
    ktotal(a(i,4),a(i,4)) = ktotal(a(i,4),a(i,4))+ k{i}(4,4)
end

%build altered penalty method element stiffness matrix
ktotal_pm = ktotal
kp = 10^30
for i = 1:1:ndisp
    n = dispdata(i,1)*2
    if dispdata(i,2) == 1
        n = n-1
    end
    ktotal_pm(n,n) = ktotal(n,n) + kp
end
for i = 1:1:2*nnode
    p(i,1) = 0
end

%build altered penalty method force matrix p
for i = 1:1:nforce
    if forcedata(i,3) ~= 0
        n = forcedata(i,1)*2
        if forcedata(i,2) == 1
            n = n-1
        end
        p(n) = forcedata(i,3)
    end
end

%ufinal are all displacements
ufinal = inv(ktotal_pm)*p;

%pfinal are all forces
pfinal = (ktotal)*ufinal;

%find internal forces and stresses for each element
for i = 1:1:nelem
    kf = [cos(angle(i))^2, cos(angle(i))*sin(angle(i)); cos(angle(i))*sin(angle(i)), sin(angle(i))^2];
    uf = [-ufinal(2*elemdata(i,1)-1)+ufinal(2*elemdata(i,2)-1);-ufinal(2*elemdata(i,1))+ufinal(2*elemdata(i,2))];
    fxy = E(i)*area(i)/length(i)*kf*uf;
    f(i) = sqrt(fxy(1)^2 + fxy(2)^2);
    x_1 = node_coor(elemdata(i,1),1);
    x_2 = node_coor(elemdata(i,2),1);
    y_1 = node_coor(elemdata(i,1),2);
    y_2 = node_coor(elemdata(i,2),2);
    newlength = sqrt((x_2 - x_1 + uf(1))^2 + (y_2 - y_1 + uf(2))^2)*12;
    if newlength < length(i)
        f(i) = -f(i);
    end
    stress(i) = f(i)/area(i);
end

%Up and Pp, displacements and forces easier to read
for i=1:nnode
   Up(i,1)=i;
   Up(i,2)=ufinal(2*i-1);
   Up(i,3)=ufinal(2*i);
   Pp(i,1)=i;
   Pp(i,2)=pfinal(2*i-1);
   Pp(i,3)=pfinal(2*i);
end

%Sp, internal forces and stresses, easier to read
for i=1:nelem
   Sp(i,1)=i;
   Sp(i,2)=f(i);
   Sp(i,3)=stress(i);
end

disp('     DISPLACEMENT RESULTS (inches)       ')
disp('     Node        x-dir(u)        y-dir(v)')
disp('                                         ')
fprintf('     %i %18.3E %15.3E\n',transpose(Up))
%***********************

disp('                                         ')
disp('                                         ')
disp('     REACTION RESULTS (lbs)         ')
disp('     Node       x-dir(u)        y-dir(v)')
disp('                                         ')
fprintf('     %i %18.3f %15.3f\n',transpose(Pp))
%**********************

disp('                                         ')
disp('                                         ')
disp('     MEMBER FORCES AND STRESSES          ')
disp('     Elem.       Force(lbs)   Stress(psi)')
disp('                                         ')
fprintf('     %i %18.3f %15.3f\n',transpose(Sp))
%**********************

mag = 1000;
limits=[min(node_coor(:,1))-1, max(node_coor(:,1))+1, min(node_coor(:,2))-1, max(node_coor(:,2))+1];
Coord1 = zeros(nelem,4);
for i=1:nelem
   Coord1(i,1:2) = node_coor(elemdata(i,1),1:2);
   Coord1(i,3:4) = node_coor(elemdata(i,2),1:2);
end
X1 = [Coord1(:,1) Coord1(:,3)];
Y1 = [Coord1(:,2) Coord1(:,4)];

node_coor = node_coor+mag*Up(:,2:3)/12;
for i=1:nelem
   Coord2(i,1:2) = node_coor(elemdata(i,1),1:2);
   Coord2(i,3:4) = node_coor(elemdata(i,2),1:2);
end
X2 = [Coord2(:,1) Coord2(:,3)];
Y2 = [Coord2(:,2) Coord2(:,4)];

plot(X1(1,:),Y1(1,:),'b','LineWidth',5),hold on
axis(limits)
plot(X2(1,:),Y2(1,:),'r--','LineWidth',2)
for i=2:nelem
    plot(X1(i,:),Y1(i,:),'b','LineWidth',5)
    plot(X2(i,:),Y2(i,:),'r--','LineWidth',2)
end
stitle = sprintf('Undeformed & Deformed Structure (magnification factor: %d)',mag);
title(stitle)
    