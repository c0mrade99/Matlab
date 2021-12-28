%       MATLAB PROGRAM FOR TRUSS ELEMENT
%                   BY
%   MT21CDM012     PRATIK TAMBE

clear all
E = input("Enter Modulus of Elasticity (E in Pascals) : ");
A = input("Enter Area: ");
EA=E*A;
% element Nodes
elementNodes = input("Enter Connectivity Matrix: ");

nodeCoordinates = input("Enter Nodes Coordinate: ");

numberElements=size(elementNodes,1);

numberNodes=size(nodeCoordinates,1);
xx=nodeCoordinates(:,1);
yy=nodeCoordinates(:,2);

GDof=2*numberNodes;
U=zeros(GDof,1);
force=zeros(GDof,1);
a = int2str(width(force));
b = int2str(height(force));

fprintf('\n\nSize of force matrix must be in the form:')
disp(['Number of rows: ' a '     Number of Columns: ' b]) ;
disp('Example:');
disp(force') ;
force=input("Enter the array for loading conditions: ")
force = force';

% applied load at nodes
force

% computation of the stiffness matrix
[stiffness]=formStiffness2Dtruss(GDof,numberElements,elementNodes,numberNodes,nodeCoordinates,xx,yy,EA);
% boundary conditions and solution
prescribedDof=input("Enter Matrix for Boundary Conditions: ");
prescribedDof = prescribedDof'
stiffness
% solution
displacements=solution(GDof,prescribedDof,stiffness,force);
us=1:2:2*numberNodes-1;
vs=2:2:2*numberNodes;

% output displacements/reactions
outputDisplacementsReactions(displacements,stiffness,GDof,prescribedDof)
% stresses at elements
stresses2Dtruss(numberElements,elementNodes,xx,yy,displacements,E)

function [stiffness]=formStiffness2Dtruss(GDof,numberElements,elementNodes,numberNodes,nodeCoordinates,xx,yy,EA);
stiffness=zeros(GDof);
% computation of the system stiffness matrix
for e=1:numberElements;
% elementDof: element degrees of freedom (Dof)
indice=elementNodes(e,:) ;
elementDof=[indice(1)*2-1 indice(1)*2 indice(2)*2-1 indice(2)*2] ;
xa=xx(indice(2))-xx(indice(1));
ya=yy(indice(2))-yy(indice(1));
length_element=sqrt(xa*xa+ya*ya);
C=xa/length_element;
S=ya/length_element;
k1=EA/length_element*...
[C*C C*S -C*C -C*S; C*S S*S -C*S -S*S;
-C*C -C*S C*C C*S;-C*S -S*S C*S S*S];
stiffness(elementDof,elementDof)= stiffness(elementDof,elementDof)+k1;
end
end

function stresses2Dtruss(numberElements,elementNodes,xx,yy,displacements,E)
% stresses at elements
for e=1:numberElements
indice=elementNodes(e,:);
elementDof=[ indice(1)*2-1 indice(1)*2 indice(2)*2-1 indice(2)*2] ;
xa=xx(indice(2))-xx(indice(1));
ya=yy(indice(2))-yy(indice(1));
length_element=sqrt(xa*xa+ya*ya);
C=xa/length_element;
S=ya/length_element;
sigma(e)=E/length_element*[-C -S C S]*displacements(elementDof);
end
disp('stresses')
sigma'
end

function displacements=solution(GDof,prescribedDof,stiffness, force)
% function to find solution in terms of global displacements 
activeDof=setdiff([1:GDof]',[prescribedDof]);
U=stiffness(activeDof,activeDof)\force(activeDof); 
displacements=zeros(GDof,1); 
displacements(activeDof)=U;

end

function outputDisplacementsReactions(displacements,stiffness,GDof,prescribedDof)
disp('Displacements')
%displacements=displacements1;
jj=1:GDof; format
[jj' displacements]
% reactions
F=stiffness*displacements;
reactions=F(prescribedDof);
disp('reactions')
[prescribedDof reactions]
end