clear all

E = input('Enter the modulus of elasticity(E) in Pascal :');
L = input('Enter the length of a rod in mm:');
numberElements=input('Enter the number of Elements :');
A1=input('Enter Area 1 in mm^2 :');
A2 = input('Enter Area 2 in mm^2 :');
fixedDof = input('Input matrix with fixed node points: ');

fixedDof = fixedDof'
temp = (A1-A2)/(numberElements-1);
A = zeros(1,(numberElements-1));
t=0;
for k=1:numberElements;
    A(1,k) = A1-A(1,k);
    t = temp+t;
    A(1,k+1) = t;    

end
A(end)=[];
A
% generation of coordinates and connectivities
% numberElements: number of elements
% generation equal spaced coordinates in xx
nodeCoordinates=linspace(0,L,numberElements+1);
xx=nodeCoordinates;
xx
% numberNodes: number of nodes
numberNodes=size(nodeCoordinates,2);
% elementNodes: connections at elements
ii=1:numberElements;
elementNodes(:,1)=ii;
elementNodes(:,2)=ii+1;
% for structure:
% displacements: displacement vector
% force : force vector
% stiffness: stiffness matrix
displacements=zeros(numberNodes,1);
force=zeros(numberNodes,1);
stiffness=zeros(numberNodes,numberNodes);
% applied load at node 2
force(2)=250000;
% computation of the system stiffness matrix
for e=1:numberElements;
% elementDof: element degrees of freedom (Dof)
elementDof=elementNodes(e,:) ;
nn=length(elementDof);
length_element=nodeCoordinates(elementDof(2))...
-nodeCoordinates(elementDof(1));
detJacobian=length_element/2;invJacobian=1/detJacobian;
% central Gauss point (xi=0, weight W=2)
[shape,naturalDerivatives]=shapeFunctionL2(0.0);
Xderivatives=naturalDerivatives*invJacobian;
% B matrix
B=zeros(1,nn); B(1:nn) = Xderivatives(:);
stiffness(elementDof,elementDof)=stiffness(elementDof,elementDof)+B'*B*2*detJacobian*E*A(1,e);

end
stiffness
% boundary conditions and solution
% prescribed dofs
%fixedDof=find(xx==min(nodeCoordinates(:)) | xx==max(nodeCoordinates(:)))';

prescribedDof=[fixedDof];
% free Dof : activeDof
activeDof=setdiff([1:numberNodes]',[prescribedDof]);
% solution
GDof=numberNodes;
displacements=solution(GDof,prescribedDof,stiffness,force);
% output displacements/reactions
outputDisplacementsReactions(displacements,stiffness,numberNodes,prescribedDof)

sigma = zeros(1,numberElements);
for i=1:numberElements;
    aa = displacements(i+1,1)- displacements(i,1);
    sigma(1,i) = aa*E/(xx(1,i+1)-xx(1,i)  );
end
sigma
%-------------------------------------------------------------
%PLOT
%-------------------------------------------------------------
figure(1)
plot(xx',displacements);
xlabel('Length')
ylabel('Displacement')
title('Node VS Displacement')
grid
numberofElements = (1:numberElements);
figure(2)
plot(numberofElements',sigma');
xlabel('Element')
ylabel('Stress')
title('Element VS Stress')
grid
%................................................................
%%%%%%%%Functions%%%%%%%%%%%%
%................................................................
function displacements=solution(GDof,prescribedDof,stiffness,force)
% function to find solution in terms of global displacements
activeDof=setdiff([1:GDof]',[prescribedDof]);
U=stiffness(activeDof,activeDof)\force(activeDof);
displacements=zeros(GDof,1);
displacements(activeDof)=U;
end

function outputDisplacementsReactions(displacements,stiffness,GDof,prescribedDof)
% output of displacements and reactions in
% tabular form
% GDof: total number of degrees of freedom of
% the problem
% displacements
disp('Displacements')
%displacements=displacements1;
jj=1:GDof; format
[jj' displacements]

% reactions
F=stiffness*displacements;
F
reactions=F(prescribedDof);
disp('reactions')
[prescribedDof reactions]
end

function [shape,naturalDerivatives]=shapeFunctionL2(xi)
% shape function and derivatives for L2 elements
% shape : Shape functions
% naturalDerivatives: derivatives w.r.t. xi
% xi: natural coordinates (-1 ... +1)
shape=([1-xi,1+xi]/2)';
naturalDerivatives=[-1;1]/2;
end % end function shapeFunctionL2
