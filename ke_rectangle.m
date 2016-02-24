function [ KE ] = ke_rectangle(a,b,h,E,nu)
% DESCRIPTION:
%   this function computes the element stiffness matrix of an arbitrary 2D
%   rectangle element (a special case of a quadrilaterial element). the 
%	element has four nodes with two d.o.f. on each node. the material is 
%   assumed isotropic.
%
%   integrations are solved analytically using a symbolic manipulation
%   software package (MuPAD).
%
% INPUT:
%   ======================
%   |                    |
%   |                    |
%   |                    2b
%   |                    |
%   |                    |
%   ========= 2a =========
%   h - thickness
%   E - Young's modulus
%   nu - Poisson's ratio
%
% OUTPUT:
%   KE: element stiffness matrix (8*8)
%
% REFERENCE:
%   https://engineering.purdue.edu/~ahvarma/CE595/
%
% LOG:
%   160222: created by Yuqing Zhou @ Univ. of Michigan;
%   please contact yuqingz@umich.edu if you have any questions.

KE_p = zeros(8,8);

KE_p(1,1) = (8/3 - (8*nu)/3)*a^2 + (16*b^2)/3;
KE_p(1,2) = 2*a*b*(nu + 1);
KE_p(1,3) = (4/3 - (4*nu)/3)*a^2 - (16*b^2)/3;
KE_p(1,4) = 2*a*b*(3*nu - 1);
KE_p(1,5) = ((4*nu)/3 - 4/3)*a^2 - (8*b^2)/3;
KE_p(1,6) = -2*a*b*(nu + 1);
KE_p(1,7) = ((8*nu)/3 - 8/3)*a^2 + (8*b^2)/3;
KE_p(1,8) = -2*a*b*(3*nu - 1);

KE_p(2,1) = 2*a*b*(nu + 1);
KE_p(2,2) = (16*a^2)/3 - (8*b^2*(nu - 1))/3;
KE_p(2,3) = -2*a*b*(3*nu - 1);
KE_p(2,4) = (8*a^2)/3 + (8*b^2*(nu - 1))/3;
KE_p(2,5) = -2*a*b*(nu + 1);
KE_p(2,6) = ((4*nu)/3 - 4/3)*b^2 - (8*a^2)/3;
KE_p(2,7) = 2*a*b*(3*nu - 1);
KE_p(2,8) = (4/3 - (4*nu)/3)*b^2 - (16*a^2)/3;

KE_p(3,1) = (4/3 - (4*nu)/3)*a^2 - (16*b^2)/3;
KE_p(3,2) = -2*a*b*(3*nu - 1);
KE_p(3,3) = (8/3 - (8*nu)/3)*a^2 + (16*b^2)/3;
KE_p(3,4) = -2*a*b*(nu + 1);
KE_p(3,5) = ((8*nu)/3 - 8/3)*a^2 + (8*b^2)/3;
KE_p(3,6) = 2*a*b*(3*nu - 1);
KE_p(3,7) = ((4*nu)/3 - 4/3)*a^2 - (8*b^2)/3;
KE_p(3,8) = 2*a*b*(nu + 1);

KE_p(4,1) = 2*a*b*(3*nu - 1);
KE_p(4,2) = (8*a^2)/3 + (8*b^2*(nu - 1))/3;
KE_p(4,3) = -2*a*b*(nu + 1);
KE_p(4,4) = (16*a^2)/3 - (8*b^2*(nu - 1))/3;
KE_p(4,5) = -2*a*b*(3*nu - 1);
KE_p(4,6) = (4/3 - (4*nu)/3)*b^2 - (16*a^2)/3;
KE_p(4,7) = 2*a*b*(nu + 1);
KE_p(4,8) = ((4*nu)/3 - 4/3)*b^2 - (8*a^2)/3;

KE_p(5,1) = ((4*nu)/3 - 4/3)*a^2 - (8*b^2)/3;
KE_p(5,2) = -2*a*b*(nu + 1);
KE_p(5,3) = ((8*nu)/3 - 8/3)*a^2 + (8*b^2)/3;
KE_p(5,4) = -2*a*b*(3*nu - 1);
KE_p(5,5) = (8/3 - (8*nu)/3)*a^2 + (16*b^2)/3;
KE_p(5,6) = 2*a*b*(nu + 1);
KE_p(5,7) = (4/3 - (4*nu)/3)*a^2 - (16*b^2)/3;
KE_p(5,8) = 2*a*b*(3*nu - 1);

KE_p(6,1) = -2*a*b*(nu + 1);
KE_p(6,2) = ((4*nu)/3 - 4/3)*b^2 - (8*a^2)/3;
KE_p(6,3) = 2*a*b*(3*nu - 1);
KE_p(6,4) = (4/3 - (4*nu)/3)*b^2 - (16*a^2)/3;
KE_p(6,5) = 2*a*b*(nu + 1);
KE_p(6,6) = (16*a^2)/3 - (8*b^2*(nu - 1))/3;
KE_p(6,7) = -2*a*b*(3*nu - 1);
KE_p(6,8) = (8*a^2)/3 + (8*b^2*(nu - 1))/3;

KE_p(7,1) = ((8*nu)/3 - 8/3)*a^2 + (8*b^2)/3;
KE_p(7,2) = 2*a*b*(3*nu - 1);
KE_p(7,3) = ((4*nu)/3 - 4/3)*a^2 - (8*b^2)/3;
KE_p(7,4) = 2*a*b*(nu + 1);
KE_p(7,5) = (4/3 - (4*nu)/3)*a^2 - (16*b^2)/3;
KE_p(7,6) = -2*a*b*(3*nu - 1);
KE_p(7,7) = (8/3 - (8*nu)/3)*a^2 + (16*b^2)/3;
KE_p(7,8) = -2*a*b*(nu + 1);

KE_p(8,1) = -2*a*b*(3*nu - 1);
KE_p(8,2) = (4/3 - (4*nu)/3)*b^2 - (16*a^2)/3;
KE_p(8,3) = 2*a*b*(nu + 1);
KE_p(8,4) = ((4*nu)/3 - 4/3)*b^2 - (8*a^2)/3;
KE_p(8,5) = 2*a*b*(3*nu - 1);
KE_p(8,6) = (8*a^2)/3 + (8*b^2*(nu - 1))/3;
KE_p(8,7) = -2*a*b*(nu + 1);
KE_p(8,8) = (16*a^2)/3 - (8*b^2*(nu - 1))/3;

KE = (E/(1-nu^2))*(h/(16*a*b))*KE_p;

end
