function V = simplexvolume(A, varargin)
%% SIMPLEXVOLUME - Volume of an n-simplex
%     SIMPLEXVOLUME(A) returns the volume of a n-simplex defined by n
%     verticies in A plus the origin (total of n+1 points). A is an
%     n-by-n, square matrix that specify n veritcies in columns.
%
%     V = SIMPLEXVOLUME(A, p) returns the same volume as above, with the
%     anchor point specified by p. p may be an n-by-m matrix with m
%     verticies in columns with n points each. V is an m-by-1 vector
%     of volumes
%
%     Examples: 
%          Area of a triangle (2-simplex) defined by the standard basis
%          A = eye(2);
%          V = simplexvolume(A)
%
%          Area of 5 3-simplexes defined by a random point in the standard
%          basis
%          A = eye(3);
%          p = rand(3, 5);
%          V = simplexvolume(A, p);
% 
%     Author: Eric Kercher, Northeastern University
%     Date:   18DEC2018

if ndims(A) ~= 2 %#ok<ISMAT>
    error('Input must be a 2-D matrix');
end

% System size, n = dimension of simplex, m = number of verticies 
[n, m] = size(A);

% A must be square
if m ~= n
    error('Input A must be a square matrix');
end

% Check for p input
if ~isempty(varargin)
    P = varargin{1};
else
    P = zeros(n, 1);
end

[n2, m2] = size(P);

% Check dimension of P
if n2 ~= n
    error('P must have the same dimension as A.');
end

V = zeros(m2, 1);

% Loop through anchor points
for i = 1:m2
    
    % "Anchor" the simplex verticies
    B = A - P(:, i);
    
    % Caluclate determinant
    V(i) = det(B);
end

% Conversion factor to compute volume
V = abs(V / factorial(n));
end