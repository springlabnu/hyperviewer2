function C = pcat3(A, P)
%pcat Patterned Concatenation. 
%   pcat(A, P) concatentes each page of 3D array A according to the pattern
%   specified in P. 

% Check A is at least 3D
if ndims(A) ~= 3
    error('A must be at least 3 dimensions');
end

% Check P is at most 2D
if ndims(P) > 2
    error('P must be a vector or matrix.');
end

x = size(A, 1);
y = size(A, 2);
N = size(A, 3);

% Restrict to square matricies for now
if x ~= y
    error('A must contain square matricies.');
end

% Number of pages in A must equal the numel in cat pattern
if numel(P) ~= N
    error('Dimsion mismatch.')
end

% Dimensions of P
[a, b] = size(P);

% Size of result
C = zeros(x*a, y*b);

% Loop through the pattern array
for i = 1:a
    for j = 1:b
        
        % Pattern index 
        idx = P(i, j);
     
        % Indicies of C to add page of A
        i1 = i*x - 1; iend = i1 + (x - 1);
        j1 = j*y - 1; jend = j1 + (y - 1);
        
        % Add to result
        C(i1:iend, j1:jend) = A(:, :, idx);
        
    end
end
end

% To upgrade for arbirary dimension of A
% function B = getDimA(A, dim)
% 
% % Pick out the n'th dimensional page of A specified by dim
% 
% 
% inds = repmat({1},1,ndims(A));
% inds{1} = 1:size(A, 1);
% inds{2} = 1:size(A, 2);
% inds{dim} = dim;
% B = A(inds{:});
% 
% 
% end
