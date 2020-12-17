function [accuracy_table] = unmix_accuracy(basis_im, th_level, members, Name)
%unmix_accuracy Calculate accuracy of unmixing for a basis image
%   Basis image must be width by height by n_basis
%   members specifies correct unmixing endmembers
%   if members is not given, default will be identity matrix matching size
%   of basis_im

n_basis_members = size(basis_im, 3);

if ~exist('members','var')
    members = eye(n_basis_members);
end

if ~exist('Name','var')
    Name = cell(n_basis_members,1);
    for i = 1:n_basis_members
        Name{i} = sprintf('Basis %d', i);
    end
end

if ~exist('th_level','var')
    th_level = graythresh(basis_im);
end

% Take maximum along basis dimension, then gate threshold on that image
% Then get id's of positive pixels
mask = max(basis_im, [], 3) > th_level;
%th_basis_im = basis_im.*mask; % Masked image of basis members

th_basis_im_flat = basis_im(repmat(mask,1,1,n_basis_members));
th_basis_im_flat = reshape(th_basis_im_flat, [], n_basis_members);
im_flat_diff = zeros(size(th_basis_im_flat));
ground_truth_all = repmat(members, [1 1 size(th_basis_im_flat, 1)]);
ground_truth_all = permute(ground_truth_all, [3 1 2]);
% Pixel wise difference maps between thresholded image and members
for i = 1:size(members, 1)
    im_flat_diff(:, i) = sum(abs(squeeze(ground_truth_all(:,i,:)) - th_basis_im_flat), 2);
end

% For each pixel find which channel of im_flat_diff is the minumum
% This will NOT work for multi-colored cells, need to re-think
[~, minIDs] = min(im_flat_diff, [], 2); % ID of closest member
indS = sub2ind(size(th_basis_im_flat), (1:size(th_basis_im_flat,1))', minIDs);

Accuracy = zeros(n_basis_members, 1);
N_pixels = zeros(n_basis_members, 1);
Std_Dev = zeros(n_basis_members, 1);
unmix_sum = sum(th_basis_im_flat, 2);

% Stats by endmember
for i = 1:size(members, 1)
    unmix_correct = th_basis_im_flat(indS(minIDs==i));
    unmix_signal = unmix_sum(minIDs==i);
    Accuracy(i) = mean(unmix_correct./unmix_signal);
    N_pixels(i) = sum(minIDs==i);
    Std_Dev(i) = std(unmix_correct./unmix_signal);
end 

accuracy_table = table(Name, Accuracy, N_pixels, Std_Dev); 
end

