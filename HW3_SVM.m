%HW #3 Classifying digits using SVM instead of Eigendigits
% 1. Load the data
close all; clear all; clc
%loads in the digits.mat file from the data folder
load ('../data/digits.mat');
whos;

%% Preview an example image from the training set
i = 500;
M=trainImages(:,:,1,i);
figure;
imagesc(M)
title(['Known digit=', num2str(trainLabels(i))])


%% Make the x by k matrix A where x is the total number of pixels in an
% image and k is the number of training images


k = 400;
for i = 1:k
    % make the column vector A_k
    M = trainImages(:,:,1,i); % M is square matrix
    x= size(M,1)*size(M,2);
    Acol = reshape(M, [x,1]);
    A(:,i) = Acol;
end
