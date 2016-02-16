function [ matout ] = flatmat( matin )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

matout = reshape(matin,1,size(matin,1)*size(matin,2));

end

