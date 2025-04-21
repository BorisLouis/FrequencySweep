function [ mov ] = getframes( path2file, frames )
%GETFRAMES Summary of this function goes here
%   Detailed explanation goes here

assert(min(size(frames))==1,'frames must be a vector of positive integers')

f_n = length(frames);

tObj = Tiff(path2file,'r');
l    = tObj.getTag(256);
w    = tObj.getTag(257);
tObj.setDirectory(frames(1));

im1  = tObj.read;
nClass = class(im1);
mov = zeros(w,l,f_n,nClass);
convert = false;

if length(size(im1))~=2
    warning('Colored images received, conversion to grayscale is performed')
    convert = true;
end

try 
    for i = 1:f_n
        f_i = frames(i);
        tObj.setDirectory(f_i)
        movTmp = tObj.read;  
        if convert
            movTmp = rgb2gray(movTmp);
        end
        mov(:,:,i) = movTmp;    
    end
catch e
    if strcmpi(e.identifier,'MATLAB:imagesci:Tiff:unableToChangeDir')
        warning('on')
        warning('Too many frames requested to be analyzed, skipping extra frames');
        warning('off')
        disp(['Found a total of ' num2str(i-1) ' frames'])
        mov= mov(:,:,1:i-1);
    end

end
tObj.close


end

