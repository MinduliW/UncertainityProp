function [ status ] = downloadSW( filename, varargin )
%DOWNLOADSW Downloads Space Weather files from CelesTrack
%   
if length(varargin)==1
    fileURL = varargin{1};
else
    fileURL = 'http://www.celestrak.com/SpaceData/sw19571001.txt';
end

[~, status] = urlwrite(fileURL, filename);

end

