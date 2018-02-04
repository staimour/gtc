function A = gtcReadProfile(profile)
% function gtcReadProfile reads in a gtc input profile
% and stores in struct A.  
% Header information is parsed and used for field names.

% provide a default file name, if one is not given by user
if nargin < 1
	profile = 'profile.dat'
end

% read the file
fid = fopen(profile);
	head = textscan(fid, '%s', 15);
	data = textscan(fid, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
fclose(fid);

% parse headers for invalid characters, and replace
for i=1:length(head)
	head{i}=strrep(head{i},'-','');
	head{i}=strrep(head{i},'+','_');
end

% store our data into a convenient struct
clear A;
A = cell2struct(data,head{1},2);
