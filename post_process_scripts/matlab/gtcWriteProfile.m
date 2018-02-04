function gtcWriteProfile(profileData,outputFileName)
% function gtcWriteProfile(profileData,outputFileName)
%
% function gtcWriteProfile takes user inputted data, profileData, and uses it
% to write a GTC input profile.  The argument profileData may be an array or
% a struct.  Arrays must be 2D, and have 15 columns.  Structs must have 15
% fields which contain numerical vectors of equal length.

% provide default output filename.
	if nargin < 2
		outputFileName = 'profile.dat';
	end

disp(['Output file: ' outputFileName]);
	
% check whether our data is a matrix or a struct.
	if isstruct(profileData) % convert to a matrix
		tmp = profileData;
		clear profileData
		profileData = struct2array(tmp);
		clear tmp
	end
		
% make sure our matrix has the proper number of columns
	[nrow ncol] = size(profileData);
	if ncol ~= 15
		error('Input data must be a either a struct with 15 fields containing numerical vectors of equal length or an array with 15 columns.');
	end

% open the file to write to.
	findex = fopen(outputFileName,'w');	

% write the header
	fprintf(findex,'   Pol-Flux            x              r              ');
		fprintf(findex,'R               R+r             Te              ');
		fprintf(findex,'ne              Ti            Zeff          ');
		fprintf(findex,'omega-tor           Er              ni           ');
		fprintf(findex,'nimp             nf              Tf\n');

% write the data
	for i=1:nrow
		for j=1:ncol
			fprintf(findex,'%16.7e',profileData(i,j));
		end
		fprintf(findex,'\n');
	end
	
% write the footer
	fprintf(findex,'\n');
	fprintf(findex,'Poloidal Flux in Webers\n');
	fprintf(findex,'x is square root of normalized toroidal flux\n');
	fprintf(findex,'r is midplane radius in cm (midplane halfwidth of flux ');
		fprintf(findex,'surface width at midplane)\n');
	fprintf(findex,'R is flux surface center in cm\n');
	fprintf(findex,'R+r is local major radius in cm on the outer midplane\n');
	fprintf(findex,'Te is electron temperature in eV\n');
	fprintf(findex,'ne is electron density in m^-3\n');
	fprintf(findex,'Ti is ion temperature in eV; last two points are artificial\n');
	fprintf(findex,'Zeff is from Carbon density profile measurement\n');
	fprintf(findex,'omega-tor is measured angular velocity from carbon ');
		fprintf(findex,'rotation in radians/sec;to get votr, multiply omega ');
		fprintf(findex,'by local R from equilibrium (or from R+r)\n');
	fprintf(findex,'Er is radial electric field in volts/cm\n');
	fprintf(findex,'ni is ion density in m^-3\n');
	fprintf(findex,'nimp is impurity density in 10^19 m^-3\n');
	fprintf(findex,'nf is fast ion density in m^-3\n');
	fprintf(findex,'Tf is fast ion temperature in eV\n');

% close profile.dat file
	fclose(findex);

% done
