function eqp = gtcReadEq(equil)
% Reads GTC equilibrium output files.
%
% Syntax
%      function [eqp eqs] = gtcReadEq(equil)
%
% The only argument, 'equil', should be a file system path to the desired
% equilibrium.out file.  If no argument is supplied, the default value of 
% the argument is 'equilibrium.out' in Matlab's current directory.
%
% Two structs are returned from gtcReadEq.  The first, eqp, contains 1D
% profile data.  The second struct, eqs, contains 2D spline mesh data.
% Fields in the struct are named descriptively, according to the convention
% in the IDL script 'equilibrium.pro', as of GTC v3, June 29, 2011.
%
if nargin < 1
	equil = 'equilibrium.out';  % default equilibrium file name
end

% clear variables before we use them
clear n1D n2D A tmp %tmpsp tmppr

tmp = load(equil);

n1D = (1+tmp(1))*tmp(2);
n2D = (tmp(n1D+3)+2)*tmp(n1D+4)*tmp(n1D+5);

%A.pdata = reshape(tmp(3:n1D+2), [tmp(2) tmp(1)+1]);
%A.spdat = reshape(tmp(n1D+6:n1D+5+n2D), [tmp(n1D+4) tmp(n1D+5) tmp(n1D+3)+2]);

tmppr = reshape(tmp(3:n1D+2), [tmp(2) tmp(1)+1]);
tmpsp = reshape(tmp(n1D+6:n1D+5+n2D), [tmp(n1D+4) tmp(n1D+5) tmp(n1D+3)+2]);

eqp.psi = tmppr(:,1);
eqp.torpsi = tmppr(:,2);
eqp.r = tmppr(:,3);
eqp.R = tmppr(:,4);
eqp.Te = tmppr(:,5);
eqp.dlnTedpsi = tmppr(:,6);
eqp.ne = tmppr(:,7);
eqp.dlnnedpsi = tmppr(:,8);
eqp.Ti = tmppr(:,9);
eqp.dlnTidpsi = tmppr(:,10);
eqp.ni = tmppr(:,11);
eqp.dlnnidpsi = tmppr(:,12);
eqp.Tf = tmppr(:,13);
eqp.dlnTfdpsi = tmppr(:,14);
eqp.nf = tmppr(:,15);
eqp.dlnnfdpsi = tmppr(:,16);
eqp.zeff = tmppr(:,17);
eqp.torrot = tmppr(:,18);
eqp.Er = tmppr(:,19);
eqp.q = tmppr(:,20);
eqp.dlnqdpsi = tmppr(:,21);
eqp.g = tmppr(:,22);
eqp.p = tmppr(:,23);
eqp.rpsi = tmppr(:,24);
eqp.torpsi = tmppr(:,25);
eqp.rgpsi = tmppr(:,26);
eqp.psitor = tmppr(:,27);
eqp.psirg = tmppr(:,28);
eqp.errspcos = tmppr(:,29);
eqp.errspsin = tmppr(:,30);

eqp.x = tmpsp(:,:,1);
eqp.z = tmpsp(:,:,2);
eqp.b = tmpsp(:,:,3);
eqp.J = tmpsp(:,:,4);
eqp.i = tmpsp(:,:,5);
eqp.zeta2phi = tmpsp(:,:,6);
eqp.del = tmpsp(:,:,7);


%A.r = tmppd(:,3);
%A.R = tmppd(:,4);

