function kout = fidall_vap_rotation(kin,phi,theta)
%FIDALL_VAP_ROTATION 3D Rotation equivalent to vap_phi/theta in fidall
%  kout = fidall_vap_rotation(kin,phi,theta)
%   kin   k-space (or grad) input [#pts/interleave,#interleaves,#groups]
%   phi   Rotation(s) around z                [deg]
% theta   Rotation(s) around y, then around z [deg]                     (0)
%  kout   Rotated k-space (or grad)
%
% 12/2022  Rolf Schulte
if (nargin<1), help(mfilename); return; end

if ~exist('theta','var'), theta = []; end
if isempty(theta),        theta = 0; end


%% pre-compute rotations
cosphi = cosd(phi(:)).';
sinphi = sind(phi(:)).';
costheta = cosd(theta(:)).';
sintheta = sind(theta(:)).';


%% actual rotations
% [#pts/interleave,#interleaves,#groups]
kout(:,:,1) =  ...
    + bsxfun(@times,kin(:,:,1),(costheta.*cosphi)) ...
    + bsxfun(@times,kin(:,:,2),(sinphi)) ...
    - bsxfun(@times,kin(:,:,3),(sintheta.*cosphi));
kout(:,:,2) =  ...
    - bsxfun(@times,kin(:,:,1),(costheta.*sinphi)) ...
    + bsxfun(@times,kin(:,:,2),(cosphi)) ...
    + bsxfun(@times,kin(:,:,3),(sintheta.*sinphi));
kout(:,:,3) = ...
    + bsxfun(@times,kin(:,:,1),(sintheta)) ...
    + bsxfun(@times,kin(:,:,3),(costheta));

end      % fidall_vap_rotation.m


% 
% if isempty(theta)
%     % 2D rotations around z: right-handed, i.e. counter-clockwise
%     rot(1,:) = [cosphi , -sinphi , 0];
%     rot(2,:) = [sinphi ,  cosphi , 0];
%     rot(3,:) = [    0  ,       0 , 1];
% 
% else
%     % 3D rotations
%     costheta = cosd(theta);
%     sintheta = sind(theta);
% 
%     % rotation of theta around y then of phi around z
%     rot(1,:) = [ costheta*cosphi , -costheta*sinphi , sintheta];
%     rot(2,:) = [          sinphi ,           cosphi ,        0];
%     rot(3,:) = [-sintheta*cosphi ,  sintheta*sinphi , costheta];
% end
