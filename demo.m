clc; clf; clear;
addpath(genpath('utils/'));

%%
mesh_dir =  './data/';
s_name = 'tr_reg_000.off';

if_output = true;
output_dir = './desc/';

% the number of eigenfunctions for LPS
K = 16; 
radius = 0.05;
% the dimension of the output desc = length(coff) * K
coff = [3,2,1]; % coff =  [4.5, 4, 3.5, 3, 2.5, 2, 1.5, 1]; 
radius_list = radius*coff;

ref = false;  % generating desc for partial matching
if ref   
    ref_name = 'cuts_victoria_shape_0.off'; % for example, you can direct set S_diameter manually
    ref_fullname = [mesh_dir, '/', ref_name];
    if ref_name(end-2:end) == 'off'
        shape_ref=read_shape(ref_fullname);
    elseif ref_name(end-2:end) == 'obj'
        shape_ref=read_shape_obj(ref_fullname);
    elseif ref_name(end-2:end) == 'ply'
        shape_ref=read_shape_ply(ref_fullname);
    end
    S_diameter=shape_diameter(shape_ref.VERT,shape_ref.TRIV);
else
    S_diameter = 0;
end

[LPS, LPS_global] = localPointSignature(mesh_dir, s_name, K, radius_list, ref, S_diameter, if_output, output_dir);