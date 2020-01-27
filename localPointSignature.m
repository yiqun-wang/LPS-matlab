function [descriptor, Cf] = localPointSignature(shape_dir, shape_name, numEigs, radius_list, ref, S_diameter, if_output, output_dir)
%% 
addpath(genpath('utils'));

shape_fullname = [shape_dir, '/', shape_name];
if shape_name(end-2:end) == 'off'
    shape=read_shape(shape_fullname);
elseif shape_name(end-2:end) == 'obj'
    shape=read_shape_obj(shape_fullname);
elseif shape_name(end-2:end) == 'ply'
    shape=read_shape_ply(shape_fullname);
end

V = shape.VERT'; F = shape.TRIV';
shape.X = shape.VERT(:,1);
shape.Y = shape.VERT(:,2);
shape.Z = shape.VERT(:,3);
shape.n = size(shape.X,1);

if ref   
    diameter=S_diameter;
else
    diameter=shape_diameter(V',F');
end
radius_list = radius_list * diameter;
idxs  = [1:length(shape.X)]';
geods = compute_geods(shape,idxs);
nface = size(shape.TRIV, 1);
options.symmetrize = 1;
options.normalize = 0;
type = 'conformal';
stitching = false;
K = numEigs;
descriptor = zeros(shape.n, K*size(radius_list,2));
parfor i = 1:shape.n
    geo = geods(i,:);
    desc = [];
    for m = 1:size(radius_list,2)
        radius = radius_list(m);
        vertices = idxs(geo <= radius);
        faces = [];
        for j = 1:nface
            if ~ismember(0, ismember(shape.TRIV(j,:),vertices))
                faces = [faces; shape.TRIV(j,:)];
            end
        end
        vertices = unique(faces);
        re_faces = faces;
        for j = 1:size(faces,1)
            for k = 1:3
                re_faces(j,k) = find(vertices==faces(j,k));
            end
        end
        sub_shape_VERT = shape.VERT(vertices,:);
        sub_shape_TRIV = re_faces;
        sub_shape_n = size(sub_shape_VERT,1);
        [L,A] = compute_mesh_laplacian_plusA_half(sub_shape_VERT',sub_shape_TRIV',type,options);
        if stitching
            boundary_edge = compute_boundary_all(sub_shape_TRIV');
            W=diag(diag(L))-L;
            inner = 1:sub_shape_n;
            inner(boundary_edge) = [];
            boundary = sort(boundary_edge);
            bs = size(boundary, 2);
            is = size(inner, 2);
            AA = zeros(size(A) + is);
            AA(1:sub_shape_n, 1:sub_shape_n) = A;
            AA(boundary,boundary) = 2*AA(boundary,boundary);
            AA(sub_shape_n+1:end, sub_shape_n+1:end) = A(inner,inner);
            LL = zeros(size(W) + is);
            LL(1:sub_shape_n, 1:sub_shape_n) = W;
            LL(boundary,boundary) = 2*LL(boundary,boundary);
            LL(sub_shape_n+1:end,boundary)=W(inner,boundary);
            LL(boundary, sub_shape_n+1:end)=W(boundary, inner);
            LL(sub_shape_n+1:end, sub_shape_n+1:end)=W(inner,inner);
            LL = diag(sum(LL,2)) - LL;
            A=sparse(AA);
            L=sparse(LL);
            VERT = zeros(sub_shape_n+is, 3);
            VERT(1:sub_shape_n, :) = sub_shape_VERT;
            VERT(sub_shape_n+1:end, :) = sub_shape_VERT(inner, :);
            sub_shape_VERT = VERT;
            sub_shape_n = size(LL, 1);
        end
        try
            [V,D] = eigs(L, A, K+1, -1);
            V=V(:,2:end); D=D(2:end,2:end);
        catch
            % In case of trouble make the laplacian definite
            [V,D] = eigs(L , A+ 1e-20*speye(sub_shape_n), K+1, -1);
            V=V(:,2:end); D=D(2:end,2:end);
        end
        C = V' * A* sub_shape_VERT;
        Cf = D * sqrt(sum(C.^2,2));
        desc = [desc; Cf];
    end
    descriptor(i,:) = desc';
end
if if_output
    dlmwrite([output_dir, '/', shape_name(1:end-3), 'txt'], descriptor,'delimiter',' ', 'precision','%4.6e');
end

options.symmetrize = 1;
options.normalize = 0;
type = 'conformal';
method = 'single'; %single_hist %batch %single (others)
use_curv = false;
K = 300; % for the whole
nfirst = K; %for single
bin = 32; %for single_hist
dim = 64; batch = K/dim; % for batch

[L,A] = compute_mesh_laplacian_plusA_half(shape.VERT',shape.TRIV',type,options);

try
    [V,D] = eigs(L, A, K+1, -1);
catch
    % In case of trouble make the laplacian definite
    [V,D] = eigs(L , A+ 1e-20*speye(shape.n), K+1, -1);
end
V=V(:,2:end); D=D(2:end,2:end);


if use_curv
    curvature = [Cmax',Cmin'];
    flen = size(curvature, 2); 
    C = V' * A* curvature; 
else
    flen =size(shape.VERT, 2);
    C = V' * A* shape.VERT;
end


Cf = C;

Cf =  D *  sqrt(sum(Cf.^2,2));

end
