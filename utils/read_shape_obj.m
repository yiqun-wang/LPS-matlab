function shape = read_shape_obj( shape_path )

[vertex, face] = read_obj(shape_path);
shape.VERT=vertex';
shape.TRIV=face';
shape.n = length(vertex);
shape.m = length(face);

end

