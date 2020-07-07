function mesh = meshinv(mesh)
% function mesh = meshinv({vertices,faces})
% function faces = meshinv(faces)

% Thomas Deneux
% Copyright 2005-2017

if iscell(mesh), faces = mesh{2}; else faces = mesh; end
if size(faces,1)~=3, faces=faces'; if size(faces,1)~=3, error('faces component should have 3 rows'), end, end
if iscell(mesh), mesh{2} = faces([1 3 2],:); else mesh = faces([1 3 2],:); end


