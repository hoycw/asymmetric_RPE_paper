function aligned_elec = fn_align_elec2mesh(elec, mesh)
% Function to reposition an electrode to the closest point in a mesh.
% elec = fieldtrip structure with electrode info. elec.chanpos field required.
% mesh = cell array with fieldtrip mesh structures. mesh.pos field required.

% Unpack if several meshes
mesh_pos = [];
for mp = 1:numel(mesh)
    mesh_pos = [mesh_pos;mesh{mp}.pos];
end

%Find closest mesh position
new_pos = NaN(size(elec.chanpos));
for e = 1:size(elec.chanpos,1)
    
    % minimum euclidian distance
    [~,mix] = min(sqrt(sum((mesh_pos - elec.chanpos(e,:)).^2,2)));
    
    % find and store new position
    new_pos(e,:) = mesh_pos(mix(1),:); % if closest to several, selec first
    aligned_elec = elec;
    aligned_elec.chanpos = new_pos;
end

end