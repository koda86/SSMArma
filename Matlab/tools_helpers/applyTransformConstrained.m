function out = applyTransformConstrained(vertices,T)
    if size(vertices,2) ~= 3
        vertices = vertices';
    end

    % Horizontal scaling only
    vertices(:,[1,3]) = vertices(:,[1,3]) * T.ScaleHorizontal;
    % Vertical (y-axis) unchanged

    % Apply rotation
    vertices = (T.Rotation * vertices')';

    % Apply translation
    vertices = vertices + T.Translation';

    out = vertices;
end
