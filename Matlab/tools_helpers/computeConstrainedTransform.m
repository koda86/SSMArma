function T = computeConstrainedTransform(floating, target, scale, sole_indices)

    assert(size(floating,2) == 3);
    assert(size(target,2) == 3);

    % Remove vertical offset (sole at y=0)
    floating(:,2) = floating(:,2) - mean(floating(sole_indices,2));
    target(:,2) = target(:,2) - mean(target(sole_indices,2));

    % Compute horizontal centroids
    centroid_floating = mean(floating(:,[1,3]),1);
    centroid_target = mean(target(:,[1,3]),1);
    
    % Horizontal translation
    translation_horizontal = centroid_target - centroid_floating;
    translation = [translation_horizontal(1), 0, translation_horizontal(2)];
    floating_translated = floating + translation;

    % Compute separate scaling (horizontal only)
    if scale
        % Horizontal scale factor (based on x,z-plane only)
        dist_floating = sqrt(sum(floating_translated(:,[1,3]).^2,'all'));
        dist_target = sqrt(sum(target(:,[1,3]).^2,'all'));
        scale_factor_horizontal = dist_target / dist_floating;
    else
        scale_factor_horizontal = 1;
    end

    % Apply horizontal scaling
    floating_scaled = floating_translated;
    floating_scaled(:,[1,3]) = floating_scaled(:,[1,3]) * scale_factor_horizontal;
    % Vertical scaling is constrained (factor = 1), so Y is untouched

    % Compute rotation around vertical (y-axis) only
    p = floating_scaled(:,[1,3]) - mean(floating_scaled(:,[1,3]));
    q = target(:,[1,3]) - mean(target(:,[1,3]));
    M = p'*q;
    [U,~,V] = svd(M);
    R2D = V * U';
    if det(R2D) < 0
        V(:,end) = -V(:,end);
        R2D = V * U';
    end
    Rotation = eye(3);
    Rotation([1,3],[1,3]) = R2D;

    % Final transformation
    T.ScaleHorizontal = scale_factor_horizontal;
    T.ScaleVertical = 1; % constrained
    T.Rotation = Rotation;
    T.Translation = translation';
end
