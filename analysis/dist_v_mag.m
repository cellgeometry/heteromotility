%% Calculate distance and vector magnitude at each point in a divergence matrix

function [result] = dist_v_mag(div, m_u, m_v, m_s);
% div = N x M divergence matrix
% m_u = N x M matrix of vector field x-component
% m_v = N x M matrix of vector field y-component
% m_s = 2 x 1 vector i, j index for minimum divergence (metastable state)
result = zeros(length(div(:)), 2);

for k = 1:length(div(:));
    [i,j] = ind2sub(size(div), k);
    u = [i j];
    d = sqrt( sum((m_s-u).^2) );
    flux = [m_u(i,j), m_v(i,j)];
    mag = sqrt( dot(flux,flux) );
    result(k,:) = [d mag];
end
end