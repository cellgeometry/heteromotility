function [result] = point_distmag(p, m_u, m_v, metastable);
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%
% p : point in question, [i,j] indices 
% m_u : matrix of vector i components 
% m_v : matrix of vector j components
% metastable : N x 1 matrix of metastable state linear indices
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Returns
%%%%%%%%%%%%%%%%%%%%%%%%%%
% result : 2 x 1 matrix containing minimum distance to metastable state
%          and magnitude of the vector at that point

distances = zeros(length(metastable), 1);

for k = 1:length(metastable);
    [i,j] = ind2sub(size(m_u), k);
    m = [i, j];
    d = sqrt( sum((p-m).^2) );
    distances(k,1) = d;
end

flux = [m_u(p(1),p(2)), m_v(p(1),p(2))];
mag = sqrt( dot(flux,flux) );

result = [min(distances), mag];

end