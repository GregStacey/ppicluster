function ARI = ari(predComplex, refComplex)

% Adjusted random index
%
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6124389/

Na = length(refComplex);
Nb = length(predComplex);

T = zeros(Na,Nb);

Ni = nan(size(refComplex));
if size(Ni,1)==1 && size(Ni,2)>1
  Ni = Ni';
end
for ii = 1:Na
  Ni(ii) = length(refComplex{ii});
  for jj = 1:Nb
    T(ii,jj) = length(intersect(refComplex{ii},predComplex{jj}));
  end
end

t1 = nan(Na,1);
for ii = 1:Na
    t1(ii) = nchoosek(length(refComplex{ii}), 2);
end
t1 = sum(t1);
t2 = nan(Nb,1);
for ii = 1:Nb
    t2(ii) = nchoosek(length(predComplex{ii}), 2);
end
t2 = sum(t2);

unqprotsa = horzcat(refComplex{:});
unqprotsb = horzcat(predComplex{:});
unqprots = unique([unqprotsa unqprotsb]);
nn = length(unqprots);
t3 = 2*t1*t2 / nn / (nn-1);

tmp = zeros(size(T));
for ii = 1:Na
    for jj = 1:Nb
        if T(ii,jj)<2; continue; end
        tmp(ii,jj) = nchoosek(T(ii,jj), 2);
    end
end

ARI = (nansum(tmp(:)) - t3) / (.5*(t1+t2) - t3);
