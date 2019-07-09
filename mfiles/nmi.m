function NMI = nmi(predComplex, refComplex)

% Normalized mutual information
%
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6124389/

Na = length(refComplex);
Nb = length(predComplex);

unqprotsa = horzcat(refComplex{:});
unqprotsb = horzcat(predComplex{:});
unqprots = unique([unqprotsa unqprotsb]);
nn = length(unqprots);
Pi = nan(1, Na);
for ii = 1:Na
    Pi(ii) = length(refComplex{ii}) / nn;
end
Pj = nan(1, Nb);
for ii = 1:Nb
    Pj(ii) = length(predComplex{ii}) / nn;
end


% 1. calculate I(C, C')
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
II = nan(size(T));
for ii = 1:Na
    for jj = 1:Nb
        II(ii,jj) = T(ii,jj) / nn * log2(T(ii,jj) / nn / Pi(ii) / Pj(jj));
    end
end
II = nansum(II(:));


% 2. calculate H(C) = I(C,C)
T = zeros(Na,Na);
Ni = nan(size(refComplex));
if size(Ni,1)==1 && size(Ni,2)>1
  Ni = Ni';
end
for ii = 1:Na
  Ni(ii) = length(refComplex{ii});
  for jj = 1:Na
    T(ii,jj) = length(intersect(refComplex{ii},refComplex{jj}));
  end
end
H1 = nan(size(T));
for ii = 1:Na
    for jj = 1:Na
        H1(ii,jj) = T(ii,jj) / nn * log2(T(ii,jj) / nn / Pi(ii) / Pi(jj));
    end
end
H1 = nansum(H1(:));


% 2. calculate H(C') = I(C',C')
T = zeros(Nb,Nb);
Ni = nan(size(predComplex));
if size(Ni,1)==1 && size(Ni,2)>1
  Ni = Ni';
end
for ii = 1:Nb
  Ni(ii) = length(predComplex{ii});
  for jj = 1:Nb
    T(ii,jj) = length(intersect(predComplex{ii},predComplex{jj}));
  end
end
H2 = nan(size(T));
for ii = 1:Nb
    for jj = 1:Nb
        H2(ii,jj) = T(ii,jj) / nn * log2(T(ii,jj) / nn / Pj(ii) / Pj(jj));
    end
end
H2 = nansum(H2(:));



NMI = II / sqrt(H1*H2);

