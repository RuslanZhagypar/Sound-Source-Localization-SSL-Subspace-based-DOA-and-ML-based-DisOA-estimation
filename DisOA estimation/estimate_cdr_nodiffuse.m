%ESTIMATE_CDR_NODIFFUSE
% Unbiased estimation of the Coherent-to-Diffuse Ratio (CDR) from the complex
% coherence of a mixed (noisy) signal, without assuming knowledge of the noise
% coherence. Equivalent to CDRprop4 in [1].
%
% CDR = estimate_cdr_nodiffuse(X, NaN, S)
%       X: complex coherence of mixed (noisy) signal
%       NaN: second argument is unused
%       S: coherence of signal component (magnitude one)
%
% Reference:
% Andreas Schwarz, Walter Kellermann, "Coherent-to-Diffuse Power Ratio
% Estimation for Dereverberation", IEEE/ACM Trans. on Audio, Speech and
% Lang. Proc., 2015 (under review); preprint available: arXiv:1502.03784
% PDF: http://arxiv.org/pdf/1502.03784
%
% Andreas Schwarz (schwarz@lnt.de)
% Multimedia Communications and Signal Processing
% Friedrich-Alexander-Universitaet Erlangen-Nuernberg (FAU)
% Cauerstr. 7, 91058 Erlangen, Germany
function CDR = estimate_cdr_nodiffuse(Cxx,~,Css)
Css = bsxfun(@times, ones(size(Cxx)), Css);

% limit the magnitude of Cxx to prevent numerical problems
magnitude_threshold = 1-1e-10;
critical = abs(Cxx)>magnitude_threshold;
Cxx(critical) = magnitude_threshold .* Cxx(critical) ./ abs(Cxx(critical));

CDR = imag(Cxx)./(imag(Css) - imag(Cxx));
CDR(imag(Css)./imag(Cxx)<=1) = Inf;
CDR(imag(Css)./imag(Cxx)<=0) = 0;

% Ensure we don't get any negative or complex results due to numerical effects
CDR = max(real(CDR),0);
end
