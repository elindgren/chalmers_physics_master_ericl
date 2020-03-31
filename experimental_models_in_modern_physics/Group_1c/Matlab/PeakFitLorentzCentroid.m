function [peak,C,ex,FWHM] = PeakFitLorentzCentroid(data,lambda,peak_guess,ran,span,plotta)
% PeakFitLorentzCentroid: Fits a Lorentzian and centroid and outputs peak position,
% centroid, extinction at peak, and FWHM (extinction from 0 to 100).
%       Usage:
%           [peak, C, ex, FWHM] = peakfitlorentz(data,lambda,peak_guess,ran,span)
%               data = matrix of spectra
%               lambda = wavelength vector
%               peak_guess = guess for peak position
%               ran = +- range around peak_guess to include in fit    
%               span = span for centroid   

war = warning('off','MATLAB:polyfit:RepeatedPointsOrRescale'); %turns off warnings
    
nr = size(data,2);
range1=peak_guess-ran;
range2=peak_guess+ran;
ranget=[range1 range2];
range=[find(ranget(1)-lambda<0,1,'first'):find(lambda-ranget(2)<0,1,'last')];

hc = 1.2398e3;
E = hc./lambda(range)';

s = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0,0,0],...
               'Upper',[3,1,5],...
               'Startpoint',[1 .5  hc/peak_guess]);
f = fittype('I*(HWHM^2/((x-p)^2+HWHM^2))','options',s);
%f = fittype('I*exp(-(x-p)^2/(2*HWHM^2))','options',s); %gaussian
hw = waitbar(0,'Please wait...');

ex = zeros(1,nr);
peak = zeros(1,nr);
FWHM = zeros(1,nr);
C = zeros(1,nr);


n = 20; % fit-polynom degree
lambdanew=lambda(range);
m=mean(lambdanew);
s=std(lambdanew);
lambdahat=(lambdanew-m)/s;
lambdahat = lambdahat';
left=peak_guess-span/2;

for i = 1:nr
tic
    waitbar(i / nr)
    
    a = fit(E,data(range,i),f);
    
    ex(i) = a.I;
    peak(i) = hc./a.p;
    FWHM(i) = hc./(a.p-a.HWHM)-hc./(a.p+a.HWHM);
    
    poly = polyfit(lambdahat,data(range,i),n);
    left=fzero(@(in) polyval(poly,(in-m)/s)-polyval(poly,(in+span-m)/s),left);
    ebase=polyval(poly,(left-m)/s);
    
    FI=polyval(poly,lambdahat);

    
    if (plotta == 1)&&(i==1)
        figure
        plot(lambda,data(:,i),'o',lambdanew,FI,'r')
    end
    if (plotta == 1)&&(i==1)
        figure
        plot(lambda(range),data(range,i),'o',lambda(range),a.I*(a.HWHM^2./((hc./lambda(range)-a.p).^2+a.HWHM^2)))
        drawnow
    end
    C1=0;
    C2=0;
    for b=1:n+1
        C1=C1+poly(b)/(n+3-b)*(((left+span-m)/s)^(n+3-b)-((left-m)/s)^(n+3-b));
        C2=C2+poly(b)/(n+2-b)*(((left+span-m)/s)^(n+2-b)-((left-m)/s)^(n+2-b));
    end
    C1=C1-ebase/2*(((left+span-m)/s)^2-((left-m)/s)^2);
    C2=C2-ebase*(((left+span-m)/s)-(left-m)/s);
    C(i)=C1/C2*s+m;
%   toc  
end

close(hw)

end

