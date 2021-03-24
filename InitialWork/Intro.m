numSamples=1000;
mu = 0;
alpha = 1;
x0=0;
HWHM = 1;

x = linspace(-5,5,numSamples);
xVoigt = linspace(-5,5,(numSamples*2)-1);
yGaus = GaussianPDF(x,mu,alpha);
yLorent = LorentzianPDF(x,x0,HWHM);
% yLorent = cauchypdf(x);
% yVoigt = conv(yGaus,yLorent);
% yVoigt = conv(yGaus,yLorent);
sigma = alpha/sqrt(2*log(2));
yVoigt = real(fadf((x + 1i*HWHM)/sigma/sqrt(2)))./ sigma./sqrt(2*pi);


figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,2,1)
plot(x,yGaus,x,yLorent);
grid on
title('Guassian and Lorentzian PDFs');
legend('Gaussian: \mu = 0, \sigma = 1','Lorentzian: x_0 = 0, \gamma = 1');
xlabel('x')
ylabel('P(x)')

subplot(1,2,2)
plot(x,yGaus);
grid on
title('Convolution of Guassian and Lorentzian PDFs')
hold on
plot(x,yLorent);
hold on
plot(x,yVoigt);
legend('Gaussian: \mu = 0, \sigma = 1','Lorentzian: x_0 = 0, \gamma = 1',...
    'Convolution');
xlabel('x')
ylabel('P(x)')