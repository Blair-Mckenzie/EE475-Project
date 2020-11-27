numSamples=200;
mu = 0;
standardDev = 1;
x0=0;
HWHM = 1;

x = linspace(-6,6,numSamples);
xVoigt = linspace(-6,6,(numSamples*2)-1);
yGaus = GaussianPDF(x,mu,standardDev);
yLorent = LorentzianPDF(x,x0,HWHM);
yVoigt = conv(yGaus,yLorent);

figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,2,1)
plot(x,yGaus);
grid on
title('Guassian and Lorentzian PDFs');
hold on
lorent = plot(x,yLorent);
legend('Gaussian','Lorentzian');

subplot(1,2,2)
plot(x,yGaus);
grid on
title('Convolution of Guassian and Lorentzian PDFs')
hold on
plot(x,yLorent);
plot(xVoigt,yVoigt);
legend('Guassian','Lorentzian','Convolution');