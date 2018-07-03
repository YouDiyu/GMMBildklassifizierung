clear all
close all
imgPath = '../../Sonnenlicht/2/';
fileNames = dir([imgPath '*.jpg']);
N = length(fileNames);
fileNo = [];
for i = 1:N
    temp = regexp(fileNames(i).name(1:end-4),'_','split');
    fileNo = [fileNo;str2double(temp(2))*1000+str2double(temp(3))];
end
[~,id] = sort(fileNo);
D = [];
imSize = [150,200];
for i = 1:50
    im = im2double(rgb2gray(imread([imgPath fileNames(id(i)).name])));
    if size(im,1) ~= imSize(1) & size(im,2) ~= imSize(2)
        im = imresize(im,imSize);
    end
    [row,colum,chan] = size(im);
    D(:,i) = im(:);
    if  rem(i,10) == 0
        fprintf('%d\n',i);
    end
end
pn = row * colum;           
K = 3;                     
mu = zeros(K,1);           
sigma2 = zeros(K,1);         
omg = zeros(K,1);        
mean_gm = zeros(pn,K);      
std_gm = zeros(pn,K);     
omg_gm = zeros(pn,K);     
for i=1:pn
    im=D(i,:)';
    t=0;
    for k=1:K
        mu(k)=im(5*k)+(k-2)*0.1000;
        sigma2(k)=0.1;
    end
    omg(:)=sigma2(:)./repmat(sum(sigma2(:)),K,1);
    E=Egama(im,omg,mu,sigma2,K);
    oldQ=0;
    newQ=Qfunc(E,im,omg,mu,sigma2,K);
    while(abs(sum(oldQ-newQ)) >= 0.01)
        for k=1:K
            mu(k)=sum(E(:,k).*im(:))/sum(E(:,k));
            sigma2(k)=sum(E(:,k).*((im(:)-repmat(mu(k),length(im),1)).^2))/sum(E(:,k));
            omg(k)=sum(E(:,k))/length(im);
        end
        E=Egama(im,omg,mu,sigma2,K);
        oldQ=newQ;
        newQ=Qfunc(E,im,omg,mu,sigma2,K);
        t=t+1;
    end
    mean_gm(i,:) = mu(:);              
    std_gm(i,:) = sqrt(sigma2(:));     
    omg_gm(i,:) = omg(:);           
end


thFeature=2000;
dir_bg = './output/method4/bg/';
dir_fg = './output/method4/fg/';
if exist(dir_bg,'dir') == 0
    mkdir(dir_bg);
end
if exist(dir_fg,'dir') == 0
    mkdir(dir_fg);
end
delete([dir_bg '/*'])
delete([dir_fg '/*'])
BeM=zeros(pn,1);
for i = 1:N
    im0 = im2double(imread([imgPath fileNames(id(i)).name]));
    im = rgb2gray(im0);
    if size(im,1) ~= imSize(1) & size(im,2) ~= imSize(2)
        im = imresize(im,imSize);
    end
    tmp=zeros(pn,K);
    Bew=zeros(pn,1);
    for k=1:K
        tmp(:,k)=im(:);
    end
    log3sigma=log(abs(tmp-mean_gm)-std_gm*3);
    for u=1:pn
        Bew(u)=isreal(sum(log3sigma(u,:)));
    end
    BeM(i)=sum(Bew(:));
    if BeM(i)>thFeature
        title([num2str(i) ': Besetzt  s:' num2str(BeM(i))]);
        imwrite(im0,[dir_fg  fileNames(id(i)).name]);
    else
        title([num2str(i) ': Frei  s:' num2str(BeM(i))]);
        imwrite(im0,[dir_bg  fileNames(id(i)).name]);
    end
    drawnow
%     pause(0.2)
end


function E=Egama(im,omg,mu,sigma2,K)
lang=length(im);
E=zeros(lang,K);
a=zeros(K,1);
for j=1:lang
    for k=1:K
        gauss=exp(-((im(j)-mu(k))^2)/(2*(sigma2(k))))/sqrt(2*pi*sigma2(k));
        a(k)=omg(k)*gauss;
    end
    E(j,:)=a(:)./sum(a(:))+0.0000001;
end
end

function Q=Qfunc(E,im,omg,mu,sigma2,K)
term1=zeros(K);
term2=zeros(K);

for k=1:K
    term1(k)=sum(E(:,k).*repmat(log(omg(k)),length(im),1));
    logwert=log(1/sqrt(2*pi))-log(sqrt(sigma2(k)));
    term2(k)=sum(E(:,k ).*repmat(logwert,length(im),1)-((im(:)-repmat(mu(k),length(im),1)).^2)./(2/sigma2(k)));
end 
Q=sum(term1(:))+sum(term2(:));
end 


