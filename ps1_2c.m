% Uncertainy Quantification: PS1

close all; clear all; clc

% Write a program simulating the process of tossing the die N times for
% each of the possible values of n where n= the number of faces with one
% dot on them out of six.

% Start by setting N
N=1000;
n=1;
switch n

    case 0 %n=0, no faces contain 1 dot
    r0 = randi([2 6],1,N);
    m0=find(r0==1);
    if isempty(m0)
        m0=0;
    end
    evidence = m0/N;
    prior= (5/6).^6;
    likelihood = 0;

    case 1 %: n=1, 1 face contains single dot
    r1= randi([1 6], 1,N);
    m1=find(r1==1);
    if isempty(m1)
        m1=0;
    end
    m1=length(m1);
    evidence = m1/N;
    prior= 6*(1/6)*(5/6).^5;
    likelihood = 1/6;

    case 2 %: n=2, 2 facea contains single dot
    r2= randi([1 6], 1,N);
    % replace all 2s with a 1
    r2(r2==2)=1;
    m2=find(r2==1)
    if isempty(m2)
        m2=0;
    end
    m2=length(m2);

    case 3 %: n=3, 3 faces contains single dot
    r3= randi([1 6], 1,N);
    % replace all 2s and 3s with a 1
    r3(r3==2)=1;
    r4(r3==3)=1;
    m3=find(r3==1);
    if isempty(m3)
        m3=0;
    end
    m3 =length(m3);


    case 4 %: n=4, 4 faces contains single dot
    r4= randi([1 6], 1,N);
    % replace all 2s, 3s and 4s with a 1
    r4(r4==2)=1;
    r4(r4==3)=1;
    r4(r4==4)=1;
    m4=find(r4==1);
    if isempty(m4)
        m4=0;
    end
    m4 =length(m4);

    case 5 %: n=5, 5 faces contains single dot
    r5= randi([1 6], 1,N);
    % replace all 2s, 3s, 4s, and 5s with a 1
    r5(r5==2)=1;
    r5(r5==3)=1;
    r5(r5==4)=1;
    m5=find(r5==1);
    if isempty(m5)
        m5=0;
    end
    m5 = length(m5);

    case 6 %: n=6, 6 faces contains single dot
    r6= randi([1 6], 1,N);
    % replace all 2s, 3s, 4s, 5s, and 6s with a 1
    r6(r6==2)=1;
    r6(r6==3)=1;
    r6(r6==4)=1;
    m6=find(r6==1);
    if isempty(m6)
        m6=0;
    end
    m6 = length(m6);
end

