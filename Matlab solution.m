global A B;
A = Read('t1_icbm_normal_1mm_pn3_rf0.rawb',[181 217 181],'uint8');
B = Read('t2_icbm_normal_1mm_pn3_rf0.rawb',[181 217 181],'uint8');

PremRotCas = zeros(1, 7);

for i =-10: 0.1: 10
%% vektorji za transformacijo
t = [i 0 0];
r = [0 0 0];

transformVector = [t r];

%figure(1); imagesc(B(:,:,91)); colormap gray; axis equal; axis tight; axis xy;

%% priprava vhodne slike B
isize=size(B); % velikost slike v vokslih
vsize=[1 1 1]; % velikost voksla v mm

originx = floor(isize(1)/2);
originy = floor(isize(2)/2);
originz = floor(isize(3)/2);

origin=[originx originy originz];

Ind=1:prod(size(A));
[I1,I2,I3] = ind2sub(size(A),Ind');

xa=(I1-1)*vsize(1)-origin(1);
ya=(I2-1)*vsize(2)-origin(2);
za=(I3-1)*vsize(3)-origin(3);


%% priprava transformacije
t=[transformVector(1) transformVector(2) transformVector(3)]; % premik
r=[transformVector(4) transformVector(5) transformVector(6)]; % rotacija
%%


T=[ cos(r(2))*cos(r(3)), cos(r(3))*sin(r(1))*sin(r(2))-cos(r(1))*sin(r(3)), cos(r(1))*cos(r(3))*sin(r(2))+sin(r(1))*sin(r(3)) t(1);
    cos(r(2))*sin(r(3)), cos(r(1))*cos(r(3))+sin(r(1))*sin(r(2))*sin(r(3)), -cos(r(3))*sin(r(1))+cos(r(1))*sin(r(2))*sin(r(3)) t(2);
    -sin(r(2)), cos(r(2))*sin(r(1)), cos(r(1))*cos(r(2)), t(3);
    0 0 0 1];

%% interpolacija

Xa=[xa' ; ya' ; za' ;  ones(size(xa')) ];
Xb= T * Xa;
xb = (Xb(1,:) + origin(1))/(vsize(1));
yb = (Xb(2,:) + origin(2))/(vsize(2));
zb = (Xb(3,:) + origin(3))/(vsize(3));

ix = floor(xb)+1;
iy = floor(yb)+1;
iz = floor(zb)+1;
jx = ix+1;
jy = iy+1;
jz = iz+1;

wx = xb - ix +1;
wy = yb - iy +1;
wz = zb - iz +1;
w1 = (1-wx).*(1-wy).*(1-wz);
w2 = (wx).*(1-wy).*(1-wz);
w3 = (1-wx).*(wy).*(1-wz);
w4 = (wx).*(wy).*(1-wz);
w5 = (1-wx).*(1-wy).*(wz);
w6 = (wx).*(1-wy).*(wz);
w7 = (1-wx).*(wy).*(wz);
w8 = (wx).*(wy).*(wz);
w1=w1(:);
w2=w2(:);
w3=w3(:);
w4=w4(:);
w5=w5(:);
w6=w6(:);
w7=w7(:);
w8=w8(:);

t1 = [ix(:), iy(:), iz(:)];
t2 = [jx(:), iy(:), iz(:)];
t3 = [ix(:), jy(:), iz(:)];
t4 = [jx(:), jy(:), iz(:)];
t5 = [ix(:), iy(:), jz(:)];
t6 = [jx(:), iy(:), jz(:)];
t7 = [ix(:), jy(:), jz(:)];
t8 = [jx(:), jy(:), jz(:)];

%inteziteta1
TB=zeros(size(A));
H=zeros(256,256);
for v=1:prod(size(A)),
    if ((t1(v,1)<size(B,1) && t1(v,1)>0) && (t1(v,2)<size(B,2) && t1(v,2)>0) && (t1(v,3)<size(B,3) && t1(v,3)>0))      
        i1=B(t1(v,1),t1(v,2),t1(v,3))+1;
        i2=B(t2(v,1),t2(v,2),t2(v,3))+1;
        i3=B(t3(v,1),t3(v,2),t3(v,3))+1;
        i4=B(t4(v,1),t4(v,2),t4(v,3))+1;
        i5=B(t5(v,1),t5(v,2),t5(v,3))+1;
        i6=B(t6(v,1),t6(v,2),t6(v,3))+1;
        i7=B(t7(v,1),t7(v,2),t7(v,3))+1;
        i8=B(t8(v,1),t8(v,2),t8(v,3))+1;
        
       
        TB(v)=(i1*w1(v)) + (i2*w2(v)) + (i3*w3(v)) + (i4*w4(v)) + (i5*w5(v)) + (i6*w6(v)) + (i7*w7(v)) + (i8*w8(v))-1;
        H(i1,A(v)+1)= H(i1,A(v)+1)+w1(v);
        H(i2,A(v)+1)= H(i2,A(v)+1)+w2(v);
        H(i3,A(v)+1)= H(i3,A(v)+1)+w3(v);
        H(i4,A(v)+1)= H(i4,A(v)+1)+w4(v);
        H(i5,A(v)+1)= H(i5,A(v)+1)+w5(v);
        H(i6,A(v)+1)= H(i6,A(v)+1)+w6(v);
        H(i7,A(v)+1)= H(i7,A(v)+1)+w7(v);
        H(i8,A(v)+1)= H(i8,A(v)+1)+w8(v);

    end
end

figure(2); imagesc( log(H+1) ); axis equal; axis tight; colormap gray;

e0=Entropija(H);

eA0=Entropija( sum(H,1) );
eB0=Entropija( sum(H,2) );

MI=eA0+eB0-e0


%{
I=[A(:),TB(:)];
jointh = hist3(I,{0:255, 0:255});

%figure(5); imagesc( log(jointh+1) ); axis equal; axis tight; colormap gray;

sum(jointh(:))
e0=Entropija(jointh);

eA0=Entropija( sum(jointh,1) );
eB0=Entropija( sum(jointh,2) );

MI=eA0+eB0-e0
%}
end
% % return;