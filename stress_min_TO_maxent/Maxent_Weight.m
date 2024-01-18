function [w  dwdx dwdy]=Maxent_Weight(xi, gpos, dm, v)
L=length(v);
dif=gpos*ones(1,L)-xi;

for i=1:L
    drdx=sign(dif(1,i))/dm(1,v(i));
    drdy=sign(dif(2,i))/dm(2,v(i));
    rx=abs(dif(1,i))/dm(1,v(i));
    ry=abs(dif(2,i))/dm(2,v(i));
    if rx>0.5
        wx=(4/3)-4*rx+4*rx^2-(4/3)*rx^3;
        dwx=(-4+8*rx-4*rx^2)*drdx;
    elseif rx<=0.5
        wx=(2/3)-4*rx^2+4*rx^3;
        dwx=(-8*rx+12*rx^2)*drdx;
    end

    if ry>0.5
        wy=(4/3)-4*ry+4*ry^2-(4/3)*ry^3;
        dwy=(-4+8*ry-4*ry^2)*drdy;
    elseif ry<=0.5
        wy=(2/3)-4*ry^2+4*ry^3;
        dwy=(-8*ry+12*ry^2)*drdy;
    end

    w(i)=wx*wy;
    dwdx(i)=wy*dwx;
    dwdy(i)=wx*dwy;
end
