%To find node in neighboring of gauss point
function [v,L]=nodes_in_support(numnode, x, gpos, dm)
one1=ones(1,numnode);
gv=gpos*one1;
dif=gv-x;
c=all(abs(dif)<=dm);
v=find(c);
L=length(v);