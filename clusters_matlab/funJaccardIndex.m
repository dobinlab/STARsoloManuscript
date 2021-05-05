function [JI, UmI]=funJaccardIndex(a,b)

u1=length(union(a,b));
i1=length(intersect(a,b));

JI=i1/u1;
UmI=u1-i1;