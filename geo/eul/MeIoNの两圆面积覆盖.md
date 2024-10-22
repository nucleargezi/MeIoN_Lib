#include<bits/z.h>
    point<ll>p1,p2;
    ll r,rr;
void sol(){
    ld dis = length(p1 - p2);
    ld al=2.0l*std::acos((r*r+dis*dis-rr*rr)/(2.0l*r*dis));
    ld R2=0.5l*al*r*r-0.5l*r*r*sin(al);
    ld beta=2.0l*std::acos((rr*rr+dis*dis-r*r)/(2.0l*rr*dis));
    ld R1=0.5l*beta*rr*rr-0.5l*rr*rr*sin(beta);
    ld R=R1+R2;
    std::cout<<R;
}