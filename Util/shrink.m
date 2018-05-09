function shrink=shrink(x,r)

shrink=sign(x).*max((abs(x)-r),0);