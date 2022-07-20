cubic_spline <- function(tp, knots)
{
  n=length(knots)
   f=rep(0,times=length(tp)*(n-1))
   dim(f)=c(length(tp),(n-1))
   #f1(x)=x#
   f[,1]=tp
   for(j in 1:(n-2))
   {
     for (i in 1:length(tp))
    {
      f[i,(j+1)]=max(0,(tp[i]-knots[j])^3, na.rm=TRUE)-max(0,(tp[i]-knots[n-1])^3, na.rm=TRUE)*(knots[n]-knots[j])/(knots[n]-knots[n-1])+
        max(0,(tp[i]-knots[n])^3, na.rm=TRUE)*(knots[n-1]-knots[j])/(knots[n]-knots[n-1])
          }
      }
    return(f)
  }



