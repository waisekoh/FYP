# Beta Binomial Test #
# ?bb.test
#to install ibb locally, install.packages('ibbfile', lib='/mnt/gtklab01/waise/lib')
#to call it library('ibb', lib.loc='/mnt/gtklab01/waise/lib')
library('ibb')
library('MASS')
args<-commandArgs(TRUE)
group <- c(rep("CTR", as.numeric(args[1])), rep("cKO", as.numeric(args[2])))
nc<-as.numeric(args[1])+as.numeric(args[2])
data<-read.table(args[3])
data<- as.matrix(data)
jc<-data[,1] #jc is now only containing entire column 1
x<-data[,2:(nc+1)] # x now contains columns 2 to 7
y<-data[,(nc+2):((2*nc)+1)] #y now contains columns 8 to 13
Class<-data[,ncol(data)] #contains the last column, which is the JN Gain, Loss, Diff
if(class(y[1,1]) != "numeric"){
	x.n <- matrix(as.numeric(x),ncol=ncol(x))
	y.n <- matrix(as.numeric(y),ncol=ncol(y))
	x <- x.n
	y <- y.n
}
out<-bb.test(x, y, group, n.threads = nc) 
head(y)
p<-as.numeric(out$p.value)
for(j in 1:length(p))
{
   	if(is.nan(p[j]))
   	{
   		p[j]=min(p)
   	}
}

# Multiple testing correction using BH #
bp<-p.adjust(p, "bonferroni",length(p))
bhp<-p.adjust(p, "BH",length(p))
byp<-p.adjust(p, "BY",length(p))

# Writing to file #
data<-cbind(jc,x,y,p,bp,bhp,byp,Class)
write.matrix(data, file = args[4], sep = "\t")