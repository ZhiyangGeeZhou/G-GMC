setwd('E:/R/Working Directory') #设置工作目录路径
library(MASS)
S=c(6,2,rep(3,7))
D=as.matrix(read.csv("(18,6^12^13^7).csv",head=FALSE))
contrast=list(
	level_2=array(c(1/sqrt(2),-1/sqrt(2)),dim=c(1,2)),
	level_3=array(c(1/sqrt(2),1/sqrt(6),0,-2/sqrt(6),-1/sqrt(2),1/sqrt(6)),dim=c(2,3)),
	level_4=array(0,dim=c(3,4)),
	level_5=array(0,dim=c(4,5)),
	level_6=array(c(1/sqrt(2),1/sqrt(6),1/sqrt(12),1/sqrt(20),1/sqrt(30),0,0,0,0,-5/sqrt(30),0,0,0,-4/sqrt(20),1/sqrt(30),0,0,-3/sqrt(12),1/sqrt(20),1/sqrt(30),0,-2/sqrt(6),1/sqrt(12),1/sqrt(20),1/sqrt(30),-1/sqrt(2),1/sqrt(6),1/sqrt(12),1/sqrt(20),1/sqrt(30)),dim=c(5,6))
)

L=function(x,subdesign){
	s=S[subdesign];d=D[,subdesign];size=length(subdesign)
	temp_Xd=NULL
	temp_Xd[size]=1
	for (i in (size-1):1){
		temp_Xd[i]=temp_Xd[i+1]*s[i+1]
	}
	temp_Xd=array(temp_Xd,dim=c(size,1))
	Xd=d%*%temp_Xd+1
	temp_s=s[1]
	Px=switch(x[1]+1,array(1/sqrt(temp_s),dim=c(1,temp_s)),contrast[[temp_s-1]])
	for (i in 2:length(s)){
		temp_s=s[i]
		Px=kronecker(Px,switch(x[i]+1,array(1/sqrt(temp_s),dim=c(1,temp_s)),contrast[[temp_s-1]]))
	}
	temp=Px[,Xd]
	if (is.matrix(temp)==FALSE) Lx=as.matrix(temp)
	else Lx=t(temp)
	return(Lx)
}

LL=function(j,subdesign){ #返回Lj以及每个Ly的最后一列对应的列号
	size=length(subdesign)
	y=c(rep(1,j),rep(0,size-j))
	if (j==0){
		Ly=L(y,subdesign)
		col_No.=1
	}else{
		location_y=1:j
		end=(size-j+1):size
		Ly=NULL
		col_No.=NULL
		repeat{	
			Ly=cbind(Ly,L(y,subdesign))
			col_No.=c(col_No.,dim(Ly)[2])
			if (sum(location_y==end)==j) break 
			else{
				if (location_y[j]!=size){
					location_y[j]=location_y[j]+1
					y=rep(0,size)	
					y[location_y]=1
				}
				else{
					temp_diff=which(location_y[-1]-location_y[-j]>1)
					location_y_move=temp_diff[length(temp_diff)]
					location_y[location_y_move]=location_y[location_y_move]+1
					location_y[(location_y_move+1):j]=location_y[location_y_move]+1:(j-location_y_move)
					y=rep(0,size)
					y[location_y]=1
				}
			}
		}
	}
	return(list(Ly,c(0,col_No.)))
}

GWLP=function(subdesign){ #返回GWLP
	s=S[subdesign];size=length(subdesign)
	A=NULL
	for (j in 1:size){
		Lj=LL(j,subdesign)[[1]]
		A=c(A,prod(s)/dim(D)[1]^2*sum(apply(Lj,2,sum)^2))
	}
	return(round(c(1,A),6)) #用四舍五入去除计算误差的影响
}

MMA_moments=function(subdesign){
	s=S[subdesign];size=length(subdesign);d=D[,subdesign]
	N=dim(d)[1]
	delta=array(0,dim=c(N,N))
	Kp=NULL
	for (p in 1:size){
		for (i in 1:(N-1)){
			for (j in (i+1):N){
				deltaijk=NULL
				for (k in 1:size){
					deltaijk=c(deltaijk,1-(d[i,k]!=d[j,k]))
				}
				delta[i,j]=sum(s*deltaijk)
			}
		}
		Kp=c(Kp,sum(delta^p))
	}
	return(round(Kp*2/N/(N-1),6)) #用四舍五入去除计算误差的影响
}

GAENP_C=function(i,j,subdesign){
	s=S[subdesign];size=length(subdesign)
	Li=LL(i,subdesign)[[1]]
	Lj=LL(j,subdesign)[[1]]
	Li=t(t(Li)/sqrt(apply(Li^2,2,sum))) #列单位化
	Lj=t(t(Lj)/sqrt(apply(Lj^2,2,sum))) #列单位化
	Cij=crossprod(Lj,Li)^2
	Cij=round(Cij,6) #用四舍五入去除计算误差的影响
	temp_order=1:dim(Cij)[2]
	for (k in temp_order){
		Cij[,k]=Cij[,k][order(Cij[,k])]
	}
	if (length(temp_order)>1){
		for (k in 2:dim(Cij)[2]){
			for (l in (k-1):1){
				temp_diff=Cij[,temp_order[l]]-Cij[,k]
				if (sum(abs(temp_diff))>0 & temp_diff[which(temp_diff!=0)[1]]>0){
					temp_order[l+1]=temp_order[l]
					temp_order[l]=k
				}else break
			}	
		}
	}
	return(as.vector(Cij[,temp_order]))
}

G2AENP_C=function(i,j,subdesign){
	s=S[subdesign];size=length(subdesign)
	Cij=NULL
	Lj=LL(j,subdesign)
	x=c(rep(1,i),rep(0,size-i))
	if (i==0){
		Lx=L(x,subdesign)
		Rxj=NULL
		inverse=ginv(crossprod(Lx))
		for (k in 1:(length(Lj[[2]])-1)) Rxj=c(Rxj,sum(((inverse%*%crossprod(Lx,Lj[[1]]))^2)[,(Lj[[2]][k]+1):Lj[[2]][k+1]]))
		Cij=cbind(Cij,Rxj[order(Rxj)])
	}else{
		location_x=1:i
		end=(size-i+1):size
		repeat{	
			Lx=L(x,subdesign)
			Rxj=NULL
			inverse=ginv(crossprod(Lx))
			for (k in 1:(length(Lj[[2]])-1)) Rxj=c(Rxj,sum(((inverse%*%crossprod(Lx,Lj[[1]]))^2)[,(Lj[[2]][k]+1):Lj[[2]][k+1]]))
			Cij=cbind(Cij,Rxj[order(Rxj)])
			if (sum(location_x==end)==i) break
			else{
				if (location_x[i]!=size){
					location_x[i]=location_x[i]+1
					x=rep(0,size)	
					x[location_x]=1
				}
				else{
					temp_diff=which(location_x[-1]-location_x[-i]>1)
					location_x_move=temp_diff[length(temp_diff)]
					location_x[location_x_move]=location_x[location_x_move]+1
					location_x[(location_x_move+1):i]=location_x[location_x_move]+1:(i-location_x_move)
					x=rep(0,size)
					x[location_x]=1
				}
			}
		}
	}
	Cij=round(Cij,6) #用四舍五入去除计算误差的影响
	temp_order=1:dim(Cij)[2]
	if (length(temp_order)>1){
		for (k in 2:(dim(Cij)[2])){
			for (l in (k-1):1){
				temp_diff=Cij[,temp_order[l]]-Cij[,k]
				if (sum(abs(temp_diff))>0 & temp_diff[which(temp_diff!=0)[1]]>0){
					temp_order[l+1]=temp_order[l]
					temp_order[l]=k
				}else break
			}	
		}
	}
	return(as.vector(Cij[,temp_order]))
}

GAENP=function(order,subdesign){
	AENP=NULL
	for (i in 1:dim(order)[2]){
		AENP=c(AENP,GAENP_C(order[1,i],order[2,i],subdesign))
	}
	return(AENP)
}

G2AENP=function(order,subdesign){
	AENP=NULL
	for (i in 1:dim(order)[2]){
		AENP=c(AENP,G2AENP_C(order[1,i],order[2,i],subdesign))
	}
	return(AENP)
}

classify_GWLP=function(column_stay,column_initial){ #column_stay为必须包含的列，column_initial为遍历的起始
	size=length(column_initial)
	column_traversal=column_initial
	subdesign=c(column_stay,column_traversal)
	end=(length(S)-size+1):length(S)
	result=list(subdesign)
	GWLP_compare=t(as.matrix(GWLP(subdesign)))
	while(sum(column_traversal==end)!=size){
		if (column_traversal[size]!=length(S)) {
			column_traversal[size]=column_traversal[size]+1
		}else{
			temp_diff=which(column_traversal[-1]-column_traversal[-size]>1)
			location_move=temp_diff[length(temp_diff)]
			column_traversal[location_move]=column_traversal[location_move]+1
			column_traversal[(location_move+1):size]=column_traversal[location_move]+1:(size-location_move)
		}
		subdesign=c(column_stay,column_traversal)
		GWLP=GWLP(subdesign)
		k=1
		while (k<=length(result)){
			if (sum(GWLP_compare[k,]==GWLP)==length(GWLP)){
				result[[k]]=rbind(result[[k]],subdesign);break
			}
			k=k+1
		}
		if (k==length(result)+1){
			result=c(result,list(subdesign))
			GWLP_compare=rbind(GWLP_compare,GWLP)
		}
	}
	temp_order=1:dim(GWLP_compare)[1]
	if (length(temp_order)>1){
		for (k in 2:dim(GWLP_compare)[1]){
			for (l in (k-1):1){
				temp_diff=GWLP_compare[temp_order[l],]-GWLP_compare[k,]
				if (sum(abs(temp_diff))>0 & temp_diff[which(temp_diff!=0)[1]]>0){
					temp_order[l+1]=temp_order[l]
					temp_order[l]=k
				}else break
			}	
		}
	}
	return(c(result,list(GWLP_compare),list(temp_order)))
}

classify_MMA=function(column_stay,column_initial){ #column_stay为必须包含的列，column_initial为遍历的起始
	size=length(column_initial)
	column_traversal=column_initial
	subdesign=c(column_stay,column_traversal)
	end=(length(S)-size+1):length(S)
	result=list(subdesign)
	MMA_compare=t(as.matrix(MMA_moments(subdesign)))
	while(sum(column_traversal==end)!=size){
		if (column_traversal[size]!=length(S)) {
			column_traversal[size]=column_traversal[size]+1
		}else{
			temp_diff=which(column_traversal[-1]-column_traversal[-size]>1)
			location_move=temp_diff[length(temp_diff)]
			column_traversal[location_move]=column_traversal[location_move]+1
			column_traversal[(location_move+1):size]=column_traversal[location_move]+1:(size-location_move)
		}
		subdesign=c(column_stay,column_traversal)
		MMA=MMA_moments(subdesign)
		k=1
		while (k<=length(result)){
			if (sum(MMA_compare[k,]==MMA)==length(MMA)){
				result[[k]]=rbind(result[[k]],subdesign);break
			}
			k=k+1
		}
		if (k==length(result)+1){
			result=c(result,list(subdesign))
			MMA_compare=rbind(MMA_compare,MMA)
		}
	}
	temp_order=1:dim(MMA_compare)[1]
	if (length(temp_order)>1){
		for (k in 2:dim(MMA_compare)[1]){
			for (l in (k-1):1){
				temp_diff=MMA_compare[temp_order[l],]-MMA_compare[k,]
				if (sum(abs(temp_diff))>0 & temp_diff[which(temp_diff!=0)[1]]>0){
					temp_order[l+1]=temp_order[l]
					temp_order[l]=k
				}else break
			}	
		}
	}
	return(c(result,list(MMA_compare),list(temp_order)))
}

classify_GAENP=function(order,subdesign_initial){
	size=length(subdesign_initial)
	subdesign=subdesign_initial
	end=(length(S)-size+1):length(S)
	result=list(subdesign)
	AENP_compare=t(as.matrix(GAENP(order,subdesign)))
	while(sum(subdesign==end)!=size){
		if (subdesign[size]!=length(S)){
			subdesign[size]=subdesign[size]+1
		}else{
			temp_diff=which(subdesign[-1]-subdesign[-size]>1)
			location_move=temp_diff[length(temp_diff)]
			subdesign[location_move]=subdesign[location_move]+1
			subdesign[(location_move+1):size]=subdesign[location_move]+1:(size-location_move)
		}
		AENP=GAENP(order,subdesign)
		k=1
		while (k<=length(result)){
			if (sum(AENP_compare[k,]==AENP)==length(AENP)){
				result[[k]]=rbind(result[[k]],subdesign);break
			}
			k=k+1
		}
		if (k==(length(result)+1)){
			result=c(result,list(subdesign))
			AENP_compare=rbind(AENP_compare,AENP)
		}
	}
	temp_order=1:dim(AENP_compare)[1]
	if (length(temp_order)>1){
		for (k in 2:dim(AENP_compare)[1]){
			for (l in (k-1):1){
				temp_diff=AENP_compare[temp_order[l],]-AENP_compare[k,]
				if (sum(abs(temp_diff))>0 & temp_diff[which(temp_diff!=0)[1]]>0){
					temp_order[l+1]=temp_order[l]
					temp_order[l]=k
				}else break
			}	
		}
	}
	return(c(result,list(AENP_compare),list(temp_order)))
}

classify_G2AENP=function(order,column_stay,column_initial){
	size=length(column_initial)
	column_traversal=column_initial
	subdesign=c(column_stay,column_traversal)
	end=(length(S)-size+1):length(S)
	result=list(subdesign)
	AENP_compare=t(as.matrix(G2AENP(order,subdesign)))
	while(sum(column_traversal==end)!=size){
		if (column_traversal[size]!=length(S)) {
			column_traversal[size]=column_traversal[size]+1
		}else{
			temp_diff=which(column_traversal[-1]-column_traversal[-size]>1)
			location_move=temp_diff[length(temp_diff)]
			column_traversal[location_move]=column_traversal[location_move]+1
			column_traversal[(location_move+1):size]=column_traversal[location_move]+1:(size-location_move)
		}
		subdesign=c(column_stay,column_traversal)
		AENP=G2AENP(order,subdesign)
		k=1
		while (k<=length(result)){
			if (sum(AENP_compare[k,]==AENP)==length(AENP)){
				result[[k]]=rbind(result[[k]],subdesign);break
			}
			k=k+1
		}
		if (k==length(result)+1){
			result=c(result,list(subdesign))
			AENP_compare=rbind(AENP_compare,AENP)
		}
	}
	temp_order=1:dim(AENP_compare)[1]
	if (length(temp_order)>1){
		for (k in 2:dim(AENP_compare)[1]){
			for (l in (k-1):1){
				temp_diff=AENP_compare[temp_order[l],]-AENP_compare[k,]
				if (sum(abs(temp_diff))>0 & temp_diff[which(temp_diff!=0)[1]]>0){
					temp_order[l+1]=temp_order[l]
					temp_order[l]=k
				}else break
			}	
		}
	}
	return(c(result,list(AENP_compare),list(temp_order)))
}

order1=matrix(c(1,1),nrow=2)
order2=matrix(c(0,2,1,2,2,0,2,1,2,2),nrow=2)
order3=matrix(c(0,3,1,3,2,3,3,0,3,1,3,2,3,3),nrow=2)
order4=matrix(c(0,4,1,4,2,4,3,4,4,0,4,1,4,2,4,3,4,4),nrow=2)
order5=matrix(c(0,5,1,5,2,5,3,5,4,5,5,0,5,1,5,2,5,3,5,4,5,5),nrow=2)
order6=matrix(c(0,6,1,6,2,6,3,6,4,6,5,6,6,0,6,1,6,2,6,3,6,4,6,5,6,6),nrow=2)
order7=matrix(c(0,7,1,7,2,7,3,7,4,7,5,7,6,7,7,0,7,1,7,2,7,3,7,4,7,5,7,6,7,7),nrow=2)
order8=matrix(c(0,8,1,8,2,8,3,8,4,8,5,8,6,8,7,8,8,0,8,1,8,2,8,3,8,4,8,5,8,6,8,7,8,8),nrow=2)

order=cbind(order1,order2,order3,order4,order5,order6,order7)
a=classify_G2AENP(order,c(1,2),c(3,4,5,6,7))

b=classify_GWLP(c(1,2),c(3))
c=classify_MMA(c(1,2),c(3,4,5))


i=1;j=1;subdesign=c(1,2,3,4,5,7,9);G2AENP_C(i,j,subdesign)
GWLP(c(1,2,3,4,5,7,8,9))
MMA_moments(c(1,2,4,5,6,7,8))