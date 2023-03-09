


DiffMatrix2=function(SMat, cVec){
  
  dimS = dim(SMat)
  n = dimS[1]
  
  cLen=length(cVec)
  
  ans = matrix(0, (n-cLen),2)
  
  nInc=1
  nIndex=1
  for(i in 1:n){
    
    if(nIndex<=cLen){
      Ind = cVec[nIndex]
    }
    
    
    if(i!=Ind){
      ans[nInc,1] = SMat[i,1]
      ans[nInc,2] = SMat[i,2]
      nInc=nInc+1
      
    }else{
      nIndex=nIndex+1
    }
    
  }
  
  return(ans)
  
}








GenerateS12 = function(nx, ny, n1, bFix=TRUE, Type=1){
  
  Totaln = nx*ny
  n2 = Totaln - n1
  
  TS = matrix(0, Totaln, 2)
  
  for(i in 1:nx){
    for(j in 1:ny){
      nIndex = (i-1) * ny + j  
      TS[nIndex, 1] = i
      TS[nIndex, 2] = j
    }
    
    
  }
  

  if(bFix == TRUE){
    if(Type==1){
      IndexVec = 1:n1
    }else if(Type == 2){
      nl = floor(ny/2)
      N1 = 1
      N2 = ny+1
      N3 = 2*ny+1
      N4 = 3*ny+1
      N5 = 4*ny+1
      
      v1 = N1:(N1+nl-1)
      v2 = N2:(N2+nl-1)
      v3 = N3:(N3+nl-1)
      v4 = N4:(N4+nl-1)
      v5 = N5:(N5+nl-1)
      
      IndexVec = cbind(t(v1), t(v2), t(v3),t(v4),t(v5) )
      IndexVec = t(IndexVec)
      #print(IndexVec)
    }else if(Type==3){
      
      nl = floor(ny/2)
      nq = floor(ny/3)
      N1 = nq*ny + nq
      N2 = N1+ny
      N3 = N2+ny
      N4 = N3+ny
      N5 = N4+ny
      
      v1 = N1:(N1+nl-1)
      v2 = N2:(N2+nl-1)
      v3 = N3:(N3+nl-1)
      v4 = N4:(N4+nl-1)
      v5 = N5:(N5+nl-1)
      
      IndexVec = cbind(t(v1), t(v2), t(v3),t(v4),t(v5) )
      IndexVec = t(IndexVec)
      
      
    }else if(Type==4){
      sqn1 = floor(sqrt(n1))
      
      XCenter = floor(nx/2)+1
      YCenter = floor(ny/2)+1
      
      x1 = XCenter - floor(sqn1/2)
      y1 = YCenter - floor(sqn1/2)
      
      V = matrix(0, sqn1, 2)
      
      for(i in 1:sqn1){
        v1 = (x1+i-2)*ny + y1
        v2 = (x1+i-2)*ny + y1 + sqn1-1
        
        vec = v1:v2
        nLen = length(vec)
        
        if(i==1){
          IndexVec = vec
        }else{
          
          OldLen = length(IndexVec)
          NewIndexVec = rep(0, times=(nLen+OldLen))
          NewIndexVec[1:OldLen] = IndexVec
          NewIndexVec[(OldLen+1):(OldLen+nLen)] = vec
          IndexVec = NewIndexVec
        }
        
       
      }

    }else if(Type == 5){
      
      IndexVec = 1:ny
    }
    
  }else{
    IndexVec = sample.int(Totaln, size=n1, replace=FALSE)  
  }
  
  S1 = matrix(0, n1, 2)
  S2 = matrix(0, n2, 2)
  
  S1 = TS[IndexVec,]
  
  nInc=1
  
  
  for(i in 1:Totaln){
    

    if(length(which(IndexVec==i) )==0){
      
      S2[nInc,] = TS[i,]
  
      nInc=nInc+1
      

    }
  }
  
  lst = list()
  lst[[1]] = S1
  lst[[2]] = S2
  lst[[3]] = TS
  return(lst)
}

CheckStatus = function(x, y, xc, yc, r){
  
  r2 = r^2
  val = (x-xc)^2+(y-yc)^2
  if(val<r2){
    ans = 1
  }else{
    ans = 0
  }
  
  return(ans)
}

GenerateCircle = function(nx, ny, r){
  
  Totaln = nx*ny
  
  XCenter = floor(nx/2)+1
  YCenter = floor(ny/2)+1
  
  TS = matrix(0, Totaln, 2)
  
  nS1Inc = 1
  nS2Inc = 1
  
  S1 = matrix(0, Totaln, 2)
  S2 = matrix(0, Totaln, 2)
  
  for(i in 1:nx){
    for(j in 1:ny){
      nIndex = (i-1) * ny + j  
      TS[nIndex, 1] = i
      TS[nIndex, 2] = j
      
      bStatus = CheckStatus(i, j, XCenter, YCenter, r)
      if(bStatus == 1){
        S1[nS1Inc, 1] = i
        S1[nS1Inc, 2] = j
        nS1Inc = nS1Inc+1
      }else{
        S2[nS2Inc, 1] = i
        S2[nS2Inc, 2] = j
        nS2Inc = nS2Inc+1
      }
      
    }
  }
  
  nS1Inc = nS1Inc-1
  nS2Inc = nS2Inc-1
  
  lst = list()
  lst[[1]] = S1[1:nS1Inc, ]
  lst[[2]] = S2[1:nS2Inc, ]
  
  TS = rbind(S1[1:nS1Inc, ], S2[1:nS2Inc, ])
  
  lst[[3]] = TS
  return(lst)
}





GenImg = function(nx, ny, Type=1, bNoise=FALSE, sig_noise=0.1){
  
  Totaln = nx*ny
  n1 = floor(Totaln/2)
  n2 = Totaln - n1
  
  TS = matrix(0, Totaln, 2)
  
  for(i in 1:nx){
    for(j in 1:ny){
      nIndex = (i-1) * ny + j
      TS[nIndex, 1] = i
      TS[nIndex, 2] = j
    }
    
    
  }
  
  
  
  if(Type==1){
    
    IndexVec = c()
    
    nl = floor(ny/3)
    nq = floor(nx/3)
    
    for(i in 1:nl){
      idx = (nq+i-1)*ny + nl
      for(j in 1:nq){
        idx = idx+1
        IndexVec = c(IndexVec, idx)
      }
    }
    
    
  }else if(Type==3){
    
    IndexVec = sample.int(Totaln, size=n1, replace=FALSE)
    IndexVec = sort(IndexVec)
    
  }
  
  
  if(Type!=2){
    S1 = TS[IndexVec,]
    S2 = DiffMatrix2(TS, IndexVec)
    
  }else{
    lst = GenerateCircle(nx, ny, floor(min(nx,ny)/3))
    S1 = lst[[1]]
    S2 = lst[[2]]
  }
  
  n1=dim(S1)[1]
  
  p1 = 1
  p2 = 0
  
  
  TrueImgMat = matrix(p2, nx, ny)
  for(i in 1:n1){
    xi = S1[i,1]; yi=S1[i,2]
    TrueImgMat[xi,yi]=p1
  }
  
  Original_TrueImgMat = TrueImgMat
  
  EpsMat = matrix(rnorm(nx*ny, 0, sig_noise), nx, ny)
  
  if(bNoise==TRUE){
    TrueImgMat = TrueImgMat+EpsMat
  }
  
  
  TSMat = matrix(0, Totaln, 3)
  
  for(i in 1:nx){
    for(j in 1:ny){
      nIndex = (i-1) * ny + j
      TSMat[nIndex, 1] = i
      TSMat[nIndex, 2] = j
      TSMat[nIndex, 3] = TrueImgMat[i,j]
    }
    
  }
  
  
  ans = list(S1=S1, S2=S2, ImgMat = TrueImgMat, OriginalImgMat=Original_TrueImgMat, TSMat=TSMat)
  return(ans)
}










GenImg2 = function(nx, ny, Type=1, bNoise=FALSE, sig_noise=0.1){
  
  
  Totaln = nx*ny
  n1 = floor(Totaln/2)
  n2 = Totaln - n1
  
  TS = matrix(0, Totaln, 2)
  
  for(i in 1:nx){
    for(j in 1:ny){
      nIndex = (i-1) * ny + j
      TS[nIndex, 1] = i
      TS[nIndex, 2] = j
    }
  }
  
  IndVec=c()
  
  p1 = 1
  p2 = 0
  
  TrueImgMat = matrix(p2, nx, ny)
  
  
  c1 = floor( (nx+1)/2)
  c2 = floor( (ny+1)/2)
  
  
  if(Type == 1){
    
    lx = nx/5
    ly = ny/5
    
    
    p1x = floor(c1 - 1.5*lx)
    p1y = floor(c2 - 0.5*ly)
    
    p2x = floor(c1 - 1.5*lx)
    p2y = floor(c2 + 0.5*ly)
    
    p3x = floor(c1 - 0.5*lx)
    p3y = floor(c2 - 1.5*ly)
    
    p4x = floor(c1 - 0.5*lx)
    p4y = floor(c2 - 0.5*ly)
    
    p5x = floor(c1 - 0.5*lx)
    p5y = floor(c2 + 0.5*ly)
    
    p6x = floor(c1 - 0.5*lx)
    p6y = floor(c2 + 1.5*ly)
    
    p7x = floor(c1 + 0.5*lx)
    p7y = floor(c2 - 1.5*ly)
    
    p8x = floor(c1 + 0.5*lx)
    p8y = floor(c2 - 0.5*ly)
    
    p9x = floor(c1 + 0.5*lx)
    p9y = floor(c2 + 0.5*ly)
    
    p10x = floor(c1 + 0.5*lx)
    p10y = floor(c2 + 1.5*ly)
    
    p11x = floor(c1 + 1.5*lx)
    p11y = floor(c2 - 0.5*ly)
    
    p12x = floor(c1 + 1.5*lx)
    p12y = floor(c2 + 0.5*ly)
    
    
    for(i in p1x: (p3x-1)){
      
      for(j in p1y:p2y){
        TrueImgMat[i,j]=p1
        
        nIndex = (i-1) * ny + j
        IndVec=c(IndVec, nIndex)
        
      }
      
    }
    
    for(i in p3x:p7x){
      
      for(j in p3y:p6y){
        TrueImgMat[i,j]=p1
        
        nIndex = (i-1) * ny + j
        IndVec=c(IndVec, nIndex)
        
        
      }
      
    }
    
    for(i in (p7x+1):p11x){
      
      for(j in p11y:p12y){
        TrueImgMat[i,j]=p1
        
        nIndex = (i-1) * ny + j
        IndVec=c(IndVec, nIndex)
        
        
      }
      
    }
    
    
    
  }else if(Type==2){
    
    r=nx/3
    
    p1x= floor(c1-r)
    p1y= c2
    
    p2x = floor(c1+r/2)
    p2y = floor(c2-sqrt(3)*r/2)
    
    p3x = floor(c1+r/2)
    p3y = floor(c2+sqrt(3)*r/2)
    
    
    for(i in p1x:p2x){
      
      tp2y = floor(c2-(i-p1x)/sqrt(3))
      tp3y = floor(c2+(i-p1x)/sqrt(3))
      
      for(j in tp2y:tp3y){
        
        TrueImgMat[i,j]=p1
        
        nIndex = (i-1) * ny + j
        IndVec=c(IndVec, nIndex)
        
        
      }
      
    }
    
  }else if(Type==3){
    
    r=nx/3
    
    p1x= floor(c1-r)
    p1y= c2
    
    p2x = floor(c1+r/2)
    p2y = floor(c2-sqrt(3)*r/2)
    
    p3x = floor(c1+r/2)
    p3y = floor(c2+sqrt(3)*r/2)
    
    
    for(i in p1x:p2x){
      
      tp2y = floor(c2-(i-p1x)/sqrt(3))
      tp3y = floor(c2+(i-p1x)/sqrt(3))
      
      for(j in tp2y:tp3y){
        TrueImgMat[j,nx-i+1]=p1
        TrueImgMat[j,i]=p1
        
        nIndex = (i-1) * ny + j
        IndVec=c(IndVec, nIndex)
        
        
      }
      
    }
    
  }
  
  
  S1 = TS[IndVec,]
  S2 = DiffMatrix2(TS, IndVec)
  
  Original_TrueImgMat = TrueImgMat
  
  EpsMat = matrix(rnorm(nx*ny, 0, sig_noise), nx, ny)
  
  if(bNoise==TRUE){
    TrueImgMat = TrueImgMat+EpsMat
  }
  
  
  TSMat = matrix(0, Totaln, 3)
  
  
  for(i in 1:nx){
    for(j in 1:ny){
      nIndex = (i-1) * ny + j
      TSMat[nIndex, 1] = i
      TSMat[nIndex, 2] = j
      TSMat[nIndex, 3] = TrueImgMat[i,j]
    }
    
  }
  
  
  
  ans = list(S1=S1, S2=S2, ImgMat = TrueImgMat, 
             OriginalMat = Original_TrueImgMat, TSMat=TSMat)
  return(ans)
  
}



GenerateImg = function(nx, ny, Type=1, bNoise=FALSE, sig_noise=0.1){
  
  if(Type<=3){
    
    lst = GenImg(nX, nY, Type=Type, bNoise=bNoise, sig_noise=sig_noise)
    
  }else{
    
    Type=Type-3 
    lst = GenImg2(nX, nY, Type=Type, bNoise=bNoise, sig_noise=sig_noise)
  }
  
  return(lst)
  
}







