



PaddingMat = function(OrigMat, xPad, yPad){
  
  dimm = dim(OrigMat)
  on1 = dimm[1]
  on2 = dimm[2]
  
  n1 = on1+2*xPad
  n2 = on2+2*yPad
  
  OutMat = matrix(0, n1, n2)
  
  OutMat[(xPad+1): (xPad+on1), (yPad+1): (yPad+on2)] = OrigMat
  
  return(OutMat)
}



RemovePaddingMat = function(OrigMat, xPad, yPad){
  
  dimm = dim(OrigMat)
  on1 = dimm[1]
  on2 = dimm[2]
  
  n1 = on1-2*xPad
  n2 = on2-2*yPad
  
  OutMat = matrix(0, n1, n2)
  
  OutMat = OrigMat[(xPad+1): (xPad+n1), (yPad+1): (yPad+n2)]
  
  return(OutMat)
}










Image2Mat = function(Img, hdex, vert1, vert2, b3D = TRUE){
  
  A = img[, , hdex]
  
  dimm = dim(A)
  nX = dimm[1]
  nY = dimm[2]
  
  
  
  xPos1 = vert1[1]; yPos1 = vert1[2]
  xPos2 = vert2[1]; yPos2 = vert2[2]
  
  #print(xPos1);print(xPos2)
  
  if(b3D == FALSE){
    nRow1 = nX-xPos1+1; nCol1 = yPos1 
    nRow2 = nX-xPos2+1; nCol2 = yPos2 
  }else{
    nRow1 = xPos1; nCol1 = yPos1 
    nRow2 = xPos2; nCol2 = yPos2  
  }
  
   
  minR = min(nRow1, nRow2);maxR = max(nRow1, nRow2)
  minC = min(nCol1, nCol2);maxC = max(nCol1, nCol2)
  
  print(nRow1);print(nRow2)
  
  
  Extracted_Image = A[minR:maxR, minC:maxC]
  
  return(Extracted_Image)
  
  
  
}






TSMat2Img = function(TSMat, nx, ny){
  
  totaln = dim(TSMat)[1]
  
  nInc = 1
  ImgMat = matrix(0, nx, ny)
  
  for(i in 1:nx){
    
    for(j in 1:ny){
      
      ImgMat[i,j] = TSMat[nInc, 3]
      nInc=nInc+1
    }
    
  }
  return(ImgMat)
}



CreateMatFromImg = function(ImgMat){
  
  
  dimm = dim(ImgMat)
  nx = dimm[1]
  ny = dimm[2]
  
  Totaln = nx*ny
  TSMat = matrix(0, Totaln, 3)
  
  for(i in 1:nx){
    for(j in 1:ny){
      nIndex = (i-1) * ny + j
      TSMat[nIndex, 1] = i
      TSMat[nIndex, 2] = j
    }
    
    
  }
  
  nInc=1
  
  for(i in 1:nx){
    
    for(j in 1:ny){
      val = ImgMat[i,j]
      TSMat[nInc, 3]=val
      nInc=nInc+1
    }
  }
  
  return(TSMat)
}




Img2TSMat = function(ImgMat, TSMat){
  
  dimm = dim(ImgMat)
  nx = dimm[1]
  ny = dimm[2]
  
  nInc=1
  
  for(i in 1:nx){
    
    for(j in 1:ny){
      val = ImgMat[i,j]
      TSMat[nInc, 3]=val
      nInc=nInc+1
    }
  }
  
  return(TSMat)
}



Back2Mat = function(TSMat2, OrderVec){
  
  dimm = dim(TSMat2)
  
  TSMat = matrix(0, dimm[1], dimm[2])
  
  for(i in 1:dimm[1]){
    
    idx = OrderVec[i]
    
    TSMat[idx,] = TSMat2[i,]
  }
  
  return(TSMat)
  
}


GetImage = function(img, p1, p2, nXb = 20, nYb=20, bFiltering=FALSE){
  
  
  dimImg = dim(img)
  xStr=2; yStr=2
  
  nPadding = 1;
  
  xPad = nPadding*nXb; yPad = nPadding*nYb
  
  nX = dimImg[1]+2*xPad; nY = dimImg[2]+2*yPad
  nOriginX = dimImg[1]; nOriginY = dimImg[2]
  
  DiceA = matrix(0, nOriginX, nOriginY)
  
  xQ = nXb/xStr; yQ = nYb/yStr 
  xyQ = xQ*yQ
  
  xIter = (nX-nXb)/xStr+1 ; yIter = (nY-nYb)/yStr+1
  
  SJ=1
  EJ=yIter
  
  SI=1
  EI=xIter
  
  Sy = (SJ-1)*yStr+1; Ey = (EJ-1)*yStr+nYb
  Sx = (SI-1)*xStr+1; Ex = (EI-1)*xStr+nXb
  
  
  nInc = 1
  
  TR1 = Sx;TC1 = Sy
  TR2 = Ex;TC2 = Ey
  
  TminR = min(TR1, TR2); TmaxR = max(TR1, TR2)
  TminC = min(TC1, TC2); TmaxC = max(TC1, TC2)
  
  TImg = img[1:nOriginX, 1:nOriginY]
  
  TImg_Padding = PaddingMat(TImg, xPad, yPad)
  
  tmpVec = c(TImg_Padding)
  
  StrInfo = c(xStr, yStr); PatchInfo = c(nXb, nYb)
  
  
  ImgMat = cppGetAllImgXY(TImg_Padding, StrInfo, PatchInfo, SI, EI, SJ, EJ, p1, p2)
  cVec = c(ImageMat)
  
  if(bFiltering){
    
    medVal = 0.5
    print(medVal)
    for(i in 1:nOriginX){
      
      for(j in 1:nOriginY){
        val = ImgMat[i,j]
        
        if(val < medVal){
          ImgMat[i,j]=0
        }
      }
    }
    
    
  }
  
  SegImg = 1 - RemovePaddingMat(ImgMat, xPad, yPad)/xyQ

  return(SegImg)
  
  
}


SegTogether = function(img, TSMat, p1, p2, bFiltering=FALSE){
  
  OrderVec = order(TSMat[,3])
  TSMat2 = TSMat[OrderVec,]
  TransImg = TSMat2Img(TSMat2, nX, nY)
  
  Intermediate_Img = GetImage(TransImg, p1, p2, bFiltering=bFiltering)
  
    
  TSMat3 = Img2TSMat(Intermediate_Img, TSMat2)
  BackTSMat = Back2Mat(TSMat3, OrderVec)
  BackImg = TSMat2Img(BackTSMat, nX, nY)
  
  return(BackImg)
  
  
}






