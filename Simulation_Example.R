rm(list = ls())

if(require("Rcpp")==FALSE){
  install.packages("Rcpp")
  library("Rcpp")
}

if(require("RcppArmadillo")==FALSE){
  install.packages("RcppArmadillo")
  library("RcppArmadillo")
}

if(require("RcppParallel")==FALSE){
  install.packages("RcppParallel")
  library("RcppParallel")
}

CurDir = dirname(rstudioapi::getActiveDocumentContext()$path)

cppFuncPath = paste(CurDir, "/src/SegToHer.cpp", sep="")
sourceCpp(cppFuncPath)



RFuncPath = paste(CurDir, "/src/FuncLib.R", sep="")
source(RFuncPath)

ImgPath = paste(CurDir, "/Image/ImgLib.R", sep="")
source(ImgPath)


############################ Create a 200x200 Pseudo QR image with a noise
nX=200; nY=200; r=3;        
idx=3

lst = GenerateImg(nX, nY, Type=idx, bNoise=TRUE, sig_noise=0.8)
S1 = lst[[1]]
S2 = lst[[2]]
ImageMat = lst[[3]]
TrueImageMat = lst[[4]]
TSMat = lst[[5]]

image(ImageMat, axes = FALSE, col = grey(seq(0, 1, length = 256)))      ### Plot the Pseudo QR image with a noise

image(TrueImageMat, axes = FALSE, col = grey(seq(0, 1, length = 256)))  ### Plot the Pseudo QR image without a noise
 

p1=1; p2=0
SegImg1 = GetImage(ImageMat, p1, p2)    ##### Resulting image without Segmenting-Together strategy
image(SegImg1, axes = FALSE, col = grey(seq(0, 1, length = 256)))


SegImg2 = SegTogether(ImageMat, TSMat, p1, p2)    ##### Resulting image with Segmenting-Together strategy
image(SegImg2, axes = FALSE, col = grey(seq(0, 1, length = 256)))








