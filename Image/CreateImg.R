
rm(list=ls())

CurDir = dirname(rstudioapi::getActiveDocumentContext()$path)

FuncPath = paste(CurDir, "/ImgLib.R", sep="")
source(FuncPath)


nX=200; nY=200  ##### Size of a generated image      


idx=1            ##### 1: square, 2: circle, 3: QR-image, 4: cross, 5: star, 6: triangle
bNoise = TRUE    ##### TRUE: Image with noises, FALSE: Image without noises
sig_noise = 0.6  ##### When bNoise is TRUE, specify the strength of noises. Noises will be generated from N(0, sig_noise).

lst = GenerateImg(nX, nY, Type=idx, bNoise=bNoise, sig_noise=sig_noise)     

S1 = lst[[1]]                #### an n1-by-2 matrix of entries of n1 white pixels: n1+n2 = nX x nY
S2 = lst[[2]]                #### an n2-by-2 matrix of entries of n2 black pixels: n1+n2 = nX x nY
ImageMat = lst[[3]]          #### an nX-by-nY matrix of the pixel values of the generated image with noises
TrueImageMat = lst[[4]]      #### an nX-by-nY matrix of the pixel values of the generated image without noises
TSMat = lst[[5]]             #### an (nX x nY)-by-3 matrix: the first and second columns are entries of pixels 
                             #### while the third column is the corresponding pixel value.   


#### png of true Image 
image(ImageMat, axes = FALSE, col = grey(seq(0, 1, length = 256)))

#### png of image with noises
image(TrueImageMat, axes = FALSE, col = grey(seq(0, 1, length = 256)))
