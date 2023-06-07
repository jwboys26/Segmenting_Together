#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <RcppParallel.h>

using namespace RcppParallel;

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//


// [[Rcpp::export]]
double AbsVal(double x){
  
  double out = 0;
  
  if(x<0){
    out = 0-x;
  }else{
    out = x;
  }
  return out;
}

// [[Rcpp::export]]
double MinVal(double x, double y){
  
  double out=0;
  
  if(x<y){
    out = x;
  }else{
    out = y;
  }
  
  return out;
}


// [[Rcpp::export]]
double MaxVal(double x, double y){
  
  double out=0;
  
  if(x>y){
    out = x;
  }else{
    out = y;
  }
  
  return out;
}


// [[Rcpp::export]]
IntegerVector RandInts(int nMany, int ceiling) {
  
  bool bMode = FALSE;
  if(nMany > ceiling){
    bMode = TRUE;
  }
  IntegerVector results(nMany) ;
  
  IntegerVector frame = seq_len(ceiling) ;
  
  IntegerVector candidate(nMany) ;
  int maxx=ceiling+1;
  
  while (maxx > ceiling) {
    
    candidate = RcppArmadillo::sample(frame, nMany, bMode, NumericVector::create() ) ;
    
    maxx = max(candidate);
    results = candidate;
    
  }
  
  return results;
}


// [[Rcpp::export]]
List cppGenerateS12(int nx, int ny, int n1, int bFix){
  
  int Totaln = nx*ny;
  int n2 = Totaln - n1;
  
  NumericMatrix TS(Totaln, 2);
  
  int nIndex=1;
  
  for(int i=1;i<=nx;i++){
    for(int j=1;j<=ny;j++){
      nIndex = (i-1)*ny+j;
      TS(nIndex-1,0) = i;
      TS(nIndex-1,1) = j;
    }
  }
  
  
  IntegerVector IndexVec(n1);
  
  
  if(bFix == 1){
    for(int i=1;i<=n1;i++){
      IndexVec[i-1] = i;
    }
  }else{
    IndexVec = RandInts(n1, Totaln);
  }
  
  NumericMatrix S1(n1,2);
  NumericMatrix S2(n2,2);
  
  //S1.zeros();
  //S2.zeros();
  
  int nInd = 1;
  
  
  int bCheck = 0;
  int nInc = 1;
  for(int i=1;i<=Totaln;i++){
    
    for(int j=1;j<=n1;j++){
      nInd = IndexVec[j-1];
      if(i==nInd){
        bCheck=1;
        S1.row(j-1) = TS.row(nInd-1);
        break;
      }
    }
    
    if(bCheck==0){
      S2.row(nInc-1) = TS.row(i-1);
      nInc += 1;
    }
    bCheck = 0;  
  }
  
  List lst(3);
  lst[0] = S1;
  lst[1] = S2;
  lst[2] = TS;
  
  return lst;
  
  
  
}


// [[Rcpp::export]]
double cppLossFunc(arma::mat S1, arma::mat S2, arma::mat zMat, double p1, double p2){
  
  int n1 = S1.n_rows;
  int n2 = S2.n_rows;
  
  double out = 0;
  
  int xIndi, yIndi, xIndj, yIndj;
  double tmpVali, tmpValj, tmpVal;
  
  for(int i=1;i<=n1;i++){ 
    xIndi = S1(i-1, 0);
    yIndi = S1(i-1, 1);
    tmpVali = zMat(xIndi-1, yIndi-1);
    
    for(int j=1;j<=n1;j++){
      xIndj = S1(j-1,0);
      yIndj = S1(j-1,1);
      tmpValj = zMat(xIndj-1, yIndj-1);
      tmpVal = AbsVal(tmpVali+tmpValj-2*p1) - AbsVal(tmpVali-tmpValj);
      out += tmpVal;
    }
    
  }
  
  
  for(int i=1;i<=n2;i++){
    xIndi = S2(i-1, 0);
    yIndi = S2(i-1, 1);
    
    tmpVali = zMat(xIndi-1, yIndi-1);
    
    for(int j=1;j<=n2;j++){
      xIndj = S2(j-1,0);
      yIndj = S2(j-1,1);
      tmpValj = zMat(xIndj-1, yIndj-1);
      tmpVal = AbsVal(tmpVali+tmpValj-2*p2) - AbsVal(tmpVali-tmpValj);
      out += tmpVal;
    }
    
  }
  
    
  
  for(int i=1;i<=n1;i++){
    xIndi = S1(i-1, 0);
    yIndi = S1(i-1, 1);
    
    tmpVali = zMat(xIndi-1, yIndi-1);
    
    for(int j=1;j<=n2;j++){
      xIndj = S2(j-1,0);
      yIndj = S2(j-1,1);
      tmpValj = zMat(xIndj-1, yIndj-1);
      tmpVal = AbsVal(tmpVali+tmpValj-(p1+p2)) - AbsVal(tmpVali-tmpValj + (p2-p1));
      out += 2*tmpVal;
    }
    
  }
  
    
  
  return out;
  
}



struct Drudgery : public Worker{
  
  const RVector<double> G1;
  const RVector<double> G2;
  
  const std::size_t n1;
  const std::size_t n2;
  
  const RMatrix<double> zMat;
  
  const double p1;
  const double p2;
  
  RVector<double> out;
  
  
  Drudgery(const NumericVector G1, const NumericVector G2, const std::size_t n1, const std::size_t n2, const NumericMatrix zMat, const double p1, const double p2, NumericVector out)
    : G1(G1), G2(G2), n1(n1), n2(n2), zMat(zMat), p1(p1), p2(p2), out(out) {}
  
  void operator()(std::size_t begin, std::size_t end){
    
    
    for(std::size_t k = begin; k<= end; k++){
      
      double tmpVal=0;
      double tmpVal1=0;
      double tmpVal2=0;
      double gii=0;
      
      double gn1 = G1[k-1];
      
      tmpVal = 0;
      
      for(std::size_t i=1;i<=n1;i++){
        if(i!=k){
          gii = G1[i-1];
          tmpVal1 = AbsVal(gii+gn1-2*p1)-AbsVal(gii-gn1);
          tmpVal2 = AbsVal(gii+gn1-(p1+p2)) - AbsVal(gii-gn1+(p2-p1));
          tmpVal += -2*tmpVal1 + 2*tmpVal2;
          
        }
      }
      
      tmpVal += 0 - 2*AbsVal(gn1-p1);
      
      
      for(std::size_t i=1;i<=n2;i++){
        gii = G2[i-1];
        tmpVal1 = AbsVal(gii+gn1-2*p2)-AbsVal(gii-gn1);
        tmpVal2 = AbsVal(gii+gn1-(p1+p2)) - AbsVal(gn1-gii+(p2-p1));
        tmpVal += 2*tmpVal1 - 2*tmpVal2;
        
      }
      
      tmpVal += 2*AbsVal(gn1-p2);
      
      
      if(tmpVal<0){
        out[k-1] = 1;
      }else{
        out[k-1] = 0;
        
      }  
      
    }
    
    
    
  }
  
  
  
};


//[[Rcpp::export]]
NumericVector parallelRcppGrade(NumericMatrix S1, NumericMatrix S2, NumericMatrix zMat, double p1, double p2){
  
  int n1 = S1.rows();
  int n2 = S2.rows();
  
  NumericVector out(n1);
  
  NumericVector G1(n1);
  NumericVector G2(n2);
  
  int xi=0;
  int yi=0;
  
  for(int i=1;i<=n1;i++){
    xi = S1(i-1,0);
    yi = S1(i-1,1);
    G1[i-1] = zMat(xi-1,yi-1);
  }
  
  for(int i=1;i<=n2;i++){
    xi = S2(i-1,0);
    yi = S2(i-1,1);
    G2[i-1] = zMat(xi-1,yi-1);
  }
  
  Drudgery drudgery(G1, G2, n1, n2, zMat, p1, p2, out);
  
  parallelFor(1, n1, drudgery);
  
  return out;
  
}




// [[Rcpp::export]]
arma::vec WH(NumericVector xVec, double tarval){
  
  int n = xVec.length();
  
  arma::vec IndexVec(n);
  IndexVec.zeros();
  
  int nCount = 0;
  double val = 0;
  
  
  for(int i=1;i<=n;i++){
    val = xVec[i-1];
    if(val == tarval){
      IndexVec[nCount] = i;
      nCount += 1;
    }
  }
  
  if(nCount==0){
    IndexVec[0] = -1;
    return IndexVec;
  }
  
  arma::vec out = IndexVec.subvec(0, nCount-1);
  
  return out;
  
}





// [[Rcpp::export]]
arma::mat cppGet_Estimated_Img_S1(NumericMatrix S1, NumericMatrix S2, NumericMatrix zMat, double p1, double p2){
  
  
  int n1 = S1.rows();
  int n2 = S2.rows();
  
  NumericVector TF1 = parallelRcppGrade(S1, S2, zMat, p1, p2);
  NumericVector TF2 = parallelRcppGrade(S2, S1, zMat, p2, p1);
  
  int nT1 = sum(TF1);
  int nT2 = sum(TF2);
  
  int new_n1 = (n1-nT1+nT2);
  int new_n2 = (n2-nT2+nT1);
  
  int n12 = n1+n2;
  
  arma::mat out(n12+1, 2);
  out.zeros();
  
  arma::mat fakeout(1, 2);
  
  int nSuccess = 0;
  
  if(new_n1 == 0){
    nSuccess = -1;
    //out(0,0) = nSuccess;
    //out(0,1) = new_n1;
    
    fakeout(0,0) = nSuccess;
    fakeout(0,1) = new_n1;
    
    return fakeout;
  }
  if(new_n2 == 0){
    nSuccess = 0;
    // out(0,0) = nSuccess;
    // out(0,1) = n12;
    // 
    // out.rows(1, n1) = S1;
    // out.rows((n1+1), n12) = S2;
    // 
    
    fakeout(0,0) = nSuccess;
    fakeout(0,1) = new_n1;
    
    
    return fakeout;
  }
  
  nSuccess = 1;
  out(0,0) = nSuccess;
  out(0,1) = new_n1;
  
  arma::vec wh10 = WH(TF1, 0);
  //arma::vec wh11 = WH(TF1, 1);
  //arma::vec wh20 = WH(TF2, 0);
  arma::vec wh21 = WH(TF2, 1);
  
  int nInc = 1;
  int nInd = 1;
  
  if(nT1 != n1){
    
    for(int i=1;i<=(n1-nT1); i++){
      nInd = wh10[i-1];
      //out.row(nInc) = S1.row(nInd-1);
      out(nInc, 0) = S1(nInd-1,0);
      out(nInc, 1) = S1(nInd-1,1);
      
      nInc += 1;
    }
    
  }
  
  if(nT2!=0){
    
    for(int i=1;i<=nT2;i++){
      nInd = wh21[i-1];
      //out.row(nInc) = S2.row(nInd-1);
      out(nInc, 0) = S2(nInd-1,0);
      out(nInc, 1) = S2(nInd-1,1);
      
      
      nInc += 1;
      
    }
  }
  
  return out.rows(0, new_n1);
  
}





// [[Rcpp::export]]
arma::mat cppGetAllImgXY(NumericMatrix TImgMat, arma::vec StrInfo, arma::vec PatchInfo, int SI, int EI, int SJ, int EJ, double p1, double p2){
  
  int nX = TImgMat.rows();
  int nY = TImgMat.cols();
  
  arma::mat ImgMat(nX, nY);
  ImgMat.zeros();
  
  int xStr = StrInfo[0];
  int yStr = StrInfo[1];
  
  int nXb = PatchInfo[0];
  int nYb = PatchInfo[1];
  
  int xPos1, xPos2, yPos1, yPos2;
  
  int minR, minC=0;
  
  int n1 = floor(nXb*nYb/2);
  int n2 = nXb*nYb-n1;
  
  NumericMatrix zMat(nXb, nYb);
  NumericMatrix PartialzMat(nXb, nY);
  
  int bFix = 0;
  
  List lst = cppGenerateS12(nXb, nYb, n1, bFix);
  
  NumericMatrix InitS1 = lst[0];
  NumericMatrix InitS2 = lst[1];
  
  arma::mat OptOut((n1+n2+1),2);
  
  int nOptS1=0;
  int OptXpos, OptYpos;
  
  int nSuccess = 0;
  
  
  for(int i=SI;i<=EI;i++){
    xPos1 = (i-1)*xStr+1;
    xPos2 = xPos1 + nXb-1;
    
    minR = MinVal(xPos1, xPos2);
    //maxR = MaxVal(xPos1, xPos2);
    
    //PartialzMat = TImgMat.rows(minR-1, maxR-1) ; 
    
    for(int k1 = 1; k1<=nXb; k1++){
      PartialzMat.row(k1-1) = TImgMat.row(minR-1+k1-1);
    }
      
    
    for(int j=SJ;j<=EJ;j++){
      
      yPos1 = (j-1)*yStr+1;
      yPos2 = yPos1 + nYb -1;
      
      minC = MinVal(yPos1, yPos2);
      //maxC = MaxVal(yPos1, yPos2);
      
      for(int k2=1;k2<=nYb;k2++){
        zMat.column(k2-1) = PartialzMat.column(minC-1+k2-1);
      }
      
      //zMat = PartialzMat.cols(minC-1, maxC-1);
      
      OptOut = cppGet_Estimated_Img_S1(InitS1, InitS2, zMat, p1, p2);
      nSuccess = OptOut(0,0);
      
      if(nSuccess == 1){
        nOptS1 = OptOut(0,1);
        
        for(int k=1;k<=nOptS1;k++){
          OptXpos = OptOut(k,0) + minR-1;
          OptYpos = OptOut(k,1) + minC-1;
          ImgMat(OptXpos-1, OptYpos-1) += 1;
        }
        
        
      }else if(nSuccess == 0){
        
        for(int k=1;k<=n1;k++){
          OptXpos = InitS1(k-1,0) + minR-1;
          OptYpos = InitS1(k-1,1) + minC-1;
          ImgMat(OptXpos-1, OptYpos-1) += 1;
        }
        
        for(int k=1;k<=n2;k++){
          OptXpos = InitS2(k-1,0) + minR-1;
          OptYpos = InitS2(k-1,1) + minC-1;
          ImgMat(OptXpos-1, OptYpos-1) += 1;
        }
        
      }
      
      
      
    }
    
    
  }
  
  
  return ImgMat;
}







