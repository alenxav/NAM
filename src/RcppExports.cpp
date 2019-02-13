// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// emBA
SEXP emBA(NumericVector y, NumericMatrix gen, double df, double R2);
RcppExport SEXP _NAM_emBA(SEXP ySEXP, SEXP genSEXP, SEXP dfSEXP, SEXP R2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type gen(genSEXP);
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    Rcpp::traits::input_parameter< double >::type R2(R2SEXP);
    rcpp_result_gen = Rcpp::wrap(emBA(y, gen, df, R2));
    return rcpp_result_gen;
END_RCPP
}
// emBB
SEXP emBB(NumericVector y, NumericMatrix gen, double df, double R2, double Pi);
RcppExport SEXP _NAM_emBB(SEXP ySEXP, SEXP genSEXP, SEXP dfSEXP, SEXP R2SEXP, SEXP PiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type gen(genSEXP);
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    Rcpp::traits::input_parameter< double >::type R2(R2SEXP);
    Rcpp::traits::input_parameter< double >::type Pi(PiSEXP);
    rcpp_result_gen = Rcpp::wrap(emBB(y, gen, df, R2, Pi));
    return rcpp_result_gen;
END_RCPP
}
// emBC
SEXP emBC(NumericVector y, NumericMatrix gen, double df, double R2, double Pi);
RcppExport SEXP _NAM_emBC(SEXP ySEXP, SEXP genSEXP, SEXP dfSEXP, SEXP R2SEXP, SEXP PiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type gen(genSEXP);
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    Rcpp::traits::input_parameter< double >::type R2(R2SEXP);
    Rcpp::traits::input_parameter< double >::type Pi(PiSEXP);
    rcpp_result_gen = Rcpp::wrap(emBC(y, gen, df, R2, Pi));
    return rcpp_result_gen;
END_RCPP
}
// emRR
SEXP emRR(NumericVector y, NumericMatrix gen, double df, double R2);
RcppExport SEXP _NAM_emRR(SEXP ySEXP, SEXP genSEXP, SEXP dfSEXP, SEXP R2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type gen(genSEXP);
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    Rcpp::traits::input_parameter< double >::type R2(R2SEXP);
    rcpp_result_gen = Rcpp::wrap(emRR(y, gen, df, R2));
    return rcpp_result_gen;
END_RCPP
}
// emBL
SEXP emBL(NumericVector y, NumericMatrix gen, double R2, double alpha);
RcppExport SEXP _NAM_emBL(SEXP ySEXP, SEXP genSEXP, SEXP R2SEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type gen(genSEXP);
    Rcpp::traits::input_parameter< double >::type R2(R2SEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(emBL(y, gen, R2, alpha));
    return rcpp_result_gen;
END_RCPP
}
// emDE
SEXP emDE(NumericVector y, NumericMatrix gen, double R2);
RcppExport SEXP _NAM_emDE(SEXP ySEXP, SEXP genSEXP, SEXP R2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type gen(genSEXP);
    Rcpp::traits::input_parameter< double >::type R2(R2SEXP);
    rcpp_result_gen = Rcpp::wrap(emDE(y, gen, R2));
    return rcpp_result_gen;
END_RCPP
}
// emEN
SEXP emEN(NumericVector y, NumericMatrix gen, double R2, double alpha);
RcppExport SEXP _NAM_emEN(SEXP ySEXP, SEXP genSEXP, SEXP R2SEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type gen(genSEXP);
    Rcpp::traits::input_parameter< double >::type R2(R2SEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(emEN(y, gen, R2, alpha));
    return rcpp_result_gen;
END_RCPP
}
// emML
SEXP emML(NumericVector y, NumericMatrix gen, Rcpp::Nullable<Rcpp::NumericVector> D);
RcppExport SEXP _NAM_emML(SEXP ySEXP, SEXP genSEXP, SEXP DSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type gen(genSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type D(DSEXP);
    rcpp_result_gen = Rcpp::wrap(emML(y, gen, D));
    return rcpp_result_gen;
END_RCPP
}
// emMX
SEXP emMX(NumericVector y, NumericMatrix gen, double R2);
RcppExport SEXP _NAM_emMX(SEXP ySEXP, SEXP genSEXP, SEXP R2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type gen(genSEXP);
    Rcpp::traits::input_parameter< double >::type R2(R2SEXP);
    rcpp_result_gen = Rcpp::wrap(emMX(y, gen, R2));
    return rcpp_result_gen;
END_RCPP
}
// calcSize
int calcSize(NumericVector col, NumericVector fam);
RcppExport SEXP _NAM_calcSize(SEXP colSEXP, SEXP famSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type col(colSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type fam(famSEXP);
    rcpp_result_gen = Rcpp::wrap(calcSize(col, fam));
    return rcpp_result_gen;
END_RCPP
}
// funI
NumericVector funI(NumericVector col, int fam, int finsiz, int f);
RcppExport SEXP _NAM_funI(SEXP colSEXP, SEXP famSEXP, SEXP finsizSEXP, SEXP fSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type col(colSEXP);
    Rcpp::traits::input_parameter< int >::type fam(famSEXP);
    Rcpp::traits::input_parameter< int >::type finsiz(finsizSEXP);
    Rcpp::traits::input_parameter< int >::type f(fSEXP);
    rcpp_result_gen = Rcpp::wrap(funI(col, fam, finsiz, f));
    return rcpp_result_gen;
END_RCPP
}
// funX
NumericVector funX(NumericVector col, int finsiz);
RcppExport SEXP _NAM_funX(SEXP colSEXP, SEXP finsizSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type col(colSEXP);
    Rcpp::traits::input_parameter< int >::type finsiz(finsizSEXP);
    rcpp_result_gen = Rcpp::wrap(funX(col, finsiz));
    return rcpp_result_gen;
END_RCPP
}
// gs
void gs(NumericMatrix C, NumericVector g, NumericVector r, int N);
RcppExport SEXP _NAM_gs(SEXP CSEXP, SEXP gSEXP, SEXP rSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type C(CSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type g(gSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    gs(C, g, r, N);
    return R_NilValue;
END_RCPP
}
// inputRow
NumericVector inputRow(NumericVector x, int exp, int n);
RcppExport SEXP _NAM_inputRow(SEXP xSEXP, SEXP expSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type exp(expSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(inputRow(x, exp, n));
    return rcpp_result_gen;
END_RCPP
}
// KMUP
SEXP KMUP(NumericMatrix X, NumericVector b, NumericVector d, NumericVector xx, NumericVector e, NumericVector L, double Ve, double pi);
RcppExport SEXP _NAM_KMUP(SEXP XSEXP, SEXP bSEXP, SEXP dSEXP, SEXP xxSEXP, SEXP eSEXP, SEXP LSEXP, SEXP VeSEXP, SEXP piSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type d(dSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type e(eSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type L(LSEXP);
    Rcpp::traits::input_parameter< double >::type Ve(VeSEXP);
    Rcpp::traits::input_parameter< double >::type pi(piSEXP);
    rcpp_result_gen = Rcpp::wrap(KMUP(X, b, d, xx, e, L, Ve, pi));
    return rcpp_result_gen;
END_RCPP
}
// KMUP2
SEXP KMUP2(NumericMatrix X, NumericVector Use, NumericVector b, NumericVector d, NumericVector xx, NumericVector E, NumericVector L, double Ve, double pi);
RcppExport SEXP _NAM_KMUP2(SEXP XSEXP, SEXP UseSEXP, SEXP bSEXP, SEXP dSEXP, SEXP xxSEXP, SEXP ESEXP, SEXP LSEXP, SEXP VeSEXP, SEXP piSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Use(UseSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type d(dSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type E(ESEXP);
    Rcpp::traits::input_parameter< NumericVector >::type L(LSEXP);
    Rcpp::traits::input_parameter< double >::type Ve(VeSEXP);
    Rcpp::traits::input_parameter< double >::type pi(piSEXP);
    rcpp_result_gen = Rcpp::wrap(KMUP2(X, Use, b, d, xx, E, L, Ve, pi));
    return rcpp_result_gen;
END_RCPP
}
// SAMP
void SAMP(NumericMatrix C, NumericVector g, NumericVector r, int N, double Ve);
RcppExport SEXP _NAM_SAMP(SEXP CSEXP, SEXP gSEXP, SEXP rSEXP, SEXP NSEXP, SEXP VeSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type C(CSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type g(gSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type Ve(VeSEXP);
    SAMP(C, g, r, N, Ve);
    return R_NilValue;
END_RCPP
}
// SAMP2
void SAMP2(NumericMatrix X, NumericVector g, NumericVector y, NumericVector xx, NumericVector E, NumericVector L, int N, double Ve);
RcppExport SEXP _NAM_SAMP2(SEXP XSEXP, SEXP gSEXP, SEXP ySEXP, SEXP xxSEXP, SEXP ESEXP, SEXP LSEXP, SEXP NSEXP, SEXP VeSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type g(gSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type E(ESEXP);
    Rcpp::traits::input_parameter< NumericVector >::type L(LSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type Ve(VeSEXP);
    SAMP2(X, g, y, xx, E, L, N, Ve);
    return R_NilValue;
END_RCPP
}
// timesMatrix
NumericMatrix timesMatrix(NumericMatrix ma1, NumericVector h, NumericMatrix ma2, int rows, int cols);
RcppExport SEXP _NAM_timesMatrix(SEXP ma1SEXP, SEXP hSEXP, SEXP ma2SEXP, SEXP rowsSEXP, SEXP colsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type ma1(ma1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type h(hSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type ma2(ma2SEXP);
    Rcpp::traits::input_parameter< int >::type rows(rowsSEXP);
    Rcpp::traits::input_parameter< int >::type cols(colsSEXP);
    rcpp_result_gen = Rcpp::wrap(timesMatrix(ma1, h, ma2, rows, cols));
    return rcpp_result_gen;
END_RCPP
}
// timesVec
NumericMatrix timesVec(NumericVector aa, NumericVector h, NumericMatrix bb, int q);
RcppExport SEXP _NAM_timesVec(SEXP aaSEXP, SEXP hSEXP, SEXP bbSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type aa(aaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type h(hSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type bb(bbSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(timesVec(aa, h, bb, q));
    return rcpp_result_gen;
END_RCPP
}
// CNT
NumericMatrix CNT(NumericMatrix X);
RcppExport SEXP _NAM_CNT(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(CNT(X));
    return rcpp_result_gen;
END_RCPP
}
// MSX
SEXP MSX(NumericMatrix X);
RcppExport SEXP _NAM_MSX(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(MSX(X));
    return rcpp_result_gen;
END_RCPP
}
// IMP
NumericMatrix IMP(NumericMatrix X);
RcppExport SEXP _NAM_IMP(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(IMP(X));
    return rcpp_result_gen;
END_RCPP
}
// NOR
SEXP NOR(NumericVector y, NumericMatrix X, double cxx, NumericVector xx, int maxit, double tol);
RcppExport SEXP _NAM_NOR(SEXP ySEXP, SEXP XSEXP, SEXP cxxSEXP, SEXP xxSEXP, SEXP maxitSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type cxx(cxxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(NOR(y, X, cxx, xx, maxit, tol));
    return rcpp_result_gen;
END_RCPP
}
// GAU
NumericMatrix GAU(NumericMatrix X);
RcppExport SEXP _NAM_GAU(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(GAU(X));
    return rcpp_result_gen;
END_RCPP
}
// GRM
NumericMatrix GRM(NumericMatrix X, bool Code012);
RcppExport SEXP _NAM_GRM(SEXP XSEXP, SEXP Code012SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< bool >::type Code012(Code012SEXP);
    rcpp_result_gen = Rcpp::wrap(GRM(X, Code012));
    return rcpp_result_gen;
END_RCPP
}
// SPC
NumericVector SPC(NumericVector y, NumericVector blk, NumericVector row, NumericVector col, double rN, double cN);
RcppExport SEXP _NAM_SPC(SEXP ySEXP, SEXP blkSEXP, SEXP rowSEXP, SEXP colSEXP, SEXP rNSEXP, SEXP cNSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type blk(blkSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type row(rowSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type col(colSEXP);
    Rcpp::traits::input_parameter< double >::type rN(rNSEXP);
    Rcpp::traits::input_parameter< double >::type cN(cNSEXP);
    rcpp_result_gen = Rcpp::wrap(SPC(y, blk, row, col, rN, cN));
    return rcpp_result_gen;
END_RCPP
}
// SPM
NumericMatrix SPM(NumericVector blk, NumericVector row, NumericVector col, double rN, double cN);
RcppExport SEXP _NAM_SPM(SEXP blkSEXP, SEXP rowSEXP, SEXP colSEXP, SEXP rNSEXP, SEXP cNSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type blk(blkSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type row(rowSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type col(colSEXP);
    Rcpp::traits::input_parameter< double >::type rN(rNSEXP);
    Rcpp::traits::input_parameter< double >::type cN(cNSEXP);
    rcpp_result_gen = Rcpp::wrap(SPM(blk, row, col, rN, cN));
    return rcpp_result_gen;
END_RCPP
}
// BRR2
SEXP BRR2(NumericVector y, NumericMatrix X1, NumericMatrix X2, double it, double bi, double df, double R2);
RcppExport SEXP _NAM_BRR2(SEXP ySEXP, SEXP X1SEXP, SEXP X2SEXP, SEXP itSEXP, SEXP biSEXP, SEXP dfSEXP, SEXP R2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X1(X1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X2(X2SEXP);
    Rcpp::traits::input_parameter< double >::type it(itSEXP);
    Rcpp::traits::input_parameter< double >::type bi(biSEXP);
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    Rcpp::traits::input_parameter< double >::type R2(R2SEXP);
    rcpp_result_gen = Rcpp::wrap(BRR2(y, X1, X2, it, bi, df, R2));
    return rcpp_result_gen;
END_RCPP
}
// emGWA
SEXP emGWA(NumericVector y, NumericMatrix gen);
RcppExport SEXP _NAM_emGWA(SEXP ySEXP, SEXP genSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type gen(genSEXP);
    rcpp_result_gen = Rcpp::wrap(emGWA(y, gen));
    return rcpp_result_gen;
END_RCPP
}
// BCpi
SEXP BCpi(NumericVector y, NumericMatrix X, double it, double bi, double df, double R2);
RcppExport SEXP _NAM_BCpi(SEXP ySEXP, SEXP XSEXP, SEXP itSEXP, SEXP biSEXP, SEXP dfSEXP, SEXP R2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type it(itSEXP);
    Rcpp::traits::input_parameter< double >::type bi(biSEXP);
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    Rcpp::traits::input_parameter< double >::type R2(R2SEXP);
    rcpp_result_gen = Rcpp::wrap(BCpi(y, X, it, bi, df, R2));
    return rcpp_result_gen;
END_RCPP
}
// emML2
SEXP emML2(NumericVector y, NumericMatrix X1, NumericMatrix X2, Rcpp::Nullable<Rcpp::NumericVector> D1, Rcpp::Nullable<Rcpp::NumericVector> D2);
RcppExport SEXP _NAM_emML2(SEXP ySEXP, SEXP X1SEXP, SEXP X2SEXP, SEXP D1SEXP, SEXP D2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X1(X1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X2(X2SEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type D1(D1SEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type D2(D2SEXP);
    rcpp_result_gen = Rcpp::wrap(emML2(y, X1, X2, D1, D2));
    return rcpp_result_gen;
END_RCPP
}
// mrr
SEXP mrr(NumericMatrix Y, NumericMatrix X);
RcppExport SEXP _NAM_mrr(SEXP YSEXP, SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(mrr(Y, X));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP NAM_calcSize(SEXP, SEXP);
RcppExport SEXP NAM_CNT(SEXP);
RcppExport SEXP NAM_emBA(SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP NAM_emBB(SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP NAM_emBC(SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP NAM_emBL(SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP NAM_emDE(SEXP, SEXP, SEXP);
RcppExport SEXP NAM_emEN(SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP NAM_emRR(SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP NAM_funI(SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP NAM_funX(SEXP, SEXP);
RcppExport SEXP NAM_GAU(SEXP);
RcppExport SEXP NAM_gs(SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP NAM_IMP(SEXP);
RcppExport SEXP NAM_inputRow(SEXP, SEXP, SEXP);
RcppExport SEXP NAM_KMUP(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP NAM_KMUP2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP NAM_MSX(SEXP);
RcppExport SEXP NAM_NOR(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP NAM_SAMP(SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP NAM_SAMP2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP NAM_SPC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP NAM_timesMatrix(SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP NAM_timesVec(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_NAM_emBA", (DL_FUNC) &_NAM_emBA, 4},
    {"_NAM_emBB", (DL_FUNC) &_NAM_emBB, 5},
    {"_NAM_emBC", (DL_FUNC) &_NAM_emBC, 5},
    {"_NAM_emRR", (DL_FUNC) &_NAM_emRR, 4},
    {"_NAM_emBL", (DL_FUNC) &_NAM_emBL, 4},
    {"_NAM_emDE", (DL_FUNC) &_NAM_emDE, 3},
    {"_NAM_emEN", (DL_FUNC) &_NAM_emEN, 4},
    {"_NAM_emML", (DL_FUNC) &_NAM_emML, 3},
    {"_NAM_emMX", (DL_FUNC) &_NAM_emMX, 3},
    {"_NAM_calcSize", (DL_FUNC) &_NAM_calcSize, 2},
    {"_NAM_funI", (DL_FUNC) &_NAM_funI, 4},
    {"_NAM_funX", (DL_FUNC) &_NAM_funX, 2},
    {"_NAM_gs", (DL_FUNC) &_NAM_gs, 4},
    {"_NAM_inputRow", (DL_FUNC) &_NAM_inputRow, 3},
    {"_NAM_KMUP", (DL_FUNC) &_NAM_KMUP, 8},
    {"_NAM_KMUP2", (DL_FUNC) &_NAM_KMUP2, 9},
    {"_NAM_SAMP", (DL_FUNC) &_NAM_SAMP, 5},
    {"_NAM_SAMP2", (DL_FUNC) &_NAM_SAMP2, 8},
    {"_NAM_timesMatrix", (DL_FUNC) &_NAM_timesMatrix, 5},
    {"_NAM_timesVec", (DL_FUNC) &_NAM_timesVec, 4},
    {"_NAM_CNT", (DL_FUNC) &_NAM_CNT, 1},
    {"_NAM_MSX", (DL_FUNC) &_NAM_MSX, 1},
    {"_NAM_IMP", (DL_FUNC) &_NAM_IMP, 1},
    {"_NAM_NOR", (DL_FUNC) &_NAM_NOR, 6},
    {"_NAM_GAU", (DL_FUNC) &_NAM_GAU, 1},
    {"_NAM_GRM", (DL_FUNC) &_NAM_GRM, 2},
    {"_NAM_SPC", (DL_FUNC) &_NAM_SPC, 6},
    {"_NAM_SPM", (DL_FUNC) &_NAM_SPM, 5},
    {"_NAM_BRR2", (DL_FUNC) &_NAM_BRR2, 7},
    {"_NAM_emGWA", (DL_FUNC) &_NAM_emGWA, 2},
    {"_NAM_BCpi", (DL_FUNC) &_NAM_BCpi, 6},
    {"_NAM_emML2", (DL_FUNC) &_NAM_emML2, 5},
    {"_NAM_mrr", (DL_FUNC) &_NAM_mrr, 2},
    {"NAM_calcSize",    (DL_FUNC) &NAM_calcSize,    2},
    {"NAM_CNT",         (DL_FUNC) &NAM_CNT,         1},
    {"NAM_emBA",        (DL_FUNC) &NAM_emBA,        4},
    {"NAM_emBB",        (DL_FUNC) &NAM_emBB,        5},
    {"NAM_emBC",        (DL_FUNC) &NAM_emBC,        5},
    {"NAM_emBL",        (DL_FUNC) &NAM_emBL,        4},
    {"NAM_emDE",        (DL_FUNC) &NAM_emDE,        3},
    {"NAM_emEN",        (DL_FUNC) &NAM_emEN,        4},
    {"NAM_emRR",        (DL_FUNC) &NAM_emRR,        4},
    {"NAM_funI",        (DL_FUNC) &NAM_funI,        4},
    {"NAM_funX",        (DL_FUNC) &NAM_funX,        2},
    {"NAM_GAU",         (DL_FUNC) &NAM_GAU,         1},
    {"NAM_gs",          (DL_FUNC) &NAM_gs,          4},
    {"NAM_IMP",         (DL_FUNC) &NAM_IMP,         1},
    {"NAM_inputRow",    (DL_FUNC) &NAM_inputRow,    3},
    {"NAM_KMUP",        (DL_FUNC) &NAM_KMUP,        8},
    {"NAM_KMUP2",       (DL_FUNC) &NAM_KMUP2,       9},
    {"NAM_MSX",         (DL_FUNC) &NAM_MSX,         1},
    {"NAM_NOR",         (DL_FUNC) &NAM_NOR,         6},
    {"NAM_SAMP",        (DL_FUNC) &NAM_SAMP,        5},
    {"NAM_SAMP2",       (DL_FUNC) &NAM_SAMP2,       8},
    {"NAM_SPC",         (DL_FUNC) &NAM_SPC,         6},
    {"NAM_timesMatrix", (DL_FUNC) &NAM_timesMatrix, 5},
    {"NAM_timesVec",    (DL_FUNC) &NAM_timesVec,    4},
    {NULL, NULL, 0}
};

RcppExport void R_init_NAM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
