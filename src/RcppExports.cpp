// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// HiddenBASiCS_MCMCcpp
Rcpp::List HiddenBASiCS_MCMCcpp(int N, int thin, int burn, NumericMatrix Counts, NumericVector mu0, NumericVector delta0, NumericVector phi0, NumericVector s0, NumericVector nu0, double theta0, double s2mu, double adelta, double bdelta, NumericVector p_Phi, double as, double bs, double atheta, double btheta, double ar, NumericVector LSmu0, NumericVector LSdelta0, double LSphi0, NumericVector LSnu0, double LStheta0, NumericVector sumByCellAll, NumericVector sumByCellBio, NumericVector sumByGeneAll, NumericVector sumByGeneBio, int StoreAdapt, int EndAdapt, int PrintProgress, double s2_delta, double prior_delta);
RcppExport SEXP BASiCS_HiddenBASiCS_MCMCcpp(SEXP NSEXP, SEXP thinSEXP, SEXP burnSEXP, SEXP CountsSEXP, SEXP mu0SEXP, SEXP delta0SEXP, SEXP phi0SEXP, SEXP s0SEXP, SEXP nu0SEXP, SEXP theta0SEXP, SEXP s2muSEXP, SEXP adeltaSEXP, SEXP bdeltaSEXP, SEXP p_PhiSEXP, SEXP asSEXP, SEXP bsSEXP, SEXP athetaSEXP, SEXP bthetaSEXP, SEXP arSEXP, SEXP LSmu0SEXP, SEXP LSdelta0SEXP, SEXP LSphi0SEXP, SEXP LSnu0SEXP, SEXP LStheta0SEXP, SEXP sumByCellAllSEXP, SEXP sumByCellBioSEXP, SEXP sumByGeneAllSEXP, SEXP sumByGeneBioSEXP, SEXP StoreAdaptSEXP, SEXP EndAdaptSEXP, SEXP PrintProgressSEXP, SEXP s2_deltaSEXP, SEXP prior_deltaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< int >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Counts(CountsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu0(mu0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type delta0(delta0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi0(phi0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s0(s0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nu0(nu0SEXP);
    Rcpp::traits::input_parameter< double >::type theta0(theta0SEXP);
    Rcpp::traits::input_parameter< double >::type s2mu(s2muSEXP);
    Rcpp::traits::input_parameter< double >::type adelta(adeltaSEXP);
    Rcpp::traits::input_parameter< double >::type bdelta(bdeltaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p_Phi(p_PhiSEXP);
    Rcpp::traits::input_parameter< double >::type as(asSEXP);
    Rcpp::traits::input_parameter< double >::type bs(bsSEXP);
    Rcpp::traits::input_parameter< double >::type atheta(athetaSEXP);
    Rcpp::traits::input_parameter< double >::type btheta(bthetaSEXP);
    Rcpp::traits::input_parameter< double >::type ar(arSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type LSmu0(LSmu0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type LSdelta0(LSdelta0SEXP);
    Rcpp::traits::input_parameter< double >::type LSphi0(LSphi0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type LSnu0(LSnu0SEXP);
    Rcpp::traits::input_parameter< double >::type LStheta0(LStheta0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sumByCellAll(sumByCellAllSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sumByCellBio(sumByCellBioSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sumByGeneAll(sumByGeneAllSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sumByGeneBio(sumByGeneBioSEXP);
    Rcpp::traits::input_parameter< int >::type StoreAdapt(StoreAdaptSEXP);
    Rcpp::traits::input_parameter< int >::type EndAdapt(EndAdaptSEXP);
    Rcpp::traits::input_parameter< int >::type PrintProgress(PrintProgressSEXP);
    Rcpp::traits::input_parameter< double >::type s2_delta(s2_deltaSEXP);
    Rcpp::traits::input_parameter< double >::type prior_delta(prior_deltaSEXP);
    __result = Rcpp::wrap(HiddenBASiCS_MCMCcpp(N, thin, burn, Counts, mu0, delta0, phi0, s0, nu0, theta0, s2mu, adelta, bdelta, p_Phi, as, bs, atheta, btheta, ar, LSmu0, LSdelta0, LSphi0, LSnu0, LStheta0, sumByCellAll, sumByCellBio, sumByGeneAll, sumByGeneBio, StoreAdapt, EndAdapt, PrintProgress, s2_delta, prior_delta));
    return __result;
END_RCPP
}
// HiddenBASiCS_MCMCcppBatch
Rcpp::List HiddenBASiCS_MCMCcppBatch(int N, int thin, int burn, NumericMatrix Counts, NumericMatrix BatchDesign, NumericVector mu0, NumericVector delta0, NumericVector phi0, NumericVector s0, NumericVector nu0, double theta0, double s2mu, double adelta, double bdelta, NumericVector p_Phi, double as, double bs, double atheta, double btheta, double ar, NumericVector LSmu0, NumericVector LSdelta0, double LSphi0, NumericVector LSnu0, double LStheta0, NumericVector sumByCellAll, NumericVector sumByCellBio, NumericVector sumByGeneAll, NumericVector sumByGeneBio, int StoreAdapt, int EndAdapt, int PrintProgress, double s2_delta, double prior_delta);
RcppExport SEXP BASiCS_HiddenBASiCS_MCMCcppBatch(SEXP NSEXP, SEXP thinSEXP, SEXP burnSEXP, SEXP CountsSEXP, SEXP BatchDesignSEXP, SEXP mu0SEXP, SEXP delta0SEXP, SEXP phi0SEXP, SEXP s0SEXP, SEXP nu0SEXP, SEXP theta0SEXP, SEXP s2muSEXP, SEXP adeltaSEXP, SEXP bdeltaSEXP, SEXP p_PhiSEXP, SEXP asSEXP, SEXP bsSEXP, SEXP athetaSEXP, SEXP bthetaSEXP, SEXP arSEXP, SEXP LSmu0SEXP, SEXP LSdelta0SEXP, SEXP LSphi0SEXP, SEXP LSnu0SEXP, SEXP LStheta0SEXP, SEXP sumByCellAllSEXP, SEXP sumByCellBioSEXP, SEXP sumByGeneAllSEXP, SEXP sumByGeneBioSEXP, SEXP StoreAdaptSEXP, SEXP EndAdaptSEXP, SEXP PrintProgressSEXP, SEXP s2_deltaSEXP, SEXP prior_deltaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< int >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Counts(CountsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type BatchDesign(BatchDesignSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu0(mu0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type delta0(delta0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi0(phi0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s0(s0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nu0(nu0SEXP);
    Rcpp::traits::input_parameter< double >::type theta0(theta0SEXP);
    Rcpp::traits::input_parameter< double >::type s2mu(s2muSEXP);
    Rcpp::traits::input_parameter< double >::type adelta(adeltaSEXP);
    Rcpp::traits::input_parameter< double >::type bdelta(bdeltaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p_Phi(p_PhiSEXP);
    Rcpp::traits::input_parameter< double >::type as(asSEXP);
    Rcpp::traits::input_parameter< double >::type bs(bsSEXP);
    Rcpp::traits::input_parameter< double >::type atheta(athetaSEXP);
    Rcpp::traits::input_parameter< double >::type btheta(bthetaSEXP);
    Rcpp::traits::input_parameter< double >::type ar(arSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type LSmu0(LSmu0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type LSdelta0(LSdelta0SEXP);
    Rcpp::traits::input_parameter< double >::type LSphi0(LSphi0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type LSnu0(LSnu0SEXP);
    Rcpp::traits::input_parameter< double >::type LStheta0(LStheta0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sumByCellAll(sumByCellAllSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sumByCellBio(sumByCellBioSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sumByGeneAll(sumByGeneAllSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sumByGeneBio(sumByGeneBioSEXP);
    Rcpp::traits::input_parameter< int >::type StoreAdapt(StoreAdaptSEXP);
    Rcpp::traits::input_parameter< int >::type EndAdapt(EndAdaptSEXP);
    Rcpp::traits::input_parameter< int >::type PrintProgress(PrintProgressSEXP);
    Rcpp::traits::input_parameter< double >::type s2_delta(s2_deltaSEXP);
    Rcpp::traits::input_parameter< double >::type prior_delta(prior_deltaSEXP);
    __result = Rcpp::wrap(HiddenBASiCS_MCMCcppBatch(N, thin, burn, Counts, BatchDesign, mu0, delta0, phi0, s0, nu0, theta0, s2mu, adelta, bdelta, p_Phi, as, bs, atheta, btheta, ar, LSmu0, LSdelta0, LSphi0, LSnu0, LStheta0, sumByCellAll, sumByCellBio, sumByGeneAll, sumByGeneBio, StoreAdapt, EndAdapt, PrintProgress, s2_delta, prior_delta));
    return __result;
END_RCPP
}
// HiddenBASiCS_MCMCcppNoSpikes
Rcpp::List HiddenBASiCS_MCMCcppNoSpikes(int N, int thin, int burn, NumericMatrix Counts, NumericMatrix BatchDesign, NumericVector mu0, NumericVector delta0, NumericVector phi0, NumericVector nu0, double theta0, double s2mu, double adelta, double bdelta, NumericVector p_Phi, double aphi, double bphi, double atheta, double btheta, double ar, NumericVector LSmu0, NumericVector LSdelta0, NumericVector LSphi0, NumericVector LSnu0, double LStheta0, NumericVector sumByCellAll, NumericVector sumByGeneAll, int StoreAdapt, int EndAdapt, int PrintProgress, double s2_delta, double prior_delta, NumericVector BatchInfo, NumericVector BatchIds, NumericVector BatchSizes, NumericVector BatchOffSet, double Constrain, NumericMatrix InvCovMu, NumericVector Index, int ref, int ConstrainType, NumericVector ExpGene, NumericVector NotExpGene);
RcppExport SEXP BASiCS_HiddenBASiCS_MCMCcppNoSpikes(SEXP NSEXP, SEXP thinSEXP, SEXP burnSEXP, SEXP CountsSEXP, SEXP BatchDesignSEXP, SEXP mu0SEXP, SEXP delta0SEXP, SEXP phi0SEXP, SEXP nu0SEXP, SEXP theta0SEXP, SEXP s2muSEXP, SEXP adeltaSEXP, SEXP bdeltaSEXP, SEXP p_PhiSEXP, SEXP aphiSEXP, SEXP bphiSEXP, SEXP athetaSEXP, SEXP bthetaSEXP, SEXP arSEXP, SEXP LSmu0SEXP, SEXP LSdelta0SEXP, SEXP LSphi0SEXP, SEXP LSnu0SEXP, SEXP LStheta0SEXP, SEXP sumByCellAllSEXP, SEXP sumByGeneAllSEXP, SEXP StoreAdaptSEXP, SEXP EndAdaptSEXP, SEXP PrintProgressSEXP, SEXP s2_deltaSEXP, SEXP prior_deltaSEXP, SEXP BatchInfoSEXP, SEXP BatchIdsSEXP, SEXP BatchSizesSEXP, SEXP BatchOffSetSEXP, SEXP ConstrainSEXP, SEXP InvCovMuSEXP, SEXP IndexSEXP, SEXP refSEXP, SEXP ConstrainTypeSEXP, SEXP ExpGeneSEXP, SEXP NotExpGeneSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< int >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Counts(CountsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type BatchDesign(BatchDesignSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu0(mu0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type delta0(delta0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi0(phi0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nu0(nu0SEXP);
    Rcpp::traits::input_parameter< double >::type theta0(theta0SEXP);
    Rcpp::traits::input_parameter< double >::type s2mu(s2muSEXP);
    Rcpp::traits::input_parameter< double >::type adelta(adeltaSEXP);
    Rcpp::traits::input_parameter< double >::type bdelta(bdeltaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p_Phi(p_PhiSEXP);
    Rcpp::traits::input_parameter< double >::type aphi(aphiSEXP);
    Rcpp::traits::input_parameter< double >::type bphi(bphiSEXP);
    Rcpp::traits::input_parameter< double >::type atheta(athetaSEXP);
    Rcpp::traits::input_parameter< double >::type btheta(bthetaSEXP);
    Rcpp::traits::input_parameter< double >::type ar(arSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type LSmu0(LSmu0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type LSdelta0(LSdelta0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type LSphi0(LSphi0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type LSnu0(LSnu0SEXP);
    Rcpp::traits::input_parameter< double >::type LStheta0(LStheta0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sumByCellAll(sumByCellAllSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sumByGeneAll(sumByGeneAllSEXP);
    Rcpp::traits::input_parameter< int >::type StoreAdapt(StoreAdaptSEXP);
    Rcpp::traits::input_parameter< int >::type EndAdapt(EndAdaptSEXP);
    Rcpp::traits::input_parameter< int >::type PrintProgress(PrintProgressSEXP);
    Rcpp::traits::input_parameter< double >::type s2_delta(s2_deltaSEXP);
    Rcpp::traits::input_parameter< double >::type prior_delta(prior_deltaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type BatchInfo(BatchInfoSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type BatchIds(BatchIdsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type BatchSizes(BatchSizesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type BatchOffSet(BatchOffSetSEXP);
    Rcpp::traits::input_parameter< double >::type Constrain(ConstrainSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type InvCovMu(InvCovMuSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Index(IndexSEXP);
    Rcpp::traits::input_parameter< int >::type ref(refSEXP);
    Rcpp::traits::input_parameter< int >::type ConstrainType(ConstrainTypeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ExpGene(ExpGeneSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type NotExpGene(NotExpGeneSEXP);
    __result = Rcpp::wrap(HiddenBASiCS_MCMCcppNoSpikes(N, thin, burn, Counts, BatchDesign, mu0, delta0, phi0, nu0, theta0, s2mu, adelta, bdelta, p_Phi, aphi, bphi, atheta, btheta, ar, LSmu0, LSdelta0, LSphi0, LSnu0, LStheta0, sumByCellAll, sumByGeneAll, StoreAdapt, EndAdapt, PrintProgress, s2_delta, prior_delta, BatchInfo, BatchIds, BatchSizes, BatchOffSet, Constrain, InvCovMu, Index, ref, ConstrainType, ExpGene, NotExpGene));
    return __result;
END_RCPP
}
