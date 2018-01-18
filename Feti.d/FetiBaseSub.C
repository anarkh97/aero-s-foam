//
// Created by Michel Lesoinne on 1/18/18.
//

#include "FetiSub.h"

void
FetiBaseSub::sendMatProps(FSCommPattern<double> *matPat)
{
	for(int i = 0; i < scomm->numNeighb; ++i) {
		FSSubRecInfo<double> sInfo = matPat->getSendBuffer(subNum(), scomm->subNums[i]);
		sInfo.data[0] = Ymod;
		sInfo.data[1] = Prat;
		sInfo.data[2] = Dens;
		sInfo.data[3] = Thih;
		sInfo.data[4] = Sspe;
	}
}

void
FetiBaseSub::collectMatProps(FSCommPattern<double> *matPat)
{
	if(!neighbYmod) neighbYmod = new double[scomm->numNeighb];
	if(!neighbPrat) neighbPrat = new double[scomm->numNeighb];
	if(!neighbDens) neighbDens = new double[scomm->numNeighb];
	if(!neighbThih) neighbThih = new double[scomm->numNeighb];
	if(!neighbSspe) neighbSspe = new double[scomm->numNeighb];
	for(int iSub = 0; iSub < scomm->numNeighb; ++iSub) {
		FSSubRecInfo<double> rInfo = matPat->recData(scomm->subNums[iSub], subNum());
		neighbYmod[iSub] = (Ymod + rInfo.data[0])/2.0;
		neighbPrat[iSub] = (Prat + rInfo.data[1])/2.0;
		neighbDens[iSub] = (Dens + rInfo.data[2])/2.0;
		neighbThih[iSub] = (Thih + rInfo.data[3])/2.0;
		neighbSspe[iSub] = (Sspe + rInfo.data[4])/2.0;
	}
}


void
FetiBaseSub::sendWaveNumbers(FSCommPattern<double> *kPat)
{
  for(int i = 0; i < scomm->numNeighb; ++i) {
    FSSubRecInfo<double> sInfo = kPat->getSendBuffer(subNum(), scomm->subNums[i]);
    sInfo.data[0] = k_p;
    sInfo.data[1] = k_s;
    sInfo.data[2] = k_s2;
    sInfo.data[3] = k_f;
  }
}

void
FetiBaseSub::collectWaveNumbers(FSCommPattern<double> *kPat)
{
  if(!neighbK_p) neighbK_p = new double[scomm->numNeighb];
  if(!neighbK_s) neighbK_s = new double[scomm->numNeighb];
  if(!neighbK_s2) neighbK_s2 = new double[scomm->numNeighb];
  if(!neighbK_f) neighbK_f = new double[scomm->numNeighb];
  for(int i = 0; i < scomm->numNeighb; ++i) {
    FSSubRecInfo<double> rInfo = kPat->recData(scomm->subNums[i], subNum());
    neighbK_p[i] = (k_p + rInfo.data[0])/2.0;
    neighbK_s[i] = (k_s + rInfo.data[1])/2.0;
    neighbK_s2[i] = (k_s2 + rInfo.data[2])/2.0;
    neighbK_f[i] = (k_f + rInfo.data[3])/2.0;
  }
}
