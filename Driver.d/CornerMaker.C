#include <Utils.d/dbg_alloca.h>
#include <Driver.d/CornerMaker.h>
#include <Driver.d/SubDomain.h>
#include <Utils.d/SolverInfo.h>

#if defined(WINDOWS) || defined(MACOSX)
#include <cfloat>
#else
#include <climits>
#endif
#include <float.h>
#include <algorithm>

//#define DEBUG_CORNER

extern SolverInfo &solInfo;
namespace {
void crossprod(double [3], double [3], double [3]);

double magnitude(double[3]);

double magnitude(double v[3]) {
	return (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

// variable used to select corner nodes for the 2D-case
DofSet XYDofs = DofSet::Xdisp | DofSet::Ydisp;
}

static double checkArea(int np, int *pt, double pxyz[][3], int dim =3)
{
	double area = 0.0;
	if(dim < 3) {
		for(int i = 1; i < np; ++i) {
			double dx1 = pxyz[pt[i-1]][0] - pxyz[pt[0]][0];
			double dy1 = pxyz[pt[i-1]][1] - pxyz[pt[0]][1];
			double dz1 = pxyz[pt[i-1]][2] - pxyz[pt[0]][2];
			area += fabs(sqrt(dx1*dx1+dy1*dy1+dz1*dz1));
		}
		return area;
	}
	for(int i = 2; i < np; ++i)
	{
		double dx1 = pxyz[pt[i-1]][0] - pxyz[pt[0]][0];
		double dy1 = pxyz[pt[i-1]][1] - pxyz[pt[0]][1];
		double dz1 = pxyz[pt[i-1]][2] - pxyz[pt[0]][2];
		double dx2 = pxyz[pt[i]][0] - pxyz[pt[0]][0];
		double dy2 = pxyz[pt[i]][1] - pxyz[pt[0]][1];
		double dz2 = pxyz[pt[i]][2] - pxyz[pt[0]][2];
		area += fabs(sqrt(
				(dy1*dz2-dy2*dz1)*(dy1*dz2-dy2*dz1) +
				(dz1*dx2-dz2*dx1)*(dz1*dx2-dz2*dx1) +
				(dx1*dy2-dx2*dy1)*(dx1*dy2-dx2*dy1)));
	}
	return area;
}

CornerMaker::CornerMaker(int _glNumSub, int _nSub, FetiSubCornerHandler **_cornerHandler,
                         FSCommPattern<int> *_cpat, FSCommunicator *_communicator)
{
	glNumSub = _glNumSub;
	nSub = _nSub;
	communicator = _communicator;
	cornerHandler = _cornerHandler;
	cpat = _cpat;

	for(int i=0; i<4; ++i) dims[i] = 0;
	paralApply(nSub, cornerHandler, &FetiSubCornerHandler::markDims, dims);
#ifdef DISTRIBUTED
	communicator->globalMax(4, dims);
#endif
	int sdim = dims[0]+dims[1]+dims[2]; // structure
	int fdim = dims[3]; // fluid
	bool mixed = (fdim && sdim);
	dim = (mixed || (fdim==0)) ? sdim : fdim; // use the structure dimension for a mixed subdomain XXXX
#ifdef DEBUG_CORNER
	cerr << "global dim = " << dim << endl;
#endif
}

CornerMaker::~CornerMaker()
{
	if(cornerHandler) {
		for(int i=0; i<nSub; ++i)
			if(cornerHandler[i]) { delete cornerHandler[i]; cornerHandler[i] = 0; }
		delete [] cornerHandler; cornerHandler = 0;
	}
	delete [] glSubGroup;
	// don't delete grToSub (GenDecDomain)
}

int
CornerMaker::makeCorners()
{
	int i, iSub;

	// First locate ``unsafe'' nodes. These nodes must be corner
	// nodes (to eliminate subdomain ZEMs) but do not work in tying two subdomains together
	paralApply(nSub, cornerHandler, &FetiSubCornerHandler::dispatchSafeNodes, cpat);
	cpat->exchange();
	paralApply(nSub, cornerHandler, &FetiSubCornerHandler::markSafeNodes, cpat);

	// Take care of the multi degreed nodes that could upset the coarse problem
	paralApply(nSub, cornerHandler, &FetiSubCornerHandler::markMultiDegNodes);

	// Take care of the corners for 4th order problems
	paralApply(nSub, cornerHandler, &FetiSubCornerHandler::dispatchRotCorners, cpat);
	cpat->exchange();
	paralApply(nSub, cornerHandler, &FetiSubCornerHandler::markRotCorners, cpat);

	// Build the list of potential corner and mark those already certainly chosen
	// Phase 1, count how many corners each sub has and to how many subdomains these corners connect
	int *cPerSub = new int[glNumSub+1];
	int *nTot = new int[glNumSub];
	for(i = 0; i < glNumSub; ++i) cPerSub[i] = nTot[i] = 0;

	// countAndMarkCornerCand internaly marks corner candidates and returns the numbers of corners
	// for which this subdomain is master and the total number of corners connectivities
	// weight[inode] contains the local numbering of corner candidates for which this sub is a master
	paralApply(nSub, cornerHandler, &FetiSubCornerHandler::countAndMarkCornerCand, cPerSub, nTot);
	communicator->globalSum(glNumSub, cPerSub);
	communicator->globalSum(glNumSub, nTot);

	int totNC, tmpNC;
	totNC = 0;
	// start creating the corner candidate to sub (master) connectivity
	for(i = 0; i < glNumSub; ++i) {
		tmpNC = cPerSub[i];
		cPerSub[i] = totNC;
		totNC += tmpNC;
	}
	cPerSub[glNumSub] = totNC;

#ifdef DEBUG_CORNER
	cerr << "Number of corner candidates = " << totNC << endl;
#endif

	int *cPtr = new int[totNC+1];
	int tot = 0;
	for(i = 0; i < glNumSub; ++i) {
		for(int j = cPerSub[i]; j < cPerSub[i+1]; ++j) cPtr[j] = tot;
		tot += nTot[i];
	}
	cPtr[cPerSub[glNumSub]] = tot;

	int *cTg = new int[tot];
	double (*xyz)[3] = new double[tot][3];
	int *essentPtr = new int[glNumSub+1];
	for(i = 0; i < tot; ++i) {
		cTg[i] = 0;
		xyz[i][0] = xyz[i][1] = xyz[i][2] = 0.0;
	}
	for(i = 0; i < glNumSub; ++i) essentPtr[i] = 0;

	char *essential = new char[totNC];
	for(i = 0; i < totNC; ++i) essential[i] = 0;

	// getCornerXYZ collects the coordinates of the corner candidates and in addition lets the
	// caller know which ones are already essential corners and finally fills the corner to sub
	// connectivity data and count the number of corner with rot dofs
	paralApply(nSub, cornerHandler, &FetiSubCornerHandler::getCornerXYZ,
	           cPerSub, xyz, essential, cPtr, cTg);
	communicator->globalMax(totNC, cPtr);
	communicator->globalSum(3*tot, (double *) xyz);
	communicator->globalSum(tot, cTg);
#if defined(LAM_MPI) || defined(SUN10) || defined(OPEN_MPI)
	// temporary fix, LAM 7.1.1 MPI_Allreduce isn't working with character so converting to int
  int *essential_tmp = new int[totNC];
  for(i = 0; i < totNC; ++i) essential_tmp[i] = int(essential[i]);
  communicator->globalSum(totNC, essential_tmp);
  for(i = 0; i < totNC; ++i) essential[i] = char(essential_tmp[i]);
  delete [] essential_tmp;
#else
	communicator->globalSum(totNC, essential);
#endif
	// now we have the "corner candidates" to subdomain connectivity
	Connectivity *cToSub = new Connectivity(totNC, cPtr, cTg);

	// exchange initial global numbering of the corner candidates (because we need
	// it to build the "rotational corner" connectivity
	paralApply(nSub, cornerHandler, &FetiSubCornerHandler::dispatchInitialNumbering, cpat, cPerSub);
	cpat->exchange();
	paralApply(nSub, cornerHandler, &FetiSubCornerHandler::recInitialNumbering, cpat, essentPtr);
	communicator->globalSum(glNumSub, essentPtr);

	int cNum = 0;
	for(iSub = 0; iSub < glNumSub; ++iSub) {
		int tmp = cNum;
		cNum += essentPtr[iSub];
		essentPtr[iSub] = tmp;
	}
	essentPtr[glNumSub] = cNum;

	int *essentTg = new int[cNum];
	for(i = 0; i < cNum; ++i) essentTg[i] = 0;

	paralApply(nSub, cornerHandler,&FetiSubCornerHandler::listRotCorners, essentPtr, essentTg);
	communicator->globalSum(cNum, essentTg);

	Connectivity subToRotCrn(glNumSub, essentPtr, essentTg);

	// initialize glCrnGroup and glSubGroup
	int *glCrnGroup = new int[totNC];
	for(i = 0; i < totNC; ++i) glCrnGroup[i] = -1;
	glSubGroup = new int[glNumSub];
	for(i = 0; i < glNumSub; ++i) glSubGroup[i] = -1;

	// Now chose the corners
	chooseCorners(essential, xyz, *cToSub, subToRotCrn, glCrnGroup);

#ifdef DEBUG_CORNER
	int xx = 0;
  for(i = 0; i < totNC; ++i) if(essential[i]) xx++;
  fprintf(stderr, "Found %d essential safe corners\n", xx);
#endif

	// Count contact corners
	int *numCnt = new int[glNumSub];
	for(iSub = 0; iSub < glNumSub; ++iSub)
		numCnt[iSub] = 0;

	paralApply(nSub, cornerHandler,&FetiSubCornerHandler::countContact,
	           numCnt, essential);
	communicator->globalSum(glNumSub, numCnt);

	int *newCPerSub = new int[glNumSub];
	int cCount = 0;
	for(iSub = 0; iSub < glNumSub; ++iSub) {
		newCPerSub[iSub] = cCount;
		for(int jj = cPerSub[iSub]; jj < cPerSub[iSub+1]; ++jj)
			if(essential[jj]) cCount++;
	}

	int cntCount = cCount;
	for(iSub = 0; iSub < glNumSub; ++iSub) {
		int tmpCntCount = numCnt[iSub];
		numCnt[iSub] = cntCount;
		cntCount += tmpCntCount;
	}
#ifdef DEBUG_CORNER
	fprintf(stderr, "Found %d essential safe+unsafe corners\n", cntCount);
#endif
	paralApply(nSub, cornerHandler, &FetiSubCornerHandler::dispatchNumbering,
	           cpat, essential, cPerSub, newCPerSub, cCount, numCnt);
	cpat->exchange();

	paralApply(nSub, cornerHandler, &FetiSubCornerHandler::recNumbering,
	           cpat, newCPerSub);

	paralApply(nSub, cornerHandler, &FetiSubCornerHandler::resendNumbers, cpat);
	cpat->exchange();
	paralApply(nSub, cornerHandler, &FetiSubCornerHandler::checkNumbers, cpat);

	delete [] essential;
	delete [] numCnt;
	delete [] newCPerSub;
	delete [] nTot;
	delete [] cPerSub;
	delete [] glCrnGroup;
	delete cToSub;

	// return the total number of corners including the ones for contact
	return cntCount;
}

void
CornerMaker::chooseCorners(char *glCornerList, double (*crnXYZ)[3],
                           Connectivity &cNConnect, Connectivity &subToRotCrn,
                           int *glCrnGroup)
{
	Connectivity *subToCrn = cNConnect.reverse();
	Connectivity *rotCrnToSub = subToRotCrn.reverse();
	// fprintf(stderr, "We have %d rot corners\n", rotCrnToSub->csize());
	Connectivity *rotSubToSub = subToRotCrn.transcon(rotCrnToSub);
	grToSub = rotSubToSub->collapse();

	double aamin = DBL_MAX, aamax = -DBL_MAX;

	// guarantee we will loop at least once
	int lastGrSize = grToSub->csize()+1;
	while(grToSub->csize() > 1 && grToSub->csize() != lastGrSize) {
		lastGrSize = grToSub->csize();
		int iGr;
		// fprintf(stderr, "Group has size %d\n", grToSub->csize());
		Connectivity *grToCrn = grToSub->transcon(subToCrn);
		Connectivity *crnToGr = grToCrn->reverse();
		// Eliminate the corners that are not corners anymore
		//grToCrn = grToCrn->trim(crnToGr);
		grToCrn->sortTargets();
		Connectivity *grToGr = grToCrn->transcon(crnToGr);
		// Now examine each group and establish a set of connection
		// priorities
		int (*fav)[2] = new int[grToSub->csize()][2];
		double (*bamax)[2] = new double[grToSub->csize()][2];
		int (*choices)[2][3] = new int[grToSub->csize()][2][3];
		int *finalFav = new int[grToSub->csize()];
		bool *tied = new bool[grToSub->csize()];

		for(iGr = 0; iGr < grToCrn->csize(); ++iGr) {
			tied[iGr] = false;
			int jGr;
			finalFav[iGr] = fav[iGr][0] = fav[iGr][1]  = -1;
			double tamax = 0, grAmax = 0.0;
			int ntmax = -1;

			//double amax = 0, amax2 = 0;
			bamax[iGr][0] = bamax[iGr][1] = 0.0;
			for(jGr = 0; jGr < grToGr->num(iGr); ++jGr) {
				int grJ = (*grToGr)[iGr][jGr];
				if(grJ == iGr) continue;
				// Find the comon nodes
				int maxCm = std::min( grToCrn->num(iGr), grToCrn->num(grJ) );
				// XML this needs to be changed
				int *cnode = (int *) dbg_alloca(maxCm*sizeof(int));
				int nc = 0;
				int n1 = 0, n2 = 0;
				while(n1 < grToCrn->num(iGr) && n2 < grToCrn->num(grJ)) {
					if( (*grToCrn)[iGr][n1] == (*grToCrn)[grJ][n2] ) {
						cnode[nc++] = (*grToCrn)[iGr][n1];
						n1++; n2++;
					} else {
						if( (*grToCrn)[iGr][n1] < (*grToCrn)[grJ][n2] )
							n1++;
						else
							n2++;
					}
				}
				if(nc < dim) continue;
				int iC;
				// pick the candidate with most number of groups connected to it
				// Or a candidate already chosen
				int bc = cnode[0];
				for(iC = 1; iC < nc; ++iC)
					if(glCornerList[cnode[iC]]) {
						if(glCornerList[bc] == 0 ||
						   crnToGr->num(cnode[iC]) > crnToGr->num(bc))
							bc = cnode[iC];
					} else if(glCornerList[bc] == 0
					          && crnToGr->num(cnode[iC]) > crnToGr->num(bc))
						bc = cnode[iC];

				// Now identify the furthest nodes from one another
				double mxDist2 = 0.0, mxDist2b = 0.0;
				int c2 = -1;
				int c2b = -1;
				for(iC = 0; iC < nc; ++iC) {
					int ndI = cnode[iC];
					int ndB = bc;
					double dst2 =
							(crnXYZ[ndB][0]-crnXYZ[ndI][0])*
							(crnXYZ[ndB][0]-crnXYZ[ndI][0]) +
							(crnXYZ[ndB][1]-crnXYZ[ndI][1])*
							(crnXYZ[ndB][1]-crnXYZ[ndI][1]) +
							(crnXYZ[ndB][2]-crnXYZ[ndI][2])*
							(crnXYZ[ndB][2]-crnXYZ[ndI][2]);
					if(dst2 > mxDist2) {
						mxDist2 = dst2;
						c2 = cnode[iC];
					}
					if(glCornerList[ndI] && dst2 > mxDist2b) {
						mxDist2b = dst2;
						c2b = ndI;
					}
				}
				// XML put it back
				if(mxDist2b >= 1.6 * mxDist2)
					c2 = c2b;
				if(c2 < 0) {
					if(nc > 1)
						fprintf(stderr, "no line %d %e %e, %e %e\n", nc,
						        crnXYZ[cnode[0]][0],
						        crnXYZ[cnode[1]][0], crnXYZ[cnode[0]][1],
						        crnXYZ[cnode[1]][1]);
					else
						fprintf(stderr,"no line %d\n",nc);
					continue; }
				// finally, find the node that maximizes the area with nodes 1 & 2
				int c3 = -1;
				if(dim >= 3) {
					double maxArea=0.0, area;
					double dx = crnXYZ[bc][0] - crnXYZ[c2][0];
					double dy = crnXYZ[bc][1] - crnXYZ[c2][1];
					double dz = crnXYZ[bc][2] - crnXYZ[c2][2];
					for(iC = 0; iC < nc; ++iC) {
						int ndI = cnode[iC];
						int ndC = c2;
						double dx2= crnXYZ[ndC][0] - crnXYZ[ndI][0];
						double dy2= crnXYZ[ndC][1] - crnXYZ[ndI][1];
						double dz2= crnXYZ[ndC][2] - crnXYZ[ndI][2];
						double cross[3] = {dy*dz2 - dz*dy2, dx2*dz-dz2*dx, dx*dy2 - dx2*dy};
						area = (cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]);
						if(area > maxArea) { maxArea = area; c3 = cnode[iC]; }
						//if(area > maxAreab) { maxAreab = area; c3b = cnode[iC]; }
					}

					if(maxArea > bamax[iGr][0]) {
						fav[iGr][1] = fav[iGr][0];
						bamax[iGr][1] = bamax[iGr][0];
						choices[iGr][1][0] = choices[iGr][0][0];
						choices[iGr][1][1] = choices[iGr][0][1];
						choices[iGr][1][2] = choices[iGr][0][2];
						fav[iGr][0] = grJ;
						choices[iGr][0][0] = bc;
						choices[iGr][0][1] = c2;
						choices[iGr][0][2] = c3;
						bamax[iGr][0] = maxArea;
					} else if(maxArea > bamax[iGr][1]) {
						fav[iGr][1] = grJ;
						choices[iGr][1][0] = bc;
						choices[iGr][1][1] = c2;
						choices[iGr][1][2] = c3;
						bamax[iGr][1] = maxArea;
					}
				} else {
					if(mxDist2 > bamax[iGr][0]) {
						fav[iGr][1] = fav[iGr][0];
						bamax[iGr][1] = bamax[iGr][0];
						choices[iGr][1][0] = choices[iGr][0][0];
						choices[iGr][1][1] = choices[iGr][0][1];
						choices[iGr][1][2] = choices[iGr][0][2];
						bamax[iGr][0] = mxDist2;
						fav[iGr][0] = grJ;
						choices[iGr][0][0] = bc;
						choices[iGr][0][1] = c2;
						choices[iGr][0][2] = -1;
					}
				}

				// Now double check if we already had the domain tied
				int ncc = 0;
				for(iC = 0; iC < nc; ++iC)
					if(glCornerList[cnode[iC]])
						cnode[ncc++] = cnode[iC];
				if(ncc >= 3) {
					double area = checkArea(ncc, cnode, crnXYZ);
					area *= area;
					//fprintf(stderr, "Area pretied is %e vs %e %e\n", area, amax, maxArea);
					if(area > tamax) {
						ntmax = grJ;
						tamax = area;
					}
					// Keeping track of the largest area for this group
				}
				if(bamax[iGr][0] > grAmax)
					grAmax = bamax[iGr][0];
			}
			if(tamax > 0.4*grAmax) {
				//fprintf(stderr, "Doing a tie\n");
				finalFav[iGr] = ntmax;
				tied[iGr] = true;
				tied[ntmax] = true;
			} //else
			//fprintf(stderr, "No tie %e %e\n", tamax, grAmax);
			if(grAmax > aamax)
				aamax = grAmax;
			if(0.6*grAmax < aamin)
				aamin = 0.6*grAmax;
		}
		// Make pairs of groups
		for(iGr = 0; iGr < grToCrn->csize(); ++iGr) {
			if(tied[iGr]) continue;
			if(fav[iGr][0] < 0) {
				//fprintf(stderr, "Warning: Group %d cannot be tied\n", iGr);
				int in;
				for(in = 0; in < grToSub->num(iGr); ++in) {
//           fprintf(stderr, "Sub %d\n", (*grToSub)[iGr][in]);
					glSubGroup[(*grToSub)[iGr][in]] = iGr;
				}
				for(in = 0; in < grToCrn->num(iGr); ++in) {
//           fprintf(stderr, "Corner %d\n", (*grToCrn)[iGr][in]);
					glCrnGroup[(*grToCrn)[iGr][in]] = iGr;
				}
			} else {
				// check if one of our favorites is not tied yet
				int pref = -1;
				if(tied[ fav[iGr][0] ] == false)
					pref = 0;
				else if(fav[iGr][1] >= 0 && tied[ fav[iGr][1] ] == false) {
					if(bamax[iGr][1] > 0.25*bamax[iGr][0])
						pref = 1;
					else
						pref = 0;
				}
				if(pref >= 0) {
					glCornerList[choices[iGr][pref][0]] = 1;
					if(choices[iGr][pref][1] > -1) glCornerList[choices[iGr][pref][1]] = 1; // PJSA
					if(choices[iGr][pref][2] > -1) glCornerList[choices[iGr][pref][2]] = 1; // PJSA 1-16-07
					tied[ fav[iGr][pref] ] = true;
					tied[ iGr ] = true;
					finalFav[iGr] = fav[iGr][pref];
				}
			}
		}
		// now attach the groups that had not been attached yet
		for(iGr = 0; iGr < grToCrn->csize(); ++iGr) {
			if(tied[iGr]) continue;
			if(fav[iGr][0] >= 0) {
				glCornerList[choices[iGr][0][0]] = 1;
				if(choices[iGr][0][1] > -1) glCornerList[choices[iGr][0][1]] = 1;
				if(choices[iGr][0][2] > -1) glCornerList[choices[iGr][0][2]] = 1; // PJSA 1-16-07
				finalFav[iGr] = fav[iGr][0];
			}
		}
		// Now create a new group to group connectivity
		int *count = new int[grToCrn->csize()+1];
		for(iGr = 0; iGr < grToCrn->csize(); ++iGr)
			count[iGr] = 1;
		count[grToCrn->csize()] = 0;
		for(iGr = 0; iGr < grToCrn->csize(); ++iGr)
			if(finalFav[iGr] >= 0)
				count[finalFav[iGr]]++;
		for(iGr = 0; iGr < grToCrn->csize(); ++iGr)
			count[iGr+1] += count[iGr];
		int *target = new int[count[grToCrn->csize()] ];
		for(iGr = 0; iGr < grToCrn->csize(); ++iGr) {
			target[--count[iGr]] = finalFav[iGr];
			if(finalFav[iGr] >= 0)
				target[--count[ finalFav[iGr] ] ] = iGr;
		}
		Connectivity *gr2Gr = new Connectivity(grToCrn->csize(), count, target);

		// collapse it
		Connectivity *grToPrevGr = gr2Gr->collapse();
		// get the group To Sub graph
		Connectivity *oldGrToSub = grToSub;
		grToSub = grToPrevGr->transcon(grToSub);
		delete [] finalFav;
		delete [] tied;
		delete [] choices;
		delete [] fav;
		delete [] bamax;
		delete oldGrToSub;
		delete grToPrevGr;
		delete gr2Gr;

		delete grToCrn;
		delete crnToGr;
		delete grToGr;
	}
	//communicator->sync();
	//fprintf(stderr, "Possible extremes: %e %e\n",aamax, aamin);
	delete [] crnXYZ;

	delete subToCrn;
	delete rotCrnToSub;
	delete rotSubToSub;
}

