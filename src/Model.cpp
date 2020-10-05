/************************************************************************
 * SnipSnip, version 1.1
 * Copyright 2013
 * Richard Howey
 * Institute of Genetic Medicine, Newcastle University
 *
 * richard.howey@ncl.ac.uk
 * http://www.staff.ncl.ac.uk/richard.howey/
 *
 * This file is part of SnipSnip, the pedigree file processing program for EMIM.
 *
 * SnipSnip is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SnipSnip is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SnipSnip.  If not, see <http://www.gnu.org/licenses/>.
 ************************************************************************/


/*! \file Model.cpp
    \brief This file contains the source of models used for logistic regression.
    
*/

#include <map>
#include <iostream>
#include <math.h>

using namespace std; // initiates the "std" or "standard" namespace

#include "Model.h"
#include "Data.h"
#include "ModelFitting.h"

//! Sets the values of all of the parameters.
void Model::setNewParameters(map<unsigned int, double> & paras)
{
	for(map<unsigned int, double>::const_iterator p = paras.begin(); p != paras.end(); ++p)
	{
		setParameter(p->first, p->second);
	};
};

//! Returns a value from the cache.
pair<double, double> Model1SNPLogReg::getNegLogBeta0(const unsigned int & noCases, const unsigned int & noControls)
{
	//set number of cases and controls without missing data
	numCases = (double)noCases;
	numControls = (double)noControls;

	//not applicable when covariates are present
	map<pair<unsigned int, unsigned int>, pair<double, double> >::const_iterator i = cache.find(make_pair(noCases, noControls));
	if(i != cache.end())
	{		
			return i->second;		
	};

	double beta0 = log(numCases/numControls);
	double negLog = (numCases + numControls)*log(1 + numCases/numControls) - beta0*(numCases);	

	cache[make_pair(noCases, noControls)] = make_pair(negLog, beta0);

	return make_pair(negLog, beta0);
};

void Model1SNPLogReg::setGenotypeCounts(unsigned int & ca0, unsigned int & ca1, unsigned int & ca2, unsigned int & co0, unsigned int & co1, unsigned int & co2)
{
		genotypeCountsCases[0] = ca0;
		genotypeCountsCases[1] = ca1;
		genotypeCountsCases[2] = ca2;
		genotypeCountsControls[0] = co0;
		genotypeCountsControls[1] = co1;
		genotypeCountsControls[2] = co2;
};

//! Sets parameters.
void Model1SNPLogReg::setParameter(const unsigned int & no, const double & value)
{
	parameters[no] = value;

	switch(no)
	{
		case 1:
			beta0 = value; return;
		case 2:
			beta1 = value; return;
	};

};

//! Returns negative log likelihood for single-SNP log regression.
double Model1SNPLogReg::negLogLikelihood()
{
	if(useCovariates) return negLogLikelihoodCovar();
	double ans = 0;
	double effects, noCases, noControls;

	//loop thro' the number of minor alleles for SNP1, i.e. start SNP
	for(unsigned int snp1 = 0; snp1 <= 2; ++snp1)
	{					
		effects = beta0 + snp1*beta1;
		noCases = (double)genotypeCountsCases[snp1];
		noControls = (double)genotypeCountsControls[snp1];
		
		ans += -noCases*effects + (noCases+noControls)*log(1 + exp(effects));
	};
	
	return ans;
};

//! Get gradient vector for 2 parameter model.
void Model1SNPLogReg::getGradientVector(map<unsigned int, double> & gradientVector, const bool & fitBeta1, const bool & fitBeta2)  const
{
	if(useCovariates) {getGradientVectorCovar(gradientVector, fitBeta1, fitBeta2); return;};
	//Always fitting beta0 and beta1 for 1 SNP model
	double ans[2] = {0, 0};

	//loop thro' the number of minor alleles for SNP1, i.e. start SNP
	for(unsigned int snp1 = 0; snp1 <= 2; ++snp1)
	{
				
		effects = beta0 + snp1*beta1;
		noCases = (double)genotypeCountsCases[snp1];
		noControls = (double)genotypeCountsControls[snp1];

		aNumber = -noCases + (noCases + noControls)/(exp(-effects) + 1);
			
		ans[0] += aNumber;		
		ans[1] += snp1*aNumber;					
	};

	gradientVector[1] = ans[0];	
	gradientVector[2] = ans[1];		
};

//! Returns the 2nd derivative w.r.t. chosen parameters of the negative log likelihood.
void Model1SNPLogReg::getHessianMatrix(map<unsigned int, map<unsigned int, double> > & hessianMatrix, const bool & fitBeta1, const bool & fitBeta2) const
{
	if(useCovariates) {getHessianMatrixCovar(hessianMatrix, fitBeta1, fitBeta2); return;};

	//Always fitting beta0 and beta1 for 1 SNP model
	double ans[2][2] = {{0, 0}, {0, 0}}; //col, row
	
	//loop thro' the number of minor alleles for SNP1, i.e. anchor SNP
	for(unsigned int snp1 = 0; snp1 <= 2; ++snp1)
	{		
		effects = beta0 + snp1*beta1;
		expEffects = exp(effects);
		noCases = (double)genotypeCountsCases[snp1];
		noControls = (double)genotypeCountsControls[snp1];
		aNumber = ((noCases+noControls)*expEffects)/((expEffects + 1)*(expEffects + 1));

		ans[0][0] += aNumber;		
		ans[0][1] += snp1*aNumber;	
		ans[1][1] += snp1*snp1*aNumber;				
	};

	//The Hessian matrix is symetric so only calculate half and then copy	
	ans[1][0] = ans[0][1];	
	
	//setup the matrix with calculated values
	map<unsigned int, double> aCol;	

	for(unsigned int col = 0; col < 2; ++col)
	{		
		for(unsigned int row = 0; row < 2; ++row)
		{
			aCol[row + 1] = ans[col][row];	
		};
		hessianMatrix[col + 1] = aCol;
	};
};


//! Returns negative log likelihood for single-SNP log regression.
double Model1SNPLogReg::negLogLikelihoodCovar()
{
	double ans = 0;
	
	//loop thro' each subject
	list<unsigned char>::const_iterator anc = anchorSNP->noMinorAllelesAllSubjects.begin();	
	list<list<double> >::const_iterator covar = covariateData->covariateDataAllSubjects.begin();
	list<unsigned char>::const_iterator cc = caseControls.begin();

	unsigned int covariateParameterNo;

	for( ; anc != anchorSNP->noMinorAllelesAllSubjects.end(); ++anc, ++covar, ++cc)
	{
		if(*anc != 3 && *cc != 0) //skip missing data
		{
			effects = beta0 + (*anc)*beta1;

			covariateParameterNo = 5;
			//add covariates effects
			for(list<double>::const_iterator cv = covar->begin(); cv != covar->end(); ++cv)
			{
				effects += (*cv)*parameters[covariateParameterNo];
				covariateParameterNo++;
			};

			if(*cc == 2) ans -= effects;
			ans += log(1 + exp(effects));
		};
	};

	return ans;
};

//! Get gradient vector for 2 parameter model.
void Model1SNPLogReg::getGradientVectorCovar(map<unsigned int, double> & gradientVector, const bool & fitBeta1, const bool & fitBeta2)  const
{
	//Fitting beta0 or beta0 & beta1 for 1 SNP model
	double ans[2] = {0, 0};
		
	//loop thro' each subject
	list<unsigned char>::const_iterator anc = anchorSNP->noMinorAllelesAllSubjects.begin();	
	list<list<double> >::const_iterator covar = covariateData->covariateDataAllSubjects.begin();
	list<unsigned char>::const_iterator cc = caseControls.begin();

	unsigned int covariateParameterNo;
	unsigned int cp = 2; //first parameter no. for covariate in ans vector
	if(fitBeta1) cp++;
	unsigned int cpi = cp;

	for(list<double>::const_iterator cv = covar->begin(); cv != covar->end(); ++cv, ++cpi) gradientVector[cpi] = 0;	

	expEffectsAllIndiv.clear();

	//add covariates parameter derivs
	for( ; anc != anchorSNP->noMinorAllelesAllSubjects.end(); ++anc, ++covar, ++cc)
	{
		if(*anc != 3 && *cc != 0) //skip missing data
		{
			effects = beta0 + (*anc)*beta1;

			covariateParameterNo = 5;
			//add covariates effects
			for(list<double>::const_iterator cv = covar->begin(); cv != covar->end(); ++cv)
			{
				effects += (*cv)*getParameter(covariateParameterNo);
				covariateParameterNo++;
			};

			//cache exp(effects) for use with the hessian matrix
			expEffects = exp(effects); 
			expEffectsAllIndiv.push_back(expEffects);

			if(*cc == 2) aNumber = -1 + expEffects/(expEffects + 1);
			else aNumber = expEffects/(expEffects + 1);		

			ans[0] += aNumber;			
			if(fitBeta1) ans[1] += (*anc)*aNumber;					

			cpi = cp;
			//add covariates parameter derivs
			for(list<double>::const_iterator cv = covar->begin(); cv != covar->end(); ++cv, ++cpi)
			{
				gradientVector[cpi] += (*cv)*aNumber;			
			};
		};
	};

	gradientVector[1] = ans[0];	
	if(fitBeta1) gradientVector[2] = ans[1];		
};

//! Returns the 2nd derivative w.r.t. chosen parameters of the negative log likelihood.
void Model1SNPLogReg::getHessianMatrixCovar(map<unsigned int, map<unsigned int, double> > & hessianMatrix, const bool & fitBeta1, const bool & fitBeta2) const
{
	//Always fitting beta0 and beta1 for 1 SNP model
	double ans[2][2] = {{0, 0}, {0, 0}}; //col, row
	map<pair<unsigned int, unsigned int>, double> covarAns; //starting count from 0, covariates start after other parameters, either 2 or 3

	//loop thro' each subject
	list<unsigned char>::const_iterator anc = anchorSNP->noMinorAllelesAllSubjects.begin();	
	list<list<double> >::const_iterator covar = covariateData->covariateDataAllSubjects.begin();
	list<unsigned char>::const_iterator cc = caseControls.begin();
	list<double>::const_iterator cv1, cv2;

	//unsigned int covariateParameterNo;
	unsigned int cp = 1; //first parameter no. for covariate in covarAns matrix, starts at 0
	if(fitBeta1) ++cp;

	double val1;

	unsigned int noParas = cp + covar->size();

	list<double>::const_iterator expEffs = expEffectsAllIndiv.begin();//set up in Vector

	//loop thro' individuals
	for( ; anc != anchorSNP->noMinorAllelesAllSubjects.end(); ++anc, ++covar, ++cc)
	{
		if(*anc != 3 && *cc != 0) //skip missing data
		{
			aNumber = (*expEffs)/((*expEffs + 1)*(*expEffs + 1));
			++expEffs;

			//only calculate for the parameters we are fitting over
			ans[0][0] += aNumber;
			
			if(fitBeta1)
			{
				ans[0][1] += (*anc)*aNumber;	
				ans[1][1] += (*anc)*(*anc)*aNumber;		
			};

			cv1 = covar->begin();

			//add covariate 2nd derivs
			for(unsigned int i1 = 0; i1 < noParas; )
			{
				cv2 = covar->begin();

				for(unsigned int i2 = 0; i2 < noParas; )
				{					
					//update any bit of matrix that is part of covariate in top half of matrix 
					if(i2 >= cp && i2 >= i1)
					{
						if(i1 < cp) //covariate and snp parameter 2nd deriv
						{
							if(i1 == 0) val1 = 1;
							else if(i1 == 1)
							{
								val1 = (*anc);								
							};							

							covarAns[make_pair(i1, i2)] += val1*(*cv2)*aNumber;
						}
						else //covariate and covariate 2nd deriv
						{
							covarAns[make_pair(i1, i2)] += (*cv1)*(*cv2)*aNumber;	
						};
					};

					++i2;
					if(i2 > cp) ++cv2; //move to next covariate value
				};

				++i1;
				if(i1 > cp) ++cv1; //move to next covariate value
			};

		};//end of skipping missing data

	};//end of subject loop


	//The Hessian matrix is symetric so only calculate half and then copy	
	if(fitBeta1) ans[1][0] = ans[0][1];	

	//setup the matrix with calculated values
	map<unsigned int, double> aCol;

	//no parameters not including covariates is given by cp
	for(unsigned int col0 = 0; col0 < cp; ++col0)
	{		
		for(unsigned int row0 = 0; row0 < cp; ++row0)
		{
			aCol[row0 + 1] = ans[col0][row0];	
		};
		hessianMatrix[col0 + 1] = aCol;
	};

	//add calculated covariate 2nd derivs
	for(unsigned int col = 0; col < noParas; ++col)
	{	
		if(col < cp) //get predefined column
		{
			aCol = hessianMatrix[col + 1];
		};

		for(unsigned int row = 0; row < noParas; ++row)
		{
			if(col >= cp && col >= row)	//top right half - except for non-covariate
			{
				aCol[row + 1] = covarAns[make_pair(row, col)];
			}
			else if(row >= cp && col < row)	//bottom left half - except for non-covariate
			{
				aCol[row + 1] = covarAns[make_pair(col, row)];
			};
		};

		hessianMatrix[col + 1] = aCol;
	};

};

//! Returns the value of the negative log likelihood for the given parameters and joint genotype counts.
double TwoSNPLogRegModel::negLogLikelihood()
{
	if(useCovariates) return negLogLikelihoodCovar();
	double ans = 0;
	
	//loop thro' the number of minor alleles for SNP1, i.e. start SNP
	for(unsigned int snp1 = 0; snp1 <= 2; ++snp1)
	{
		//loop thro' the number of minor alleles for SNP2, i.e. partner SNP
		for(unsigned int snp2 = 0; snp2 <= 2; ++snp2)
		{			
			effects = beta0 + snp1*beta1 + snp2*beta2;
			noCases = (double)partnerJointGenotypeCountsCases->counts[snp1][snp2];
			noControls = (double)partnerJointGenotypeCountsControls->counts[snp1][snp2];
		
			ans += -noCases*effects + (noCases+noControls)*log(1 + exp(effects));
		};

	};
	
	return ans;
};


//! Returns the value of the negative log likelihood for the given parameters and joint genotype counts.
double TwoSNPLogRegModel::negLogLikelihoodCovar()
{
	double ans = 0;
	
	//loop thro' each subject
	list<unsigned char>::const_iterator anc = anchorSNP->noMinorAllelesAllSubjects.begin();
	list<unsigned char>::const_iterator part = partnerSNP->noMinorAllelesAllSubjects.begin();
	list<list<double> >::const_iterator covar = covariateData->covariateDataAllSubjects.begin();
	list<unsigned char>::const_iterator cc = caseControls.begin();

	unsigned int covariateParameterNo;

	for( ; anc != anchorSNP->noMinorAllelesAllSubjects.end(); ++anc, ++part, ++covar, ++cc)
	{
		if(*anc != 3 && *part != 3 && *cc != 0) //skip missing data
		{
			effects = beta0 + (*anc)*beta1 + (*part)*beta2;

			covariateParameterNo = 5;
			//add covariates effects
			for(list<double>::const_iterator cv = covar->begin(); cv != covar->end(); ++cv)
			{
				effects += (*cv)*parameters[covariateParameterNo];
				covariateParameterNo++;
			};

			if(*cc == 2) ans -= effects; //a case
			ans += log(1 + exp(effects));
		};
	};

	return ans;
};

//! Sets the value of a given parameter for the model.
void Model2SNPs::setParameter(const unsigned int & no, const double & value)
{
	parameters[no] = value;

	switch(no)
	{
		case 1:
			beta0 = value; return;
		case 2:
			beta1 = value; return;
		case 3:
			beta2 = value; return;
		case 4:
			score = value; return;	
	};
};	

//! Checks if there are any minor alleles for the anchor SNP 
bool TwoSNPLogRegModel::checkAnchorSNPData() const
{
	
	return (partnerJointGenotypeCountsCases->counts[1][0] > 0 || partnerJointGenotypeCountsCases->counts[1][1] > 0 || partnerJointGenotypeCountsCases->counts[1][2] > 0 ||
	partnerJointGenotypeCountsCases->counts[2][0] > 0 || partnerJointGenotypeCountsCases->counts[2][1] > 0 || partnerJointGenotypeCountsCases->counts[2][2] > 0 ||
	partnerJointGenotypeCountsControls->counts[1][0] > 0 || partnerJointGenotypeCountsControls->counts[1][1] > 0 || partnerJointGenotypeCountsControls->counts[1][2] > 0 ||
	partnerJointGenotypeCountsControls->counts[2][0] > 0 || partnerJointGenotypeCountsControls->counts[2][1] > 0 || partnerJointGenotypeCountsControls->counts[2][2] > 0 );

};

//! Checks if there are any minor alleles for the partner SNP 
bool TwoSNPLogRegModel::checkPartnerSNPData() const
{
	 return (partnerJointGenotypeCountsCases->counts[0][1] > 0 || partnerJointGenotypeCountsCases->counts[0][2] > 0 ||
		 partnerJointGenotypeCountsCases->counts[1][1] > 0 || partnerJointGenotypeCountsCases->counts[1][2] > 0 ||
		 partnerJointGenotypeCountsCases->counts[2][1] > 0 || partnerJointGenotypeCountsCases->counts[2][2] > 0 ||
		 partnerJointGenotypeCountsControls->counts[0][1] > 0 || partnerJointGenotypeCountsControls->counts[0][2] > 0 ||
		 partnerJointGenotypeCountsControls->counts[1][1] > 0 || partnerJointGenotypeCountsControls->counts[1][2] > 0 ||
		 partnerJointGenotypeCountsControls->counts[2][1] > 0 || partnerJointGenotypeCountsControls->counts[2][2] > 0);

};

//! Checks if there are any minor allele data for the anchor and partner is different.
bool TwoSNPLogRegModel::checkAnchorDiffPartnerSNPData() const
{
		
	return (partnerJointGenotypeCountsCases->counts[0][1] > 0 || partnerJointGenotypeCountsCases->counts[0][2] > 0 ||
	partnerJointGenotypeCountsCases->counts[1][0] > 0 || partnerJointGenotypeCountsCases->counts[1][2] > 0 ||
	partnerJointGenotypeCountsCases->counts[2][0] > 0 || partnerJointGenotypeCountsCases->counts[2][1] > 0 ||
	partnerJointGenotypeCountsControls->counts[0][1] > 0 || partnerJointGenotypeCountsControls->counts[0][2] > 0 ||
	partnerJointGenotypeCountsControls->counts[1][0] > 0 || partnerJointGenotypeCountsControls->counts[1][2] > 0 ||
	partnerJointGenotypeCountsControls->counts[2][0] > 0 || partnerJointGenotypeCountsControls->counts[2][1] > 0);
	
};

//! Fits 2 SNP model to quantitive traits.
bool TwoSNPLinearRegModel::fitModel(double & rss, const bool & fitAnc, const bool & fitPart, const bool & lr)
{
	if(useCovariates) return fitModelCovar(rss, fitAnc, fitPart, lr);

	//Solve equn X^T y = (X^T X)betaHat, to find betaHat, where X is the design matrix
	bool fitOK = false;

	//construct vector X^T y and matrix (X^T X)
	map<unsigned int, double> vectorXTy;
	map<unsigned int, map<unsigned int, double> > matrixXTX;
	list<unsigned char>::const_iterator anc = anchorSNP->noMinorAllelesAllSubjects.begin();
	list<unsigned char>::const_iterator part = partnerSNP->noMinorAllelesAllSubjects.begin();
	list<unsigned char>::const_iterator cc = caseControls.begin();

	double totalQT = 0, ancQT = 0, partQT = 0;
	double sumAnc = 0, sumPart = 0, sumAncPart = 0;
	double sumAncSq = 0, sumPartSq = 0;
	
	unsigned int totalNotMissing = 0;
	
	for(map<unsigned int, double>::const_iterator i = quantitiveTraits->values.begin(); i != quantitiveTraits->values.end(); ++i)
	{
		if((*anc) != 3 && ((*part) != 3 || lr) && i->second != missingQTValue && *cc != 0) // miss out indivs where there is missing SNP data or missing QT data, for standard LR not intersted in missing partner SNP data
		{
			totalQT += i->second;

			totalNotMissing++;			

			if(fitAnc)
			{
				ancQT += (i->second)*(*anc);
				sumAnc += *anc;
				sumAncSq += (*anc)*(*anc);
			};

			if(fitPart)
			{
				partQT += (i->second)*(*part);
				sumPart += *part;
				sumPartSq += (*part)*(*part);
				if(fitAnc) sumAncPart += (*anc)*(*part);			
			};			

			fitOK = true;
		};
		++anc;
		++part;
		++cc;
	};

	if(!fitOK) return false;

	//construct vector XTy
	vectorXTy[1] = totalQT;
	if(fitAnc) vectorXTy[2] = ancQT;
	if(fitPart) vectorXTy[3] = partQT;	

	//construct matrix XTX
	//row 1
	map<unsigned int, double> rowOne;
	rowOne[1] = totalNotMissing;
	if(fitAnc) rowOne[2] = sumAnc;
	if(fitPart) rowOne[3] = sumPart;	
	matrixXTX[1] = rowOne;

	if(fitAnc)
	{
		//row 2
		map<unsigned int, double> rowTwo;
		rowTwo[1] = sumAnc;
		rowTwo[2] = sumAncSq;
		if(fitPart) rowTwo[3] = sumAncPart;	
		matrixXTX[2] = rowTwo;
	};

	if(fitPart)
	{
		//row 3
		map<unsigned int, double> rowThree;
		rowThree[1] = sumPart;
		if(fitAnc) rowThree[2] = sumAncPart;
		rowThree[3] = sumPartSq;	
		matrixXTX[3] = rowThree;
	};

	//solve matrix equation to get betaHat
	map<unsigned int, double> ans = getSolnMatrixEqun(matrixXTX, vectorXTy);

	beta0 = 0; beta1= 0; beta2 = 0;

	//set parameters
	for(map<unsigned int, double>::const_iterator par = ans.begin(); par != ans.end(); ++par)
	{
		setParameter(par->first, par->second);	
	};

	//set residual sum of squares
	anc = anchorSNP->noMinorAllelesAllSubjects.begin();
	part = partnerSNP->noMinorAllelesAllSubjects.begin();
	cc = caseControls.begin();

	double diff;
	double xiTbetaHat; //model fit

	rss = 0;

	//loop thro' subjects and calc difference of model fit with the observed data.
	for(map<unsigned int, double>::const_iterator i = quantitiveTraits->values.begin(); i != quantitiveTraits->values.end(); ++i)
	{	
		if((*anc) != 3 && ((*part) != 3 || lr) && i->second != missingQTValue && *cc != 0) // miss out indivs where there is missing SNP data or missing QT data
		{
			xiTbetaHat = beta0;
			if(fitAnc) xiTbetaHat += beta1*(*anc);
			if(fitPart) xiTbetaHat += beta2*(*part);		

			diff = i->second - xiTbetaHat;
			rss += diff*diff;
		};
		++anc;
		++part;
		++cc;
	};

	if(lr) totalNotMissingSNPs = totalNotMissing;
	else totalNotMissingBetweenSNPs = totalNotMissing;

	return fitOK;
};

//! Fits 2 SNP model to quantitive traits.
bool TwoSNPLinearRegModel::fitModelCovar(double & rss, const bool & fitAnc, const bool & fitPart, const bool & lr)
{
	//Solve equn X^T y = (X^T X)betaHat, to find betaHat, where X is the design matrix
	bool fitOK = false;

	//construct vector X^T y and matrix (X^T X)
	map<unsigned int, double> vectorXTy;
	map<unsigned int, map<unsigned int, double> > matrixXTX;
	list<unsigned char>::const_iterator anc = anchorSNP->noMinorAllelesAllSubjects.begin();
	list<unsigned char>::const_iterator part = partnerSNP->noMinorAllelesAllSubjects.begin();
	list<unsigned char>::const_iterator cc = caseControls.begin();
	list<list<double> >::const_iterator covar = covariateData->covariateDataAllSubjects.begin();
	list<double>::const_iterator cv, cv2;

	double totalQT = 0, ancQT = 0, partQT = 0;
	list<double> covarQT, sumCV, covarAnc, covarPart;
	list<double>::iterator cvQT, sCV, cvAn, cvPa;
	for(cv = covar->begin(); cv != covar->end(); ++cv)
	{
		covarQT.push_back(0); //set initial totals to 0
		sumCV.push_back(0);
		covarAnc.push_back(0);
		covarPart.push_back(0);		
	};

	list<list<double> > covarCovarSums;
	list<double> aCovarCovarSum;
	for(cv = covar->begin(); cv != covar->end(); ++cv)
	{
		aCovarCovarSum.push_back(0); //set initial totals to 0			
	};

	for(cv = covar->begin(); cv != covar->end(); ++cv)
	{
		covarCovarSums.push_back(aCovarCovarSum);
	};

	list<list<double> >::iterator cvcvSums;
	list<double>::iterator cvcvSums2;

	double sumAnc = 0, sumPart = 0, sumAncPart = 0;
	
	double sumAncSq = 0, sumPartSq = 0;
	unsigned int covariateParameterNo;

	unsigned int totalNotMissing = 0;

	for(map<unsigned int, double>::const_iterator i = quantitiveTraits->values.begin(); i != quantitiveTraits->values.end(); ++i, ++anc, ++part, ++cc, ++covar)
	{
		if((*anc) != 3 && ((*part) != 3 || lr) && i->second != missingQTValue && *cc != 0) // miss out indivs where there is missing SNP data or missing QT data, for standard LR not intersted in missing partner SNP data
		{
			totalQT += i->second;

			totalNotMissing++;			

			if(fitAnc)
			{
				ancQT += (i->second)*(*anc);
				sumAnc += *anc;
				sumAncSq += (*anc)*(*anc);
			};

			if(fitPart)
			{
				partQT += (i->second)*(*part);
				sumPart += *part;
				sumPartSq += (*part)*(*part);
				if(fitAnc) sumAncPart += (*anc)*(*part);			
			};			
		
			//do covariate sums		
			cvQT = covarQT.begin();
			sCV = sumCV.begin();
			cvAn = covarAnc.begin();
			cvPa = covarPart.begin();	
			cvcvSums = covarCovarSums.begin();	
			for(cv = covar->begin(); cv != covar->end(); ++cv)
			{
				*cvQT += (i->second)*(*cv);
				*sCV += (*cv);
				*cvAn += (*anc)*(*cv);
				*cvPa += (*part)*(*cv);
			
				//loop thro' covar*covar sums
				cvcvSums2 = cvcvSums->begin();
				for(list<double>::const_iterator cv2 = covar->begin(); cv2 != covar->end(); ++cv2, ++cvcvSums2)
				{
					*cvcvSums2 += (*cv)*(*cv2);
				};

				++cvQT; ++sCV; ++cvAn; ++cvPa; ++cvcvSums;
			};

			fitOK = true;
		};		
	};

	if(!fitOK) return false;

	//construct vector XTy
	vectorXTy[1] = totalQT;
	if(fitAnc) vectorXTy[2] = ancQT;
	if(fitPart) vectorXTy[3] = partQT;	

	//do covariate elements, row 5 onwards
	for(covariateParameterNo = 5, cvQT = covarQT.begin(); cvQT != covarQT.end(); ++cvQT, ++covariateParameterNo)
	{
		vectorXTy[covariateParameterNo] = *cvQT;		
	};

	//construct matrix XTX
	//row 1
	map<unsigned int, double> rowOne;
	rowOne[1] = totalNotMissing;
	if(fitAnc) rowOne[2] = sumAnc;
	if(fitPart) rowOne[3] = sumPart;	
	//do covariate elements, row 5 onwards
	for(covariateParameterNo = 5, sCV = sumCV.begin(); sCV != sumCV.end(); ++sCV, ++covariateParameterNo)
	{
		rowOne[covariateParameterNo] = *sCV;		
	};

	matrixXTX[1] = rowOne;

	if(fitAnc)
	{
		//row 2
		map<unsigned int, double> rowTwo;
		rowTwo[1] = sumAnc;
		rowTwo[2] = sumAncSq;
		if(fitPart) rowTwo[3] = sumAncPart;	
		//do covariate elements, row 5 onwards
		for(covariateParameterNo = 5, cvAn = covarAnc.begin(); cvAn != covarAnc.end(); ++cvAn, ++covariateParameterNo)
		{
			rowTwo[covariateParameterNo] = *cvAn;		
		};
		matrixXTX[2] = rowTwo;
	};

	if(fitPart)
	{
		//row 3
		map<unsigned int, double> rowThree;
		rowThree[1] = sumPart;
		if(fitAnc) rowThree[2] = sumAncPart;
		rowThree[3] = sumPartSq;
		//do covariate elements, row 5 onwards
		for(covariateParameterNo = 5, cvPa = covarPart.begin(); cvPa != covarPart.end(); ++cvPa, ++covariateParameterNo)
		{
			rowThree[covariateParameterNo] = *cvPa;		
		};
		matrixXTX[3] = rowThree;
	};

	//do covariate rows, row 5 onwards
	sCV = sumCV.begin();
	cvAn = covarAnc.begin();
	cvPa = covarPart.begin();
	map<unsigned int, double> rowCovariate;
	unsigned int covariateParameterNo2;
	cvcvSums = covarCovarSums.begin();

	for(covariateParameterNo = 5; sCV != sumCV.end(); ++covariateParameterNo, ++cvcvSums)
	{
		rowCovariate[1] = *sCV; 
		if(fitAnc) rowCovariate[2] = *cvAn; 
		if(fitPart) rowCovariate[3] = *cvPa;

		
		for(covariateParameterNo2 = 5, cvcvSums2 = cvcvSums->begin(); cvcvSums2 != cvcvSums->end(); ++cvcvSums2, ++covariateParameterNo2)
		{
			rowCovariate[covariateParameterNo2] = *cvcvSums2;
		};

		matrixXTX[covariateParameterNo] = rowCovariate;
		++sCV; ++cvAn; ++cvPa;
	};

	//solve matrix equation to get betaHat
	map<unsigned int, double> ans = getSolnMatrixEqun(matrixXTX, vectorXTy);

	beta0 = 0; beta1= 0; beta2 = 0;

	//set parameters
	for(map<unsigned int, double>::const_iterator par = ans.begin(); par != ans.end(); ++par)
	{
		setParameter(par->first, par->second);	
	};

	//set residual sum of squares
	anc = anchorSNP->noMinorAllelesAllSubjects.begin();
	part = partnerSNP->noMinorAllelesAllSubjects.begin();
	cc = caseControls.begin();
	covar = covariateData->covariateDataAllSubjects.begin();

	double diff;
	double xiTbetaHat; //model fit

	rss = 0;

	//loop thro' subjects and calc difference of model fit with the observed data.
	for(map<unsigned int, double>::const_iterator i = quantitiveTraits->values.begin(); i != quantitiveTraits->values.end(); ++i, ++anc, ++part, ++cc, ++covar)
	{	
		if((*anc) != 3 && ((*part) != 3 || lr) && i->second != missingQTValue && *cc != 0) // miss out indivs where there is missing SNP data or missing QT data
		{
			xiTbetaHat = beta0;
			if(fitAnc) xiTbetaHat += beta1*(*anc);
			if(fitPart) xiTbetaHat += beta2*(*part);	
			covariateParameterNo = 5;
			for(cv = covar->begin(); cv != covar->end(); ++cv)
			{
				xiTbetaHat += (*cv)*(getParameter(covariateParameterNo));
				covariateParameterNo++;
			};

			diff = i->second - xiTbetaHat;
			rss += diff*diff;
		};
		
	};

	if(lr) totalNotMissingSNPs = totalNotMissing;
	else totalNotMissingBetweenSNPs = totalNotMissing;

	return fitOK;
};

void TwoSNPLogRegModel::getGradientVector(map<unsigned int, double> & gradientVector, const bool & fitBeta1, const bool & fitBeta2)  const
{
	if(useCovariates) {getGradientVectorCovar(gradientVector, fitBeta1, fitBeta2); return;};

	double ans[3] = {0, 0, 0};

	//loop thro' the number of minor alleles for SNP1, i.e. start SNP
	for(unsigned int snp1 = 0; snp1 <= 2; ++snp1)
	{
		//loop thro' the number of minor alleles for SNP2, i.e. partner SNP
		for(unsigned int snp2 = 0; snp2 <= 2; ++snp2)
		{			
			effects = beta0 + snp1*beta1 + snp2*beta2;
			noCases = (double)partnerJointGenotypeCountsCases->counts[snp1][snp2];
			noControls = (double)partnerJointGenotypeCountsControls->counts[snp1][snp2];
			aNumber = -noCases + (noCases + noControls)/(exp(-effects) + 1);
			
			ans[0] += aNumber;
			if(fitBeta1)
			{
				ans[1] += snp1*aNumber;
				if(fitBeta2) ans[2] += snp2*aNumber;
			}
			else if(fitBeta2) ans[1] += snp2*aNumber;			
				
		};

	};

	gradientVector[1] = ans[0];
	if(fitBeta1)
	{
			gradientVector[2] = ans[1];
			if(fitBeta2) gradientVector[3] = ans[2];
	}
	else if(fitBeta2) gradientVector[2] = ans[1];
	
};


//! Returns the 2nd derivative w.r.t. chosen parameters of the negative log likelihood.
void TwoSNPLogRegModel::getHessianMatrix(map<unsigned int, map<unsigned int, double> > & hessianMatrix, const bool & fitBeta1, const bool & fitBeta2) const
{
	if(useCovariates) {getHessianMatrixCovar(hessianMatrix, fitBeta1, fitBeta2); return;};

	double ans[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}; //col, row
	
	//loop thro' the number of minor alleles for SNP1, i.e. anchor SNP
	for(unsigned int snp1 = 0; snp1 <= 2; ++snp1)
	{
		//loop thro' the number of minor alleles for SNP2, i.e. partner SNP
		for(unsigned int snp2 = 0; snp2 <= 2; ++snp2)
		{			
			effects = beta0 + snp1*beta1 + snp2*beta2;
			expEffects = exp(effects);
			noCases = (double)partnerJointGenotypeCountsCases->counts[snp1][snp2];
			noControls = (double)partnerJointGenotypeCountsControls->counts[snp1][snp2];
			aNumber = ((noCases+noControls)*expEffects)/((expEffects + 1)*(expEffects + 1));

			//only calculate for the parameters we are fitting over
			ans[0][0] += aNumber;

			if(fitBeta1)
			{
				ans[0][1] += snp1*aNumber;	
				ans[1][1] += snp1*snp1*aNumber;

				if(fitBeta2)
				{
					ans[0][2] += snp2*aNumber;
					ans[2][2] += snp2*snp2*aNumber;
					ans[1][2] += snp1*snp2*aNumber;			
				};
			}
			else if(fitBeta2)
			{
				ans[0][1] += snp2*aNumber;	
				ans[1][1] += snp2*snp2*aNumber;
			};

		};
	};

	//The Hessian matrix is symetric so only calculate half and then copy
	if(fitBeta1)
	{
		ans[1][0] = ans[0][1];

		if(fitBeta2)
		{
			ans[2][0] = ans[0][2];

			ans[2][1] = ans[1][2];			
		};
	}
	else if(fitBeta2)
	{
		ans[1][0] = ans[0][1];
	};

	//setup the matrix with calculated values
	map<unsigned int, double> aCol;
	unsigned int noParas = 2;

	if(fitBeta2 && fitBeta1) noParas = 3;

	for(unsigned int col = 0; col < noParas; ++col)
	{		
		for(unsigned int row = 0; row < noParas; ++row)
		{
			aCol[row + 1] = ans[col][row];	
		};
		hessianMatrix[col + 1] = aCol;
	};
};

//! Returns the 2nd derivative w.r.t. chosen parameters of the negative log likelihood with covariates.
void TwoSNPLogRegModel::getHessianMatrixCovar(map<unsigned int, map<unsigned int, double> > & hessianMatrix, const bool & fitBeta1, const bool & fitBeta2) const
{
	double ans[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}; //col, row
	map<pair<unsigned int, unsigned int>, double> covarAns; //starting count from 0, covariates start after other parameters, either 2 or 3

	//loop thro' each subject
	list<unsigned char>::const_iterator anc = anchorSNP->noMinorAllelesAllSubjects.begin();
	list<unsigned char>::const_iterator part = partnerSNP->noMinorAllelesAllSubjects.begin();
	list<list<double> >::const_iterator covar = covariateData->covariateDataAllSubjects.begin();
	list<unsigned char>::const_iterator cc = caseControls.begin();
	list<double>::const_iterator cv1, cv2;

	//unsigned int covariateParameterNo;
	unsigned int cp = 1; //first parameter no. for covariate in covarAns matrix
	if(fitBeta1) cp++;
	if(fitBeta2) cp++;
	//unsigned int cpi1, cpi2;
	double val1;

	unsigned int noParas = cp + covar->size();

	list<double>::const_iterator expEffs = expEffectsAllIndiv.begin();//set up in Vector

	//set 2nd dervis to 0 for covariates
	//for(unsigned int i1 = 1; i1 <= noParas; ++i1)
	//{
	//	for(unsigned int i2 = 1; i2 <= noParas; ++i2)
	//	{
	//		covarAns[make_pair(cpi1, cpi2)] = 0;	
	//	};
	//};

	//loop thro' individuals
	for( ; anc != anchorSNP->noMinorAllelesAllSubjects.end(); ++anc, ++part, ++covar, ++cc)
	{
		if(*anc != 3 && *part != 3 && *cc != 0) //skip missing data
		{

			//effects = beta0 + (*anc)*beta1 + (*part)*beta2;

			//covariateParameterNo = 5;
			////add covariates effects
			//for(list<double>::const_iterator cv = covar->begin(); cv != covar->end(); ++cv)
			//{
			//	effects += (*cv)*(getParameter(covariateParameterNo));
			//	covariateParameterNo++;
			//};

			//expEffects = exp(effects);			

			aNumber = (*expEffs)/((*expEffs + 1)*(*expEffs + 1));
			++expEffs;

			//only calculate for the parameters we are fitting over
			ans[0][0] += aNumber;

			if(fitBeta1)
			{
				ans[0][1] += (*anc)*aNumber;	
				ans[1][1] += (*anc)*(*anc)*aNumber;

				if(fitBeta2)
				{
					ans[0][2] += (*part)*aNumber;
					ans[2][2] += (*part)*(*part)*aNumber;
					ans[1][2] += (*anc)*(*part)*aNumber;			
				};
			}
			else if(fitBeta2)
			{
				ans[0][1] += (*part)*aNumber;	
				ans[1][1] += (*part)*(*part)*aNumber;
			};

			cv1 = covar->begin();

			//add covariate 2nd derivs
			for(unsigned int i1 = 0; i1 < noParas; )
			{
				cv2 = covar->begin();

				for(unsigned int i2 = 0; i2 < noParas; )
				{					
					//update any bit of matrix that is part of covariate in top half of matrix 
					if(i2 >= cp && i2 >= i1)
					{
						if(i1 < cp) //covariate and snp parameter 2nd deriv
						{
							if(i1 == 0) val1 = 1;
							else if(i1 == 1)
							{
								if(fitBeta1) val1 = (*anc);
								else val1 = (*part);
							}
							else if(i1 == 2)
							{							
								val1 = (*part);
							};

							covarAns[make_pair(i1, i2)] += val1*(*cv2)*aNumber;
						}
						else //covariate and covariate 2nd deriv
						{
							covarAns[make_pair(i1, i2)] += (*cv1)*(*cv2)*aNumber;	
						};
					};

					++i2;
					if(i2 > cp) ++cv2; //move to next covariate value
				};

				++i1;
				if(i1 > cp) ++cv1; //move to next covariate value
			};

		};//end of skipping missing data

	};//end of subject loop

	//The Hessian matrix is symetric so only calculate half and then copy
	if(fitBeta1)
	{
		ans[1][0] = ans[0][1];

		if(fitBeta2)
		{
			ans[2][0] = ans[0][2];

			ans[2][1] = ans[1][2];			
		};
	}
	else if(fitBeta2)
	{
		ans[1][0] = ans[0][1];
	};

	//setup the matrix with calculated values
	map<unsigned int, double> aCol;

	//no parameters not including covariates is given by cp
	for(unsigned int col0 = 0; col0 < cp; ++col0)
	{		
		for(unsigned int row0 = 0; row0 < cp; ++row0)
		{
			aCol[row0 + 1] = ans[col0][row0];	
		};
		hessianMatrix[col0 + 1] = aCol;
	};

	//add calculated covariate 2nd derivs
	for(unsigned int col = 0; col < noParas; ++col)
	{	
		if(col < cp) //get predefined column
		{
			aCol = hessianMatrix[col + 1];
		};

		for(unsigned int row = 0; row < noParas; ++row)
		{
			if(col >= cp && col >= row)	//top right half - except for non-covariate
			{
				aCol[row + 1] = covarAns[make_pair(row, col)];
			}
			else if(row >= cp && col < row)	//bottom left half - except for non-covariate
			{
				aCol[row + 1] = covarAns[make_pair(col, row)];
			};
		};

		hessianMatrix[col + 1] = aCol;
	};

};

//! Gradient vector with covariates.
void TwoSNPLogRegModel::getGradientVectorCovar(map<unsigned int, double> & gradientVector, const bool & fitBeta1, const bool & fitBeta2)  const
{
	double ans[3] = {0, 0, 0};
		
	//loop thro' each subject
	list<unsigned char>::const_iterator anc = anchorSNP->noMinorAllelesAllSubjects.begin();
	list<unsigned char>::const_iterator part = partnerSNP->noMinorAllelesAllSubjects.begin();
	list<list<double> >::const_iterator covar = covariateData->covariateDataAllSubjects.begin();
	list<unsigned char>::const_iterator cc = caseControls.begin();

	unsigned int covariateParameterNo;
	unsigned int cp = 2; //first parameter no. for covariate in ans vector
	if(fitBeta1) ++cp;
	if(fitBeta2) ++cp;
	unsigned int cpi = cp;

	for(list<double>::const_iterator cv = covar->begin(); cv != covar->end(); ++cv, ++cpi) gradientVector[cpi] = 0;	

	expEffectsAllIndiv.clear();

	//add covariates parameter derivs
	for( ; anc != anchorSNP->noMinorAllelesAllSubjects.end(); ++anc, ++part, ++covar, ++cc)
	{
		if(*anc != 3 && *part != 3 && *cc != 0) //skip missing data
		{
			effects = beta0 + (*anc)*beta1 + (*part)*beta2;

			covariateParameterNo = 5;
			//add covariates effects
			for(list<double>::const_iterator cv = covar->begin(); cv != covar->end(); ++cv)
			{
				effects += (*cv)*getParameter(covariateParameterNo);
				covariateParameterNo++;
			};

			//cache exp(effects) for use with the hessian matrix
			expEffects = exp(effects); 
			expEffectsAllIndiv.push_back(expEffects);

			if(*cc == 2) aNumber = -1 + expEffects/(expEffects + 1);
			else aNumber = expEffects/(expEffects + 1);		

			ans[0] += aNumber;

			if(fitBeta1)
			{
				ans[1] += (*anc)*aNumber;
				if(fitBeta2) ans[2] += (*part)*aNumber;
			}
			else if(fitBeta2) ans[1] += (*part)*aNumber;	

			cpi = cp;
			//add covariates parameter derivs
			for(list<double>::const_iterator cv = covar->begin(); cv != covar->end(); ++cv, ++cpi)
			{
				gradientVector[cpi] += (*cv)*aNumber;			
			};
		};
	};

	gradientVector[1] = ans[0];
	if(fitBeta1)
	{
			gradientVector[2] = ans[1];
			if(fitBeta2) gradientVector[3] = ans[2];
	}
	else if(fitBeta2) gradientVector[2] = ans[1];
	
};

