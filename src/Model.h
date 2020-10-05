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


/*! \file Model.h
    \brief This file contains the models used for logistic regression.
    
*/

#ifndef __MODEL
#define __MODEL


#include <map>

using namespace std;

struct JointGenotypeCounts;
struct QuantitiveTraits;
struct SNPData;
struct CovariateData;

#include "Data.h"

//! General Class for Models.
class Model
{
protected:

	map<unsigned int, double> parameters;
	
	//used for fitting covariates and/or linear model
	bool useCovariates;
	CovariateData * covariateData;	
	SNPData * anchorSNP;		
	list<unsigned char> caseControls;

public:

	Model() : parameters(), useCovariates(false) {};
	//setup parameter values and initValues of variables, to be over written in subclass of each model
	//where the variables created and added to the vector for the variables
	Model(map<unsigned int, double> paras) : parameters(paras), useCovariates(false) {};
	
	virtual ~Model() {};

	double getParameter(const unsigned int & no) const {return parameters.find(no)->second;};
	map<unsigned int, double> getParameters() const {return parameters;};
	void setNewParameters(map<unsigned int, double> & paras);
	void setCovariateData(CovariateData * co, bool & uc) {covariateData = co; useCovariates = uc;};
	void setAnchorData(SNPData * anc) {anchorSNP = anc;};	
	void setCaseControls(const list<unsigned char> & ccs) {caseControls = ccs;};

	virtual void setQuantitiveTraits(QuantitiveTraits * qts, double & mQTV) {};

	virtual void setParameter(const unsigned int & no, const double & value) {parameters[no] = value;};	
	virtual double negLogLikelihood() {return 0;};
	virtual void getGradientVector(map<unsigned int, double> & gradientVector, const bool & fitBeta1, const bool & fitBeta2) const {};
	virtual void getHessianMatrix(map<unsigned int, map<unsigned int, double> > & hessianMatrix, const bool & fitBeta1, const bool & fitBeta2) const {};

};

//! Cache of negative log values for beta0 for no. cases and no. controls.
class NegLogLikeCache
{
private:
	

public:
	NegLogLikeCache() {};

	~NegLogLikeCache() {};

	void setValue(const unsigned int & noCases, const unsigned int & noControls, const double & negLog, const double & beta0);
	pair<double, double> getValue(const unsigned int & noCases, const unsigned int & noControls);
};

//! Model for logistic regression on anchor SNP only.
class Model1SNPLogReg : public Model
{
private:

	double beta0, beta1, numCases, numControls;
	
	map<pair<unsigned int, unsigned int>, pair<double, double> > cache; //no. cases, no. controls, neg. log value, beta0
	unsigned int genotypeCountsCases[3];
	unsigned int genotypeCountsControls[3];

	//varibles used for fitting the model
	mutable double effects;
	mutable double expEffects;
	mutable double noCases, noControls;
	mutable double aNumber;	
	mutable list<double> expEffectsAllIndiv;

public:

	Model1SNPLogReg() {};
	
	~Model1SNPLogReg() {};

	void setParameter(const unsigned int & no, const double & value);	
	double negLogLikelihood();
	
	void getGradientVector(map<unsigned int, double> & gradientVector, const bool & fitBeta1, const bool & fitBeta2) const;
	void getHessianMatrix(map<unsigned int, map<unsigned int, double> > & hessianMatrix, const bool & fitBeta1, const bool & fitBeta2) const;

	pair<double, double> getNegLogBeta0(const unsigned int & noCases, const unsigned int & noControls);
	void setGenotypeCounts(unsigned int & ca0, unsigned int & ca1, unsigned int & ca2, unsigned int & co0, unsigned int & co1, unsigned int & co2);

	double negLogLikelihoodCovar();
	void getGradientVectorCovar(map<unsigned int, double> & gradientVector, const bool & fitBeta1, const bool & fitBeta2) const;
	void getHessianMatrixCovar(map<unsigned int, map<unsigned int, double> > & hessianMatrix, const bool & fitBeta1, const bool & fitBeta2) const;
};

//! General Class for Models with 2 SNPs.
class Model2SNPs : public Model
{
protected:
	
	double beta0, beta1, beta2, score; // base, SNP1 effect, SNP2 effect, SNP1xSNP2 effect, correlation of the two SNPs
	
	//used for fitting covariates and/or linear model	
	SNPData * partnerSNP;

public:

	Model2SNPs() {};
	//setup parameter values and initValues of variables, to be over written in subclass of each model
	//where the variables created and added to the vector for the variables
		
	virtual ~Model2SNPs() {};

	void setParameter(const unsigned int & no, const double & value);	
	void setPartnerData(SNPData * par) {partnerSNP = par;};	

	virtual void setJointGenoTypeCounts(JointGenotypeCounts * gjcca, JointGenotypeCounts * gjcco) {};
	
	virtual unsigned int getNoNonMissingSubjects() const {return 0;};
	virtual unsigned int getNoNonMissingSubjectsLR() const {return 0;};
	virtual bool checkAnchorSNPData() const {return true;};
	virtual bool checkPartnerSNPData() const {return true;};
	virtual bool checkAnchorDiffPartnerSNPData() const {return true;};
	virtual bool fitModel(double & rss, const bool & fitAnc, const bool & fitPart, const bool & lr = false) {return 0;};
};

//! Two SNP logistic regression model.
class TwoSNPLogRegModel : public Model2SNPs
{
private:

	JointGenotypeCounts * partnerJointGenotypeCountsCases; //points to the same partner object in the SNP window object
	JointGenotypeCounts * partnerJointGenotypeCountsControls;

protected:

	//varibles used for fitting the model
	mutable double effects;
	mutable double expEffects;
	mutable double noCases, noControls;
	mutable double aNumber;	
	mutable list<double> expEffectsAllIndiv;

public:

	TwoSNPLogRegModel() {};
	
	~TwoSNPLogRegModel() {};
	
	void setJointGenoTypeCounts(JointGenotypeCounts * gjcca, JointGenotypeCounts * gjcco) {partnerJointGenotypeCountsCases = gjcca; partnerJointGenotypeCountsControls = gjcco;};

	double negLogLikelihood();
	
	void getGradientVector(map<unsigned int, double> & gradientVector, const bool & fitBeta1, const bool & fitBeta2) const;
	void getHessianMatrix(map<unsigned int, map<unsigned int, double> > & hessianMatrix, const bool & fitBeta1, const bool & fitBeta2) const;

	bool checkAnchorSNPData() const;
	bool checkPartnerSNPData() const;
	bool checkAnchorDiffPartnerSNPData() const;

	double negLogLikelihoodCovar();
	void getGradientVectorCovar(map<unsigned int, double> & gradientVector, const bool & fitBeta1, const bool & fitBeta2) const;
	void getHessianMatrixCovar(map<unsigned int, map<unsigned int, double> > & hessianMatrix, const bool & fitBeta1, const bool & fitBeta2) const;
};

//! Two SNP linear regression model.
class TwoSNPLinearRegModel : public Model2SNPs
{
private:

	QuantitiveTraits * quantitiveTraits;
	double missingQTValue;
	unsigned int totalNotMissingBetweenSNPs, totalNotMissingSNPs;

public:

	TwoSNPLinearRegModel() {};
	
	~TwoSNPLinearRegModel() {};

	void setQuantitiveTraits(QuantitiveTraits * qts, double & mQTV) {quantitiveTraits = qts; missingQTValue = mQTV;};
	unsigned int getNoNonMissingSubjects() const {return totalNotMissingBetweenSNPs;};
	unsigned int getNoNonMissingSubjectsLR() const {return totalNotMissingSNPs;};

	bool fitModel(double & rss, const bool & fitAnc, const bool & fitPart, const bool & lr = false); 
	bool fitModelCovar(double & rss, const bool & fitAnc, const bool & fitPart, const bool & lr);
};

#endif
