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


/*! \file Data.h
    \brief This file contains classes for manipulating SNP data.
    
*/

#ifndef __DATA
#define __DATA

#include <set>
#include <map>
#include <list>
#include <stack>
#include <fstream>

using namespace std;

class SNPWindow;
class Model2SNPs;
class Model1SNPLogReg;

#include "Model.h"

//! Stores SNP data for all subjects for a given SNP.
struct SNPData
{
	
	list<unsigned char> noMinorAllelesAllSubjects; //ordered list of the number of minor alleles for the subjects	
	string name;
	unsigned int bpPos;
	unsigned int chromosome;

	SNPData(SNPWindow * snpWindow, const string & nm, const unsigned int & bp, const unsigned int & chr);	

	~SNPData()
	{


	};

	void readInNewData(SNPWindow * snpWindow, const string & nm, const unsigned int & bp, const unsigned int & chr);
};

//! Covariate data for all subjects
struct CovariateData
{
	
	//map<unsigned int, list<double> > covariateDataAllSubjects; //covariate number, covariate value for each individual in order
	list<list<double> > covariateDataAllSubjects; //ordered individuals, each with ordered list of covariates
	list<unsigned char> caseControls; 

	CovariateData() {};
	CovariateData(string & covariateFilename, string & covariates, string & fname, double & missingQTValue);	

	~CovariateData()
	{

	};

};

//! Joint Genotype Counts
struct JointGenotypeCounts
{
	unsigned int total;
	unsigned int counts[3][3];

	JointGenotypeCounts()
	{
		total = 0;
		counts[0][0] = 0;
		counts[0][1] = 0;
		counts[0][2] = 0;
		counts[1][0] = 0;
		counts[1][1] = 0;
		counts[1][2] = 0;
		counts[2][0] = 0;
		counts[2][1] = 0;
		counts[2][2] = 0;		
	};

	void addCount(const unsigned int & count1, const unsigned int & count2);
	
};

void setCounts(JointGenotypeCounts * jgc1, JointGenotypeCounts & jgc2);

//! Joint Genotype Counts
struct QuantitiveTraits
{
	map<unsigned int, double> values; //subject no, quantitive trait

	QuantitiveTraits() : values()
	{
		
	};	

};

//! General class for defining how to test the potential partners 
class ScorePotentialPartners
{
private:

protected:

	list<unsigned char> caseControls; //list of values to whether a subject is a case or not: 0=missing, 1=control, 2=case
	Model2SNPs * model2SNPs;
	
	JointGenotypeCounts * partnerJointGenotypeCountsCases;
	JointGenotypeCounts * partnerJointGenotypeCountsControls;
	
	string anchorSNPName, partnerSNPName;
	unsigned int totalNoSubjects;
	SNPWindow * snpWindow;

public:

	ScorePotentialPartners(Model2SNPs *& m2);

	virtual ~ScorePotentialPartners()
	{		
		delete partnerJointGenotypeCountsCases;
		delete partnerJointGenotypeCountsControls;		
	};

	void setupCaseControls(string & fname);
	void updateCaseControlWithMissing(list<unsigned char> & caseControls, bool & linear);
	list<unsigned char> getCaseControls() {return caseControls;};
	void setupQuantitiveTraits(QuantitiveTraits * quantitiveTraits, string & fname, double & missingQTValue);
	void setSNPWindow(SNPWindow * sw) {snpWindow = sw;};
	int checkSNPData(SNPData *& someSNPData) const;
	unsigned int getTotalNoSubjects() {return totalNoSubjects;};

	JointGenotypeCounts * getCaseCounts() {return partnerJointGenotypeCountsCases;};
	JointGenotypeCounts * getControlCounts() {return partnerJointGenotypeCountsControls;};
	string getAnchorSNPName() const {return anchorSNPName;};
	string getPartnerSNPName() const {return partnerSNPName;};
	void setAnchorSNPName(const string & a) {anchorSNPName = a;};
	void setPartnerSNPName(const string & p) {partnerSNPName = p;};

	void setPartnerJointGenotypeCounts(SNPData *& anchorSNPData, SNPData * partnerSNPData, const double & correlation);
	virtual bool calculatePPScore(SNPData *& startSNPData, SNPData * otherSNPData, double & corr, double & score, double & partnerScore) {return 0;};
	virtual void outputPartnerResults(unsigned int & snpID, ofstream & resultsFile) {};
};

//! Calculates the score used to a pick a partner based on cached values
class PPScoreCorrelation : public ScorePotentialPartners
{
protected:

public:

	PPScoreCorrelation(Model2SNPs *& m2) : ScorePotentialPartners(m2)
	{
	
	};

	~PPScoreCorrelation() {};

	bool calculatePPScore(SNPData *& startSNPData, SNPData * otherSNPData, double & corr, double & score, double & partnerScore);
	virtual double calculateCorrScore(const double & corr);
};


//! Calculates the score used to a pick a partner based on correlation evaluations
class PPScoreCorrMulti : public PPScoreCorrelation
{
private:	
	
public:

	PPScoreCorrMulti(Model2SNPs *& m2);
	double calculateCorrScore(const double & corr);

	~PPScoreCorrMulti() {};
};

//! Calculates the score used to a pick a partner based on correlation evaluations with a dominant LR model
class PPScoreCorrDom : public PPScoreCorrelation
{
private:	
	
public:

	PPScoreCorrDom(Model2SNPs *& m2);
	double calculateCorrScore(const double & corr);

	~PPScoreCorrDom() {};
};

//! Calculates the score used to a pick a partner based on correlation evaluations with a recessive LR model
class PPScoreCorrReces : public PPScoreCorrelation
{
private:	
	
public:

	PPScoreCorrReces(Model2SNPs *& m2);
	double calculateCorrScore(const double & corr);

	~PPScoreCorrReces() {};

};

//! Stores SNP data for the SNP window in question.
class SNPWindow
{
private:

	list<unsigned int> basePairs; //base-pair positions of all SNPs in order
	list<unsigned int> chromosomes; //chromosomes of all SNPs in order
	list<string> namesOfSNPs;	
	map<unsigned int, SNPData *> window; //SNP no, SNP Data
	stack<SNPData *> spareSNPs;
	unsigned int windowBPSize; //in base-pair position units (bp units)
	unsigned int windowSNPSize; //if this is set then the window is a fixed size of so many SNPs
	unsigned int totalNoSubjects;
	list<unsigned char> caseControls; //list of bool values to whether a subject is a case or not	

	list<double> partnerScores;
	ScorePotentialPartners * scorePotentialPartners;
	SNPData * anchorSNP;
	SNPData * partnerSNP;
	unsigned int anchorSNPNo;
	map<pair<unsigned int, unsigned int>, pair<double, double> > cachePotentialScores; //snp1, snp2, correlation, joint score

	list<unsigned int> windowSizes;
	list<unsigned int>::const_iterator anchorBPSNP; //will be at the start of the SNP window unless double window is chosen
	list<unsigned int>::const_iterator endBPSNP; //the endSNP iterator is the SNP after the SNP window (i.e. is not in window)
	list<string>::const_iterator nameSNPiter; //used for adding new SNP data
	list<unsigned int>::const_iterator bpSNPiter; //used for adding new SNP data
	list<unsigned int>::const_iterator chrSNPiter; //used for adding new SNP data
	unsigned int addingSNPNo; //used for adding new SNP data
	bool newSNPObject;	
	unsigned int excludeSNP;
	bool excludingSNP;
	double lastAnchorPartnerCorr;
	bool linear;

	//for reading in data from binary files
	ifstream readSNPData;
	unsigned int bitCount;
	int one;
	int aBit;	
	char buffer[1];

public:

	SNPWindow(string & filename, unsigned int & ws, unsigned int & wss, unsigned int & startSNPno, ScorePotentialPartners *& spp, unsigned int & excSNP, bool & lin, bool & lr);

	~SNPWindow()
	{
		delete scorePotentialPartners;		
		readSNPData.close();		

		while(!spareSNPs.empty())
		{
			delete spareSNPs.top();
			spareSNPs.pop();
		};
	};
	
	//methods regarding the SNP descriptions
	void setUpSNPDesciptionData(string & filename);	
	
	void updateEndSNP(); 
	void updateStartSNP();
	unsigned int getNoSNPs() {return basePairs.size();};
	double getLastPartnerScore() {return *(partnerScores.rbegin());};

	void setUpFirstWindow(string & filename, unsigned int & startSNPno);
	void setUpStartAndEndSNPs(unsigned int & startSNPno);
	bool nextWindow();	
	void addNextSNPDataToWindow();
	bool calculatePartnerScores();
	void outputPartnerResults(unsigned int & snpID, ofstream & resultsFile) {scorePotentialPartners->outputPartnerResults(snpID, resultsFile);};

	SNPData * getAnchorSNPData() const {return anchorSNP;};
	void setPartnerSNP(SNPData *& psd) {partnerSNP = psd;};
	SNPData * getPartnerSNPData() const {return partnerSNP;};

	void displayWindowStats();
	unsigned char getNextNoOfMinorAlleles();
	unsigned int getTotalNoSubjects() {return totalNoSubjects;};
	unsigned int getStartBPWindowSNP(unsigned int & startSNPno);
	void advanceToFirstWindow(unsigned int & startSNPno);
	void advanceSNPData();
	SNPData * getNewSNPData();
	
	unsigned int getAnchorSNPBP() const {return anchorSNP->bpPos;};
	unsigned int getAnchorChromosome() const {return anchorSNP->chromosome;};
	unsigned int getPartnerSNPBP() const {return partnerSNP->bpPos;};
	unsigned int getPartnerChromosome() const {return partnerSNP->chromosome;};
	string getAnchorSNPName() const {return scorePotentialPartners->getAnchorSNPName();};
	string getPartnerSNPName() const {return scorePotentialPartners->getPartnerSNPName();};
	string getLastAnchorSNPName() const {return *namesOfSNPs.rbegin();};
	pair<unsigned int, unsigned int> getNoCasesControlsAnchor(Model1SNPLogReg * model1SNPLogReg) const;

	JointGenotypeCounts * getCaseCounts() {return scorePotentialPartners->getCaseCounts();};
	JointGenotypeCounts * getControlCounts() {return scorePotentialPartners->getControlCounts();};
	bool calculatePPScore(SNPData *& startSNPData, SNPData * otherSNPData, double & corr, double & score, double & partnerScore) {return scorePotentialPartners->calculatePPScore(startSNPData, otherSNPData, corr, score, partnerScore);}; 
	string getSNPName(const unsigned int & snpNo) const;

	void setLastAnchorPartnerCorr(const double & c) {lastAnchorPartnerCorr = c;};
	double getLastAnchorPartnerCorr() {return lastAnchorPartnerCorr;};
	
	void writeSNPWindowBPs(ofstream & outFile, const unsigned int & startBP, const unsigned int & endBP, unsigned int & noSNPsInHap);	
	
};



#endif
