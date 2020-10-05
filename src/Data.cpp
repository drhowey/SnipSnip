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


/*! \file Data.cpp
    \brief This file contains the source for manipulating SNP data.
    
*/

#include <string>
#include <map>
#include <set>
#include <list>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>

using namespace std; // initiates the "std" or "standard" namespace

#include "main.h"
#include "Data.h"
#include "ModelFitting.h"

//! Constructor for scoring potential partners
ScorePotentialPartners::ScorePotentialPartners(Model2SNPs *& m2)
	: caseControls(), model2SNPs(m2)
{
		//setup classes for counting cases and control joint genotypes
		partnerJointGenotypeCountsCases = new JointGenotypeCounts();
		partnerJointGenotypeCountsControls = new JointGenotypeCounts();

		model2SNPs->setJointGenoTypeCounts(partnerJointGenotypeCountsCases, partnerJointGenotypeCountsControls);
};

//! Constructor
PPScoreCorrMulti::PPScoreCorrMulti(Model2SNPs *& m2) : PPScoreCorrelation(m2)
{

};

//! Constructor
PPScoreCorrDom::PPScoreCorrDom(Model2SNPs *& m2) : PPScoreCorrelation(m2)
{

};

//! Constructor
PPScoreCorrReces::PPScoreCorrReces(Model2SNPs *& m2) : PPScoreCorrelation(m2)
{

};

//! Sets up family case/control data from the .fam file
void ScorePotentialPartners::setupCaseControls(string & fname)
{
	//try and find the family file and read in data
	unsigned int length = fname.length();
	string famFilename = fname.substr(0, length-4) + ".fam";

	ifstream readFamilyFile;
	readFamilyFile.open(famFilename.c_str());
	if(!readFamilyFile.is_open())
	{
		outErr("Cannot read family file: "); outErr(famFilename); outErr("!\n");
		exit(0);
	};

	string famID, indivID, FatherId, MotherID, sexID, famIndivID;
	string prevFamIndivID = "";
	int phenoType;

	//loop thro' subjects and store the cases
	do{		
		readFamilyFile >> famID >> indivID >> FatherId >> MotherID >> sexID >> phenoType;
		famIndivID = famID + "-" + indivID;

		//do not duplicate the last row
		if(famIndivID != prevFamIndivID) 
		{
			if(phenoType == 2)
			{
					caseControls.push_back(2);					
			}			
			else if(phenoType == 1) caseControls.push_back(1);
			else caseControls.push_back(0);
		};

		prevFamIndivID = famIndivID;
	}while(!readFamilyFile.eof());

	readFamilyFile.close();
};

//!Updates caseControl status with more missing indivs (with missing covariates)
void ScorePotentialPartners::updateCaseControlWithMissing(list<unsigned char> & otherCaseControls, bool & linear)
{
	unsigned int noCases = 0;
	unsigned int noMissing = 0;
	list<unsigned char>::iterator cc = caseControls.begin();
//unsigned int test = 1;
	if(otherCaseControls.size() > 0)
	{
		for(list<unsigned char>::const_iterator oc = otherCaseControls.begin(); oc != otherCaseControls.end(); ++oc, ++cc)
		{	
			if(*cc == 0) ++noMissing;
			else if(*oc == 0)
			{
				*cc = 0;
				++noMissing;	
			}
			else if(*cc == 2) ++noCases;	
		};
	}
	else
	{
		for( ; cc != caseControls.end(); ++cc)
		{	
			if(*cc == 0) ++noMissing; else if(*cc == 2) ++noCases;	
		};
	};

	totalNoSubjects = caseControls.size();

	double totalNoSubjectsNoMiss = caseControls.size() - noMissing;

	double casesPercent = (((double)(noCases))/((double)(totalNoSubjectsNoMiss)))*100;
	double controlsPercent = (((double)(totalNoSubjectsNoMiss - noCases))/((double)(totalNoSubjectsNoMiss)))*100;

	out("Data summary statistics:\n");
	out("Number of subjects: "); out(totalNoSubjectsNoMiss);
	if(noMissing > 0) {out(" (from "); out(totalNoSubjects); out(" in file)"); };
	out("\n");
	if(!linear)
	{
		out("Number of cases: "); out(noCases); out(" ("); out(casesPercent); out("%)\n");
		out("Number of controls: "); out(totalNoSubjectsNoMiss - noCases); out(" ("); out(controlsPercent); out("%)\n");
	};
	out("Number of missing: "); out(noMissing);  out("\n");

};

//! Sets up family case/control data from the .fam file
void ScorePotentialPartners::setupQuantitiveTraits(QuantitiveTraits * quantitiveTraits, string & fname, double & missingQTValue)
{
	//try and find the family file and read in data
	unsigned int length = fname.length();
	string famFilename = fname.substr(0, length-4) + ".fam";

	ifstream readFamilyFile;
	readFamilyFile.open(famFilename.c_str());
	if(!readFamilyFile.is_open())
	{
		outErr("Cannot read family file: "); outErr(famFilename); outErr("!\n");
		exit(0);
	};

	string famID, indivID, FatherId, MotherID, sexID, famIndivID;
	string prevFamIndivID = "";
	double phenoType;
	unsigned int noCases = 0;
	unsigned int subjectNo = 1;
	unsigned int noMissingPhenos = 0;

	//loop thro' subjects and store the cases
	do{		
		readFamilyFile >> famID >> indivID >> FatherId >> MotherID >> sexID >> phenoType;
		famIndivID = famID + "-" + indivID;
	
		//do not duplicate the last row
		if(famIndivID != prevFamIndivID) 
		{
			
			quantitiveTraits->values[subjectNo] = phenoType;
			subjectNo++;
			if(phenoType == missingQTValue)
			{
				++noMissingPhenos; 
				caseControls.push_back(0); //add as missing
			}
			else caseControls.push_back(1); //add all as controls for convience
		};

		prevFamIndivID = famIndivID;
	}while(!readFamilyFile.eof());

	totalNoSubjects = caseControls.size();

	double totalNoSubjectsNoMiss = caseControls.size() - noMissingPhenos;

	readFamilyFile.close();
};

//! Adds a count to the joint genotypes between two SNPs.
void JointGenotypeCounts::addCount(const unsigned int & count1, const unsigned int & count2)
{
	if(count1 != 3 && count2 != 3)
	{
		counts[count1][count2]++;
		++total;
	};
};

void setCounts(JointGenotypeCounts * jgc1, JointGenotypeCounts & jgc2)
{
	jgc1->counts[0][0] = jgc2.counts[0][0];
	jgc1->counts[0][1] = jgc2.counts[0][1];
	jgc1->counts[0][2] = jgc2.counts[0][2];
	jgc1->counts[1][0] = jgc2.counts[1][0];
	jgc1->counts[1][1] = jgc2.counts[1][1];
	jgc1->counts[1][2] = jgc2.counts[1][2];
	jgc1->counts[2][0] = jgc2.counts[2][0];
	jgc1->counts[2][1] = jgc2.counts[2][1];
	jgc1->counts[2][2] = jgc2.counts[2][2];
	jgc1->total = jgc2.total;
};

//! Read in SNP data.
SNPData::SNPData(SNPWindow * snpWindow, const string & nm, const unsigned int & bp, const unsigned int & chr)
{
	unsigned int totalNoSubjects = snpWindow->getTotalNoSubjects();

	for(unsigned int subjectNo = 1; subjectNo <= totalNoSubjects; ++subjectNo)
	{
		noMinorAllelesAllSubjects.push_back(snpWindow->getNextNoOfMinorAlleles());
	};

	name = nm;
	bpPos = bp;	
	chromosome = chr;
};

//! Constucts covariate data and reads in data from file.
CovariateData::CovariateData(string & covariateFilename, string & covariates, string & fname, double & missingQTValue)
{
	//try and find the family file and read in data to create map
	map<string, unsigned int> idToIndivNo; //id name, position of indiv in file
	map<unsigned int, string> indivNoToId; //position of indiv in file, id name

	unsigned int length = fname.length();
	string famFilename = fname.substr(0, length-4) + ".fam";

	ifstream readFamilyFile;
	readFamilyFile.open(famFilename.c_str());
	if(!readFamilyFile.is_open())
	{
		outErr("Cannot read family file: "); outErr(famFilename); outErr("!\n");
		exit(0);
	};

	string famID, indivID, FatherId, MotherID, sexID, famIndivID;
	string prevFamIndivID = "";
	double phenoType;	
	unsigned int subjectNo = 1;

	//loop thro' subjects and store the cases
	do{		
		readFamilyFile >> famID >> indivID >> FatherId >> MotherID >> sexID >> phenoType;
		famIndivID = famID + " " + indivID;

		//do not duplicate the last row
		if(famIndivID != prevFamIndivID) 
		{		
			idToIndivNo[famIndivID] = subjectNo;
			indivNoToId[subjectNo] = famIndivID;
			subjectNo++;
		};

		prevFamIndivID = famIndivID;
	}while(!readFamilyFile.eof());

	unsigned int totalSubjects = subjectNo - 1;
	readFamilyFile.close();

	ifstream readCovariateFile;
	readCovariateFile.open(covariateFilename.c_str());
	if(!readCovariateFile.is_open())
	{
		outErr("Cannot read covariate file: "); outErr(covariateFilename); outErr("!\n");
		exit(0);
	};

	//determine number of covariates
	unsigned int noCovars = 3;
	ifstream readCovariateFileNoCovar(covariateFilename.c_str());
	string firstLine;
	getline(readCovariateFileNoCovar, firstLine);
	unsigned int noCols = 0;
	bool lastWasSpace = true;	

	for(unsigned int ch = 0; ch < firstLine.length(); ++ch)
	{
		if((firstLine.substr(ch, 1) != " " && firstLine.substr(ch, 1) != "\t") && lastWasSpace)
		{
			++noCols;
			lastWasSpace = false;
		};

		if(firstLine.substr(ch, 1) == " " || firstLine.substr(ch, 1) == "\t" ) lastWasSpace = true;
	};

	if(noCols <= 2)
	{
		outErr("Not enough columns in covariate file: "); outErr(covariateFilename); outErr("!\n");
		exit(0);
	};

	noCovars = noCols - 2;
	readCovariateFileNoCovar.close();

	//determine if header present
	ifstream readCovariateFileHeader(covariateFilename.c_str());
	string aStr;

	readCovariateFileHeader >> aStr >> aStr; //first two are expected strings, even if no header

	bool isANumber = false;
	double aNum;

	//if others are numbers assume there is no header
	for(unsigned int h = 3; h <= noCols; ++h)
	{
		readCovariateFileHeader >> aStr;
		aNum = atof(aStr.c_str());
		if(aNum != 0) isANumber = true;
		else if(aNum == 0 && (aStr.substr(0, 1) == "0" || (aStr.length() >= 2 && aStr.substr(0, 2) == "-0"))) isANumber = true; 
	};

	readCovariateFileHeader.close();

	bool headerExists = !isANumber;

	map<string, unsigned int> headerNames;
	map<string, unsigned int>::const_iterator hn;

	if(headerExists)
	{
		ifstream readCovariateFileHeaderNames(covariateFilename.c_str());

		readCovariateFileHeaderNames >> aStr >> aStr; //fam id, id
		readCovariateFile >> aStr >> aStr; //fam id, id

		for(unsigned int h = 1; h <= noCovars; ++h)
		{
			readCovariateFileHeaderNames >> aStr;
			headerNames[aStr] = h;
			readCovariateFile >> aStr;
		};

		readCovariateFileHeaderNames.close();
	};

	//pick covariates to analyse
	set<unsigned int> covariatesToUse;
	set<unsigned int>::const_iterator ctu;
	set<string> someCovariates;
	unsigned int aCol, bCol;

	if(covariates == "") for(unsigned int co = 1; co <= noCovars; ++co)	covariatesToUse.insert(co);
	else
	{		
		//parse covariate string to find covariates to use		
		char * pch;

		pch = strtok(const_cast<char *>(covariates.c_str()),",");
		while(pch != NULL)
		{
			someCovariates.insert(pch);
			pch = strtok(NULL, ",");
		};

		unsigned int dashChr = 0;
		string str1, str2;

		//covariates may still be of form 1-3, so check for this
		for(set<string>::const_iterator sc = someCovariates.begin(); sc != someCovariates.end(); ++sc)
		{
			dashChr = 0;
			
			for(unsigned int cch = 0; cch < (*sc).length(); ++cch)
			{
				if((*sc).substr(cch, 1) == "-") dashChr = cch;
			};

			if(dashChr == 0)
			{
				str1 = *sc;
				aCol = atoi(str1.c_str());				

				if(aCol == 0)
				{
					hn = headerNames.find(str1);
					if(hn != headerNames.end()) aCol = hn->second; 
					else
					{
						outErr("Cannot find covariate column: "); outErr(str1); outErr("!\n");
						exit(0);
					};
				};
				
				covariatesToUse.insert(aCol);
			}
			else
			{
				//find covariates a-b, a to b
				//pch1 = strtok(const_cast<char *>((*sc).c_str()),"-");
				//pch2 = strtok(NULL, "-");
				str1 = (*sc).substr(0, dashChr);
				str2 = (*sc).substr((dashChr+1), ((*sc).length()-dashChr-1));
				aCol = atoi(str1.c_str());
				bCol = atoi(str2.c_str());

				if(aCol == 0)
				{
					hn = headerNames.find(str1);
					if(hn != headerNames.end()) aCol = hn->second; 
					else
					{
						outErr("Cannot find covariate column: "); outErr(str1); outErr("!\n");
						exit(0);
					};
				};

				if(bCol == 0)
				{
					hn = headerNames.find(str2);
					if(hn != headerNames.end()) bCol = hn->second; 
					else
					{
						outErr("Cannot find covariate column: "); outErr(str2); outErr("!\n");
						exit(0);
					};
				};				

				//add all covariates from a to b
				for(unsigned int cab = aCol; cab <= bCol; ++cab) covariatesToUse.insert(cab);
			};
		};
	};

	prevFamIndivID = "";
	list<double> covariateValues;
	double aCovarValue;
	unsigned int coToUse;
	map<string, unsigned int>::const_iterator itin;
	map<unsigned int, list<double> > covariateDataSubjects;
	list<double> missingCovars;
	for(unsigned int m = 1; m <= noCovars; ++m)	missingCovars.push_back(missingQTValue);

	//loop thro' subjects and store the covariate values
	do{
		
		readCovariateFile >> famID >> indivID;
		famIndivID = famID + " " + indivID;

		covariateValues.clear();
		coToUse = 1;
		for(unsigned int i = 1; i <= noCovars; ++i)
		{
			readCovariateFile >> aCovarValue;
			ctu = covariatesToUse.find(i);
			if(ctu != covariatesToUse.end())
			{
				covariateValues.push_back(aCovarValue);
			};
		};

		//add covariate values to subject
		if(famIndivID != prevFamIndivID) 
		{		
			itin = idToIndivNo.find(famIndivID);
			if(itin != idToIndivNo.end())
			{
				//existMissingCovar = false;
				////check for missing phenotypes, if one is missing add all as missing
				//for(list<double>::const_iterator chkc = covariateValues.begin(); chkc != covariateValues.end(); ++chkc)
				//{
				//	if(*chkc == missingQTValue) existMissingCovar = true;
				//};

				//if(existMissingCovar) covariateDataSubjects[itin->second] = missingCovars;
				//else
					
				covariateDataSubjects[itin->second] = covariateValues;
			}
			else
			{
				outErr("Subject \""); outErr(famID); outErr(" "); outErr(indivID); outErr("\" in covariate file "); outErr(covariateFilename); outErr("\n");
				outErr(" cannot be found in family file "); outErr(famFilename); outErr("!\n");
				exit(0);
			};
		};

		prevFamIndivID = famIndivID;
	}while(!readCovariateFile.eof());

	//add covariate data for use later and set missing
	subjectNo = 1;
	bool existMissingCovar;	
	unsigned int noMissingIndivPhenos = 0;
	
	for(map<unsigned int, list<double> >::const_iterator cds = covariateDataSubjects.begin(); cds != covariateDataSubjects.end(); ++cds)
	{
		if(subjectNo != cds->first)
		{
			//numbers "subject No" to cds->first-1 missing
			for(unsigned int j = subjectNo; j < cds->first; ++j)
			{
				covariateDataAllSubjects.push_back(missingCovars);
				caseControls.push_back(0); //add missing individual
				++noMissingIndivPhenos;
			};
		};

		covariateDataAllSubjects.push_back(cds->second);
		//check for missing covariates, if one is missing add all as missing
		existMissingCovar = false;				
		for(list<double>::const_iterator chkc = cds->second.begin(); chkc != cds->second.end(); ++chkc)
		{
			if(*chkc == missingQTValue) existMissingCovar = true; 
		};

		if(existMissingCovar) caseControls.push_back(0); 
		else caseControls.push_back(1); //add all as controls (or missing) for now will be updated later

		subjectNo = cds->first + 1;
	};

	//add missing indivs at end
	for(unsigned int k = subjectNo; k <= totalSubjects; ++k)
	{
		covariateDataAllSubjects.push_back(missingCovars);
		caseControls.push_back(0); //add missing individual
		++noMissingIndivPhenos;
	};

	readCovariateFile.close();
};

//! Replaces SNP data of existing SNP object
void SNPData::readInNewData(SNPWindow * snpWindow, const string & nm, const unsigned int & bp, const unsigned int & chr)
{	
	for(list<unsigned char>::iterator i = noMinorAllelesAllSubjects.begin(); i != noMinorAllelesAllSubjects.end(); ++i)
	{
		*i = snpWindow->getNextNoOfMinorAlleles();
	};

	name = nm;
	bpPos = bp;
	chromosome = chr;
};

//! Gets the next number of minor alleles from file.
unsigned char SNPWindow::getNextNoOfMinorAlleles()
{
	int allele1, allele2;
	unsigned char noMinorAlleles = 0;

	//read in the next piece of data
	if(bitCount == 9)
	{
		
		readSNPData.read(buffer, 1);
		if(readSNPData.eof())
		{			
			outErr("Error: binary SNP file (.bed) is incomplete!\n");
			exit(0);
		};
			
		aBit = buffer[0];
			
		bitCount = 1;
	};

	allele1 = aBit & one; //read the least significant bit				
	aBit = aBit >> 1; //shift bits to the right
	allele2 = aBit & one; //read the new least significant bit				
	aBit = aBit >> 1; //shift bits to the right for next time

	bitCount += 2;	

	//if genotype is encoded 1/0 then the genotype is missing so do not add it
	if(allele1 == 1 && allele2 == 1)
	{	
		noMinorAlleles = 0;
	}
	else if(allele1 == 0 && allele2 == 1)
	{	
		noMinorAlleles = 1;
	}
	else if(allele1 == 0 && allele2 == 0)
	{	
		noMinorAlleles = 2;
	}
	else
		noMinorAlleles = 3; //denotes missing genotype

	return noMinorAlleles;
};

//! Sets up SNP data from the .bim file
void SNPWindow::setUpSNPDesciptionData(string & filename)
{
	//try and find the binary map file, .bim, and read in data
	unsigned int length = filename.length();
	string mapFilename = filename.substr(0, length-4) + ".bim";

	ifstream readMapFile;
	readMapFile.open(mapFilename.c_str());
	if(!readMapFile.is_open())
	{
		outErr("Cannot read map file: "); outErr(mapFilename); outErr("!\n");
		exit(0);
	};

	//assign each subject an ID in the order they appear in the .fam file corresponding to data in the .bed file
	unsigned int snpID = 1;

	string snpIdentifier;
	string prevSnpIdentifier = "";
	unsigned int basePairPosition, chromosome;
	string geneticDistance, alleleName1, alleleName2;

	//loop thro' subjects and store the cases
	do{
		
		readMapFile >> chromosome >> snpIdentifier >> geneticDistance >> basePairPosition >> alleleName1 >> alleleName2;
		
		if(snpIdentifier != prevSnpIdentifier)
		{
				chromosomes.push_back(chromosome);
				basePairs.push_back(basePairPosition);
				namesOfSNPs.push_back(snpIdentifier);				
		};
		
		prevSnpIdentifier = snpIdentifier;
		snpID++;		
	}while(!readMapFile.eof());

	readMapFile.close();

	out("Number of SNPs: "); out(basePairs.size()); out("\n");
};

//! Checks there is sufficient data for cases and controls.
int ScorePotentialPartners::checkSNPData(SNPData *& someSNPData) const
{
	list<unsigned char>::const_iterator nma = someSNPData->noMinorAllelesAllSubjects.begin();
	list<unsigned char>::const_iterator cc = caseControls.begin();
	bool existsCase = false;
	bool existsControl = false;

	while(cc != caseControls.end())
	{

		//check data exists
		if(*nma != 3)
		{
			if(*cc == 2) existsCase = true;
			else if(*cc == 1) existsControl = true;
		};

		if(existsCase && existsControl) return 0;
		++cc;
		++nma;
	};

	if(!existsCase && !existsControl) return -1;
	else if(!existsCase) return -2;
	else return -3;
};

//! Adds the data for the next SNP to the SNP window.
void SNPWindow::addNextSNPDataToWindow()
{
	//create object to store data for the one SNP
	bitCount = 9; //move reading data onto the next byte
	SNPData * someSNPData;
	
	someSNPData = getNewSNPData();
	
	if(!newSNPObject)
	{
		someSNPData->readInNewData(this, *nameSNPiter, *bpSNPiter, *chrSNPiter);
		nameSNPiter++;
		bpSNPiter++;
		chrSNPiter++;
	};

	window[addingSNPNo++] = someSNPData;

	//make sure SNP window size is increased by one to account for extra ignored SNP when using fixed window size
	//do not add extra if excluded SNP is the first SNP
	if(windowSNPSize != 0 && excludingSNP && someSNPData->bpPos == excludeSNP && window.size() > 1)
	{
		endBPSNP++;
		addNextSNPDataToWindow();		
	};
};

//! Returns SNP name for a given SNP no.
string SNPWindow::getSNPName(const unsigned int & snpNo) const
{
	//this is used just for testing purposes and is a list and not a map due to being iterated thro' in main code
	unsigned int snpCount = 1;
	list<string>::const_iterator n = namesOfSNPs.begin();

	while(snpCount < snpNo)
	{		
		snpCount++;
		n++;
	};

	return *n;
};

//! Update the LHS of the SNP window based on BP position.
void SNPWindow::updateStartSNP()
{		
	unsigned int startBasePairPos = 0;
	if(anchorSNP->bpPos > windowBPSize) startBasePairPos = *anchorBPSNP - windowBPSize;	

	//remove all SNPs in window	that are outside the new LHS BP SNP window boundary
	map<unsigned int, SNPData *>::iterator aSNP = window.begin();

	while(aSNP != window.end() && aSNP->second->bpPos < startBasePairPos)
	{
		spareSNPs.push(aSNP->second);
		window.erase(aSNP);
		aSNP = window.begin();
	};	
};

//! Add SNPs to the equired window based on BP size until all the RHS is filled with SNP data. 
void SNPWindow::updateEndSNP()
{		
	unsigned int endBasePairPos = *anchorBPSNP + windowBPSize;	

	//run iterator thro' SNPs until we go past the end of the window or come to the end of the file
	do{

		if(endBPSNP == basePairs.end() || *endBPSNP > endBasePairPos) break;
		if(endBPSNP != basePairs.end()) endBPSNP++;			
			
		if(endBPSNP != basePairs.end()) addNextSNPDataToWindow();
	}while(true);
};

//! Creates the SNP window.
SNPWindow::SNPWindow(string & filename, unsigned int & ws, unsigned int & wss, unsigned int & startSNPno, ScorePotentialPartners *& spp, unsigned int & excSNP, bool & lin, bool & lr)
	: basePairs(), window(), windowBPSize(ws),  windowSNPSize(wss), partnerScores(), windowSizes(), scorePotentialPartners(spp), anchorSNPNo(1), addingSNPNo(1), excludeSNP(excSNP), linear(lin)
{	
	excludingSNP = excludeSNP != 0;

	scorePotentialPartners->setSNPWindow(this);

	//if(!linear)
	//{

	if(lr) caseControls = scorePotentialPartners->getCaseControls();

	//};

	totalNoSubjects = scorePotentialPartners->getTotalNoSubjects();

	//used for reading in SNP binary data
	one = '\1';
	bitCount = 9;
	
	setUpSNPDesciptionData(filename); 
	setUpFirstWindow(filename, startSNPno);
};

//! Moves through SNP data by one SNP.
void SNPWindow::advanceSNPData()
{
	//a new byte is started after each SNP, 4 subject genotypes per byte,
	// so the no. of bytes is rounded up when divided by 4
	unsigned int bufferSize;
	if(totalNoSubjects%4 == 0) bufferSize = totalNoSubjects/4;
	else bufferSize = (unsigned int)(totalNoSubjects/4) + 1;

	char buffer[1];
	for(unsigned int i = 1; i <= bufferSize; ++i)
	{
		readSNPData.read(buffer, 1);
	};

};

//! Gets initial SNP of the pot. part. window when using BP pos and anchor BP. 
unsigned int SNPWindow::getStartBPWindowSNP(unsigned int & startSNPno)
{
	//get BP position of start SNP
	unsigned int posStartBP;
	unsigned int snpNo = 1;
	for(list<unsigned int>::const_iterator i = basePairs.begin(); i != basePairs.end(); ++i)
	{
		if(startSNPno == snpNo) {posStartBP = *i; break;};
		snpNo++;
	};

	unsigned int startBPWindow; //LHS SNP window BP start position
	
	if(posStartBP < windowBPSize) return 1;
	else startBPWindow = posStartBP - windowBPSize;

	snpNo = 1;
	for(list<unsigned int>::const_iterator i = basePairs.begin(); i != basePairs.end(); ++i)
	{
		if(*i >= startBPWindow) return snpNo;
		snpNo++;
	};

	return 1;//this should not be possible
};

//! Moves SNP window onto the chosen first SNP.
void SNPWindow::advanceToFirstWindow(unsigned int & startSNPno)
{
	unsigned int snpCount = 1;
	unsigned int snpToStopAt = startSNPno; //data for this SNP is needed

	if(windowSNPSize == 0)
	{
		snpToStopAt = getStartBPWindowSNP(startSNPno);
	}
	else
	{
		if(windowSNPSize >= startSNPno + 1) snpToStopAt = 1;
		else snpToStopAt = startSNPno - windowSNPSize + 1; //for fixed SNP window size		
	};
		
	do{
		if(snpCount == snpToStopAt) break;

		//read SNP data in for the previous SNP that will not be used
		advanceSNPData();
		addingSNPNo++; //keep this in sync as if the SNP data were added

		snpCount++;

		nameSNPiter++;
		bpSNPiter++;
		chrSNPiter++;
		anchorBPSNP++;
		endBPSNP++;
	}while(bpSNPiter !=  basePairs.end());

	unsigned int noTimesToAdvAnchorBP = startSNPno - snpToStopAt;
	for(unsigned int i = 1; i <= noTimesToAdvAnchorBP; ++i) anchorBPSNP++;
};

//! Sets up the SNP IDs for the start and end SNPs of the first SNP window.
void SNPWindow::setUpStartAndEndSNPs(unsigned int & startSNPno)
{
	if(startSNPno > basePairs.size())
	{
		outErr("The start SNP must not be greater than the total number of SNPs!\n");
		exit(0);
	};

	anchorBPSNP = basePairs.begin();
	anchorSNPNo = startSNPno; //SNP to start analysis

	endBPSNP = basePairs.begin();
	nameSNPiter = namesOfSNPs.begin();
	bpSNPiter = basePairs.begin();
	chrSNPiter = chromosomes.begin();

	advanceToFirstWindow(startSNPno);

	//add SNP data for the first SNP in window
	addNextSNPDataToWindow();

	if(windowSNPSize != 0)
	{	
		unsigned int noSNPsToAdd = windowSNPSize;
		
		//no of SNPs in LHS window also to be added			
		if(startSNPno < windowSNPSize) noSNPsToAdd += startSNPno - 1; 
		else noSNPsToAdd += windowSNPSize - 1; 
		
		//move end SNP on a fixed no of times, the end SNP is not included in the window
		for(unsigned int i = 1; i <= noSNPsToAdd; ++i)
		{		
			if(endBPSNP != basePairs.end())
			{
				endBPSNP++;
				if(endBPSNP == basePairs.end()) break;
				if(i < noSNPsToAdd) addNextSNPDataToWindow(); //do not add data for the last SNP as it is not in the window				
			};
		};	
	}
	else
	{
		//move on end SNP according to the base-pairs window size
		updateEndSNP();
	};

	anchorSNP = window.find(anchorSNPNo)->second;
};

//! Creates the first SNP window, opens SNP file and reads in data for the first window
void SNPWindow::setUpFirstWindow(string & filename, unsigned int & startSNPno)
{
	//try and find the binary pedigree file, .bed, and read in data for the first window
	readSNPData.open(filename.c_str(), ios::binary);
	
	if(!readSNPData.is_open())
	{
		outErr("Cannot read binary pedigree file: "); outErr(filename); outErr("!\n");
		exit(0);
	};

	char buffer[3];
	readSNPData.read(buffer, 3);

	//check the plink magic numbers for the file type
	//3rd number indicates format of genotype data, 1 => subjects x SNPs, 0 => SNPs x subjects
	unsigned int magicNumber1 = buffer[0];
	unsigned int magicNumber2 = buffer[1];

	if(magicNumber1 != 108 || magicNumber2 != 27)
	{
		outErr( "Detected an old version .bed file!\n");
		outErr("Please use PLINK to update the .bed file.\n");
			
		readSNPData.close();		
		exit(0);
	};

	//determine binary file type
	unsigned int mode = buffer[2];
	if(mode == 0)
	{
		outErr("The binary pedigree file must be in SNP-major mode!\n");
		outErr("Please use PLINK to update the .bed file.\n");
			
		readSNPData.close();		
		exit(0);
	};
	
	setUpStartAndEndSNPs(startSNPno);
 
	//record size of window not including the anchor
	windowSizes.push_back(window.size()-1);
	
	//set up joint counts with the partner SNP by calculating LD between start SNP and others
	calculatePartnerScores(); 
};

//! Gets new SNP data from spare SNP data objects or creates a new SNP Data object.
SNPData * SNPWindow::getNewSNPData()
{
	if(spareSNPs.size() == 0)
	{	
		newSNPObject = true;
		return new SNPData(this, *(nameSNPiter++), *(bpSNPiter++), *(chrSNPiter++));
	};

	SNPData * snpData = spareSNPs.top();
	spareSNPs.pop();
	newSNPObject = false;

	return snpData;
};

//! Moves the window on for the next SNP.
bool SNPWindow::nextWindow()
{

	if(windowSNPSize != 0 && anchorSNPNo + 1 > windowSNPSize)
	{
		spareSNPs.push(window.begin()->second);

		//remove first element from the window		
		window.erase(window.begin());
	};

	//move the window along one SNP in the list of all SNPs
	anchorBPSNP++;

	if(windowSNPSize != 0)
	{			
		if(endBPSNP != basePairs.end() && (!excludingSNP || (window.begin()->second)->bpPos != excludeSNP))
		{			
			endBPSNP++;		
			addNextSNPDataToWindow();
		};
	}
	else
	{
		updateEndSNP();
	};

	//record size of window not including the anchor
	windowSizes.push_back(window.size()-1);
	
	//setup new anchor SNP
	anchorSNPNo++;
	anchorSNP = window.find(anchorSNPNo)->second;

	if(windowSNPSize == 0)
	{
		updateStartSNP(); //set LHS snip window based on BP pos
	};

	return calculatePartnerScores();
};

bool SNPWindow::calculatePartnerScores()
{	

	scorePotentialPartners->setAnchorSNPName(anchorSNP->name);
	scorePotentialPartners->setPartnerSNPName("NA");

	//setup genotype counts for the first SNP, in order to compare single SNP analysis with two-SNP
	if(!linear)
	{
		int checkData = scorePotentialPartners->checkSNPData(anchorSNP);

		if(checkData != 0)
		{
			partnerScores.push_back(checkData); //magic number to say SNP had either no cases or controls for this SNP		
			return false;
		};
	};

	map<unsigned int, SNPData *>::iterator i = window.begin();

	double partnerScore = -1;
	double partnerScoreCache = -1;
	double ppCorr; //potential partner correlation
	double ppScore; //potenial partner score
	
	bool partnerSNPFound = false;
	bool partnerSNPFoundLHS = false;
	bool partnerSNPFoundLHSCache = false;
	map<pair<unsigned int, unsigned int>, pair<double, double> >::iterator pp;
	bool firstPPCache = false;
	map<pair<unsigned int, unsigned int>, pair<double, double> >::iterator fpp;
	map<pair<unsigned int, unsigned int>, pair<double, double> >::iterator bestLHSPP;
	SNPData * bestLHSPPSNPData;

	//calculate potential partner score of every SNP in window against anchor SNP in window 
	for( ; i != window.end(); ++i)
	{
		if((!excludingSNP || i->second->bpPos != excludeSNP) && anchorSNPNo != i->first) 
		{
			if(i->first < anchorSNPNo)
			{
				//left hand side window
				pp = cachePotentialScores.find(make_pair(i->first, anchorSNPNo));
				if(!firstPPCache) {fpp = pp; firstPPCache = true;};

				if(pp != cachePotentialScores.end())
				{  
					//found partner in cached values
					if(pp->second.second > partnerScore)
					{
						partnerScore = pp->second.second;
						partnerScoreCache = partnerScore;
						bestLHSPP = pp;
						bestLHSPPSNPData = i->second;					
						partnerSNPFoundLHSCache = true;
					};
				}
				else
				{
					partnerSNPFoundLHS = scorePotentialPartners->calculatePPScore(anchorSNP, i->second, ppCorr, ppScore, partnerScore) || partnerSNPFoundLHS;
				};
			}
			else
			{
				//right hand side window
				partnerSNPFound = scorePotentialPartners->calculatePPScore(anchorSNP, i->second, ppCorr, ppScore, partnerScore) || partnerSNPFound;
				cachePotentialScores[make_pair(anchorSNPNo, i->first)] = make_pair(ppCorr, ppScore);
			};
			
		};

	};

	//update cache of scores
	
	//better partner not found in RHS window so set
	if(!partnerSNPFound && partnerSNPFoundLHSCache && partnerScoreCache == partnerScore)
	{
		scorePotentialPartners->setPartnerJointGenotypeCounts(anchorSNP, bestLHSPPSNPData, bestLHSPP->second.first);
	};

	//remove first pp
	if(cachePotentialScores.size() > 0)
	{
		fpp = cachePotentialScores.begin();
		cachePotentialScores.erase(fpp); //need to erase more or less if a bp window
	};
	
	if(partnerSNPFound || partnerSNPFoundLHSCache || partnerSNPFoundLHS)
	{
		//record partner score
		partnerScores.push_back(partnerScore);			
		return true;
	}
	else
	{
		if(window.size() < 2) partnerScores.push_back(-5); //magic number to say SNP window was too small to find a partner SNP
		else partnerScores.push_back(-4); //magic number to say no suitable partner SNP could be found in window
		return false;
	};
};

double PPScoreCorrelation::calculateCorrScore(const double & corr)
{
	return 0;
};

//! Returns the score of a given correlation
double PPScoreCorrMulti::calculateCorrScore(const double & corr)
{
	double corr2 = corr*corr;
	double corr3 = corr2*corr;
	double corr4 = corr3*corr;
	double corr5 = corr4*corr;
	
	if(corr <= 0.025) return 26.47733*corr; 
	else if(corr <= 0.075) return -120.870*corr + 5893.893*corr2; 
    else if(corr <= 0.975) return -58.86439 + 1406.49125*corr - 4464.18906*corr2 + 6465.82903*corr3 - 4896.66871*corr4 + 1548.48253*corr5; 
    else return 36.40664 - 36.40664*corr; 
};

//! Returns the score of a given correlation using a dominant model
double PPScoreCorrDom::calculateCorrScore(const double & corr)
{
	double corr2 = corr*corr;
	double corr3 = corr2*corr;

	if(corr <= 0.075) return 7.99162*corr + 651.54268*corr2; 
	else if(corr <= 0.925) return  -51.84352 +  851.88206*corr - 1430.52744*corr2 + 624.31772*corr3;
	else return 992.4562  -1974.9955*corr + 982.5393*corr2; 
};

//! Returns the score of a given correlation using a recessive model
double PPScoreCorrReces::calculateCorrScore(const double & corr)
{
	double corr2 = corr*corr;
	double corr3 = corr2*corr;

	if(corr <= 0.025) return 11.02800*corr;
	else if(corr <= 0.075) return -58.94837*corr + 2799.05496*corr2; 
	else if(corr <= 0.925) return -36.55695 + 719.91359*corr - 1118.65846*corr2 + 425.33700*corr3;
	else return 1889.211 - 3804.103*corr + 1914.891*corr2;
};

//! Calculates the LD using Pearson's corellation coeffient, r^2, as given by Welleck and Ziegler.
bool PPScoreCorrelation::calculatePPScore(SNPData *& startSNPData, SNPData * otherSNPData, double & corr, double & score, double & partnerScore)
{
	JointGenotypeCounts jointGenotypeCountsCases;
	JointGenotypeCounts jointGenotypeCountsControls;

	//loop thro' all the subjects for each SNP
	list<unsigned char>::const_iterator snp1 = startSNPData->noMinorAllelesAllSubjects.begin();
	list<unsigned char>::const_iterator snp = otherSNPData->noMinorAllelesAllSubjects.begin();
	list<unsigned char>::const_iterator subs = caseControls.begin();

	//count the minor alleles for each subject
	do{

		//count cases and controls separately
		if(*subs == 2) jointGenotypeCountsCases.addCount(*snp1, *snp);
		else if(*subs == 1) jointGenotypeCountsControls.addCount(*snp1, *snp);

		++snp1; ++snp; ++subs;
	}while(snp1 != startSNPData->noMinorAllelesAllSubjects.end());

	double jointFreqs[3][3];
	double total = (double)(jointGenotypeCountsCases.total + jointGenotypeCountsControls.total);

	//too much missing data to calculate LD
	if(total == 0)
	{	
		return false;
	};

	//convert counts to freqs
	jointFreqs[0][0] = (double)(jointGenotypeCountsCases.counts[0][0] + jointGenotypeCountsControls.counts[0][0])/total;
	jointFreqs[0][1] = (double)(jointGenotypeCountsCases.counts[0][1] + jointGenotypeCountsControls.counts[0][1])/total;
	jointFreqs[0][2] = (double)(jointGenotypeCountsCases.counts[0][2] + jointGenotypeCountsControls.counts[0][2])/total;
	jointFreqs[1][0] = (double)(jointGenotypeCountsCases.counts[1][0] + jointGenotypeCountsControls.counts[1][0])/total;
	jointFreqs[1][1] = (double)(jointGenotypeCountsCases.counts[1][1] + jointGenotypeCountsControls.counts[1][1])/total;
	jointFreqs[1][2] = (double)(jointGenotypeCountsCases.counts[1][2] + jointGenotypeCountsControls.counts[1][2])/total;
	jointFreqs[2][0] = (double)(jointGenotypeCountsCases.counts[2][0] + jointGenotypeCountsControls.counts[2][0])/total;
	jointFreqs[2][1] = (double)(jointGenotypeCountsCases.counts[2][1] + jointGenotypeCountsControls.counts[2][1])/total;
	jointFreqs[2][2] = (double)(jointGenotypeCountsCases.counts[2][2] + jointGenotypeCountsControls.counts[2][2])/total;

	//plug in the frequencies into the formula to calculate the correlation
	double sum1 = jointFreqs[1][0] + jointFreqs[1][1] + jointFreqs[1][2];
	double sum2 = jointFreqs[2][0] + jointFreqs[2][1] + jointFreqs[2][2];
	double sum3 = jointFreqs[0][1] + jointFreqs[1][1] + jointFreqs[2][1];
	double sum4 = jointFreqs[0][2] + jointFreqs[1][2] + jointFreqs[2][2];

	double e1 = sum1 + 2*sum2;
	double e2 = sum3 + 2*sum4;
	double c = jointFreqs[1][1] + 2*(jointFreqs[1][2] + jointFreqs[2][1]) + 4*jointFreqs[2][2] - e1*e2;
	double v1 = sum1 + 4*sum2 - e1*e1;
	double v2 = sum3 + 4*sum4 - e2*e2;

	if(v1 == 0 || v2 == 0) corr = 1; //avoid nans
	else corr = (c*c)/(v1*v2);

	//get corresponding score for the calculated correlation
	score =	calculateCorrScore(corr);

	//set new partner SNP if on the same chromosome
	if(score > partnerScore && startSNPData->chromosome == otherSNPData->chromosome)
	{		
		partnerScore = score;
		partnerSNPName = otherSNPData->name;
		snpWindow->setLastAnchorPartnerCorr(corr);
		snpWindow->setPartnerSNP(otherSNPData);

		setCounts(partnerJointGenotypeCountsCases, jointGenotypeCountsCases); 
		setCounts(partnerJointGenotypeCountsControls, jointGenotypeCountsControls);

		return true;
	};

	return false;
};

//! Sets the joint genotype data for the given anchor and partner SNP 
void ScorePotentialPartners::setPartnerJointGenotypeCounts(SNPData *& anchorSNPData, SNPData * partnerSNPData, const double & correlation)
{
	JointGenotypeCounts jointGenotypeCountsCases;
	JointGenotypeCounts jointGenotypeCountsControls;

	//loop thro' all the subjects for each SNP
	list<unsigned char>::const_iterator snp1 = anchorSNPData->noMinorAllelesAllSubjects.begin();
	list<unsigned char>::const_iterator snp = partnerSNPData->noMinorAllelesAllSubjects.begin();
	list<unsigned char>::const_iterator subs = caseControls.begin();

	//count the minor alleles for each subject
	do{

		//count cases and controls separately
		if(*subs == 2) jointGenotypeCountsCases.addCount(*snp1, *snp);
		else if(*subs == 1) jointGenotypeCountsControls.addCount(*snp1, *snp);

		++snp1; ++snp; ++subs;
	}while(snp1 != anchorSNPData->noMinorAllelesAllSubjects.end());

	setPartnerSNPName(partnerSNPData->name);
	snpWindow->setLastAnchorPartnerCorr(correlation);
	snpWindow->setPartnerSNP(partnerSNPData);

	setCounts(partnerJointGenotypeCountsCases, jointGenotypeCountsCases); 
	setCounts(partnerJointGenotypeCountsControls, jointGenotypeCountsControls);
};

//! Returns no. of cases and controls without missing Anchor data, and sets genotype counts.
pair<unsigned int, unsigned int> SNPWindow::getNoCasesControlsAnchor(Model1SNPLogReg * model1SNPLogReg) const
{
	unsigned int noCases = 0;
	unsigned int noControls = 0;
	unsigned int noCases0 = 0;
	unsigned int noControls0 = 0;
	unsigned int noCases1 = 0;
	unsigned int noControls1 = 0;
	unsigned int noCases2 = 0;
	unsigned int noControls2 = 0;

	list<unsigned char>::const_iterator snp = anchorSNP->noMinorAllelesAllSubjects.begin();

	for(list<unsigned char>::const_iterator cc = caseControls.begin(); cc != caseControls.end(); ++cc, ++snp)
	{
		//count non missing
		if(*snp != 3)
		{
			if(*cc == 2) noCases++;
			else if(*cc == 1) noControls++;
		};

		//set genotype counts
		if(*snp == 0)
		{
			if(*cc == 2) noCases0++;
			else if(*cc == 1) noControls0++;
		}
		else if(*snp == 1)
		{
			if(*cc == 2) noCases1++;
			else if(*cc == 1) noControls1++;
		}
		else if(*snp == 2)
		{
			if(*cc == 2) noCases2++;
			else if(*cc == 1) noControls2++;
		};
	};

	model1SNPLogReg->setGenotypeCounts(noCases0, noCases1, noCases2, noControls0, noControls1, noControls2);

	return make_pair(noCases, noControls);
};

//! Write BPs between and including given BPs
void SNPWindow::writeSNPWindowBPs(ofstream & outFile, const unsigned int & startBP, const unsigned int & endBP, unsigned int & noSNPsInHap)
{
	for(list<unsigned int>::const_iterator i = basePairs.begin(); i != basePairs.end(); ++i)
	{
		if(*i >= startBP && *i <= endBP)
		{
			outFile << *i << "\n";
			noSNPsInHap++;
		}
		else if(*i > endBP) break;
	};
};

//! Displays stats of the window size, mean, sd etc.
void SNPWindow::displayWindowStats()
{
	double scoreMean = 0;
	double scoreSD = 0;
	double scoreTotal = 0;
	double minScore = 0;
	double maxScore = 0;

	for(list<double>::const_iterator pc = partnerScores.begin(); pc != partnerScores.end(); ++pc)
	{
		if(*pc >= 0)
		{
			scoreMean += *pc;
			scoreTotal++;
			if(*pc < minScore || scoreTotal == 1) minScore = *pc;
			if(*pc > maxScore) maxScore = *pc;
		};
	};

	scoreMean /= scoreTotal;

	for(list<double>::const_iterator pc = partnerScores.begin(); pc != partnerScores.end(); ++pc)
	{
		if(*pc >= 0)
		{
			scoreSD += (*pc - scoreMean)*(*pc - scoreMean);
		};
	};

	scoreSD /= scoreTotal;
	scoreSD = sqrt(scoreSD);

	double winMean = 0;
	double winSD = 0;
	double minWin = (*(windowSizes.begin()))*2;
	double maxWin = (*(windowSizes.begin()))*2;

	for(list<unsigned int>::const_iterator ws = windowSizes.begin(); ws != windowSizes.end(); ++ws)
	{		
		winMean += *ws;
		if(*ws < minWin) minWin = *ws;
		if(*ws > maxWin) maxWin = *ws;
	};

	winMean /= (double)(windowSizes.size());
	
	for(list<unsigned int>::const_iterator ws = windowSizes.begin(); ws != windowSizes.end(); ++ws)
	{		
		winSD += (*ws - winMean)*(*ws - winMean);		
	};

	winSD /= (double)(windowSizes.size());
	winSD = sqrt(winSD);

	out("\nSNP window summary statistics:\n");
	out("Mean: "); out(winMean); out("\n");
	out("Standard deviation: "); out(winSD); out("\n");
	out("Range: ("); out(minWin); out(", "); out(maxWin); out(")\n\n");

	out("Score of anchor SNP and partner SNP summary statistics:\n");
	out("Mean: "); out(scoreMean); out("\n");
	out("Standard Deviation: "); out(scoreSD); out("\n");
	out("Range: ("); out(minScore); out(", "); out(maxScore); out(")\n\n");

};

