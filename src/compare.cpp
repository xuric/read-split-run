/*
  File:       compare.cpp

  Copyright 2015 Jeff Kinne, Yongheng Bai, Brandon Donham.
  Permission to use for academic, non-profit purposes is granted, in which case this 
  copyright notice should be maintained and original authors acknowledged.

  Author:     Jeff Kinne, jkinne@cs.indstate.edu

  Contents:   Helper program to take multiple RSR results files and match up junctions
              that are in more than one results file, and display them together.  Useful
	      for determining junctions that are present in one version of an experiment
	      and not in another.

  To compile: g++ compare.cpp -o comp -O4

  To run:     ./comp linesToSkip supportPosTolerance outputBasename file1 file2 file3 ...
              The files are assumed to be results files produced
	      by the splitPairs.cpp program.  See that file for 
	      the file format.  

	      linesToSkip is how many leading lines to skip (which
	      the splitPairs program puts onto the beginning of its results files - 
	      statistics mostly).  

	      supportPosTolerance is the "buffer boundary" parameter
	      for determining whether two junctions are "the same" - if they are the 
	      same length and differ by at most supportPosTolerance amount they are
	      considered "the same".

	      outPutBasename is the base filename for output files - the output
	      summary files are outPutBasename.comparedResults.txt and outPutBasename.comparisonSummary.txt

  Version history...  Note - executable saved as comp, comp2, comp3, ... (+1 each time testing/debugging a new version)

    8/15/2015 - RSR version 1.0.0 on github

    comp2 - 8/6/2015

    8/6/2015 - add another column to summary that says the max support of any junction.

    8/6/2015 - add another column to summary that shows the total number of known+novel junctions in each file.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <iostream>
#include <limits.h>
#include <vector>
#include <algorithm>
#include <time.h>
#include <unordered_map>
#include <set>

using namespace std;

#include "RSW.h"

vector<RSW_result *> *results;

// counts of the number of novel/known junctions, and max support for any junction.
// pointers to arrays - these are counts for all files.
int *myNovelCount, *myKnownCount, *maxSupportNum;

// parameters that come off the command line
int numResultsFiles, 
  supportPosTolerance;

// keep track of junctions that have already been printed - used in lining up 
// junctions betwen results files
set<RSW_result *> alreadyPrinted;


char sLine[MAX_LINE+1];

char temp[MAX_LINE+1];

/*
  Function:   read_results, reads a results file from disk into memory.

  Parameters: lineToSkip - off the command line, filename - file to open, data_results - 
              vector to store results records into.

  If can't open file, print error message and abort program.
 */
void read_results(int linesToSkip, const char *filename, vector<RSW_result *> &data_results) {
  FILE * f = fopen(filename, "r");
  if (f == NULL) {printf("Error reading from file %s\n", filename); exit(0); }
  
  int result=1; int lineNo=0;
  while (result > 0) {
    int i; char dir;

    RSW_result * rs = new RSW_result; // new record to store a line into

    // get a line
    result = get_line(f, sLine, MAX_LINE);
    if (result < 0) {
      printf("Error reading results file %s, line exceeded %i characters.\n", filename, MAX_LINE);
      break;
    }
    lineNo++;
    if (lineNo <= linesToSkip) continue; // skip past a certain number of lines that aren't records

    // loop through line, split into fields separated by tabs
    char *temp = strtok(sLine, "\t");
    char *temp1; char *sTemp;
    i=0;
    while (temp != NULL) {
      switch (i) {
      case 0: // geneName
	sTemp = (char *)malloc((strlen(temp)+1)*sizeof(char));
	strcpy(sTemp, temp);
	rs->geneName = sTemp;
	break;
      case 1: // chromosome
	sTemp = (char *)malloc((strlen(temp)+1)*sizeof(char));
	strcpy(sTemp, temp);
	rs->chromosome = sTemp;
	break;
      case 2: // supportCount
	rs->supportCount = atoi(temp);
	break;
      case 3: // spliceLength
	rs->spliceLength = atoi(temp);
	break;
      case 4: // min/max positions of boundaries of supporting reads
	temp1 = strstr(temp, "--");
	if (temp1 == NULL) {
	  rs->minSmallSupport = rs->maxLargeSupport = 0;
	}
	else {
	  temp1[0] = '\0';
	  rs->minSmallSupport = atol(temp);
	  rs->maxLargeSupport = atol(temp1+2);
	}
	break;
      case 5: // novel or not
	if (strcmp(temp, "Novel") == 0)
	  rs->novel = true;
	else 
	  rs->novel = false;
	break;
      }
      i++;
      temp = strtok(NULL,"\t");
      if (temp == NULL) break;
    }
    if (i < 6) break; // should be at least 6 fields in each line
    
    // put new record into data_results vector
    data_results.push_back(rs);
  }
  
  fclose(f);
}

/*
  Function: compare_data_results, compare two split junctions, used in sorting them
         
  Sort them based on chromosome name, gene, splice length, and # of supporting reads
 */
bool compare_data_results(RSW_result const *aa, RSW_result const *bb) {
  int temp = strcmp(aa->chromosome, bb->chromosome);
  if (temp < 0) return true;
  else if (temp > 0) return false;

  temp = strcmp(aa->geneName, bb->geneName);
  if (temp < 0) return true;
  else if (temp > 0) return false;

  if (aa->spliceLength < bb->spliceLength) 
    return true;
  else if (bb->spliceLength < aa->spliceLength)
    return false;

  if (aa->minSmallSupport < bb->minSmallSupport) 
    return true;

  return false;
}

/*
  Function: overlapResult, return true/false whether the two splices overlap
  
  For true, must be on same chromosome and gene, have same splice length, and 
  be within supportPosTolerance of each other.
 */
bool overlapResult(RSW_result const *aa, RSW_result const *bb) {
  if (strcmp(aa->chromosome, bb->chromosome) != 0) return false;
  if (strcmp(aa->geneName, bb->geneName) != 0) return false;
  if (aa->spliceLength != bb->spliceLength) return false;

  if (abs(aa->minSmallSupport-bb->minSmallSupport) > supportPosTolerance) return false; 

  return true;
}


// record of all results, used for keeping track of whether we've printed a given
// junction yet or not.
vector<RSW_result *> allResults;



int main(int argc, char *argv[]) {

  if (argc < 5) {
    printf("Usage: ./comp linesToSkip supportPosTolerance outputBasename file1 file2 file3 ...\n"
	   "See splitPairs.cpp file file format of results files.\n");
    exit(0);
  }
  
  numResultsFiles = argc-4;

  supportPosTolerance = atoi(argv[2]);

  // allocate storage for all the results
  results = new vector<RSW_result *>[numResultsFiles];
  myNovelCount = new int[numResultsFiles];
  myKnownCount = new int[numResultsFiles];
  maxSupportNum = new int[numResultsFiles];
  int ii[numResultsFiles];

  // read results from each of the results files given as command-line parameter
  int i;
  for(i=4; i < argc; i++) {
    read_results(atoi(argv[1]), argv[i], results[i-4]);

    sort(results[i-4].begin(), results[i-4].end(), compare_data_results);
  }


  // put all the results into a single array as well.  We go through that 
  // array below.
  for(i=0; i < numResultsFiles; i++) {
    int j;
    ii[i] = 0;
    maxSupportNum[i] = 0; // initialized to 0
    for(j=0; j < results[i].size(); j++) {
      allResults.push_back(results[i][j]);
      if (results[i][j]->supportCount > maxSupportNum[i])
	maxSupportNum[i] = results[i][j]->supportCount;
    }
    myNovelCount[i] = myKnownCount[i] = 0; // initialized to 0
  }


  // sort them, so can match up ones that overlap
  sort(allResults.begin(), allResults.end(), compare_data_results);


  // open output files...
  string sResultsName = argv[3]; sResultsName += ".comparedResults.txt";
  FILE * fResults = fopen(sResultsName.c_str(),"w");
  if (fResults == NULL) {
    printf("Unable to open file %s for writing.\n", sResultsName.c_str());
    exit(1);
  }

  string sSummaryName = argv[3]; sSummaryName += ".comparisonSummary.txt";
  FILE * fSummary = fopen(sSummaryName.c_str(),"w");
  if (fSummary == NULL) {
    printf("Unable to open file %s for writing.\n", sSummaryName.c_str());
    exit(1);
  }

  fprintf(fResults,"Comparison run as: ");
  fprintf(fSummary,"Comparison run as: ");
  for(i=0; i < argc; i++) {
    fprintf(fResults,"%s ",argv[i]);
    fprintf(fSummary,"%s ",argv[i]);
  }
  fprintf(fResults,"\n");
  fprintf(fSummary,"\n");


  // go through the allResults array.  For each one, we check each of the 
  // results files to see if it gets printed with that result.
  int k;
  for(k=0; k < allResults.size(); k++) {
    // if already printed this result, continue
    if (alreadyPrinted.find(allResults[k]) != alreadyPrinted.end()) continue;

    int i_first = -1, i_last = -1;

    fprintf(fResults,"%s\t%s\t", allResults[k]->geneName, allResults[k]->chromosome);

    // check each file and the current spot we are at in that file ...
    for(i=0; i < numResultsFiles; i++) {
      bool printIt = false;
      if (ii[i] < results[i].size() && overlapResult(allResults[k], results[i][ii[i]])) {
	printIt = true;
      }
      if (printIt) {
	fprintf(fResults,"%li\t%li\t%li--%li\t%s\t",
	       results[i][ii[i]]->supportCount,
	       results[i][ii[i]]->spliceLength,
	       results[i][ii[i]]->minSmallSupport,
	       results[i][ii[i]]->maxLargeSupport,
	       results[i][ii[i]]->novel ? "Novel" : "*");
	alreadyPrinted.insert(results[i][ii[i]]);
	ii[i]++; // if printed, move one spot further in that file.  note that it will get printed at some point, so it's safe to only increment when it gets printed.
	if (i_first == -1) i_first = i;
	else i_last = i;
      }
      else {
	fprintf(fResults,"-\t-\t-\t-\t");
      }
    }
    fprintf(fResults,"\n");

    // keep track of number of novel/known that are in exactly one of 
    // the results files.
    if (i_first != -1 && i_last == -1) {
      if (allResults[k]->novel) myNovelCount[i_first]++;
      else myKnownCount[i_first]++;
    }
  }

  printf("Finished running comparison.  Results written to %s.  Summary writen to %s.\n",
	 sResultsName.c_str(), sSummaryName.c_str());
  printf("Summary of results...\n");

  fprintf(fSummary,"file \t#known unique to file \t#novel unique to file\ttotal known+novel in file\tmax support for any read\n");
  printf("file \t#known unique to file \t#novel unique to file\ttotal known+novel in file\tmax support for any read\n");

  for(i=0; i < numResultsFiles; i++) {
    fprintf(fSummary,"%s\t%i\t%i\t%li\t%i\n", argv[i+4], myKnownCount[i], myNovelCount[i], results[i].size(), maxSupportNum[i]);
    printf("%s\t%i\t%i\t%li\t%i\n", argv[i+4], myKnownCount[i], myNovelCount[i], results[i].size(), maxSupportNum[i]);
  }


  fclose(fResults);
  fclose(fSummary);
  

  return 0;
}
