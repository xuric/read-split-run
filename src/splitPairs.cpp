/*
  File:        splitPairs.cpp

  Copyright 2015 Jeff Kinne, Yongheng Bai, Brandon Donham.
  Permission to use for academic, non-profit purposes is granted, in which case this
  copyright notice should be maintained and original authors acknowledged.

  Author:      Jeff Kinne, jkinne@cs.indstate.edu

  Contents:    Program to take alignments of pieces of unaligned reads and determine
               which could be "split" or "matched" pairs - indicating that the
               unaligned read resulted from a splice.  Also, determine which
               matched pairs support each other (resulted from the same splice junction).

  To compile: g++ splitPairs.cpp -o sp -O4 -std=c++11

  To run:     ./sp options.txt

              Where options.txt is an options file.  If the program is run
              with no command-line arguments it by default processes
              RSW_test.txt with the same parameters as the scripts downloaded
              from Yongsheng Bai's website.  When you run the program it
              prints which files output is written to.  See readOptionsFromFile
              function for the order of the parameters in the options file.

Modification history...  

11/24/2015 - Modify RSW.h and this file to print out the actual sequence for split pairs.
             This will appear in the output in the .results files.  
           - Update comments in this file and RSW.h.  Replace %i with %li where required
             for 64 bit systems (when printing results of time(NULL) or vector.size())

8/15/2015 - version 1.0.0 of RSR on github.

* sp4 - 7/30/2015
* 7/30/2015 - fix bug in computing supporting reads that appeared in sp3.
* 7/30/2015 - fixed bug where compare_data was being used in main loop that computes
              matched pairs; compare_data was written to be used in sorting, but was
	      being used looking for exact match between lines.  Just stopped using
	      it in the main loop and put in the if tests there (to avoid repeating
	      this mistake in the future).
* 7/30/2015 - had been working on using openmp to parallelize the main loops.  was
              getting seg fault for unknown reason, so commented out all #pragma
	      lines that were part of that effort.
* sp3 - 7/29/2015
* 7/29/2015 - Jeff Kinne - fixed bug that was reporting the number of supporting reads
              as one two small for each junction.
* 7/9/2015 - Jeff Kinne - store strings in stringTable to avoid duplicating the 
             strings.  Should save memory and speed up (because compare string pointers
	     rather than strcmp).
* 7/8/2015 - started tracking mod history...
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <iostream>
#include <limits.h>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <time.h>
#include <sys/resource.h>
//#include <omp.h>
using namespace std;

#include "RSW.h"


// parameters input from options file
int maxDistance;  // max difference in aligned pieces to be considered matched pair
int sampleLength; // length of reads in this data set
int minSpliceLength; // minimum size of splice to consider
int supportPosTolerance; // "buffer boundary" parameter for determining supporting reads
int minSupportingReads;  // only output junctions with at least this many supporting reads
char const *sampleDataFile; // input file name
char const *refFlatFile;    // refFlat file of gene locations
char const *refFlatBoundaryFile; // refFlat file of known intron/extron boundaries
char const *resultsBaseName;     // base file name used for output file names

char buff[MAX_STR_LEN];

char options[MAX_STR_LEN];

time_t beginTime, endTime;  // for keeping track of running time of program

// string table is used to reduce memory usage of program.  for any string
// we need we store it in the string table and then only use the pointer to the
// string any time we need it.  So two records from the input file with the same
// read id won't take up twice as much memory.  See reading of data file for
// how this works.
unordered_set<string> stringTable;

// used in printing of results to avoid printing multiple junction sites
// multiple times
unordered_set<const char *> readIdsReported;

// data read into the program
vector<struct RSW> data; // from input file of alignments of pieces of unaligned reads
vector<struct RSW_Known> data_known; // from refFlat
vector<struct RSW_Boundaries> data_boundaries; // from refFlat intron/extron boundaries

const char unfound_string[100] = "UNFOUND_";   // gene name of any junction outside of genes

vector<RSW_splice *> data_splice; // used to store possible jucntions, see RSW.h for RSW_splice definition

int numDifferentReads; // counter...

// function not currently used
int compute_hash(RSW *d) {
  int h = 1, i;
  for(i=0; d->id[i] != '\0'; i++) h = h * d->id[i] % SHRT_MAX;
  for(i=0; d->chromosome[i] != '\0'; i++) h = h * d->chromosome[i] % SHRT_MAX;
  h = h * d->direction % SHRT_MAX;
  return h;
}


char sLine[MAX_LINE+1];
char temp[MAX_LINE+1];

/*
  Function: read_data, read in data file into data vector

  Parameters: filename - file to open and read

  Note: if file is .gz or .lrz then attempt to unzip before reading.  This will
  only work if gunzip and/or lrunzip can be run from the current directory.
*/
void read_data(const char *filename) {
  // open file for reading (from pipe if trying to unzip)
  FILE *f;
  int len = strlen(filename);
  if (len > 3 && strcmp(filename+len-3,".gz")==0) {
    sprintf(temp,"gunzip -c %s", filename);
    f = popen(temp, "r");
  }
  else if (len > 4 && strcmp(filename+len-4,".lrz")==0) {
    sprintf(temp,"cat %s | ./lrunzip", filename);
    f = popen(temp, "r");
  }
  else f = fopen(filename, "r");
  if (f == NULL) {printf("Error reading from file %s\n", filename); exit(0); }

  // read data file one line at a time.
  int result=1;
  while (result > 0) {
    int i; char dir;
    RSW *r = new RSW;
    result = get_line(f, sLine, MAX_LINE);
    if (result < 0) {
      printf("Error reading data file %s, line exceeded %i characters.\n", filename, MAX_LINE);
      delete r;
      break;
    }

    // break line into fields, separated by tab
    char *tempA = strtok(sLine, "\t");
    i=0;
    while (tempA != NULL) {
      string temp = tempA; temp.shrink_to_fit();
      pair<unordered_set<string>::iterator,bool> result;
      switch (i) {
      case 0: // id of read
        // put into string table and store pointer to string.  note
        // that insert just returns a pointer if the string already was in the string table.
	result = stringTable.insert(temp);
	r->id = (result.first)->c_str();
	break;
      case 1: // side
	r->side = temp[0];
	break;
      case 2: // length of piece
	r->length = atoi(temp.c_str());
	break;
      case 3: // direction
	r->direction = temp[0];
	break;
      case 4: // chromosome
	result = stringTable.insert(temp);
	r->chromosome = (result.first)->c_str();
	break;
      case 5: // position
	r->position = atol(temp.c_str());
	break;
      case 6: // sequence 
        result = stringTable.insert(temp);
        r->sequence = (result.first)->c_str();
        break;
      case 8: // count, unused currently
	r->count = atoi(temp.c_str());
	break;
      }
      i++;
      tempA = strtok(NULL,"\t");
      if (tempA == NULL) break;
    }
    if (i < 9) {delete r; break;}

    r->hash = compute_hash(r); // not currently used
    data.push_back(*r); // save into vector
    delete r;

    //if (data.size() >= 15000000) break; // cut off early, for debugging to prevent program from running for too long.
  }

  fclose(f);
}

/*
  Function: read_knownGene, similar to read_data but read the format of the
            refFlat file of gene locations
*/
void read_knownGene(const char *filename) {
  FILE * f = fopen(filename, "r");
  if (f == NULL) {printf("Error reading from file %s\n", filename); exit(0); }

  int result=1;
  while (result > 0) {
    int i; char dir;
    RSW_Known * rk = new RSW_Known;
    result = get_line(f, sLine, MAX_LINE);
    if (result < 0) {
      printf("Error reading data file %s, line exceeded %i characters.\n", filename, MAX_LINE);
      delete rk;
      break;
    }
    char *tempA = strtok(sLine, "\t");
    i=0;
    while (tempA != NULL) {
      string temp = tempA; temp.shrink_to_fit();
      pair<unordered_set<string>::iterator,bool> result;
      switch (i) {
      case 0:
	result = stringTable.insert(temp);
	rk->id1 = (result.first)->c_str();
	break;
      case 1:
	result = stringTable.insert(temp);
	rk->id2 = (result.first)->c_str();
	break;
      case 2:
	result = stringTable.insert(temp);
	rk->chromosome = (result.first)->c_str();
	break;
      case 3:
	rk->direction = temp[0];
	break;
      case 4:
	rk->position1 = atol(temp.c_str());
	break;
      case 5:
	rk->position2 = atol(temp.c_str());
	break;
      }
      i++;
      tempA = strtok(NULL,"\t");
      if (tempA == NULL) break;
    }
    if (i < 11) {delete rk; break;}

    if (rk->position1 > rk->position2) {
      int t = rk->position1;
      rk->position1 = rk->position2;
      rk->position2 = t;
    }
    data_known.push_back(*rk);
    delete rk;
  }

  fclose(f);
}

/*
  Function: read_boundaries, similar to read_data but read the format of the
            refFlat file of intron/extron boundaries.
*/
void read_boundaries(const char *filename) {
  FILE * f = fopen(filename, "r");
  if (f == NULL) {printf("Error reading from file %s\n", filename); exit(0); }

  int result=1;
  while (result > 0) {
    int i; char dir;
    RSW_Boundaries * rk = new RSW_Boundaries;
    result = get_line(f, sLine, MAX_LINE);
    if (result < 0) {
      printf("Error reading data file %s, line exceeded %i characters.\n", filename, MAX_LINE);
      delete rk;
      break;
    }
    char *tempA = strtok(sLine, "\t");
    i=0;
    while (tempA != NULL) {
      string temp = tempA; temp.shrink_to_fit();
      pair<unordered_set<string>::iterator,bool> result;
      switch (i) {
      case 0:
	result = stringTable.insert(temp);
	rk->id1 = (result.first)->c_str();
	break;
      case 1:
	result = stringTable.insert(temp);
	rk->id2 = (result.first)->c_str();
	break;
      case 2:
	result = stringTable.insert(temp);
	rk->chromosome = (result.first)->c_str();
	break;
      case 3:
	rk->direction = temp[0];
	break;
      case 11:
	rk->length = atoi(temp.c_str());
	break;
      case 12:
	char *tempB = (char *) malloc(sizeof(char)* (temp.size()+1));
	strcpy(tempB,temp.c_str());
	char *temp1 = strstr(tempB, "--");
	if (temp1 == NULL) {
	  rk->position1 = rk->position2 = 0;
	}
	else {
	  temp1[0] = '\0';
	  rk->position1 = atol(tempB);
	  rk->position2 = atol(temp1+2);
	}
	free(tempB);
	break;
      }
      i++;
      tempA = strtok(NULL,"\t");
      if (tempA == NULL) break;
    }
    if (i < 13) {delete rk; break;}

    data_boundaries.push_back(*rk);
    delete rk;
  }

  fclose(f);
}

/*
  Function:   compare_data, used for sorting input data

  Sorts based on id, direction, chromosome - all must be the same
  if two pieces were from the same read
*/
bool compare_data(RSW const &aa, RSW const &bb) {
  int temp = aa.id-bb.id;
  if (temp < 0) return true;
  else if (temp > 0) return false;

  if (aa.direction < bb.direction) 
    return true;
  else if (bb.direction < aa.direction)
    return false;

  temp = aa.chromosome-bb.chromosome;
  if (temp < 0) return true;
  else if (temp > 0) return false;

  return false;
}

/*
  Function:   compare_dataToSort, used for sorting input data

  Sorts based on id, direction, chromosome, position
*/
bool compare_dataToSort(RSW const &aa, RSW const &bb) {
  int temp = aa.id-bb.id;
  if (temp < 0) return true;
  else if (temp > 0) return false;

  if (aa.direction < bb.direction) 
    return true;
  else if (bb.direction < aa.direction)
    return false;

  temp = aa.chromosome-bb.chromosome;
  if (temp < 0) return true;
  else if (temp > 0) return false;

  if (aa.position < bb.position) return true;
  else return false;
}


/*
  Function:  compare_data_known, used for sorting results from refFlat file

  Sort based on chromosome and position.
*/
bool compare_data_known(RSW_Known const &aa, RSW_Known const &bb) {
  int temp = aa.chromosome-bb.chromosome;
  if (temp < 0) return true;
  else if (temp > 0) return false;

  return aa.position1 < bb.position1;
}

/*
  Function:   compare_spliceBySupport, used for sorting junctions

  Sorts based on number of supporting reads - consider higher ones first.
*/
bool compare_spliceBySupport(RSW_splice *aa, RSW_splice *bb) {
  if (bb->supported_reads.size() < aa->supported_reads.size()) return true;
  else return false;
}

/*
  Function:   compare_spliceByChromLen, used for sorting junctions

  Sorts by chromosome, length of splice, number of supporting reads - used in deciding
  which splices to print.
*/
bool compare_spliceByChromLen(RSW_splice *aa, RSW_splice *bb) {
  int temp = aa->chromosome-bb->chromosome;
  if (temp < 0) return true;
  else if (temp > 0) return false;

  temp = (aa->positionLarger - aa->positionSmaller) - 
    (bb->positionLarger - bb->positionSmaller);
  if (temp < 0) return true;
  else if (temp > 0) return false;

  if (bb->supported_reads.size() < aa->supported_reads.size()) return true;
  return false;
}

/*
  Function: compare_spliceByChromPos, used for sorting junctions

  Sorts based on chromosome, splice length, position - used in sorting
  before computing supporting reads.
*/
bool compare_spliceByChromPos(RSW_splice *aa, RSW_splice *bb) {
  int temp = aa->chromosome-bb->chromosome;
  if (temp < 0) return true;
  else if (temp > 0) return false;

  temp = (aa->positionLarger - aa->positionSmaller) - 
    (bb->positionLarger - bb->positionSmaller);
  if (temp < 0) return true;
  else if (temp > 0) return false;

  if (aa->positionSmaller < bb->positionSmaller) return true;
  return false;
}


/*
  Function: readOptionsFromFile, reads parameters for program

  Options assumed to be one per line in a predefined order.
*/
void readOptionsFromFile(const char *filename) {
  FILE * fOptions = fopen(filename, "r");
  if (fOptions == NULL) { printf("Error opening file %s\n", filename); exit(0); }

  int pos = 0; int count = 0;
  char *fields[9]; fields[0] = &options[0];
  int ch;
  while ((ch = fgetc(fOptions)) != EOF) {
    if (pos >= MAX_STR_LEN) { printf("Options file %s is more than the max of %i bytes.\n", filename, MAX_STR_LEN); exit(0); }
    if (ch == '\n') {
      options[pos++] = '\0'; count++;
      if (count < 9)
	fields[count] = &options[pos];
      else break;
    }
    else options[pos++] = ch;
  }

  fclose(fOptions);

  if (count < 9) { printf("Not enough lines in options file.\n"); exit(0); }

  sampleDataFile = fields[0];
  maxDistance = atoi(fields[1]);
  sampleLength = atoi(fields[2]);
  refFlatFile = fields[3];
  refFlatBoundaryFile = fields[4];
  minSpliceLength = atoi(fields[5]);
  supportPosTolerance = atoi(fields[6]);
  resultsBaseName = fields[7];
  minSupportingReads = atoi(fields[8]);
}

/*
  Function:  setDefaultOptions, sets some default options if now options
             file is given on the command-line.

  Useful for debugging - set the default options to be whatever you are testing,
  so don't have to type in name of options file each time you run the program.
*/
void setDefaultOptions() {
  sampleDataFile = "RSW_test.txt";
  maxDistance = 40000;
  sampleLength = 33;
  refFlatFile = "refFlat.txt";
  refFlatBoundaryFile = "refFlat.txt.intronBoundary.exonsgaps";
  minSpliceLength = 2;
  supportPosTolerance = 5;
  resultsBaseName = "RSW_tst";
  minSupportingReads = 2;

  printf("Not enough arguments given, using default values.\n");
  printf("Usage is to load options from file: ./splitPairs optionsFile.txt \n");
  printf("And make sure options file has options in order, each on their own line with no extra lines.\n");
}

/*
  Function:  printCurrentOptions, print options to given file pointer.

  Print to file pointer, so can print the options that were used to stdout and/or
  to output files with results.
*/
void printCurrentOptions(FILE *f) {
  fprintf(f, "Running with options...\n"
	  "  file with read data          %s\n"
	  "  max distance between matches %i\n"
	  "  length of samples            %i\n"
	  "  refFlat file                 %s\n"
	  "  refFlat intron boundary file %s\n"
	  "  minimum splice length        %i\n"
	  "  tolerance of difference in position for supporting reads  %i\n"
	  "  base of file name for writing results                     %s\n"
	  "  minimum number of supporting reads                        %i\n\n",
	  sampleDataFile, maxDistance, sampleLength, refFlatFile, refFlatBoundaryFile, minSpliceLength, supportPosTolerance, resultsBaseName, minSupportingReads);

  if (strlen(sampleDataFile) > MAX_STR_LEN - 100) {
    fprintf(f,"Error, filename %s is too long.\n", sampleDataFile); exit(0);
  }
}

// files we will write out to
FILE * fKnown, *fUnknown, *fKnownFull, *fUnknownFull, *fSplitPairs;

// open output files to be ready to write out to them.
void openOutputFiles() {
  sprintf(buff,"%s.results", resultsBaseName);
  fKnown = fopen(buff, "w");
  if (fKnown == NULL) { printf("Error opening file %s for writing.\n", buff); exit(0); }
  printf("Will write summary results that match in known genes to file\n"
	 "   %s\n", buff);

  sprintf(buff,"%s.results.full", resultsBaseName);
  fKnownFull = fopen(buff, "w");
  if (fKnownFull == NULL) { printf("Error opening file %s for writing.\n", buff); exit(0); }
  printf("Will write full results that match in known genes to file\n"
	 "   %s\n", buff);

  sprintf(buff,"%s.results.unknown", resultsBaseName);
  fUnknown = fopen(buff, "w");
  if (fUnknown == NULL) { printf("Error opening file %s for writing.\n", buff); exit(0); }
  printf("Will write summary results that do NOT match in known genes to file\n"
	 "   %s\n", buff);

  sprintf(buff,"%s.results.unknown.full", resultsBaseName);
  fUnknownFull = fopen(buff, "w");
  if (fUnknownFull == NULL) { printf("Error opening file %s for writing.\n", buff); exit(0); }
  printf("Will write full results that do NOT match in known genes to file\n"
	 "   %s\n", buff);

  sprintf(buff,"%s.results.splitPairs", resultsBaseName);
  fSplitPairs = fopen(buff, "w");
  if (fSplitPairs == NULL) { printf("Error opening file %s for writing.\n", buff); exit(0); }
  printf("Will write split pairs to file\n"
	 "   %s\n", buff);
}

/*
  Function:  printSplice, print a given junction to the given opened file

*/
void printSplice(FILE *f, RSW_splice *sp) {
  fprintf(f,
          "%s\t%s\t%li\t%li\t%li--%li\t%s\t%s--%s", 
          sp->geneName, sp->chromosome,
          sp->supported_reads.size(),
          sp->positionLarger-sp->positionSmaller,
          sp->minSmallSupport,sp->maxLargeSupport,
          sp->novel ? "Novel" : "*",
	  sp->sequenceSmaller, sp->sequenceLarger 

          );
}


/*
  Function:   updateSpliceSupport - compare the two junctions and put
              each other into each others list of supporting reads
              if they support each other and haven't been put into
              each other's list of supporting reads yet.
*/
void updateSpliceSupport(RSW_splice *splice1, RSW_splice *splice2) {

  for(int i=0; i < 2; i++) {
    RSW_splice * sp1, *sp2;
    if (i == 0) {sp1 = splice1; sp2 = splice2;}
    else {sp1 = splice2; sp2 = splice1;}

    //#pragma omp critical
    {
      if (sp1->supported_reads.size() <= 0 ||
	  sp2->positionSmaller < sp1->minSmallSupport)
	sp1->minSmallSupport = sp2->positionSmaller;
      if (sp1->supported_reads.size() <= 0 ||
	  sp2->positionLarger > sp1->maxLargeSupport)
	sp1->maxLargeSupport = sp2->positionLarger;
      sp1->supported_reads.insert({sp2->id, sp2});
    }
  }
}

/*
  Function:  printStats, prints statistics gathered so far to the
             opened file - useful for debugging to see some partial
             information as each phase of the program finishes.
*/
void printStats(FILE *f) {
  FILE *fStatus = fopen("/proc/self/status","r");
  char s[1000], mem[20]="", units[20]="";
  while (fscanf(fStatus,"%999s",s) == 1) {
    if (strcmp(s,"VmRSS:") == 0) {
      fscanf(fStatus,"%19s %19s",mem, units);
      break;
    }
  }
  fclose(fStatus);

  endTime = time(NULL);
  fprintf(f, "Finished processing data, results written to files.\n");
  fprintf(f, "Number of entries in data file:             %li\n", data.size());
  fprintf(f, "Number of different reads:                  %i\n", numDifferentReads);
  fprintf(f, "Number of entries in refFlat file:          %li\n", data_known.size());
  fprintf(f, "Number of entries in refFlat boundary file: %li\n", data_boundaries.size());
  fprintf(f, "Number of matches:                          %li\n", data_splice.size());
  fprintf(f, "String table size:                          %li\n", stringTable.size());
  fprintf(f, "VmRSS, memory resident set size:            %s %s\n", mem, units);
  fprintf(f, "Total time to process:                      %li seconds\n", endTime-beginTime);
  fprintf(f, "\n");
}

int main(int argc, char *argv[]) {
  setpriority(0, 0, 20); // so other processes get priority over this one

  beginTime = time(NULL);

  // read options, from file or default options
  if (argc > 1) 
    readOptionsFromFile(argv[1]);
  else 
    setDefaultOptions();

  // write out options to all output files and stdout

  openOutputFiles();

  printCurrentOptions(stdout);
  printCurrentOptions(fKnown);
  printCurrentOptions(fKnownFull);
  printCurrentOptions(fUnknown);
  printCurrentOptions(fUnknownFull);
  printCurrentOptions(fSplitPairs);


  // read from refFlat file into data_known array, 
  read_knownGene(refFlatFile);
  sort(data_known.begin(), data_known.end(), compare_data_known);
  printf("Done reading/sorting refFlat, total time elapsed %li seconds\n", time(NULL)-beginTime);
  printStats(stdout);

  // read from refFlat boundary file into data_boundaries array, 
  read_boundaries(refFlatBoundaryFile);
  printf("Done reading refFlat intron/exon boundaries, total time elapsed %li seconds\n", time(NULL)-beginTime);
  printStats(stdout);

  // read the read data
  read_data(sampleDataFile);
  printf("Done reading read data, total time elapsed %li seconds\n", time(NULL)-beginTime);
  printStats(stdout);

  // sort the read data
  sort(data.begin(), data.end(), compare_dataToSort);
  printf("Done sorting read data, total time elapsed %li seconds\n", time(NULL)-beginTime);
  printStats(stdout);

  // look at all pairs of read segments, looking for matches
  numDifferentReads = 0;
  const int inputSize = data.size();
  //#pragma omp parallel for
  for(int left=0; left < inputSize; left++) {
    int right;
    if (left == 0 || //strcmp(data[left].id,data[left-1].id) != 0) 
	(data[left].id  != data[left-1].id) ) {
      //#pragma omp atomic
      numDifferentReads++;
    }

    for(right=left+1; right < inputSize; right++) {
      // read segments are ordered by id/chromosome/strand, so if 
      // there isn't a match we can skip the rest of the read segments
      // for "right", and go to the next iteration of the "left" loop
      if (data[left].id != data[right].id ||
	  data[left].direction != data[right].direction ||
	  data[left].chromosome != data[right].chromosome) break;

      if (data[right].position - data[left].position > maxDistance) break;

      // want it to be from two sides of the same segment
      if (data[left].side == data[right].side) continue;

      // total length should be correct
      if (data[left].length + data[right].length != sampleLength) continue;

      // calculate the end of the segments, since what is given
      // in the data is the beginning of the segments
      int endSmaller, endLarger; // splice is between endSmaller and endLarger
      int first, second;
      if (data[left].side == 'L' && data[left].direction == '+' ||
	  data[left].side == 'R' && data[left].direction == '-') { first = left; second = right; }
      else { first = right; second = left; }
      endSmaller = data[first].position + data[first].length;
      endLarger = data[second].position;

      // splice length, and check that it is within specified bounds
      int spliceLength = endLarger - endSmaller;
      if (spliceLength > maxDistance) continue;
      if (spliceLength < minSpliceLength) continue;

      // check if we already have this splice from this read...
      int i;
      for(i=0; i < data_splice.size(); i++) {
	if (data_splice[i]->positionSmaller != endSmaller ||
	    data_splice[i]->positionLarger != endLarger ||
	    data_splice[i]->id != data[left].id ||
	    data_splice[i]->chromosome != data[left].chromosome) 
	  continue;
	break;
      }
      // if already have this exact splice for this chromosome from this read, don't include it again.
      if (i < data_splice.size()) continue;

      // note: could print this match here, step 5 done.

      // look for this in the known gene...
      int k; int foundInGene = 0;
      for(k=0; k < data_known.size(); k++) {
	if ((((data[left].chromosome == data_known[k].chromosome) &&
	      (data[left].position >= data_known[k].position1 &&
	       data[right].position >= data_known[k].position1) &&
	      (data[left].position <= data_known[k].position2 &&
	       data[right].position <= data_known[k].position2)))) {
	  if (foundInGene) ;//printf("DUPLICATE_"); // duplicate
	  foundInGene = 1;
	  break; // just cut off search, don't look for duplicates
	}
      }
      int geneIndex = k;

      // make a new splice record and put into vector of splices
      RSW_splice *sp = new RSW_splice;
      sp;
      if (foundInGene) {
	sp->geneName = data_known[geneIndex].id1;
	sp->geneUnknown = 0;
      }	 
      else {
	sp->geneName = unfound_string;
	sp->geneUnknown = 1;
      }
      sp->id = data[left].id;
      sp->chromosome = data[left].chromosome;
      sp->direction = data[left].direction;
      sp->positionSmaller = sp->minSmallSupport = endSmaller;
      sp->positionLarger = sp->maxLargeSupport = endLarger;
      sp->supported_reads.insert({sp->id, sp});

      if (data[left].position < data[right].position) {
        sp->sequenceSmaller = data[left].sequence;
        sp->sequenceLarger = data[right].sequence;
      }
      else {
        sp->sequenceSmaller = data[right].sequence;
        sp->sequenceLarger = data[left].sequence;
      }


      //#pragma omp critical
      {
	data_splice.push_back(sp);
      }
    }
  }

  printf("Done finding matched pairs, total time elapsed %li seconds\n", time(NULL)-beginTime);
  printStats(stdout);

  // sort splices by chromosome and length
  sort(data_splice.begin(), data_splice.end(), compare_spliceByChromPos);

  printf("Done sorting matched pairs, total time elapsed %li seconds\n", time(NULL)-beginTime);
  printStats(stdout);


  // compute supporting reads

  //#pragma omp parallel for
  for(int sp1=0; sp1 < data_splice.size(); sp1++) {
    int sp2;
    for(sp2=sp1; sp2 < data_splice.size(); sp2++) {
      // see if these reads support each other

      // because splices are sorted by chromosome and length, if these don't match then 
      // can break sp2 loop.
      if (data_splice[sp1]->chromosome != data_splice[sp2]->chromosome) break;
      if (abs(data_splice[sp1]->positionLarger - data_splice[sp1]->positionSmaller) != 
	  abs(data_splice[sp2]->positionLarger - data_splice[sp2]->positionSmaller)) break; // splice length must be the same
      if (data_splice[sp1]->geneName != data_splice[sp2]->geneName) break; // ???

      // difference in where the splice is, 0 for exact match.  if doesn't match, 
      // continue to next sp2.
      if (abs(data_splice[sp1]->positionSmaller-data_splice[sp2]->positionSmaller) > supportPosTolerance) 
	continue;
	//break;
      if (data_splice[sp1]->id == data_splice[sp2]->id) continue; // read can't support itself

      // they are matches for each other, so add each to each other's list, as 
      // long as not a match already there for the same read id.
      updateSpliceSupport(data_splice[sp1], data_splice[sp2]);
    }
  }

  printf("Done computing supporting reads, total time elapsed %li seconds\n", time(NULL)-beginTime);
  printStats(stdout);

  // sort splices by chromosome and length, and # supporting reads
  sort(data_splice.begin(), data_splice.end(), compare_spliceByChromLen);

  printf("Done sorting matched pairs again, total time elapsed %li seconds\n", time(NULL)-beginTime);
  printStats(stdout);

  printf("Done filtering matched pairs, total time elapsed %li seconds\n", time(NULL)-beginTime);
  printStats(stdout);

  endTime = time(NULL);


  printf("Done calculating/filtering supporting reads, total time elapsed %li seconds\n", time(NULL)-beginTime);
  printStats(stdout);

  // print the splice results.  start with the statistics
  printStats(stdout);
  printStats(fKnown);
  printStats(fKnownFull);
  printStats(fUnknown);
  printStats(fUnknownFull);
  printStats(fSplitPairs);

  int k;


  // save all the splices, and reads that support each one into .splitPairs file - this is
  // the full results, with many duplicates of splices.  The file also can be HUGE, so
  // generally this file will often get deleted unless needed for debugging.
  fprintf(fSplitPairs, "Line No\tId\tGene\tChr\tStrand\tSplice region\tSupporting reads\tSupporting splice range\tBracketed sequence\n"); 
  for(k=0; k < data_splice.size(); k++) {
    fprintf(fSplitPairs, "%s\t%s\t%s\t%c\t%li-%li\t%li\t%li-%li\t%s-%s\n", 
            data_splice[k]->id, data_splice[k]->geneName,
            data_splice[k]->chromosome, data_splice[k]->direction,
            data_splice[k]->positionSmaller,data_splice[k]->positionLarger,
            data_splice[k]->supported_reads.size(),
            data_splice[k]->minSmallSupport,data_splice[k]->maxLargeSupport,
	    data_splice[k]->sequenceSmaller,data_splice[k]->sequenceLarger
            );

    for(auto &x: data_splice[k]->supported_reads)
      fprintf(fSplitPairs, "%s %li--%li %s--%s, ", x.second->id, x.second->positionSmaller, x.second->positionLarger, x.second->sequenceSmaller, x.second->sequenceLarger); 
    fprintf(fSplitPairs, "\n");
  }

  // save the tabulated results.  the fKnown file is the only one normally looked at.
  // fKnownFull has the same information as fKnown, but also has the list of supporting reads for each junction.
  // fUnknown and fUnknownFull are for junctions not within genes.
  fprintf(fKnown, "GeneName\tChromosome\t# supporting reads\tsplice length\trange of supporting reads\tNovel or not (*)\tBracketed sequence\n"); 
  fprintf(fUnknown, "GeneName\tChromosome\t# supporting reads\tsplice length\trange of supporting reads\tNovel or not (*)\tBracketed sequence\n");
  fprintf(fKnownFull, "GeneName\tChromosome\t# supporting reads\tsplice length\trange of supporting reads\tNovel or not (*)\tBracketed sequence\n");
  fprintf(fUnknownFull, "GeneName\tChromosome\t# supporting reads\tsplice length\trange of supporting reads\tNovel or not (*)\tBracketed sequence\n");
  FILE * f, *fFull;
  for(k=0; k < data_splice.size(); k++) {
    // if read already reported, don't report it again.
    if (readIdsReported.find(data_splice[k]->id) != readIdsReported.end()) {continue;}

    // if doesn't have enough supporting reads, don't print it
    if (data_splice[k]->supported_reads.size() < minSupportingReads) continue;

    // will be reported, so log that.
    readIdsReported.insert(data_splice[k]->id);

    // check if novel or not.
    int j;
    data_splice[k]->novel = true;
    for(j=0; j < data_boundaries.size(); j++) {
      if (data_splice[k]->positionLarger-data_splice[k]->positionSmaller != 
	  data_boundaries[j].length) continue;
      if (abs(data_splice[k]->minSmallSupport-data_boundaries[j].position1) <= supportPosTolerance &&
	  abs(data_splice[k]->maxLargeSupport-data_boundaries[j].position2) <= supportPosTolerance)
	break;
    }
    if (j < data_boundaries.size())  // if found in the boundaries data, not novel
      data_splice[k]->novel = false;

    // going into known file or unknown
    if (data_splice[k]->geneUnknown) {
      f = fUnknown; fFull = fUnknownFull;
    }
    else {
      f = fKnown; fFull = fKnownFull;
    }

    // print out results
    printSplice(f, data_splice[k]);
    fprintf(f, "\n");

    // print out full results, that includes id's of supporting reads
    printSplice(fFull, data_splice[k]);

    // print out full results, that includes id's of supporting reads
    fprintf(fFull, "\t");
    for(auto& x: data_splice[k]->supported_reads) {
      fprintf(fFull,", %s", x.second->id);
      readIdsReported.insert(x.second->id);
    }
    fprintf(fFull, "\n");
  }

  fclose(fKnown); fclose(fKnownFull); fclose(fUnknown); fclose(fUnknownFull); fclose(fSplitPairs);

  printf("Done saving results, total time elapsed %li seconds\n", time(NULL)-beginTime);
  printStats(stdout);

  // free memory.  good to do so we can run a memory checker and verify
  // we don't have any memory leaks.
  while (data_splice.size() > 0) {
    delete data_splice.back(); data_splice.pop_back();
  }
  
  return 0;
}
