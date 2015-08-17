/*
  To compile: g++ splitPairs.cpp -o sp -O4 -std=c++11
  To run:     ./sp options.txt
  where options.txt is an options file.  If the program is run
  with no command-line arguments it by default processes 
  RSW_test.txt with the same parameters as the scripts downloaded
  from Yongsheng Bai's website.  When you run the program it
  prints which files output is written to.

Modification history...
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


int maxDistance;
int sampleLength;
int minSpliceLength;
int supportPosTolerance;
char const *sampleDataFile;
char const *refFlatFile;
char const *refFlatBoundaryFile;
char const *resultsBaseName;
int minSupportingReads;

char buff[MAX_STR_LEN];
char options[MAX_STR_LEN];

time_t beginTime, endTime;

unordered_set<string> stringTable;
unordered_set<const char *> readIdsReported;

vector<struct RSW> data;
vector<struct RSW_Known> data_known;
vector<struct RSW_Boundaries> data_boundaries;
const char unfound_string[100] = "UNFOUND_";
vector<RSW_splice *> data_splice;
int numDifferentReads;

int compute_hash(RSW *d) {
  int h = 1, i;
  for(i=0; d->id[i] != '\0'; i++) h = h * d->id[i] % SHRT_MAX;
  for(i=0; d->chromosome[i] != '\0'; i++) h = h * d->chromosome[i] % SHRT_MAX;
  h = h * d->direction % SHRT_MAX;
  return h;
}


char sLine[MAX_LINE+1];
char temp[MAX_LINE+1];

void read_data(const char *filename) {
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
    char *tempA = strtok(sLine, "\t");
    i=0;
    while (tempA != NULL) {
      string temp = tempA; temp.shrink_to_fit();
      pair<unordered_set<string>::iterator,bool> result;
      switch (i) {
      case 0:
	result = stringTable.insert(temp);
	r->id = (result.first)->c_str();

	//r->id = (char *)malloc((strlen(temp)+1)*sizeof(char));
	//strcpy(r->id, temp);
	break;
      case 1:
	r->side = temp[0];
	break;
      case 2:
	r->length = atoi(temp.c_str());
	break;
      case 3:
	r->direction = temp[0];
	break;
      case 4:
	result = stringTable.insert(temp);
	r->chromosome = (result.first)->c_str();
	//r->chromosome = (char *) malloc(sizeof(char) * (strlen(temp)+1));
	//strcpy(r->chromosome, temp);
	break;
      case 5:
	r->position = atol(temp.c_str());
	break;
      case 8:
	r->count = atoi(temp.c_str());
	break;
      }
      i++;
      tempA = strtok(NULL,"\t");
      if (tempA == NULL) break;
    }
    if (i < 9) {delete r; break;}

    r->hash = compute_hash(r);
    data.push_back(*r);
    delete r;

    //if (data.size() >= 15000000) break; // cut off early, for debugging.
  }

  fclose(f);
}

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

	//rk->id1 = (char *)malloc((strlen(temp)+1)*sizeof(char));
	//strcpy(rk->id1, temp);
	break;
      case 1:
	result = stringTable.insert(temp);
	rk->id2 = (result.first)->c_str();

	//rk->id2 = (char *)malloc((strlen(temp)+1)*sizeof(char));
	//strcpy(rk->id2, temp);
	break;
      case 2:
	result = stringTable.insert(temp);
	rk->chromosome = (result.first)->c_str();

	//rk->chromosome = (char *) malloc(sizeof(char) * (strlen(temp)+1));
	//strcpy(rk->chromosome, temp);
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

	//rk->id1 = (char *)malloc((strlen(temp)+1)*sizeof(char));
	//strcpy(rk->id1, temp);
	break;
      case 1:
	result = stringTable.insert(temp);
	rk->id2 = (result.first)->c_str();

	//rk->id2 = (char *)malloc((strlen(temp)+1)*sizeof(char));
	//strcpy(rk->id2, temp);
	break;
      case 2:
	result = stringTable.insert(temp);
	rk->chromosome = (result.first)->c_str();

	//rk->chromosome = (char *) malloc(sizeof(char) * (strlen(temp)+1));
	//strcpy(rk->chromosome, temp);
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

bool compare_data(RSW const &aa, RSW const &bb) {
  //int temp = strcmp(aa.id, bb.id);
  int temp = aa.id-bb.id;
  if (temp < 0) return true;
  else if (temp > 0) return false;

  if (aa.direction < bb.direction) 
    return true;
  else if (bb.direction < aa.direction)
    return false;

  //temp = strcmp(aa.chromosome, bb.chromosome);
  temp = aa.chromosome-bb.chromosome;
  if (temp < 0) return true;
  else if (temp > 0) return false;

  return false;
}

bool compare_dataToSort(RSW const &aa, RSW const &bb) {
  //int temp = strcmp(aa.id, bb.id);
  int temp = aa.id-bb.id;
  if (temp < 0) return true;
  else if (temp > 0) return false;

  if (aa.direction < bb.direction) 
    return true;
  else if (bb.direction < aa.direction)
    return false;

  //temp = strcmp(aa.chromosome, bb.chromosome);
  temp = aa.chromosome-bb.chromosome;
  if (temp < 0) return true;
  else if (temp > 0) return false;

  if (aa.position < bb.position) return true;
  else return false;
}


bool compare_data_known(RSW_Known const &aa, RSW_Known const &bb) {
  //int temp = strcmp(aa.chromosome, bb.chromosome);
  int temp = aa.chromosome-bb.chromosome;
  if (temp < 0) return true;
  else if (temp > 0) return false;

  return aa.position1 < bb.position1;
}

bool compare_spliceBySupport(RSW_splice *aa, RSW_splice *bb) {
  if (bb->supported_reads.size() < aa->supported_reads.size()) return true;
  else return false;
}

// sort by chromosome, length fo splice, # supporting reads - used in deciding which splices to print
bool compare_spliceByChromLen(RSW_splice *aa, RSW_splice *bb) {
  //int temp = strcmp(aa->chromosome, bb->chromosome);
  int temp = aa->chromosome-bb->chromosome;
  if (temp < 0) return true;
  else if (temp > 0) return false;

  temp = (aa->positionLarger - aa->positionSmaller) - 
    (bb->positionLarger - bb->positionSmaller);
  if (temp < 0) return true;
  else if (temp > 0) return false;

  if (bb->supported_reads.size() < aa->supported_reads.size()) return true;
  //if (aa->positionSmaller < bb->positionSmaller) return true;
  return false;
}

// sort by chromosome, length of splice, position - used for calculating supporting reads
bool compare_spliceByChromPos(RSW_splice *aa, RSW_splice *bb) {
  //int temp = strcmp(aa->chromosome, bb->chromosome);
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

FILE * fKnown, *fUnknown, *fKnownFull, *fUnknownFull, *fSplitPairs;

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

void printSplice(FILE *f, RSW_splice *sp) {
  fprintf(f,
	  "%s\t%s\t%li\t%li\t%li--%li\t%s",
	  sp->geneName, sp->chromosome, 
	  sp->supported_reads.size(),
	  sp->positionLarger-sp->positionSmaller, 
	  sp->minSmallSupport,sp->maxLargeSupport,
	  //sp->positionSmaller, sp->positionLarger, 
	  sp->novel ? "Novel" : "*"
	  );
}


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
  fprintf(f, "Number of entries in data file:             %i\n", data.size());
  fprintf(f, "Number of different reads:                  %i\n", numDifferentReads);
  fprintf(f, "Number of entries in refFlat file:          %i\n", data_known.size());
  fprintf(f, "Number of entries in refFlat boundary file: %i\n", data_boundaries.size());
  fprintf(f, "Number of matches:                          %i\n", data_splice.size());
  fprintf(f, "String table size:                          %i\n", stringTable.size());
  fprintf(f, "VmRSS, memory resident set size:            %s %s\n", mem, units);
  fprintf(f, "Total time to process:                      %i seconds\n", endTime-beginTime);
  fprintf(f, "\n");
}

int main(int argc, char *argv[]) {
  setpriority(0, 0, 20); // so other processes get priority over this one

  beginTime = time(NULL);

  if (argc > 1) 
    readOptionsFromFile(argv[1]);
  else 
    setDefaultOptions();

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
  printf("Done reading/sorting refFlat, total time elapsed %i seconds\n", time(NULL)-beginTime);
  printStats(stdout);

  // read from refFlat boundary file into data_boundaries array, 
  read_boundaries(refFlatBoundaryFile);
  printf("Done reading refFlat intron/exon boundaries, total time elapsed %i seconds\n", time(NULL)-beginTime);
  printStats(stdout);

  // read the read data
  read_data(sampleDataFile);
  printf("Done reading read data, total time elapsed %i seconds\n", time(NULL)-beginTime);
  printStats(stdout);

  // sort the read data
  sort(data.begin(), data.end(), compare_dataToSort);
  printf("Done sorting read data, total time elapsed %i seconds\n", time(NULL)-beginTime);
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
	    //strcmp(data_splice[i]->id, data[left].id) != 0 ||
	    data_splice[i]->chromosome != data[left].chromosome) 
	  //strcmp(data_splice[i]->chromosome, data[left].chromosome) != 0) 
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
	      //if ((((strcmp(data[left].chromosome, data_known[k].chromosome) == 0) &&
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
      //if (!found) printf("UNFOUND_\t"); // not found in any known gene.
      //else printf("%s\t", data_known[k].id1);

      // add this splice match to the data_splice vector.
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
      //sp->reported = false;
      sp->supported_reads.insert({sp->id, sp});
      //#pragma omp critical
      {
	data_splice.push_back(sp);
      }
    }
  }

  printf("Done finding matched pairs, total time elapsed %i seconds\n", time(NULL)-beginTime);
  printStats(stdout);

  //fclose(fKnown); fclose(fKnownFull); fclose(fUnknown); fclose(fUnknownFull); fclose(fSplitPairs);
  //exit(0);


  // sort splices by chromosome and length
  sort(data_splice.begin(), data_splice.end(), compare_spliceByChromPos);

  printf("Done sorting matched pairs, total time elapsed %i seconds\n", time(NULL)-beginTime);
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
      //if (strcmp(data_splice[sp1]->chromosome, data_splice[sp2]->chromosome) != 0) break;
      if (abs(data_splice[sp1]->positionLarger - data_splice[sp1]->positionSmaller) != 
	  abs(data_splice[sp2]->positionLarger - data_splice[sp2]->positionSmaller)) break; // splice length must be the same
      if (data_splice[sp1]->geneName != data_splice[sp2]->geneName) break; // ???
      //if (strcmp(data_splice[sp1]->geneName, data_splice[sp2]->geneName) != 0) break;

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

  printf("Done computing supporting reads, total time elapsed %i seconds\n", time(NULL)-beginTime);
  printStats(stdout);

  // sort splices by chromosome and length, and # supporting reads
  sort(data_splice.begin(), data_splice.end(), compare_spliceByChromLen);

  printf("Done sorting matched pairs again, total time elapsed %i seconds\n", time(NULL)-beginTime);
  printStats(stdout);

  // set anything that is included within another splice region, or the same splice region,
  // as already reported.  if saying something is reported, take one with larger supporting reads as
  // one to print.
  /*  for(sp1=0; sp1 < data_splice.size(); sp1++) {
    if (data_splice[sp1]->reported) continue; // already reported sp1, so go to next.

    for(sp2=sp1+1; sp2 < data_splice.size(); sp2++) {
      // sorted by chromosome and length, if those don't match then not supporting each other, so 
      // go to next sp1.
      if (data_splice[sp1]->chromosome != data_splice[sp2]->chromosome) break;
      //if (strcmp(data_splice[sp1]->chromosome, data_splice[sp2]->chromosome) != 0) break;
      if (abs(data_splice[sp1]->positionLarger - data_splice[sp1]->positionSmaller) != 
	  abs(data_splice[sp2]->positionLarger - data_splice[sp2]->positionSmaller)) break; // splice length must be the sam

      if (data_splice[sp2]->reported) continue; // already not reporting sp2, so go to next.

      int m1 = data_splice[sp1]->minSmallSupport - data_splice[sp2]->minSmallSupport;
      int m2 = data_splice[sp1]->maxLargeSupport - data_splice[sp2]->maxLargeSupport;
      if (m1 <= 0 && m2 >= 0) {// splice sp2 region included in splice sp1 region
	// and # reads of sp2 should be <= # reads of sp1 since sorted by that ordering
	// so choose not to report sp2.
	data_splice[sp2]->reported = true;
      }
    }
    }*/

  printf("Done filtering matched pairs, total time elapsed %i seconds\n", time(NULL)-beginTime);
  printStats(stdout);

  endTime = time(NULL);


  printf("Done calculating/filtering supporting reads, total time elapsed %i seconds\n", time(NULL)-beginTime);
  printStats(stdout);

  // print the splice results.  start with the statistics
  printStats(stdout);
  printStats(fKnown);
  printStats(fKnownFull);
  printStats(fUnknown);
  printStats(fUnknownFull);
  printStats(fSplitPairs);

  int k;


  // save all the splices, and reads that support each one.
  fprintf(fSplitPairs, "Line No\tId\tGene\tChr\tStrand\tSplice region\tSupporting reads\tSupporting splice range\n");
  for(k=0; k < data_splice.size(); k++) {
    fprintf(fSplitPairs, "%s\t%s\t%s\t%c\t%i-%i\t%i\t%i-%i\n",
	    data_splice[k]->id, data_splice[k]->geneName, 
	    data_splice[k]->chromosome, data_splice[k]->direction,
	    data_splice[k]->positionSmaller,data_splice[k]->positionLarger,
	    data_splice[k]->supported_reads.size(),
	    data_splice[k]->minSmallSupport,data_splice[k]->maxLargeSupport);

    for(auto &x: data_splice[k]->supported_reads)
      fprintf(fSplitPairs, "%s, ", x.second->id);
    fprintf(fSplitPairs, "\n");
  }

  // save the tabulated results.
  fprintf(fKnown, "GeneName\tChromosome\t# supporting reads\tsplice length\trange of supporting reads\tNovel or not (*)\n");
  fprintf(fUnknown, "GeneName\tChromosome\t# supporting reads\tsplice length\trange of supporting reads\tNovel or not (*)\n");
  fprintf(fKnownFull, "GeneName\tChromosome\t# supporting reads\tsplice length\trange of supporting reads\tNovel or not (*)\n");
  fprintf(fUnknownFull, "GeneName\tChromosome\t# supporting reads\tsplice length\trange of supporting reads\tNovel or not (*)\n");
  FILE * f, *fFull;
  for(k=0; k < data_splice.size(); k++) {
    if (readIdsReported.find(data_splice[k]->id) != readIdsReported.end()) {continue;}

    /*    if (data_splice[k]->reported) {continue;}  // if already reported as supporting another read, don't print it.
	  data_splice[k]->reported = true;*/

    // if doesn't have enough supporting reads, don't print it
    if (data_splice[k]->supported_reads.size() < minSupportingReads) continue;

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
    if (j < data_boundaries.size()) 
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

    fprintf(fFull, "\t");
    for(auto& x: data_splice[k]->supported_reads) {
      fprintf(fFull,", %s", x.second->id);
      //x.second->reported = true; // already reported this splice as supporting another read that has higher support, so don't report it again
      readIdsReported.insert(x.second->id);
    }
    fprintf(fFull, "\n");
  }

  fclose(fKnown); fclose(fKnownFull); fclose(fUnknown); fclose(fUnknownFull); fclose(fSplitPairs);

  printf("Done saving results, total time elapsed %i seconds\n", time(NULL)-beginTime);
  printStats(stdout);

  // free memory.  good to do so we can run a memory checker and verify
  // we don't have any memory leaks.
  while (data_splice.size() > 0) {
    delete data_splice.back(); data_splice.pop_back();
  }
  
  return 0;
}
