// old file, apparently no longer used...

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>


#define SUFFIX ".split1stcolumn"

//D2FC08P1:143:D0KHCACXX:6:1101:1563:1953 2:N:0:-R-68/2	

void split_field(char *readName, int len) {
    int i, offset;

    if (!readName || strlen(readName) <= len) {
        fprintf(stderr,"No field supplied to split_field. aborting\n");
        exit(1);
    }
//actual output
    if (readName[len-2] == '/' && (readName[len-1] == '1' || readName[len-1] == '2')) { //paired
        //before:
        //D2FC08P1:143:D0KHCACXX:6:1101:1563:1953 2:N:0:-R-68/2	
        //let's do some horse-trading...
        //put the /1 or /2 at the beginning of the substring...
        offset = 7;
        char t1,t2;
        t1 = readName[len - 2];
        t2 = readName[len - 1];
        //shift everything right
        for (i = (len - 3);i > (len - offset);i--) {
             readName[i+2] = readName[i];
        }
        readName[i] = t1;
        readName[i+1] = t2;
    }
    //after, OR, the /2 is irrelevant
    //D2FC08P1:143:D0KHCACXX:6:1101:1563:1953 2:N:0:/2-R-68
    i = len - 5;
    readName[i] = '\t';
    readName[i+2] = '\t';
}

void parse(FILE *in, FILE *out) {
    char *line = 0;
    char *readName,*tmp;
    size_t n;
    int pos;

    tmp = 0;
    while((!feof(in)) && (getline(&line,&n,in) != -1)) {
        if (tmp) { free(tmp); tmp=0; }
        tmp = strdup(line);
        char *found = strchr(line, '\t');
        pos = (found - line);
        split_field(tmp, pos);
        fputs(tmp,out);
    }
}


int main(int argc, char *argv[]) {
    FILE *fr, *fw;
    char *outname;
    size_t newLen, as, ss;

    if (argc < 2) {
        fprintf(stderr, "usage -- %s <file to split>\n",argv[0]);
        return 1;
    }

    if (!(fr = fopen(argv[1],"r"))) {
        fprintf(stderr, "Could not open %s for reading", argv[1]);
        perror("");
        return 1;
    }
    as = strlen(argv[1]) *sizeof(char);
    ss = strlen(SUFFIX) *sizeof(char);
    outname = (char *)malloc(as+ss + 1); // bug fix 6 jan 2016, added +1 to include space for terminating NULL character
    memset(outname, 0, as+ss+1); // bug fix 6 jan 2016, added +1 to include space for terminating NULL character
    strcpy(outname, argv[1]);
    strcat(outname, SUFFIX);
    outname[as+ss] = '\0';
    if (!(fw = fopen(outname, "w"))) {
        fprintf(stderr,"Could not open %s for writing", outname);
        perror("");
        close(fr);
        return 1;
    }

    parse(fr,fw);
    fclose(fr);
    fclose(fw);
    return 0;
}
