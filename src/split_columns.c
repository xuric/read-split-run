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

void parse(char *filename) {
    char *line = 0;
    char *readName,*tmp;
    size_t n;
    int pos;
    FILE *in, *out;
    char *outname;

    
    if (!(in = fopen(filename,"r"))) {
        fprintf(stderr, "Could not open %s for reading", filename);
        perror("");
        return;
    }
    
    outname = (char *)malloc(sizeof(char) * (strlen(filename)+strlen(SUFFIX)));
    memset(outname, 0, sizeof(char) * (strlen(filename)+strlen(SUFFIX)));
    sprintf(outname, "%s%s",filename,SUFFIX);
    outname[strlen(filename)+strlen(SUFFIX)] = 0;
    if (!(out = fopen(outname, "w"))) {
        fprintf(stderr,"Could not open %s for writing", outname);
        perror("");
        close(in);
        return;
    }

    tmp = 0;
    while((!feof(in)) && (getline(&line,&n,in) != -1)) {
        if (tmp) { free(tmp); tmp=0; }
        tmp = strdup(line);
        char *found = strchr(line, '\t');
        pos = (found - line);
        split_field(tmp, pos);
        fputs(tmp,out);
    }

    fclose(in);
    fclose(out);
}


int main(int argc, char *argv[]) {
    int i;

    if (argc < 2) {
        fprintf(stderr, "usage -- %s <file to split> [additonal files...]\n",argv[0]);
        return 1;
    }

    for (i = 1;i < argc;i++)
        parse(argv[i]);

    return 0;
}
