#include <config.h>

#include <string.h>
#include <stdlib.h>

#include "htslib/sam.h"
#include "htslib/bgzf.h"

int main(int argc, char **argv)
{
	bam_hdr_t *header;
	bam1_t *aln = bam_init1();

	int n_threads = 0;
	char *samin = NULL;
	char *samout = NULL;

	int i;
	for(i = 1; i < argc; i++) {
		if(strcmp(argv[i], "-t") == 0 && argc - i > 2) {
			n_threads = atoi(argv[++i]);
			continue;
		}
		if(samin == NULL) {
			samin = argv[i];
			continue;
		}
		if(samout == NULL){
			samout = argv[i];
			break;
		}
	}
	if(samin == NULL){
		fprintf(stderr, "Usage:\n\t%s [-t threads] SAMIN [SAMOUT]\n", argv[0]);
		return;
	}
	if(samout == NULL) samout = "-";

	// OPEN SAMFIEL FOR READ
	samFile *in = sam_open(samin, "r");
	if(in == NULL) {
		perror(samin);
		exit(1);
	}

	// OPEN SAMFILE FOR WRITE
	samFile *out = sam_open(samout, "wb");
	if(out == NULL) {
		perror(samout);
		exit(1);
	}

	// PARALLEL
	if(n_threads > 1) {
		hts_set_threads(in, n_threads);
		bgzf_set_cache_size(in->fp.bgzf, BGZF_MAX_BLOCK_SIZE * n_threads * 256);
		hts_set_threads(out, n_threads);
	}

	// PROCESS
	if((header = sam_hdr_read(in)) == NULL) {
		fprintf(stderr, "fail to read the header from \"%s\".\n", argv[i]);
		exit(2);
	}
	if(sam_hdr_write(out, header) < 0) {
		fprintf(stderr, "fail to write header to \"%s\".\n", samout);
		exit(3);
	}
	
	char *qname, *umi;
	size_t l_umi;
	int n_colon = 0;
	char tag[2] = "RX";
	while(sam_read1(in, header, aln) >= 0) {
		qname = strdup(bam_get_qname(aln));
		for(i = 0; i < strlen(qname); i++) {
			if(qname[i] == ':')
				n_colon++;
			if(n_colon == 7) {
				n_colon = 0;
				l_umi = strlen(qname+i+1);
				umi = malloc(l_umi + 1);
				strncpy(umi, qname+i+1, l_umi);
				break;
			}
		}
		if(bam_aux_update_str(aln, tag, l_umi+1, umi) < 0) {
			fprintf(stderr, "fail to add \"RX\" tag to \"%s\".\n", samout);
			exit(3);
		}
		uint8_t *s = bam_aux_get(aln, tag);
		if(sam_write1(out, header, aln) < 0) {
			fprintf(stderr, "Fail to Write \"%s\" to \"%s\".\n", qname, samout);
			exit(3);
		}
		free(qname);
		free(umi);
	}
	
	sam_close(in);
	sam_close(out);
	bam_destroy1(aln);

	return 0;
}
			
