#include <iostream>
#include <fstream>
#include <cstring>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "MersenneTwister.h"

using namespace std;

#define MAXSITE 500000

int main (int argc, char* argv[]) {

	if (argc<5 || argc >6) {
		cerr << "\nratchblock: simple PAUP parsimony ratchet program\n\n";
		cerr << "usage:\n";
		cerr << "  ratchblock <outfile> <nchar> <# replicates> <% upweight> <multiplier> <blockf>\n\n";
		cerr << "  outfile      = name used for output files\n";
		cerr << "    (ouput files include the ratchet block, log, and treefiles)\n";
		cerr << "  nchar        = number of characters in the data file\n";
		cerr << "  # replicates = number of ratchet searches\n";
		cerr << "  % upweight   = percentage (1-100) of sites upweighted\n";
		cerr << "  multiplier   = maxtrees multiplier\n";
		cerr << "    (final maxtrees value for swapping will be # replicates * multiplier)\n";
		cerr << "  blockf       = name of file with additional commands\n";
		cerr << "    NOTE: blockf is optional; if used, it should contain a set of\n";
		cerr << "          commands executed at the beginning of the analysis\n";
		return 0;
	}
	bool useblockf=false;
	char line[2048];
	if (argc==7) {
		useblockf=true;
	}

    char file[256];
	strcpy(file,argv[1]);
	char treefile[256];
	strcpy(treefile,file);
	strcat(treefile,".tre");
	
	unsigned nchar=atoi(argv[2]);
	unsigned ratrep=atoi(argv[3]);
	float percentupwt=atof(argv[4]);
	float nwtchar=nchar*(percentupwt/100);
	float multiplier=atof(argv[5]);
	if ( multiplier < 2 ) { multiplier = 2; }
	
	unsigned maxtrees = 1000;
	if ( ratrep > 500 ) {
		maxtrees = multiplier * ratrep;
	}
	
	unsigned i, j, k;
	bool upwtsite[MAXSITE];
	bool changed;

	char logfile[256];
	strcpy(logfile,file);
	strcat(logfile,".log");
	char outfile[256];
	strcpy(outfile,file);
	strcat(outfile,".ratchblock");
	ofstream outf(outfile);
	
	// Start the ratchblock file
	outf << "#NEXUS\n\n";
	outf << "BEGIN PAUP;\n\tlog file=" << logfile << " start replace;\n";
	// Add the commands in blockf if blockf name is passed
	if (useblockf) {
		ifstream bf(argv[6]);
		line[0]='-';
		while (line[0]!='\0') {
			line[0]='\0';
			bf.getline(line,2047,'\n');
			if (line[0]!='\0') { outf << line << '\n'; }
		}
	}
	// Write the remainder of the ratchblock file
	outf << "[! commandfile created by ratchblock program: ]\n";
	outf << "[!   % sites upweighted in each replicate = " << percentupwt << " ]\n";
	outf << "[!   number of ratchet replicates = " << ratrep << " ]\n";
	outf << "[!   total number of sites = " << nchar << " ]\n";
	outf << "[! ]";
	outf << "\tset autoclose=yes warntree=no warnreset=no dstatus=60 status=no root=outgroup monitor=yes notifybeep=no taxlabels=full [errorstop=no] warnreset=no warntree=no warntsave=no showtaxnum=no [visnotif=none] [checkevts=yes] [background=yes] MaxTrees=100 increase=auto autoinc=100;\n";
	outf << "\tset criterion=parsimony;\n\tcleartrees;\n\thsearch start=stepwise multrees=no;\n\ttime;\n";
	outf << "\tsavetrees file=" << treefile << " brlens=yes replace=yes;";
		
	MTRand mtrand1;
	
	for (i=0; i<ratrep; i++) {
		cout << "Writing ratchet iteration # " << i+1 << '\n';
		outf << "[!\n*Ratchet iteration " << i+1 << "\n]\n";
		outf << "\twts 2: ";
		for (j=0; j<nchar; j++) {
			upwtsite[j]=false;
		}
		for (j=0; j<nwtchar; j++) {
			changed=false;
			while (!changed) {
				k=mtrand1.randInt(nchar-1);
				if (!upwtsite[k]) {
					changed=true;
					upwtsite[k]=true;
				}
			}
		}
		for (j=0; j<nchar; j++) {
			if (upwtsite[j]) outf << " " << j+1;
		}
		outf << ";\n";
		outf << "\thsearch start=current;\n\twts 1: all;\n\thsearch start=current;\n\ttime;\n";
		outf << "\tsavetrees file=" << treefile << " brlens=yes append=yes;\n\n";
	}
	
	outf << "\t[!\nGenerate consensus tree:\n]\n";
	outf << "\tgettrees file=" << treefile << " allblocks=yes mode=3 storebrlens=yes duptrees=eliminate;\n";
	outf << "\tfilter best;\n\tcontree all /strict=yes root=outgroup outroot=monophyl showtree=yes treefile=";
	outf << file << ".cons.tre replace=yes;\n";
	outf << "\tpscores /TL=yes CI=yes RI=yes RC=Yes HI=Yes;\n";
	outf << "\t[!\nSwap further on shortest tree:\n]\n";
	outf << "\tgettrees file=" << treefile << " allblocks=yes mode=3 storebrlens=yes duptrees=eliminate;\n";
	outf << "\tfilter best;\n";
	outf << "\tset maxtrees=" << maxtrees << " increase=no;\n";
	outf << "\tpset collapse=minBrlen;\n";
	outf << "\thsearch start=current swap=TBR steepest=no multrees=yes;\n";
	outf << "\tfilter best;\n";
	outf << "\tsavetrees file=" << file << ".finalswap.tre brlens=yes format=altnexus replace=yes;\n";
	outf << "\tcontree all /strict=yes root=outgroup outroot=monophyl showtree=yes treefile=" << file << ".finaltreecons.tre replace=yes;\n";
	outf << "\tpscores /TL=yes CI=yes RI=yes RC=Yes HI=Yes scorefile=" << file << ".scores.txt replace;\n";
	outf << "\n\tlog stop;\n";
	outf << "\tquit;\n";
	outf << "\nEND;\n";

	return 0;
}
