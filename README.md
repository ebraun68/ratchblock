# ratchblock
Code to implement a simple parsimony ratchet search in PAUP

This simple C++ program uses Mersenne twister random generator in MersenneTwister.h

It is straightforward to compile. Simply place ratchblock.cpp and MersenneTwister.h in the
same directory and use a standard compiler. For example:

```
g++ -c ratchblock.cpp
g++ -o ratchblock ratchblock.o
```

If you execute the ratchblock file without any arguments the program will output information
regarding its use.

```
ratchblock: simple PAUP parsimony ratchet program

 usage:
  ratchblock <outfile> <nchar> <# replicates> <% upweight> <multiplier> <blockf>

  outfile      = name used for output files
    (ouput files include the ratchet block, log, and treefiles)
  nchar        = number of characters in the data file
  # replicates = number of ratchet searches
  % upweight   = percentage (1-100) of sites upweighted
  multiplier   = maxtrees multiplier
    (final maxtrees value for swapping will be # replicates * multiplier)
  blockf       = name of file with additional commands
    NOTE: blockf is optional; if used, it should contain a set of
          commands executed at the beginning of the analysis
```

Typical use would be:

```
ratchblock outputfile 5000 100 25 10
```

This will write a file with a PAUP block that can be appended to a data matrix (in this case,
I have assumed the data matrix has 5000 characters). The block generated by this command will
conduct an analysis with 100 replicates. 25% of the characters will be upweighted per replicate.

Typically, I run a series of searches with 15%, 20% 25%, 30%, and 35% upweighting.

I have uploaded ratchblock.cpp and MersenneTwister.h to this repository. I did not write
MersenneTwister.h (it was written by Richard J. Wagner) and you will note that I have left the
copyright information intact, as requested by the author, who states:

```
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (given details provided in MersenneTwister.h file)
```
