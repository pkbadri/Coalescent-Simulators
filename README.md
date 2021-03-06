# Standard Coalescent Simulator


https://github.com/pkbadri/Coalescent-Simulators/blob/main/arg.c

https://github.com/pkbadri/Coalescent-Simulators/blob/main/mtrand.h

https://github.com/pkbadri/Coalescent-Simulators/blob/main/mtrand.cpp



Download the 3 files in the links above. arg.c is a C++ program that can efficiently simulate datasets under the standard coalescent with uniform mutation rate, constant population size and with either non-uniform or uniform crossing-over and conversion rates. Parameters can be changed at the top of the code. Can be compiled using: g++ -O3 -o arg arg.c mtrand.cpp -Wno-deprecated and run by typing ./arg. Outputs SNP count and SNP locations with the alleles at those locations.












# Fixed Segregating Sites Simulator


https://github.com/pkbadri/Coalescent-Simulators/blob/main/fixedseg.c

https://github.com/pkbadri/Coalescent-Simulators/blob/main/mtrand.h

https://github.com/pkbadri/Coalescent-Simulators/blob/main/mtrand.cpp



Download the 3 files in the links above. fixedseg.c allows the user to fix the positions of segregating sites in the sequence. Parameters can be changed at the top of the code. Can be compiled by typing: g++ -O3 -o fixedseg fixedseg.c mtrand.cpp -Wno-deprecated and run by typing ./fixedseg. Outputs SNP count, their positions and the haplotypes in 0-1 format where 0 denotes ancestral allele and 1 denotes derived allele.
