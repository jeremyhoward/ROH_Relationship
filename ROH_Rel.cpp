#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <string>
#include <vector>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <omp.h>

/**********************************************************/
/*      ROH Index Object (Keeps track of ROH Windows)     */
/**********************************************************/
class ROH_Index
{
private:
    int Chromosome;                 /* Which Chromosome belongs to */
    int StartPosition;              /* Start position (Mb) of roh window */
    int EndPosition;                /* End position (Mb) of roh window */
    int StartIndex;                 /* Start position (column number) of roh window */
    int EndIndex;                   /* End position (column number) of roh window */
    int NumberSNP;                  /* Number of SNP in ROH Window */
public:
    // Constructors
    ROH_Index();
    ROH_Index(int chr = 0, int stpos = 0, int enpos = 0, int stind = 0, int enind = 0, int numsnp = 0);
    // Destructors
    ~ROH_Index();
    // Functions to grab parts of object
    int getChr(){return Chromosome;}
    int getStPos(){return StartPosition;}
    int getEnPos(){return EndPosition;}
    int getStInd(){return StartIndex;}
    int getEnInd(){return EndIndex;}
    int getNumSNP(){return NumberSNP;}
};
// Constructors ROH_Index
ROH_Index::ROH_Index(){Chromosome = 0;  StartPosition = 0; EndPosition = 0; StartIndex = 0; EndIndex = 0; NumberSNP = 0;}
ROH_Index::ROH_Index(int chr, int stpos, int enpos, int stind, int enind, int numsnp)
{
    Chromosome = chr;  StartPosition = stpos; EndPosition = enpos; StartIndex = stind; EndIndex = enind; NumberSNP = numsnp;
}
// Destructors
ROH_Index::~ROH_Index(){}

using namespace std;

int main(int argc, char* argv[])
{
    cout<<"\n###############################################################\n";
    cout<<"####################################################   ########\n";
    cout<<"################################################   /~|   ######\n";
    cout<<"#############################################   _- `~~~', #####\n";
    cout<<"###########################################  _-~       )  #####\n";
    cout<<"########################################  _-~          |  #####\n";
    cout<<"#####################################  _-~            ;  ######\n";
    cout<<"###########################  __---___-~              |   ######\n";
    cout<<"########################   _~   ,,                  ;  `,,  ###\n";
    cout<<"######################  _-~    ;'                  |  ,'  ; ###\n";
    cout<<"####################  _~      '                    `~'   ; ####\n";
    cout<<"#############   __---;                                 ,' #####\n";
    cout<<"#########   __~~  ___                                ,' #######\n";
    cout<<"######  _-~~   -~~ _          N                    ,' #########\n";
    cout<<"###### `-_         _           C                  ; ###########\n";
    cout<<"########  ~~----~~~   ;         S                ; ############\n";
    cout<<"##########  /          ;         U               ; ############\n";
    cout<<"########  /             ;                      ; ##############\n";
    cout<<"######  /                `                    ; ###############\n";
    cout<<"####  /                                      ; ################\n";
    cout<<"##                                            #################\n";
    cout <<"+-----+-----+-----+-----+-----+-----+-----+------+-----+-----+ " <<endl;
    cout <<"| This program is free software: you can redistribute it and |" << endl;
    cout <<"| /or modify it under the terms of the GNU General Public    |" << endl;
    cout <<"| License as published by the Free Software Foundation,      |" << endl;
    cout <<"| either version 3 of the License, or (at your option) any   |" << endl;
    cout <<"| later version.                                             |" << endl;
    cout <<"+-----+-----+-----+-----+-----+-----+-----+------+-----+-----+ " <<endl;
    cout <<"|  Optimized Run-of-Homozygosity                             |" << endl;
    cout <<"|  Genomic Relationship Creator                              |" << endl;
    cout <<"|  Author: Jeremy T. Howard                                  |" << endl;
    cout <<"|  Updated: Date: May 5, 2018                                |" << endl;
    cout <<"+-----+-----+-----+-----+-----+-----+-----+------+-----+-----+ " <<endl;
    if(argc != 7){cout << "Wrong nummber of parameters given" << endl; exit (EXIT_FAILURE);}
    time_t begin_time = time(0);
    cout << "- Parameter File Information:" << endl;
    /* Genotype File: ID then string of genotypes */
    string genofile = argv[1];
    /* Map File: Chromosome Mb position with a space dilemeter */
    string mapfile = argv[2];
    /* Mb cutoff for ROH */
    int roh_cutoff = atoi(argv[3]);
    roh_cutoff = roh_cutoff * 1000000;
    /* Number times the mean to cutoff if too low of SNP in ROH */
    int roh_threshold = atoi(argv[4]);
    /* Number of threads */
    int cores = atoi(argv[5]);
    omp_set_num_threads(cores);
    /* output file name */
    string genooutfile = argv[6];
    cout << "   - Genotype File: " << genofile << endl;
    cout << "   - Map File: " << mapfile << endl;
    cout << "   - ROH Nucleotide Cutoff: " << roh_cutoff << endl;
    cout << "   - ROH Threshold: " << roh_threshold << endl;
    cout << "   - Number of Threads: " << cores << endl;
    cout << "   - Output File Name: " << genooutfile << endl << endl;
    cout << "- Reading in Map file:" << endl;
    /* Read in map file don't need to know how many SNP are in the file */
    vector <string> numbers;
    string line;
    /* Import file and put each row into a vector then grab chromsome and position */
    ifstream infile;
    infile.open(mapfile.c_str());
    while (getline(infile,line)){numbers.push_back(line);}
    /* Total number of SNP */
    int rows = numbers.size();
    cout << "   - Total Number of SNP in Map file is " << rows << endl;
    /* stores chromosome in vector */
    vector < int > chr(rows,0);
    /* stores position in Mb */
    vector < int > positionMb(rows,0);
    /* position SNP is referring to in genotype string it is on */
    int index[rows];
    for(int i = 0; i < numbers.size(); i++)
    {
        string temp = numbers[i];                       // grab a line
        /* Check to see if tab delimeter exists */
        size_t pos = temp.find("\t", 0);                 // Find position where delimiter is at
        if(pos != std::string::npos){cout << endl << "Delimeter in Map file is a tab!! Needs to be a space." << endl; exit (EXIT_FAILURE);}
        else{
            vector < string > solvervariables(31,"");
            for(int j = 0; j < 31; j++)
            {
                size_t pos = temp.find(" ",0);
                solvervariables[j] = temp.substr(0,pos);
                solvervariables[j].erase(remove(solvervariables[j].begin(), solvervariables[j].end(), ' '),solvervariables[j].end());
                if(pos != std::string::npos){temp.erase(0, pos + 1);}
                if(pos == std::string::npos){temp.clear(); j = 31;}
            }
            int start = 0;
            while(start < solvervariables.size())
            {
                if(solvervariables[start] == ""){solvervariables.erase(solvervariables.begin()+start);}
                if(solvervariables[start] != ""){start++;}
            }
            if(solvervariables.size() != 2)
            {
                cout<<endl<<"Map file did not get parsed out to chr and pos correctly (Line"<< i+1 << ")!!"<<endl;
                exit (EXIT_FAILURE);
            }
            chr[i] = atoi(solvervariables[0].c_str());          // Convert it to an integer
            positionMb[i] = atoi(solvervariables[1].c_str());                   // Convert it to an integer
            index[i] = i;
        }
    }
    cout << "       - First line got parsed out to: " << endl;
    cout << "           - Chromosome: '" << chr[0] << "'" << endl;
    cout << "           - Position: '" << positionMb[0] << "'" << endl;
    numbers.clear();                                    // clear vector that holds each row
    /* Double check to make sure it is ordered chr and pos */
    for(int i = 1; i < chr.size(); i++)
    {
        if(chr[i] == chr[i-1])
        {
            if(positionMb[i] < positionMb[i-1])
            {
                cout<<endl<<"Map file is not sorted chromosome then position (Look at line "<<i<<" "<<i-1<<")!!"<<endl;
                exit (EXIT_FAILURE);
            }
        }
    }
    numbers.clear();
    cout << "   - Indexing Positions for ROH:" << endl;
    /* Fill ROH Index Object with each ROH window */
    vector < ROH_Index > roh_index;
    /* Create index to grab correct columns from genotype file when constructing ROH and Autozygosity */
    for(int i = 0; i < rows; i++)
    {
        int diff[rows];
        for(int j = 0; j < rows; j++)
        {
            /* if within the same chromosome then can be in an roh */
            if(chr[i] == chr[j]){diff[j] = positionMb[j] - positionMb[i];}
            /* if in another chromosome give it a 0 */
            if(chr[i] != chr[j]){diff[j] = 0;}
        }
        /* determine where to grab */
        int j = 0;
        while(j < rows)
        {
            /* Determine where first one that is over threshold is reached */
            if(diff[j] > roh_cutoff){break;}
            j++;
        }
        /* Tabulate number of SNP in roh window */
        int numsnp = index[j] - index[i] + 1;
        if(j < rows)
        {
            /* fill ROH index with the window */
            ROH_Index roh_temp(chr[i],positionMb[i],positionMb[j],index[i],index[j],numsnp);
            roh_index.push_back(roh_temp);                                      /* store in vector of roh_index objects */
        }
    }
    cout << "       - Total number of ROH windows prior to removal " << roh_index.size() << endl;
    /* Figure out mean and Standard Deviation of number of SNP within ROH windows */
    double sum = 0.0; double sq_sum = 0.0;
    for(int i = 0; i < roh_index.size(); i++)
    {
        sq_sum += roh_index[i].getNumSNP() * roh_index[i].getNumSNP();
        sum += roh_index[i].getNumSNP();
    }
    double mean = sum / double(roh_index.size());
    double stdev = sqrt(sq_sum / roh_index.size() - mean * mean);
    cout << "           - Mean +/- S.D. SNP size for an ROH: " << mean << " " << stdev << endl;
    int SNPSizeCutoff = mean - (roh_threshold * stdev + 0.5);     /* Ensures it rounds up */
    cout << "           - Any ROH window with SNP size below: " << SNPSizeCutoff << " removed." << endl;
    /* Remove ROH windows that fall below threshold */
    int row = roh_index.size();                                                     /* Figure out number of rows */
    int i = 0;                                                                      /* counter to determine where you are at */
    /* Keep ROH that are greater than a given threshold */
    while(i < row)
    {
        while(1)
        {
            if(roh_index[i].getNumSNP() < SNPSizeCutoff)
            {
                roh_index.erase(roh_index.begin()+i);   /* Remove from allfreq */
                row = row - 1;                          /* One less column in G and row in allfreq */
                break;
            }
            else
            {
                i++;
                break;
            }
        }
    }
    cout << "       - Total number of ROH windows after removal " << roh_index.size() << endl << endl;
    cout << "- Reading in Genotype file:" << endl;
    /* Stores Individual ID */
    vector < string > ID;
    /* Stores genotypes as a string */
    vector < string > geno;
    /* Read in Genotype file */
    ifstream infile1;
    infile1.open(genofile.c_str());
    while (getline(infile1,line))
    {
        size_t pos = line.find("\t", 0);                 // Find position where delimiter is at
        if(pos != std::string::npos){cout << endl << "Delimeter in Genotype file is a tab!! Needs to be a space." << endl; exit (EXIT_FAILURE);}
        /* Check to make sure it has two column ID and genotype */
        vector < string > solvervariables(31,"");
        for(int j = 0; j < 31; j++)
        {
            size_t pos = line.find(" ",0);
            solvervariables[j] = line.substr(0,pos);
            solvervariables[j].erase(remove(solvervariables[j].begin(), solvervariables[j].end(), ' '),solvervariables[j].end());
            if(pos != std::string::npos){line.erase(0, pos + 1);}
            if(pos == std::string::npos){line.clear(); j = 31;}
        }
        int start = 0;
        while(start < solvervariables.size())
        {
            if(solvervariables[start] == ""){solvervariables.erase(solvervariables.begin()+start);}
            if(solvervariables[start] != ""){start++;}
        }
        if(solvervariables.size() != 2)
        {
            cout<<endl<<"Genotype file did not get parsed out to chr and pos correctly (Line"<< i+1 << ")!!"<<endl;
            exit (EXIT_FAILURE);
        }
        ID.push_back(solvervariables[0]);                                /* Save ID */
        geno.push_back(solvervariables[1]);                              /* Save Genotype File */
    }
    cout << "   - Number of animals with genotypes: " << geno.size() << "." << endl;
    cout << "   - Number of Genotypes per animal: " << geno[0].size() << "." << endl << endl;
    cout << "- Constructing ROH Relationship Matrix:" << endl;
    /* store relationship values in 2-D vector and build first and initialize to 0.0 */
    double* ROH_Relationship = new double[geno.size()*geno.size()];
    for(int i = 0; i < geno.size(); i++)
    {
        for(int j = 0; j < geno.size(); j++)
        {
            ROH_Relationship[(i*geno.size())+j] = 0.0;
        }
    }
    cout << "   - ROH Relationship Matrix (Size = " << geno.size() << " x " << geno.size() << ")." << endl;
    /* Loop across all ROH Windows  */
    for(int r = 0; r < roh_index.size(); r++)
    {
        /* First step is to get number of unique haplotypes to set dimension of Haplotype similarity matrix (H) */
        /* stores each unique ROH roh */
        vector < string > haplotypes(geno.size()*2,"");
        vector < string > animalpatstring(geno.size(),"");
        vector < string > animalmatstring(geno.size(),"");
        /* Index for paternal haplotype number of a given animal */
        vector < int > AnimalPatHap(geno.size(),-5);
        /* Index for maternal haplotype number of a given animal */
        vector < int > AnimalMatHap(geno.size(),-5);
        /* loop across animals */
        #pragma omp parralel for
        for(int i = 0; i < geno.size(); i++)
        {
            /* Paternal haplotypes */
            string homo1 = geno[i].substr(roh_index[r].getStInd(),roh_index[r].getNumSNP());
            string homo2 = homo1;                           /* Maternal haplotypes */
            for(int g = 0; g < homo1.size(); g++)
            {
                if(homo1[g] == '0'){homo1[g] = '1'; homo2[g] = '1';}       /* a1a1 genotype */
                if(homo1[g] == '2'){homo1[g] = '2'; homo2[g] = '2';}       /* a2a2 genotype */
                if(homo1[g] == '3'){homo1[g] = '1'; homo2[g] = '2';}       /* a1a2 genotype */
                if(homo1[g] == '4'){homo1[g] = '2'; homo2[g] = '1';}       /* a2a1 genotype */
            }
            animalpatstring[i] = homo1;
            animalmatstring[i] = homo2;
            haplotypes[(i*2)+0] = homo1;
            haplotypes[(i*2)+1] = homo2;
        }
        /* Now Sort them and only keep unique ones */
        sort(haplotypes.begin(),haplotypes.end());
        haplotypes.erase(unique(haplotypes.begin(),haplotypes.end()),haplotypes.end());
        #pragma omp parralel for
        for(int i = 0; i < geno.size(); i++)
        {
            int j = 0;
            string foundpat = "nope"; string foundmat = "nope";
            /* assign paternal and maternal string to a numeric value */
            while(foundpat == "nope" || foundmat == "nope")
            {
                if(animalpatstring[i] == haplotypes[j]){AnimalPatHap[i] = j; foundpat = "yes";}
                if(animalmatstring[i] == haplotypes[j]){AnimalMatHap[i] = j; foundmat = "yes";}
                j++;
            }
            if(AnimalPatHap[i] == -5){cout << "error 1" << endl; exit (EXIT_FAILURE);}
        }
        animalpatstring.clear(); animalmatstring.clear();
        /* all unique haplotypes are tabulated now create Haplotype similarity matrix (H) in an array format */
        double* H_matrix = new double[haplotypes.size()*haplotypes.size()];
        int i, j;
        #pragma omp parallel for private(j)
        for(i = 0; i < haplotypes.size(); i++)
        {
            for(j = 0; j < haplotypes.size(); j++)
            {
                if(i == j){H_matrix[(i*haplotypes.size())+j] = 1.0;}
                if(i != j){H_matrix[(i*haplotypes.size())+j] = 0.0;}
            }
        }
        /* Haplotype similarity matrix (H) created now construct add result to ROH_Relationship Matrix */
        /* This results in an animal being related only if it has the exact same haplotype, if not the */
        /* exact same assume unrelated */
        #pragma omp parallel for private(j)
        for(i = 0; i < geno.size(); i++)
        {
            for(j = i; j < geno.size(); j++)
            {
                ROH_Relationship[(i*geno.size())+j] += (H_matrix[(AnimalPatHap[i]*haplotypes.size())+AnimalPatHap[j]] +
                                                        H_matrix[(AnimalPatHap[i]*haplotypes.size())+AnimalMatHap[j]] +
                                                        H_matrix[(AnimalMatHap[i]*haplotypes.size())+AnimalPatHap[j]] +
                                                        H_matrix[(AnimalMatHap[i]*haplotypes.size())+AnimalMatHap[j]]) / double(2);
                ROH_Relationship[(j*geno.size())+i] = ROH_Relationship[(i*geno.size())+j];
            }
        }
        if(r % 1000 == 0 && r != 0){cout << "       - " << r << endl;}
    }
    /* Once finished computing the relationship matrix now divide by total number of ROH windows */
    int j;
    #pragma omp parallel for private(j)
    for(i = 0; i < geno.size(); i++)
    {
        for(j = i; j < geno.size(); j++)
        {
            ROH_Relationship[(i*geno.size())+j] = ROH_Relationship[(i*geno.size())+j] / roh_index.size();
            ROH_Relationship[(j*geno.size())+i] = ROH_Relationship[(i*geno.size())+j];
        } /* Finish loop across ind2 */
    } /* Finish loop across ind1 */
    cout << "   - Finished Constructing ROH Relationship Matrix (Size = " << geno.size() << " x ";
    cout << geno.size() << ")." << endl << endl;
    cout << "- Output Matrix :" << endl << endl;
    /* Output Matrix to file */
    fstream checkRel; checkRel.open(genooutfile.c_str(), std::fstream::out | std::fstream::trunc); checkRel.close();
    std::ofstream output9(genooutfile.c_str(), std::ios_base::app | std::ios_base::out);
    for(int ind1 = 0; ind1 < geno.size(); ind1++)
    {
        for(int ind2 = 0; ind2 < geno.size(); ind2++)
        {
             if(ind2 >= ind1)
             {
                output9 << ID[ind1] << " "<< ID[ind2] <<" ";
                output9 << ROH_Relationship[(ind1*geno.size())+ind2];
                output9<<endl;
             }
        } /* Finish loop across ind2 */
    } /* Finish loop across ind1 */
    time_t end_time = time(0);
    cout << "- Finished with Program (Took: " << difftime(end_time,begin_time) << " seconds)" << endl << endl;
}

