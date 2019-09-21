//
//  main.cpp
//  reverseSA
//
//  Created by muhammed oguzhan kulekci on 31.03.2019.
//  Copyright © 2019 muhammed oguzhan kulekci. All rights reserved.
//


#include <sdsl/suffix_arrays.hpp>
#include <sdsl/suffix_trees.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/lcp.hpp>
#include <string>
#include <iostream>
#include <algorithm>
#include <iomanip>

#include <fstream>
#include <math.h>
#include <sstream>
#include <stdlib.h>
#include <ctime>
#include <iomanip>
#include <divsufsort.h>
#include <chrono>
#include <random>

//default sdsl values are SA_SAMPLE_FREQ=32 ISA_SAMPLE_FREQ=64
#define SA_SAMPLE_FREQ  32
#define ISA_SAMPLE_FREQ 32

#define BLCD

using namespace sdsl;
using namespace std;

//#define SA_SAMPLING_STRAT sa_order_sa_sampling<>
#define SA_SAMPLING_STRAT text_order_sa_sampling<>


#ifdef BLCD
typedef csa_wt<wt_blcd<rrr_vector<127>>,SA_SAMPLE_FREQ,ISA_SAMPLE_FREQ,SA_SAMPLING_STRAT, isa_sampling<>> fmindexType;
#else
typedef csa_wt<wt_hutu<rrr_vector<127>>,SA_SAMPLE_FREQ,ISA_SAMPLE_FREQ,SA_SAMPLING_STRAT, isa_sampling<>> fmindexType;
#endif

fmindexType fm_index;
fmindexType Rfm_index;

char *buffer, *reversebuffer;

std::ifstream::pos_type file_size(const char* filename)
{
    std::ifstream in(filename, std::ifstream::ate | std::ifstream::binary);
    return in.tellg();
}


unsigned long long SuffixArray(unsigned long long k, unsigned int* steps, fmindexType FM_index, std::chrono::duration<long,std::micro>* durationSUS, std::chrono::duration<long,std::micro>* durationSAaccess)
//returns SA[k] on the reverse text, with the number of times while loop executed (steps), which is the length of SUS on the queried position.
//Besides, the time spent on SUS extraction and regular SA access via the FM-index are returned separately in microseconds
{
    std::chrono::high_resolution_clock::time_point start,end;
    
    //------------------------------------------------------------------------------------------
    //Phase-I: SUS extraction phase
    start = std::chrono::high_resolution_clock::now();
    unsigned long long leftBoundary, rightBoundary;
    unsigned char symbol;
    leftBoundary  = 0;
    rightBoundary = FM_index.size()-1;
    unsigned int m=0;
    while (rightBoundary>leftBoundary){
        auto res1 = quantile_freq(FM_index.wavelet_tree, leftBoundary, rightBoundary, k); //! Returns the k-th smallest element and its frequency in wt[leftBoundary..rightBoundary].
        symbol = get<0>(res1);
        auto res = FM_index.wavelet_tree.lex_count(leftBoundary, rightBoundary+1, symbol); //! How many symbols are lexicographic smaller/greater than c in [i..j-1].
                                                                                              //Returns A triple containing: rank(i,c) #symbols smaller than c in [i..j-1] #symbols greater than c in [i..j-1]
        k -= get<1>(res);
        backward_search(FM_index, leftBoundary, rightBoundary, symbol, leftBoundary, rightBoundary);
        m++;
    }
    end = std::chrono::high_resolution_clock::now();
    *durationSUS = std::chrono::duration_cast<std::chrono::microseconds>( end - start );
    //----------------------------------------------------------------------------------
    
    //----------------------------------------------------------------------------------
    //Phase-II: SA access on the Forward FM-Index
    start = std::chrono::high_resolution_clock::now();
    unsigned long long retVal =(FM_index[leftBoundary]-1+m) % FM_index.size();
    end = std::chrono::high_resolution_clock::now();
    *durationSAaccess = std::chrono::duration_cast<std::chrono::microseconds>( end - start );
    //----------------------------------------------------------------------------------
    retVal = FM_index.size()-retVal-2; // following two lines, two subtractions and one assignment operations are atomic, that run in less than milliseconds and creates a negligible overhead in measuring SAaccesstime in  Phase-II
    *steps = m;

    return retVal;
}


unsigned long long InverseSuffixArray(unsigned long long position, unsigned int* steps, fmindexType FM_index, std::chrono::duration<long,std::micro>* durationSUS, std::chrono::duration<long,std::micro>* durationISAaccess)
//returns ISA[position] on the reverse text, with the number of times while loop executed (steps), which is the length of SUS on the queried position.
//Besides, the time spent on SUS extraction and regular ISA access via the FM-index are returned separately in microseconds
{
    std::chrono::high_resolution_clock::time_point start,end;
    unsigned long long leftBoundary=0, rightBoundary = FM_index.size()-1, retVal=0, m=0;
    
    //------------------------------------------------------------------------------------------
    //Phase-I: Find the position on the reversed text, which requires to access to ISA on the FM-index
    start = std::chrono::high_resolution_clock::now();
    position = FM_index.size() - position - 2;
    unsigned long long k = FM_index.isa[position+1]; // Usually when we have the BWT, we throw away the text, thus we need to find k where T[position] is mapped to BWT[k]
    end = std::chrono::high_resolution_clock::now();
    *durationISAaccess = std::chrono::duration_cast<std::chrono::microseconds>( end - start );
    //------------------------------------------------------------------------------------------

    
    //------------------------------------------------------------------------------------------
    //Phase-II: Keep backwards search on the FM-index till reaching the SUS
    start = std::chrono::high_resolution_clock::now();
    unsigned char symbol = FM_index.bwt[k];
    while (rightBoundary>leftBoundary){
        auto res = FM_index.wavelet_tree.lex_count(leftBoundary, rightBoundary+1, symbol);
        retVal += get<1>(res);
        backward_search(FM_index, leftBoundary, rightBoundary, symbol, leftBoundary, rightBoundary);
        k = FM_index.lf[k];
        symbol = FM_index.bwt[k];
        m++;
    }
    end = std::chrono::high_resolution_clock::now();
    *durationSUS = std::chrono::duration_cast<std::chrono::microseconds>( end - start );

    //------------------------------------------------------------------------------------------
    *steps = m;
    return retVal;
}

/*
 // this is the version where we assume the forward text is available besides the BWT. So no need for LF mapping and ISA to detect the next symbol !!!
unsigned  long long InverseSuffixArrayWITHForwardText(unsigned long long position)
{
    position = fm_index.size() - position - 2;
    unsigned long long leftBoundary, rightBoundary;
    leftBoundary  = 0;
    rightBoundary = fm_index.size()-1;
    unsigned char* symbol;
    symbol = (unsigned char*)&buffer[position];
    unsigned long long retVal = 0;
    while (rightBoundary>leftBoundary){
        auto res = fm_index.wavelet_tree.lex_count(leftBoundary, rightBoundary+1, *symbol);
        retVal += get<1>(res);
        backward_search(fm_index, leftBoundary, rightBoundary, *symbol, leftBoundary, rightBoundary);
        symbol--;
    }
    return retVal;
}
*/


int main(int argc, char** argv)// arguments: inputfilename samplingRatio
{

    std::chrono::high_resolution_clock::time_point t1,t2;
    std::chrono::duration<long,std::micro> SUSduration,SAaccessduration,duration1,duration2;
    
    unsigned long int filesize = file_size(argv[1]);

    ifstream in(argv[1]);
    buffer = new char[filesize+1];
    memset(buffer, 0, filesize+1);
    in.read(&buffer[0],filesize);
    in.close();
    //cout << "\nFile Name and File Size in bytes\t:" << argv[1] << " " << filesize << endl;
   
    
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//construct the fm_index on forward text, this computes everything SA,ISA,BWT,LCP,....
    t1 = std::chrono::high_resolution_clock::now();
    construct_im(fm_index, (const char*)buffer,1);
    t2 = std::chrono::high_resolution_clock::now();
    //cout << "Forward FM-index creation elapsed time in milliseconds:\t"<<   std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count() << endl;
    /*
    cout << " i SA ISA PSI LF BWT   T[SA[i]..SA[i]-1]" << endl;
    csXprintf(cout, "%2I %2S %3s %3P %2p %3B   %:3T", fm_index);
    */
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//reverse the file
    reversebuffer = new char[filesize+1];
    memset(reversebuffer, 0, filesize+1);
    for(unsigned  int i=0;i<filesize;i++)    reversebuffer[i]= buffer[filesize-i-1];
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// compute the reverse FM-Index
    cst_sct3<fmindexType, lcp_vlc<>> cst;
    t1 = std::chrono::high_resolution_clock::now();
    construct_im(cst, (const char*)reversebuffer,1);
    Rfm_index = cst.csa;
    t2 = std::chrono::high_resolution_clock::now();
    //cout << "Reverse FM-index creation via CST elapsed time in milliseconds:\t"<<  std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count() << endl ;
    /*
    cout << " i SA ISA PSI LF BWT   T[SA[i]..SA[i]-1]" << endl;
    csXprintf(cout, "%2I %2S %3s %3P %2p %3B   %:3T", Rfm_index);
    */
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// generate random queries
    unsigned int sampleSize = stoi(argv[2]); //( samplingRatio * (unsigned int)filesize / 100 ; // 10% of the file size
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist(2,(unsigned int)filesize-100); // distribution in range [1, filesize-2]
    unsigned int* querySequence = new unsigned int[sampleSize];
    for(unsigned int s=0; s< sampleSize;s++)  querySequence[s]= dist(rng);
    
    unsigned int ss=0;
    for(unsigned int g=1;g<filesize-1;g++)
        if (((unsigned int)cst.lcp[g] > (unsigned int)SA_SAMPLE_FREQ) && ((unsigned int)cst.lcp[g+1] > (unsigned int)SA_SAMPLE_FREQ)) ss++;
    cout << 100.00 - 100.00 * (double)ss/(double)filesize  << "\t"; // number of positions in the file that has an SUS shorter than the SA_SAMPLE_FREQ = ISA_SMAPLE_FREQ
    
    ss=0;
    while(ss<sampleSize){ // randomly select among those eligible ones
        querySequence[ss] = dist(rng);
        if (((unsigned int)cst.lcp[querySequence[ss]] <= (unsigned int)SA_SAMPLE_FREQ) && ((unsigned int)cst.lcp[querySequence[ss]+1] <= (unsigned int)SA_SAMPLE_FREQ)) ss++;
    }
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------/


                                        /*************************  SUFFIX ARRAY TESTS ******************************************/

//------------------------------------------------------------------------------------------------------------------------------------------------------------------------/
//GRAND AVERAGE SA ACCESS TIME TEST: Compute the average access time for reverse SA[i] on the !!!REGULAR!!!  FM-index created over the reversed text
    unsigned long long *tempres = new unsigned long long[sampleSize];
    t1 = std::chrono::high_resolution_clock::now();
    duration1 = std::chrono::duration_cast<std::chrono::microseconds>(t1-t1); //set duration to zero
    for(unsigned int s=0; s< sampleSize;s++){
        t1 = std::chrono::high_resolution_clock::now();
        tempres[s] = Rfm_index[querySequence[s]];
        t2 = std::chrono::high_resolution_clock::now();
        duration1 += std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 );
    }
    double avg = (double)duration1.count() / (double)sampleSize;
    cout << avg << "\t";
   //cout << "Average Access Time to reverse SA[i] for a random position i on FM-Index with SA sampling frequency " << fm_index.sa_sample_dens << " is:\t " << (double)duration1.count() / (double)sampleSize << " microseconds" << endl;
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    
//GRAND AVERAGE SA ACCESS TIME TEST: Compute the average access time for reverse SA[i] by the !!!PROPOSED METHOD!!! using the forward FM-index  ------------------------------------------------------------------
    unsigned int lessThan[10]={0,0,0,0,0,0,0,0,0,0};
    unsigned int SUSlength,totalSUSlength=0;
    t1 = std::chrono::high_resolution_clock::now();
    duration1 = std::chrono::duration_cast<std::chrono::microseconds>(t1-t1); //set duration1 to zero
    duration2 = std::chrono::duration_cast<std::chrono::microseconds>(t1-t1); //set duration2 to zero
    for(unsigned int s=0; s< sampleSize;s++){
        tempres[s] = SuffixArray(querySequence[s], &SUSlength, fm_index, &SUSduration, &SAaccessduration);
        totalSUSlength+=SUSlength;
        duration1 += SUSduration;
        duration2 += SAaccessduration;
        switch ((SUSduration.count()+SAaccessduration.count()) / (unsigned int)avg) {
            case 0: lessThan[0]++; break; //proposed method achieves better performence in accessing the reverse SA than regular FM-index over reverse T
            case 1: lessThan[1]++; break; //performance is more than 1 regular SA access but better than 2
            case 2: lessThan[2]++; break; //performance is more than 2 regular SA access but better than 3
            case 3: lessThan[3]++; break; //performance is more than 3 regular SA access but better than 4
            case 4: lessThan[4]++; break; //performance is more than 4 regular SA access but better than 5
            case 5: lessThan[5]++; break; //performance is more than 5 regular SA access but better than 6
            case 6: lessThan[6]++; break; //performance is more than 6 regular SA access but better than 7
            case 7: lessThan[7]++; break; //performance is more than 7 regular SA access but better than 8
            case 8: lessThan[8]++; break; //performance is more than 8 regular SA access but better than 9
            default: lessThan[9]++; break; //performance is more than 9 regular SA access
        }
    }
    //cout << "Average Access Time to reverse SA[i] for a random position i by using forward FM-Index with SA sampling frequency " << fm_index.sa_sample_dens << " is:\t " << (double)(duration1.count())/(double)sampleSize  << " microseconds" << "\t" << (double)(duration2.count())/(double)sampleSize << "\t" << (double)(duration1.count() + duration2.count())/(double)sampleSize << "\t" << (double)totalSUSlength/(double)sampleSize << endl;
    cout << (double)(duration2.count())/(double)sampleSize << "\t" << (double)(duration1.count())/(double)sampleSize << "\t" << (double)(duration1.count()+duration2.count())/(double)sampleSize << "\t" << (double)totalSUSlength/(double)sampleSize << "\t"<< 100.00 *((double)(lessThan[0]) / (double)(sampleSize)) << "\t" << 100.00 *((double)(lessThan[1]) / (double)(sampleSize)) << "\t" << 100.00 *((double)(lessThan[2]) / (double)(sampleSize)) << "\t" << 100.00 *((double)(lessThan[3]) / (double)(sampleSize)) << "\t" <<  100.00 *((double)(lessThan[4]) / (double)(sampleSize)) << "\t" << 100.00 *((double)(lessThan[5]) / (double)(sampleSize)) << "\t" << 100.00 *((double)(lessThan[6]) / (double)(sampleSize)) << "\t" << 100.00 *((double)(lessThan[7]) / (double)(sampleSize)) << "\t" << 100.00 *((double)(lessThan[8]) / (double)(sampleSize)) << "\t" << 100.00 *((double)(lessThan[9]) / (double)(sampleSize)) << "\t\t";
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------


                                                        /************************* INVERSE SUFFIX ARRAY TESTS *******************************************/
    
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------/
//GRAND AVERAGE ISA ACCESS TIME TEST: Compute the average access time for reverse ISA[i] on the !!!REGULAR!!!  FM-index created over the reversed text
    t1 = std::chrono::high_resolution_clock::now();
    duration1 = std::chrono::duration_cast<std::chrono::microseconds>(t1-t1); //set duration to zero
    for(unsigned int s=0; s< sampleSize;s++){
        t1 = std::chrono::high_resolution_clock::now();
        tempres[s] = Rfm_index.isa[querySequence[s]];
        t2 = std::chrono::high_resolution_clock::now();
        duration1 += std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 );
    }
    avg = (double)duration1.count() / (double)sampleSize;
    cout << avg << "\t";
    //cout << "Average Access Time to reverse ISA[i] for a random position i on FM-Index with SA sampling frequency " << fm_index.isa_sample_dens << " is:\t " << (double)duration1.count() / (double)sampleSize << " microseconds" << endl;
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    

//------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//GRAND AVERAGE ISA ACCESS TIME TEST: Compute the average access time for reverse  ISA[i] by the !!!PROPOSED METHOD!!! using the forward FM-index
    
    //first take the inverse of the querySequence
    for(unsigned int i=0;i<sampleSize;i++) querySequence[i] = Rfm_index[querySequence[i]];

    lessThan[0]=lessThan[1]=lessThan[2]=lessThan[3]=lessThan[4]=lessThan[5]=lessThan[6]=lessThan[7]=lessThan[8]=lessThan[9]=0;
    totalSUSlength=0;
    t1 = std::chrono::high_resolution_clock::now();
    duration1 = std::chrono::duration_cast<std::chrono::microseconds>(t1-t1); //set duration1 to zero
    duration2 = std::chrono::duration_cast<std::chrono::microseconds>(t1-t1); //set duration2 to zero
    for(unsigned int s=0; s< sampleSize;s++){
        tempres[s] = InverseSuffixArray(querySequence[s], &SUSlength, fm_index, &SUSduration, &SAaccessduration);
        totalSUSlength+=SUSlength;
        duration1 += SUSduration;
        duration2 += SAaccessduration;
        switch ((SUSduration.count()+SAaccessduration.count()) / (unsigned int)avg) {
            case 0: lessThan[0]++; break; //proposed method achieves better performence in accessing the reverse SA than regular FM-index over reverse T
            case 1: lessThan[1]++; break; //performance is more than 1 regular SA access but better than 2
            case 2: lessThan[2]++; break; //performance is more than 2 regular SA access but better than 3
            case 3: lessThan[3]++; break; //performance is more than 3 regular SA access but better than 4
            case 4: lessThan[4]++; break; //performance is more than 4 regular SA access but better than 5
            case 5: lessThan[5]++; break; //performance is more than 5 regular SA access but better than 6
            case 6: lessThan[6]++; break; //performance is more than 6 regular SA access but better than 7
            case 7: lessThan[7]++; break; //performance is more than 7 regular SA access but better than 8
            case 8: lessThan[8]++; break; //performance is more than 8 regular SA access but better than 9
            default: lessThan[9]++; break; //performance is more than 9 regular SA access
        }
    }
    //cout << "Average Access Time to reverse SA[i] for a random position i by using forward FM-Index with SA sampling frequency " << fm_index.sa_sample_dens << " is:\t " << (double)(duration1.count())/(double)sampleSize  << " microseconds" << "\t" << (double)(duration2.count())/(double)sampleSize << "\t" << (double)(duration1.count() + duration2.count())/(double)sampleSize << "\t" << (double)totalSUSlength/(double)sampleSize << endl;
    cout << (double)(duration2.count())/(double)sampleSize << "\t" << (double)(duration1.count())/(double)sampleSize << "\t" << (double)(duration1.count()+duration2.count())/(double)sampleSize << "\t" << (double)totalSUSlength/(double)sampleSize << "\t"<< 100.00 *((double)(lessThan[0]) / (double)(sampleSize)) << "\t" << 100.00 *((double)(lessThan[1]) / (double)(sampleSize)) << "\t" << 100.00 *((double)(lessThan[2]) / (double)(sampleSize)) << "\t" << 100.00 *((double)(lessThan[3]) / (double)(sampleSize)) << "\t" <<  100.00 *((double)(lessThan[4]) / (double)(sampleSize)) << "\t" << 100.00 *((double)(lessThan[5]) / (double)(sampleSize)) << "\t" << 100.00 *((double)(lessThan[6]) / (double)(sampleSize)) << "\t" << 100.00 *((double)(lessThan[7]) / (double)(sampleSize)) << "\t" << 100.00 *((double)(lessThan[8]) / (double)(sampleSize)) << "\t" << 100.00 *((double)(lessThan[9]) / (double)(sampleSize)) << "\t" << endl;
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    
    
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// CORRECTNESS TEST: Test whether the proposed and the FM-Index compute the same SA/ISA values
    for(unsigned int s=0; s< sampleSize;s++){
        if (InverseSuffixArray(querySequence[s], &SUSlength, fm_index, &SUSduration, &SAaccessduration) !=Rfm_index.isa[querySequence[s]] )
            cout << "Error " << s << "\t" << InverseSuffixArray(querySequence[s], &SUSlength, fm_index, &SUSduration, &SAaccessduration) << "\t" << Rfm_index.isa[querySequence[s]] << endl;
        if (SuffixArray(querySequence[s], &SUSlength, fm_index, &SUSduration, &SAaccessduration) != Rfm_index[querySequence[s]])
            cout << "Error " << s << "\t" << SuffixArray(querySequence[s], &SUSlength, fm_index, &SUSduration, &SAaccessduration) << "\t" << Rfm_index[querySequence[s]] << endl;
        
    }
    cout << "Correctness test completed." << endl;
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    return 1;

}
