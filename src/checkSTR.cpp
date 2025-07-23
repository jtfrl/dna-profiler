#include "readData.hpp"
#include <cmath>
#include <utility>




/*
//verifica as partições úteis para verificar os STRs
std::vector<int> seekForDivsinSTR(size_t seqtam){
    std::vector<int> divs_inSTR;

    if(seqtam<2){
        std::cout<<"Impossible to get STRs (too short)\n"; return divs_inSTR;
    } 

    const int MAX_PART=seqtam/2;
    const int MIN_PART=2;


    for(int len=MIN_PART; len<MAX_PART; ++len){
        int v=seqtam/len;
        divs_inSTR.push_back(v);
    }

    return divs_inSTR;
}
*/


int countMaxRepeats(const std::string& dnasq, const std::string& str){
    int max_repeats=0;
    size_t str_len=str.length();


    if(str_len==0){
        return 0;
    }


    for(size_t i=0; i<dnasq.length(); ++i){
        int repCurrent=0;
        while(i+(repCurrent+1)*str_len<=dnasq.length() && 
        dnasq.substr(i+repCurrent*str_len, str_len) == str){ 
            //compara todas as sequências contíguas
            repCurrent++;
        }
    
        if(repCurrent>0){
            if(repCurrent>max_repeats){
                max_repeats=repCurrent;
            }
            i+=(repCurrent*str_len)-1; //move-se para uma nova contagem
        }
    
    }

    return max_repeats;
}


void DNAt::buildDNAProfile(const std::vector<std::string>& findSTRs){
    sq_freqSTR.clear(); //dados já usados são apagados

    for(const auto& str: findSTRs){ //para cada STR no vetor
        sq_freqSTR[str]=countMaxRepeats(DNAsq, str);
        //atualiza a contagem de STRs.
    }

}

