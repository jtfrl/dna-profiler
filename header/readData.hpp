#ifndef READ_DATA
#define READ_DATA

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map> 
#include <algorithm>
#include <set>


class DNAbase{
    
    std::string nomebase;
    std::vector<std::string> base_STR; //cada um dos STRs no header do csv
    std::map<std::string, int> freqSTR;
    inline static std::set<std::string> missingSTR;

    public:
    DNAbase(std::string u_nomebase, std::vector<std::string> base_STRu): 
    nomebase(u_nomebase), base_STR(base_STRu){
        for (const std::string& str : base_STRu){
            freqSTR[str]=0; //frequencia padrão de ocorrências
        }
    }

    std::string getbNome() const {
        return nomebase;
    }

    std::vector<std::string> getbaseSTR(){
        return base_STR;
    }

    void setSTRFreq(const std::string& str_name, int freq){
        if(freqSTR.count(str_name)){ //checa o str
            freqSTR[str_name]=freq; //associa cada STR com uma frequencia
            //importante ver que aqui se trata de uma associação individual
        }else{
            std::cerr<<"Alert: STR '"<<str_name<<" was not found.\n"<<std::endl;

        }
    }

    //frequencia de um STR
    int getFreqSTR(const std::string& str_name) const{
        if(freqSTR.count(str_name)){
            return freqSTR.at(str_name);
        }
        if(missingSTR.insert(str_name).second){
            std::cerr<<"The following STR was not found: "<<str_name<<std::endl;
        //o STR não foi encontrado
        }
        return -1;
    }

    //mostrar todas as frequencias
    const std::map<std::string, int>& getAllFreq() const{
        return freqSTR; 
    }


    void feedFreqFromCSVRows(const std::vector<std::string> & freq_values){
        //long long int tam=freqSTR.size();
        if(freqSTR.size()!=base_STR.size()){
            std::cerr<<"Error: Different sizes for freq. and expeceted STRs"<<std::endl;
            //"Frequências e tamanho de container para STRs diferentes."
            return;
        }

        for(size_t i=0; i<freqSTR.size(); i++){
            try{
                int f=std::stoi(freq_values[i]);
                freqSTR[base_STR[i]]=f;
            }catch(const std::invalid_argument& e){
                std::cerr<<"Error: invalid freq. value. for a given STR\n | \n "<<"freq.:"<<freq_values[i]<<" | "<<base_STR[i]<<e.what();
                return;
            }catch(const std::out_of_range& e){
                std::cerr<<"Error: freq. value "<<freq_values[i]<<" out of range for STR"<<e.what();

            }
        }
    }

    //leitura de linhas da base 
    static std::vector<DNAbase> readingRowsFromDNAbase(
        std::ifstream& fin, const std::vector<std::string>& STRsfromBase){

        std::vector<DNAbase> dna_db;
        std::string row;

        //laço para o resto do arquivo
        while(std::getline(fin, row)){ 
            std::stringstream ss(row);
            std::string name;
            std::vector<std::string> values; //cada valor das linhas

            if(!std::getline(ss, name, ',')){
                std::cerr<<"Could not read names from CSV file. Problem on: "<<row<<std::endl;
                continue;
            }

            std::string cell;
            while(std::getline(ss, cell, ',')){ //cada célula é lida
                values.push_back(cell);
            }


            if(values.size() != STRsfromBase.size()){
                std::cerr<<"Different sizes: could not count properly STRs. Problem on:  "<<row<<std::endl;
                continue;
            }

            DNAbase person(name, STRsfromBase);
            person.feedFreqFromCSVRows(values); //frequencias "entram" no objeto
            dna_db.push_back(person);

        }

        return dna_db;
    }

    
};


class DNAt{

    std::string DNAsq;
    std::map<std::string, int> sq_freqSTR; //frequencia de STR para cada padrão repetido
    std::vector<std::string> tSTR; //os STRs presentes no arquivo
    size_t seqtam;

    public:

    DNAt(std::string u_DNAsq, std::map<std::string, int> sq_freqSTRu={}): DNAsq(u_DNAsq), 
    sq_freqSTR(sq_freqSTRu), seqtam(u_DNAsq.length()){
    } 


    std::string getDNAsq(){
        return DNAsq;
    }

    std::map<std::string, int> getsq_freqSTR(){
        return sq_freqSTR;
    }

    size_t getSeqTam() const{
        return seqtam;
    }


    // verifica quantos STR estão presentes, assim como a sua freq.
    int getFreqTestDNA(const std::map<std::string, int>& sq_freqSTR, std::string STRtoBeVer) const{
        if(sq_freqSTR.count(STRtoBeVer)){
            return sq_freqSTR.at(STRtoBeVer);
        }

        std::cerr<<"The following STR was not found at the tested DNA sequence: "<<STRtoBeVer<<std::endl;
        return -1;
    }
    //std::map<std::vector<std::string>, std::pair<int, int> > possibleSTR() const;
   // int contaVizinhoRpt(const std::string& str,const std::map<std::vector<std::string>, 
    //    std::pair<int, int> >& possible_STR) const;


    //leitura de string (DNA a ser testado); fin aqui -> fin2 (dnaprofiler.cpp)
    static DNAt readingLikeOneUniqueRow(std::ifstream& fin){
        std::string DNA2betested;
        std::string uniqueLine;

        //le todo o arquivo como uma única linha
        while(std::getline(fin, uniqueLine)){
            if(uniqueLine.empty()|| uniqueLine[0]=='>') continue;
            DNA2betested+=uniqueLine;
        }
        
        for(char c:DNA2betested){
            c=toupper(c);
            if(c!='A' && c!='C' && c!='T' && c!='G') throw std::runtime_error("Invalid: " + 
                std::string(1,c));
        }
          
        std::map<std::string, int> freqMap;
        return DNAt(DNA2betested, freqMap);
    }


    void buildDNAProfile(const std::vector<std::string>& findSTRs);

    //retorna o perfil
    const std::map<std::string, int>& getProfile() const { return sq_freqSTR;}
    
    
};



class DNAMatcher {

    public:

    //comparação com a base e o perfil de DNA gerado
    static std::vector<std::string> findMatches( const DNAt& test_dna,
        const std::vector<DNAbase>& dna_db){

        std::vector<std::string> matches;
        const std::map<std::string, int>& profileInTest = test_dna.getProfile();

        for(const auto& person:dna_db){
            if(isMatch(profileInTest, person)){
                matches.push_back(person.getbNome());
            }
        }

        return matches;
    }

    private:

    //test_profile é justamente a variável que vai estar como DNAt logo acima
    static bool isMatch(const std::map<std::string, int>& test_profile, 
        const DNAbase& person){
    
        const std::map<std::string, int>& person_profile = person.getAllFreq();

        if(test_profile.size()!=person_profile.size()){
            return false;
        }

        for(const auto& [str, person_count]:person_profile){
            auto it=test_profile.find(str);
            if(it==test_profile.end()){
                return false; //STR não encontrado dentro do perfil de DNA
            }
    

            if(test_profile.at(str)!=person_count){
                return false;
            }
        }

        return true;


    }

};


#endif
