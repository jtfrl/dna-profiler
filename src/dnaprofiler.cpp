#include "readData.hpp"
#include "visuals.hpp"

int main(int argc, char* argv[]){
    if(argc!=5){
        std::cerr<<"Program usage: "<<argv[0]<<
        "-d <database_file> -s <dna_sequence_file>"<<std::endl;
        return 1;
    }
    //cada um dos parâmetros
    std::string dbFlag(argv[1]), dbFile(argv[2]), seqFlag(argv[3]), seqFile(argv[4]);
    
    if ((dbFlag != "-d" || seqFlag != "-s")) {
        std::cerr << "Error. Expected format: " << argv[0] << 
        " -d <database_file> -s <dna_sequence_file>" << std::endl;
        return 1; 
    }

    /* ::::::::: DNA DATABASE ::::::::: */

    std::ifstream fin(dbFile); //base como file input

    if(!fin.is_open()){
        std::cerr<<"Error: Could not open .csv file: "<<dbFile<<std::endl;
        return 1;
    }

    std::string header; //toda a linha de cabeçalho (dbFile)
    std::vector<std::string> STRsinBaseFile; //variável a ser usada na classe DNAbase

    if(std::getline(fin, header)){
        std::stringstream ss(header);
        std::string column;
        bool col1=true;

        while(std::getline(ss, column, ',')){
            if(col1){
                col1=false;
                continue;
            }
            STRsinBaseFile.push_back(column); 
        }
    }

    //uso da classe para tomar os dados do CSV
    std::vector<DNAbase> dna_database= DNAbase::readingRowsFromDNAbase(fin, STRsinBaseFile);


    /* ::::::::: DNA SEQUENCE & GETTING PROFILE ::::::::: */

    std::ifstream fin2(seqFile); //file input: DNA a ser testado

    if(!fin2){
        std::cerr<<"DNA file could not be opened. Try again. \n";
        return 1;
    }

    /* ############### FINAL CHECK ################ */

    try{

        UI::printProgramHeader("JEFFERSON");
        UI::showLoadingAnimation("DNA database file [" + dbFile + "]");
        UI::showLoadingAnimation("Allegedly unknown DNA sequence file [" + seqFile + "]");
        UI::showProgressBar("Searching the database for a match...");

        DNAt test_dna=DNAt::readingLikeOneUniqueRow(fin2);
        test_dna.buildDNAProfile(STRsinBaseFile);

        std::vector<std::string> matches = DNAMatcher::findMatches(test_dna, dna_database);


        if(!matches.empty()){
            UI::printMatchResult(matches[0]);

            const std::string& dna_seq=test_dna.getDNAsq();
            const auto& profile=test_dna.getProfile();
            //mapa das posições dos STRs
           // std::map<size_t, std::pair<std::string, int>> poSTR; 
            UI::printHighlightedDna(dna_seq, profile);
            UI::printStrProfileSummary(profile);
        }else{
             UI::printNoMatchFound(); 
        }
    }catch(const std::exception& e){
        UI::printError(e.what()); 
        return 1;
    }    
    return 0; 
    
}