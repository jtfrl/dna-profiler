#ifndef VISUALS_HPP
#define VISUALS_HPP

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <string>
#include <chrono>
#include <thread>
#include <cstdlib> 
#include <cstddef> 

size_t findLongestRunPos(const std::string& dnasq, const std::string str){
    int maxRp=0;
    size_t best_pos=std::string::npos; //avalia toda a sequência
    size_t str_len=str.length();

    if(str_len==0) return best_pos;


    for (size_t i = 0; i < dnasq.length(); ++i) {
        int current_repeats = 0;
        while (i + (current_repeats + 1) * str_len <= dnasq.length() &&
               dnasq.substr(i + current_repeats * str_len, str_len) == str) {
            current_repeats++;
        }
        
        if (current_repeats > maxRp) {
            maxRp = current_repeats;
            best_pos = i;
        }

        if (current_repeats > 0) {
            i += (current_repeats * str_len) - 1;
        }
    }
    return best_pos;

}

namespace UI{

    //cores em uso (ANSI)
    const std::string RESET   = "\033[0m";
    const std::string BOLD    = "\033[1m";
    const std::string RED     = "\033[1;31m";
    const std::string GREEN   = "\033[1;32m";
    const std::string YELLOW  = "\033[1;33m";
    const std::string BLUE    = "\033[1;34m";
    const std::string MAGENTA = "\033[1;35m";
    const std::string CYAN    = "\033[0;36m";
    const std::string ORANGE  = "\033[38;5;208m"; 



    void clearscreen(); //limpa a tela para próximos usos
    void printProgramHeader(const std::string& authorName); //nome do programa e outros detalhes
    void showLoadingAnimation(const std::string& message); //*.*.*.*.*.*. ANIMAÇÕES *.*.*.*.*.*.*
    void showProgressBar(const std::string& message);
    void printMatchResult(const std::string& matchName);
    void printNoMatchFound();

    //impl para DNA (exibição)

    void printHighlightedDna(const std::string& dna_seq, 
        const std::map<size_t, std::pair<std::string, int>>& str_positions);
    void printStrProfileSummary(const std::map<std::string, int>& profile);


    void printError(const std::string& errorMessage);

   //*.*.*.*.*.*. funções: *.*.*.*.*.*.*

    void clearScreen() {
    #ifdef _WIN32
        system("cls");
    #else
        system("clear");
    #endif
    }

    void printProgramHeader(const std::string& authorName) { 
        clearScreen();
        std::cout << BLUE << "Welcome to the C++ DNA Profiler, v1.0\n";
        std::cout << "Copyright (C) 2025, " << authorName << "\n\n" << RESET;
        std::cout << "This program loads a DNA database and an unknown DNA sequence\n";
        std::cout << "and tries to find a match between the input DNA and the database.\n\n";
    }


    void showLoadingAnimation(const std::string& message) {
        std::cout << MAGENTA << "[+] Loading " << message;
        for (int i = 0; i < 3; ++i) {
            std::cout << "." << std::flush;
            std::this_thread::sleep_for(std::chrono::milliseconds(300));
        }
        std::cout << " " << GREEN << "[OK]" << RESET << "\n";
   }

   
    void showProgressBar(const std::string& message) {
        std::cout << MAGENTA << "[+] " << message << RESET << "\n";
        std::cout << "[";
        for (int i = 0; i < 50; ++i) {
            std::cout << "=" << std::flush;
            std::this_thread::sleep_for(std::chrono::milliseconds(50));
        }
        std::cout << "]\n\n";
    }


    void printMatchResult(const std::string& matchName) {
    std::cout << GREEN << "Match ID (99.9%): " << matchName << "\n\n" << RESET;
    }

    void printNoMatchFound() {
        std::cout << RED << "NO MATCH FOUND FOR THE PROVIDED DNA." << RESET << "\n";
    }

    void printHighlightedDna(const std::string& dna_seq, const std::map<std::string, int>& profile) {
        const int context=20; 
        //caracteres (qtd.) para mostrar antes e depois das cadeias da seq. de char
        std::cout<<"Showing up the STRs\n";

            
        for (const auto& [str, count] : profile) {
            if (count == 0) continue;

            size_t pos = findLongestRunPos(dna_seq, str);
            if (pos == std::string::npos) continue;

            // Calcula limites de trechos
            size_t run_length = str.length() * count;
            size_t start = (pos > context) ? pos - context : 0;
            size_t end = (pos + run_length + context < dna_seq.length()) ? pos + run_length + context : dna_seq.length();


            //efeitos
            std::cout << "..." << dna_seq.substr(start, pos - start);
            std::cout << MAGENTA; // Highlight color
            std::cout << dna_seq.substr(pos, run_length);
            std::cout << RESET;
            std::cout << " [x" << count << "]";
            std::cout << dna_seq.substr(pos + run_length, end - (pos + run_length)) << "...\n";
        }
        std::cout << "\n\n\n";
    }


    void printStrProfileSummary(const std::map<std::string, int>& profile) {
        std::cout << "STR Profile Summary:\n";
        for (const auto& [str, count] : profile) {
            std::cout << ORANGE << str << RESET << " [x" << count << "]\t";
        }
        std::cout << "\n\n\n";
    }

    void printError(const std::string& errorMessage) {
        std::cerr << RED << "Error: " << errorMessage << RESET << std::endl;
    }


}


#endif //VISUALS_HPP