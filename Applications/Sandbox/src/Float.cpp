#include "Float.h"
#include "NKLogger/NkLog.h"
#include "NKContainers/String/NkStringUtils.h"
#include <bitset>

using namespace nkentseu::string;

// Sommation de Kahan — O(N) comme std::accumulate mais BEAUCOUP plus précis 
float NkMath::kahanSum(std::vector<double>& data) {
    float sum  = 0.0f; 
    float comp = 0.0f;  // compensation des erreurs perdues 
 
    for(int i = 0; i < data.size(); i++) { 
        float y = data[i] - comp;    // compenser l'erreur précédente 
        float t = sum + y;           // t est grand, y est petit → perte de bits 
        comp    = (t - sum) - y;     // capture les bits perdus dans y 
        sum     = t; 
    } 
    return sum; 
}
float NkMath::kahanSum(std::vector<float>& data) {
    float sum  = 0.0f; 
    float comp = 0.0f;  // compensation des erreurs perdues 
 
    for(int i = 0; i < data.size(); i++) { 
        float y = data[i] - comp;    // compenser l'erreur précédente 
        float t = sum + y;           // t est grand, y est petit → perte de bits 
        comp    = (t - sum) - y;     // capture les bits perdus dans y 
        sum     = t; 
    } 
    return sum; 
}

// Inspecter un float, écrit en binaire
// Norme IEEE 754 : valeur = (-1)^signe × 1.mantisse × 2^(exposant - 127)
void NkMath::inspectFloat(float x) {
    // Interpréter les bits du float comme un entier 32 bits
    uint32_t bits;
    std::memcpy(&bits, &x, sizeof(bits));
    
    // Utiliser bitset pour la représentation binaire
    std::bitset<32> binary(bits);
    
    // Extraire les champs pour l'affichage séparé
    uint32_t sign = (bits >> 31) & 0x1;
    uint32_t exponent = (bits >> 23) & 0xFF;
    uint32_t mantissa = bits & 0x7FFFFF;
    
    // Affichage détaillé
    logger.Info("Inspection du float {0} :\n"
                "Représentation binaire complète : {1}\n"
                "Signe (1 bit)     : {2} ({3})\n"
                "Exposant (8 bits) : {4} (valeur : {5})\n"
                "Mantisse (23 bits): {6}\n"
                "Valeur décodée    : (-1)^{7} × 1.{8} × 2^({9} - 127)",
        x,
        binary.to_string(),
        sign, (sign == 0 ? "positif" : "négatif"),
        std::bitset<8>(exponent).to_string(), exponent,
        std::bitset<23>(mantissa).to_string(),
        sign,
        mantissa,
        exponent
    );
}
// Norme IEEE 754 : valeur = (-1)^s × 1.mantisse × 2^(exposant - 1023)
void NkMath::inspectDouble(double x) {
    // Interpréter les bits du double comme un entier 64 bits
    uint64_t bits;
    std::memcpy(&bits, &x, sizeof(bits));

    uint64_t sign = (bits >> 63) & 0x1;         // Le dernier bit après décalage
    uint64_t exponent = (bits >> 52) & 0x7FF;   // Les 11 derniers bits après décalage
    uint64_t mantissa = bits & 0xFFFFFFFFFFFFF; // Les 52 derniers bits après décalage
    
    // Affichage en binaire
    logger.Info("Affichage du double {0} en format binaire :\n- Signe    : {1}\n- Exposant : {2}\n- Mantisse : {3}",
        x,
        sign,        
        exponent,
        mantissa
    );
}

// Variance naïve : moyenne des carrés - carré de la moyenne
float NkMath::varianceNaive(const std::vector<float>& data) {
    float sum = 0.0f, sumSq = 0.0f;

    for (float x : data) {
        sum += x;
        sumSq += x * x;
    }

    float mean = sum / data.size();
    return (sumSq / data.size()) - (mean * mean);
}
// Formule de Welford
float NkMath::varianceWelford(const std::vector<float>& data) {
    float mean = 0.0f;
    float M2 = 0.0f;
    int n = 0;

    for (float x : data) {
        n++;
        float delta = x - mean;
        mean += delta / n;
        float delta2 = x - mean;
        M2 += delta * delta2;
    }

    return M2 / n;
}

// Epsilon machine par boucle
float NkMath::epsilonMachine() {
    float eps = 1.0f;

    while ((1.0f + eps / 2.0f) > 1.0f) {
        eps /= 2.0f;
    }

    return eps;
}