#include <string>
#include <cmath>
using namespace std;

std::string amino_check(int num); 

int check_amino(char amino){
    int result;
    
    switch (amino) {
		case 'A':
            result = 0;
			break;
		case 'C':
			result = 1;
            break;
        case 'D':
			result = 2;
            break;
        case 'E':
			result = 3;
            break;
        case 'F':
			result = 4;
            break;
        case 'G':
			result = 5;
            break;
        case 'H':
			result = 6;
            break;
        case 'I':
			result = 7;
            break;
        case 'K':
			result = 8;
            break;
        case 'L':
			result = 9;
            break;
        case 'M':
			result = 10;   
            break;
        case 'N':
			result = 11;
            break;
        case 'P':
			result = 12;
            break;
        case 'Q':
			result = 13;
            break;
        case 'R':
			result = 14;
            break;
        case 'S':
			result = 15;
            break;
        case 'T':
			result = 16;
            break;
        case 'U':
			result = 17;
            break;
        case 'V':
			result = 18;
            break;
        case 'W':
			result = 19;
            break;
        case 'X':
			result = 20;
            break;
        case 'Y':
			result = 21;
            break;
        default:
			break;
	}
    return result;
}
int transform_scs_10(int n, string amino5){
    int result = 0;      
    int multiplier = n-1; 
    for (int i = 0; i < n; i++){
        char amino = amino5[i];

        switch (amino) {
		case 'A':
			break;
		case 'C':
			result += 1*std::pow(20, multiplier);   
            break;
        case 'D':
			result += 2*std::pow(20, multiplier);   
            break;
        case 'E':
			result += 3*std::pow(20, multiplier);   
            break;
        case 'F':
			result += 4*std::pow(20, multiplier);   
            break;
        case 'G':
			result += 5*std::pow(20, multiplier);   
            break;
        case 'H':
			result += 6*std::pow(20, multiplier);   
            break;
        case 'I':
			result += 7*std::pow(20, multiplier);   
            break;
        case 'K':
			result += 8*std::pow(20, multiplier);   
            break;
        case 'L':
			result += 9*std::pow(20, multiplier);   
            break;
        case 'M':
			result += 10*std::pow(20, multiplier);   
            break;
        case 'N':
			result += 11*std::pow(20, multiplier);   
            break;
        case 'P':
			result += 12*std::pow(20, multiplier);   
            break;
        case 'Q':
			result += 13*std::pow(20, multiplier);   
            break;
        case 'R':
			result += 14*std::pow(20, multiplier);   
            break;
        case 'S':
			result += 15*std::pow(20, multiplier);   
            break;
        case 'T':
			result += 16*std::pow(20, multiplier);   
            break;
        case 'U':
			result += 17*std::pow(20, multiplier);   
            break;
        case 'V':
			result += 17*std::pow(20, multiplier);   
            break;
        case 'W':
			result += 18*std::pow(20, multiplier);   
            break;
        case 'X':
			result += 20*std::pow(20, multiplier);   
            break;
        case 'Y':
			result += 19*std::pow(20, multiplier);   
            break;
        default:
			break;
	}
    multiplier--;
    }
    return result;
}
std::string transform_10_scs(int element, int n){ 
    std::string result; 

    switch (n) {
	case 3:
        if (element / 22 < 1){
            result = "AA" + amino_check(element);
        }else if (element / 484 < 1){
            result = "A" + amino_check(element / 22);
            result += amino_check(element % 22);
        }else{
            result = amino_check(element / 484) + amino_check((element % 484) / 22);;
            result += amino_check((element % 22) % 22);
        }
		break;
    case 4:
        if (element / 22 < 1){
            result = "AAA" + amino_check(element);
        }else if (element / 484 < 1){
            result = "AA" + amino_check(element / 22);
            result += amino_check(element % 22);
        }else if (element / 10648 < 1){
            result = "A" + amino_check(element / 484);
            result += amino_check((element % 484) / 22);
            result += amino_check((element % 22) % 22);
        }else{
            result = amino_check(element / 10648) + amino_check((element % 10648) / 484);
            result += amino_check((element % 484) / 22);
            result += amino_check((element % 22) % 22);
        }
        break;
    case 5:
        if (element / 22 < 1){
            result = "AAAA" + amino_check(element);
        }else if (element / 484 < 1){
            result = "AAA" + amino_check(element / 22);
            result += amino_check(element % 22);
        }else if (element / 10648 < 1){
            result = "AA" + amino_check(element / 484);
            result += amino_check((element % 484) / 22);
            result += amino_check((element % 22) % 22);
        }else if (element / 234256 < 1){
            result = "A" + amino_check(element / 10648);
            result += amino_check((element % 10648) / 484);
            result += amino_check((element % 484) / 22);
            result += amino_check((element % 22) % 22);
        }else{
            result = amino_check(element / 234256) + amino_check((element % 234256) / 10648);
            result += amino_check((element % 10648) / 484);
            result += amino_check((element % 484) / 22);
            result += amino_check((element % 22) % 22);
        }
    }
    return result;
}
std::string amino_check(int num){
    std::string result;
    
    switch (num) {
	case 0:
        result = "A";
		break;
	case 1:
        result = "C";
		break;
    case 2:
        result = "D";
		break;
    case 3:
        result = "E";
		break;
    case 4:
        result = "F";
		break;
    case 5:
        result = "G";
		break;
    case 6:
        result = "H";
		break;
    case 7:
        result = "I";
		break;
    case 8:
        result = "K";
		break;
    case 9:
        result = "L";
		break;
    case 10:
        result = "M";
		break;
    case 11:
        result = "N";
		break;
    case 12:
        result = "P";
		break;
    case 13:
        result = "Q";
		break;
    case 14:
        result = "R";
		break;
    case 15:
        result = "S";
		break;
    case 16:
        result = "T";
		break;
    case 17:
        result = "U";
		break;
    case 18:
        result = "V";
		break;
    case 19:
        result = "W";
		break;
    case 20:
        result = "X";
		break;
    case 21:
        result = "Y";
		break;
    default:
		break;
	}
    return result;
}