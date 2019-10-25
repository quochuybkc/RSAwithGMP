////////////////////////////////////////////
// ASSIGMENT MAT MA & AN NINH MANG
// 
// IMPLEMENT RSA C++ (512 bits key)
//
// USE LIBRARY: GNU MP
//
// quochuy
////////////////////////////////////////////
#include <iostream>
#include <gmp.h>
#include <string>
#include <fstream>
#include <gmpxx.h>
#include <cmath>

using namespace std;

class rsa{
private:
	mpz_class p; // prime p
	mpz_class q; // prime q
public:
	rsa();
	~rsa();
	void power(mpz_t, mpz_t, unsigned long int);
	void power(mpz_t, unsigned long int , unsigned long int);
	void mypowm(mpz_t, mpz_t, mpz_t, mpz_t);
	void random(mpz_t, mp_bitcnt_t);
	void random(mpz_t, mpz_t);
	bool miller_prime_test(mpz_class);
	void find_p_q(mp_bitcnt_t);
	bool genkey(mp_bitcnt_t);
	void gcd(mpz_t,mpz_t, mpz_t);
	void encrypt(mpz_t, mpz_t, mpz_t, mpz_t);
	bool encrypt_str(string, string);
	void decrypt(mpz_t, mpz_t, mpz_t, mpz_t);
	bool decrypt_str(string, string);
	mpz_class mygcd(mpz_class, mpz_class);
	mpz_class mygcd_ext(mpz_class, mpz_class);
	mpz_class string2num(string);
};
// constructor
rsa::rsa(void){
	cout<<"Executing...\n";
}

// destructor
rsa::~rsa(void){
	cout<<"Done\n";
}

// calculate a = b^n , type b is mpz_t, O(n)
void rsa::power(mpz_t a, mpz_t b, unsigned long int n){
	mpz_set(a,b);
	unsigned long int i;
	for(i = 1; i < n; i++){
		mpz_mul(a,a,b);
	}
}

// calculate a = b^n, type b is unsigned long int, O(n)
void rsa::power(mpz_t a, unsigned long int b, unsigned long int n){
	mpz_set_ui(a,b);
	unsigned long int i;
	for(i = 1; i < n; i++){
		mpz_mul_ui(a,a,b);
	}
}

// cal gcd(a,b)
mpz_class rsa::mygcd(mpz_class a, mpz_class b){
	mpz_class tmp;
    while(b != 0) {
        tmp = a % b;
        a = b;
        b = tmp;
    }
    return a;
}

// extend gcd, use to find d: ed mod phi(n) = 1 or ed + xphi(n) = 1
mpz_class rsa::mygcd_ext(mpz_class e, mpz_class phi){
	mpz_class s, t, r, old_s, old_t, old_r, q, t1, t2, t3;
	s = 0;
	t = 1;
	r = e;
	old_s = 1;
	old_t = 0;
	old_r = phi;
	while(1){
		if(r == 1){
			return t;
		}
		q = old_r / r;
		t1 = old_r - q*r;
		t2 = old_s - q*s;
		t3 = old_t - q*t;
		old_r = r;
		old_s = s;
		old_t = t;
		r = t1;
		s = t2;
		t = t3;
	}
}

// calculate a = b^c mod n, O(log2(n)), square and multiply algorithm
void rsa::mypowm(mpz_t a, mpz_t b, mpz_t c, mpz_t n){
	mpz_mod(b,b,n);// b = b mod n
	mpz_set_ui(a,1);
	mpz_t tmp, t_c;
	mpz_init(tmp);mpz_init(t_c);
	while(mpz_cmp_ui(c,0) > 0){
		mpz_set_ui(tmp,1);
		mpz_and(t_c,c,tmp);//c = c & 1
		if(mpz_cmp_ui(t_c,1) == 0){// c == 1
			mpz_mul(tmp,a,b); // tmp = a*b
			mpz_mod(tmp,tmp,n);// tmp =a*b mod n
			mpz_set(a,tmp);
		}
		mpz_mul(tmp,b,b); //tmp = b*b
		mpz_mod(tmp,tmp,n);// tmp = b*b mod n
		mpz_set(b,tmp);
		mpz_fdiv_q_2exp(c,c,1);// c = c>>1
	}
	mpz_clear(tmp);mpz_clear(t_c);
}

// test prime with miller rabin algorithm
bool rsa::miller_prime_test(mpz_class n){
	// find k, q so that 2^k*q = n-1
	mpz_t q, temp, temp_k;
	unsigned long int k = 0;
	mpz_init(q);
	mpz_init(temp_k);
	mpz_init(temp);
	mpz_sub_ui(temp,n.get_mpz_t(),1);// temp = n -1
	mpz_set(temp_k, temp);// k = n - 1
	mpz_t s, r; // k = s*2 + r;
	mpz_init(s);
	mpz_init(r);
	mpz_cdiv_qr_ui(s,r,temp_k,2);
	while(mpz_cmp_ui(r, 0) == 0){// r = 0
		mpz_set(temp_k, s);
		mpz_cdiv_qr_ui(s,r,temp_k,2);
		k++;
	}
	mpz_set(q,temp_k);
	// random a integer a: 1 < a < n-1
	mpz_class a;
	random(a.get_mpz_t(),n.get_mpz_t());
	while(a < 2){
		random(a.get_mpz_t(),n.get_mpz_t());
	}
	// check res = a^q mod n = 1 ?
	mpz_t res;
	mpz_init(res);
	mpz_powm(res,a.get_mpz_t(),q,n.get_mpz_t());
	if(mpz_cmp_ui(res,1) == 0){
		mpz_clear(q);mpz_clear(temp_k);mpz_clear(temp);mpz_clear(s);
		mpz_clear(r);mpz_clear(res);
		return true; // maybe prime
	}
	// check a^(2^j*q) mod n = n -1 ?
	mpz_t tj;
	mpz_init(tj);
	unsigned long int j;
	for(j = 0; j < k; j++){
		power(tj,2,j);// tj = 2^j
		mpz_mul(tj,tj,q);// tj = 2^j*q
		mpz_powm(tj,a.get_mpz_t(),tj,n.get_mpz_t());
		if(mpz_cmp(tj,n.get_mpz_t() - 1) == 0){// if a^(2^j*q) mod n = n -1
			mpz_clear(q);mpz_clear(temp_k);mpz_clear(temp);mpz_clear(s);
			mpz_clear(r);mpz_clear(res);mpz_clear(tj);
			return true; // maybe prime
		}
	}
	//
	mpz_clear(q);mpz_clear(temp_k);mpz_clear(temp);mpz_clear(s);
	mpz_clear(r);mpz_clear(res);mpz_clear(tj);
	return false; // composite
}

// find 2 prime number(n bits)
void rsa::find_p_q(mp_bitcnt_t n){
	mpz_class p, q;
	random(p.get_mpz_t(),n);
	random(q.get_mpz_t(),n);
	// loop until find p and q
	while(!miller_prime_test(p)){
		random(p.get_mpz_t(),n);
	}
	while(!miller_prime_test(q)){
		random(q.get_mpz_t(),n);
	}
	this->p = p;
	this->q = q;
}

// generate public and private key and save
bool rsa::genkey(mp_bitcnt_t n_bits_key){
	//find 2 prime number p and q
	find_p_q(n_bits_key);
	// n = p*q
	mpz_class n;
	n = p*q;
	// phi(n) = (p-1)*(q-1)
	mpz_class phi;
	phi = (p - 1)*(q - 1);
	// generate randome e, 1 < e < phi(n), gcd(e,phi(n)) = 1
	mpz_class e,gcd_e;
	random(e.get_mpz_t(),phi.get_mpz_t());
	gcd_e = mygcd(e,phi);
	while(gcd_e != 1 || e < 0){	
		random(e.get_mpz_t(),phi.get_mpz_t());
		gcd_e = mygcd(e,phi);
	}
	// find d from e, e*d mod phi(n) = 1, 0 <= d <= phi(n)
	// e*d + k * phi(n) = 1 , gcd(e,phi(n)) = 1 => extend euclidean algorithm
	mpz_class d;
	d = mygcd_ext(e,phi); 
	// write public key and private key in a file
	ofstream pub, pvt;
	pub.open("pub_key.txt");
	pvt.open("pvt_key.txt");
	pub << e << " " << n << endl;
	pvt << d << " " << n << endl;
	pub.close();
	pvt.close();
	return true;
}

// greatest commom divisor of a and b
void rsa::gcd(mpz_t r, mpz_t a, mpz_t b){
	mpz_gcd(r,a,b);
}

// encrypt, C is cipher text, P is plain text
void rsa::encrypt(mpz_t C, mpz_t P, mpz_t e, mpz_t n){
	mpz_powm(C,P,e,n);
}

//decrypt, D is decrypt text, C is cipher text
void rsa::decrypt(mpz_t D, mpz_t C, mpz_t d, mpz_t n){
	mpz_powm(D,C,d,n);
}

// random range: 2^(n-1) to 2^n -1, use random func in gmp
void rsa::random(mpz_t result, mp_bitcnt_t n){
	unsigned long int seed;
	seed = rand();
	gmp_randstate_t state;
	gmp_randinit_mt(state);
	gmp_randseed_ui(state, seed);
	mpz_rrandomb(result, state, n);
	gmp_randclear(state);
}

// random range: 0 to n-1, use random func in gmp
void rsa::random(mpz_t result, mpz_t n){
	unsigned long int seed;
	seed = rand();
	gmp_randstate_t state;
	gmp_randinit_mt(state);
	gmp_randseed_ui(state, seed);
	mpz_urandomm(result, state, n);
	gmp_randclear(state);
}

// encrypt string
bool rsa::encrypt_str(string plain_text_file, string pub_key_file){
	ofstream enc;
	ifstream pub_file;
	ifstream pl_file(plain_text_file, ios::binary | ios::ate);
	enc.open("enc.txt");
	// open public key file
	pub_file.open(pub_key_file);
	// open plain text file
	if(!pub_file || !pl_file)
		return false;
	mpz_t enc_data, e, n;
	mpz_init(enc_data);
	mpz_init(e);mpz_init(n);
	// read pub key
	pub_file >> e >> n;
	// encrypt plain text, 4bytes
	int f_size = pl_file.tellg(); // length of file
	f_size--;
    char bytes[4];
    string tmp;
    mpz_class M;
    int i, k;
    k = 0;
    pl_file.seekg(0); // set read point
    while(k <= f_size/4){
    	i = 0;
    	if(k == f_size/4 ){
    		pl_file.read(bytes, f_size%4);    		
    	}
    	else{
	    	pl_file.read(bytes, 4);   
		}
	    while(i < 4){
	    	tmp += bytes[i];
	    	i++;
    	} 
    	k++;
    	pl_file.seekg(4*k);
    	for(int j = 0; j < 4;j++)
    		bytes[j] = 32;
    	M = string2num(tmp);// 4 bytes plain text
    	encrypt(enc_data, M.get_mpz_t(), e, n);
    	enc << enc_data << " ";
    	tmp = "";
	}
	enc << endl;
	mpz_clear(enc_data);
	mpz_clear(e);mpz_clear(n);
	enc.close();
	pub_file.close();
	pl_file.close();
	return true;
}

// string to number, e.g: abc -> 197198199
mpz_class rsa::string2num(string in){
	mpz_class tmp;
	mpz_class tmp_data;
	int i = 0;
	while(i < in.length()){
		if(in[i] == 0)
			break;
		tmp = in[i];
		tmp += 100;
		tmp_data = tmp_data*1000 + tmp;
		i++;
	}
	return tmp_data;
}

//decrypt string
bool rsa::decrypt_str(string enc_file, string pvt_key_file){
	ifstream dec;
	ifstream pvt_file;
	// open pvt key file
	pvt_file.open(pvt_key_file);
	dec.open(enc_file);
	if(!pvt_file || !dec)
		return false;
	mpz_class tmp, num;
	mpz_t d, n;
	mpz_init(d);mpz_init(n);
	//read pvt key
	ofstream dec_file;
	dec_file.open("plaintext.txt");
	pvt_file >> d >> n;
	// decrypt file
	int t;
	int i;
	string plain_text;
	while(dec >> tmp){	
		decrypt(num.get_mpz_t(), tmp.get_mpz_t(), d, n);
		i = 1;
		while(i < 5){
			tmp = num / pow(10,12 - 3*i);
			num = num % pow(10,12 - 3*i);
			tmp -= 100;
			t = mpz_get_ui(tmp.get_mpz_t());
			plain_text += (char)t;
			i++;
		}
		dec_file << plain_text;
		plain_text = "";
	}
	dec_file << endl;
	dec.close();
	dec_file.close();
	pvt_file.close();
	mpz_clear(d);mpz_clear(n);
	return true;
}

// main function
int main(int argc, char **argv){
	if(argc == 1){
		cout << "Few arguments\n" << "Usage: ./RSA arguments\n";
		cout << "Arguments: genkey or encrypt plaintextfile pubkeyfile or decrypt encryptfile pvtkeyfile\n";
		cout << "Attention: must be generate key before encrypt or decrypt\n";
		return 1;
	}
	if(argc > 4){
		cout << "Too much arguments\n" << "Usage: ./RSA arguments\n";
		cout << "Arguments: genkey or encrypt plaintextfile pubkeyfile or decrypt encryptfile pvtkeyfile\n";
		cout << "Attention: must be generate key before encrypt or decrypt\n";
		return 1;
	}
	if(argc == 2 && argv[1] == string("genkey")){
		srand(time(NULL));
		rsa RSA;
		if(RSA.genkey(512)) 
			cout << "Generate key successfully. Check at pub_key.txt and pvt_key.txt\n";
		else
			cout << "Can't generate keypair\n";
	}
	else if(argc == 4 && argv[1] == string("encrypt")){
		rsa RSA;
		if(!RSA.encrypt_str(argv[2],argv[3])){
			cout << "Try again.Can't read file\n";
		}
		else{
			RSA.encrypt_str(argv[2],argv[3]);
		}
	}
	else if(argc == 4 && argv[1] == string("decrypt")){
		rsa RSA;
		if(!RSA.decrypt_str(argv[2],argv[3])){
			cout << "Try again.Can't read file\n";
		}
		else{
			RSA.decrypt_str(argv[2],argv[3]);
		}
	}
	else if(argc == 2 && argv[1] == string("--help")){
		cout << "Usage: ./RSA arguments\n\n";
		cout << "Arguments: genkey or encrypt plaintextfile pubkeyfile or decrypt encryptfile pvtkeyfile\n\n";
		cout << "Example:\n\tGenerate keypair: ./RSA genkey\n\tEncrpt: ./RSA encrypt testRSA.txt pub_key.txt\n\tDecrypt: ./RSA decrypt enc.txt pvt_key.txt\n\n";
		cout << "Attention: must be generate key before encrypt or decrypt\n";
	}
	else{
		cout << "Unknow command.Try again\n";
		cout << "Use: ./RSA --help for details\n";
	}
	return 0;
}