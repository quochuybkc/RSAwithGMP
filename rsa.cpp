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

using namespace std;

class rsa{
private:
	mpz_t p; // prime p
	mpz_t q; // prime q
	mpz_t n; // n = p*q
	mpz_t phi; // phi(n) = (p-1)*(q-1)
	mpz_t e; // e,n is public key use encrypt
	mpz_t d; // d,n is private key use decrypt
public:
	rsa();
	~rsa();
	void power(mpz_t, mpz_t, unsigned long int);
	void power(mpz_t, unsigned long int , unsigned long int);
	void mypowm(mpz_t, mpz_t, mpz_t, mpz_t);
	void random(mpz_t, mp_bitcnt_t);
	void random(mpz_t, mpz_t);
	void myrandom(mpz_t, int);
	bool miller_prime_test(mpz_t);
	void find_p_q(mp_bitcnt_t);
	bool genkey(mp_bitcnt_t);
	void gcd(mpz_t,mpz_t, mpz_t);
	void encrypt(mpz_t, mpz_t, mpz_t, mpz_t);
	bool encrypt_str(string, string);
	void decrypt(mpz_t, mpz_t, mpz_t, mpz_t);
	bool decrypt_str(string, string);
};
// constructor
rsa::rsa(void){
	mpz_init(p);mpz_init(q);mpz_init(n);
	mpz_init(phi);mpz_init(e);mpz_init(d);
}

// destructor
rsa::~rsa(void){
	mpz_clear(p);mpz_clear(q);mpz_clear(n);
	mpz_clear(phi);mpz_clear(e);mpz_clear(d);
}

// calculate a = b^n , type b is mpz_t
void rsa::power(mpz_t a, mpz_t b, unsigned long int n){
	mpz_set(a,b);
	unsigned long int i;
	for(i = 1; i < n; i++){
		mpz_mul(a,a,b);
	}
}

// calculate a = b^n, type b is unsigned long int
void rsa::power(mpz_t a, unsigned long int b, unsigned long int n){
	mpz_set_ui(a,b);
	unsigned long int i;
	for(i = 1; i < n; i++){
		mpz_mul_ui(a,a,b);
	}
}

// calculate a = b^c mod n
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
bool rsa::miller_prime_test(mpz_t n){
	// find k, q so that 2^k*q = n-1
	mpz_t q, temp, temp_k;
	unsigned long int k = 0;
	mpz_init(q);
	mpz_init(temp_k);
	mpz_init(temp);
	mpz_sub_ui(temp,n,1);// temp = n -1
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
	mpz_t a;
	mpz_init(a);
	random(a,n);
	while(mpz_cmp_ui(a,2) < 0){
		random(a,n);
	}
	// check res = a^q mod n = 1 ?
	mpz_t res;
	mpz_init(res);
	mpz_powm(res,a,q,n);
	if(mpz_cmp_ui(res,1) == 0){
		mpz_clear(q);mpz_clear(temp_k);mpz_clear(temp);mpz_clear(s);
		mpz_clear(r);mpz_clear(a);mpz_clear(res);
		return true; // maybe prime
	}
	// check a^(2^j*q) mod n = n -1 ?
	mpz_t tj;
	mpz_init(tj);
	unsigned long int j;
	for(j = 0; j < k; j++){
		power(tj,2,j);// tj = 2^j
		mpz_mul(tj,tj,q);// tj = 2^j*q
		mpz_powm(tj,a,tj,n);
		if(mpz_cmp(tj,n-1) == 0){// if a^(2^j*q) mod n = n -1
			mpz_clear(q);mpz_clear(temp_k);mpz_clear(temp);mpz_clear(s);
			mpz_clear(r);mpz_clear(a);mpz_clear(res);mpz_clear(tj);
			return true; // maybe prime
		}
	}
	//
	mpz_clear(q);mpz_clear(temp_k);mpz_clear(temp);mpz_clear(s);
	mpz_clear(r);mpz_clear(a);mpz_clear(res);mpz_clear(tj);
	return false; // composite
}

// find 2 prime number(n bits)
void rsa::find_p_q(mp_bitcnt_t n){
	mpz_t p, q;
	mpz_init(p);mpz_init(q);
	myrandom(p,n);
	myrandom(q,n);
	// loop until find p and q
	while(!miller_prime_test(p)){
		random(p,n);
	}
	while(!miller_prime_test(q)){
		random(q,n);
	}
	mpz_set(this->p, p);
	mpz_set(this->q, q);
	mpz_clear(p);mpz_clear(q);
}

// generate public and private key and save
bool rsa::genkey(mp_bitcnt_t n_bits_key){
	//find 2 prime number p and q
	find_p_q(n_bits_key);
	// n = p*q
	mpz_t n;
	mpz_init(n);
	mpz_mul(n,p,q);
	mpz_set(this->n,n);
	// phi(n) = (p-1)*(q-1)
	mpz_t phi, t_p, t_q;
	mpz_init(phi);mpz_init(t_p);mpz_init(t_q);
	mpz_sub_ui(t_p,p,1);
	mpz_sub_ui(t_q,q,1);
	mpz_mul(phi,t_p,t_q);
	mpz_set(this->phi,phi);
	// generate randome e, 1 < e < phi(n), gcd(e,phi(n)) = 1
	mpz_t e,gcd_e;
	mpz_init(e);
	mpz_init(gcd_e);
	random(e,phi);
	gcd(gcd_e,e,phi);
	while((mpz_cmp_ui(gcd_e,1) != 0) || (mpz_cmp_ui(e,2) < 0)){	
		random(e,phi);
		gcd(gcd_e,e,phi);
	}
	mpz_set(this->e,e);// set e
	// find d from e, e*d mod phi(n) = 1, 0 <= d <= phi(n)
	// e*d + k * phi(n) = 1 , gcd(e,phi(n)) = 1 => extend euclidean algorithm
	mpz_t d, k, g;
	mpz_init(d);mpz_init(k);mpz_init(g);
	mpz_set_ui(g,1); // set gcd of phi(n) and e is 1
	mpz_gcdext(g,d,k,e,phi);// find d and k satifies: e*d + k * phi(n) = 1 
	mpz_set(this->d,d);// set d
	// write public key and private key in a file
	ofstream pub, pvt;
	pub.open("pub_key.txt");
	pvt.open("pvt_key.txt");
	pub << e << " " << n << endl;
	pvt << d << " " << n << endl;
	pub.close();
	pvt.close();
	mpz_clear(n);mpz_clear(phi);mpz_clear(t_p);mpz_clear(t_q);
	mpz_clear(e);mpz_clear(gcd_e);mpz_clear(d);mpz_clear(k);mpz_clear(g);
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

// random range: 2^(n-1) to 2^n -1, not use random func in gmp
void rsa::myrandom(mpz_t result, int n){
	mpz_t base, off;
	mpz_init(base);mpz_init(off);
	power(base,2,n-1); // base = 2^(n-1)
	unsigned int tmp1, tmp2;
	tmp1 = rand()%31;
	tmp2 = rand();
	power(off,tmp2,tmp1);
	mpz_add(result,base,off);
	mpz_clear(base);mpz_clear(off);
}

// encrypt string
bool rsa::encrypt_str(string plain_text_file, string pub_key_file){
	ofstream enc;
	ifstream pub_file, pl_file;
	enc.open("enc.txt");
	// open public key file
	pub_file.open(pub_key_file);
	// open plain text file
	pl_file.open(plain_text_file);
	if(!pub_file || !pl_file)
		return false;
	mpz_t str_data, enc_data, e, n;
	mpz_init(str_data);mpz_init(enc_data);
	mpz_init(e);mpz_init(n);
	// read pub key
	pub_file >> e >> n;
	// read plain text
	string pl, tmp;
	int i, len;
	while(pl_file >> tmp){
		pl +=  tmp;
		pl += " ";
	}
	// encrypt plain text
	i = 0;
	len = pl.length();
	while(i != len){
		mpz_set_ui(str_data,pl[i]);
		encrypt(enc_data, str_data, e, n);
		enc << enc_data << " ";
		i++;
	}
	enc << endl;
	mpz_clear(str_data);mpz_clear(enc_data);
	mpz_clear(e);mpz_clear(n);
	enc.close();
	pub_file.close();
	pl_file.close();
	return true;
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
	mpz_t tmp, pt, d, n;
	mpz_init(tmp);mpz_init(pt);
	mpz_init(d);mpz_init(n);
	int t;
	string plain_text;
	//read pvt key
	ofstream dec_file;
	dec_file.open("plaintext.txt");
	pvt_file >> d >> n;
	while(dec >> tmp){	
		decrypt(pt, tmp, d, n);
		t = mpz_get_ui(pt);
		plain_text += (char)t;
	}
	dec_file << plain_text << endl;
	dec.close();
	dec_file.close();
	pvt_file.close();
	mpz_clear(tmp);mpz_clear(pt);
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