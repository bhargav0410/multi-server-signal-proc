//FFTW library 
#include <fftw3.h>
//Shared Memory 
#include "CSharedMemSimple.hpp"
#include "ShMemSymBuff.hpp"
#include <csignal>
#include <fstream>
#include <boost/asio.hpp>
#include <boost/lexical_cast.hpp>
#include <sys/socket.h>
#define mode 0
#define fileNameForX "Pilots.dat"
#ifndef addr
	#define addr "10.10.21.1"
#endif

#ifndef remote_addr
	#define remote_addr "10.10.0.10"
#endif

#ifndef port
	#define port 1001
#endif

/*
	mode:
		= 1 -> master -> creates shared memory 
		= 0 -> slave -> doesn't create the shared memory
		
	Waits to read dimension vector then does fft on it and then divides by 1+i 
*/

//! Install dependencies: apt-get -y install libboost-program-options-dev libfftw3-dev 
//!How to Compile:   g++ -o cpu ../../examples/cpuLS.cpp -lfftw3f -lrt
// ./cpu

//LS
//Y = 16 x 1024
//X = 1 x 1023
//H = 16 x 1023
ShMemSymBuff* buffPtr;

using namespace std;
using boost::asio::ip::udp;

std::string file = "Output_cpu.dat";
//std::ofstream outfile;
/*
if (not file.empty()) {
	outfile.open(file.c_str(), std::ofstream::binary);
}
*/

static bool stop_signal_called = false;
void sig_int_handler(int){stop_signal_called = true;}

//Reads in Vector X from file -> 1xcols
void matrix_readX(complexF* X, int cols){
	ifstream inFile;
	inFile.open(fileNameForX, std::ifstream::binary);
	if (!inFile) {
		cerr << "Unable to open file data file, filling in 1+i for x\n";
		float c=1.0f;
		for (int col = 0; col <  cols; col++){
			X[col].real=c;
			X[col].imag=c;
		}
		return;
	}
	
	inFile.read((char*)X, (cols)*sizeof(*X));
	//printOutArr(X, 1, cols);
	/*
	float c=0;
	for (int col = 0; col <  cols; col++){
		inFile >> c;
		X[col].real=c;
		inFile >> c;
		X[col].imag=c;
	}
	*/
	
	complexF* temp = 0;
	temp=(complexF*)malloc ((cols-1)/2* sizeof (*temp));
	//copy second half to temp
	memmove(temp, &X[(cols+1)/2], (cols-1)/2* sizeof (*X));
	//copy first half to second half
	memmove(&X[(cols-1)/2], X, (cols+1)/2* sizeof (*X));
	//copy temp to first half
	memmove(X, temp, (cols-1)/2* sizeof (*X));
	
	free(temp);
	
	inFile.close();
}

//Shifts first and second half of vector fft(Y)
void shiftOneRow(complexF* Y, int cols, int row){
	complexF* Yf = &Y[row*cols];
	//std::cout << "Here...\n";
	complexF* temp = 0;
	temp=(complexF*)malloc ((cols+1)/2* sizeof (*temp));
	//copy second half to temp
	memmove(temp, &Yf[(cols-1)/2], (cols+1)/2* sizeof (*Yf));
	//copy first half to second half
	memmove(&Yf[(cols+1)/2], Yf, (cols-1)/2* sizeof (*Yf));
	//copy temp to first half
	memmove(Yf, temp, (cols+1)/2* sizeof (*Yf));
	
	free(temp);
	
}

//FFT on one vector of Y in row row
void fftOneRow(complexF* Y, int cols, int row){
	
	fftwf_complex* org = (fftwf_complex*)&Y[row*cols];
	fftwf_complex* after = (fftwf_complex*)&Y[row*cols];

	fftwf_plan fft_p = fftwf_plan_dft_1d(cols, org, after, FFTW_FORWARD, FFTW_MEASURE /*FFTW_ESTIMATE*/);
	fftwf_execute(fft_p);
	fftwf_destroy_plan(fft_p);
	
	
}

//Element by element multiplication
void matrixMultThenSum(complexF* Y, complexF* Hconj, complexF* Yf, int rows, int cols ){
	//Y x conj(H) -> then sum all rows into elements in Yf
	//Y = 16x1023
	//conjH = 16x1023
	for (int i = 0; i<rows; i++){
		for(int j=0; j<cols-1; j++){
			float Yreal = Y[i*(cols-1)+j].real;
			float Yimag = Y[i*(cols-1)+j].imag;
			float Hreal = Hconj[i*(cols-1)+j].real;
			float Himag = Hconj[i*(cols-1)+j].imag;
			//(a+bi)(c+di) = a*c - b*d + (bc + ad)i
			if(i==0){
				Yf[j].real = 0;
				Yf[j].imag=0;
			}
			
			Yf[j].real=Yf[j].real+(Yreal*Hreal - Yimag*Himag);
			Yf[j].imag=Yf[j].imag+(Yreal*Himag + Yimag*Hreal);	
		}
	}
	
} 

// Gives Sum of Hsqrd = |H|^2 as a 1x1023 vector
void findDistSqrd(complexF* H, complexF* Hsqrd, int rows, int cols){
	//initialize first row since Hsqrd currently holds X
	for (int j = 0; j<cols; j++){
		//|H|^2 = real^2 + imag^2
		//Sum of |H|^2 is summing all elements in col j
		Hsqrd[j].real = (H[j].real*H[j].real)+ (H[j].imag*H[j].imag);
		Hsqrd[j].imag =0;
	}
	
	for (int i = 1; i<rows; i++){  
		for (int j = 0; j<cols; j++){
			//|H|^2 = real^2 + imag^2
			//Sum of |H|^2 is summing all elements in col j
			Hsqrd[j].real = Hsqrd[j].real+ (H[i*cols + j].real*H[i*cols + j].real)+ (H[i*cols + j].imag*H[i*cols + j].imag);
		}
	}
	
}



//Divide matrix A by vector B -> element division
void divideOneRow(complexF * A, complexF * B, int cols, int row){
	int i=row;
	for(int j=0; j<cols; j++){
		float fxa = A[i*cols+j].real;
		float fxb = A[i*cols+j].imag;
		float fya = B[j].real;
		float fyb = B[j].imag;
		A[i*cols+j].real=((fxa*fya + fxb*fyb)/(fya*fya+fyb*fyb));
		A[i*cols+j].imag=((fxb*fya - fxa*fyb)/(fya*fya + fyb*fyb));	
	}
	
}

void doOneSymbol(complexF* Y, complexF* Hconj, complexF* Hsqrd,int rows, int cols, int it, udp::socket &socket, udp::endpoint &server_endpoint){
	
	if(it==numberOfSymbolsToTest-1){
		//read in 16x1024
		buffPtr->readLastSymbol(Y);
	}
	else{
		//read in 16x1024
		buffPtr->readNextSymbol(Y, it);
	}
	//printOutArr(Y, 1, cols);
	complexF* Yf = 0;
	Yf = (complexF*)malloc((cols-1)* sizeof (*Yf));
	complexF* Ytemp = 0;
	Ytemp = (complexF*)malloc(rows*(cols-1)*sizeof(*Ytemp));
	
	
	
	clock_t start, finish;
	if(timerEn){
		start = clock();
	}
	
	for(int row=0; row<rows; row++){
		//FFT one row 
		fftOneRow(Y, cols, row);
	}
	if(timerEn){
		finish = clock();
		fft[it] = ((float)(finish - start))/(float)CLOCKS_PER_SEC;
	}
	
	if(timerEn){
		start = clock();
	}
	for(int row=0; row<rows; row++){
		memcpy(&Ytemp[row*(cols-1)], &Y[row*cols+1], (cols-1)* sizeof (*Y));
		//shiftOneRow(Ytemp, cols-1, row);
	}
	
	//Find sum of YH* -> 1x1023
	matrixMultThenSum(Ytemp,Hconj,Yf, rows, cols);
	/*
	std::string Sym;
	for (int j = 0; j < cols-1; j++) {
		Sym += boost::lexical_cast<std::string>(Yf[j].real);
		Sym += boost::lexical_cast<std::string>(Yf[j].imag);
	}
	*/
//	std::cout << "Symbol:" << it << std::endl;
	int len = 0;
	while (len == 0 and not stop_signal_called) {
		len = socket.send_to(boost::asio::buffer((char*)&Yf[0], (cols-1)*sizeof(*Yf)), server_endpoint);
	}
	/*
	//Divide YH* / |H|^2
	//divideOneRow(Yf, Hsqrd, cols-1, 0);
	for (int j = 0; j < cols - 1; j++) {
		Yf[j].real = Yf[j].real/Hsqrd[j].real;
		Yf[j].imag = Yf[j].imag/Hsqrd[j].real;
	}
	shiftOneRow(Yf, cols-1, 0);
	if(timerEn){
		finish = clock();
		decode[it] = ((float)(finish - start))/(float)CLOCKS_PER_SEC;
	}
	
	if (it <= 1) {
		outfile.open(file.c_str(), std::ofstream::binary | std::ofstream::trunc);
	} else {
		outfile.open(file.c_str(), std::ofstream::binary | std::ofstream::app);
	}
	outfile.write((const char*)Yf, (cols-1)*sizeof(*Yf));
	outfile.close();
	*/
	/*
	if(testEn){
		printOutArr(Yf, 1, cols-1);
	}
	*/
	
	free(Yf);
}

//Finds |H|^2 and H*=Hconj, rows=16 cols=1024
void firstVector(complexF* Y, complexF* Hconj, complexF* X, int rows, int cols, udp::socket &socket, udp::endpoint &server_endpoint){
	//Read in X vector -> 1x1023
	matrix_readX(X, cols-1);
	//printOutArr(X, 1, cols-1);
	for (int i = 0; i<rows; i++){  
		for (int j = 0; j<cols; j++){
			Y[i*cols + j].real=0;
			Y[i*cols + j].imag=0;
			if(j<cols-1){
				Hconj[i*(cols-1)+j].real=0;
				Hconj[i*(cols-1)+j].imag=0;
			}
		}
	}
	
	//Do pre FFT bc library doesn't work first time
	fftOneRow(Y, cols, 0);
	
	//Read in Y (get rid of prefix)
	buffPtr->readNextSymbol(Y, 0);
	
	clock_t start, finish;
	if(timerEn){
		start = clock();
	}
	
	for(int row=0; row<rows; row++){
		//FFT one row 
		fftOneRow(Y, cols, row);
	}
	if(timerEn){
		finish = clock();
		fft[0] = ((float)(finish - start))/(float)CLOCKS_PER_SEC;
	}
	
	if(timerEn){
		start = clock();
	}
	for(int row=0; row<rows; row++){
		//Drop first element and copy it into Hconj
		memcpy(&Hconj[row*(cols-1)], &Y[row*cols+1], (cols-1)* sizeof (*Y));
		
		//shift the row
		//shiftOneRow(Hconj, cols-1, row);
		
		//Divide FFT(Y) by X
		divideOneRow(Hconj, X, cols-1, row);
	}
	
	//take conjugate of H
	
	for (int i = 0; i<rows; i++){  
		for (int j = 0; j<cols-1; j++){
			Hconj[i*(cols-1) + j].imag = -1*Hconj[i*(cols-1) + j].imag;
		}
	}
	if(timerEn){
		finish = clock();
		decode[0] = ((float)(finish - start))/(float)CLOCKS_PER_SEC;
	}
	
	//Now Hconj holds H
	//Save |H|^2 into X
	findDistSqrd(Hconj,X,rows, cols-1);
	/*
	std::string chanEstSym;
	for (int j = 0; j < cols-1; j++) {
		chanEstSym += boost::lexical_cast<std::string>(X[j].real);
		chanEstSym += boost::lexical_cast<std::string>(X[j].imag);
	}
	*/
	int len = 0;
	while (len == 0 and not stop_signal_called) {
		len = socket.send_to(boost::asio::buffer((char*)&X[0], (cols-1)*sizeof(*X)), server_endpoint);
	}
	
}

int main(){
	int rows = numOfRows; // number of vectors -> 16
	int cols = dimension;//dimension -> 1024
	string shm_uid = shmemID;
	
	boost::asio::io_service io_service;
	udp::endpoint local_endpoint(boost::asio::ip::address::from_string(addr), port);
	udp::endpoint server_endpoint(boost::asio::ip::address::from_string(remote_addr), port);
	udp::socket socket(io_service, local_endpoint);
	if (!socket.is_open()) {socket.open(udp::v4());}
	//printf("CPU LS: \n");
	//printInfo();
	
	//Y = 16x1024
	complexF* Y = 0;
	Y = (complexF*)malloc(rows*cols*sizeof(*Y));
	//H (and Hconj) = 16x1023
	complexF* Hconj = 0;
	Hconj = (complexF *)malloc(rows*(cols-1)* sizeof (*Hconj));
	//X = 1x1023 -> later can become |H|^2
	complexF* X = 0;
	X = (complexF *)malloc((cols-1)* sizeof(*X));
	
	// Create shared memory space, return pointer and set as master.
	buffPtr=new ShMemSymBuff(shm_uid, mode);
	std::signal(SIGINT, &sig_int_handler);
	
	//Find H* (H conjugate) ->16x1023 and |H|^2 -> 1x1023
	//while (not stop_signal_called) {
		if(testEn){
			printf("Symbol #0:\n");
		}
		firstVector(Y, Hconj, X, rows, cols, socket, server_endpoint);
	
		for(int i=1; i<numberOfSymbolsToTest; i++){
		
		if(testEn){
			printf("Symbol #%d:\n", i);
		}
		
		
			doOneSymbol(Y, Hconj, X, rows, cols, i, socket, server_endpoint);
			buffIter = i;
		}
	//}
	
	free(Y);
	free(Hconj);
	free(X);
	//delete buffPtr;
	if(timerEn) {
//		printTimes(true);
		storeTimes(true);
}
	
	return 0;

}
