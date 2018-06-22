#include <fftw3.h>
//Shared Memory 
#include "CSharedMemSimple.hpp"
#include "ShMemSymBuff.hpp"
#include <csignal>
#include <fstream>
#include <cstring>
#include <boost/lexical_cast.hpp>
#include <complex>
#include <vector>
#include <cstdlib>

#define mode 1
#define FFT_size dimension
#define cp_size prefix
#define numSymbols lenOfBuffer
#define chan numOfRows
ShMemSymBuff* buffPtr;

static bool stop_signal_called = false;
void sig_int_handler(int){stop_signal_called = true;}

std::string file = "corr_rec_ch_";
std::ifstream infile;
std::vector<std::vector<std::complex<float> > > copy_buff;

int main() {
	std::cout << "File to shared memory:\n";
	printInfo();
	std::signal(SIGINT, &sig_int_handler);
	int iter = 1;
	std::complex<float>* copy_to_mem = 0;
	std::string shm_uid = shmemID;
	buffPtr=new ShMemSymBuff(shm_uid, mode);
	
	std::string File = file + "0_binary";
	infile.open(File.c_str(), std::ifstream::binary);
	infile.seekg(0, infile.end);
	size_t num_tx_samps = infile.tellg()/sizeof(std::complex<float>);
	infile.seekg(0, infile.beg);
	infile.close();
	
	std::cout << "Number of samples = " << num_tx_samps << std::endl;
	
	copy_buff.resize(
        chan, std::vector<std::complex<float> >(num_tx_samps)
    );
	
	std::cout << "Copying from file...\n";
	for (int i = 0; i < chan; i++) {
		//std::cout << "File " << i << std::endl;
		std::string File = file + boost::lexical_cast<std::string>(i) + "_binary";
		infile.open(File.c_str(), std::ifstream::binary);
		infile.read((char*)&copy_buff[i].front(), num_tx_samps*sizeof(std::complex<float>));
		std::cout << "From file \"" << File << "\" - ";
		std::cout << copy_buff[i][0] << ",";
		infile.close();
	}
	
	copy_to_mem = (std::complex<float>*)malloc((chan*(FFT_size+cp_size)*sizeof(*copy_to_mem)));
	std::cout << "Copying to shared memory...\n";
	//while(not stop_signal_called) {
		for (int i = 0; i < numSymbols; i++) {
			//std::cout << "Symbol " << i << std::endl;
			for (int j = 0; j < chan; j++) {
				memcpy(&copy_to_mem[j*(FFT_size+cp_size)], &copy_buff[j][i*(FFT_size+cp_size)], (FFT_size+cp_size)*sizeof(*copy_to_mem));
				/*
				for (int k = 0; k < FFT_size+cp_size; k++) {
					copy_to_mem[j*(FFT_size+cp_size) + k].real = copy_buff[j][i*(FFT_size+cp_size) + k].real();
					copy_to_mem[j*(FFT_size+cp_size) + k].imag = copy_buff[j][i*(FFT_size+cp_size) + k].imag();
				}
				*/
			}
			buffPtr->writeNextSymbolWithWait(copy_to_mem);
		}
	//	printOutArr(copy_to_mem,1,cp_size+FFT_size);
		iter++;
	//}
	
	
	delete buffPtr;
	return 0;
}