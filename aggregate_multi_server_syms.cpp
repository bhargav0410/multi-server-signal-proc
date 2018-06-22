#include "ShMemSymBuff.hpp"
#include <csignal>
#include <fstream>
#include <boost/asio.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/thread.hpp>
#ifndef numOfServers
	#define numOfServers 4
#endif

#ifndef addr
	#define addr "10.10.0.10"
#endif

using namespace std;
using boost::asio::ip::udp;

boost::asio::io_service io_service;
//boost::array<complexF, numOfServers*(dimension-1)> rec;
complexF rec[numOfServers][(numOfRows/numOfServers)*dimension-1];


void sync_receiver(std::string addr, std::string remote_addr, int port, int server) {
	udp::endpoint local_endpoint(boost::asio::ip::address::from_string(addr), port);
	udp::endpoint server_endpoint(boost::asio::ip::address::from_string(remote_addr), port);
	udp::socket socket(io_service, local_endpoint);
	if (!socket.is_open()) {socket.open(udp::v4());}
	int len = 0;
	while (len <= 0) {
		len = socket.receive_from(boost::asio::buffer(&rec[server]), server_endpoint);
	}
	udp_check.clear();
	for (int i = 0; i < len; i++) {
		udp_check.push_back(rec[i]);
	}
	//std::cout << udp_check << " " << len;
	//len = 0;
}

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

int main() {
	int rows = numOfRows;
	int cols = dimension;
	
	return 0;
}