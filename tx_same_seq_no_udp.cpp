#include <uhd/utils/thread_priority.hpp>
#include <uhd/utils/safe_main.hpp>
#include <uhd/utils/static.hpp>
#include <uhd/usrp/multi_usrp.hpp>
#include <uhd/exception.hpp>
#include <uhd/types/time_spec.hpp>
#include <uhd/transport/udp_simple.hpp>
#include <boost/program_options.hpp>
#include <boost/math/special_functions/round.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/thread.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/asio.hpp>
#include <boost/array.hpp>
#include <iostream>
#include <csignal>
#include <ctime>
#include <fstream>
#include <cmath>
#include <sys/socket.h>
#include <string>

#define pi 3.141592654
#define c 299792458.0

namespace po = boost::program_options;
using boost::asio::ip::udp;
/***********************************************************************
 * Signal handlers
 **********************************************************************/
static bool stop_signal_called = false;
void sig_int_handler(int){stop_signal_called = true;}
size_t len = 0;
void handler(
  const boost::system::error_code& error, // Result of operation.
  size_t length          // Number of bytes received.
  ) {
	  len = length;
  } 

void udp_receiver(udp::socket &socket, udp::endpoint &server_endpoint) {
	boost::array<char, uhd::transport::udp_simple::mtu> rec;
	while (len <= 0) {
		len = socket.receive_from(boost::asio::buffer(rec), server_endpoint);
	}
	/*
	std::cout << "Rec: ";
	std::cout.write(rec.data(),len);
	std::cout << std::endl;
	*/
}
  
/***********************************************************************
 * Main function
 **********************************************************************/
int UHD_SAFE_MAIN(int argc, char *argv[]){
    uhd::set_thread_priority_safe();
	
    //variables to be set by po
    std::string args, ant, subdev, ref, pps, otw, cpu, channel_list, file, addr, remote_addr;
    double rate, freq, gain, bw;
    float delay;
	int port, start ,end, num_times;
	size_t nsamps; 
    //setup the program options
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "help message")
        ("args", po::value<std::string>(&args)->default_value(""), "single uhd device address args")
        ("rate", po::value<double>(&rate), "rate of outgoing samples")
        ("freq", po::value<double>(&freq), "RF center frequency in Hz")
        ("gain", po::value<double>(&gain), "gain for the RF chain")
		("nsamps", po::value<size_t>(&nsamps)->default_value(4096), "Number of samples of the packet")
        ("ant", po::value<std::string>(&ant), "antenna selection")
        ("subdev", po::value<std::string>(&subdev), "subdevice specification")
        ("bw", po::value<double>(&bw), "analog frontend filter bandwidth in Hz")
        ("ref", po::value<std::string>(&ref)->default_value("internal"), "clock reference (internal, external, mimo, gpsdo)")
        ("pps", po::value<std::string>(&pps), "PPS source (internal, external, mimo, gpsdo)")
        ("otw", po::value<std::string>(&otw)->default_value("sc16"), "specify the over-the-wire sample mode (sc16 or sc8)")
		("cpu", po::value<std::string>(&otw)->default_value("fc32"), "specify the cpu sample mode (fc32 or sc16)")
        ("channels", po::value<std::string>(&channel_list)->default_value("0"), "which channels to use (specify \"0\", \"1\", \"0,1\", etc)")
        ("int-n", "tune USRP with integer-N tuning")
		("delay", po::value<float>(&delay)->default_value(0), "Delay value in stream time spec")
		("file", po::value<std::string>(&file), "name of the file to take samples from")
		("verbose","shows the sync messages sent and received via UDP")		
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    //print the help message
    if (vm.count("help")){
        std::cout << boost::format("UHD TX Waveforms %s") % desc << std::endl;
        return ~0;
	}
    //create a usrp device
    std::cout << std::endl;
    std::cout << boost::format("Creating the usrp device with: %s...") % args << std::endl;
    uhd::usrp::multi_usrp::sptr usrp = uhd::usrp::multi_usrp::make(args);

    //detect which channels to use
    std::vector<std::string> channel_strings;
    std::vector<size_t> channel_nums;
    boost::split(channel_strings, channel_list, boost::is_any_of("\"',"));
    for(size_t ch = 0; ch < channel_strings.size(); ch++){
        size_t chan = boost::lexical_cast<int>(channel_strings[ch]);
        if(chan >= usrp->get_tx_num_channels())
            throw std::runtime_error("Invalid channel(s) specified.");
        else
            channel_nums.push_back(boost::lexical_cast<int>(channel_strings[ch]));
    }


    //Lock mboard clocks
    usrp->set_clock_source(ref);

    //always select the subdevice first, the channel mapping affects the other settings
    if (vm.count("subdev")) usrp->set_tx_subdev_spec(subdev);

    std::cout << boost::format("Using Device: %s") % usrp->get_pp_string() << std::endl;

    //set the sample rate
    if (not vm.count("rate")){
        std::cerr << "Please specify the sample rate with --rate" << std::endl;
        return ~0;
    }
    std::cout << boost::format("Setting TX Rate: %f Msps...") % (rate/1e6) << std::endl;
    usrp->set_tx_rate(rate);
    std::cout << boost::format("Actual TX Rate: %f Msps...") % (usrp->get_tx_rate()/1e6) << std::endl << std::endl;

    //set the center frequency
    if (not vm.count("freq")){
        std::cerr << "Please specify the center frequency with --freq" << std::endl;
        return ~0;
    }

    for(size_t ch = 0; ch < channel_nums.size(); ch++) {
        std::cout << boost::format("Setting TX Freq: %f MHz...") % (freq/1e6) << std::endl;
        uhd::tune_request_t tune_request(freq);
        if(vm.count("int-n")) tune_request.args = uhd::device_addr_t("mode_n=integer");
        usrp->set_tx_freq(tune_request, channel_nums[ch]);
        std::cout << boost::format("Actual TX Freq: %f MHz...") % (usrp->get_tx_freq(channel_nums[ch])/1e6) << std::endl << std::endl;

        //set the rf gain
        if (vm.count("gain")){
            std::cout << boost::format("Setting TX Gain: %f dB...") % gain << std::endl;
            usrp->set_tx_gain(gain, channel_nums[ch]);
            std::cout << boost::format("Actual TX Gain: %f dB...") % usrp->get_tx_gain(channel_nums[ch]) << std::endl << std::endl;
        }

        //set the analog frontend filter bandwidth
        if (vm.count("bw")){
            std::cout << boost::format("Setting TX Bandwidth: %f MHz...") % bw << std::endl;
            usrp->set_tx_bandwidth(bw, channel_nums[ch]);
            std::cout << boost::format("Actual TX Bandwidth: %f MHz...") % usrp->get_tx_bandwidth(channel_nums[ch]) << std::endl << std::endl;
        }

        //set the antenna
        if (vm.count("ant")) usrp->set_tx_antenna(ant, channel_nums[ch]);
    }

    boost::this_thread::sleep(boost::posix_time::seconds(1)); //allow for some setup time

	//create a transmit streamer
	//linearly map channels (index0 = channel0, index1 = channel1, ...)
	uhd::stream_args_t stream_args("fc32");
	stream_args.channels = channel_nums;
	uhd::tx_streamer::sptr tx_stream = usrp->get_tx_stream(stream_args);

	std::cout << boost::format("Setting device timestamp to 0...") << std::endl;
	if (channel_nums.size() > 1)
	{
		// Sync times
		if (pps == "mimo")
		{
			UHD_ASSERT_THROW(usrp->get_num_mboards() == 2);

			//make mboard 1 a slave over the MIMO Cable
			usrp->set_time_source("mimo", 1);

			//set time on the master (mboard 0)
			usrp->set_time_now(uhd::time_spec_t(0.0), 0);

			//sleep a bit while the slave locks its time to the master
			boost::this_thread::sleep(boost::posix_time::milliseconds(100));
		}
		else
		{
			if (pps == "internal" or pps == "external" or pps == "gpsdo")
				usrp->set_time_source(pps);
			usrp->set_time_unknown_pps(uhd::time_spec_t(0.0));
			boost::this_thread::sleep(boost::posix_time::seconds(1)); //wait for pps sync pulse
		}
	}
	else
	{
		usrp->set_time_now(0.0);
	}

	//Check Ref and LO Lock detect
	std::vector<std::string> sensor_names;
	const size_t tx_sensor_chan = channel_list.empty() ? 0 : boost::lexical_cast<size_t>(channel_list[0]);
	sensor_names = usrp->get_tx_sensor_names(tx_sensor_chan);
	if (std::find(sensor_names.begin(), sensor_names.end(), "lo_locked") != sensor_names.end()) {
		uhd::sensor_value_t lo_locked = usrp->get_tx_sensor("lo_locked", tx_sensor_chan);
		std::cout << boost::format("Checking TX: %s ...") % lo_locked.to_pp_string() << std::endl;
		UHD_ASSERT_THROW(lo_locked.to_bool());
	}
	const size_t mboard_sensor_idx = 0;
	sensor_names = usrp->get_mboard_sensor_names(mboard_sensor_idx);
	if ((ref == "mimo") and (std::find(sensor_names.begin(), sensor_names.end(), "mimo_locked") != sensor_names.end())) {
		uhd::sensor_value_t mimo_locked = usrp->get_mboard_sensor("mimo_locked", mboard_sensor_idx);
		std::cout << boost::format("Checking TX: %s ...") % mimo_locked.to_pp_string() << std::endl;
		UHD_ASSERT_THROW(mimo_locked.to_bool());
	}
	if ((ref == "external") and (std::find(sensor_names.begin(), sensor_names.end(), "ref_locked") != sensor_names.end())) {
		uhd::sensor_value_t ref_locked = usrp->get_mboard_sensor("ref_locked", mboard_sensor_idx);
		std::cout << boost::format("Checking TX: %s ...") % ref_locked.to_pp_string() << std::endl;
		UHD_ASSERT_THROW(ref_locked.to_bool());
	}



	//Copying data from file to buffers
	std::ifstream infile(file.c_str(), std::ifstream::binary);
	infile.seekg(0, infile.end);
	size_t num_tx_samps = infile.tellg()/sizeof(std::complex<float>);
	infile.seekg(0, infile.beg);
	if (!infile.is_open()) return ~0;
	std::vector<std::vector<std::complex<float> > > temp_buff(channel_nums.size(), std::vector<std::complex<float> > (num_tx_samps*channel_nums.size()));
	for (int i = 0; i < channel_nums.size(); i++) {
		if (i != 0) {
			for (int n = 0; n < num_tx_samps*(i-1); n++) {
				temp_buff.at(i).at(n) = 0;
			}
		}
		infile.read((char*)&temp_buff.at(i).at(num_tx_samps*i), num_tx_samps*sizeof(std::complex<float>));
		//if (i != channel_nums.size()) {
			/*
		for (int n = num_tx_samps*(i+1); n < temp_buff.at(i).size(); n++) {
			temp_buff.at(i).at(n) = 0;
		}
		//}
		infile.seekg(0, infile.beg);
		std::cout << "Pos: " << infile.tellg() << std::endl;
		*/
	}
	//allocate a buffer which we re-use for each channel
	//if (!infile.is_open()) {std::cout << "File not read...\n";}
	std::cout << std::endl << "Tx samps: " << num_tx_samps << std::endl;
	std::vector<std::vector<std::complex<float> > > buff(channel_nums.size(), std::vector<std::complex<float> > (num_tx_samps*channel_nums.size()));
	std::vector<std::complex<float> *> buffs;
	for (int ch = 0; ch < channel_nums.size(); ch++) buffs.push_back(&buff[ch].front());

	std::signal(SIGINT, &sig_int_handler);
	std::cout << "Press Ctrl + C to stop streaming..." << std::endl;

	
	//UDP connection for sync
	/* std::string flag = "TX";
	boost::asio::io_service io_service;
	udp::endpoint local_endpoint(boost::asio::ip::address::from_string(addr), port);
	udp::endpoint server_endpoint(boost::asio::ip::address::from_string(remote_addr), port);
	udp::socket socket(io_service, local_endpoint);
	if (!socket.is_open()) {socket.open(udp::v4());}
	while (len == 0 and not stop_signal_called) {
		len = socket.send_to(boost::asio::buffer(flag.c_str(), flag.size()), server_endpoint);
	}
	flag.clear();
	if (vm.count("verbose")) {
	std::cout << "Starting transmission...\n";
	}
	boost::array<char, uhd::transport::udp_simple::mtu> rec;
	//boost::array<char, 5> temp = ;
	len = 0;
	//rec.resize(uhd::transport::udp_simple::mtu);
	while (len <= 0 and not stop_signal_called) {
		len = socket.receive_from(boost::asio::buffer(rec), server_endpoint);
	}
	if (vm.count("verbose")) {
	std::cout << "Rec: ";
	std::cout.write(rec.data(),len);
	std::cout << std::endl;
	}
	len = 0; */
	
	
	
	
	
	uhd::tx_metadata_t md;

	//for (int j = 0; j < num_times; j++) {
		//for (int i  = 0; i <= channel_nums.size(); i++){
	if (pps == "internal" or pps == "external" or pps == "gpsdo") {	
		usrp->set_time_unknown_pps(uhd::time_spec_t(0.0));
		boost::this_thread::sleep(boost::posix_time::seconds(2)); //wait for pps sync pulse
	}
	//if (stop_signal_called) break;
	
	
	
	/*
	if (i != channel_nums.size()) {
		for (size_t ch = 0; ch < channel_nums.size(); ch++) {
			for (size_t n = 0; n < buff.at(ch).size(); n++){
				if (ch == i) { 
					buff.at(ch)[n] = temp_buff.at(ch)[n];
				}
				else {
					buff.at(ch)[n] = 0;
				}
			}
		}
	}
	else {
		for (size_t ch = 0; ch < channel_nums.size(); ch++) {
			for (size_t n = 0; n < buff.at(ch).size(); n++){
				buff.at(ch)[n] = temp_buff.at(ch)[n];
			}
		}
	}
	*/
	for (size_t ch = 0; ch < channel_nums.size(); ch++) {
		for (size_t n = 0; n < buff.at(ch).size(); n++){
			buff.at(ch)[n] = temp_buff.at(ch)[n];
		}
	}

	
	/* while (len <= 0 and not stop_signal_called) {
		len = socket.receive_from(boost::asio::buffer(rec), server_endpoint);
	}
	if (vm.count("verbose")) {
	std::cout << "Rec: ";
	std::cout.write(rec.data(),len);
	std::cout << std::endl;
	}
	len = 0; */

	
	
	md.start_of_burst = false;
	md.end_of_burst   = false;
	md.has_time_spec = true;
	std::cout << "Sending samples...\n";
	/* flag += "Start RX";
	if (vm.count("verbose")) {
	std::cout << flag << std::endl;
	}
	while (len == 0 and not stop_signal_called) {
		len = socket.send_to(boost::asio::buffer(flag.c_str(), flag.size()), server_endpoint);
	}
	flag.clear();
	len = 0; */
	
	//socket.non_blocking(true);
	//boost::thread t(udp_receiver, boost::ref(socket), boost::ref(server_endpoint));
	md.time_spec = usrp->get_time_now() + uhd::time_spec_t(delay);
	while(not stop_signal_called){
		//send the entire contents of the buffer
		tx_stream->send(buffs, buff.at(0).size(), md);
		md.has_time_spec = false;
	}
	md.end_of_burst = true;
	tx_stream->send("", 0, md);
	//len = 0;
//	boost::this_thread::sleep(boost::posix_time::milliseconds(1000));
	/* flag += "Stopped";
	if (vm.count("verbose")) {
	std::cout << flag << std::endl;
	}
	while (len == 0 and not stop_signal_called) {
		len = socket.send_to(boost::asio::buffer(flag.c_str(), flag.size()), server_endpoint);
	}
	flag.clear();
	len = 0; */
	
	std::cout << "Samples sent...\n";
			
	//	}
//	}


	//send a mini EOB packet
	//md.end_of_burst = true;
    //finished
    std::cout << std::endl << "Done!" << std::endl << std::endl;
    return EXIT_SUCCESS;
}


