# inmarsatc
C++ library with functions to receive Inmarsat-C signals

Mostly based on Scytale-C source code(https://bitbucket.org/scytalec/scytalec/), rewritten in C++ with some optimizations

Original code license: GNU General Public License 3, microp11 2017

From my information, Scytale-C is based on the Tekmanoid's java app source code, so thanks him for it!


Related projects:

    stdcdec: set of programs to receive inmarsat-c signals
    https://github.com/cropinghigh/stdcdec

    sdrpp-inmarsatc-demodulator: SDR++ module, which replaces stdcdec_demod and stdcdec_decoder
    https://github.com/cropinghigh/sdrpp-inmarsatc-demodulator

    qstdcdec: qt version of stdcdec_parser
    https://github.com/cropinghigh/qstdcdec

Building:

  0.  If you have an arch-like system, just install libinmarsatc-git

  1.  Build

          mkdir build
          cd build
          cmake ..
          make
          sudo make install


Header files:

        /usr/include/inmarsatc_demodulator.h
        /usr/include/inmarsatc_decoder.h
        /usr/include/inmarsatc_parser.h

Library files:

        /usr/lib/libinmarsatc_demodulator.so
        /usr/lib/libinmarsatc_decoder.so
        /usr/lib/libinmarsatc_parser.so

Classes:

  1.  inmarsatc::demodulator::Demodulator

    Class used to demodulate BPSK symbols from the audio samples. Currently the worst working part of the library.

        struct demodulator_result {
            double meanMagnitude; //Mean magnitude of the input samples
            uint8_t bitsDemodulated[DEMODULATOR_SYMBOLSPERCHUNK]; //size=5000, hard symbols, 1byte=1bit
        };
        Demodulator();  // Basic constructor
        bool isCmaEnabled(); //Check, if CMA equalizer is enabled
        bool isAgcEnabled(); //Check, if AGC is enabled
        int getLowFreq(); //Get the low audio frequency limit
        int getHighFreq(); //Get the high audio frequency limit
        double getCenterFreq(); //Get the current audio frequency
        bool getIsInSync(); //Check, if demodulator is in sync with the signal(may be wrong)
        int getNoSyncCount(); //Get count of demodulate() calls resulted in no sync
        std::complex<double> getScatterPoint(); //Get current constellation point
        void setCmaEnabled(bool cmaEnabled); //Enable/disable CMA equalizer
        void setAgcEnabled(bool agcEnabled); //Enable/disable AGC
        void setLowFreq(int lowFreq); //Set the low audio frequency limit
        void setHighFreq(int highFreq); //Set the high audio frequency limit
        void setCenterFreq(double centerFreq); //Force tune demodulator to specified frequency
        void cmaReset(); //Reset the CMA equalizer
        std::vector<demodulator_result> demodulate(std::complex<double> samples[], int length); //Main demodulation function. Accepts array of complex samples. Sample rate=48kS/s. Return vector of demodulated symbols.

   Usage example: stdc_demod

  2.  inmarsatc::decoder::Decoder

    Class used to decode frames from BPSK symbols.

        struct decoder_result {
            uint8_t decodedFrame[DESCRAMBLER_FRAME_LENGTH]; //Actual decoded frame, soft packed, 1byte=8bits, size=640
            int length; //Length of decoded frame, always equals to 640
            int frameNumber; //Number of decoded frame
            std::chrono::time_point<std::chrono::high_resolution_clock> timestamp; //Time point, when the frame was decoded
            bool isHardDecision; //always true
            bool isReversedPolarity; //Are input symbols in reverse polarity
            bool isMidStreamReversePolarity; //Is the polarity reversed in the middle of the frame
            bool isUncertain; //If true, frame may be corrupted or non-complete
            int BER; //Bit Error Rate
        };
        Decoder(int tolerance); // Basic constructor, tolerance=max acceptable BER
        std::vector<decoder_result> decode(uint8_t inputBits[DEMODULATOR_SYMBOLSPERCHUNK]); //Main decoding function. Accepts array of symbols after demodulator, size=640. Returns vector of decoded frames(if any)

   Usage example: stdc_decoder & sdrpp-inmarsatc-demodulator

  3.  inmarsatc::frameParser::FrameParser

    Class used to extract packets from decoded frames

        struct packetDecoder_payload {
            int presentation = -1; //Format of payload(PACKETDECODER_PRESENTATION_IA5, PACKETDECODER_PRESENTATION_ITA2 or PACKETDECODER_PRESENTATION_BINARY)
            std::vector<uint8_t> data8Bit; //Actual data, 1byte=8bits
        };
        struct packetDecoder_result {
            bool isDecodedPacket; //always true
            int frameNumber; //Number of input frame
            std::chrono::time_point<std::chrono::high_resolution_clock> timestamp; //Time point, when the frame was decoded
            uint8_t packetDescriptor; //Byte, represents packet type. Check example or source for the available values
            int packetLength; //Length of packet in bytes
            int decodingStage; //Is this packet fully decoded? PACKETDECODER_DECODING_STAGE_NONE=not, PACKETDECODER_DECODING_STAGE_PARTIAL=not fully, PACKETDECODER_DECODING_STAGE_COMPLETE=yes
            bool isCrc; //If crc of the packet is valid, always true
            packetDecoder_payload payload; //Packet payload, like message
            std::map<std::string, std::string> packetVars; //Variables of decoded packet, check example or source code fo the available fields for any descriptor
            packetDecoder_multiFramePacket mfp; //Variable for internal usage
        };
        struct frameParser_result {
            PacketDecoder::packetDecoder_result decoding_result;
        };
        FrameParser(); // Basic constructor
        std::vector<frameParser_result> parseFrame(decoder::Decoder::decoder_result inputFrame); //Main parsing function. Accepts decoded frame. Returns vector of parsed packets(if any)

   Usage example: stdc_parser & qstdcdec

   WARNING! All messages are directed to their recipients! If you're not the recipient, you should delete received message!
