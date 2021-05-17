//DSP AND PACKET DECODER CODE TAKEN FROM SCYTALE-C, WHICH IS THE OPEN-SOURCE INMARSAT-C DECODER WRITTEN ON C#. I JUST PORTED IT TO C++(with some optimizations).
//Code license: GNU GPL

#ifndef INMARSATC_H
#define INMARSATC_H

#if defined(_MSC_VER) || defined(WIN64) || defined(_WIN64) || defined(__WIN64__) || defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
#  define Q_DECL_EXPORT __declspec(dllexport)
#  define Q_DECL_IMPORT __declspec(dllimport)
#else
#  define Q_DECL_EXPORT     __attribute__((visibility("default")))
#  define Q_DECL_IMPORT     __attribute__((visibility("default")))
#endif

#if defined(INMARSATC_LIBRARY)
#  define INMARSATC_EXPORT Q_DECL_EXPORT
#else
#  define INMARSATC_EXPORT Q_DECL_IMPORT
#endif

#include <complex>
#include <vector>
#include <algorithm>
#include <math.h>
#include <mutex>

namespace inmarsatc {
    namespace demodulator {
        class FIR {
            public:
                FIR(double b[], int len);
                double filter(double sample);
                double filter();
                void resetFilter();
            private:
                int psLength;
                double* b;
                double* prevSamples;
        };

        class RRC {
            public:
                RRC(double alpha, int firstSize, double sampleRate, double symbolRate);
                double filter(double sample);
            private:
                double* points;
                FIR* rrcFIR;
        };

        class Gardner {
            public:
                Gardner(double sampleRate, double symbolRate);
                bool step(double re, double im, std::complex<double>* output, bool* onPoint);
            private:
                double ts;
                double symbol2xTimer;
                int samplesAgo;
                std::complex<double> maxErrorxAggression;
                int isOnPoint;
                std::complex<double> currentSample;
                std::complex<double> sampleThisOn;
                std::complex<double> error;
                std::complex<double> sampleLastOn;
                std::complex<double> sampleOff;
                double aggression;
                bool result;
                std::complex<double> dummyOutput;
        };

        class CMA {
            public:
                #define CMA_SZ 9
                CMA();
                bool step(std::complex<double> sample, bool isOnPoint, std::complex<double>* output);
                void CMAReset();
                bool result;
            private:
                std::complex<double> dummyOutput;
                std::complex<double> cmaW[CMA_SZ];
                std::complex<double> cmaX[CMA_SZ];
                double beta;
                double stepSize;
                double mean;
                std::complex<double> cmaEqualizerOut;
                std::complex<double> error;
        };

        class AGC {
            public:
                AGC();
                void apply(double* i, double* q);
                void apply(std::complex<double>* value);
                void reset();
            private:
                double mAGC;
        };

        class INMARSATC_EXPORT Demodulator {
            public:
                #define DEMODULATOR_SYMBOLSPERCHUNK 5000
                struct demodulator_result {
                    double meanMagnitude;
                    uint8_t bitsDemodulated[DEMODULATOR_SYMBOLSPERCHUNK];
                };

                Demodulator();
                bool isCmaEnabled();
                bool isAgcEnabled();
                int getLowFreq();
                int getHighFreq();
                double getCenterFreq();
                bool getIsInSync();
                int getNoSyncCount();
                std::complex<double> getScatterPoint();
                void setCmaEnabled(bool cmaEnabled);
                void setAgcEnabled(bool agcEnabled);
                void setLowFreq(int lowFreq);
                void setHighFreq(int highFreq);
                void setCenterFreq(double centerFreq);
                void cmaReset();
                std::vector<demodulator_result> demodulate(std::complex<double> samples[], int length);

            private:
                #define DEMODULATOR_SYMBOLRATE 1200.0
                #define DEMODULATOR_SAMPLERATE 48000.0
                // As the carrier gets locked, the positive and the negative error are equally distributed and
                // their sum gets closer and closer to zero. This is the absolute number I came up with.
                #define DEMODULATOR_MAXABSERRSUM 12
                #define DEMODULATOR_ALPHA 0.0065
                // this implementation of a Costas loop tracks easier to the left, so we set the
                // center frequency on purpuse a bit higher that where we expect to find the sync
                #define DEMODULATOR_DEFAULT_FREQ 2600.0
                #define DEMODULATOR_LOW_FREQ 500.0
                #define DEMODULATOR_HIGH_FREQ 4500.0
                double freq = DEMODULATOR_DEFAULT_FREQ;
                std::mutex freq_mtx;
                const double alpha = DEMODULATOR_ALPHA;
                // see literature
                const double beta = alpha * alpha / 4;
                double omega;
                double phase = 0.0;
                int noSyncCounter = 0;
                double error;
                bool isInSync = false;
                int flagcounter;
                std::complex<double> scatterPoint;
                std::complex<double> lastScatterPoint;
                // LPF1, LPF2
                // The literature suggestes that the square-shaped LPF filter should have a
                // bandwidth of no less than half the symbol rate. This will remove the most
                // noise possible without reducing the amplitude of the desired signal at
                // its sampling instants.
                double blpf12[2] = {0.9, 0.9};
                FIR* lpf1;
                FIR* lpf2;
                // The literature indicates that the loop filter should have a
                // response far outside the LPF1/LPF2.The fastest achievable settle
                // time is one in which the VCO has a gain 8x that of the LPF1/LPF2 pole
                // frequency.For LPF3, a factor of six times K or 12 times the LPF1/LPF2
                // pole is a better choice
                double blpf3[2] = {0.48, 0.48};
                FIR* lpf3;
                RRC* rrcRe;
                RRC* rrcIm;
                RRC* rrcRe2;
                RRC* rrcIm2;
                Gardner* gardner;
                CMA* cma;
                bool cmaEnabled = false;
                std::complex<double> gardnerOutputSample;
                std::complex<double> cmaOutputSample;
                uint8_t symbolBuffer[DEMODULATOR_SYMBOLSPERCHUNK];
                AGC* agc;
                bool agcEnabled = false;
                int loFreq = DEMODULATOR_LOW_FREQ;
                int hiFreq = DEMODULATOR_HIGH_FREQ;
        };

        class complexMath {
            public:
                static bool lessThan(std::complex<double> left, std::complex<double> right);
                static bool biggerThan(std::complex<double> left, std::complex<double> right);
        };
    }

    namespace decoder {

        class UWFinder {
            /// The UW finder works with a matrix of interleaved symbols as described above.
            /// The BPSK demodulator has phase ambiguity, and also there is the possibility
            /// that mid-stream the polarity of the symbols will change if the demodulator
            /// steps over a symbol.
            ///
            /// This also represents the frame synchronization.
            public:
                #define UWFINDER_FRAME_LENGTH 10368
                #define UWFINDER_MAXTOLERANCE 30
                /// From this moment on we are working on "frames" of data, the stream has been cut into parts of given length.
                /// The symbols in the uwfinder_result are still in the order they came out of the demodulator.
                struct uwfinder_result {
                    int symbolCount;
                    int length;
                    bool isHardDecision;
                    bool isReversedPolarity;
                    bool isMidStreamReversePolarity;
                    bool isUncertain;
                    int BER;
                    uint8_t uwFrame[UWFINDER_FRAME_LENGTH]; //array of bits
                };
                /// We set the highest allowed number of incorrect UW words in a frame.
                /// If a frame has more incorrect UW than the tolerance, we skip the frame.
                void SetTolerance(int tolerance);
                /// One full data frame has 10240 symbols.
                /// Add to this 2 x 64 UW symbols and we get 10368 symbols per received frame.
                /// We use a buffer of 2 full frames to be able to store the symbols that come before a frame is detected.
                /// Beyound 2 frames, what we miss is discarded.
                ///
                /// We store the symbols into the symbol buffer, one symbol at a time.
                /// The first symbol received is always at the highest index in the symbols array parameter.
                /// As new symbols come, the buffer content is rotated to the right to make space for the newcomer.
                /// As the newcomer arrives, we run the unique word finding sequence.
                /// When the current frame is detected, we check the number of symbols pushed beyound the frame size and if a good
                /// amount, say 80% of a frame is in there, we attempt to decode that as well.
                /// If we manage to recover that frame, we send it down the processing chain before the current frame.
                std::vector<uwfinder_result> Decode(uint8_t bitsDemodulated[DEMODULATOR_SYMBOLSPERCHUNK], bool isHardDecision);
            private:

                /// This is the UW pattern at the receiver, 64 symbols for Inmarsat-C, see table in specs.
                std::vector<uint8_t> nrmPolUwPattern = {
                        0, 0, 0, 0,    0, 1, 1, 1,   1, 1, 1, 0,   1, 0, 1, 0,
                        1, 1, 0, 0,    1, 1, 0, 1,   1, 1, 0, 1,   1, 0, 1, 0,
                        0, 1, 0, 0,    1, 1, 1, 0,   0, 0, 1, 0,   1, 1, 1, 1,
                        0, 0, 1, 0,    1, 0, 0, 0,   1, 1, 0, 0,   0, 0, 1, 0
                };
                /// We must use a reverse polarity search pattern as well due to phase ambiguity
                std::vector<uint8_t> revPolUwPattern = {
                    1, 1, 1, 1,    1, 0, 0, 0,   0, 0, 0, 1,   0, 1, 0, 1,
                    0, 0, 1, 1,    0, 0, 1, 0,   0, 0, 1, 0,   0, 1, 0, 1,
                    1, 0, 1, 1,    0, 0, 0, 1,   1, 1, 0, 1,   0, 0, 0, 0,
                    1, 1, 0, 1,    0, 1, 1, 1,   0, 0, 1, 1,   1, 1, 0, 1
                };
                uint8_t symbolRegister[UWFINDER_FRAME_LENGTH*2];
                uint8_t symbolRegisterUncertain[UWFINDER_FRAME_LENGTH];
                // There is no reason why this should not start at zero.
                // however if the input is a file, by setting the counter higher we force the evaluation of
                // less than full frames as well.
                // Of course, to make this really nice, we should also on close flush our register, just in
                // case there is one almost complete frame not getting pushed through fully.
                int symbolCount = 10239;
                int tolerance = UWFINDER_MAXTOLERANCE;
                /// Frame detector method. It computes the distribution of unique words inside a packet of known length.
                /// If we perform this operation every time a new symbol comes out of the demodulator, there will be a
                /// point where the distribution will match the encoded stream and in that moment we have detected the frame.
                /// The array of symbols we are using has 2 * UWFINDER_FRAMELENGTH, that is 2 * 10368 symbols.
                /// We are parting the symbol array in two parts, | lower | higher | where the last arrived symbol sits at the
                /// leftmost position in the lower part.
                /// We always compute 10368 symbols, either the lower or the higher part only.
                /// lowestFrame: If true, compute the lower part, otherwise the higher part.
                /// nUW: Returns the number of times the normal UW was found.
                /// rUW: returns the number of times the reverse polarity UW was found.
                /// isReversedPolarity: Indicates reversed polarity. To decode this frame we will have to reverse all the symbols later on.
                /// isMidStreamReversePolarity: This indicates a frame that has symbols that were one polarity for a while and then another.
                /// This could happen id the demodulator loses sync, or the audio input is losing samples, or the source is of low quality, or...
                /// isReversedFirst: It indicates the polarity of the first computed symbols.
                bool IsFrameDetected(bool lowestFrame, int* nUW, int* rUW, bool* isReversedPolarity, bool* isMidStreamReversePolarity, bool* isReversedFirst);

        };


        class Depermuter {
            /**
             *   Bibliography:
             *
             *   Calcutt, D. M., & Tetley, L. (2004). Satellite communications: principles and applications. Oxford: Elsevier.
             *   https://www.amazon.com/Satellite-Communications-Applications-David-Calcutt/dp/034061448X
             *
             *   Nera. (2015) Nera Inmarsat-C Service Manual. Billingstad: Nera ASA.
             *   https://www.manualslib.com/manual/1201514/Nera-Inmarsat-C.html
             *
             * The interleave matrix (block) consists of 64 rows by 162 columns.
             * The matrix holds 10368 symbols.
             * The block is transmitted on a row by row basis.
             * The symbols in a row are transmitted in ascending order of column position, UW words first.
             * The rows are not transmitted in sequential order.
             * They are transmitted according to a permutted sequence.
             * If the rows in the interleaved block are numbered from i = 0 to i = 63 sequentially and the
             * transmitted order is from j = 0 sequentially to j = 63, then i and j are related as it follows:
             * i = (j x 39) modulo 64
             * j = (i x 23) modulo 64
             * 1st output j = 0, i = (0 x 39) modulo 64 = 0
             * 2nd output j = 1, i = (1 x 39) modulo 64 = 39
             * 3rd output j = 2, i = (2 x 39) modulo 64 = 14
             *
             * For any of the results, both i an j can only have values from 0 to 63.
             * Any bits representing values of 64, 128 and above can be ignored.
             * The final value is given by the right-hand six bits only.
             *
             * When the frame gets into the Depermuter on the receiver side, the original matrix is an
             * array of symbols, where the first arriving line (j above) occupies the highest index in the array.
             * The highest position in the array is occupied by the two UW symbols of the first line.
             *
             * We first reverse this symbol array aka frame aka block aka matrix.
             * After reversing the array, we obtain the original matrix, permutted, where i=j=0
             * We depermute it by using the i derived from j above:
             * 1st output i = 0, j = (0 x 23) modulo 64 = 0
             * 2nd output i = 1, j = (1 x 23) modulo 64 = 23
             * 3rd output i = 2, j = (2 x 23) modulo 64 = 46
             *
             * We also know that each row has 162 columns. It means we can precompute a depermutting array that
             * would help us depermute in a simple loop the entire frame. We can precompute this array either in
             * the constructor or when processing the first packet. It is more economical than calculating the
             * values on-the-fly.
             *
             */
            public:
                #define DEPERMUTER_FRAME_LENGTH 10368
                struct depermuter_result {
                    int length;
                    uint8_t depermutedFrame[DEPERMUTER_FRAME_LENGTH];
                    bool isHardDecision;
                };
                Depermuter();
                depermuter_result depermute(uint8_t uwFrame[UWFINDER_FRAME_LENGTH], bool isHardDecision);
            private:
                std::vector<int> depermuttingArray;
        };


        class Deinterleaver {
            /*
            * The interleave matrix (block) consists of 64 rows by 162 columns.
            * The matrix holds 10368 symbols.
            *
            * The deinterleaver matrix contains 64 rows and only 160 columns as the UW are removed.
            * The deinterleaver matrix contains 10240 symbols.
            *
            * At the end of processing, the output will contain 10240 data symbols representing the received
            * output of the transmitter's convolutional encoder
            */
            public:
                #define DEINTERLEAVER_FRAME_LENGTH 10240
                struct deinterleaver_result {
                    int length;
                    uint8_t deinterleavedFrame[DEINTERLEAVER_FRAME_LENGTH];
                    bool isHardDecision;
                };
                deinterleaver_result deinterleave(uint8_t depermutedFrame[DEPERMUTER_FRAME_LENGTH], bool isHardDecision);
            private:
                //create an interleaving matrix
                //we do not need space for the UW as we will be descarding them at this stage
                uint8_t deinterleverMatrix[64][160];
        };


        class ViterbiDecoder {
            /*
             * Bibliography:
             *
             * ***************************************************************************************
             * Title: KA9Q Viterbi decoder V1.0 source code
             * Author: Karn, P
             * Date: 1995, March 18
             * Code version: V1.0
             * Availability: ftp://ftp.ucsd.edu/hamradio/dsp/viterbi.txt
             *               ftp://ftp.ucsd.edu/hamradio/dsp/viterbi.zip
             * Copyright 1995 Phil Karn, KA9Q
             *
             * Viterbi decoder for the NASA standard rate 1/2 constraint length 7 convolutional code.
             *
             * ***************************************************************************************
             *
             * Nera. (2015) Nera Saturn C Technical Manual.Billingstad: Nera ASA.
             * ftp://161.87.213.193.static.cust.telenor.com/Manuals/Saturn%20C/SatC_Marine_Tech_Manual_A.pdf
             *
             * Another C# Viterbi decoder:
             * https://github.com/n8ohu/HamModem/blob/master/HamModem/Viterbi.cs
             *
             */
            public:
                #define VITERBIDECODER_FRAME_LENGTH 640
                struct viterbidecoder_result {
                    int length;
                    uint8_t viterbiFrame[VITERBIDECODER_FRAME_LENGTH];
                };
                // Generate metric tables for a soft-decision convolutional decoder
                // assuming gaussian noise on a PSK channel.
                //
                // Works from "first principles" by evaluating the normal probability
                // function and then computing the log-likelihood function
                // for every possible received symbol value
                //
                // Copyright 1995 Phil Karn, KA9Q
                //
                // Symbols are offset-binary, with 128 corresponding to an erased (no
                // information) symbol
                //
                // Normal function integrated from -Inf to x. Range: 0-1
                //
                // Generate log-likelihood metrics for 8-bit soft quantized channel
                // assuming AWGN and BPSK
                ViterbiDecoder();
                // According to the documentation, the encoding scrambler produces 639 bytes of information and one flush byte.
                // These 640 bytes are then sent to the convolutional encoder.
                // As such we do not need be concerned with padding in any way the symbols we input into the viterbi decoder.
                //
                // Viterbi decoder for K=7 rate=1/2 convolutional code
                // Copyright 1995 Phil Karn, KA9Q
                //
                //
                // The basic Viterbi decoder operation, called a "butterfly"
                // operation because of the way it looks on a trellis diagram. Each
                // butterfly involves an Add-Compare-Select (ACS) operation on the two nodes
                // where the 0 and 1 paths from the current node merge at the next step of
                // the trellis.
                //
                // The code polynomials are assumed to have 1's on both ends. Given a
                // function encode_state() that returns the two symbols for a given
                // encoder state in the low two bits, such a code will have the following
                // identities for even 'n' less than 64:
                //
                // encode_state(n) = encode_state(n+65)
                // encode_state(n+1) = encode_state(n+64) = (3 ^ encode_state(n))
                //
                // Any convolutional code you would actually want to use will have
                // these properties, so these assumptions aren't too limiting.
                //
                // Doing this as a macro lets the compiler evaluate at compile time the
                // many expressions that depend on the loop index and encoder state and
                // emit them as immediate arguments.
                // This makes an enormous difference on register-starved machines such
                // as the Intel x86 family where evaluating these expressions at runtime
                // would spill over into memory.
                viterbidecoder_result decode(uint8_t deinterleavedFrame[DEINTERLEAVER_FRAME_LENGTH], bool isHardDecision);
            private:
                // 8-bit parity lookup table, generated by partab.c
                uint8_t partab[256] = {
                        0, 1, 1, 0, 1, 0, 0, 1,   1, 0, 0, 1, 0, 1, 1, 0,   1, 0, 0, 1, 0, 1, 1, 0,   0, 1, 1, 0, 1, 0, 0, 1,
                        1, 0, 0, 1, 0, 1, 1, 0,   0, 1, 1, 0, 1, 0, 0, 1,   0, 1, 1, 0, 1, 0, 0, 1,   1, 0, 0, 1, 0, 1, 1, 0,
                        1, 0, 0, 1, 0, 1, 1, 0,   0, 1, 1, 0, 1, 0, 0, 1,   0, 1, 1, 0, 1, 0, 0, 1,   1, 0, 0, 1, 0, 1, 1, 0,
                        0, 1, 1, 0, 1, 0, 0, 1,   1, 0, 0, 1, 0, 1, 1, 0,   1, 0, 0, 1, 0, 1, 1, 0,   0, 1, 1, 0, 1, 0, 0, 1,
                        1, 0, 0, 1, 0, 1, 1, 0,   0, 1, 1, 0, 1, 0, 0, 1,   0, 1, 1, 0, 1, 0, 0, 1,   1, 0, 0, 1, 0, 1, 1, 0,
                        0, 1, 1, 0, 1, 0, 0, 1,   1, 0, 0, 1, 0, 1, 1, 0,   1, 0, 0, 1, 0, 1, 1, 0,   0, 1, 1, 0, 1, 0, 0, 1,
                        0, 1, 1, 0, 1, 0, 0, 1,   1, 0, 0, 1, 0, 1, 1, 0,   1, 0, 0, 1, 0, 1, 1, 0,   0, 1, 1, 0, 1, 0, 0, 1,
                        1, 0, 0, 1, 0, 1, 1, 0,   0, 1, 1, 0, 1, 0, 0, 1,   0, 1, 1, 0, 1, 0, 0, 1,   1, 0, 0, 1, 0, 1, 1, 0
                };
                // Index into partab
                uint8_t partabIdx[32] = {
                    0, 1, 3, 2, 3, 2, 0, 1, 0, 1, 3, 2, 3, 2, 0, 1, 2, 3, 1, 0, 1, 0, 2, 3, 2, 3, 1, 0, 1, 0, 2, 3
                };
                const int offset = 128;
                int mettab[2][256];
                double sqrt2 = sqrt(2);
                double log2 = log(2);
                // The path memory for each state is 32 bits. This is slightly shorter
                // than we'd like for K=7, especially since we chain back every 8 bits.
                // But it fits so nicely into a 32-bit machine word...
                struct State {
                    ulong path; // Decoded path to this state
                    long metric; // Cumulative metric to this state
                };
        };


        class Descrambler {
            /*
             * Bibliography:
             *
             *   Calcutt, D. M., & Tetley, L. (2004). Satellite communications: principles and applications. Oxford: Elsevier.
             *   https://www.amazon.com/Satellite-Communications-Applications-David-Calcutt/dp/034061448X
             *
             *   Nera. (2015) Nera Inmarsat-C Service Manual. Billingstad: Nera ASA.
             *   https://www.manualslib.com/manual/1201514/Nera-Inmarsat-C.html
             *
             *   Nera. (2015) Nera Saturn C Technical Manual. Billingstad: Nera ASA.
             *   ftp://161.87.213.193.static.cust.telenor.com/Manuals/Saturn%20C/SatC_Marine_Tech_Manual_A.pdf
             *
             *
             * Scrambling prevents 0 or 1 from continuing excessively; if 0 or 1 continues, clock
             * recovery would be reduced at the BPSK modulator. For scrambling, the output is
             * gained by inputting to the output from the scramble generator and the modulo-2 adder.
             *
             * Descrambling is the reverse of scrambling.
             * The 640 bytes frame is split into 160 groups, each of 4 consecutive bytes.
             * Each group either has all its bytes inverted, or is left as is depending of the
             * 0 or 1 output of the descramble generator, G = X3 + X4 + X5 + X7
             *
             *Instead of recalculating on the fly, we populate a descrambler array with the output of the
             *descramble generator.
             *
             */
            public:
                #define DESCRAMBLER_GROUP_COUNT 160
                #define DESCRAMBLER_FRAME_LENGTH 640
                struct descrambler_result {
                    int length;
                    uint8_t descramblerFrame[DESCRAMBLER_FRAME_LENGTH];
                    int frameNumber;
                    bool isBadBulletinBoard;
                    unsigned long timestamp;
                };
                Descrambler();
                descrambler_result decode(uint8_t viterbiFrame[VITERBIDECODER_FRAME_LENGTH]);
            private:
                uint8_t descramblerArray[DESCRAMBLER_GROUP_COUNT];
                // bit 7 = bit 0
                // bit 6 = bit 1
                // bit 5 = bit 2
                // bit 4 = bit 3
                uint8_t invertBits(uint8_t input);
        };

        class INMARSATC_EXPORT Decoder {
            public:
                Decoder(int tolerance);
                std::vector<Descrambler::descrambler_result> decode(uint8_t inputBits[DEMODULATOR_SYMBOLSPERCHUNK]);
            private:
                UWFinder* uwFinder;
                Depermuter* depermuter;
                Deinterleaver* deinterleaver;
                ViterbiDecoder* viterbiDecoder;
                Descrambler* descrambler;
        };
    }
}

#endif // INMARSATC_H
